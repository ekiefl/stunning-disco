import numpy as np
import os
import anvio
import anvio.terminal as terminal
from anvio import constants
import pandas as pd
from anvio.errors import ConfigError
import anvio.filesnpaths as filesnpaths

pp = terminal.pretty_print
progress = terminal.Progress()
run = terminal.Run(width=62)
pd.options.display.max_rows = 100
pd.options.display.max_columns = 500


class SubstitutionMatrix:
    def __init__(self, args):
        self.args = args
        A = lambda x, t: t(args.__dict__[x]) if x in args.__dict__ else None
        null = lambda x: x
        self.output = A("output", null)
        self.min_cov = A("min_cov", null)
        self.genes_only = A('genes_only', null)
        self.variant_type = A("variant_type", null)
        self.variability_table = A("variability_table", null)
        self.samples_of_interest = A('samples_of_interest', null)
        self.filter_by_coverage_outliers = A("filter_by_coverage_outliers", null)
        self.filter_by_coverage_outliers_std = A("filter_by_coverage_outliers_std", null)

        self.variant_items = {"nt": constants.nucleotides, "codon": constants.codons, "aa": constants.amino_acids}
        self.items = sorted(self.variant_items[self.variant_type])

        if not self.output:
            raise ConfigError(" you must supply an output")

        if not self.filter_by_coverage_outliers:
            run.warning(" you didn't supply the file to filter by coverage outliers. Probably you should to avoid nonspecific mapping")

        self.master = pd.read_csv(args.variability_table, sep="\t")


    def filter(self):

        self.master = pd.read_csv(self.variability_table, sep="\t")

        # gets rid of samples not in our samples of interest
        if self.samples_of_interest:
            samples_of_interest = [n.strip() for n in open(self.samples_of_interest).readlines()]
            self.master = self.master[self.master["sample_id"].isin(samples_of_interest)]

        # gets rid of variation not lying in genes
        if self.genes_only:
            self.master = self.master[self.master["corresponding_gene_call"] >= 0]

        # Gets rid of all positions in which all samples do not have at least 20X
        if self.min_cov:
            groupby_pos = self.master.groupby("unique_pos_identifier")
            def func(x):
                if np.all(x["coverage"] >= 20):
                    return x
                else:
                    return
            self.master = self.master.groupby("unique_pos_identifier").apply(func).reset_index(drop=True)

        # filters according to coverage lying outside the non-outlier average
        if self.filter_by_coverage_outliers:
            # add non-outlier coverages and stds
            nonoutlier_cov = pd.read_csv(self.filter_by_coverage_outliers, sep="\t")
            nonoutlier_cov = pd.melt(nonoutlier_cov, id_vars="gene_callers_id", var_name="sample_id", value_name="nonoutlier_gene_cov")

            nonoutlier_std = pd.read_csv(self.filter_by_coverage_outliers_std, sep="\t")
            nonoutlier_std = pd.melt(nonoutlier_std, id_vars="gene_callers_id", var_name="sample_id", value_name="nonoutlier_gene_std")

            nonoutlier = pd.merge(nonoutlier_cov, nonoutlier_std)
            nonoutlier["corresponding_gene_call"] = nonoutlier["gene_callers_id"]
            self.master = pd.merge(self.master, nonoutlier)

            #Filter if the coverage value is outside one standard deviation
            self.master = self.master[self.master.nonoutlier_gene_cov != 0]
            self.master = self.master[(self.master.coverage > (self.master.nonoutlier_gene_cov - 1*self.master.nonoutlier_gene_std)) & \
                            (self.master.coverage < (self.master.nonoutlier_gene_cov + 1*self.master.nonoutlier_gene_std))]


    def collapse_samples(self):
        not_coverage_cols = [x for x in self.master.columns if x not in ["coverage","sample_id"] + self.items]
        def func(x):
            y = x.loc[:, ["coverage"] + self.items].sum(axis=0)
            y = y.append(x.iloc[0,:][not_coverage_cols])
            return y
        self.master = self.master.groupby("unique_pos_identifier").apply(func).reset_index(drop=True)


    def calculate_BLOSUM(self):
        '''
        I'm using amino acid substitution matrices from  protein blocks,  that is, the original
        BLOSUM paper.  since coverage varies from amino acid to amino acid,  each amino acid is
        considered its own block within the language used in the paper.  this means w = 1 the number
        of blocks equals the number of  unique position identifiers within the variability table.
        https://www.ncbi.nlm.nih.gov/pmc/articles/PMC50453/pdf/pnas01096-0363.pdf, and a good
        explanation can be found here:
        http://www.cs.columbia.edu/4761/assignments/assignment1/reference1.pdf

        Q is the matrix of qij elements in the first equation of the paper
        P is the matrix of pij elements in the second equation of the paper
        '''

        # for debugging
        #self.items = ["A", "S"]
        #self.master = pd.DataFrame(np.array([[9, 1], [0, 0]]), columns=self.items, index=self.items)

        # CALCUTING Q
        # ===========

        for first in self.items:
            for second in self.items:
                # this condition ensures no double counting, i.e j <= i in equation 1
                if self.items.index(second) > self.items.index(first):
                    self.master[first+second] = 0
                    continue
                if first == second:
                    self.master[first+second] = self.master[first]*(self.master[first]-1)/2.0
                else:
                    self.master[first+second] = self.master[first]*self.master[second]

        # The frequency table has redundancy, i.e. accounts for AB and BA
        frequency_table = np.asarray(self.master.loc[:, self.items[0]+self.items[0]:self.items[-1]+self.items[-1]].sum(axis=0))

        # just what equation 1 is
        self.Q = frequency_table.reshape((len(self.items),len(self.items)))
        self.Q = self.Q / np.sum(self.Q)

        # CALCUTING P
        # ===========

        # this one's a thinker
        self.P = np.sum(self.Q, axis=1)/2 + np.sum(self.Q, axis=0)/2

        # CALCUTING E
        # ===========

        # this one's even more of a thinker
        outer_product = np.outer(self.P, self.P)
        self.E = np.tril(2*outer_product - np.identity(len(self.P)) * self.P**2)

        # CALCUTING S
        # ===========

        self.S = np.round(2 * np.log2(self.Q / self.E))
        for i in range(len(self.items)):
            for j in range(i):
                self.S[j, i] = self.S[i, j]
        self.S = pd.DataFrame(self.S, columns=self.items, index=self.items)


    def save_BLOSUM(self):
        if os.path.exists(self.output):
            pass
        else:
            os.mkdir(self.output)
        self.S.to_csv(os.path.join(self.output, "BLOSUM.txt"), sep="\t")


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="")

    parser.add_argument('variability_table', help="variability table")
    parser.add_argument('-c', '--min-cov', type=float, required=False)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-t', '--variant-type', required=True, help="either aa or nt")
    parser.add_argument('-s', '--samples-of-interest', required=False, help="file of samples of interest")
    parser.add_argument('-g', '--genes-only', action='store_true', help="uses only positions lying within a gene call")
    parser.add_argument('-u', '--filter-by-coverage-outliers', required=False, help="filter coverages deemed outliers from anvi-summarize. Give outliers txt file per gene")
    parser.add_argument('-w', '--filter-by-coverage-outliers-std', required=False, help="filter coverages deemed outliers from anvi-summarize. Give outliers std txt file per gene")

    args = parser.parse_args()

    try:
        SSM = SubstitutionMatrix(args)
        SSM.filter()
        SSM.collapse_samples()
        SSM.calculate_BLOSUM()
        SSM.save_BLOSUM()

    except ConfigError as e:
        print(e)
        sys.exit(-1)







