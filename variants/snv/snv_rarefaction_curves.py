import variants.snv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import utils as ut

def snv_rarefaction_curves(df, c, genome_wide=True, c_gene=None, genome_len=None):
    """
    Currently only plots. Calculating k-cutoff vector functionality will be added
    when genome_wide = False soon.

    Plots SNV density rarefaction curves.

    Plots SNV density vs sample where samples are ordered from low to high coverage.

    ATTRIBUTES
    ----------
    df : pandas DataFrame
        SNV table
    c : str 
        file path to {collection}-mean_coverage.txt or {collection}-mean_coverage_Q2Q3.txt
        coverage matrix generated by anvi-summarize. Required whether or not
        genome_wide = True.
    genome_wide: Boolean (= True)
        if True, plots SNV density (calculated over whole genome) vs sample--in this 
        case, c and genome_len must be supplied. If False, plots SNV density vs sample for
        each gene, one at a time--in this case, c_gene argument must be supplied.
    c_gene : str (= None)
        file path to {collection}-gene_coverages.txt generated by anvi-summarize when
        --init-gene-coverages flag is supplied. Required if gene_wide = False.
    genome_len : Int (= None)
        Length of genome. Required when genome_wide = True.
    """

    sanity_check(genome_wide,c,c_gene,genome_len)
    if genome_wide:
        df = genome_wide_routine(df,c,genome_len)
    else:
        df = gene_by_gene_routine(df,c,c_gene)
    return df

def gene_by_gene_routine(df,c,c_gene):
#   first, we append the density_per_sample_per_gene column:
    df = snv.snv_density_by_gene(df,boxplot=False,cohort=False)
#   next we append gene cov_per_sample_per_gene
#   gc is a gene table (indices are sample, columns are genes)
    gc = pd.read_csv(c_gene,header=0,index_col=0,sep='\t').transpose()
    gc = ut.melt(gc,rowname="sample_id",colname="corresponding_gene_call",valname="cov_per_gene_per_sample")
    df = df.merge(gc)

#   col_keep is the columns of interest
    col_keep = ["sample_id",
                "cohort",
                "corresponding_gene_call",
                "density_per_gene_per_sample",
                "cov_per_gene_per_sample"]
#   for each gene sorted by number
    for gene in np.sort(df["corresponding_gene_call"].unique()):
        df_temp = df[df["corresponding_gene_call"]==gene][col_keep]
#       drop duplicates (there are as many duplicates as snv's in gene), 
#       then sort by coverage
        df_temp = df_temp.drop_duplicates().sort_values(by="cov_per_gene_per_sample")

        fig = plt.figure()
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        ax2.set_yscale("log")
        df_temp.groupby("cohort").plot(x="sample_id",y="density_per_gene_per_sample",ax=ax1)
        df_temp.groupby("cohort").plot(x="sample_id",y="cov_per_gene_per_sample",ax=ax2)
                                                                                             
        plt.show()
        print(df_temp.head())
    import sys; sys.exit()
#   

def genome_wide_routine(df,c,genome_len):
#   first, we append the density_per_sample column:
    df = snv.snv_density_by_gene(df,boxplot=False,cohort=False,genome_length=genome_len)
#   next we append coverage by gene cov_per_sample
    mc = pd.read_csv(c,header=None, index_col = 0, sep = "\t")
    mc = mc.transpose().reset_index(drop=True)
    cov_name = mc.columns.values[-1]
    mc[cov_name] = mc[cov_name].astype(float)
    mc = mc.rename(columns={"bin":"sample_id",cov_name:"cov_by_sample"})
    df = pd.merge(df,mc)

    df_temp = df.drop_duplicates(subset=["sample_id","density_per_sample"])
    df_temp = df_temp.sort_values(by="cov_by_sample")

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    df_temp.groupby("cohort").plot(x="sample_id",y="density_per_sample",ax=ax1)
    df_temp.groupby("cohort").plot(x="sample_id",y="cov_by_sample",ax=ax2)

    plt.show()

    return df

def sanity_check(genome_wide,c,c_gene,genome_len):
    if c == None:
        raise ValueError("c, the mean coverage matrix file must be supplied.")
    if genome_wide:
        if genome_len == None:
            raise ValueError("genome_len must be provided")
        if c_gene != None:
            print("c_gene is provided but will not be used")
    else:
        if genome_len != None:
            print("genome_len is provided but will not be used")
        if c_gene == None:
            raise ValueError("c_gene must be provided")
