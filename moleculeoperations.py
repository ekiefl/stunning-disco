import numpy as np
import shutil
import pandas as pd
import os
import SNV as snv
import glob
import random
import utils as ut

import pymol
from pymol import cmd

class MoleculeOperations():
    
    def __init__(self, args):

        self.args = args
        #A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        A = lambda x: args[x] if x in args else None
        self.saav_table_fname = A("saav-table")
        self.gene_list_fname = A("gene-list")
        self.output_dir = A("output-dir")
        self.raptor_repo = A("raptor-repo")
        self.simplify_sample_id_method = A("simplify-sample-id-method")

        # loads saav table
        if os.path.isfile("{}_quick_load.pkl".format(self.saav_table_fname)):
            self.saav_table = pd.read_pickle("{}_quick_load.pkl".format(self.saav_table_fname))
        else:
            self.saav_table = snv.load_gene_table(self.saav_table_fname)
            self.saav_table = snv.filter_by_snv(self.saav_table, filt_quince=True, filt_incgene=False)
            self.saav_table = self.saav_table.dropna()
            if self.simplify_sample_id_method:
                self.simplify_sample_ids()
            self.saav_table.to_pickle("{}_quick_load.pkl".format(self.saav_table_fname))

        print(self.saav_table.columns)

        # create the output directory if it doesn't already exist
        ut.mkdirp(self.output_dir)

        # get sample names that have saavs in this protein
        self.samples = list(self.saav_table["sample_id"].unique())

        # load gene list (list for which raptorx solved structures for)
        self.genes = np.loadtxt(self.gene_list_fname).astype(int)

    def simplify_sample_ids(self):
        """
        If ad hoc sample_id modifications need to be made for whatever reason, you can define a method
        here for such a purpose. For an example, see the case when simplify_sample_id_method == "SAR11",
        for which this method was first defined.
        """
        if self.simplify_sample_id_method == "SAR11":
            self.saav_table["sample_id"] = self.saav_table["sample_id"].str.replace("SAR11_","")
            self.saav_table["sample_id"] = self.saav_table["sample_id"].str.replace("_BOWTIE2","")
            self.saav_table["cohort"] = self.saav_table["sample_id"].str.split("_",expand=True)[0]
        elif self.simplify_sample_id_method == "placeholder":
            pass
        else:
            raise ValueError("simplify_sample_ids subroutine for simplify_sample_id_method = {} not defined".\
                              format(self.simplify_sample_id_method))

    def main(self):
    
    #   loop through each gene
        for gene in self.genes:
            print ("##################################### GENE {}".format(gene))
    
    #       create subfolder for the gene if it doesn't exist
            gene_dir = os.path.join(self.output_dir, str(gene))
            ut.mkdirp(gene_dir)
    
            pdb_files = glob.glob(os.path.join(self.raptor_repo, "{}.all_in_one".format(gene), "*.pdb"))

            if len(pdb_files) != 1:
                raise ValueError("Expecting 1 pdb file but found {}".format(len(pdb_files)))
            pdb_file = pdb_files[0]

            colormap = self.generate_colormap(gene, pdb_file)
    
    #       create PSE file for new gene
            self.create_protein_pse_file(gene, pdb_file)
    
    #       invoke the visualization routine, which makes a pse file for each sample
            self.visualize_routine(gene, pdb_file, colormap)

    def create_protein_pse_file(self, gene, pdb_file):
    
    #   load and create scaffold and surface objects
        cmd.reinitialize()
        cmd.load(pdb_file,"scaffold")
        cmd.hide() # hides default sticks
        cmd.copy("surface","scaffold")
    
    #   specify common features between genes and saavs
        self.set_common_pse_attributes(kind="protein")
    
    #   save it in the gene subfolder
        cmd.save(os.path.join(self.output_dir, str(gene), "00_{}.pse".format(gene)))

    def visualize_routine(self, gene, pdb_file, colormap):

    #   subset the saav table to only include gene 
        gene_saav_table = self.saav_table[self.saav_table["corresponding_gene_call"]==gene]
    
        for sample in self.samples:
    #       subset to only include on sample from the cohort
    
            gene_saav_table_sample = gene_saav_table[gene_saav_table["sample_id"]==sample]
    
    #       create saav_data to be passed to visualization subroutine
            codon = list(gene_saav_table_sample["codon_order_in_gene"].values)
            codon_score = list(gene_saav_table_sample["competing_aas"].values)
            saav_data = dict(zip(codon, codon_score))
    
            if len(saav_data) == 0:
                # this creates an empty file if there are no SAAVs
                cmd.save(os.path.join(self.output_dir, str(gene), "{}.pse".format(sample)))
            else:
    #           start visualization procedure
                self.create_saav_pse_file(gene, sample, saav_data, colormap, pdb_file)
    
    ###########################################################################
    
    
    def create_saav_pse_file(self, gene, sample_name, saav_data, colormap, pdb_file):
    
    #   we have to load the pdb so we know where the amino acids are in space
        cmd.load(pdb_file)
        cmd.hide()
    
    #   set commonalities between pse files
        self.set_common_pse_attributes(kind="saav")
    
    #   create the colors for each SAAV
        pymol.saav_colors = {}
        for saav in saav_data.keys():
            pymol.saav_colors[str(saav)] = colormap[saav_data[saav]]
    
    #   define the selection of the SAAVs
        sites = "+".join([str(x) for x in list(saav_data.keys())])
        cmd.select(sample_name+"_sel","resi {} & name ca".format(sites))
    
    #   make the SAAV selection its own object
        cmd.create(sample_name,sample_name+"_sel")
    
    #   change the color for each saav according to the saav_colors dict
        cmd.alter(sample_name,"color = pymol.saav_colors.get(resi, color)")
        cmd.alter(sample_name+"_sel","color = pymol.saav_colors.get(resi, color)")
    
    #   displays the spheres
        cmd.show("spheres",sample_name)
    
        cmd.save(os.path.join(self.output_dir, str(gene), "{}.pse".format(sample_name)), sample_name)
        cmd.reinitialize()
    
    ###########################################################################
    
    def set_common_pse_attributes(self, kind=None):
        """
        This function sets commond parameters I want shared between PSEs
        
        INPUT
        -----
        kind : str
            accepted inputs: "saav" or "protein". Some of the settings are specific
            to whether we are talking about a protein or a collection of saavs and 
            so the type must be specified.
        """

        if not kind:
            raise ValueError("must be kind")
    
    #   make it pretty initial settings (can edit later interactively)
        cmd.bg_color("white")
        cmd.set("sphere_scale","2")
        cmd.set("fog","off")
    
        if kind == "saav":
            pass
        if kind == "protein":
            cmd.set(name="transparency",
                    value=0.60)
            cmd.show('cartoon','scaffold')
            #cmd.show('surface','surface')
            cmd.color('wheat','scaffold')
            cmd.set("ray_trace_mode","0")
            cmd.set("ray_opaque_background","off")
            cmd.color("gray90","surface")
    #       orient to a 'best position', then save position to return to visualize_routine
            cmd.orient()
    
    ###########################################################################
    
    def generate_colormap(self, gene, pdb_file):
        """
        Gets color mapping for codon_score scores of each saav. Since each gene will
        procur different amino acid substitutions, the color range is defined
        elative to each gene.
        
        INPUT
        -----
        pdb_file : str
            the protein structure file generated by raptorX
        OUTPUT
        ------
        colormap : dict
            a dictionary with the keys being codon_score values and the values
            being color indices for pymol
        """

        gene_saav_table = self.saav_table[self.saav_table["corresponding_gene_call"]==gene]   

        cmd.reinitialize()
        cmd.load(pdb_file)
    #   There are 15 possible codon_score values. temp is a 15 atom selection
    #   which the cmd.spectrum command is applied to in order to get a 
    #   pretty spectrum of colors for each codon_score value. Actually, since
    #   it would be very rare to see all 15 mutations, I make the length
    #   whatever the range of data is. For example, if the minimum and
    #   maximum BLOSUM scores are -3 and 4, temp is length 7.
    
        values = gene_saav_table["competing_aas"].unique()

        minimum = 0
        maximum = len(values) + 1
        values = random.sample(values, len(values))
    
        color_range = range(minimum, maximum)
    
        cmd.select("temp","resi 1-{} & name ca".format(maximum-minimum))
        cmd.spectrum("count","gcbmry","temp")
    #   codon_score_colors is a list of color indices for the colors
        pymol.codon_score_colors = []
        cmd.iterate("temp","codon_score_colors.append(color)")
    #   colormap is a dictionary. The words are codon_score values
    #   and the definitions are color indices
        colormap = dict(zip(values, pymol.codon_score_colors))
    
        cmd.reinitialize()
    
        return colormap
    
    ###########################################################################
    
    
    ################################################################################
    
# visualizes SAAVs in the following SAAV table sample by sample 
#gene_by_gene(f="HIMB083-74stations-AA-20x-v15-799genes-dep-10",gene_list="proteins_to_keep_thresh60")

args = {"saav-table":"HIMB083-74stations-AA-20x-v15-799genes-dep-10",
        "gene-list":"00_filter_goi/proteins_to_keep_thresh60",
        "output-dir":"03_pymol",
        "raptor-repo":"02_raptor_out",
        "simplify-sample-id-method":"SAR11"}


MoleculeOperations(args).main()
