import numpy as np
import shutil
import pandas as pd
import os
import variants.snv as snv
import glob
import random
import utils as ut

#import pymol
#from pymol import cmd

class AddRaptorXProperty():
    """
    This class holds methods for adding columns to the SAAV table that are
    specifically related to the output of RaptorX Structure Prediction:
    http://raptorx.uchicago.edu/. An input directory is required for the
    instantiation of this class, and should be formatted as follows:

        path/to/input_dir
            gene_id1.all_in_one
                ...
            gene_id2.all_in_one
                ...
            gene_id3.all_in_one
                ...
            ...
    Each folder gene_id#.all_in_one is the standard output from each of the 
    RaptorX job submissions.

    If you want to add another method that adds columns to the SAAV table, you
    should add it here and it should start with "append_"
    """

    def __init__(self, table, input_dir, method_list=None):

        self.table = table
        self.input_dir = input_dir
    #   make sure table is VariantsTable object
        self.validate_table()
    #   make sure the input_dir is in good shape
        self.validate_input_dir()
    #   define a list of all append methods in this class
        self.all_append_methods = self.get_append_methods()

    #   validate append_methods if methods_list provided
        if method_list:
            self.methods_list = methods_list
            self.validate_method_list()
    #   otherwise run all append methods
        else:
            self.methods_list = self.all_append_methods

    #   add all the columns associated with the append methods in method_list
        self.add_columns()
    

    def add_columns(self):
        """
        This function calls all the append methods in self.methods_list. If any of the
        methods have optional parameters, they should be passed in the line where
        append is defined. For example, append = Append(parameter = 4)
        """
        append = self.Append(self)
        for method in self.methods_list:
            append_method = getattr(append, method)
            append_method()
        self.table.saav_table = append.saav_table

    def return_new_table(self):
        """
        Once the new columns have been added, this method should be called by the 
        VariantsTable object.
        """
        return self.table.saav_table


    def validate_table(self):
        class_name = self.table.__class__.__name__
        if not class_name == "VariantsTable":
            raise ValueError("You have passed an object of class '{}' to the variable \
                              table. It must be an object of class VariantsTable".\
                              format(class_name))


    def validate_input_dir(self):
        """
        checks that the directory exists and that all the gene names match the
        gene names in the saav table
        """
    #   check that folder exists
        if not os.path.isdir(self.input_dir):
            raise ValueError("you're so bad at this.")

    #   if there are zero folders matching the gene_id2.all_in_one format, raise hell
        subdir_list = glob.glob(os.path.join(self.input_dir, "*.all_in_one"))
        if len(subdir_list) == 0:
            raise ValueError("what the fuck man")

    #   if genes in the table aren't found in the raptorx folder, no
        raptor_genes = [int(os.path.splitext(os.path.basename(x))[0]) for x in subdir_list]
        in_both = [gene for gene in self.table.genes if gene in raptor_genes]
        if not self.table.genes == in_both:
            missing_in_raptorx = [gene for gene in self.table.genes if gene not in in_both]
            raise ValueError("You have genes in your table that are missing in your raptorX \
                              structure repository. Here are some that are missing: {}".\
                              format(missing_in_raptorx))


    def get_append_methods(self):
        """
        This function defines a list of all the append methods. It searches for all the methods
        in the embedded Append class.
        """
        append_methods = [func for func in dir(self.Append) \
                          if callable(getattr(self.Append, func)) and "__" not in func] 
        return append_methods


    class Append():
        """
        This simple class is simply a repository for different methods used by AddRaptorXProperty
        """
        def __init__(self, RaptorX, ss3_confidence=0.66, ss8_confidence=0.5, solvent_acc_confidence=0.66):
            """
            RaptorX is an instance of the AddRaptorXProperty class.
            """
            self.saav_table = RaptorX.table.saav_table
            self.input_dir = RaptorX.input_dir

        #   parameters for all methods in class
            self.ss3_confidence = ss3_confidence
            self.ss8_confidence = ss8_confidence
            self.solvent_acc_confidence = solvent_acc_confidence

        def append_8_class_structure(self):
            """
            RaptorX outputs a {}.all_in_one folder for each sequence it predicts structure for. All
            of these need to be placed into a folder self.input_dir. The secondary structure predictions
            for RaptorX are found in self.input_dir/{}.all_in_one/labeling/{}.ss8. This function uses
            this information from every gene in self.saav_table to construct 9 new columns:

            ss8 : ["H", "G", "I", "E", "B", "T", "S", "L", "U"]
                The most likely secondary structure. "H" = alpha helix, "E" = beta sheet, "C" = loop, ..., 
                and "U" = unknown. We classify the SAAV as "U" whenever the highest confidence score 
                attributed to the 3 structures is less than self.ss8_confidence.
            ss8_genewide_X : integer
                The number of amino acids (not SAAVs--you can count this yourself using the table) in the
                protein that are classified as X
            """
            
            columns = ("codon_order_in_gene","AA","ss8","prob_H","prob_G","prob_I","prob_E","prob_B","prob_T","prob_S","prob_L")

            def calc_ss8_data_in_1gene(x):
            #   get gene id of this groupby object, then find path of .ss8 file for that gene
                gene = x["corresponding_gene_call"].values[0]
                ss8_path = glob.glob(os.path.join(self.input_dir,"{}.all_in_one".format(gene),"labeling","*.ss8"))[0]
            #   load ss8 data for gene as pandas DataFrame
                ss8 = pd.read_csv(ss8_path, skiprows=3, names=columns, delim_whitespace=True)
            #   0 index the raptor results to conform to anvio convention :\
                ss8["codon_order_in_gene"] -= 1
            #   add gene column (used with 'reference' to uniquely map ss8 entries to saav entries)
                ss8["corresponding_gene_call"] = gene
            #   if the highest confidence for the 8 classes < self.prob_confidence, the SAAV is classified as "U" for unknown
                l = [col for col in ss8.columns if "prob_" in col]
                ss8["ss8"] = ss8.apply(lambda row: row["ss8"] if any(row[l] > self.ss8_confidence) else "U", axis = 1)
            #   ss8_genewide_X is the total number of AAs in the gene with secondary structure X
                ss8["ss8_genewide_H"] = len(ss8[ss8["ss8"] == "H"])
                ss8["ss8_genewide_G"] = len(ss8[ss8["ss8"] == "G"])
                ss8["ss8_genewide_I"] = len(ss8[ss8["ss8"] == "I"])
                ss8["ss8_genewide_E"] = len(ss8[ss8["ss8"] == "E"])
                ss8["ss8_genewide_B"] = len(ss8[ss8["ss8"] == "B"])
                ss8["ss8_genewide_T"] = len(ss8[ss8["ss8"] == "T"])
                ss8["ss8_genewide_S"] = len(ss8[ss8["ss8"] == "S"])
                ss8["ss8_genewide_L"] = len(ss8[ss8["ss8"] == "L"])
                ss8["ss8_genewide_U"] = len(ss8[ss8["ss8"] == "U"])
            #   The `AA` column (1-letter AA code) is redundant, we already have `reference` (3-letter AA code)
            #   ss8_X are also dropped to cut down on saav_table_size, but code could be modified to retain them
                ss8 = ss8.drop(["AA"]+l, axis=1)
            #   merge with original dataframe
                return pd.merge(x,ss8)
            
            saav_table_grouped = self.saav_table.groupby("corresponding_gene_call")
            self.saav_table = saav_table_grouped.apply(calc_ss8_data_in_1gene)


        def append_3_class_structure(self):
            """ 
            The 3-class secondary structure predictions for
            RaptorX are found in self.input_dir/{}.all_in_one/labeling/{}.ss3.
            This function uses this information from every gene in self.saav_table.genes
            to construct 4 new columns:

            ss3 : ["H", "E", "C", "U"] The most likely secondary structure. "H" =
            alpha helix, "E" = beta sheet, "C" = loop, and "U" = unknown. We
            classify the SAAV as "U" whenever the highest confidence score
            attributed to the 3 structures is less than self.prob_confidence.

            prob_genewide_X : integer The number of amino acids (not SAAVs--you can
            count this yourself using the table) in the protein that are classified
            as X.  
            """
            
            columns = ("codon_order_in_gene","AA","ss3","prob_H","prob_E","prob_C")

            def calc_ss3_data_in_1gene(x):
            #   get gene id of this groupby object, then find path of .ss3 file for that gene
                gene = x["corresponding_gene_call"].values[0]
                ss3_path = glob.glob(os.path.join(self.input_dir,"{}.all_in_one".format(gene),"labeling","*.ss3"))[0]
            #   load ss3 data for gene as pandas DataFrame
                ss3 = pd.read_csv(ss3_path, skiprows=3, names=columns, delim_whitespace=True)
            #   0 index the raptor results to conform to anvio convention :\
                ss3["codon_order_in_gene"] -= 1
            #   add gene column (used with 'reference' to uniquely map ss3 entries to saav entries)
                ss3["corresponding_gene_call"] = gene
            #   if the highest confidence for the 3 classes < self.ss3_confidence, the SAAV is classified as "U" for unknown
                l = [col for col in ss3.columns if "prob_" in col]
                ss3["ss3"] = ss3.apply(lambda row: row["ss3"] if any(row[l] > self.ss3_confidence) else "U", axis = 1)
            #   prob_genewide_X is the total number of AAs in the gene with secondary structure X
                ss3["ss3_genewide_H"] = len(ss3[ss3["ss3"] == "H"])
                ss3["ss3_genewide_E"] = len(ss3[ss3["ss3"] == "E"])
                ss3["ss3_genewide_C"] = len(ss3[ss3["ss3"] == "C"])
                ss3["ss3_genewide_U"] = len(ss3[ss3["ss3"] == "U"])
            #   The `AA` column (1-letter AA code) is redundant, we already have `reference` (3-letter AA code)
            #   prob_X are also dropped to cut down on table_size, but code could be modified to retain them
                ss3 = ss3.drop(["AA"]+l, axis=1)
            #   merge with original dataframe
                return pd.merge(x,ss3)

            table_grouped = self.saav_table.groupby("corresponding_gene_call")
            self.saav_table = table_grouped.apply(calc_ss3_data_in_1gene)


        def append_solvent_accessibility(self):
            """
            RaptorX outputs a {}.all_in_one folder for each sequence it predicts structure for. All
            of these need to be placed into a folder self.input_dir. This function incorporates the
            solvent accesibility predictions found in self.input_dir/{}.all_in_one/labeling/{}.acc
            by adding one column, "solvent_acc".

            solvent_acc : ["B"=Buried(pACC=1-10), "M"=Medium(pACC=11-40), "E"=Exposed(pACC=41-100), "U"=Unknown]
                pACC is equal to the relative solvent accessibility calculated by DSSP. If the highest confidence
                for the classifications B, M, and E is less than self.solvent_acc_confidence,
                the SAAV is considered U.
            """

            columns = ("codon_order_in_gene","AA","solvent_acc","prob_B","prob_M","prob_E")

            def calc_solv_acc_in_1gene(x):
            #   get gene id of this groupby object, then find path of .ss3 file for that gene
                gene = x["corresponding_gene_call"].values[0]
                solvent_acc_path = glob.glob(os.path.join(self.input_dir,"{}.all_in_one".format(gene),"labeling","*.acc"))[0]
            #   load solvent_acc data for gene as pandas DataFrame
                acc = pd.read_csv(solvent_acc_path, skiprows=3, names=columns, delim_whitespace=True)
            #   0 index the raptor results to conform to anvio convention :\
                acc["codon_order_in_gene"] -= 1
            #   add gene column (used to uniquely map acc entries to saav entries)
                acc["corresponding_gene_call"] = gene
            #   if the highest confidence for the 3 classes < self.solvent_acc_confidence, the SAAV is classified as "U" for unknown
                l = [col for col in acc.columns if "prob_" in col]
                acc["solvent_acc"] = acc.apply(lambda row: row["solvent_acc"] if any(row[l] > self.solvent_acc_confidence) else "U", axis = 1)
            #   prob_genewide_X is the total number of AAs in the gene with secondary structure X
                acc["solvent_acc_genewide_B"] = len(acc[acc["solvent_acc"] == "B"])
                acc["solvent_acc_genewide_M"] = len(acc[acc["solvent_acc"] == "M"])
                acc["solvent_acc_genewide_E"] = len(acc[acc["solvent_acc"] == "E"])
                acc["solvent_acc_genewide_U"] = len(acc[acc["solvent_acc"] == "U"])
            #   The `AA` column (1-letter AA code) is redundant, we already have `reference` (3-letter AA code)
            #   prob_X are also dropped to cut down on saav_table_size, but code could be modified to retain them
                acc = acc.drop(["AA"]+l, axis=1)
            #   merge with original dataframe
                return pd.merge(x,acc)

            saav_table_grouped = self.saav_table.groupby("corresponding_gene_call")
            self.saav_table = saav_table_grouped.apply(calc_solv_acc_in_1gene)

#==================================================================================================

class VariantsTable():
    
    def __init__(self, args):
        
        self.args = args
        A = lambda x: args[x] if x in args else None

    #   converting input to self variables
        self.saav_table_fname = A("saav-table")
        self.genes_list_fname = A("gene-list")
        self.samples_list_fname = A("samples-list")
        self.simplify_sample_id_method = A("simplify-sample-id-method")
        self.output_dir = A("output-dir")
        self.input_dir = A("raptor-repo")

    #   load the saav table
        self.load()

    #   simplify sample id 
        if self.simplify_sample_id_method:
            self.simplify_sample_ids()

    #   get genes and samples lists
        self.genes = self.load_genes_file_as_list()
        self.samples = self.load_samples_file_as_list()

    #   filter by genes and samples
        self.filter_table("sample_id", self.samples)
        self.filter_table("corresponding_gene_call", self.genes)

    #   get array of columns that exist in this table
        self.columns = self.get_columns()

    def get_columns(self):
        """
        Defines an array of the columns in self.saav_table
        """
        return self.saav_table.columns.values

    def load_genes_file_as_list(self):
        """
        You really don't know?
        """
    #   if there is no file name provided, return all genes in SAAV table
        if self.genes_list_fname is None:
            return list(self.saav_table["corresponding_gene_call"].unique())
    #   otherwise return genes_list
        else:
            if not os.path.isfile(self.genes_list_fname):
                raise ValueError("Your genes-list path sucks")
            return [int(x.strip()) for x in open(self.genes_list_fname).readlines()]

    def load_samples_file_as_list(self):
        """
        You really don't know?
        """
    #   if there is no file name provided, return all genes in SAAV table
        if not self.samples_list_fname:
            return list(self.saav_table["sample_id"].unique())
    #   otherwise return samples_list
        else:
            if not os.path.isfile(self.samples_list_fname):
                raise ValueError("Your samples-list path sucks")
            return [int(x.strip()) for x in open(self.samples_list_fname).readlines()]

    def filter_table(self, column, elements):
        """
        This function filters the SAAV table according to single column criteria,
        where column is the column you are filtering by and elements is a list of 
        accepted elements in that column. All other rows are filtered
        """
        self.saav_table = self.saav_table[self.saav_table[column].isin(elements)]

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
            raise ValueError("fuck you, man")

    def load(self):
        """
        Load the SAAV table output from gen-variability-profile.
        """
        if not self.saav_table_fname:
            raise ValueError("fuck you... you ass")
        self.saav_table = pd.read_csv(self.saav_table_fname, sep='\t', header=0, index_col=False)

    def save(self):
        """
        Saves any operations performed on the table.
        """
    #   sanity checks
        if not self.table_out_fname:
            raise ValueError("fuck you")
        if os.path.isfile(self.table_out_fname):
            raise ValueError("fuck you")
    #   save the table in same format as the output of gen-variability-profile
        self.saav_table.to_csv(self.table_out_fname, sep='\t', index=False)

#class MoleculeOperations():
#    
#    def __init__(self, args):
#
#    #   transform input variables into `self` variables
#        self.args = args
#        A = lambda x: args[x] if x in args else None
#        self.coloring_hierarchy = A("coloring-hierarchy")
#        self.saav_table_fname = A("saav-table")
#        self.genes_list_fname = A("gene-list")
#        self.samples_list_fname = A("samples-list")
#        self.output_dir = A("output-dir")
#        self.input_dir = A("raptor-repo")
#        self.simplify_sample_id_method = A("simplify-sample-id-method")
#
#
#
#        self.saav_table = Table(args)
#
#    #   loads and cleans up saav table
#        self.saav_table = snv.load_snv_table(self.saav_table_fname)
#        self.clean_up_saav_table()
#    
#    #   create the output directory if it doesn't already exist
#        ut.mkdirp(self.output_dir)
#
#    #   get sample names that have saavs in this protein
#        self.samples = list(self.saav_table["sample_id"].unique())
#
#        # load gene list (list for which raptorx solved structures for)
#        self.genes_of_interest = np.loadtxt(self.genes_list_fname).astype(int)
#
#    def clean_up_saav_table(self):
#        """
#        This function performs any required or desired touch ups on the SAAV table. It is
#        not for adding column names, but could be use for little stuff, like simplifying
#        sample IDs and removing entries with invalid entries.
#        """
#    #   must filter quince mode entries (otherwise practically everywhere is a SAAV)
#        self.saav_table = snv.filter_by_snv(self.saav_table, filt_quince=True, filt_incgene=False)
#    #   entries in the SAAV table are rarely N/A (like BLOSUM scores) and are removed here
#        self.saav_table = self.saav_table.dropna()
#    #   simplifies sample_id according to self.simplify_sample_ids() if simplify-sample-id-methods is provided
#        if self.simplify_sample_id_method:
#            self.simplify_sample_ids()
#
#    def simplify_sample_ids(self):
#        """
#        If ad hoc sample_id modifications need to be made for whatever reason, you can define a method
#        here for such a purpose. For an example, see the case when simplify_sample_id_method == "SAR11",
#        for which this method was first defined.
#        """
#        if self.simplify_sample_id_method == "SAR11":
#            self.saav_table["sample_id"] = self.saav_table["sample_id"].str.replace("SAR11_","")
#            self.saav_table["sample_id"] = self.saav_table["sample_id"].str.replace("_BOWTIE2","")
#            self.saav_table["cohort"] = self.saav_table["sample_id"].str.split("_",expand=True)[0]
#        elif self.simplify_sample_id_method == "placeholder":
#            pass
#        else:
#            raise ValueError("simplify_sample_ids subroutine for simplify_sample_id_method = {} not defined".\
#                              format(self.simplify_sample_id_method))
#
#    def main(self):
#    
#    #   loop through each gene
#        for gene in self.genes_of_interest:
#            print ("##################################### GENE {}".format(gene))
#    
#    #       create subfolder for the gene if it doesn't exist
#            gene_dir = os.path.join(self.output_dir, str(gene))
#            ut.mkdirp(gene_dir)
#    
#            pdb_files = glob.glob(os.path.join(self.input_dir, "{}.all_in_one".format(gene), "*.pdb"))
#
#            if len(pdb_files) != 1:
#                raise ValueError("Expecting 1 pdb file but found {}".format(len(pdb_files)))
#            pdb_file = pdb_files[0]
#
#            colormap = self.generate_colormap(gene, pdb_file)
#    
#    #       create PSE file for new gene
#            self.create_protein_pse_file(gene, pdb_file)
#    
#    #       invoke the visualization routine, which makes a pse file for each sample
#            self.visualize_routine(gene, pdb_file, colormap)
#
#    def create_protein_pse_file(self, gene, pdb_file):
#    
#    #   load and create scaffold and surface objects
#        cmd.reinitialize()
#        cmd.load(pdb_file,"scaffold")
#        cmd.hide() # hides default sticks
#        cmd.copy("surface","scaffold")
#    
#    #   specify common features between genes and saavs
#        self.set_common_pse_attributes(kind="protein")
#    
#    #   save it in the gene subfolder
#        cmd.save(os.path.join(self.output_dir, str(gene), "00_{}.pse".format(gene)))
#
#    def visualize_routine(self, gene, pdb_file, colormap):
#
#    #   subset the saav table to only include gene 
#        gene_saav_table = self.saav_table[self.saav_table["corresponding_gene_call"]==gene]
#    
#        for sample in self.samples:
#    #       subset to only include on sample from the cohort
#    
#            gene_saav_table_sample = gene_saav_table[gene_saav_table["sample_id"]==sample]
#    
#    #       create saav_data to be passed to visualization subroutine
#            codon = list(gene_saav_table_sample["codon_order_in_gene"].values)
#            codon_score = list(gene_saav_table_sample["competing_aas"].values)
#            saav_data = dict(zip(codon, codon_score))
#    
#            if len(saav_data) == 0:
#                # this creates an empty file if there are no SAAVs
#                cmd.save(os.path.join(self.output_dir, str(gene), "{}.pse".format(sample)))
#            else:
#    #           start visualization procedure
#                self.create_saav_pse_file(gene, sample, saav_data, colormap, pdb_file)
#    
#    ###########################################################################
#    
#    
#    def create_saav_pse_file(self, gene, sample_name, saav_data, colormap, pdb_file):
#    
#    #   we have to load the pdb so we know where the amino acids are in space
#        cmd.load(pdb_file)
#        cmd.hide()
#    
#    #   set commonalities between pse files
#        self.set_common_pse_attributes(kind="saav")
#    
#    #   create the colors for each SAAV
#        pymol.saav_colors = {}
#        for saav in saav_data.keys():
#            pymol.saav_colors[str(saav)] = colormap[saav_data[saav]]
#    
#    #   define the selection of the SAAVs
#        sites = "+".join([str(x) for x in list(saav_data.keys())])
#        cmd.select(sample_name+"_sel","resi {} & name ca".format(sites))
#    
#    #   make the SAAV selection its own object
#        cmd.create(sample_name,sample_name+"_sel")
#    
#    #   change the color for each saav according to the saav_colors dict
#        cmd.alter(sample_name,"color = pymol.saav_colors.get(resi, color)")
#        cmd.alter(sample_name+"_sel","color = pymol.saav_colors.get(resi, color)")
#    
#    #   displays the spheres
#        cmd.show("spheres",sample_name)
#    
#        cmd.save(os.path.join(self.output_dir, str(gene), "{}.pse".format(sample_name)), sample_name)
#        cmd.reinitialize()
#    
#    ###########################################################################
#    
#    def set_common_pse_attributes(self, kind=None):
#        """
#        This function sets commond parameters I want shared between PSEs
#        
#        INPUT
#        -----
#        kind : str
#            accepted inputs: "saav" or "protein". Some of the settings are specific
#            to whether we are talking about a protein or a collection of saavs and 
#            so the type must be specified.
#        """
#
#        if not kind:
#            raise ValueError("must be kind")
#    
#    #   make it pretty initial settings (can edit later interactively)
#        cmd.bg_color("white")
#        cmd.set("sphere_scale","2")
#        cmd.set("fog","off")
#    
#        if kind == "saav":
#            pass
#        if kind == "protein":
#            cmd.set(name="transparency",
#                    value=0.60)
#            cmd.show('cartoon','scaffold')
#            #cmd.show('surface','surface')
#            cmd.color('wheat','scaffold')
#            cmd.set("ray_trace_mode","0")
#            cmd.set("ray_opaque_background","off")
#            cmd.color("gray90","surface")
#    #       orient to a 'best position', then save position to return to visualize_routine
#            cmd.orient()
#    
#    ###########################################################################
#    
#    def generate_colormap(self, gene, pdb_file):
#        """
#        Gets color mapping for codon_score scores of each saav. Since each gene will
#        procur different amino acid substitutions, the color range is defined
#        elative to each gene.
#        
#        INPUT
#        -----
#        pdb_file : str
#            the protein structure file generated by raptorX
#        OUTPUT
#        ------
#        colormap : dict
#            a dictionary with the keys being codon_score values and the values
#            being color indices for pymol
#        """
#
#        gene_saav_table = self.saav_table[self.saav_table["corresponding_gene_call"]==gene]   
#
#        cmd.reinitialize()
#        cmd.load(pdb_file)
#    #   There are 15 possible codon_score values. temp is a 15 atom selection
#    #   which the cmd.spectrum command is applied to in order to get a 
#    #   pretty spectrum of colors for each codon_score value. Actually, since
#    #   it would be very rare to see all 15 mutations, I make the length
#    #   whatever the range of data is. For example, if the minimum and
#    #   maximum BLOSUM scores are -3 and 4, temp is length 7.
#    
#        values = gene_saav_table["competing_aas"].unique()
#
#        minimum = 0
#        maximum = len(values) + 1
#        values = random.sample(values, len(values))
#    
#        color_range = range(minimum, maximum)
#    
#        cmd.select("temp","resi 1-{} & name ca".format(maximum-minimum))
#        cmd.spectrum("count","gcbmry","temp")
#    #   codon_score_colors is a list of color indices for the colors
#        pymol.codon_score_colors = []
#        cmd.iterate("temp","codon_score_colors.append(color)")
#    #   colormap is a dictionary. The words are codon_score values
#    #   and the definitions are color indices
#        colormap = dict(zip(values, pymol.codon_score_colors))
#    
#        cmd.reinitialize()
#    
#        return colormap
#    
