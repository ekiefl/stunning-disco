import ConfigParser
import numpy as np
import shutil
import pandas as pd
import os
import variants.snv as snv
import glob
import random

import pymol
from pymol import cmd

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

    def __init__(self, table, args):

        self.args = args
        A = lambda x: args[x] if x in args else None

    #   making input variables attributes
        self.input_dir = A("input-dir")
        self.method_list = A("method-list")
        self.ss3_confidence = A("ss3-confidence")
        self.ss8_confidence = A("ss8-confidence")
        self.solvent_acc_confidence = A("solvent-acc-confidence")

    #   make sure table is VariantsTable object
        self.table = table
        MoleculeOperations.validate_table(self.table)
    #   make sure the input_dir is in good shape
        MoleculeOperations.validate_input_dir(self.input_dir, self.table)
    #   define a list of all append methods in this class
        self.all_append_methods = self.get_append_methods()

    #   validate append methods if methods_list provided
        if self.method_list:
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
        for method in self.methods_list:
            append_method = getattr(self, method)
            append_method()


    def return_new_table(self):
        """
        Once the new columns have been added, this method should be called by the 
        VariantsTable object.
        """
        return self.table.saav_table



    def get_append_methods(self):
        """
        This function defines a list of all the append methods. It searches for all the methods
        in this class that start with "append_"
        """
        append_methods = [func for func in dir(self) \
                          if callable(getattr(self, func)) and func.startswith("append")] 
        return append_methods


    def append_8_class_structure(self):
        """
        The 8-class secondary structure predictions for RaptorX are found in
        self.input_dir/{}.all_in_one/labeling/{}.ss8. This function uses this
        information from every gene in self.table.saav_table to construct 9 new
        columns:

        ss8 : ["H", "G", "I", "E", "B", "T", "S", "L", "U"]
            The most likely secondary structure. "H" = alpha helix, "E" = beta
            sheet, "C" = loop, ..., and "U" = unknown. We classify the SAAV as
            "U" whenever the highest confidence score attributed to the 3
            structures is less than self.ss8_confidence.
        ss8_genewide_X : integer
            The number of amino acids (not SAAVs--you can count this yourself
            using the table) in the protein that are classified as X
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
        
        saav_table_grouped = self.table.saav_table.groupby("corresponding_gene_call")
        self.table.saav_table = saav_table_grouped.apply(calc_ss8_data_in_1gene)


    def append_3_class_structure(self):
        """ 
        The 3-class secondary structure predictions for RaptorX are found in
        self.input_dir/{}.all_in_one/labeling/{}.ss3.  This function uses this
        information from every gene in self.table.saav_table.genes to construct
        4 new columns:

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

        table_grouped = self.table.saav_table.groupby("corresponding_gene_call")
        self.table.saav_table = table_grouped.apply(calc_ss3_data_in_1gene)


    def append_solvent_accessibility(self):
        """
        This function incorporates the solvent accesibility predictions found
        in self.input_dir/{}.all_in_one/labeling/{}.acc by adding one column,
        "solvent_acc".

        solvent_acc : ["B"=Buried(pACC=1-10), "M"=Medium(pACC=11-40),
        "E"=Exposed(pACC=41-100), "U"=Unknown] pACC is equal to the relative
        solvent accessibility calculated by DSSP. If the highest confidence for
        the classifications B, M, and E is less than
        self.solvent_acc_confidence, the SAAV is considered U.
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

        saav_table_grouped = self.table.saav_table.groupby("corresponding_gene_call")
        self.table.saav_table = saav_table_grouped.apply(calc_solv_acc_in_1gene)

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
    #   save the table in same format style as the output of gen-variability-profile
        self.saav_table.to_csv(self.table_out_fname, sep='\t', index=False)

#==================================================================================================

class MoleculeOperations():
    
    def __init__(self,  args):
        
    #   define inputs as attributes
        self.args = args
        A = lambda x: args[x] if x in args else None
        self.output_dir = A("output-dir")
        self.input_dir = A("input-dir")
        self.color_vars = A("color-vars")
        self.pymol_config_fname = A("pymol-config")
        self.saav_table_fname = A("saav-table")

    #   initialize a VariantsTable object
        self.table = VariantsTable(args)

    #   make sure the input_dir is in good shape
        self.validate_input_dir(self.input_dir, self.table)

    #   load pymol_config file as ConfigParse object and validate its logic
        """ For now, all I do is load the settings. Later I will make sure all
        of the logic works out. data type in each column (e.g. sidechain is
        Boolean), radius and transparency should deal with scalar data, if
        group_by should be qualitative. """
        self.load_and_validate_pymol_config_file()
        
    #   create columns if they don't exist in saav-table already
        """ For now, I will assume any columns mentioned in self.pymol_config
        already exist in self.table. Later I will have to make a list of
        columns found in pymol_config, and tag to each a class and
        corresponding method responsible for appending them to saav_table. Then
        I will call each of those classes with a method_list parameter, e.g.
        AddRaptorXProperty(self.table, methods_list, args). Each Add class should
        have an attribute called "columns_i_add" """

    #   create the output directory if it doesn't already exist
        self.mkdirp(self.output_dir)
        self.mkdirp(os.path.join(self.output_dir, os.path.splitext(self.saav_table_fname)[0]))

    #   loop_through_sections calls loop_through_genes which calls loop_through_groups
        self.loop_through_sections()


    def loop_through_sections(self):
        """
        """
    #   create pse files for all sections in pymol_config
        for section in self.pymol_config.sections():

        #   create folder for section if doesn't exist
            self.section_dir = os.path.join(self.output_dir, os.path.splitext(self.saav_table_fname)[0], section)
            self.mkdirp(self.section_dir)

        #   This is where the "color_hierarchy" attribute comes in. colormap is
        #   the colormap used to color the groups and it is redefined based on
        #   the color_hierarchy attribute. If color_hierarchy = global, then
        #   colormap is defined only once and therefore the colors are
        #   conserved across genes. If color_hierarchy = gene, then colormap is
        #   redefined for every gene. If color_hierarchy = group, then colormap
        #   is redefined for every group. To avoid colormap being redefined
        #   inappropriately, self.get_colormap is always called conditioned on
        #   self.colored = False
            self.color_hierarchy = self.pymol_config.get(section, "color_hierarchy")
            if self.color_hierarchy == "global":
                self.colormap = self.get_colormap()

        #   loop through each gene
            self.loop_through_genes(section)


    def loop_through_genes(self, section):
    
        for gene in self.table.genes:
        
        #   make directory for the gene
            self.mkdirp(os.path.join(self.section_dir, gene))
        #   make the protein .pse 
            self.do_all_protein_pse_things(gene)
        #   color/recolor if appropriate
            if self.color_hierarchy == "gene":
                self.colormap = self.get_colormap()
        #   loop through each group
            self.loop_through_groups(section, gene)

    def loop_through_groups(self, section, gene):

        groups = self.get_group_list(self, section)
        for group in groups:
        
        #   color/recolor if appropriate
            if self.color_hierarchy == "group":
                self.colormap = self.get_colormap()

        #   
        


    def get_group_list(self, gene_saav_table, section):
        """
        Returns a list of the unique elements in the SAAV table column specified
        by the "group_by" attribute in self.pymol_config_fname.
        """
        return list(saav_table[self.pymol_config.get(section, "group_by")].unique())


    def get_color_mapping(self, table):
        
        table
        pass


    def do_all_protein_pse_things(self, gene):
        """
        This function creates a directory for the gene if it doesn't exist, locates
        the .pdb file, loads in in a pymol session, applies all programmed settings,
        and saves the file.
        """
    #   get pdb_file
        pdb_file = self.get_pdb_file(gene)

    #   loads pdb file into pymol, applies default settings, and saves
        self.create_protein_pse_file(gene, pdb_file)

    def create_protein_pse_file(self, gene, pdb_file):

    #   load and create scaffold and surface objects
        cmd.reinitialize()
        cmd.load(pdb_file, "scaffold")
        cmd.copy("surface","scaffold")
        cmd.hide() # creates blank slate to work with

    #   All settings related to the protein .pse should be set here
        cmd.bg_color("white")
        cmd.set("sphere_scale", "2")
        cmd.set("fog", "off")
        cmd.color("wheat", "scaffold")
        cmd.set("ray_trace_mode", "0")
        cmd.set("ray_opaque_background", "off")
        cmd.color("gray90", "surface")
        cmd.orient()
        cmd.show("cartoon", "scaffold")

    #   save it in the gene subfolder
        cmd.save(os.path.join(self.section_dir, gene, "00_{}.pse".format(gene)))


    def load_and_validate_pymol_config_file(self):
        """
        The format of this file should be standard INI format. an example would be

            [<unique_name_1>]
            # a descriptive name for <unique_name_1> should be chosen, but could be simple, like "1"
            color            = <a column name from SAAV table, default = red>
            radius           = <a column name from SAAV table, default = 2>
            transparency     = <a column name from SAAV table, default = 1>
            sidechain        = <True, default = False>
            group_by         = <a column name from SAAV table>
            color_hierarchy  = <global, gene, group>

            [unique_name_2]
            #group_by is the only required argument
            group_by = <only requirement

            [DEFAULT]
            # if any settings are missing from the above sections, they are given the
            # default values defined here

        """

        if not os.path.isfile(self.pymol_config_fname):
            raise("{} isn't even a file".format(self.pymol_config_name))

    #   this is a list of all possible attributes
        attributes_list   = ["color","radius","transparency","sidechain","color_hierarchy","group_by"]
        required_defaults = ["color","radius","transparency","sidechain","color_hierarchy"]

    #   load file
        self.pymol_config = ConfigParser.ConfigParser()
        self.pymol_config.read(self.pymol_config_fname)

        """ IMPORTANT: The indexing syntax for configparser is fundamentally
        different between python2 and python3. Once we start running pymol
        through python 3, this syntax will have to be imported over. """

    #   make sure all attributes are present in DEFAULT section
        default_attributes = [x[0] for x in self.pymol_config.items("DEFAULT")]
        if set(default_attributes) != set(required_defaults):
            raise ValueError("DEFAULT must have exactly these attributes: {}".format(attributes_list))

        for section in self.pymol_config.sections():

        #   don't allow whitespace in section names
            if " " in section:
                raise ValueError("Please no whitespace in section names.")

            for name, value in self.pymol_config.items(section):
            #   demand all user attributes are in attributes_list
                if name not in attributes_list:
                    raise ValueError("{} in {} is not a valid attribute.".format(name, section))


    def get_pdb_file(self, gene):
        """
        Returns the pdb file for a given gene.
        """ 
        pdb_files = glob.glob(os.path.join(self.input_dir, "{}.all_in_one".format(gene), "*.pdb"))
        if len(pdb_files) != 1:
            raise ValueError("Expecting 1 pdb file but found {}".format(len(pdb_files)))
        pdb_file = pdb_files[0]
        return pdb_file


    def mkdirp(self, path):
        """
        This function makes a directory if it doesn't exist, otherwise it
        doesn't do anything.  It is a python wrapper for the "mkdir -p" command
        in bash
        """
        if os.path.exists(path):
            pass
        else:
            os.mkdir(path)


    """
    Below are a couple of variable validation methods that I've made static so
    they can be borrwed by other classes.
    """

    @staticmethod
    def validate_table(table):
        class_name = table.__class__.__name__
        if not class_name == "VariantsTable":
            raise ValueError("You have passed an object of class '{}' to the variable\
                              table. It must be an object of class VariantsTable".\
                              format(class_name))

    @staticmethod
    def validate_input_dir(input_dir, table):
        """
        checks that the directory exists and that all the gene names match the
        gene names in the saav table
        """
    #   check that folder exists
        if not os.path.isdir(input_dir):
            raise ValueError("you're so bad at this.")
 
    #   if there are zero folders matching the gene_id2.all_in_one format, raise hell
        subdir_list = glob.glob(os.path.join(input_dir, "*.all_in_one"))
        if len(subdir_list) == 0:
            raise ValueError("what the fuck man")
 
    #   if genes in the table aren't found in the raptorx folder, no
        raptor_genes = [int(os.path.splitext(os.path.basename(x))[0]) for x in subdir_list]
        in_both = [gene for gene in table.genes if gene in raptor_genes]
        if not table.genes == in_both:
            missing_in_raptorx = [gene for gene in table.genes if gene not in in_both]
            raise ValueError("You have genes in your table that are missing in your raptorX \
            structure repository. Here are some that are missing: {}".format(missing_in_raptorx))



