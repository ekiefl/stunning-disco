import ConfigParser
import time
import random
import numpy as np
import shutil
import pandas as pd
import os
import variants.snv as snv
import glob
import random
import matplotlib.cm

import pymol
from pymol import cmd
import sys; print(sys.version)

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
        self.save_file = A("save-file") 
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
            print("appending method {}".format(method))
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
        self.save_file = A("save-file")

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
        if not self.save_file:
            raise ValueError("fuck you")
        if os.path.isfile(self.save_file):
            raise ValueError("fuck you")
    #   save the table in same format style as the output of gen-variability-profile
        self.saav_table.to_csv(self.save_file, sep='\t', index=False)

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
        self.no_images = A("no-images")
        self.ray = A("ray")
        self.res = A("res")

    #   initialize a VariantsTable object
        self.table = VariantsTable(args)

    #   make sure the input_dir is in good shape
        self.validate_input_dir(self.input_dir, self.table)

    #   load pymol_config file as ConfigParse object and validate its logic
        """ For now, all I do is load the settings. Later I will make sure all
        of the logic works out. data type in each column (e.g. sidechain is
        Boolean), radii and alpha should deal with scalar data, if
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
        creates pse files for all sections in pymol_config
        """
        for section in self.pymol_config.sections():

        #   create folder for section if doesn't exist
            self.section = section
            self.section_dir = os.path.join(self.output_dir, os.path.splitext(self.saav_table_fname)[0], section)
            self.mkdirp(self.section_dir)
        #   within the section, create two subdirectories: 1 for PyMOL 1 for images
            self.section_dir_pymol = os.path.join(self.section_dir, "PyMOL")
            if not self.no_images:
                self.section_dir_images = os.path.join(self.section_dir, "Images")
                self.mkdirp(self.section_dir_images)
            self.mkdirp(self.section_dir_pymol)

        #   This is where the "color_hierarchy" attribute comes in. colormap is
        #   the colormap used to color the groups and it is redefined based on
        #   the color_hierarchy attribute. If color_hierarchy = global, then
        #   colormap is defined only once (all genes use the same colormap) If
        #   color_hierarchy = gene, then colormap is redefined for every gene.
        #   Finally, if color_hierarchy = group, then colormap is redefined for
        #   every group. The last would be useful if you you're coloring by
        #   competing_aas and there are more than, say, 40 AASs per gene. To
        #   avoid colormap being redefined inappropriately, self.get_colormap
        #   is always called conditioned by color_hierarchy being either
        #   global, gene, or group
            self.color_hierarchy = self.pymol_config.get(section, "color_hierarchy")
            if self.color_hierarchy == "global":
                print("I JUST ENTERED SECTION COLORING")
                saav_table_subset = self.get_relevant_saav_table()
                self.colorObject = Color(saav_table_subset, self.pymol_config, self.section)
                self.colorObject.export_legend(self.section_dir_pymol, "{}_legend.txt".format(section))
        #   loop through each gene
            self.loop_through_genes()


    def loop_through_genes(self):
    
        for gene in self.table.genes:

            print("Currently on SECTION: {}, GENE: {}".format(self.section, gene))
        
        #   get access to protein_pdb
            self.protein_pdb_path = self.get_protein_pdb(gene)
        #   make pymol directory for the gene
            self.gene_dir_pymol = os.path.join(self.section_dir_pymol, str(gene))
            self.mkdirp(self.gene_dir_pymol)
        #   make image directory for the gene
            if not self.no_images:
                self.gene_dir_images = os.path.join(self.section_dir_images, str(gene))
                self.mkdirp(self.gene_dir_images)

        #   make the protein .pse 
            self.do_all_protein_pse_things(gene)
        #   color/recolor if appropriate
            if self.color_hierarchy == "gene":
                print("I JUST ENTERED GENE COLORING")
                saav_table_subset = self.get_relevant_saav_table(gene=gene)
                self.colorObject = Color(saav_table_subset, self.pymol_config, self.section, gene=gene)
                self.colorObject.export_legend(self.gene_dir_pymol, "00_{}_legend.txt".format(gene))
        #   loop through each group
            self.loop_through_groups(gene)


    def loop_through_groups(self, gene):
        """
        """
        groups = self.get_group_list(self, self.section)
        for group in groups:

            print("GROUP: {}".format(group))

        #   subset the saav_table to include only group and gene
            saav_table_subset = self.get_relevant_saav_table(gene=gene, group=group)
        #   color/recolor if appropriate
            if self.color_hierarchy == "group":
                print("I JUST ENTERED GROUP COLORING")
                self.colorObject = Color(saav_table_subset, self.pymol_config, self.section, gene=gene, group=group)
                self.colorObject.export_legend(os.path.join(self.gene_dir_pymol),"{}_legend.txt".format(group))

        #   make the saav .pse
            self.do_all_saav_pse_things(saav_table_subset, gene, group)
        #   make an image for the group
            self.do_all_image_things(gene, group)


    def do_all_image_things(self, gene, group):
        """
        """
    #   make directory for the group 
        group_dir = os.path.join(self.gene_dir_images, group)
        self.mkdirp(group_dir)

    #   merge the pses
        pse_list = [self.protein_pse_path, self.saav_pse_path]
        self.join_pses(pse_list)
    #   I can't believe this, but the alpha is messed up unless I do this
    #   AFTER merging
        cmd.set("sphere_transparency", 0.0)
        cmd.set("sphere_transparency", 1.0)

    #   I could have sophisticated routines here for taking multiple images,
    #   collages, any view setting like orientation, etc. this really deserves
    #   its own class. Instead I'll just call this 1-liner that saves the image
        self.create_image_file(gene, group, group_dir) 
        

    def create_image_file(self, gene, group, group_dir):
        """
        saves a single png image
        """
        if self.ray:
        #   very costly!
            cmd.ray(self.res)

        save_path = os.path.join(group_dir,"{}.png".format(group))
        cmd.png(save_path)
        cmd.reinitialize()


    def join_pses(self, pse_list):
        """
        Merges multiple pse sessions using the settings of the first file.
    
        INPUT
        -----
        pse_list : list
            list of paths of the pse's to merge. The first one in the list 
            is the one the settings are matched to.
        save : str, default None
            If provided, the merged pse is saved with the filepath `save`
        """
        cmd.load(pse_list[0])
        for pse in pse_list[1:]:
            cmd.load(pse,partial=1)


    def do_all_saav_pse_things(self, saav_table_subset, gene, group):
        """
        This function creates creates a saav .pse file for each group, and then
        saves the file.
        """
    #   holds all the info for each SAAV like color, radii, & alpha
        self.saav_properties = self.fill_saav_properties_table(saav_table_subset)

    #   create save directory for the saav pse
        self.saav_pse_path = os.path.join(self.section_dir_pymol, str(gene), "{}.pse".format(group))

    #   create the file
        self.create_saav_pse_file(gene, group)

    def create_saav_pse_file(self, gene, group):
    
        s = time.time()
    #   set any settings specific to the SAAV .pse files here. Once a SAAV .pse
    #   is merged with its corresponding protein .pse, the settings of the
    #   merged .pse inherits the settings defined under create_protein_pse_file
        cmd.bg_color("white")
        cmd.set("fog","off")

    #   if there are no SAAVs,just save an empty file
        if len(self.saav_properties.index) == 0:
            cmd.save(os.path.join(self.section_dir_pymol, str(gene), "{}.pse".format(group)), group)
            cmd.reinitialize()

    #   otherwise do the stuff we were planning to
        else:
        #   we have to load the pdb so we know where the amino acids are in space
            cmd.load(self.protein_pdb_path)

        #   define the selection of the SAAVs and create their object
            sites = "+".join([str(resi) for resi in self.saav_properties.index])
            cmd.select(group+"_sel","resi {}".format(sites))
            cmd.create(group,group+"_sel")

        #   delete the protein
            cmd.delete(os.path.splitext(os.path.basename(self.protein_pdb_path))[0])
        
            pymol.saav_properties = self.saav_properties
        #   change the color for each saav according to the saav_colors dict
            cmd.alter(group,"color = pymol.saav_properties.loc[int(resi),'color_indices']")
            cmd.alter(group,"s.sphere_transparency = pymol.saav_properties.loc[int(resi),'alpha']")
            cmd.alter(group,"s.sphere_scale = pymol.saav_properties.loc[int(resi),'radii']")
            cmd.rebuild()

        #   displays the spheres
            cmd.hide()
            cmd.show("spheres","name ca")

        #   save the file
            cmd.save(self.saav_pse_path, group)
            print( time.time() - s)
            cmd.quit()
            #cmd.reinitialize()


    def fill_saav_properties_table(self, saav_table_subset):
        """
        This function produces a table of SAAV properties (size, color, radii,
        maybe more later) that, as an example, looks like this:

        resi   color_indices   alpha          radii
        148    5342            0.8            2.0
        244    2042            0.7            2.0
        248    4222            0.8            8.0
        ...    ...             ...            ...

        This is the table fed to PyMOL's cmd.alter function in order to change
        the appearance of the SAAVs in the SAAV .pse file. Getting this table
        requires some data massaging. For one, if your group is composed of
        multiple samples, its possible to have multiple SAAVs at a single codon
        position, and their columns which specify color_indices, alpha,
        and/or radii could be non-identical. To reconcile this, I take the
        following approach: columns that have number data are averaged (for
        example alpha and radii must be from number data), and for
        columns that have string data, the most frequent entry is the one
        chosen (e.g. if competing_aas is your color column and you observe
        AspGlu AspGlu and IleGlue at resi 148, the color_indices are calculated
        for AspGlu.
        """

        """ IMPORTANT: Currently, this only works when a color column is
        specified. Will need to be refactored to include single pymol colors.
        One way around it would be appending a new column to the saav_table
        when column_color is "False" and calling that the new column_color.
        This would actually work quite well for alpha and radii as well"""

    #   I add +1 to account for the zero-indexing anvio does
        saav_table_subset["resi"] = saav_table_subset["codon_order_in_gene"]+1

    #   rebrand the column names for style points
        color_column = self.colorObject.color_column
        alpha_column = self.pymol_config.get(self.section,"alpha")
        radii_column = self.pymol_config.get(self.section,"radii") 

    #   determine/define whether columns are string-type or number-type
        string_or_number = {color_column : self.colorObject.columntype,
                            alpha_column : "numbers",
                            radii_column : "numbers"}

    #   if number-type, take the mean. if string-type, take most frequent string
        def resolve_ambiguity(dtype):
            if dtype == "numbers":
                return np.mean
            if dtype == "strings":
                return lambda x: x.value_counts().idxmax()

        methods_dictionary  = {}
        for column, dtype in string_or_number.items():
            methods_dictionary[column] = resolve_ambiguity(dtype)

    #   this single line illustrates why pandas is so fucking good
        saav_data = saav_table_subset.groupby("resi").agg(methods_dictionary)

        """ The first part of this method created an aggregated form of the
        data for situations in which the same SAAV position was found multiple
        times within the group. The aggregated data is defined in `saav_data`.
        But pymol doesn't know what to do with data such as `AspGlu`, so the
        second part of this method transforms this data into terms that are
        directly accessible to PyMOL. For example, in place of `AspGlu` would
        be `1452`, the color index corresponding to `AspGlu`. This new data is
        the table spelled out in the above docstring and is called
        `saav_properties`. """

    #   Now I have to translate the data in saav_data to data understood by PyMOL
        saav_properties = pd.DataFrame({}, index=saav_data.index)

    #   add color_indices column
        saav_properties["color_indices"] = self.colorObject.create_color_indices_for_group(saav_data)
        """ IMPORTANT: I will eventually have alpha and radii classes that do a
        similar task to that done directly above for color. But for now, I just
        write specific and non-scalable code. In what follows I assume that the
        columns for radii and alpha are both values that very between 0 and
        1"""
        if (saav_data[alpha_column].max() or saav_data[radii_column].max()) > 1:
            print("alpha: {}".format(saav_data[alpha_column].max()))
            print("radii: {}".format(saav_data[radii_column].max()))
            raise ValueError("Expecting columns to be less than 0")
    #   add alpha column (my def of alpha is: 0 means opaque, 1 means translucent)
        saav_properties["alpha"] = 1 - saav_data[alpha_column]
    #   add radii column
        min_r = 0.5; max_r = 3.0
        saav_properties["radii"] = min_r + (max_r-min_r) * saav_data[radii_column]

        return saav_properties

    def do_all_protein_pse_things(self, gene):
        """
        This function creates a directory for the gene if it doesn't exist, locates
        the .pdb file, loads in in a pymol session, applies all programmed settings,
        and saves the file.
        """
    #   loads pdb file into pymol, applies default settings, and saves
        self.create_protein_pse_file()
    #   the protein is in the PyMOL state. now I save it
        self.protein_pse_path = os.path.join(self.section_dir_pymol, str(gene), "00_{}.pse".format(gene))
    #   save it in the gene subfolder and then reinitialize
        cmd.save(self.protein_pse_path)
        cmd.reinitialize()


    def create_protein_pse_file(self):

    #   load and create scaffold and surface objects
        cmd.reinitialize()
        cmd.load(self.protein_pdb_path, "scaffold")
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
        cmd.set("cartoon_transparency",0.2)
        cmd.set("sphere_transparency", 1.0)
        #cmd.show("surface", "surface")
        #cmd.set("transparency",0.9)
        cmd.show("cartoon", "scaffold")
        cmd.orient()


    def get_relevant_saav_table(self, gene=None, group=None):
        
        saav_table_subset = self.table.saav_table
        if gene:
            saav_table_subset = saav_table_subset[saav_table_subset["corresponding_gene_call"]==gene]
        if group:
            saav_table_subset = saav_table_subset[self.table.saav_table[self.pymol_config.get(self.section, "group_by")]==group]
        return saav_table_subset


    def get_group_list(self, gene_saav_table, section):
        """
        Returns a list of the unique elements in the SAAV table column specified
        by the "group_by" attribute in self.pymol_config_fname.
        """
        return list(self.table.saav_table[self.pymol_config.get(section, "group_by")].unique())


    def load_and_validate_pymol_config_file(self):
        """
        The format of this file should be standard INI format. an example would be

            [<unique_name_1>]
            # a descriptive name for <unique_name_1> should be chosen, but could be simple, like "1"
            color_column     = <a column name from SAAV table, default = False>
            radii           = <a column name from SAAV table, default = 2>
            alpha     = <a column name from SAAV table, default = 1>
            sidechain        = <True, default = False>
            group_by         = <a column name from SAAV table>
            color_hierarchy  = <global, gene, group>
            color_scheme     = <a color_scheme from matplotlib (for now) or a pymol color>

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
        attributes_list   = ["color_column","radii","alpha","sidechain","color_hierarchy","color_scheme","group_by"]
        required_defaults = ["color_column","radii","alpha","sidechain","color_hierarchy","color_scheme"]

    #   load file
        self.pymol_config = ConfigParser.ConfigParser()
        self.pymol_config.read(self.pymol_config_fname)

        """ IMPORTANT: The indexing syntax for configparser is fundamentally
        different between python2 and python3. Once we start running pymol
        through python 3, this syntax will have to be imported over. """

        """ IMPORTANT: make the following code: If color_column = False, ensure
        that color_scheme is a pymol_color."""

        """ IMPORTANT: color_column either takes a string OR the boolean "False". Currently I just
        test whether the string = "False". Ideally I will convert to a boolean if the user input
        is the string "False" """

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


    def get_protein_pdb(self, gene):
        """
        Returns the pdb file for a given gene.
        """ 
        protein_pdbs = glob.glob(os.path.join(self.input_dir, "{}.all_in_one".format(gene), "*.pdb"))
        if len(protein_pdbs) != 1:
            raise ValueError("Expecting 1 pdb file but found {}".format(len(protein_pdbs)))
        protein_pdb = protein_pdbs[0]
        return protein_pdb


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
    they can be borrowed by other classes.
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




class Color():
    
    def __init__(self, saav_table, pymol_config, section, gene=None, group=None):

    #   make input parameters class attributes
        self.saav_table = saav_table
        self.pymol_config = pymol_config
        self.section = section
        self.gene = gene
        self.group = group

    #   get color_scheme from pymol_config
        self.color_scheme = pymol_config.get(section,"color_scheme")
    #   get color_column from pymol_config
        self.color_column = pymol_config.get(section,"color_column")

    #   if the color_column is False, the color index for pymol color is stored and we're done
        if self.color_column == "False":
            self.pymol_colors = True
            self.pymol_color_index = [x[1] for x in cmd.get_color_indices() if x[0]==self.color_scheme][0]
            return
        else:
            self.pymol_colors = False

    #   create template data: an array of the unique entries
        self.template_data_for_colormap = self.extract_unique_color_values()
    #   determine datatype of template_data_for_colormap (is it a string or a number?)
        self.columntype = self.find_columntype()

    #   create colormap
        self.colormap = self.create_colormap()

    #   create a legend
        self.legend = self.create_legend()


    def create_colormap(self):
        """
        returns a colormap for self.color_column
        """
    #   conditional for if the column has string-like data
        if self.columntype == "strings":
        #   number of distinct colors
            colormap = matplotlib.cm.get_cmap(self.color_scheme, len(self.template_data_for_colormap))

    #   conditional for if column has number-like data 
        else:
        #   the colormap is scaled to run from the min to the max of color_map_data
            bounds = [np.min(self.template_data_for_colormap), np.max(self.template_data_for_colormap)]
            colorObject = matplotlib.cm.ScalarMappable(norm=bounds, cmap=self.color_scheme)
            colormap = colorObject.get_cmap()
        return colormap

    def access_colormap(self, value):
        """
        Because either string or number data is accepted for coloring, this
        wrapper function is used to access the RGB values in self.colormap
        depending on self.columntype. The tuple is converted into a list
        and only the first three elements are considered (the fourth is alpha)
        """
    #   string-type data
        if self.columntype == "strings":
            return list(self.colormap(self.legend[value]))[:3]
    #   number-type data
        else:
            return list(self.colormap(value))[:3]


    def create_legend(self):
        """
        creates a dictionary self.legend where each key is a value in
        self.template_data_for_colormap and each value is the corresponding RGB
        held in self.colormap. The self.legend attribute is useful for two
        reasons. The first is for creating a legend for the user. The second
        requires some explaining. The self.colormap function houses 256 RGB
        values that can be accessed by indices or by the value mapped to that
        index. For example, self.colormap(10) returns the 10th RGB value
        whereas self.colormap(10.0) returns the RGB which most closely matches
        the value 10.0. When self.columntype == "numbers" we access RGB values the
        second way, but when self.columntype == "strings" we access RGB values the
        first way. And in order to do that, we need to associate each index to
        a unique value in self.template_data_for_colormap. We do that with
        self.legend.  To spell it out fully, suppose
        self.template_data_for_colormap == [ Ile, Val, Glu ].  Then self.legend
        = {"Ile":0, "Val":1, "Glu":2} and to access the color associated with
        Ile, we call self.colormap(self.legend["Ile"]). The legend only exists
        for when column_name is provided.
        """
        legend = {}
    #   if columntype = strings, a value must be defined for each color
        if self.columntype == "strings":
            for ind, value in enumerate(self.template_data_for_colormap):
            #   index only first 3 since 4th is alpha
                legend[value] = ind

    #   if columntype = numbers, only the min and max need to be defined 
        else:
            mini = self.template_data_for_colormap.min()
            maxi = self.template_data_for_colormap.max()
            midd = (maxi - mini) / 2
            legend[mini] = 0
            legend[midd] = 127
            legend[maxi] = 255
        return legend


    def export_legend(self, path, name=None):
        """
        Exports a text file called color_legend.txt
        """

    #   nothing to report
        if self.pymol_colors:
            return

        if not name:
            name = "color_legend.txt"

    #   creates text file
        text_legend = open(os.path.join(path, name), "w")
        text_legend.write("value\tR\tG\tB\thex\n")
        for key, value in self.legend.items():

        #   get the RGBs as list
            RGB = self.access_colormap(key)
        #   matplotlib uses RGB values bounded by [0,1]. I convert to [0,255]
            RGB = [int(255*X) for X in RGB]
        #   get corresponding hex code
            hexcode = '#%02x%02x%02x' % (RGB[0], RGB[1], RGB[2])
        #   write this and then on to the next--on, on to the next one
            text_legend.write("{}\t{}\t{}\t{}\t{}\n".format(key, RGB[0], RGB[1], RGB[2], hexcode))
        text_legend.close()


    def extract_unique_color_values(self): 
        """
        Takes SAAV table (already assumed to be subsetted according to
        gene and group) and returns the unique values
        """
        template_data_for_colormap = self.saav_table[self.color_column].unique()
        return template_data_for_colormap
        

    def find_columntype(self):
        """
        Determines the data type of the column (either string or number)
        """
    #   for some reason strings come up as type == object in pandas
        return "strings" if self.template_data_for_colormap.dtype==object else "numbers"


    def create_color_indices_for_group(self, saav_data):
        """
        """
        num_saavs = len(saav_data.index)

    #   if a single pymol color was chosen by the user, we're done
        if self.pymol_colors:
            color_indices = [self.pymol_color_index for _ in range(num_saavs)]
            return color_indices

    #   otherwise, we have some work to do

############ AVOID DIRECT EYE CONTACT ############
        """IMPORTANT: PyMOL won't let me define color names that contain any
        numbers whatsoever so I have to name the color of each SAAV some random
        and unique alphabetic string. This sucks and I'm pissed.  Here's this
        piece of shit code that generates a set of random unique alphabetic
        strings. lol."""
        def generate(unique):
            chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
            while True:
                value = "".join(random.choice(chars) for _ in range(5))
                if value not in unique:
                    unique.append(value)
                    break
        color_names = []
        for _ in range(num_saavs):
            generate(color_names)
############ AVOID DIRECT EYE CONTACT ############

        color_indices = []
        for resi in saav_data.index:
            rgb = self.access_colormap(saav_data.loc[resi, self.color_column])
            cmd.set_color(color_names[0], rgb)
            color_index = [x[1] for x in cmd.get_color_indices() if x[0]==color_names[0]][0]
            color_indices.append(color_index)
            color_names.pop(0)
        return color_indices

    @staticmethod
    def is_it_a_pymol_color(color_scheme):
        pymol_colors = cmd.get_color_indices()
        pymol_colors = [x[0] for x in pymol_colors]
        return True if color_scheme in pymol_colors else False

    @staticmethod
    def is_it_a_column_color(color_scheme):
        return True if color_scheme in self.saav_table.columns.values else False

class Timer():

    def __init__(self):
        self.start = time.time()
    def timestamp(self):
        print("{} seconds".format(time.time()-self.start))
