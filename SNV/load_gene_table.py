#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import pandas as pd

def load_gene_table(gene_file):
    """
    Loads the anvi'o generated gene table {userinput}-gene_calls.txt 
    
    INPUTS
    
    gene_file:
        path to gene_table 
        
    RETURNS
    
    gene_table:
        gene table in form of pandas DataFrame object
    """
    gene_table = pd.DataFrame.from_csv(gene_file, sep='\t', header=0)   
    
    print(("gene table {} loaded".format(gene_file)))
    return gene_table