#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import pandas as pd
import variants.snv

def gen_input_for_cluster(df,new_col,sep,file_out,*cols):
    """
    Outputs file used for clustering. The clustering algorithm only uses a 
    single column so information from several columns in concatenated into 
    one. According to Tom, it iis required that the concatenated column is 
    added twice. Also returns SNV table
    
    INPUTS
    
    df:
        SNV table pandas DataFrame object
    new_col:
        name of the new column, e.g. "gene_pos_pair"
    sep:
        value separator in the concatenated column. e.g. if sep='/', a row in 
        the concatenated column would look like "str1/str2/str3..."
    file_out:
        name of output directory. convention for the file name is 
        "cluster_itep.txt" and should be placed in appropriate directory.
    *cols:
        all of the column names to be concatenated. To make it clear, when the
        function is passed, *cols should be passed as a series of column names
        e.g. gen_input_for_cluster(...,"corresponding_gene_call","pos","competing_nts")
    """
    
    df[new_col] = snv.concat_cols(df,sep,*cols)             
    df[['sample_id','sample_id',new_col]].to_csv(file_out.replace(".txt","")+\
                                                 ".txt",sep='\t',
                                                 index=False,
                                                 header=False)
    print(("\ngenerated {}\n".format(file_out)))
