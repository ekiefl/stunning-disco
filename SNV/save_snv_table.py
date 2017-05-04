#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import pandas as pd

def save_snv_table(df,file_out, pkl = False,txt = False):
    """
    saves snv table as .tsv file, .pkl file, or both
    
    INPUTS
    
    df:
        SNV table pandas DataFrame object
    file_out:
        file name. e.g. "path/to/file"
    pkl:
        if True, df is saved as a pickle file
    txt:
        if True, df is saved as a .txt tab-delimited file
    """
    
    file_out = file_out.replace(".txt","").replace(".pkl","")
    
    if (pkl==False) & (txt==False):
        raise ValueError("pkl and txt cannot both be false")
    
    if txt:    
        df.to_csv(file_out+".txt",sep='\t',index=False)
        print(("\nText file {}.txt generated".format(file_out)))
        
    if pkl:
        df.to_pickle(file_out+".pkl") 
        print(("\nPickle file {}.pkl generated".format(file_out)))