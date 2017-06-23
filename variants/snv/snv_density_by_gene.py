#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np

def snv_density_by_gene(df,boxplot=False,cohort=False,genome_length=None):
    """
    appends SNV density on a per gene per sample basis

    INPUTS    
    
    df:
        SNV table (with gene_length column appended)
    boxplot:
        plots boxplot if True
    cohort : Boolean
        if True, boxplot is plotted for each cohort (see splice_sample_id(...)).
        boxplot must also be True.
    genome_length: Int
        if provided, an additional column is added that's the SNV density on a 
        per sample basis, i.e. across all SNVs. This can be found in the output
        of anvi-summarize (bin_by_bin/{collection}/{collection}-total_length.txt)
        
    RETURNS
    
    df:
        pandas Dataframe object containing SNV density info. More specifically,
        it appends a column 'density_per_gene_per_sample' as well as a 
        'density_per_sample' column if genome_length is provided.
    """
    
    if any("gene_length" == c for c in df.columns.values.tolist()) == False:
        raise ValueError("gene_length column does not exist. See SNV.append_gene_len(...)")
    
    df2 = df.groupby(["corresponding_gene_call","sample_id"]).size().reset_index()
    df = pd.merge(df,df2)
    df["density_per_gene_per_sample"] = df[0].astype(float)/df["gene_length"]
    df = df.drop(0, 1)
    
    if genome_length != None:
        df3 = df.groupby(["sample_id"]).size().reset_index()
        df = pd.merge(df,df3)
        df["density_per_sample"] = df[0].astype(float)/genome_length
        df = df.drop(0, 1)

    if boxplot:
        import matplotlib.pyplot as plt
        import seaborn as sns
        if cohort:
            sns.boxplot(data=df,x='corresponding_gene_call',y='density_per_gene_per_sample',hue='cohort')
            #df.groupby("cohort").boxplot(column='density_per_gene_per_sample',by='corresponding_gene_call')
        else:
            df.boxplot(column='density_per_gene_per_sample',by='corresponding_gene_call')
    return df
