#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def snv_density_by_gene(df,boxplot=True,cohort=False):
    """
    Calculates SNV density on a per gene per sample basis

    INPUTS    
    
    df:
        SNV table (with gene_length column appended)
    boxplot:
        plots boxplot if True
    cohort : Boolean
        if True, boxplot is plotted for each cohort (see splice_sample_id(...)).
        boxplot must also be True.
        
    RETURNS
    
    df:
        pandas Dataframe object containing SNV density info
    """
    
    if any("gene_length" == c for c in df.columns.values.tolist()) == False:
        raise ValueError("gene_length column does not exist. See SNV.append_gene_len(...)")
    
    df2 = df.groupby(["corresponding_gene_call","sample_id"]).size().reset_index()
    df = pd.merge(df,df2)
    df["density_per_gene_per_sample"] = df[0].astype(float)/df["gene_length"]

    if boxplot:
        if cohort:
            sns.boxplot(data=df,x='corresponding_gene_call',y='density_per_gene_per_sample',hue='cohort')
            #df.groupby("cohort").boxplot(column='density_per_gene_per_sample',by='corresponding_gene_call')
        else:
            df.boxplot(column='density_per_gene_per_sample',by='corresponding_gene_call')
    return df
