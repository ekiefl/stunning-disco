#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def snv_density_by_gene(df,
					  boxplot=True):
    """
    Calculates SNV density on a per gene per sample basis

    INPUTS    
    
    df:
        SNV table with gene_length column appended.
    boxplot:
        plots boxplot if True
        
    RETURNS
    
    snv_dens:
        pandas Dataframe object containing SNV density info
    """
    
    if any("gene_length" == c for c in df.columns.values.tolist()) == False:
        raise ValueError("gene_length column does not exist. See SNV.append_gene_len(...)")
    
    snv_dens = df.groupby(['gene_length','corresponding_gene_call','sample_id'])\
               ['pos_in_contig'].aggregate(len)
    snv_dens = snv_dens.reset_index(level=['gene_length','corresponding_gene_call','sample_id'])
    snv_dens.columns = snv_dens.columns.str.replace("pos_in_contig","num_snvs")
    snv_dens['density'] = snv_dens.num_snvs / snv_dens.gene_length
    
    if boxplot == True:
        snv_dens.boxplot(column='density',by='corresponding_gene_call')
#        snv.dens.boxplot(column='density',by='sample_id')
        
    return snv_dens