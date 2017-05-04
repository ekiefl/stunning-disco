#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def filter_by_snv_pos(df,
                      req_all = False,
                      at_least = None,
                      ret_params = False):    
    """
    Filters a SNV table according to snv position properties. This is different
    from SNV.sample_by_snv. For example, SNV.filter_by_snv can be used to 
    filter out all SNVs (across samples and positions) that have a coverage
    of less than 500. However, SNV.filter_by_snv_pos filters according to 
    SNV positions. For example, SNV.filter_by_snv_pos could be used to filter
    out any SNV positions where a SNV is not present in all samples. So if you 
    have 100 samples, and at the SNV position defined by pos_in_contig = 242,
    only 99 of the samples display an SNV at this position, this SNV position
    will be removed from the table.
    
    INPUTS
    
    df:
        SNV table as pandas DataFrame object 
        
    req_all:
        if True, only SNV positions present in all samples will pass the 
        filter.
        
    at_least:
        if req_all = False, at_least defines the minimum number of samples
        that a SNV must be present in for the SNV position to pass the filter
        
    RETURNS
    
    df:
        SNV table pandas DataFrame object filtered by SNV position
    params:
        dictionary of the filtering parameters used (if they are not None),
        where the key is the string of the parameter, e.g. "mmc" and the value
        is the parameter passed to filter_by_sample.
    """
    
#   saves params used as dictionary
    params = locals()
    del params['df']
    
#   pre-filtering stats
    print("\npre-filter          \t:\t{} entries".\
    format(np.shape(df)[0]))
    print("                      \t:\t{} samples".\
    format(len(df["sample_id"].unique())))
    
#   filters SNV positions based on req_all and at_least
    if req_all:
        num_samples = len(df["sample_id"].unique())
        bysnvpos = df.groupby("pos_in_contig")        
        count = bysnvpos["sample_id"].count()
        filtered_ids = count[count == num_samples].index.values    
        df = df[df["pos_in_contig"].isin(filtered_ids)]
        print("filter req_all      \t:\t{} entries".\
        format(np.shape(df)[0]))
        print("                    \t:\t{} samples\n".\
        format(len(df["sample_id"].unique())))    
                
    else:
        if at_least != None:
            bysnvpos = df.groupby("pos_in_contig")        
            count = bysnvpos["sample_id"].count()
            filtered_ids = count[count >= at_least].index.values    
            df = df[df["pos_in_contig"].isin(filtered_ids)]
            print("filter at_least      \t:\t{} entries".\
            format(np.shape(df)[0]))
            print("                     \t:\t{} samples\n".\
            format(len(df["sample_id"].unique())))   
            
    if ret_params:    
        return df, params
    else:
        return df
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    