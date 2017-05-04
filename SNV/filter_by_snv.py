#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def filter_by_snv(df,
                  min_cov=None,
                  min_n2n1=None,
                  min_depcon=None,
                  filt_nongene=True,
                  filt_incgene=True,
                  filt_quince=True,
                  ret_params=False):
    """
    Filters a SNV table according to individual SNV properties. For example,
    this function can be used to filter out all SNVs with coverage less than
    500. For sample-based filtering, see filter_by_sample.
    
    INPUTS
        
    df:
        SNV table pandas DataFrame object
    min_cov:
        all SNVs with coverage less than min_cov are filtered
    min_n2n1:
        all n2n1 ratios less than min_n2n1 are filtered
    min_depcon:
        all SNVs with a "departure from consensus" of less than min_depcon are
        filtered. This may be useful when looking at a single metagenome and 
        you want to filter out insignificant variation.
    filt_nongene:
        filters out all SNVs in non- or unidentified-genes
    filt_incgene:
        filters out all SNVs in incomplete genes
    filt_quince:
        if --quincemode flag is provided for anvi-gen-variability-profile, some 
        entries will have zero variation and are removed here (they have
        departure_from_consensus = 0 even if flag -j is used in 
        anvi-gen-variability-profile)
        
    RETURNS
    
    df:
        SNV table pandas DataFrame object filtered by SNVs
    params:
        dictionary of the filtering parameters used (if they are not None),
        where the key is the string of the parameter, e.g. "mmc" and the value
        is the parameter passed to filter_by_sample.
    """
    
#   saves params used as dictionary
    params = locals()
    del params['df']
    
#   pre-filtering stats
    print(("\npre-filter          \t:\t{} entries".\
    format(np.shape(df)[0])))
    print(("                      \t:\t{} samples".\
    format(len(df["sample_id"].unique()))))
            
#   filters out entries with no ambiguity
    if filt_quince:
        df = df[df['departure_from_consensus']>0]
        print(("filter quince       \t:\t{} entries".\
        format(np.shape(df)[0])))
        print(("                    \t:\t{} samples".\
        format(len(df["sample_id"].unique()))))

#   filters out low coverage positions across all samples
    if min_cov != None:
        df = df[df['coverage']>=min_cov]
        print(("filter coverage < {}\t:\t{} entries".\
        format(min_cov,np.shape(df)[0])))
        print(("                    \t:\t{} samples".\
        format(len(df["sample_id"].unique()))))

#   filters out positions that fall outside of known gene calls
    if filt_nongene:
        df = df[df['corresponding_gene_call']!=-1]
        print(("filter non-gene calls\t:\t{} entries".\
        format(np.shape(df)[0])))
        print(("                     \t:\t{} samples".\
        format(len(df["sample_id"].unique()))))

#   filters out incomplete genes
    if filt_incgene:
        df = df[df['in_complete_gene_call']==1]
        print(("filter incomplete genes\t:\t{} entries".\
        format(np.shape(df)[0])))
        print(("                       \t:\t{} samples".\
        format(len(df["sample_id"].unique()))))
    
#   filters out n2n1_ratios too low
    if min_n2n1 != None:
        df = df[df['n2n1ratio']>=min_n2n1]
        print(("filter n2n1 ratio < {}\t:\t{} entries".\
        format(min_n2n1, np.shape(df)[0])))
        print(("                      \t:\t{} samples".\
        format(len(df["sample_id"].unique()))))
        
#	filters out departure from consensus values that are too low
    if min_depcon != None:
        df = df[df['departure_from_consensus']>=min_depcon]
        print(("filter dep con < {}\t:\t{} entries".\
        format(min_depcon, np.shape(df)[0])))
        print(("                   \t:\t{} samples\n".\
        format(len(df["sample_id"].unique()))))
        
    if ret_params:    
        return df, params
    else:
        return df