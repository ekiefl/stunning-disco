import pandas as pd
import numpy as np

def filter_by_sample(df,
                     mmc = None,
                     ret_params = False):    
    """
    Filters a SNV table according to sample properties. This function offers
    several parameters to filter by, where what is filtered is entire samples
    (if you want to sample individual SNV positions instead, see
    SNV.sample_by_snv). For example, we are typically only interested in 
    samples that recruit reads to the reference genome. Hence, we may want to 
    filter out all samples whose mean coverage at the SNV locations does not 
    exceed some threshold amount. This function carries out this purpose.
    
    INPUTS
    
    df:
        SNV table as pandas DataFrame object 
    mmc:
        minimum mean coverage. all samples with a mean coverage (defined only 
        at SNV positions) less than mmc will be filtered out of the SNV table.
        
    RETURNS
    
    df:
        SNV table pandas DataFrame object filtered by sample
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
    
#   filters out samples with mean coverage less than mmc
    if mmc != None:
#       group df according to column "sample_id"
        bysampleid = df.groupby("sample_id")
#       calculate mean coverage for each sample_id
        mean_cov = bysampleid["coverage"].mean()
#       filter all samples with mean coverages less than mmc
        filtered_ids = mean_cov[mean_cov > mmc].index.values       
        df = df[df["sample_id"].isin(filtered_ids)]   

        print("filter mean coverage\t:\t{} entries".\
        format(np.shape(df)[0]))
        print("                    \t:\t{} samples\n".\
        format(len(df["sample_id"].unique())))

    if ret_params:    
        return df, params
    else:
        return df
