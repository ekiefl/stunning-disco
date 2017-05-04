
import pandas as pd
import numpy as np
import SNV as snv

def load_snv_table(fname):
    """
    Loads snv_table generated by the output of anvio-gen-variability-profile
    
    INPUTS
    
    fname:
        file name.
        
    RETURNS
    
    df:
        pandas DataFrame obect
    """
    
    df = pd.read_csv(fname,sep='\t',header=0,index_col=False)
    print(("\ndf has {} entries\n".format(np.shape(df)[0])))
    return df