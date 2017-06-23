#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from functools import reduce

def concat_cols(df,sep,*cols):
    """
    Returns a pandas Series that is the string concatenation of several 
    columns in a pandas DataFrame. 
    
    INPUTS
    
    df:
        the DataFrame from which you are concatenating columns
    sep:
        value separator in the concatenated column. e.g. if sep='/', a row in 
        the concatenated column would look like "str1/str2/str3..."
    *cols:
        all of the column names to be concatenated. To make it clear, when the
        function is passed, *cols should be passed as a series of column names
        e.g. gen_input_for_cluster(df,sep,"column1","column2","column3")
        
    Example usage:
        x['new_col'] = concat_cols(x,'-','col1','col2','col3')
    """
    from functools import reduce
    return reduce(lambda x, y: x.astype(str).str.cat(y.astype(str), sep=sep),
                  [df[col] for col in cols])
