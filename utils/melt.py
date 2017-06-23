import pandas as pd

def melt(df,dont_melt=None,rowname="rows",colname="cols",valname="value"):
    """
    Melts pandas DataFrame from a matrix-style format to a 3-column format.

    The default behaviour of this program converts this pandas DataFrame:

          col1  col2
    row1  0     1
    row2  2     3
    
    to:
       rows   cols    value
    0  row1   col1    0
    1  row1   col2    1
    2  row2   col1    2
    3  row2   col2    3

    INPUTS
    ------
    df : pandas DataFrame
        The format should be a DataFrame whose columns you want to melt like the
        example provided above.
    rowname : string, default "rows"
        The column name of what used to be the rows of the matrix, i.e. the index.
    colname : string, default "cols"
        The column name of what used to be the rows
    valname : string, default "value"
        The column of the data of the matrix

    RETURNS
    -------
    df2 : pandas DataFrame
        melted dataframe. See above example.
    """
#   for some reason df gets changed globally. This prevents that.
    df2 = df.copy(deep=True)

    df2[rowname] = df2.index
    index_name = [rowname]
    df2 = pd.melt(df2, id_vars=index_name, var_name=colname, value_name=valname)
    return df2
