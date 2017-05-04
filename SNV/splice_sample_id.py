import re
import pandas as pd

error1 = \
"""Your separator does not split each row in sample_id
into 2 strings. In other words, sep is not present exactly
once in every sample id"""

def splice_sample_id(df,sep=None):
    """
    Splices sample_id in the SNV table to easily filter samples and cohorts

    This use of this function splices the sample_id column into 2 columns in
    order to easily filter or group the table according to cohort or
    environmental conditions. Obviously, this method relies on the fact that
    cohort or environmental information is embedded in your sample names.
    Ideally, you named your samples to contain some tag of your cohort,
    followed by some sample number or environmental condition. For example,
    let's suppose you collect 3 samples each for 2 patients. If your samples
    were named ABC1, ABC2, ABC3, DEF1, DEF2, and DEF3, then this method creates
    2 columns `cohort` and `envrionment` that look like:

    sample_id     cohort     environment
    ABC1          ABC        1
    ABC2          ABC        2
    ABC3          ABC        3
    DEF1          DEF        1
    DEF2          DEF        2
    DEF3          DEF        3

    This illustrates the default behaviour of this method. It splices sample_id
    at the last non-digit character. If your sample_id is more complicated,
    `sep` must be provided (see Attributes). For example, suppose your sample
    names are PATIENT1_GUT, PATIENT1_ORAL, PATIENT2_GUT, and PATIENT2_ORAL.
    Then with sep="_", the method creates 2 columns `cohort` and `environment`
    that look like:

    sample_id      cohort     environment
    PATIENT1_GUT   PATIENT1   GUT
    PATIENT1_ORAL  PATIENT1   ORAL
    PATIENT2_GUT   PATIENT2   GUT
    PATIENT2_ORAL  PATIENT2   ORAL

    NOTE: If the resultant column is ALL digit characters, the column is converted
    from a string column to an int column

    Attributes
    ----------
    df : pandas Dataframe
        SNV table that the cohort and environment columns will be appended to.
    sep : str 
        a string identifier present in sample_id that separates the
        cohort info from environment info that should be present exactly once in
        each sample id.

    Returns
    -------
    df : pandas DataFrame
        The SNV table with cohort and environment columns added.
    """

    def find_split_index(string):
        matches = [match for match in re.finditer(r"\D\d",string)]
        split_index = matches[-1].start() + 1
        return string[:split_index], string[split_index:]

    col_list = ["cohort","environment"]

#   default behaviour: split at last non-digit number
    if sep == None:
        for ind,col in enumerate(col_list):
            df[col] = df["sample_id"].apply(lambda x: find_split_index(x)[ind])

#   behaviour if sep is provided
    else:
        try:
            df[col_list] = df["sample_id"].str.split(sep,expand=True)
        except:
            raise ValueError(error1)
        
#   if either columns are all digits, they are converted from strings to ints
    for col in col_list:
        isdigit_list = [x.isdigit() for x in df[col]]
        if all(isdigit_list):
            df[col].astype(int)
            print("WARNING: '{}' converted to type int".format(col))
#   returns SNV table
    return df

