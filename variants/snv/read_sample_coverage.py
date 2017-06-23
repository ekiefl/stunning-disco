import pandas as pd

def read_sample_coverage(mean_cov_file,thresh=None,save=False):
    """
    Read in {collection}-mean_coverage(_Q2Q3).txt as python DataFrame.
    
    Returns {collection}-mean_coverage(_Q2Q3).txt file generated from
    anvi-summarize as a pandas DataFrame that has been filtered to only include
    samples that have a mean coverage above a threshold value. Can also save a
    new file that's the same format as {collection}-mean_coverage(_Q2Q3).txt
    except with the samples below the threshold filtered out.

    INPUTS
    ------
    mean_cov_file : string
        The file path to {collection}-mean_coverage.txt or 
        {collection}-mean_coverage_Q2Q3.txt 
    thresh : integer, default None
        if not None, the df is filtered to only include samples
        with >= thresh
    save : Boolean, default False
        saves output in same format as {collection}-mean_coverage.txt
        except with those with less than thresh filtered out. It is
        output to the same filepath as {collection}-mean_coverage.txt
        with the name {collection}-mean_coverage_gt{thresh}x.txt 

    RETURNS
    -------
    df : pandas DataFrame
        df with columns = ("bin", "mean_coverage(_Q2Q3)") and as many
        rows as there are samples.
    """

    df = pd.read_csv(mean_cov_file, header=None, index_col=0, sep='\t')
    df = df.transpose()
    cols = df.columns.values
    df[cols[1]] = df[cols[1]].astype(float)
#   sorts samples by coverage just because
    df = df.sort_values(by=cols[1],ascending=True).reset_index(drop=True)

#   filter samples by coverage values
    if thresh != None:
        df = df[df[cols[1]] >= thresh]

#   write new file if save = True
    if save:
        out = df.transpose()
        name = mean_cov_file.replace(".txt","")+"_gt{}x".format(thresh)+".txt"
        out.to_csv(name,index=True,sep="\t",header=False)

    return df
