import pandas as pd

def blast_to_pandas(f, cols=None):
    """
    reads the results of a BLAST search as a pandas DataFrame. Unless otherwise
    specified, this function assumes the flag `--outfmt 6` was used. Just as an
    example, if instead the flag `outfmt "6 qseqid sseqid pident qlen length
    mismatch gapopen qstart qend sstart send evalue bitscore"` was used, then a
    second argument `cols` should be passed, defined below.

    PARAMETERS
    ----------
    f : string
        The table output by BLAST.
    cols : string, default None
        The manually supplied outfmt format given to blastp. For example, if you used
        `--outfmt 6 "qseqid sseqid qlen"`, cols = "qseqid sseqid qlen"
    
    RETURNS
    -------
    df : pandas DataFrame
        The blast search as a pandas DataFrame with correct column headers
        and no columns used as indices.
    """

    if cols:
        cols = cols.split(" ") 
    else:
        cols = ["qseqid","sseqid","pident","length","mismatch","gapopen",
                "qstart","qend","sstart","send","evalue","bitscore"]

    return pd.read_csv(f, sep="\t", names=cols)

