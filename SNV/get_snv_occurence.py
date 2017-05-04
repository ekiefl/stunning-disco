import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def get_snv_occurence(snv_table,
					  plot=False):
    """
	For every SNV position, the number of samples that have non-zero 
	variability at that position is calculated.
	
	INPUTS
	
	snv_table:
		SNV table DataFrame object
	plot:
		If true, SNV occurences are plotted from high to low occurence
		
	RETURNS
	
	(snv_pos_freq, freq):
		freq and snv_pos_freq are numpy arrays, each the length of the number
		of SNV positions in the reference genome. snv_pos_freq is the
		position in the contig of each SNV position and freq is the number
		of samples that variation is observed for that SNV position, ordered
		from most observed to least observed.
    """
    
    snv_pos = snv_table['pos_in_contig']

    snv_pos_unique = snv_pos.unique()
    freq = np.zeros(len(snv_pos_unique))
    
    for ind,pos in enumerate(snv_pos_unique):
        freq[ind] = len(snv_pos[snv_pos==pos])

    order = np.argsort(freq)[::-1]

    snv_pos_freq = snv_pos_unique[order]
    freq = freq[order]
    
    if plot:
        plt.scatter(np.arange(len(freq)),freq)
        plt.show()
        plt.close()
        
    return (snv_pos_freq,freq)