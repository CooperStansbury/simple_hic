import sys
import pandas as pd
import pyranges as pr
import numpy as np
from pathlib import Path


if __name__ == "__main__":
    gtf_path = sys.argv[1]  
    outpath = sys.argv[2]  

    # read the GTF file
    gr = pr.read_gtf(gtf_path)
    gf = gr.as_df()

    # get only protien-coding genes
    genes = gf[gf['Feature'] == 'gene'].reset_index(drop=True)
    genes = genes[genes['gene_biotype'] == 'protein_coding']

    # subset the columns
    gene_cols = [
        'Chromosome',
        'Start', 
        'End', 
        'Strand',
        'gene_name'
    ]
    
    genes = genes[gene_cols]

    # write output
    genes.to_csv(outpath, index=False)