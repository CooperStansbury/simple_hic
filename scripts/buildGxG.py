import pandas as pd
import os
import numpy as np
import sys
import pyranges as pr


def read_pairs(fpath):
    """A function to read pairs file """
    
    cols = [
        'readID',
        'chrom1',
        'pos1',
        'chrom2',
        'pos2',
        'strand1',
        'strand2',
        'pair_type',
        'sam1',
        'sam2',
        'mapq1',
        'mapq2',
        'read_len1',
        'read_len2',
        'algn_read_span1',
        'algn_read_span2',
        'algn_ref_span1',
        'algn_ref_span2'
    ]
    
    usecols = [
        'readID',
        'chrom1',
        'pos1',
        'chrom2',
        'pos2',
        'strand1',
        'strand2',
        'pair_type',
        'mapq1',
        'mapq2',
        'algn_ref_span1',
        'algn_ref_span2'
    ]
    
    dtypes = {
        'readID': 'str',
        'chrom1': 'str',
        'pos1': 'int',
        'chrom2': 'str',
        'pos2': 'int',
        'strand1': 'str',
        'strand2': 'str',
        'pair_type': 'str',
        'mapq1': 'str',
        'mapq2': 'str',
        'algn_ref_span1': 'str',
        'algn_ref_span2': 'str',
    }
    
    df = pd.read_csv(fpath, 
                     sep="\t", 
                     names=cols,
                     usecols=usecols,
                     dtype=dtypes,
                     na_filter=False,
                     comment="#", 
                     compression='gzip',
                     low_memory=True,
                     header=None)

    return df
    

def makePyRanges(df, frag=1):
    """A function to produce a PyRanges frame from a pairs 
    dataframe """

    cols = [
        'readID', 
        f'chrom{frag}', 
        f'pos{frag}',
        f'strand{frag}',
        f'algn_ref_span{frag}',
        f'mapq{frag}'
    ]

    prdf = df[cols].copy()
    prdf.columns = ['readID', 'Chromosome', 'Start', 'Strand', 'Span', 'MAPQ']
    prdf['End'] = prdf['Start'] + prdf['Span']
    prdf = prdf[prdf['End'].notna()]
    
    prdf['Start'] = prdf['Start'].astype(int)
    prdf['End'] = prdf['End'].astype(int)

    prdf = prdf[['readID', 'Chromosome', 'Start', 'End', 'Strand', 'MAPQ']]
    prdf['Fragment'] = frag
    return pr.PyRanges(prdf)


def getGeneHits(prdf, gr, slack=300):
    """A function to get contact hits. Returns a list of readID """
    prdf = prdf.join(gr, slack=slack)
    prdf = prdf.as_df() # convert back to df
    prdf = prdf[prdf['gene_name'] != "-1"] # drop contacts not in gene regions
    return prdf

    
if __name__ == "__main__":
    pairs_path = sys.argv[1]  
    genes_path = sys.argv[2]  
    slack = int(sys.argv[3])
    outpath = sys.argv[4]  

    print(f"Working {pairs_path=}")

    chromlist = [
        '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', 
        '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y',
    ]

    # read pairs
    df = read_pairs(pairs_path)
    print(f"(pre-filter) {df.shape=}")

    # simple filters
    df = df[df['pair_type'] == 'UU'].reset_index(drop=True)
    df = df[df['chrom1'].isin(chromlist)].reset_index(drop=True)
    df = df[df['chrom2'].isin(chromlist)].reset_index(drop=True)
    
    df = df[df['mapq1'] != ''].reset_index(drop=True)
    df = df[df['mapq2'] != ''].reset_index(drop=True)

    df = df[df['algn_ref_span1'] != ''].reset_index(drop=True)
    df = df[df['algn_ref_span2'] != ''].reset_index(drop=True)

    # data typing
    df['mapq1'] = df['mapq1'].astype(int)
    df['mapq2'] = df['mapq2'].astype(int)
    
    df['algn_ref_span1'] = df['algn_ref_span1'].astype(int)
    df['algn_ref_span2'] = df['algn_ref_span2'].astype(int)
    
    print(f"(post-filter) {df.shape=}")

    print(df.info(memory_usage="deep"))

    # read gene positions
    genes = pd.read_csv(genes_path)
    gr = pr.PyRanges(genes)

    # build pyranges 
    p1 = makePyRanges(df, frag=1)
    p2 = makePyRanges(df, frag=2)

    print(f"Build ranges done.")

    del df

    # get gene hits
    p1 = getGeneHits(p1, gr, slack=slack)
    p2 = getGeneHits(p2, gr, slack=slack)

    print(f"Gene hits done.")

    del gr

    # filter for common
    shared = np.intersect1d(p1['readID'], p2['readID'])
    p1 = p1[p1['readID'].isin(shared)]
    p2 = p2[p2['readID'].isin(shared)]

    print(f"Read filtering done.")

    cols = [
        'readID', 
        'Chromosome', 
        'Start_b', 
        'End_b', 
        'Fragment',
        'gene_name',
        'MAPQ',
    ]

    p1 = p1[cols]
    p2 = p2[cols]

    # merge the gene x gene relationsships
    gx = pd.merge(p1, p2, 
              how='left',
              left_on='readID',
              right_on='readID',
              suffixes=['_1', '_2'])

    print(f"Merge done.")

    gx.to_csv(outpath, index=False)

    

    

    
    



    