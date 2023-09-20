import pandas as pd
import os
import sys


def getFastq(path):
    """A function to import basecalls """
    imageDf = pd.read_csv(path,  comment="#")
    return imageDf


def getRawFastqPath(fastDf, outdir):
    "A function to get a list of file names for copy"
    for idx, row in fastDf.iterrows():
        pass
        

def expandImageIds(imageDf):
    """A function to return snakemake `expand` inputs """
    iids = []
    tids = []
    for idx, row in imageDf.iterrows():
        
        imgPath = row['filePath']
        imgId = row['fileId']

        with TiffFile(imgPath) as tif:
            series = tif.series[0]
            
        nTimePoints = series.shape[1]
        for t in range(1, nTimePoints+1):
            iids.append(imgId)
            tids.append(t)
            
    return iids, tids


