import pandas as pd
import yaml
from pathlib import Path
import re
import os
import sys
from utils import utils

BASE_DIR = Path(workflow.basedir)
configfile: str(BASE_DIR) + "/config/config.yaml"

# big picture variables
OUTPUT = config['output_path']

# structure file names and get Id lists
fastqPath = os.path.abspath(config['fastq_path'])
fastqDf = utils.getFastq(fastqPath)
fids = fastqDf['fileId'].to_list()

print(fastqDf)
