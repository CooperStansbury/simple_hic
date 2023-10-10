## Adapted from: https://data.4dnucleome.org/resources/data-analysis/hi_c-processing-pipeline

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
rids = fastqDf['rep'].to_list()


rule all:
    input:
        OUTPUT + 'references/reference.fa',
        OUTPUT + 'references/reference.bwt',
        OUTPUT + 'references/reference.chrom.sizes',
        expand(f"{OUTPUT}fastq/{{fids}}.{{rids}}.raw.fastq.gz", zip, fids=fids, rids=rids),
        expand(f"{OUTPUT}alignment/{{fids}}.bam", fids=fids),
        expand(f"{OUTPUT}pairs/{{fids}}.pairs.gz", fids=fids),
        expand(f"{OUTPUT}pairs/{{fids}}.sorted.pairs.gz", fids=fids),
        expand(f"{OUTPUT}pairs/{{fids}}.final.pairs.gz", fids=fids),
        expand(f"{OUTPUT}cooler/{{fids}}.base.cool", fids=fids),
        expand(f"{OUTPUT}cooler/{{fids}}.mcool", fids=fids),

rule getReference:
    input:
        refgenome=config['ref_path'],
    output:
        OUTPUT + 'references/reference.fa.gz'
    shell:
        "cp {input} {output}"
        
        
rule prepReference:
    input:
        refgenome=OUTPUT + 'references/reference.fa.gz'
    output:
        ref=OUTPUT + 'references/reference.fa',
        flag=touch(OUTPUT + 'references/reference.done')
    shell:
        "cat {input} | gzip -d > {output.ref}"


rule bwa_index:
    input:
        OUTPUT + 'references/reference.fa',
    output:
        idx=multiext(OUTPUT + 'references/reference', 
                         ".amb", 
                         ".ann", 
                         ".bwt", 
                         ".pac", 
                         ".sa"),
    params:
        algorithm="bwtsw",
    wrapper:
        "v2.6.0/bio/bwa/index"


rule make_chromsizes:
    input:
        OUTPUT + 'references/reference.fa',
    output:
        OUTPUT + 'references/reference.chrom.sizes'
    shell:
        "faidx {input} -i chromsizes > {output}"


rule getFastq:
    input:
        fastq=fastqDf['filePath'].to_list(),
    output:
        expand(f"{OUTPUT}fastq/{{fids}}.{{rids}}.raw.fastq.gz", zip, fids=fids, rids=rids)
    wildcard_constraints:
        rids='|'.join([re.escape(x) for x in set(rids)]),
        fids='|'.join([re.escape(x) for x in set(fids)]),
    run:
        from shutil import copyfile
        for i, refPath in enumerate(input.fastq):

            outPath = output[i]
            copyfile(refPath, outPath)


rule bwa_mem:
    input:
        reads=[OUTPUT + "fastq/{fid}.r1.raw.fastq.gz", OUTPUT + "fastq/{fid}.r2.raw.fastq.gz"],
        idx=multiext(OUTPUT + 'references/reference', 
                         ".amb", 
                         ".ann", 
                         ".bwt", 
                         ".pac", 
                         ".sa"),
    output:
        OUTPUT + "alignment/{fid}.bam"
    params:
        extra=r"-5SPM -R '@RG\tID:{fid}\tSM:{fid}'",
        sorting="none", 
        sort_order="queryname", 
        sort_extra="",
    threads: 
        16
    wildcard_constraints:
        fids='|'.join([re.escape(x) for x in set(fids)]),
    wrapper:
        "v2.6.0/bio/bwa/mem"


rule samtools_index:
    input:
        OUTPUT + "alignment/{fid}.bam"
    output:
        OUTPUT + "alignment/{fid}.bam.bai"
    params:
        extra="",  # optional params string
    threads: 
        4 
    wrapper:
        "v2.6.0/bio/samtools/index"


rule pairtools_parse:
    input:
        sizes=OUTPUT + 'references/reference.chrom.sizes',
        bam=OUTPUT + "alignment/{fid}.bam",
        # bai=OUTPUT + "alignment/{fid}.bam.bai"
    output:
        pairs=OUTPUT + "pairs/{fid}.pairs.gz",
        stats=OUTPUT + "pairs/{fid}.pairs.stats",
    wildcard_constraints:
        fids='|'.join([re.escape(x) for x in set(fids)]),
    shell:
        """pairtools parse -o {output.pairs} -c {input.sizes} \
          --output-stats {output.stats} \
          --assembly hg38 \
          --add-columns mapq,read_len,algn_read_span,algn_ref_span \
          {input.bam}"""


rule pairtools_sort:
    input:
        pairs=OUTPUT + "pairs/{fid}.pairs.gz",
    output:
        sorted=OUTPUT + "pairs/{fid}.sorted.pairs.gz",
    wildcard_constraints:
        fids='|'.join([re.escape(x) for x in set(fids)]),
    shell:
        """pairtools sort --output {output.sorted} {input.pairs} """ 


rule pairtools_dedup:
    input:
        pairs=OUTPUT + "pairs/{fid}.sorted.pairs.gz",
    output:
        uniqpairs=OUTPUT + "pairs/{fid}.nodups.pairs.gz", 
        unmappedpairs=OUTPUT + "pairs/{fid}.unmapped.pairs.gz",
        dupspairs=OUTPUT + "pairs/{fid}.dups.pairs.gz",
        dupstats=OUTPUT + "pairs/{fid}.dedup.stats",
    wildcard_constraints:
            fids='|'.join([re.escape(x) for x in set(fids)]),
    shell:
        """pairtools dedup --max-mismatch 3 --mark-dups \
               --output {output.uniqpairs} \
               --output-unmapped {output.unmappedpairs} \
               --output-dups {output.dupspairs} \
               --output-stats {output.dupstats} \
                {input.pairs} """


rule pairtools_filter:
    input:
        pairs=OUTPUT + "pairs/{fid}.nodups.pairs.gz",
    output:
        OUTPUT + "pairs/{fid}.final.pairs.gz",
    wildcard_constraints:
                fids='|'.join([re.escape(x) for x in set(fids)]),
    shell:
        """pairtools select 'mapq1>0 and mapq2>0' {input.pairs} -o {output}"""


rule to_cooler:
    input:
        pairs=OUTPUT + "pairs/{fid}.final.pairs.gz",
        sizes=OUTPUT + 'references/reference.chrom.sizes',
    output:
        cooler=OUTPUT + "cooler/{fid}.base.cool"
    params:
        bins=config['bins'] 
    wildcard_constraints:
            fids='|'.join([re.escape(x) for x in set(fids)]),
    shell:
        """cooler cload pairs \
            -c1 2 -p1 3 -c2 4 -p2 5 \
            --assembly hg38 \
            {input.sizes}:{params.bins} \
            {input.pairs} \
            {output}"""


rule cooler_zoom:
    input:
        cooler=OUTPUT + "cooler/{fid}.base.cool"
    output:
        OUTPUT + "cooler/{fid}.mcool"
    threads: 
        8
    params:
        zoom=",".join(map(str, config["zoomify"]))
    shell:
        """cooler zoomify \
                --nproc {threads} \
                --out {output} \
                --resolutions {params.zoom} \
                --balance \
                {input.cooler} """