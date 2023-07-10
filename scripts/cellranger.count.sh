#!/bin/bash

### input 
# $1. --id <ID> : A unique run id and output folder name [a-zA-Z0-9_-]+
# $2. --fastqs <PATH>: Path to input FASTQ data
# $3. --transcriptome <PATH>: Path of folder containing 10x-compatible transcriptome reference
# $4. output_prefix


mkdir -p $4
cd $4

cellranger \
        count --id=$1 \
        --fastqs=$2 \
        --sample=$1 \
        --transcriptome=$3 \
        --expect-cells=10000 \
        --localcores=35 \
        --localmem=200 \
        --chemistry=auto
cd -

