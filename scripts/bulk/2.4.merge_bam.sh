#!/bin/bash

echo "### Start mapping ..."

STAR="/software/team282/download/STAR-2.7.3a/bin/Linux_x86_64/STAR"

# $1, runThreadN : ncpus
# $2, name of output bam file
# $3, List of input BAM filenames, one per line


       # STAR="/software/team282/download/STAR-2.7.3a/bin/Linux_x86_64/STAR"
        # runThreadN=4
        # index_dir="/lustre/scratch119/realdata/mdt3/teams/davenport/ws/ref/IMGT_star_index/idx/chr_HLA-A_03_01_01_01"
        # fastq_1="/lustre/scratch119/humgen/projects/gains_team282/hla/1.1.fq_for_mhc/gains8032857/gains8032857.forMHCmapping_1.fq.gz"
        # fastq_2="/lustre/scratch119/humgen/projects/gains_team282/hla/1.1.fq_for_mhc/gains8032857/gains8032857.forMHCmapping_2.fq.gz"
        # prefix="./"



# samtools merge -@ 3 merged.bam A_01_01.bam A_01_01.bam A_01_02.bam B_01_01.bam B_01_02.bam C_01_01.bam C_01_02.bam  ..
samtools merge -@ $1 $2 -b $3

#samtools sort $5.Aligned.out.bam -o $5.Aligned.sorted.bam
#samtools index $5.Aligned.sorted.bam




