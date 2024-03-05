#!/bin/bash

echo "### Start mapping ..."

STAR="/software/team282/download/STAR-2.7.3a/bin/Linux_x86_64/STAR"


# $1, runThreadN : ncpus
# $2, genomeDir : directory of genome index
# $3, fastq_1
# $4, fastq_2
# $5, outFileNamePrefix


       # STAR="/software/team282/download/STAR-2.7.3a/bin/Linux_x86_64/STAR"
        # runThreadN=4
        # index_dir="/lustre/scratch119/realdata/mdt3/teams/davenport/ws/ref/IMGT_star_index/idx/chr_HLA-A_03_01_01_01"
        # fastq_1="/lustre/scratch119/humgen/projects/gains_team282/hla/1.1.fq_for_mhc/gains8032857/gains8032857.forMHCmapping_1.fq.gz"
        # fastq_2="/lustre/scratch119/humgen/projects/gains_team282/hla/1.1.fq_for_mhc/gains8032857/gains8032857.forMHCmapping_2.fq.gz"
        # prefix="./"



$STAR --runMode alignReads --runThreadN $1 \
	--genomeDir $2 \
	--readFilesIn $3 $4 \
	--readFilesCommand zcat \
	--outFileNamePrefix $5. \
	--outSAMtype BAM Unsorted \
	--quantMode GeneCounts


picard SortSam I=$5.Aligned.out.bam \
        O=$5.queryname.sorted.bam \
        SORT_ORDER=queryname

samtools index $5.queryname.sorted.bam

#samtools sort $5.Aligned.out.bam -o $5.Aligned.sorted.bam
#samtools index $5.Aligned.sorted.bam








