#!/bin/bash


echo "### Build an index ..."


# $1, runThreadN : ncpus
# $2, genomeDir : OUT_DIR
# $3, genomeFastaFiles 
# $4, sjdbGTFfile

STAR="/software/team282/download/STAR-2.7.3a/bin/Linux_x86_64/STAR"



mkdir $2


$STAR --runThreadN $1 \
	--runMode genomeGenerate \
	--genomeSAindexNbases 8 \
	--genomeDir  $2 \
	--genomeFastaFiles $3 \
	--sjdbGTFfile $4 \
	--sjdbOverhang 99 \
	--outTmpDir $2/tmp

