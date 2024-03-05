#!/bin/bash

echo "### Sorting bam by queryname ..."

# conda activate picard 


# $1, The SAM, BAM or CRAM file to sort
# $2, The sorted SAM, BAM or CRAM output file

picard SortSam I=$1 \
	O=$2 \
	SORT_ORDER=queryname







