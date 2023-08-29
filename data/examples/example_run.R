

# add a path of DIR_personalizedHLA 
# DIR_personalizedHLA=your_PATH/personalisedHLAmapping

source( paste0(DIR_personalizedHLA, "/scripts/load_ref.R") )
source( paste0(DIR_personalizedHLA, "/scripts/align_and_adjust_annotation.R") )
source( paste0(DIR_personalizedHLA, "/scripts/make_personalized_HLA_ref.R") )


input_data <- paste0(DIR_personalizedHLA, "/data/examples/example_input3.txt")
input_alleles <- read.table(input_data, header=T)


build_personalized_HLA_ref(input_alleles, output_directory="./data/examples" )
# build_personalized_HLA_ref(input_alleles, output_directory="/lustre/scratch126/humgen/teams/davenport/ws/ref/test" )


# 'make_personalized_HLA_ref' will produce 3 files
# 1. fasta file of personalized HLA genes for each individual
# 2. gtf file of personalized HLA genes for each individual
# 3. bed file with HLA genes to mask primary fasta and gtf annotation





