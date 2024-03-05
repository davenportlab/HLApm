


# add a path of DIR_personalizedHLA 
# DIR_personalizedHLA=your_PATH/personalisedHLAmapping

# DIR_personalizedHLA="/lustre/scratch126/humgen/teams/davenport/ws/github_ws/HLApm"
source( paste0(DIR_personalizedHLA, "/scripts/load_ref.R") )
# imgt_gtf
# subject_seq : imgt seq
source( paste0(DIR_personalizedHLA, "/scripts/make_personalized_HLA_ref.R") )
source( paste0(DIR_personalizedHLA, "/scripts/align_and_adjust_annotation.R") )




input_data <- paste0(DIR_personalizedHLA, "/data/examples/example_input.txt")
input_alleles <- read.table(input_data, header=T)


# output_directory="/lustre/scratch126/humgen/teams/davenport/ws/github_ws/HLApm/scripts/bulk"

bulkRNA_build_personalized_HLA_ref(input_alleles, output_directory="./" )



