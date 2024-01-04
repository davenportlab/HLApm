


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




print("### Build an index ...")






OUT_DIR= paste0(main_dir, "idx/", allele_name)
fasta_file = paste0(output_imgt, allele_name, ".fa")
gtf_file = paste0(output_imgt, allele_name, ".gtf")


if ( !file.exists(OUT_DIR)){
  dir.create(file.path(OUT_DIR ))
}

cmd=paste0("/software/team282/download/STAR-2.7.3a/bin/Linux_x86_64/STAR --runThreadN 4 --runMode genomeGenerate ",
           "--genomeSAindexNbases 8 ",
           "--genomeDir ", OUT_DIR, " ",
           "--genomeFastaFiles ", fasta_file, " ",
           "--sjdbGTFfile ", gtf_file, " ",
           "--sjdbOverhang 99 --outTmpDir ", OUT_DIR, "/tmp")

print(cmd)
