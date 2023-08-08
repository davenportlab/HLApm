library(seqinr)
library(Biostrings)
library(stringr)
library(DECIPHER)
library(data.table)
library(dplyr)

# need to add a path of 'DIR_personalizedHLA' 

print("Loading IMGT GTF ...")

imgt_gtf_fname <- gzfile( paste0(DIR_personalizedHLA, 
                                 "/data/references/IMGT.3.52.0/two_field.imgt_allele.renamed.gtf.gz"),
                          'rt') 
imgt_gtf <- read.table(imgt_gtf_fname, sep="\t", header=F)
imgt_gtf$V2 <- "IMGT_3.52.0"
close(imgt_gtf_fname)


print("Loading IMGT allele sequences ...")

imgt_fa_fname <- gzfile( paste0(DIR_personalizedHLA,  
                                "/data/references/IMGT.3.52.0/two_field.imgt_allele.renamed.fa.gz"),
                         'rt') 
subject_seq <- suppressWarnings ( read.fasta(file = imgt_fa_fname) )
close(imgt_fa_fname)


print("Loading chr6 ...")
chr6_fname <- gzfile( paste0(DIR_personalizedHLA, 
                             "/data/references/hg38/chr6.fa.gz"),
                          'rt') 
seq_chr6 <- suppressWarnings ( read.fasta(file = chr6_fname) )
seq_chr6 <- toupper( c2s(seq_chr6[[1]]) )
#print ( paste0("length of chr6: ", nchar(seq_chr6) ) )
close(chr6_fname)
print ( paste0("reference data loaded" ) )

print("Loading IMGT exon sequences ...")
exon_seq_fname <- gzfile( paste0(DIR_personalizedHLA, 
                           "/data/references/IMGT.3.52.0/imgt_exon_seq.txt.gz"),
                    'rt') 
exon_seq <- read.table(exon_seq_fname)
exon_seq$id <- gsub("(.*\\*\\d+\\:\\d+).*", "\\1", exon_seq$V1 )
exon_seq$V1 <- gsub("\\*", "\\_", as.character(exon_seq$V1))
exon_seq$V1 <- gsub("\\:", "\\_", exon_seq$V1)
exon_seq$id <- gsub("\\*", "\\_", as.character(exon_seq$id))
exon_seq$id <- gsub("\\:", "\\_", exon_seq$id)

print("Loading alternative exon sequences ...")
exon_alt_name <- gzfile( paste0(DIR_personalizedHLA, 
                           "/data/references/IMGT.3.52.0/imgt_allele_L234.alternative_exons.txt.gz"),
                          'rt') 
exon_alt <- read.table(exon_alt_name)
exon_alt$id <- gsub("(.*\\*\\d+\\:\\d+).*", "\\1", exon_alt$V1 )
exon_alt$id <- gsub("[\\*|\\:]", "\\_", as.character(exon_alt$id))

print("Loading allele levels ...")
allele_level_name <- gzfile( paste0(DIR_personalizedHLA, 
                                "/data/references/IMGT.3.52.0/allele_levels.txt.gz"),
                         'rt') 
allele_level <- read.table(allele_level_name, header = T)


print("Loading HLA regions to mask ...")
hla_to_mask <- read.table( paste0(DIR_personalizedHLA, 
                                  "/data/references/hg38/Ensembl99.primary.12_HLA_genes.bed") )


source(paste0(DIR_personalizedHLA, "/scripts/DRB_haplotypes.R") )
print("Loading DRB genes ...")

DRB_fa <- readDNAStringSet(file = paste0(DIR_personalizedHLA,
                                         "/data/references/hg38/Ensembl99.all_drb_genes.fa.gz"), 
                           format="fasta")

DRB_gtf_fname <- gzfile( paste0(DIR_personalizedHLA, 
                               "/data/references/hg38/Ensembl99.all_drb_genes.gtf.gz"),'rt') 
DRB_gtf <- suppressWarnings ( rtracklayer::import(DRB_gtf_fname) )
close(DRB_gtf_fname)


print("Start: Roading HLA genes and their annotation from primary reference ...")
# pattern_seq : ref_seq from chr6
pattern_seq.all.name <- paste0(DIR_personalizedHLA, 
                          "/data/references/hg38/hla_seq.fa.gz")
pattern_seq.all <- suppressWarnings ( read.fasta(file = pattern_seq.all.name) )

pattern_anno.name <- gzfile(paste0(DIR_personalizedHLA,
                              "/data/references/hg38/hla_anno.txt.gz"),'rt') 
pattern_anno <- read.table(pattern_anno.name, header=T)
close(pattern_anno.name)

print("Finished: Roading IMGT seqeunces, HLA genes and their annotation from primary reference ...")


