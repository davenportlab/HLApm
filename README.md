# PersonalisedHLAmapping
Personalised HLA mapping to quantify HLA gene expression from scRNA and bulkRNA sequencing data



* [Introduction](#Introduction)
* [Requirements](#Requirements)
* [Usage](#usage)
* [Contact](#contact)



   
   
# Introduction


- HLA alleles from 20 HLA genes are available. 
"HLA-A"    "HLA-B"    "HLA-C"    "HLA-DMA"  "HLA-DMB"  "HLA-DOA"  "HLA-DOB"  "HLA-DPA1" "HLA-DPB1" "HLA-DQA1" "HLA-DQA2" "HLA-DQB1", "HLA-DRA"  "HLA-DRB1" "HLA-DRB3" "HLA-DRB4" "HLA-DRB5" "HLA-E"    "HLA-F"    "HLA-G"   



# Requirements


* R packages

```R
# packages to install to build personalized referece :

install.packages( "data.table", "dplyr", "stringr")

source("https://bioconductor.org/biocLite.R")
biocLite( c("Biostrings", "rtracklayer", "DECIPHER") )

```

* bedtools 

```bash
# To install conda package  
conda install -c bioconda bedtools

# or see here, https://bedtools.readthedocs.io/en/latest/content/installation.html 
```

* Tools for aligning reads with scRNA-seq data

cellranger-7.0.0 or later




* Primary reference and annotation: hg38 human genome 

```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/GRCh38.primary_assembly.genome.fa.gz
```
* hg38 annotation version 33 (Ensembl 99)

```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.primary_assembly.annotation.gtf.gz
```





# Usage



## 1. Prepare input data: HLA allelels 


1. HLA alleles 
* see an [example input](./data/examples/example_input.txt)  
- tab-delimited text file
- First column - individual ids : will be used as a prefix of output file name
- Second column - HLA alleles : HLA prefix with at least 2-field resolution (e.g., HLA-A*02:05:01, HLA-DQA1*03:01)

2. Raw scRNA-seq data 
- FASTQ OR BAM



## 2. Build personalized reference and annotation 


HLA alleles from the [IMGT/HLA database](https://github.com/ANHIG/IMGTHLA) were preprocessed to fill missing exons and UTR regions.

Some alleles in the IMGT/HLA database have (untranslated) exons missing. UTR annotation is absent or variable in different alleles. 
For instance, the HLA-B gene has 8 exons, which has a stop codon in the 7th exon and a 3’ UTR overlapping with the 8th exon. Reported HLA-B allele sequences in the IMGT don’t have 8th exon (UTR) annotation and some of them have even no sequences of the UTR (sequences at the 7th exon.). Without UTR sequences or UTR annotation, quantification of gene expression can be underestimated. 
To complete gene sequences and annotation, we extend both 5’ and 3’ sequences of each allele in the IMGT based on primary reference sequences. And we adjusted gene coordinates based on extended allele sequences



* [here](./data/examples/example_run.R)'s an example script to build personalized HLA reference.



### 2.1. Set a path of a package

A path name should be `DIR_personalizedHLA`

```R
# where the package installed
yourpath="./"
DIR_personalizedHLA=paste0(yourpath, "/personalisedHLAmapping") 
```



### 2.2. Gernerate personalized HLA reference and annotation for each individual

A single `build_personalized_HLA_ref` 

```R
source( paste0(DIR_personalizedHLA, "/scripts/load_ref.R") )
source( paste0(DIR_personalizedHLA, "/scripts/align_and_adjust_annotation.R") )
source( paste0(DIR_personalizedHLA, "/scripts/build_personalized_HLA_ref.R") )

# reading input alleles
input_alleles <- read.table(paste0(DIR_personalizedHLA, 
                                  "/data/examples/example_input.txt"), head=T)

# Get personalized HLA reference and annotation
build_personalized_HLA_ref(input_alleles, output_directory="./your_directory" )

```

This will produce 

- two files (per.fa and per.gtf) for each individual and 
- one bed file to mask HLA genes in primary reference and annotation;


1) fasta file of personalized HLA genes: 
__`output_directory/out/individual_ID/individual_ID.per.fa`__

2) gtf file of personalized HLA genes:
__`output_directory/out/individual_ID/individual_ID.per.gtf`__

3) bed file with HLA genes to mask primary fasta and gtf annotation: 
__`output_directory/out/mask_these_genes.bed`__






### 2.3. Mask primary reference and annotation. 


- HLA genes __only__ used for personalized mapping. 
- Use a `mask_these_genes.bed` file generated by `build_personalized_HLA_ref`. 
- <span style="color:blue">__*You can skip this step if you already have masked reference and annotation for the same set of HLA genes.*__ </span>




```bash
# primary fasta 
## use `bedtools maskfasta`

bedtools maskfasta -fi GRCh38.primary_assembly.genome.fa -bed mask_these_genes.bed -fo GRCh38.primary.HLA_masked.fa
```


```R
# primary gtf

gtf.primary <- rtracklayer::import(gzfile("gencode.v33.primary_assembly.annotation.gtf.gz"))
hla_genes_to_mask <- read.table("output_directory/mask_these_genes.bed")

HLA_genes.masked <- gtf.primary[!(gtf.primary$gene_name %in% hla_genes_to_mask$V4), ]
rtracklayer::export(HLA_genes, "GRCh38.primary.HLA_masked.gtf"))
```


### 2.4. Get final personalized genome and annotation for each individual

* hg38 masked reference + personalized HLA genes for each individual

- Using R

```R
# reference sequences

refSeq.primary.masked <- GRCh38.primary.HLA_masked.fa
refSeq.personalizedHLA <- output_directory/out/individual_ID/individual_ID.per.fa

Biostrings::writeXStringSet(refSeq.final, 
                            "output_directory/out/individual_ID/individual_ID.primaryMasked_and_HLA.fa" )

# gtf

refGTF.primary.masked <- rtracklayer::import("GRCh38.primary.HLA_masked.gtf")
refGTF.personalizedHLA <- rtracklayer::import("output_directory/out/individual_ID/individual_ID.per.gtf")
refGTFfinal <- c(refGTF.personalizedHLA, refSeq.primary.masked)
rtracklayer::export(refGTF.final, 
                    "output_directory/out/individual_ID/individual_ID.primaryMasked_and_HLA.gtf"))

```




# Mapping


## scRNA-seq mapping with Cell Ranger


### Indexing: mkref
```bash
#!/bin/bash

cellranger \
        mkref \
        --nthreads=12 \
        --genome=output_directory/ref.individual_ID \
        --fasta=output_directory/out/individual_ID/individual_ID.primaryMasked_and_HLA.fa \
        --genes=output_directory/out/individual_ID/individual_ID.primaryMasked_and_HLA.gtf

```

### Mapping: count


[here](./scripts/cellranger.count.sh)'s a script to run `cellranger count`.

```bash
./scripts/cellranger.count.sh sample_name fastq_path path_of_personalised_ref out_prefix

# ./scripts/cellranger.count.sh individual_ID fastq_path output_directory/ref.individual_ID out_prefix
```




## bulkRAN-seq







### Contact

You are welcome to:
- submit suggestions and bug-reports at: https://github.com/davenportlab/personalisedHLAmapping/issues
- compose a friendly e-mail to: wl2@sanger.ac.uk


## to do..
make directory tutorials, data (parsed imgt fa and gtf), example,


