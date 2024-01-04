
# Usage for scRNA-seq data

<br>




## 1. Prepare input data: HLA allelels 

<br>


### 1.1 If you already have HLA alleles for each individual 
<br>
Please make sure if you have a right input format

* see an [example input](../data/examples/example_input.txt)  
- tab-delimited text file
- First column - individual ids : will be used as a prefix of output file name
- Second column - HLA alleles : HLA prefix with at least 2-field resolution (e.g., HLA-A\*02:05:01, HLA-DQA1\*03:01)



### 1.2 If you need to call HLA alleles

Please see [here](./HLA_typing.md)



<br>

## 2. Build personalized reference and annotation 

<br>
HLA alleles from the [IMGT/HLA database](https://github.com/ANHIG/IMGTHLA) were preprocessed to fill missing exons and UTR regions.

Some alleles in the IMGT/HLA database have (untranslated) exons missing. UTR annotation is absent or variable in different alleles. 
For instance, the HLA-B gene has 8 exons, which has a stop codon in the 7th exon and a 3’ UTR overlapping with the 8th exon. Reported HLA-B allele sequences in the IMGT don’t have 8th exon (UTR) annotation and some of them have even no sequences of the UTR (sequences at the 7th exon.). Without UTR sequences or UTR annotation, quantification of gene expression can be underestimated. 
To complete gene sequences and annotation, we extend both 5’ and 3’ sequences of each allele in the IMGT based on primary reference sequences. And we adjusted gene coordinates based on extended allele sequences


<br>

### 2.1 simple run example

*  [here](../data/examples/example_run.R)'s an example script to build personalized HLA reference.

<br>

### 2.2 Step by step procedures
<br>

#### A. Set a path of a package

* *Set a path of this package with a name __`DIR_personalizedHLA`__*
* a directory of `scripts` should be under the `DIR_personalizedHLA`

```R
# where the package installed
yourpath="./"
DIR_personalizedHLA=paste0(yourpath, "/personalisedHLAmapping") 
```



#### B. Gernerate personalized HLA reference and annotation for each individual


Just running a single script `build_personalized_HLA_ref` can provide personalized HLA reference. 

```R
source( paste0(DIR_personalizedHLA, "/scripts/load_ref.R") )
source( paste0(DIR_personalizedHLA, "/scripts/align_and_adjust_annotation.R") )
source( paste0(DIR_personalizedHLA, "/scripts/build_personalized_HLA_ref.R") )

```

All resource for building references loaded ! 

Next, read input alleles. 

```R
# reading input alleles
input_alleles <- read.table(paste0(DIR_personalizedHLA, 
                                  "/data/examples/example_input.txt"), head=T)

# Get personalized HLA reference and annotation
build_personalized_HLA_ref(input_alleles, output_directory="./your_directory" )

```

#### C. ouput 

You will have three files:

- For each individual: ```per.fa``` and ```per.gtf``` 
- one ```bed file``` to mask HLA genes in primary reference and annotation;


1) fasta file of personalized HLA genes: 
`output_directory/out/individual_ID/individual_ID.per.fa`

2) gtf file of personalized HLA genes:
`output_directory/out/individual_ID/individual_ID.per.gtf`

3) bed file with HLA genes to mask primary fasta and gtf annotation: 
`output_directory/out/mask_these_genes.bed`



#### D. Mask primary reference and annotation. 


- Mask HLA genes __only__ used for personalized mapping. 
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
rtracklayer::export(HLA_genes.masked, “GRCh38.primary.HLA_masked.gtf”)
```


#### E. Get final personalized genome and annotation 

* __For each individual__, 
- reference genome includes __hg38 HLA masked primary reference__ and __personalized HLA genes' sequences__.
- annotation includes __gtf of hg38 annotation excluding HLA genes__ and __gtf of personalized HLA genes__.



```R
# reference sequences (final)

ref_primary <- "GRCh38.primary.HLA_masked.fa" # generated from 2.3.
ref_personalized_hla <- "output_directory/out/individual_ID/individual_ID.per.fa" # generated from 2.2

refSeq.primary.masked <- readBStringSet(ref_primary)
refSeq.personalizedHLA <- readBStringSet(ref_personalized_hla)
refSeq.final <- c(refSeq.personalizedHLA, refSeq.primary.masked)
Biostrings::writeXStringSet(refSeq.final, 
                            "output_directory/out/individual_ID/individual_ID.primaryMasked_and_HLA.fa" )


# gtf (final)

gtf_primary <- "GRCh38.primary.HLA_masked.gtf" # generated from 2.3.
gtf_personalized_hla <- "output_directory/out/individual_ID/individual_ID.per.gtf" # generated from 2.2

refGTF.primary.masked <- rtracklayer::import(gtf_primary) 
refGTF.personalizedHLA <- rtracklayer::import(gtf_personalized_hla) 
refGTFfinal <- c(refGTF.personalizedHLA, refGTF.primary.masked)
rtracklayer::export(refGTF.final, 
                    "output_directory/out/individual_ID/individual_ID.primaryMasked_and_HLA.gtf"))

```


<br>
<br>

## 3. Mapping with Cell Ranger


### 3.1 Indexing: mkref


<br>

__*For each sample*,__

```bash
#!/bin/bash

cellranger \
        mkref \
        --nthreads=12 \
        --genome=output_directory/ref.individual_ID \
        --fasta=output_directory/out/individual_ID/individual_ID.primaryMasked_and_HLA.fa \
        --genes=output_directory/out/individual_ID/individual_ID.primaryMasked_and_HLA.gtf

```

### 3.2 Mapping: count


[here](../scripts/cellranger.count.sh)'s a script to run `cellranger count`.

```bash

cellranger.count.sh sample_name fastq_path path_of_personalised_ref out_prefix

# ./scripts/cellranger.count.sh individual_ID fastq_path output_directory/ref.individual_ID out_prefix
```


