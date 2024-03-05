

# 1.Build HLA personalized reference. 


## 1.1. Set a path of a package

* *Set a path of this package with a name __`DIR_personalizedHLA`__*
* a directory of `scripts` should be under the `DIR_personalizedHLA`

```R
# where the package installed
yourpath="./"
DIR_personalizedHLA=paste0(yourpath, "/personalisedHLAmapping") 
```

<br>

## 1.2. Load functions


Just running a single script `build_personalized_HLA_ref` can provide personalized HLA reference. 

```R
source( paste0(DIR_personalizedHLA, "/scripts/load_ref.R") )
source( paste0(DIR_personalizedHLA, "/scripts/make_personalized_HLA_ref.R") )
source( paste0(DIR_personalizedHLA, "/scripts/align_and_adjust_annotation.R") )
source( paste0(DIR_personalizedHLA, "/scripts/bulk/functions_personalized_HLA_ref.R") )
source( paste0(DIR_personalizedHLA, "/scripts/bulk/functions_star.R") )


```

All resource for building references loaded ! 


<br>

## 1.3. Gernerate personalized HLA reference and annotation for each individual


### Input

- input_alleles: Alleles for each sample. Please see [example input](../data/examples/example_input.txt)
- output_directory


### Output
a pair of gtf and annotation files for each HLA allele in a sample directory 


### Example

```R
# reading input alleles
input_alleles <- read.table(paste0(DIR_personalizedHLA, 
                                  "/data/examples/example_input.txt"), head=T)

# Get personalized HLA reference and annotation
bulkRNA_build_personalized_HLA_ref(input_alleles, output_directory="./" )

```

Please see [example code](../scripts/bulk/1.build_HLApm.R)


<br>

# 2. Personalised mapping

__Three steps requred.__ 

- 1. Building star index for each allele reference
- 2. Mapping with STAR for each allele
- 3. Assigning reads to best allele for each sample  

<br>
__For each sample,__ <br>
__Run 1 and 2 for each allele, and run 3.__

## 2.1. Build star index


### Input

For each allele,<br>
 `*.fa` and `*.gtf` files (output from `bulkRNA_build_personalized_HLA_ref`) 


### Output

STAR index

### Example


Please see [example code](../scripts/bulk/2.1.build_star_index.sh) to run

__run STAR indexing__

```bash
sh 2.1.build_star_index.sh 1 output_directory A_01_01.fa A_01_01.gtf
```

<br>

## 2.2. STAR mapping


### Input

For each allele,

- $1, runThreadN : NumberOfThreads
- $2, genomeDir : /path/to/output/dir/genome_index, output from 2.1.build_star_index
- $3, fastq_1 : /path/to/read1
- $4, fastq_2 : /path/to/read2
- $5, outFileNamePrefix /path/to/output/dir/prefix


### Output

sorted bam and bam index file for each allele


### Example
Please see [example code](../scripts/bulk/2.2.star_mapping.sh) to run




<br>

## 2.3. Assign reads to allele






# delete this later

/lustre/scratch123/hgi/teams/davenport/ws/GitHub/HLA_mapping/gains/alt.4.2.mhc_mapping
/lustre/scratch126/humgen/teams/davenport/ws/github_ws/HLApm/scripts/bulk/out/individual_a
/lustre/scratch126/humgen/projects/gains_team282/hla/alt.4.2.mhc_mapping/gains8033273
/lustre/scratch126/humgen/teams/davenport/ws/github_ws/HLA_modules/HLA_mapping



