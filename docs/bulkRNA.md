


## 1. input data


<br>

## 2. Set a path of a package

* *Set a path of this package with a name __`DIR_personalizedHLA`__*
* a directory of `scripts` should be under the `DIR_personalizedHLA`

```R
# where the package installed
yourpath="./"
DIR_personalizedHLA=paste0(yourpath, "/personalisedHLAmapping") 
```

<br>

## 3. Load functions


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

## 4. Gernerate personalized HLA reference and annotation for each individual


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

<br>

## 5. build star indexing 


### Input
- sample_directory: where all `*.fa` and `*.gtf` files for alleles (output from `bulkRNA_build_personalized_HLA_ref`) for each sample
- run_STAR: a path where `STAR` installed


### Output

`run_star_index.sh` : command lines to build index for all alleles in `sample_directory`. 



### Example
```R
dir_test="/lustre/scratch126/humgen/teams/davenport/ws/github_ws/HLApm/scripts/bulk/out/individual_a"
star_link="/software/team282/download/STAR-2.7.3a/bin/Linux_x86_64/STAR"

indexing_star(sample_directory=dir_test, run_STAR=star_link)
```


__run STAR indexing__

```bash
sh run_star_indexing.sh 
```

<br>

## 6. STAR mapping


### Input

- run_STAR: a path where `STAR` installed
- runThreadN: NumberOfThreads 
- index_dir: a directory path where star index files are 
- fastq_1 and fastq_2 : /path/to/read1 and /path/to/read2
- outFileNamePrefix : /path/to/output/dir/prefi


### Output

- commands to run star mapping for each allele



### Example
`mapping_star_for_single_allele` runs for each allele. 

You can run all alleles for each sample with combined commands. 

```R
# for each sample

file_to_write <- "textFile.sh"

run_star_mapping <- c()
for(a in alleles){
    index_dir_for_allele = path/of/star_index/for_each_allele

    cmd.to.run.star <- mapping_star_for_single_allele(run_STAR="STAR", 
                                                  runThreadN=4, 
                                                  index_dir=index_dir_for_allele, 
                                                  fastq_1, fastq_2, 
                                                  outFileNamePrefix)
                                                  
    run_star_mapping <- rbind(run_star_mapping, cmd.to.run.star)
}
writeLines(run_star_mapping, file_to_write)

```

__submit bsub__
```bash
bsub  ~~~ bash textFile.sh
```


<br>

## 7. merge bam files for each sample

### Input
input_alleles: Alleles for each sample. Please see [example input](../data/examples/example_input.txt)
sample_name: single sample from input_alleles
bam_directory: where individual bam files are (output from `mapping_star_for_single_allele`)

### Output

`sample_name.initial_HLA_mapping.merged.sorted.bam`


### Example
```R
merge_bams (input_alleles, sample_name="individual_a", bam_directory="./" )
```




# delete this later

/lustre/scratch123/hgi/teams/davenport/ws/GitHub/HLA_mapping/gains/alt.4.2.mhc_mapping
/lustre/scratch126/humgen/teams/davenport/ws/github_ws/HLApm/scripts/bulk/out/individual_a
/lustre/scratch126/humgen/projects/gains_team282/hla/alt.4.2.mhc_mapping/gains8033273
/lustre/scratch126/humgen/teams/davenport/ws/github_ws/HLA_modules/HLA_mapping



