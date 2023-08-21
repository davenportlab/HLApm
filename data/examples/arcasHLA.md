### arcasHLA conda install

```bash
conda install -c bioconda arcas-hla
```


### prepare reference

You can select the version you like using the commithash from the [IMGT/HLA Github](https://github.com/ANHIG/IMGTHLA).

```bash
arcasHLA reference --version [commithash]
```


### Run arcasHLA


```bash

#!/bin/bash

## activate conda before running arcasHLA
# conda activate arcasHLA

# $1, samplename
# $2, input fastq 
# $3, output_dir
# $4 threads

arcasHLA genotype \
        --single -t $4 --drop_iterations 4 \
        -o $3arcasHLA.$1 \
        -v --log $3/arcasHLA.$1.log \
        -g A,B,C,DPA1,DPB1DQA1,DQB1,DRB1,DRB3,DRB4,DRB5 \
        $2
```
