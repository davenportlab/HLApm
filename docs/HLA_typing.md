# HLA typing 


### with array data

- [HIBAG](https://bioconductor.org/packages/release/bioc/html/HIBAG.html), there are [pre-fit models](https://hibag.s3.amazonaws.com/hlares_index.html) for different ancestry / Multi-ethnic models specific to a genotyping platform. * To generate a HIBAG ref. panel, filter SNPS that overlapping with your array data in a corresponding genome build.

- [SNP2HLA](http://software.broadinstitute.org/mpg/snp2hla/)





### with scRNA-seq


- [arcasHLA](https://github.com/RabadanLab/arcasHLA)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;__- make sure that you use the same version of IMGT/HLA data for running arcasHLA__

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;__- add `--single`flag for single-end reads (only trascript sequences as input, R2)__

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;__- here's [simple procedure](../data/examples/arcasHLA.md) how to run arcasHLA.__



### with bulk-RNA seq
- [arcasHLA](https://github.com/RabadanLab/arcasHLA) 

- [PHLAT](https://sites.google.com/site/phlatfortype) 
* It can also use DNA sequences as Input.
* To obtain access, please contact them at [regn.phlat@gmail.com](mailto:regn.phlat@gmail.com). A valid email with academic domain name is needed



### with WES, WGS

- [HLA-LA](https://github.com/DiltheyLab/HLA-LA) 

- [HLAscan](https://github.com/SyntekabioTools/HLAscan)

- [Kourami](https://github.com/Kingsford-Group/Kourami) 

- [Michigan Imputation Server](https://imputationserver.readthedocs.io/en/latest/): HLA reference panel spanning five global populations based on ~20,000 whole-genome sequencing data. * Based on my personal experience, they are probably good at G-group resolution, but may have low accuracy at 2-field resolution comparing to output from HIBAG.
