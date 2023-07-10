

## make_personalized_HLA_ref: 
## inputs:
# input_alleles : input alleles, see `./data/examples/example_input.txt`
# output_directory: directory to write output
## outputs:
# individual_ID.fa : HLA reference sequences
# individual_ID.gtf : adjusted HLA annotation

build_personalized_HLA_ref <- function(input_alleles, output_directory="./" ){
    
    colnames(input_alleles) <- c("individual_ID", "HLA_allele")

    # number of samples
    samples <- unique( gsub("(.*)\\*.*", "\\1", input_alleles$individual_ID) )
    print(paste0(length(samples),  " individuals found" ))
    HLA_genes <- unique( gsub("(.*)\\*.*", "\\1", input_alleles$HLA_allele) )
    print(paste0(length(HLA_genes),  " HLA genes found" ))
    
    print("### important: ")
    print(paste0( (length(HLA_genes)+1), " HLA genes (HLA-DRB5, ", 
                paste0(HLA_genes, collapse = ", "),  
                ") need to be masked in primary reference and annotation " ))
  
    mask_genes <- c(HLA_genes, "HLA-DRB5")
    mask_genes <- rbind( hla_to_mask[hla_to_mask$V4 %in% mask_genes, ],
                         hla_to_mask[grep("DRB", hla_to_mask$V4), ] )
    mask_genes <- unique(mask_genes)
    
    tmp_dir <- paste0(output_directory, "/out/")
    if (!file.exists(tmp_dir)){
      dir.create(file.path(output_directory, "/out/"), 
                 recursive=T)
    }
    write.table(mask_genes,
                paste0(output_directory, "/out/mask_these_genes.bed"),
                sep="\t", quote=F, row.names = F)
    
    for(s in samples){
            per_sample(input_alleles, individual_ID=s, output_directory)
    }
      
}


## get_DRB_haplotypes
# this is to make DRB haplotypes based on DRB1 alleles
# inputs
#   input_alleles: 
#   individual_ID: name of sample(individual), should be the same name from input_alleles
# outputs
#   


get_DRB_haplotypes <- function(input_alleles, individual_ID=""){
  
  if(individual_ID==""){
      print("no sample name found")
  }
  
  drb1 <- input_alleles[grep("DRB1", input_alleles$HLA_allele),]
  drb1 <- drb1[drb1$individual_ID==individual_ID, c("HLA_allele")]
  drb1 <- unique( gsub("HLA-DRB1\\*(\\d{2})\\:\\d+.*", "\\1", drb1) )
  
  other_drbs <- drb_haplotypes[drb_haplotypes$DRB1 %in% drb1,]$DRB_par
  other_drbs <- unique( unlist( strsplit( other_drbs, ";") ) )
  #print(other_drbs)
  other_drbs <- other_drbs[order(other_drbs)]
  protein_coding_drbs <- unique( other_drbs[(other_drbs %in% c("DRB3", "DRB4", "DRB5"))] )
  print( paste0( "DRB paralogous genes: ", paste0(protein_coding_drbs, collapse = ", ") ) )
  
  
  if(length(protein_coding_drbs) >0){
      per.drb_fa <- DRB_fa[names(DRB_fa) %in% c( paste0("chr_HLA-", protein_coding_drbs) ) ]
      per.drb_gtf <- DRB_gtf[DRB_gtf$gene_name %in% c( paste0("HLA-", protein_coding_drbs) ), ]
  
      out_drb <- list(per.drb_fa, per.drb_gtf)
  }
  if(length(protein_coding_drbs) ==0){
    print("no DRB paralogous genes")
    out_drb="no_DRB_paralogous_genes"
  }
  return(out_drb)
}



per_sample <- function(input_alleles, individual_ID="", output_directory="./"){
    
    if(individual_ID==""){
      print("no sample name found")
    }
  
    print(paste0("sample name (individual_ID): ", individual_ID ) )
    
    tmp_dir <- paste0(output_directory, "/out/", individual_ID)
    if (!file.exists(tmp_dir)){
        dir.create(file.path(output_directory, "/out/", individual_ID), 
                   recursive=T)
    }
    
    drb_hap <- get_DRB_haplotypes(input_alleles, individual_ID=individual_ID)
    
    per.alleles <- unique( input_alleles[input_alleles$individual_ID==individual_ID, c("HLA_allele")] )
    
    sample.all.seq <- drb_hap[[1]]
    sample.all.gtf <- list()
    sample.all.gtf[[1]] <- drb_hap[[2]]
    i=1
    for(a in per.alleles){
          out <- align_and_adjust_annotation (a, output_directory=tmp_dir, prefix= individual_ID)
          i=i+1
          tmp_gtf <- paste0(tmp_dir, "/", a, ".gtf") 
          write.table (out[[2]], tmp_gtf, sep="\t", quote=F, row.names = F, col.names=F)
          sample.all.gtf[[i]] <- rtracklayer::import( tmp_gtf )
          sample.all.gtf[[i]]$gene_id <- paste0("IMGT.", sample.all.gtf[[i]]$gene_name[1])
          sample.all.gtf[[i]]$transcript_id <- gsub("trans_", "", sample.all.gtf[[i]]$transcript_id)
          sample.all.gtf[[i]]$transcript_name <- gsub("transcript_", "", sample.all.gtf[[i]]$transcript_name)
          
          sample.all.seq <- c(sample.all.seq, out[[1]])
          file.remove(tmp_gtf)
    }
    
    # combine them all
    suppressWarnings ( rtracklayer::export ( do.call(c, sample.all.gtf),
                          paste0(tmp_dir, "/", individual_ID, ".per.gtf") ) )
    
    
    Biostrings::writeXStringSet(sample.all.seq, 
                                paste0(tmp_dir, "/", individual_ID, ".per.fa") )
    
}



