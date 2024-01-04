
bulkRNA_per_sample <- function(input_alleles, individual_ID="", output_directory="./"){
  
  if(individual_ID==""){
    print("no sample name found")
  }
  
  print(paste0("sample name (individual_ID): ", individual_ID ) )
  
  tmp_dir <- paste0(output_directory, "/out/", individual_ID)
  if (!file.exists(tmp_dir)){
    dir.create(file.path(output_directory, "/out/", individual_ID), 
               recursive=T)
  }
  
  print("##################################################")
  print(paste0("### DRB paralogous genes of sample '", individual_ID, "'"))
  drb_hap <- get_DRB_haplotypes(input_alleles, individual_ID=individual_ID)
  
  
  print("##################################################")
  print(paste0("### Generating allele sequences and adjusting gene annotation"))
  
  per.alleles <- unique( input_alleles[input_alleles$individual_ID==individual_ID, c("HLA_allele")] )
  per.alleles.2f <- gsub("(HLA\\-[A|B|C|DRB1|DQA1|DQB1|DPA1|DPB1|DOA|DOB|K|G|E|F]+\\*\\d+\\:\\d+).*", "\\1", 
                         per.alleles)
  per.alleles.2f <- unique(per.alleles.2f)
  #all_genes <- unique( gsub("(.*)\\*.*", "\\1", per.alleles.2f) )
  

  for(a in per.alleles.2f){
    print("--------------------------------------------------")
    print(paste0("### allele: ", a))
    
    allele_level_2 <- allele_level
    allele_level_2$allele <- gsub("(.*\\d+)[A-Z]$", "\\1", 
                                  allele_level_2$allele)
    
    find_level <- allele_level_2[ grep(paste0( gsub("[\\*|\\:]", "_", a), "_" ), 
                                       paste0(allele_level_2$allele, "_") ) ,]
    
    print(paste0( "Allele level: ", find_level) )
    
    if(find_level$level=="L1" ){
      out <- align_and_adjust_annotation (a, output_directory=tmp_dir, prefix= individual_ID)
    }
    if(find_level$level!="L1"){
      out <- align_and_adjust_annotation_L234 (a, output_directory=tmp_dir, prefix= individual_ID)
    }
    

    tmp_gtf <- paste0(tmp_dir, "/", a, ".gtf") 
    write.table (out[[2]], tmp_gtf, sep="\t", quote=F, row.names = F, col.names=F)
    
    sample.gtf <- rtracklayer::import( tmp_gtf )
    sample.gtf$gene_id <- paste0("IMGT.", sample.gtf$gene_name[1])
    sample.gtf$transcript_id <- gsub("trans_", "", sample.gtf$transcript_id)
    sample.gtf$transcript_name <- gsub("transcript_", "", sample.gtf$transcript_name)
    
    rename_a <- gsub("\\:|\\*", "\\_", a)
    suppressWarnings ( rtracklayer::export ( sample.gtf,
                                             paste0(tmp_dir, "/", rename_a, ".gtf") ) )
    
    file.remove(tmp_gtf)
    
    write.fasta(sequences=out[[1]], 
                names=names(out[[1]]),
                file.out=paste0(tmp_dir, "/", rename_a, ".fa") )
  }
  
  
  
  ## add DRB_paralogous
  #if(drb_hap=="no_DRB_paralogous_genes"){}
  
  #if(drb_hap!="no_DRB_paralogous_genes"){
  if(length(drb_hap)==2){
      for(drb_name in names(drb_hap[[1]]) ){
          drb_seq = drb_hap[[1]][names(drb_hap[[1]])==drb_name]
          drb_gtf = drb_hap[[2]][seqnames(drb_hap[[2]])==drb_name]
          
          rename_a <- gsub("chr_", "", drb_name)
          suppressWarnings ( rtracklayer::export ( drb_gtf,
                                                   paste0(tmp_dir, "/", rename_a, ".gtf") ) )
          
          write.fasta(sequences=drb_seq, 
                      names=drb_name,
                      file.out=paste0(tmp_dir, "/", rename_a, ".fa") )
      }

  }
  
}





bulkRNA_build_personalized_HLA_ref <- function(input_alleles, output_directory="./" ){
  
  colnames(input_alleles) <- c("individual_ID", "HLA_allele")
  
  # number of samples
  samples <- unique( gsub("(.*)\\*.*", "\\1", input_alleles$individual_ID) )
  print(paste0(length(samples),  " individuals found" ))
  HLA_genes <- unique( gsub("(.*)\\*.*", "\\1", input_alleles$HLA_allele) )
  print(paste0(length(HLA_genes),  " HLA genes found" ))
  
  
  print("##################################################")
  print("### important: ")
  print(paste0( (length(HLA_genes)+1), " HLA genes (HLA-DRB5, ", 
                paste0(HLA_genes, collapse = ", "),  
                ") need to be masked in primary reference and annotation " ))
  
  mask_genes <- c(HLA_genes, "HLA-DRB5")
  mask_genes <- rbind( hla_to_mask[hla_to_mask$V4 %in% mask_genes, ],
                       hla_to_mask[grep("DRB", hla_to_mask$V4), ] )
  mask_genes <- unique(mask_genes)
  
  
  print(paste0("### make a directory for personalized referece : ") )
  tmp_dir <- paste0(output_directory, "/out/")
  print(tmp_dir)
  
  if (!file.exists(tmp_dir)){
    dir.create(file.path(output_directory, "/out/"), 
               recursive=T)
  }
  
  write.table(mask_genes,
              paste0(tmp_dir, "/mask_these_genes.bed"),
              sep="\t", quote=F, row.names = F)
  
  
  for(s in samples){
    bulkRNA_per_sample(input_alleles, individual_ID=s, output_directory)
  }
  
}



