

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
    
    print("##################################################")
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
    
    sample.all.gtf <- list()
    i=0
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
          
          i=i+1
          tmp_gtf <- paste0(tmp_dir, "/", a, ".gtf") 
          write.table (out[[2]], tmp_gtf, sep="\t", quote=F, row.names = F, col.names=F)
          sample.all.gtf[[i]] <- rtracklayer::import( tmp_gtf )
          sample.all.gtf[[i]]$gene_id <- paste0("IMGT.", sample.all.gtf[[i]]$gene_name[1])
          sample.all.gtf[[i]]$transcript_id <- gsub("trans_", "", sample.all.gtf[[i]]$transcript_id)
          sample.all.gtf[[i]]$transcript_name <- gsub("transcript_", "", sample.all.gtf[[i]]$transcript_name)
          
          if(i==1){
                sample.all.seq <- out[[1]]
          }
          if(i>1){
                sample.all.seq <- c(sample.all.seq, out[[1]])
          }
          
          # remove gtf
          file.remove(tmp_gtf)
    }
    
    ## add DRB_paralogous
    #if(drb_hap=="no_DRB_paralogous_genes"){}
    
    #if(drb_hap!="no_DRB_paralogous_genes"){
    if(length(drb_hap)==2){
      i=i+1
      sample.all.seq <-  c(sample.all.seq, drb_hap[[1]] )
      sample.all.gtf[[i]] <- drb_hap[[2]]
    }
    
    sample.all.gtf <- do.call(c, sample.all.gtf)
    
    
    print("--------------------------------------------------")
    print("##################################################")
    print("### Combining alleles seqeunces and annotation")
    
    # combine them all
    new_seq_name=paste0("chr_personalised.", individual_ID)
    tmp <- as.data.frame(sample.all.gtf)
    
    all_pers_alleles <- names(sample.all.seq)[order(names(sample.all.seq))]
    
    comb.sample.all.seq <- sample.all.seq[names(sample.all.seq)==all_pers_alleles[1]]
    comb.sample.all.gtf <- sample.all.gtf[tmp$seqnames==all_pers_alleles[1] ]
    
    for( chr in all_pers_alleles[2:length(all_pers_alleles)] ){
          print(chr)
          s = sample.all.seq[names(sample.all.seq)==chr]
          len_fa <- nchar(comb.sample.all.seq)
          
          g = sample.all.gtf[tmp$seqnames==chr ]
          g@ranges@start <- as.integer( g@ranges@start + len_fa + 10 )
          
          comb.sample.all.seq <- paste0(comb.sample.all.seq, "NNNNNNNNNN" , s)
          comb.sample.all.gtf <- c(comb.sample.all.gtf, g)
    }
    comb.sample.all.gtf@seqnames@values <- 
      rep(new_seq_name, length(comb.sample.all.gtf@seqnames@values ) )
 

    #tmp <- as.data.frame(comb.sample.all.gtf)
    
    suppressWarnings ( rtracklayer::export ( comb.sample.all.gtf,
                          paste0(tmp_dir, "/", individual_ID, ".per.gtf") ) )
    
    ## adjust gene ranges
    gtf <- rtracklayer::import( paste0(tmp_dir, 
                                       "/", individual_ID, ".per.gtf") )
    
    for(gene in unique(gtf$gene_name)){
          #print(gene)
          ir <- gtf[gtf$gene_name==gene & gtf$type=="gene"]
          if(  length(ir) ==2 ) {
                end_pos <- max( c( (ir@ranges@start[1] + ir@ranges@width[1]),
                               (ir@ranges@start[2] + ir@ranges@width[2] ) ) )
                wid <- as.integer (end_pos - min(ir@ranges@start)  + 10  )
                ir[ir$gene_name==gene & ir$type=="gene"]@ranges@start <- 
                  rep( min(ir@ranges@start), 2)
                ir@ranges@width  <- rep(wid, 2)
                
                tmp1 <- gtf[gtf$gene_name!=gene | gtf$type!="gene"]
                gtf <- c(ir[1], tmp1)
          }
    }
    gtf <- gtf[order(gtf@ranges@start)]
    suppressWarnings ( rtracklayer::export ( gtf,
                                             paste0(tmp_dir, "/", individual_ID, ".per.gtf") ) )
    
    
    # Biostrings::writeXStringSet(sample.all.seq, 
    #                             paste0(tmp_dir, "/", individual_ID, ".per.fa") )
    write.fasta(sequences=comb.sample.all.seq, 
                names=new_seq_name,
                file.out=paste0(tmp_dir, "/", individual_ID, ".per.fa") )
    
}


