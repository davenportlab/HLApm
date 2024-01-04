

indexing_star <- function(sample_directory="./", run_STAR="STAR"){
  
      alleles.fa <- list.files(path=sample_directory, pattern="fa$")
      alleles.gtf <- list.files(path=sample_directory, pattern="gtf$")
      
      alleles <- gsub("\\.fa", "", alleles.fa)
      print( paste0(length(alleles), " alleles found"))
      
      if ( all( alleles.gtf %in% paste0(alleles, ".gtf") != TRUE) ){
          stop("fa and gtf files are not matching...")
      }
      
      run_indexing <- c()
      run_mapping <- c()
      
      for (a in alleles){
            print(a)
            a.fa=paste0(sample_directory, "/", a, ".fa")
            a.gtf=paste0(sample_directory, "/", a, ".gtf")
            out_dir = paste0(sample_directory, "/idx.", a)
            
            if ( !file.exists(out_dir)){
              dir.create(file.path(out_dir ))
            }
            
            cmd=paste0(run_STAR, 
                 " --runThreadN 4 --runMode genomeGenerate ",
                 "--genomeSAindexNbases 8 ",
                 "--genomeDir ", out_dir, " ",
                 "--genomeFastaFiles ", a.fa, " ",
                 "--sjdbGTFfile ", a.gtf, " ",
                 "--sjdbOverhang 99 --outTmpDir ", out_dir, "/tmp")
            
            cmd2=paste0(run_STAR, )
            
            run_indexing <- rbind(run_indexing, cmd)
      }

      write.table(run_index, 
                  file=paste0(sample_directory, "/run_star_indexing.sh"), 
                  quote=FALSE,
                  row.names=FALSE, 
                  col.names=FALSE)
}




mapping_star_for_single_allele <- function(run_STAR="STAR", 
                         runThreadN=4, 
                         index_dir, 
                         fastq_1, fastq_2, outFileNamePrefix){
  
  # STAR="/software/team282/download/STAR-2.7.3a/bin/Linux_x86_64/STAR"
  # runThreadN=4
  # index_dir="/lustre/scratch119/realdata/mdt3/teams/davenport/ws/ref/IMGT_star_index/idx/chr_HLA-A_03_01_01_01" 
  # fastq_1="/lustre/scratch119/humgen/projects/gains_team282/hla/1.1.fq_for_mhc/gains8032857/gains8032857.forMHCmapping_1.fq.gz"
  # fastq_2="/lustre/scratch119/humgen/projects/gains_team282/hla/1.1.fq_for_mhc/gains8032857/gains8032857.forMHCmapping_2.fq.gz"
  # prefix="./"
  
        
      cmd1 <- paste0(run_STAR, " ",
                     "--genomeDir ", index_dir, " ",
                     "--readFilesIn ", fastq_1, " ", fastq_2, " ",
                     "--readFilesCommand zcat ",
                     "--outFileNamePrefix ", outFileNamePrefix, ". ",
                     "--outSAMtype BAM Unsorted ",
                     "--quantMode GeneCounts " )
      
      cmd2 <- paste0("samtools sort ", outFileNamePrefix, ".Aligned.out.bam -o ", outFileNamePrefix, ".Aligned.sorted.bam")
      cmd3 <- paste0("samtools index ", outFileNamePrefix, ".Aligned.sorted.bam")
      cmd.to.run.star <- rbind(cmd1, cmd2, cmd3)
  
  #tmp <- paste0("rm ", outFileNamePrefix, ".Aligned.out.bam")
  #cmd.to.run.star <- rbind(cmd.to.run.star, tmp)
      
      return(cmd.to.run.star)
      
}



merge_bams <- function(input_alleles, sample_name="individual_a", bam_directory="./" ){
  colnames(input_alleles) <- c("individual_ID", "HLA_allele")
  
  tmp <- input_alleles[input_alleles$individual_ID==sample_name,]
  tmp$HLA_allele <- gsub("(HLA\\-.*\\*\\d+\\:\\d+).*", "\\1", tmp$HLA_allele)
  tmp$HLA_allele <- gsub("\\*|\\:", "\\_", tmp$HLA_allele)
  tmp <- unique(tmp)
  
  bam_list <- c()
  for(a in tmp$HLA_allele){
      bam_list <- paste0(bam_list,
                         bam_directory, "/", a, ".Aligned.sorted.bam ") 
  }

  merge_bams <- paste0(bam_directory, "/", sample_name, ".initial_HLA_mapping.merged.bam")
  cmd1 <- paste0("samtools merge ",merge_bams, " ",  bam_list)
  system(cmd1)
  cmd2 <- paste0("samtools sort ", merge_bams, " -o ", gsub("merged", "merged.sorted", merge_bams) )
  system(cmd2)
  cmd3 <- paste0("samtools index ", gsub("merged", "merged.sorted", merge_bams) )
  system(cmd3)
  cmd4 <- paste0("rm ", merge_bams)
  system(cmd4)

}

  

### picard sorting
sort.by.picard.queryname <- function(bam, out_file_name){
  sort.bam=paste0("/usr/lib/jvm/java-1.8.0-openjdk-amd64/bin/java -Xmx1000m -jar /software/team282/jar/picard.jar ",
                  "SortSam ",
                  "I=", bam, " ",
                  "O=", out_file_name, " ",
                  "SORT_ORDER=queryname" )
  system(sort.bam)
}


