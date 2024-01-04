
star_mapping_to_hla_genes <- function(STAR, runThreadN=4, index_dir, fastq_1, fastq_2, outFileNamePrefix){
	# STAR="/software/team282/download/STAR-2.7.3a/bin/Linux_x86_64/STAR"
	# runThreadN=4
	# index_dir="/lustre/scratch119/realdata/mdt3/teams/davenport/ws/ref/IMGT_star_index/idx/chr_HLA-A_03_01_01_01" 
	# fastq_1="/lustre/scratch119/humgen/projects/gains_team282/hla/1.1.fq_for_mhc/gains8032857/gains8032857.forMHCmapping_1.fq.gz"
	# fastq_2="/lustre/scratch119/humgen/projects/gains_team282/hla/1.1.fq_for_mhc/gains8032857/gains8032857.forMHCmapping_2.fq.gz"
	# prefix="./"


	cmd.to.run.star <- paste0("/software/team282/download/STAR-2.7.3a/bin/Linux_x86_64/STAR --runMode alignReads --runThreadN 4 ",
                                                  "--genomeDir ", index_dir, " ",
                                                  "--readFilesIn ", fastq_1, " ", fastq_2, " ",
                                                "--readFilesCommand zcat ",
                                                "--outFileNamePrefix ", outFileNamePrefix, ". ",
                                                "--outSAMtype BAM Unsorted ",
                                                "--quantMode GeneCounts " )

	tmp <- paste0("samtools sort ", outFileNamePrefix, ".Aligned.out.bam -o ", outFileNamePrefix, ".Aligned.sorted.bam")
	cmd.to.run.star <- rbind(cmd.to.run.star, tmp)
	tmp <- paste0("samtools index ", outFileNamePrefix, ".Aligned.sorted.bam")
	cmd.to.run.star <- rbind(cmd.to.run.star, tmp)
	#tmp <- paste0("rm ", outFileNamePrefix, ".Aligned.out.bam")
	#cmd.to.run.star <- rbind(cmd.to.run.star, tmp)

	return(cmd.to.run.star)
}



featureCounts.mhc <- function(input.bam, gtf, outFileNamePrefix){
	# gtf="/lustre/scratch119/humgen/teams/davenport/ws/ref/GRCh38.99/gencode.v33.annotation.gtf.gz"
	out.bam <- paste0(outFileNamePrefix, ".mhc.bam")

	cmd1=paste0("samtools view --threads 1 -b -o ",
            out.bam, " ",  input.bam, " 6:28500000-33400000" )
	cmd2=paste0("samtools index ", out.bam)
	cmd3= paste0("featureCounts -R BAM -T 1 -a ", 
		     gtf, " ",
		     "-g gene_name  -o ",
		     outFileNamePrefix, ".gene_name ",
		     "-p -s 2 ",  out.bam)

	cmd <- rbind(cmd1, cmd2, cmd3)
	return(cmd)
}	

featureCounts.GL000255 <- function(input.bam, gtf, outFileNamePrefix){
	gtf <- "/lustre/scratch119/realdata/mdt3/teams/davenport/ws/ref/GRCh38.99/by_chrs/GL000255.gtf"
	cmd <- paste0("featureCounts -R BAM -T 1 -a ",
		      gtf, " ",
		      "-g gene_name  -o ",
		      outFileNamePrefix, ".gene_name ",
		      "-p -s 2 ",  input.bam)
	return(cmd)
}


featureCounts.GL000256 <- function(input.bam, gtf, outFileNamePrefix){
        gtf <- "/lustre/scratch119/realdata/mdt3/teams/davenport/ws/ref/GRCh38.99/by_chrs/GL000256.gtf"
        cmd <- paste0("featureCounts -R BAM -T 1 -a ",
                      gtf, " ",
                      "-g gene_name  -o ",
                      outFileNamePrefix, ".gene_name ",
                      "-p -s 2 ",  input.bam)
        return(cmd)
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



read_featureCountsBam <- function(bam_with_featureCounts){
	library(Rsamtools)
	library(dplyr)

	# bam_with_featureCounts <- paste0(OUT_DIR, "/", s, "/", s, ".mhc.bam.featureCounts.bam")
	what <- c("qname" )
	tag=c("NH", "nM", "XS", "XT")
	param=ScanBamParam(tag=tag, what=what)
	featureCounts.bam <- scanBam(bam_with_featureCounts, param=param)
	featureCounts <- data.frame(qname=featureCounts.bam[[1]]$qname,
                            NH=featureCounts.bam[[1]]$tag$NH,
                            NM=featureCounts.bam[[1]]$tag$nM,
                            XS=featureCounts.bam[[1]]$tag$XS,
                            XT=featureCounts.bam[[1]]$tag$XT)
	dim(featureCounts)
	featureCounts.assigned <- featureCounts[featureCounts$XS=="Assigned",]
	# filtering un-paried reads
	tmp <- featureCounts.assigned  %>% count(qname, XT)
        tmp <- tmp[tmp$n==2,]
	featureCounts.assigned <- featureCounts.assigned[featureCounts.assigned$qname %in% tmp$qname,]

	return(featureCounts.assigned)
}

reads_from_bam <- function(merged_bams){
	library(Rsamtools)

	what <- c("qname" , "rname")
	tag=c("NH", "nM")
	param=ScanBamParam(tag=tag, what=what)
	bam <- scanBam(merged_bams, param=param)
	bam <- data.frame(rname=bam[[1]]$rname,
			  qname=bam[[1]]$qname,
			  NH=bam[[1]]$tag$NH,
			  NM=bam[[1]]$tag$nM)
	return(bam)
}



#java -jar picard.jar FilterSamReads \
#       I=input.bam \ 
#       O=output.bam \
#       READ_LIST_FILE=read_names.txt FILTER=filter_value

filter_bam_by_qnames <- function(input.bam, output.bam, READ_LIST_FILE){
	cmd=paste0("/usr/lib/jvm/java-1.8.0-openjdk-amd64/bin/java -Xmx1000m -jar /software/team282/jar/picard.jar ",
                        "FilterSamReads ",
                        "I=", input.bam, " ",
                        "O=", output.bam, " ",
                        "READ_LIST_FILE=", READ_LIST_FILE,
		       " FILTER=includeReadList SORT_ORDER=coordinate"	)
        system(cmd)
}



# samtools merge merged.bam HLA-A_03_01.bam HLA-DQB2.bam
# samtools sort merged.bam -o merged.sorted.bam 
# samtools index merged.sorted.bam 
# samtools view merged.sorted.bam chr_HLA-A_03_01_01_01 -b > test.bam



#find_missing_hla <- function(sample_directory, ){
#}




