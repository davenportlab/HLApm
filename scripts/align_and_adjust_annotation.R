


#######################################################
### example 1

# input data
#pattern_seq = "./ref_tx_seq/seq_HLA-B.fa"
#allele_name="HLA-B*15:01"
#allele_name="HLA-B*08:01"	#del case
#allele_name="HLA-B*57:03"
#anno_file="./ref_tx_seq/anno_HLA-B.txt"
#anno <- read.table(anno_file, sep="\t", header=T)
#output_directory <- "./imgt_ref/"



# awk '$0 ~ "^>" { match($1, /^>([^:]+)/, id); filename=id[1]} {print >> filename".fa"}' input.fa
# load chr6 sequences


############ note
##### use '2.2.extend_allele_seq_L2.R' 
#HLA-DOB*01:02
#HLA-DOB*01:04

# allele_name= "HLA-B*07:02"



align_and_adjust_annotation <- function(allele_name, output_directory="./", prefix="sample" ){
  
  print(paste0("Input allele: ", allele_name) )
	# allele_name: name of  subject_seq 
	# pattern_anno = annotation file of pattern seq (primary ref seq)

	#############################################
	###
	### 1. extend allele sequence
	###
	#############################################

  gene <- gsub("(HLA-.*)\\_\\d+\\_\\d+", "\\1", allele_name)	
  gene <- gsub("^([A-Z].*)\\*\\d+.*", "\\1", allele_name)
  print( paste0("HLA gene: ", gene) )
  
  
  anno <- pattern_anno[pattern_anno$gene_name==gene,]
    
	if(gene=="HLA-A"){
		anno <- subset(anno, transcript_name=="HLA-A-202")
	}
	if(gene=="HLA-C"){
                anno <- subset(anno, transcript_name=="HLA-C-201")
        }
	if(gene=="HLA-DMA"){
                anno <- subset(anno, transcript_name=="HLA-DMA-202")
        }
	if( gene=="HLA-DQA1"){
		# need to fix UTR region
                anno <- subset( anno, transcript_name=="HLA-DQA1-201")
        }
	if( gene=="HLA-DRA"){
		anno <- subset( anno, transcript_name=="HLA-DRA-201")
	}
	if( gene=="HLA-DQB1"){
                anno <- subset( anno, transcript_name=="HLA-DQB1-201")
        }
	
	
	strand <- as.character(anno$strand[1])
	allele_name <- gsub("\\*", "\\_", allele_name)
	allele_name <- gsub("\\:", "\\_", allele_name)
	allele_name.2f <- gsub("(HLA\\-[A|B|C|DRB1|DQA1|DQB1|DPA1|DPB1|DOA|DOB|K|G]+\\_\\d+\\_\\d+).*", "\\1", 
	                       allele_name)
	
	subject_seq_names <- gsub ( "[N|Q|L|A]$", "", names(subject_seq) )
	
	if ( length( grep( allele_name.2f, subject_seq_names ) ) >0 ){
	    allele_name <- names(subject_seq)[grep( allele_name.2f, subject_seq_names ) ][1]
	}
	print(paste0( "Subject allele matched to an input allele: ", allele_name.2f) )

	seq1 <- pattern_seq.all[grep(gene, names(pattern_seq.all) )]
	seq2<- subject_seq[names(subject_seq) %in% allele_name ]

	seq1string <- DNAString(toupper(c2s(seq1[[1]])) )
	seq2string <- DNAString(toupper(c2s(seq2[[1]])) )
# 	print(paste0("Length of pattern allele : ", length(seq1string) ))
#   print(paste0("Length of subject allele : ", length(seq2string) ))


	#######################################################
  ### n. of exons from primary ref.
  #exons_in_pri_ref <- anno[anno$type=="exon", ]
  
  
	### get the first exon sequence
	gtf_sub <- subset(imgt_gtf, V1==allele_name)
  
  subject_exons <- subset(gtf_sub, V3=="exon")
  e <- min(subject_exons$V4)
	e_seq <- substr(toupper(c2s(seq2[[1]])), e, (e+15) )
	
	# first exon seq
	first_exon <- subject_exons[subject_exons$V4==min(subject_exons$V4), ]
	first_exon_seq <- substr(seq2string, first_exon$V4, first_exon$V5)
	
	# last exon seq
	last_exon <- subject_exons[subject_exons$V4==max(subject_exons$V4), ]
	last_exon_seq <- substr(seq2string, last_exon$V4, last_exon$V5)
	

	
	seqs <- DNAStringSet( c(toupper(c2s(seq1[[1]]) ) , toupper(c2s(seq2[[1]]) ) ) )

	if(strand=="-"){
		# reversing sequences and/or complementing DNA
		#seq2string <- reverseComplement(seq2string)
		#toString(seq2string)
		#seqs <- OrientNucleotides(seqs, verbose = FALSE)		# DECIPHER
		seqs <- DNAStringSet( c(toupper(c2s(seq1[[1]]) ) , 
					 toString(reverseComplement(DNAString(toupper(c2s(seq2[[1]]) ))  )) ) )
		first_exon_seq <- reverseComplement( first_exon_seq )
		last_exon_seq <- reverseComplement( last_exon_seq )
	}

	
	
	print(paste0( "Aligning ", allele_name, "(subject allele) on primary gene (pattern allele) sequences...") )
	aligned <- DECIPHER::AlignSeqs(seqs, iterations = 5, verbose = FALSE)
	
	
	
	# seqs_1 <- DNAStringSet( c(toupper(c2s(seq1[[1]]) ) ,
	#                           toString( first_exon_seq ) ) )
	# aligne_first_exon_seq <- AlignSeqs(seqs_1, iterations = 5, verbose = TRUE)
	# # toString(aligne_first_exon_seq[1] )
	# # toString(aligne_first_exon_seq[2] )
	# 
	# seqs_2 <- DNAStringSet( c(toupper(c2s(seq1[[1]]) ) ,
	#                           toString( last_exon_seq ) ) )
	# aligne_last_exon_seq <- AlignSeqs(seqs_2, iterations = 5, verbose = TRUE)
	# tmp_seq <- gsub("", "", toString(aligne_last_exon_seq[2] ) ) 
	# get_right_extra_seq <-  toString(aligne_last_exon_seq[2] )
	# # toString(aligne_last_exon_seq[1] )
	# # toString(aligne_last_exon_seq[2] )
	

	
	
	#alignedDNA <- AlignProfiles( DNAStringSet( c(toupper(c2s(seq1[[1]]) ))),  DNAStringSet( c(toupper(c2s(seq2[[1]]) ) ) ) ) 

	### save aligned sequences
	#writeXStringSet(aligned, file="test.fa")

	
	aligned_pattern <- toString(aligned[1] )
  aligned_subject <- toString(aligned[2] )
	nchar(aligned_pattern)
	nchar(aligned_subject)
	
	aligned_subject <- gsub("(A|C|G|T)+", "0", aligned_subject)


	### check if there is deletion in aligned_pattern.
	check_del_in_pattern <- stringr::str_count(aligned_pattern, '-')
	if(check_del_in_pattern >0){
		      # print(paste0(check_del_in_pattern, " deletion(s) found in aligned_pattern... ") )
	}

	### check if there is deletion in aligned_subject.
	del_aligned_subject <-0
	check_del <- stringr::str_count(aligned_subject, '0')
	if(check_del >1){
		    # print(paste0( (check_del-1), " del found in aligned_subject... ") )
		    aligned_subject <- gsub("0.*0", "0", aligned_subject)
	}


	left <- nchar( gsub("^(.*)0(.*)" , "\\1", aligned_subject) )
	right <- nchar( gsub("^(.*)0(.*)" , "\\2", aligned_subject) )
	
	left_seq <- substr(aligned_pattern, 1, left) 
	right_seq <- substr(aligned_pattern, (nchar(aligned_pattern)-right +1) , nchar(aligned_pattern) )

	if (left!=nchar(left_seq) | right!=nchar(right_seq) ) {
		print("A length of sequence doesn't match...")
	}

	new_allele_seq <- paste0(left_seq, 
				 substr(toString(aligned[2]), (left+1), ( nchar(aligned_pattern)-right  )) ,
				 right_seq)

	del_aligned_subject <-  stringr::str_count(new_allele_seq, '-')
	new_allele_seq <- gsub("-", "", new_allele_seq)
	

	### to check if it's correct
	to_adj_gtf_cor <- ""
	if(strand=="+"){
		to_adj_gtf_cor <- substr(new_allele_seq, (left + e) , (left +e +15) ) == e_seq
	}
	if(strand=="-"){
		to_adj_gtf_cor <- substr(new_allele_seq, (nchar(new_allele_seq) - right - e -15 +1), (nchar(new_allele_seq) -right -e +1) ) == 
			toString(reverseComplement(DNAString(e_seq)))
	}


	#new_allele_seq2 <- reverseComplement(DNAString(new_allele_seq))
	#substr(new_allele_seq2, (left +e), (left + e+ 15) )



	### save new allele sequence
	# write.fasta(sequences = new_allele_seq, names = allele_name, 
	# 	    file.out = paste0(output_directory, "/", prefix, ".", allele_name, ".fa"))

	print(paste0("Fasta of ", allele_name, " generated...") )
	
	#globalAlign<- pairwiseAlignment(seq1string, seq2string, type="local-global")
	#globalAlign



	#########################################
	###
	### 2. adjust coordinates
	###
	##########################################

		# HLA-DRB1, HLA-DPA1
	#if(gene=="HLA-DRB1"){
	if(to_adj_gtf_cor!=TRUE){
		stop(paste0(allele_name, ": sequence alignemnt didn't work well... check it") )

		tmp <- strsplit(new_allele_seq, e_seq )
		len <- nchar(tmp[[1]][1])
		
		substr(new_allele_seq, (len+1) , (len +10) )
		gtf_sub$V4 <- gtf_sub$V4 + len+1 - e
		gtf_sub$V5 <- gtf_sub$V5 + len+1 - e
		new_gtf_sub <- gtf_sub

		## need to finx
		new_allele_seq2 <- reverseComplement(DNAString(new_allele_seq))
		new_gtf_sub <- data.frame(new_gtf_sub[,1:3],
					  V4= ( nchar(new_allele_seq2) - new_gtf_sub$V5 +1),
					  V5= ( nchar(new_allele_seq2) - new_gtf_sub$V4 +1),
					   new_gtf_sub[,6:9])
		new_gtf_sub$V7="-"
		new_gtf_sub <- new_gtf_sub[order(new_gtf_sub$V4),]

		write.fasta(sequences = new_allele_seq2, names = allele_name,
                    file.out = paste0(output_directory, allele_name, ".fa"))

	} else{

		# check a length sequence
		#gtf_sub <- subset(gtf, V1==allele_name)

		if( nchar(new_allele_seq) !=  max(gtf_sub$V5) + nchar(left_seq) + nchar(right_seq)  ) {
		  print( paste0("new_allele_seq:", nchar(new_allele_seq)) )
		  print( paste0("gtf seq:", max(gtf_sub$V5) + nchar(left_seq) + nchar(right_seq)) )
		  
		  stop("Check new allele sequence, deletion was not corrected...")
		  
		}

		if(strand=="-"){
			# convert coordinates + to - : all IMGT alleles are + strand
			#gtf_sub$new_start = max(gtf_sub$V5) - gtf_sub$V4 +1
			#gtf_sub$new_end = max(gtf_sub$V5) - gtf_sub$V5 +1

			gtf_sub <- data.frame(gtf_sub[,1:3],
					      V4= ( max(gtf_sub$V5) - gtf_sub$V5 +1),
					      V5= ( max(gtf_sub$V5) - gtf_sub$V4 +1),
				 		gtf_sub[,6:9])
			gtf_sub$V7="-" 
		}

		# change coordinates by adding nchar(left_seq)
		gtf_sub$V4 <- gtf_sub$V4 + nchar(left_seq)
		gtf_sub$V5 <- gtf_sub$V5 + nchar(left_seq)
		nchar(new_allele_seq)


		# get lenths of utr regions from main reference  annotation (chr6)
		# add UTR annotation 
		# change utr and transcript length annotation in the allele annotation 

		anno_utr <- subset(anno, type=="UTR")		# annotation in chr6
		seq1=toupper(c2s(seq1[[1]]))
		seq2=toupper(c2s(seq2[[1]]))
		aligned_pattern <- toString(aligned[1] )
		aligned_subject <- toString(aligned[2] )
		gtf_sub_utr <- subset(gtf_sub, V3=="UTR")	# annotation in imgt


	
		add_utr_anno <- c()

		for(i in 1:nrow( anno_utr ) ){
			
		  # print(i)
		  
		  # get utr from chr6
			#ref_utr <- substr(seq1, anno_utr$new_start[i] , anno_utr$new_end[i] )
			ref_utr_in_chr6 <- substr(seq_chr6, anno_utr$start[i] , anno_utr$end[i] )

      # get utr seq from chr6
			#ref_utr_left <- substr(seq1, anno_utr$new_start[i] , (anno_utr$new_start[i] +30) )
			ref_utr_left <- substr(seq_chr6, anno_utr$start[i] , (anno_utr$start[i] +30) )
		
			# split chr6 by ref_utr_left
			tmp <- strsplit(aligned_pattern, ref_utr_left)
			len <- nchar(tmp[[1]][1])
			
			# get aligned utr in allele seq
			get_seq_from_aligned_subject <- substr(aligned_subject, (len+1), (len+31))
			#get_seq_from_aligned_pattern <- substr(aligned_pattern, (len+1), (len+31))
			

			uniq_string <- paste0(unique(strsplit(get_seq_from_aligned_subject, "")[[1]] ) , collapse="")
			uniq_string.len <- length(unique( strsplit(get_seq_from_aligned_subject, "")[[1]]) )
			if( length(uniq_string[grep("-", uniq_string)] ) ==1 & uniq_string.len ==1){
			  get_seq_from_aligned_subject <- strsplit(new_allele_seq, ref_utr_left)
			  len_subject <- nchar(get_seq_from_aligned_subject[[1]][1])
			}
			
			if( length(uniq_string[grep("-", uniq_string)] ) ==1 & uniq_string.len >1){
				#get_seq_from_aligned_subject <- strsplit(new_allele_seq, ref_utr_left)
			  get_seq_from_aligned_subject <- gsub("-", "", get_seq_from_aligned_subject)
				get_seq_from_aligned_subject <- strsplit(new_allele_seq, get_seq_from_aligned_subject)
				len_subject <- nchar(get_seq_from_aligned_subject[[1]][1])
			}

			
			#if( uniq_string != "-" & get_seq_from_aligned_subject!=""){
			if( length(uniq_string[grep("-", uniq_string)] ) ==0 ){
				get_seq_from_aligned_subject <- strsplit(new_allele_seq, get_seq_from_aligned_subject)
				len_subject <- nchar(get_seq_from_aligned_subject[[1]][1])
			}

			add_utr_anno <- rbind(add_utr_anno,
					      data.frame(gtf_sub[1,1:2],
							 V3="UTR",
							 V4=len_subject+1,
							 V5=len_subject+nchar(ref_utr_in_chr6),
							 gtf_sub[1,6:8],
							 V9=gtf_sub_utr$V9[1]  ) )

		}
		new_gtf_sub <- rbind(subset(gtf_sub, V3!="UTR"), add_utr_anno)
		new_gtf_sub <- new_gtf_sub[order(new_gtf_sub$V4),]
		gene_start <- min(new_gtf_sub[new_gtf_sub$V3!="gene" & new_gtf_sub$V3!="transcript",]$V4)
		gene_end <- max(new_gtf_sub[new_gtf_sub$V3!="gene" & new_gtf_sub$V3!="transcript",]$V5)

		new_gtf_sub[new_gtf_sub$V3=="gene" | new_gtf_sub$V3=="transcript",]$V4 <- gene_start
		new_gtf_sub[new_gtf_sub$V3=="gene" | new_gtf_sub$V3=="transcript",]$V5 <- gene_end
	}

		add_exon <- subset(new_gtf_sub, V3=="UTR")
		add_exon$V3 <- "exon"
		new_gtf_sub <- rbind(new_gtf_sub, add_exon)
		new_gtf_sub <- new_gtf_sub[order(new_gtf_sub$V4),]


		comb_exon <- c()
		exons <- subset(new_gtf_sub, V3=="exon")
		end_pos=10000
		b="no"
		ex <- exons[1,]
		for(e in 2:nrow(exons)){
			#print(e)
			#print(exons$V4[e])
			#print(ex$V5)
			if(exons$V4[e] < ex$V5 | (exons$V4[e] -1)==ex$V5 ){
				comb_exon <- rbind(comb_exon,
						   data.frame(ex[,1:4], exons[e,5:9]) ) 
				b="ok"
			}
			else{
				if(b!="ok"){
					comb_exon <- rbind(comb_exon, ex)
				}
				b="no"
			}
			ex <- exons[e,]
		}
		if(b=="no"){
			comb_exon <- rbind(comb_exon, ex)
		}


		new_gtf_sub <- rbind(subset(new_gtf_sub, V3!="exon"),
				     comb_exon)
		
		v9_na <- new_gtf_sub[is.na(new_gtf_sub$V9),]
		if(nrow(v9_na) > 0){
			v9_na_to_replace <- new_gtf_sub[new_gtf_sub$V3=="transcript",]$V9
			v9_na$V9 <- v9_na_to_replace
			new_gtf_sub <- rbind(new_gtf_sub[!is.na(new_gtf_sub$V9),],
					     v9_na)
		}
		new_gtf_sub <- new_gtf_sub[order(new_gtf_sub$V4),]


	# write.table(new_gtf_sub, 
	# 	    paste0(output_directory, "/", prefix, ".", allele_name, ".gtf"), 
	# 	    sep="\t", quote=F, row.names=F, col.names=F)

	print(paste0("GTF of ", allele_name, " generated...") )
	
	new_allele_seq2 <- DNAStringSet(new_allele_seq)
	names(new_allele_seq2) <- allele_name
	output <- list(new_allele_seq2, new_gtf_sub)
	return(output)
}



