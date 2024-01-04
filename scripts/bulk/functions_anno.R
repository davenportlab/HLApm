

splitChr <- function(fa){
	#fa="./GRCh38.p13.genome.fa.gz"

	require(Biostrings)
	scaffolds = readAAStringSet(filepath = fa, use.names = T)
	t2 = split(scaffolds, names(scaffolds))
	sapply(names(t2), function (x) writeXStringSet(t2[[x]], filepath = paste0(gsub("(.*)\\s.*", "\\1", x),".fasta")))
}


annoHLA <- function(gtf_patch_hapl_scaff){
	# gtf_patch_hapl_scaff =  "/lustre/scratch119/realdata/mdt3/teams/davenport/ws/ref/GRCh38.99/gencode.v33.chr_patch_hapl_scaff.annotation.gtf.gz"

	require(rtracklayer)
	require(reshape2)

	gtf <- rtracklayer::import(gtf_patch_hapl_scaff)
	study_rowdata <- as.data.frame(gtf)
	study_rowdata$gene_id <- gsub("\\.\\d+", "", study_rowdata$gene_id)
	
	anno_hla <- study_rowdata[grep("HLA-", study_rowdata$gene_name), ]
	anno_hla <- anno_hla[anno_hla$type=="gene", ]

	anno_hla.wide <- dcast(anno_hla, gene_name ~ seqnames, value.var="gene_id")
	anno_hla.wide.chr6_na <- anno_hla.wide[is.na(anno_hla.wide$chr6), ]

	#tmp <- data.frame(table(anno_hla$seqnames))
	#tmp <- tmp[tmp$Freq >0,]

	return(anno_hla.wide)
}



make_star_index.hap <- function(fa, gtf_patch_hapl_scaff, hap_name, genomeDir){
	# gtf_patch_hapl_scaff =  "/lustre/scratch119/realdata/mdt3/teams/davenport/ws/ref/GRCh38.99/gencode.v33.chr_patch_hapl_scaff.annotation.gtf.gz"
	# fa="/lustre/scratch119/realdata/mdt3/teams/davenport/ws/ref/GRCh38.99/by_chrs/GL000255.2.fasta"
	# genomeDir="/lustre/scratch119/realdata/mdt3/teams/davenport/ws/ref/GRCh38.99/by_chrs"
	# hap_name= "GL000255"

	split_gtf <- paste0("zless ", gtf_patch_hapl_scaff, " |  grep ", hap_name, " > ", hap_name, ".gtf")
	print("splitting gtf")
	system(split_gtf)
	
	output_dir <- file.path(genomeDir, paste0("star_index.", hap_name))
	if (!dir.exists(output_dir)){
		dir.create(output_dir)
	}


	run_index <- paste0("/software/team282/bin/STAR --runThreadN 6 ",
			    "--runMode genomeGenerate ",
			    "--genomeDir ", output_dir, " ",
			    "--genomeFastaFiles ", fa, " ",
			    "--sjdbGTFfile ", hap_name, ".gtf ",
			    "--genomeSAindexNbases 10 ",
			    "--sjdbOverhang 99")
	write.table(run_index, paste0("tmp.index.", hap_name, ".sh"), 
		    quote=F, row.names=F, col.names=F )

	# bsub -q normal -G team282 -o ./logs/log.genomeIndex.%J -e ./logs/log.genomeIndex.error.%J -n6 -R "span[hosts=1]" -R "select[mem>42000] rusage[mem=42000]" -M42000 sh tmp.index.GL000255.sh 
}

# make_star_index.hap(fa, gtf_patch_hapl_scaff, hap_name="GL000256", genomeDir)



