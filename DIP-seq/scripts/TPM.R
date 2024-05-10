#! /usr/bin/env Rscript

#TODO:
# add filtering of features and samples
# add support for all smps libSize input

suppressPackageStartupMessages(library(argparse))
parser <- ArgumentParser(description='convert count matrix to TPM matrix (gene|length format as input ids)')
parser$add_argument('-i', '--matrix', type='character', required=TRUE,
                    help='input count matrix. Rows are genes. Columns are samples.')
parser$add_argument('-s', '--sample', type='character', default="all",
                    help='tsv/txt file: pass filter colnames/samples, default: all (keep all colnames/samples), optional: txt path of valid samples, will change mat to the same order')
parser$add_argument('--pseudo-count', type='double', default=1.0,
                    help='pseudo-count added to log2 transform in ttest')
parser$add_argument('-o', '--output-file', type='character', required=TRUE,
                    help='output file')
# parser$add_argument('-r', '--regiontype', type='character', required=F,
#                     help='region-type:gene,promoter...')
parser$add_argument('-e', '--externalLibSize', type='character', default="FALSE",
                    help='to replace colSums, default: "FALSE", optional: txt path of libsize (numbers);tsv without header first col lists valid samples same order of matrix, second col list corresponding libsize')
parser$add_argument('--geneLength', type='character', default="FALSE",
                    help='default: FALSE, optional: txt path of geneLength (int), in same order of valid genes')
args <- parser$parse_args()
sa <- (args$sample)
gl <- (args$geneLength)


# # test
# args <- list()
# args$matrix <- "/data2/lulab1/bpf/projects/WCHSU-FTC/output/SLE/call_peak_dedup/count_matrix/cfpeakCNN_b5_d50_p1.txt"
# sa <- "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/SLE/240316/sample_ids.txt"
# gl <- 'FALSE'

message('read count matrix: ', args$matrix)
mat <- read.table(args$matrix, header = TRUE, row.names=1, check.names=FALSE, sep='\t')
#region_type <- args$regiontype
matrix_tpm <- as.matrix(mat)

if(sa=="all"){
  message(paste0("keep all samples from input mat"))
  s <- colnames(matrix_tpm)
} else {
  message(paste0("only keep samples from path: ",sa))
  s <- read.delim(sa,header = F,stringsAsFactors = F)$V1
}
#count.matrix <- count.matrix[,colnames(count.matrix) %in% s,,drop=F] # do not change col order
matrix_tpm <- matrix_tpm[,s,drop=F] # change col order


if(gl=='FALSE'){
	message("try to get feature length from gene format (2nd '|' col in feature)")
	gene.len <- as.numeric(lapply(strsplit(rownames(matrix_tpm),"\\|"), function(x) x[2]))
	if(is.na(gene.len[1])){
		message("try to get feature length from peak format (7th-6th '|' col in feature)")
		gene.len <- as.numeric(lapply(strsplit(rownames(matrix_tpm),"\\|"), function(x) x[7])) - as.numeric(lapply(strsplit(rownames(matrix_tpm),"\\|"), function(x) x[6]))
		if(is.na(gene.len[1])){
			 stop("wrong length format")
		}
	}
}else{
	message(paste0("gene length provided"))
	geneLen <- read.delim(gl,header = F,stringsAsFactors = F)$V1	
}
	
#if (is.na(gene.len[2])){
#	message("length not found in gene|length, using default length file in /BioII/lulab_b/baopengfei/shared_reference/hg38/ ")
#	gene.len <- read.delim(paste0("/BioII/lulab_b/baopengfei/shared_reference/hg38/",region_type,".length"),sep="\n",stringsAsFactors=F,header=F)[,1]
#}
message("gL:",gene.len[1])

matrix_tpm <- 1000*matrix_tpm / gene.len
if(args$externalLibSize=="FALSE"){
message("use rowSums as libSize")
matrix_tpm <- t(t(matrix_tpm) * 1e6 / colSums(matrix_tpm))
}else{
message("use provided file as libSize")
#libSize <- read.delim(ex,header = F,stringsAsFactors = F)$V1
#libSize <- read.table(args$externalLibSize,header = F,stringsAsFactors = F)$V2 #V1: sample
  libSize <- read.table(args$externalLibSize,header = F,stringsAsFactors = F)
  tmp <- libSize$V1
  libSize <- libSize$V2
  names(libSize) <- tmp
  matrix_tpm <- t(t(matrix_tpm) * 1e6 / libSize[s])
}
matrix_tpm <- as.data.frame(matrix_tpm)
matrix_tpm$`gene_id` <- rownames(matrix_tpm)
matrix_tpm <- matrix_tpm[,c(ncol(matrix_tpm),1:(ncol(matrix_tpm)-1))]

message('Write results to output file: ', args$output_file)
write.table(x=matrix_tpm, file=args$output_file, sep='\t', quote=FALSE, row.names=F, col.names=T)

