#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
parser <- ArgumentParser(description='convert count matrix to TPM matrix (gene|length format as input ids)')
parser$add_argument('-i', '--matrix', type='character', required=TRUE,
                    help='input count matrix. Rows are genes. Columns are samples.')
parser$add_argument('--pseudo-count', type='double', default=1.0,
                    help='pseudo-count added to log2 transform in ttest')
parser$add_argument('-o', '--output-file', type='character', required=TRUE,
                    help='output file')
parser$add_argument('-r', '--regiontype', type='character', required=F,
                    help='region-type:gene,promoter...')
args <- parser$parse_args()


message('read count matrix: ', args$matrix)
mat <- read.table(args$matrix, header = TRUE, row.names=1, check.names=FALSE, sep='\t')
region_type <- args$regiontype
matrix_tpm <- as.matrix(mat)
gene.len <- as.numeric(lapply(strsplit(rownames(matrix_tpm),"\\|"), function(x) x[2]))
if (is.na(gene.len[2])){
	message("length not found in gene|length, using default length file in /BioII/lulab_b/baopengfei/shared_reference/hg38/ ")
	gene.len <- read.delim(paste0("/BioII/lulab_b/baopengfei/shared_reference/hg38/",region_type,".length"),sep="\n",stringsAsFactors=F,header=F)[,1]
}
matrix_tpm <- 1000*matrix_tpm / gene.len
matrix_tpm <- t(t(matrix_tpm) * 1e6 / colSums(matrix_tpm))

matrix_tpm <- as.data.frame(matrix_tpm)
matrix_tpm$`gene_id` <- rownames(matrix_tpm)
matrix_tpm <- matrix_tpm[,c(ncol(matrix_tpm),1:(ncol(matrix_tpm)-1))]

message('Write results to output file: ', args$output_file)
write.table(x=matrix_tpm, file=args$output_file, sep='\t', quote=FALSE, row.names=F, col.names=T)

