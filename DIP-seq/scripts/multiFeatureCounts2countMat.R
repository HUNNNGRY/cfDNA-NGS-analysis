# multiFeatureCounts2countMat.R
# raw multibam featurecounts result to count matrix
# last 20210512

#getwd()
#setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq")

suppressPackageStartupMessages(library(argparse))
parser <- ArgumentParser(description='get multi bam featurecounts results to count matrix')
parser$add_argument('-i', '--input', type='character', required=TRUE,help='Input multibam featurecounts result file path')
parser$add_argument('-o', '--output', type='character', required=TRUE,help='Output "gene_id|length CRC-PKU-10-me" format file path')
args <- parser$parse_args()

c <- read.table(args$input,stringsAsFactors = F,header = T,check.names = F)

colnames(c)[c(1,6)] <- c("gene_id","lenth")
c$gene_id <- paste(c$gene_id,c$lenth,sep = "|")
c <- c[,-c(2:6)]

colnames(c)[2:ncol(c)] <- basename(colnames(c)[2:ncol(c)])
colnames(c)[2:ncol(c)] <- as.character(lapply(strsplit(colnames(c)[2:ncol(c)],"\\."), function(x) x[1]))
#tail(c$gene_id)

write.table(c,args$output,row.names = F, col.names = T, quote = F,sep = "\t")
message("finished.")
