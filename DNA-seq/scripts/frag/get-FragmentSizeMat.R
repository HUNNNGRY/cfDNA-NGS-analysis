# merge multiple fragment size depth ratio to matrix
# last 220721 by pengfei

options(stringsAsFactors = F)

suppressPackageStartupMessages(library(argparse))
parser <- ArgumentParser(description='GO/KEGG/DO geneset enrichment by ORA/GSEA using diff table (network independent version)')
parser$add_argument('-s', '--small', type='character', required=TRUE,
                    help='small (100-150bp) GC-content corrected depth or CPM matrix. Rows are genes (^ENSG*, *_ENSG*). Columns are samples')
parser$add_argument('-l', '--long', type='character', required=TRUE,
                    help='long (151-220bp) GC-content corrected depth or CPM matrix. Rows are genes (^ENSG*, *_ENSG*). Columns are samples')
parser$add_argument('-o', '--outfile', type='character', required=TRUE,
                    help='output FragSizeRatio matrix. Rows are genes (^ENSG*, *_ENSG*). Columns are samples')
args <- parser$parse_args()
small_matrix <- args$small
long_matrix <- args$long
outfile <- args$outfile


s.ratio <- read.table(small_matrix,header = T,row.names = 1, check.names = F) #.not tsv
l.ratio <- read.table(long_matrix,header = T,row.names = 1, check.names = F) #.not tsv

#s.ratio <- log2(as.matrix(s.ratio)+0.01)
#l.ratio <- log2(as.matrix(l.ratio)+1)
l.ratio <- l.ratio[,gsub("-short","-long",colnames(s.ratio))]
all(gsub("-short|-long","",colnames(s.ratio))==gsub("-short|-long","",colnames(l.ratio)))  # should be true

ratio <- s.ratio/(l.ratio+1)
colnames(ratio) <- gsub("-short|-long","",colnames(s.ratio))
ratio[1:3,1:3]


## write out mat
ratio.out <- cbind(rownames(ratio),ratio)
colnames(ratio.out)[1] <- "gene_id"
#write.table(ratio.out,paste0("./output/",dst,"/matrix/DNA-FragRatio_matrix_",region,".txt"),row.names = F,col.names = T,sep = "\t",quote = F)
message("write to ",outfile)
write.table(ratio.out,outfile,row.names = F,col.names = T,sep = "\t",quote = F)


