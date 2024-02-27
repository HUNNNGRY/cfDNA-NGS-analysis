# split matrix to train and test sets
# last 20210519 by pengfei

suppressPackageStartupMessages(library(argparse))
parser <- ArgumentParser(description='split matrix to train and test sets')
parser$add_argument('-i', '--matrix', type='character', required=TRUE,
                    help='input feature matrix tsv. usually Rows are features(genes), cols are samples(CRC-*-*-*)')
parser$add_argument('-l', '--label', type='character', required=TRUE,
                    help='input label/tag txt. number equal to samples')
parser$add_argument('-a', '--addlabel', type='logical', default=FALSE,
                    help='whether to add label/tag to out files, default=FALSE')
parser$add_argument('-r', '--trainRatio', type='double', default=0.7,
                    help='trainRatio (default=0.7)')
#parser$add_argument('-s', '--sampleInRow', type='logical', default=FALSE,
#                    help='whether sampleInRow (default=FALSE, if TRUE, then mat will be transposed)')
parser$add_argument('-o', '--outdir', type='character', default="./",
                    help='outdir, result will be sampleInRow, comma sep ')
args <- parser$parse_args()
dir.create(args$outdir,showWarnings = F)


# test
# setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq")
# matrix <- read.table("/BioII/lulab_b/baopengfei/projects/methylation/output/lulab-diffbind/14CRCvs19NC-rpkm-pTop100-pairedWGS.txt",sep = "\t", header = T, row.names = 1, stringsAsFactors = F, check.names = F)
# label <- read.table("/BioII/lulab_b/baopengfei/projects/methylation/output/lulab-diffbind/14CRCvs19NC-label.txt",sep = "\t", header = F, stringsAsFactors = F, check.names = F)$V1
# addlabel <- TRUE
# trainRatio <- 0.6
# outdir <- "."

matrix <- read.table(args$matrix,sep = "\t", header = T, row.names = 1, stringsAsFactors = F, check.names = F)
label <- read.table(args$label,sep = "\t", header = F, stringsAsFactors = F ,check.names = F)$V1
addlabel <- args$addlabel
trainRatio <- args$trainRatio
outdir <- args$outdir

if(ncol(matrix)==length(label))
{
  matrix <- as.data.frame(t(matrix))   # df to matrix
  message("samples in column")
} else if (nrow(matrix)==length(label)) {
  matrix <- matrix
  message("samples in row")
} else {
  message("not equal number between label/tag and matrix")
}

if(addlabel){
  matrix$group <- label
  message("add label/tag to column")
}
message(paste0('rows: ',nrow(matrix),'   cols: ',ncol(matrix)))

indices <- sample(1:nrow(matrix))
indices <- indices[1:floor(nrow(matrix) * trainRatio)]

write.table(matrix[indices, ],paste0(outdir,"/train-mat.txt"),quote = F,sep = ",",col.names = T,row.names = F)
message("write trainset done")
write.table(matrix[-indices, ],paste0(outdir,"/test-mat.txt"),quote = F,sep = ",",col.names = T,row.names = F)
message("write testset done")