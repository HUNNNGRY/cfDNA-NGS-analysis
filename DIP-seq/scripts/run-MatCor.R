# run cor from 2 matrix
# last 2101110 by pengfei (not finished)
options(max.print = 3)

#suppressPackageStartupMessages(library())
parser <- argparse::ArgumentParser(description='calculate combination records with correlation, p and FDR from 2 mat')
parser$add_argument('-a', '--matrix1', type='character', required=TRUE,
                    help='first input matrix, tsv, usually Columns are samples/attribute/variable, default only calculates correlation among columns, or you need to set --record col')
parser$add_argument('-b', '--matrix2', type='character', required=FALSE,
                    help='second input matrix, the same requirement as first input matrix, can be skipped')
parser$add_argument('-r', '--record', type='character', default="row",
                    help='which dim should be considered as records (not attribute/variable): row or col, (default: row)')
parser$add_argument('-atag', '--matrix1tag', type='character', default="mat1",
                    help='char surfix for first input matrix, will be appended to feature name to avoid the same feature name between mat')
parser$add_argument('-btag', '--matrix2tag', type='character', default="mat2",
                    help='char surfix for second input matrix, the same as first input matrix')
parser$add_argument('-q', '--padj', type='character', default="BH",
                    help='adj P method: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none", (default=BH)')
parser$add_argument('-m', '--method', type='character', default="spearman",
                    help='cor method: "pearson","spearman", (default=spearman)')
parser$add_argument('-o', '--outfile', type='character', default="./cor",
                    help='outfile prefix (default=./cor)')
parser$add_argument('-d', '--diag', type='logical', default=TRUE,
                    help='whether to additionally output diag matrix r and p (default=TRUE)')
parser$add_argument('--omitNA', type='logical', default=TRUE,
                    help='whether to omit NA records before cor (default=TRUE)')
args <- parser$parse_args()


a.in <- args$matrix1
a.tag <- args$matrix1tag

matrix2 <- try(args$matrix2)
#b.exist <- exists("matrix2")
b.exist <- TRUE
if(b.exist){
  b.in <- args$matrix2
  b.tag <- args$matrix2tag
}


record <- args$record
adjp <- args$padj
type <- args$method
outfile <- args$outfile
diag <- args$diag
omitNA <- args$omitNA

message("b.exist: ",b.exist)
message("a.in: ",a.in)
message("a.tag: ",a.tag)
#message("b.in: ",b.in)
#message("b.tag: ",b.tag)
message("record: ",record)
message("adjp: ",adjp)
message("type: ",type)
message("outfile: ",outfile)
message("diag: ",diag)
message("omitNA: ",omitNA)


# test
# /usr/bin/Rscript scripts/run-MatCor.R \
# -a output/cnv.txt \
# -atag cnv \
# -b output/exp.txt \
# -btag exp \
# -r col \
# -m spearman \
# -o output/cor/cnv.exp.gene.spearman.omitNA

# setwd("/BioII/lulab_b/baopengfei/projects/multi-omics-explore/")
# a.in <- "output/cnv.txt"
# b.in <- "output/met.txt"
# a.tag <- "cnv"
# b.tag <- "met"
# adjp <- 'BH'
# type <- "spearman"
# record <- "col"
# outfile <- "./"
# diag <- TRUE
# b.exist <- TRUE
# omitNA <- TRUE

## read input mat
a <- read.table(a.in,sep = "\t", header = T, row.names = 1,check.names = F,stringsAsFactors = F)
if(b.exist){
  b <- read.table(b.in,sep = "\t", header = T, row.names = 1,check.names = F,stringsAsFactors = F)
}


if(record!="row"){
  message("treat input mat cols as records, cols as attributes/variables, flipping mat...")
  a <- t(a)
  if(b.exist){b <- t(b)}
} else{
  message("treat input mat rows as records, cols as attributes/variables, will not transpose...")
}

## omit before t() convert to mat
message("a.rows: ",nrow(a), "; a.cols: ",ncol(a))
if(omitNA){
  a <- na.omit(a)
  message("omit NA rows, remaining: ")
  message("a.rows: ",nrow(a), "; a.cols: ",ncol(a))
} 
if(b.exist){
  message("b.rows: ",nrow(b), "; b.cols: ",ncol(b))
  if(omitNA){
    b <- na.omit(b)
    message("omit NA rows, remaining: ")
    message("b.rows: ",nrow(b), "; b.cols: ",ncol(b))
  } 
}

## change sample id (in-house matrix only)
#change to patient name, and only keep the same row records/patients
rownames(a) <- sub("-T$|-N$|-PBMC$|-pico$","",rownames(a))
colnames(a) <- sub("-T$|-N$|-PBMC$|-pico$","",colnames(a))

## add colnames to diff tags as surfix, to prevent the same feature name
colnames(a) <- paste0(colnames(a),".",a.tag)

if(b.exist){
  a <- a[rownames(a) %in% rownames(b),]
  #all(colnames(a)==colnames(b))
  #all(rownames(a)==rownames(b))

  rownames(b) <- sub("-T$|-N$|-PBMC$|-pico$","",rownames(b))
  colnames(b) <- sub("-T$|-N$|-PBMC$|-pico$","",colnames(b))
  b <- b[rownames(b) %in% rownames(a),]
  b <- b[match(rownames(a),rownames(a)),]
  colnames(b) <- paste0(colnames(b),".",b.tag)
}




## run cor
flat_cor_mat <- function(cor_r, cor_p){
  #This function provides a simple formatting of a correlation matrix
  #into a table with 4 columns containing :
  # Column 1 : row names (variable 1 for the correlation test)
  # Column 2 : column names (variable 2 for the correlation test)
  # Column 3 : the correlation coefficients
  # Column 4 : the p-values of the correlations
  library(tidyr)
  library(tibble)
  cor_r <- rownames_to_column(as.data.frame(cor_r), var = "row")
  cor_r <- gather(cor_r, column, cor, -1)
  cor_p <- rownames_to_column(as.data.frame(cor_p), var = "row")
  cor_p <- gather(cor_p, column, p, -1)
  cor_p_matrix <- dplyr::left_join(cor_r, cor_p, by = c("row", "column"))
  cor_p_matrix
}

# if too many cols, will meet error, 
# if more than 5000 cols left, will random select 5000 cols each mat
if(b.exist){
  if(ncol(a)>=10000 | ncol(b)>=10000){
    message("more than 10000 cols, will random select 5000 cols")
    cl <- sample(1:min(ncol(a),ncol(b)),size = min(ncol(a),ncol(b),5000),replace = F)
    a <- a[,cl]
    b <- b[,cl]
  }
}else{
  if(ncol(a)>=10000){
    message("more than 10000 cols, will random select 5000 cols")
    cl <- sample(1:ncol(a),size = min(ncol(a),5000),replace = F)
    a <- a[,cl]
  }
}

#dim(a)
if(b.exist){
  tmp <- as.matrix(cbind(a,b))
} else {
  tmp <- as.matrix(a)
}


cor_3 <- Hmisc::rcorr(tmp,type=type)

my_cor_matrix <- flat_cor_mat(cor_3$r, cor_3$P)
my_cor_matrix$FDR.BH <- p.adjust(my_cor_matrix$p, method=adjp)

colnames(my_cor_matrix)[1:2] <- c("from","to")

## remove dup rows
cols = c(1,2)
newdf = my_cor_matrix[,cols]
for (i in 1:nrow(my_cor_matrix)){
  newdf[i, ] = sort(my_cor_matrix[i,cols])
}

my_cor_matrix <- my_cor_matrix[!duplicated(newdf),]
my_cor_matrix <- my_cor_matrix[my_cor_matrix$from!=my_cor_matrix$to,]  # rm self cmp

write.table(my_cor_matrix,outfile,sep = "\t",quote = F,row.names = F,col.names = T)
if(diag==T){
  write.table(cor_3$r,paste0(outfile,".r"),sep = "\t",quote = F,row.names = T,col.names = T)
  write.table(cor_3$n,paste0(outfile,".n"),sep = "\t",quote = F,row.names = T,col.names = T)
  write.table(cor_3$P,paste0(outfile,".p"),sep = "\t",quote = F,row.names = T,col.names = T)
}
message("done!")
