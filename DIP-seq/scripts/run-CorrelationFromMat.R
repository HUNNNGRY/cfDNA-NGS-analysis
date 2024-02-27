# run cor from 2 matrix
# last 210803 by pengfei (newer version: run-MatCor.R)
options(max.print = 3)

#suppressPackageStartupMessages(library())
parser <- argparse::ArgumentParser(description='cal cor record with cor, p and fdr from 2 mat')
parser$add_argument('-a', '--matrix1', type='character', required=TRUE,
                    help='first input matrix, tsv, usually Rows are sample records and Columns are genes/genesets, or you need to set --record flag as row')
parser$add_argument('-b', '--matrix2', type='character', required=TRUE,
                    help='second input matrix, the same requirement as first input matrix')
parser$add_argument('-r', '--record', type='character', default="col",
                    help='which is record: row or col, (default: col)')
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
args <- parser$parse_args()

a.in <- args$matrix1
b.in <- args$matrix2
a.tag <- args$matrix1tag
b.tag <- args$matrix2tag
record <- args$record
adjp <- args$padj
type <- args$method
outfile <- args$outfile
diag <- args$diag

# test
# a.in <- "/BioII/lulab_b/baopengfei/projects/multi-omics-explore/data/pico-matrix/gsva_matrix_ENSG_noMT.txt"
# b.in <- "/BioII/lulab_b/baopengfei/projects/multi-omics-explore/data/pico-matrix/ssgsea_matrix_ENSG_noMT.txt"
# a.tag <- "tumor"
# b.tag <- "normal"
# adjp <- 'BH'
# type <- "spearman"
# record <- "col"
# outfile <- "./test"
# diag <- True


## read input mat
a <- read.table(a.in,sep = "\t", header = T, row.names = 1,check.names = F,stringsAsFactors = F)
b <- read.table(b.in,sep = "\t", header = T, row.names = 1,check.names = F,stringsAsFactors = F)

if(record!="row"){
  message("input mat cols are records, flipping mat...")
  a <- t(a)
  b <- t(b)
} else{
  message("input mat rows are records, will not change...")
}

## change to patient name, and only keep the same row records/patients
rownames(a) <- sub("-T$|-N$|-PBMC$|-pico$","",rownames(a))
rownames(b) <- sub("-T$|-N$|-PBMC$|-pico$","",rownames(b))
a <- a[rownames(a) %in% rownames(b),]
b <- b[rownames(b) %in% rownames(a),]

b <- b[match(rownames(a),rownames(a)),]
#all(colnames(a)==colnames(b))
#all(rownames(a)==rownames(b))

## add colnames to diff tags as surfix, to prevent the same feature name
colnames(a) <- paste0(colnames(a),".",a.tag)
colnames(b) <- paste0(colnames(b),".",b.tag)

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

#dim(a)
cor_3 <- Hmisc::rcorr(as.matrix(cbind(a,b)),type=type)

my_cor_matrix <- flat_cor_mat(cor_3$r, cor_3$P)
my_cor_matrix$FDR.BH <- p.adjust(my_cor_matrix$p, method=adjp)

colnames(my_cor_matrix)[1:2] <- c("from","to")

## remove dup rows
cols = c(1,2)
newdf = my_cor_matrix[,cols]
for (i in 1:nrow(my_cor_matrix)){
  newdf[i, ] = sort(my_cor_matrix[i,cols])
}
# suppressPackageStartupMessages(library(foreach))
# suppressPackageStartupMessages(library(doParallel))
# registerDoParallel(10)  # use multicore, set to the number of our cores
# newdf <- foreach(i=1:nrow(my_cor_matrix),.verbose=F,.combine=rbind) %dopar% {
#   sort(my_cor_matrix[i,cols])
# }
# write.table(newdf,paste0(outfile,".test"),sep = "\t",quote = F,row.names = T,col.names = T)

#newdf[,1] <- apply(my_cor_matrix,1, function(x) sort(x[cols])[1])
#newdf[,2] <- apply(my_cor_matrix,1, function(x) sort(x[cols])[2])
my_cor_matrix <- my_cor_matrix[!duplicated(newdf),]
my_cor_matrix <- my_cor_matrix[my_cor_matrix$from!=my_cor_matrix$to,]  # rm self cmp

write.table(my_cor_matrix,outfile,sep = "\t",quote = F,row.names = F,col.names = T)
if(diag==T){
  write.table(cor_3$r,paste0(outfile,".r"),sep = "\t",quote = F,row.names = T,col.names = T)
  write.table(cor_3$P,paste0(outfile,".p"),sep = "\t",quote = F,row.names = T,col.names = T)
}
message("done!")

# require(ggpubr)
# require(tidyverse)
# require(Hmisc)
# require(corrplot)
# 
# head(mtcars)
# library(ggpubr)
# my_data <- mtcars
# my_data$cyl <- factor(my_data$cyl)
# str(my_data)
# ggscatter(my_data, x = "wt", y = "mpg",
#           add = "reg.line", conf.int = TRUE,
#           cor.coef = TRUE, cor.method = "pearson",
#           xlab = "Weight (1000 lbs)", ylab = "Miles/ (US) gallon")
# 
# # Shapiro-Wilk normality test for mpg
# shapiro.test(my_data$mpg) # => p = 0.1229
# 
# # Check for the normality of "mpg""
# ggpubr::ggqqplot(my_data$mpg, ylab = "MPG")
# 
# # diff cor.test methods
# res <- cor.test(my_data$wt, my_data$mpg, method = "pearson")
# res
# str(res)
# 
# res2 <- cor.test(my_data$mpg, my_data$wt, method = "kendall")
# res2
# 
# res3 <- cor.test(my_data$wt, my_data$mpg, method = "spearman")
# res3
# 
# 
# # cor
# my_data <- dplyr::select(mtcars, mpg, disp, hp, drat, wt, qsec)
# head(my_data)
# cor_1 <- round(cor(my_data), 2)
# cor_1
# 
# cor_2 <- Hmisc::rcorr(as.matrix(my_data))
# cor_2
# 
# 
# flat_cor_mat <- function(cor_r, cor_p){
#   #This function provides a simple formatting of a correlation matrix
#   #into a table with 4 columns containing :
#   # Column 1 : row names (variable 1 for the correlation test)
#   # Column 2 : column names (variable 2 for the correlation test)
#   # Column 3 : the correlation coefficients
#   # Column 4 : the p-values of the correlations
#   library(tidyr)
#   library(tibble)
#   cor_r <- rownames_to_column(as.data.frame(cor_r), var = "row")
#   cor_r <- gather(cor_r, column, cor, -1)
#   cor_p <- rownames_to_column(as.data.frame(cor_p), var = "row")
#   cor_p <- gather(cor_p, column, p, -1)
#   cor_p_matrix <- left_join(cor_r, cor_p, by = c("row", "column"))
#   cor_p_matrix
# }
