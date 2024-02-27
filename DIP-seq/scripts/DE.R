#! /usr/bin/env Rscript
#last update at 20211125 by pengfei
#211211: change norm-method choice; rm cpm(method)

suppressPackageStartupMessages(library(argparse))
parser <- ArgumentParser(description='Differential expression')
parser$add_argument('-i', '--matrix', type='character', required=TRUE,
                    help='input matrix. Rows are genes. Columns are samples.')
parser$add_argument('-d', '--datatype', type='character', default="count",
                    choices=c("count","normalized"),
                    help='input matrix type. count (raw counts) or normalized (cpm/tpm). for normalized mat, only wilcox and t-test supported. default:count')
parser$add_argument('-p', '--positive-ids', type='character', required=TRUE,
                    help='file contains ids of positive class')
parser$add_argument('-n', '--negative-ids', type='character', required=TRUE,
                    help='file contains ids of negative class')
parser$add_argument('-m', '--method', type='character', default="edger_glmlrt",
                    choices=c('deseq2_wald','deseq2_lrt','deseq2_sc', 'edger_exact', 'edger_glmqlf', 'edger_glmlrt', 'wilcox', 'limma', 'ttest'),
                    help='differential expression method to use')
#parser$add_argument('--norm-method', type='character', default='TMM',
#                    choices=c('RLE', 'CPM', 'TMM', 'upperquartile'))
parser$add_argument('--norm-method', type='character', default='TMM',
                    choices=c("TMM","TMMwsp","RLE","upperquartile","none"))
parser$add_argument('--pseudo', type='double', default=1.0,
                    help='pseudo-count added to log2 transform in ttest')
parser$add_argument('--sort', type='logical', default=FALSE,
                    help='whether to sort diff table according to p and FC, default: FALSE')
parser$add_argument('--cores', type='double', default=1,
                    help='multi-cores using for deseq2, default:1')
parser$add_argument('-o', '--output-file', type='character', required=TRUE,
                    help='output file')
args <- parser$parse_args()

# test
# setwd("/BioII/lulab_b/baopengfei/cooperation/yinjianhua")
# mat <-  read.table("output/rna-matrix/rna-matrix/count_matrix_gene.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
# positive_samples <- read.delim("meta/wes_rt_ids.txt",sep="\n",stringsAsFactors=F,header=F)[,1]
# negative_samples <- read.delim("meta/wes_pt_ids.txt",sep="\n",stringsAsFactors=F,header=F)[,1]
# samples <- c(positive_samples,negative_samples)
# mat <- as.matrix(mat[,samples])
# colnames(mat)
# table(samples %in% colnames(mat))

pseudo <- args$pseudo
#message('read matrix: ', args$matrix)
mat <- read.table(args$matrix, header = TRUE, row.names=1, check.names=FALSE, sep='\t')
#datatype <- args$datatype

#message('Read positive sample ids')
positive_samples <- read.delim(args$positive_ids,sep="\n",stringsAsFactors=F,header=F)[,1]
#message('Read negative sample ids')
negative_samples <- read.delim(args$negative_ids,sep="\n",stringsAsFactors=F,header=F)[,1]
message('Number of positive samples: ', length(positive_samples))
message('Number of negative samples: ', length(negative_samples))

samples <- c(positive_samples, negative_samples)
group <- c(rep("positive",length(positive_samples)),rep("negative",length(negative_samples)))
mat <- as.matrix(mat[,samples])
mat <- na.omit(mat) # omit rows that contain NA, or will get error in wilcox.test
norm_method <- args$norm_method
m <- args$method
d <- args$datatype

message(paste0("input is ",d," perform diff using ", m))
if(args$datatype=="count") {
# Required columns for a differential expression file: baseMean, log2FoldChange, pvalue, padj
if(grepl('^deseq2_', args$method)){
  suppressPackageStartupMessages(library(DESeq2))
  suppressPackageStartupMessages(library(BiocParallel))
  register(MulticoreParam(args$cores))
  dds <- DESeqDataSetFromMatrix(countData = mat,
                                colData = as.matrix(data.frame(group=group)),
                                design = ~group)
  if(args$method == 'deseq2_wald'){
    dds <- DESeq(dds,test="Wald",parallel=T)
    res <- results(dds, contrast=c('group', 'positive', 'negative'))
    } else if(args$method == 'deseq2_lrt'){
    dds <- DESeq(dds,test="LRT",reduced = ~ 1)
    res <- results(dds, contrast=c('group', 'positive', 'negative'))
    } else if(args$method == 'deseq2_sc'){
    dds <- DESeq(dds,test="LRT",reduced = ~ 1,useT=TRUE, minmu=1e-6, minReplicatesForReplace=Inf)
    res <- results(dds, contrast=c('group', 'positive', 'negative'))
    }
  res <- as.data.frame(res)
} else if(grepl('^edger_', args$method)) {
  suppressPackageStartupMessages(library(edgeR))
  y <- DGEList(counts=mat, samples=samples, group=group)
  y <- calcNormFactors(y, method=norm_method)
  design <- model.matrix(~group)
  y <- estimateDisp(y, design)
  if(args$method == 'edger_glmqlf'){
    fit <- glmQLFit(y, design)
    test <- glmQLFTest(fit, coef=2)
    res <- topTags(test, n=nrow(mat), sort.by='none')
  } else if(args$method == 'edger_glmlrt'){
    fit <- glmFit(y, design)
    test <- glmLRT(fit, coef=2)
    res <- topTags(test, n=nrow(mat), sort.by='none')
  } else if(args$method == 'edger_exact'){
    test <- exactTest(y)
    res <- topTags(test, n=nrow(mat), sort.by='none')
  }
  res <- cbind(res$table, baseMean=2^(res$table$logCPM))
  # rename columns
  mapped_names <- colnames(res)
  for(i in 1:ncol(res)){
    if(colnames(res)[i] == 'logFC'){
      mapped_names[i] <- 'log2FoldChange'
    }else if(colnames(res)[i] == 'PValue'){
      mapped_names[i] <- 'pvalue'
    }else if(colnames(res)[i] == 'FDR') {
      mapped_names[i] <- 'padj'
    }else{
      mapped_names[i] <- colnames(res)[i]
    }
  }
  colnames(res) <- mapped_names
  
} else if(args$method == 'wilcox') {
  suppressPackageStartupMessages(library(edgeR))
  # normalize
  test_func <- function(x){
    wilcox.test(x[group == 'negative'], x[group == 'positive'], alternative='two.sided')$p.value
  }
  #matrix_cpm <- cpm(mat, method=norm_method)
  matrix_cpm <- cpm(mat)
  pvalues <- apply(matrix_cpm, 1, test_func)
  matrix_logcpm <- log2(matrix_cpm + pseudo)
  treatMeans <- apply(matrix_logcpm[,which(group == 'positive')], 1, mean)
  ctrlMeans <- apply(matrix_logcpm[,which(group == 'negative')], 1, mean)
  logFC <- treatMeans - ctrlMeans
  res <- data.frame(log2FoldChange=logFC,
                    pvalue=pvalues, 
                    padj=p.adjust(pvalues, method='BH'),
                    baseMean=apply(matrix_cpm, 1, mean),
                    treatMean=treatMeans,
                    ctrlMean=ctrlMeans)
  #message('Write results to output file: ', args$output_file)
  
} else if(args$method == 'limma'){
  suppressPackageStartupMessages(library(limma))
  suppressPackageStartupMessages(library(edgeR))
  y <- DGEList(counts=mat, samples=samples, group=group)
  y <- calcNormFactors(y, method=norm_method)
  model <- model.matrix(~group)
  y <- voom(y, model, plot=FALSE)
  fit <- lmFit(y, model)
  fit <- eBayes(fit, robust=TRUE, trend=TRUE)
  #fit2 <- contrasts.ft(fit)
  #fit2 <- eBayes(fit2, robust=TRUE, trend=TRUE)
  #top_table <- topTable(fit2, sort.by='none', n=Inf)
  top_table <- topTable(fit, coef=2, sort.by='none', n=Inf)
  # rename columns
  mapped_names <- colnames(top_table)
  for(i in 1:ncol(top_table)){
    if(colnames(top_table)[i] == 'logFC'){
      mapped_names[i] <- 'log2FoldChange'
    }else if(colnames(top_table)[i] == 'P.Value'){
      mapped_names[i] <- 'pvalue'
    }else if(colnames(top_table)[i] == 'adj.P.Val') {
      mapped_names[i] <- 'padj'
    }else{
      mapped_names[i] <- colnames(top_table)[i]
    }
  }
  colnames(top_table) <- mapped_names
  res <- top_table
} else if(args$method == "ttest") {
  suppressPackageStartupMessages(library(genefilter))
  res <- rowttests(mat, as.factor(group))
  res$padj <- p.adjust(res$p.value, method='BH')
  res$log2FoldChange <- rowMeans(mat[, group == 'positive']) - rowMeans(mat[, group == 'negative'])
  res$treatMean <- rowMeans(mat[, group == 'positive'])
  res$ctrlMean <- rowMeans(mat[, group == 'negative'])
  
} else {stop('unknown differential expression method: ', args$method)}

} else if(args$datatype=="normalized") {
  if (max(as.matrix(mat),na.rm = T)>50 & min(as.matrix(mat),na.rm = T)>0){
    message("input seeems like not log2tpm matrix or gsva outfile, convert to log2 ")
    matrix_log <- log2(mat + pseudo)
  } else {
    message("input seeems like log2tpm matrix or other normalized matrix, no normalization performed !")
    matrix_log <- mat
  }
  
  if(args$method == 'wilcox') {
    test_func <- function(x){
      wilcox.test(x[group == 'negative'], x[group == 'positive'], alternative='two.sided')$p.value
    }
    pvalues <- apply(matrix_log, 1, test_func)
    treatMeans <- apply(matrix_log[,which(group == 'positive')], 1, mean)
    ctrlMeans <- apply(matrix_log[,which(group == 'negative')], 1, mean)
    logFC <- treatMeans - ctrlMeans
    res <- data.frame(log2FoldChange=logFC,
                      pvalue=pvalues, 
                      padj=p.adjust(pvalues, method='BH'),
                      baseMean=apply(mat, 1, mean),
                      treatMean=treatMeans,
                      ctrlMean=ctrlMeans)
    #message('Write results to output file: ', args$output_file)
    
  } else if(args$method == "ttest") {
    suppressPackageStartupMessages(library(genefilter))
    res <- rowttests(matrix_log, as.factor(group))
    res$padj <- p.adjust(res$p.value, method='BH')
    res$log2FoldChange <- rowMeans(matrix_log[, group == 'positive']) - rowMeans(matrix_log[, group == 'negative'])
    res$treatMean <- rowMeans(matrix_log[, group == 'positive'])
    res$ctrlMean <- rowMeans(matrix_log[, group == 'negative'])
    res$baseMean <- apply(mat, 1, mean)  # new adding
    
    } else {stop('unknown differential expression method (for normalized matrix input): ', args$method)}
} else (stop('unknown matrix datatype !'))

# sort diff table
if(args$sort==T){
  message('Sort results according to pvalue and log2FoldChange')
  res <- res[order(-res$pvalue,abs(res$log2FoldChange),decreasing = T),]
} else {
  message('Not sort results according to pvalue and log2FoldChange')
}

# write results to file
message('Write results to output file: ', args$output_file)
res <- cbind(rownames(res),res)
colnames(res)[1] <- "id"

write.table(res, args$output_file, sep='\t', quote=FALSE, row.names=FALSE)
