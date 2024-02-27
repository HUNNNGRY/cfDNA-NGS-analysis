#! /usr/bin/env Rscript
#last update at 20210428 by pengfei

suppressPackageStartupMessages(library(argparse))
parser <- ArgumentParser(description='Differential expression')
parser$add_argument('-i', '--matrix', type='character', required=TRUE,
                    help='input count matrix. Rows are genes. Columns are samples.')
parser$add_argument('-p', '--positive-ids', type='character', required=TRUE,
                    help='file contains ids of positive class')
parser$add_argument('-n', '--negative-ids', type='character', required=TRUE,
                    help='file contains ids of negative class')
parser$add_argument('-m', '--method', type='character', default="edger_glmlrt",
                    choices=c('deseq2_wald','deseq2_lrt', 'edger_exact', 'edger_glmqlf', 'edger_glmlrt', 'wilcox', 'limma', 'ttest'),
                    help='differential expression method to use')
parser$add_argument('--norm-method', type='character', default='TMM',
                    choices=c('RLE', 'CPM', 'TMM', 'upperquartile'))
parser$add_argument('--pseudo-count', type='double', default=1.0,
                    help='pseudo-count added to log2 transform in ttest')
parser$add_argument('--cores', type='double', default=1,
                    help='multi-cores using for deseq2, default:1')
parser$add_argument('-o', '--output-file', type='character', required=TRUE,
                    help='output file')
args <- parser$parse_args()

message('read count matrix: ', args$matrix)
mat <- read.table(args$matrix, header = TRUE, row.names=1, check.names=FALSE, sep='\t')

message('Read positive sample ids')
positive_samples <- read.delim(args$positive_ids,sep="\n",stringsAsFactors=F,header=F)[,1]
message('Read negative sample ids')
negative_samples <- read.delim(args$negative_ids,sep="\n",stringsAsFactors=F,header=F)[,1]
message('Number of positive samples: ', length(positive_samples))
message('Number of negative samples: ', length(negative_samples))


samples <- c(positive_samples, negative_samples)
group <- c(rep("positive",length(positive_samples)),rep("negative",length(negative_samples)))
mat <- as.matrix(mat[,samples])

norm_method <- args$norm_method
message('perform differential expression using ', args$method)
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
    }
  write.table(as.data.frame(res), args$output_file, sep='\t', quote=FALSE, row.names=TRUE)
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
  
  # write results to file
  message('Write results to output file: ', args$output_file)
  write.table(res, args$output_file, sep='\t', quote=FALSE, row.names=TRUE)
} else if(args$method == 'wilcox') {
  suppressPackageStartupMessages(library(edgeR))
  # normalize
  matrix_cpm <- cpm(mat, method=norm_method)
  test_func <- function(x){
    wilcox.test(x[group == 'negative'], x[group == 'positive'], alternative='two.sided')$p.value
  }
  pvalues <- apply(matrix_cpm, 1, test_func)
  matrix_logcpm = log2(matrix_cpm + 1)
  treatMeans <- apply(matrix_logcpm[,which(group == 'positive')], 1, mean)
  ctrlMeans <- apply(matrix_logcpm[,which(group == 'negative')], 1, mean)
  logFC <- treatMeans - ctrlMeans
  res <- data.frame(log2FoldChange=logFC,
                    pvalue=pvalues, 
                    padj=p.adjust(pvalues, method='BH'),
                    baseMean=apply(matrix_cpm, 1, mean),
                    treatMean=treatMeans,
                    ctrlMean=ctrlMeans)
  message('Write results to output file: ', args$output_file)
  write.table(res, args$output_file, sep='\t', quote=FALSE, row.names=TRUE)
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
  
  # write results to file
  message('Write results to output file: ', args$output_file)
  write.table(top_table, args$output_file, sep='\t', quote=FALSE, row.names=TRUE)
} else if(args$method == "ttest") {
  suppressPackageStartupMessages(library(genefilter))
  #mat <- log2(mat + args$pseudo_count)
  res <- rowttests(mat, as.factor(group))
  res$padj <- p.adjust(res$p.value, method='BH')
  res$log2FoldChange <- rowMeans(mat[, group == 'positive']) - rowMeans(mat[, group == 'negative'])
  res$treatMean <- rowMeans(mat[, group == 'positive'])
  res$ctrlMean <- rowMeans(mat[, group == 'negative'])
  write.table(res, args$output_file, sep='\t', quote=FALSE, row.names=TRUE)
} else {
  stop('unknown differential expression method: ', args$method)
}
