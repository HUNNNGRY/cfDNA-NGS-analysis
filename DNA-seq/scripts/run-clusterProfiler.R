# run clusterProfiler
# last at 20210905 by pengfei

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(DOSE))
suppressPackageStartupMessages(library(enrichplot))
suppressPackageStartupMessages(library(org.Hs.eg.db))

parser <- ArgumentParser(description='GO/KEGG/DO geneset enrichment by ORA/GSEA using diff table')
parser$add_argument('-i', '--matrix', type='character', required=TRUE,
                    help='input diff matrix. Rows are genes (^ENSG*, *_ENSG*). Columns are deseq2 tags')
parser$add_argument('-c', '--cutoff', type='double', default=0.05,
                    help='Over-Representation Analysis cutoff pvalue (default=0.05, genes with pvalue < 0.05 will be treated as top rank genes)')
parser$add_argument('-p', '--pValue', type='double', default=0.05,
                    help='pValue cutoff to be saved in output file (default=0.05, will be qvalue if adjustPval is TRUE)')
parser$add_argument('-a', '--adjustPval', type='logical', default=FALSE,
                    help='whether to adjust pValue and cutoff (default=FALSE, if TRUE, then defined pValue wil send to qValue)')
parser$add_argument('-o', '--outdir', type='character', default="./",
                    help='outdir (default=./)')
parser$add_argument('-f', '--fig', type='logical', default=TRUE,
                    help='whether to plot figures (default=TRUE)')
args <- parser$parse_args()
matrix <- args$matrix
cutoff <- args$cutoff
pValue <- args$pValue
adjustPval <- args$adjustPval
outdir <- args$outdir
fig <- args$fig

if(!adjustPval)
{
  message("not adjust pvalue")
  pAdjust <- "none"
  qValue <- 1
} else {
  message("adjust pvalue by BH")
  pAdjust <- "BH"
  qValue <- pValue
}

# test 
# matrix <- "/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/output/lulab/diff/promoter/CRCvsNC-passQC-edger_glmlrt/CRCvsNC-passQC-edger_glmlrt.txt"
# cutoff <- 0.05
# pValue <- 0.1
# adjustPval <- FALSE
# outdir <- "./run-clusterprofiler-test"
# fig <- TRUE

dir.create(outdir,recursive = T, showWarnings = F)
res <- read.table(matrix,stringsAsFactors = F)

if(!adjustPval){
  message(paste0('pvalue<=',cutoff,':  '), sum(res$pvalue<=cutoff, na.rm = T),
          " up:",sum(res$pvalue<=cutoff & res$log2FoldChange>0, na.rm = T),
          " down:",sum(res$pvalue<=cutoff & res$log2FoldChange<0, na.rm = T))
} else{
  message(paste0('padj<=',cutoff,':  '), sum(res$padj<=cutoff, na.rm = T),
          " up:",sum(res$padj<=cutoff & res$log2FoldChange>0, na.rm = T),
          " down:",sum(res$padj<=cutoff & res$log2FoldChange<0, na.rm = T))
}

resOrder <- res[order(res$log2FoldChange,-res$pvalue,decreasing = T),]  # order rows by logFC and then pVal
resOrder$name <- rownames(resOrder)
resOrder$name <- as.character(lapply(strsplit(resOrder$name,".",fixed = TRUE),function(x) x[1]))

if (length(grep("^ENSG",rownames(resOrder)))>=2) {
	message("input gene id is ^ENSG*")
	resOrder$name <- resOrder$name
} else {
        message("input gene id is not ^ENSG*")
	resOrder$name <- as.character(lapply(strsplit(resOrder$name,"_",fixed = TRUE),function(x) x[2]))
}
#resOrder$name <- as.character(lapply(strsplit(resOrder$name,"_",fixed = TRUE),function(x) x[2]))
resOrder <- resOrder[!duplicated(resOrder$name),]
gene.list.ncbiid <- bitr(resOrder$name, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
resOrder <- merge(x = resOrder, y = gene.list.ncbiid, by.x="name",by.y="ENSEMBL")
resOrder <- resOrder[!duplicated(resOrder$ENTREZID),]
resOrder <- resOrder[order(resOrder$log2FoldChange,-resOrder$pvalue,decreasing = T),]  # order rows by logFC and pVal
resOrderPval <- resOrder[order(-resOrder$pvalue,resOrder$log2FoldChange,decreasing = T),]  # order rows by pVal then logFC
# get all geneList (rank)
gene.list.ncbiid <- resOrder$log2FoldChange
names(gene.list.ncbiid) <- resOrder$ENTREZID
gene.list.ensgid <- resOrder$log2FoldChange
names(gene.list.ensgid) <- resOrder$name
# get top geneList
gene.top.ncbiid <- resOrderPval$ENTREZID[resOrderPval$pvalue<=cutoff]


# 1.Over-Representation Analysis
message("start ORA enrich...")
enrich.GO<- enrichGO(
  gene = gene.top.ncbiid, 
  universe = names(gene.list.ncbiid),
  OrgDb = org.Hs.eg.db,
  ont = "ALL", 
  keyType = "ENTREZID",  
  pvalueCutoff = pValue,
  pAdjustMethod = pAdjust, 
  qvalueCutoff = qValue
)

enrich.KEGG <- enrichKEGG(
  gene = gene.top.ncbiid, 
  keyType = "ncbi-geneid", 
  pvalueCutoff = pValue,
  pAdjustMethod = pAdjust,
  qvalueCutoff = qValue
)

enrich.DO <- enrichDO(
  gene  = gene.top.ncbiid, 
  ont = "DO",  
  pvalueCutoff = pValue,
  pAdjustMethod = pAdjust,
  qvalueCutoff = qValue,
  universe = names(gene.list.ncbiid),
  minGSSize = 5,
  maxGSSize = 500,
  readable = FALSE  
  )

# 2.Gene Set Enrichment Analysis
message("start gsea...")
gse.GO <- gseGO(
  geneList = gene.list.ncbiid, 
  OrgDb = org.Hs.eg.db,
  ont = "ALL", 
  keyType = "ENTREZID",  
  pvalueCutoff = pValue,
  pAdjustMethod = pAdjust
  #by = "fgsea", #'DOSE'
)

gse.KEGG <- gseKEGG(
  geneList = gene.list.ncbiid,
  keyType = 'kegg', 
  organism = 'hsa',
  nPerm  = 1000,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = pValue,
  pAdjustMethod = pAdjust
  #by = "fgsea",
  #seed = T  
)

gse.DO <- gseDO(
  geneList = gene.list.ncbiid,
  nPerm  = 1000,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = pValue,
  pAdjustMethod = pAdjust
  #by = "fgsea",
  #seed = T  
)


# 3.write table output
types <- c(enrich.GO,enrich.KEGG,enrich.DO,gse.GO,gse.KEGG,gse.DO)
names(types) <- c("enrich.GO","enrich.KEGG","enrich.DO","gse.GO","gse.KEGG","gse.DO")

message("start writing tables...")
for(i in 1:length(types)){
  write.table(types[[i]],paste0(outdir,"/", names(types)[i], ".txt"),row.names = F,col.names = T, quote = F,sep = "\t")
}

# plot figure output
if(fig){
# 4.plot
message("start plotting...")
dir.create(paste0(outdir,"/plot/"),showWarnings = F)
## 4.1.barplot for top gene enrich (not informative...)
#barplot(enrich.GO, showCategory = 10)  
#barplot(enrich.KEGG, showCategory = 10) 

## 4.2.dotplot for gene set enrich
try(                   # keep running when plot error exists
for (i in 1:length(types)){
  pdf(paste0(outdir,"/plot/dotplot-",names(types)[i],".pdf"))
  print(enrichplot::dotplot(types[[i]], font.size = 9, title=names(types)[i],showCategory = 20))  # must print in for loop !
  dev.off()
}
)
#cowplot::plot_grid(p21, p22, p23, p24, ncol=2, labels=paste(LETTERS[1:4],txt,sep = ": "))


## 4.3.emapplot for top gene and gene set enrich?  (pathway relations)
try(
for (i in 1:length(types)){
  pdf(paste0(outdir,"/plot/emapplot-",names(types)[i],".pdf"))
  print(enrichplot::emapplot(types[[i]],showCategory = 20))
  dev.off()
}
)
## 4.4.cnetplot top gene and gene set enrich? (genes centric relations)
try(
for (i in 1:length(types)){
  pdf(paste0(outdir,"/plot/cnetplot-",names(types)[i],".pdf"))
  print(enrichplot::cnetplot(types[[i]], showCategory = 5))
  dev.off()
}
)
}
