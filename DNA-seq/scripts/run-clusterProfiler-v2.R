# run clusterProfiler
# last at 20211022 by pengfei (warning: p and q the same)
options(stringsAsFactors = F)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(DOSE))
suppressPackageStartupMessages(library(enrichplot))
suppressPackageStartupMessages(library(org.Hs.eg.db))

parser <- ArgumentParser(description='GO/KEGG/DO geneset enrichment by ORA/GSEA using diff table (network independent version)')
parser$add_argument('-i', '--matrix', type='character', required=TRUE,
                    help='input diff matrix. Rows are genes (^ENSG*, *_ENSG*). Columns are deseq2 tags')
parser$add_argument('-c', '--cutoff', type='double', default=0.05,
                    help='Over-Representation Analysis cutoff pvalue (default=0.05, genes with pvalue < 0.05 will be treated as top rank genes)')
parser$add_argument('-p', '--pValue', type='double', default=0.05,
                    help='pValue cutoff to be saved in output file (default=0.05, will be qvalue if adjustPval is TRUE)')
parser$add_argument('-a', '--adjustPval', type='logical', default=FALSE,
                    help='whether to adjust pValue and cutoff (default=FALSE, not use this flag, if TRUE, then defined pValue wil feed to qValue)')
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

# message(cutoff)
# message(pValue)
# message(adjustPval)
# # test 
# matrix <- "/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/diff/gene/CRCvsNC-passQC-edger_glmlrt/CRCvsNC-passQC-edger_glmlrt.txt"
# cutoff <- 0.05
# pValue <- 1
# adjustPval <- FALSE
# outdir <- "./run-clusterprofiler-test"
# fig <- TRUE

dir.create(outdir,recursive = T, showWarnings = F)
res <- read.table(matrix,stringsAsFactors = F)


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


if(!adjustPval){
  message(paste0('output table contains records within pvalue<=',pValue,':  '), sum(res$pvalue<=pValue, na.rm = T),
          " up:",sum(res$pvalue<=pValue & res$log2FoldChange>0, na.rm = T),
          " down:",sum(res$pvalue<=pValue & res$log2FoldChange<0, na.rm = T))
  message(paste0('ORA enrich selects records within pvalue<=',cutoff,':  '), sum(res$pvalue<=cutoff, na.rm = T),
          " up:",sum(res$pvalue<=cutoff & res$log2FoldChange>0, na.rm = T),
          " down:",sum(res$pvalue<=cutoff & res$log2FoldChange<0, na.rm = T))
} else{
  message(paste0('output table contains records within padj<=',qValue,':  '), sum(res$padj<=pValue, na.rm = T),
          " up:",sum(res$padj<=pValue & res$log2FoldChange>0, na.rm = T),
          " down:",sum(res$padj<=pValue & res$log2FoldChange<0, na.rm = T))
  message(paste0('ORA enrich selects records within padj<=',cutoff,':  '), sum(res$padj<=cutoff, na.rm = T),
          " up:",sum(res$padj<=cutoff & res$log2FoldChange>0, na.rm = T),
          " down:",sum(res$padj<=cutoff & res$log2FoldChange<0, na.rm = T))
}


resOrder <- res[order(res$log2FoldChange,-res$pvalue,decreasing = T),]  # order rows by logFC and then pVal
resOrder$ENSG <- rownames(resOrder)
resOrder$ENSG <- as.character(lapply(strsplit(resOrder$ENSG,".",fixed = TRUE),function(x) x[1]))

if (length(grep("^ENSG",rownames(resOrder)))>=2) {
	message("input gene id is ^ENSG*")
	resOrder$ENSG <- resOrder$ENSG
} else {
        message("input gene id is not ^ENSG*")
	resOrder$ENSG <- as.character(lapply(strsplit(resOrder$ENSG,"_",fixed = TRUE),function(x) x[2]))
}
resOrder <- resOrder[!duplicated(resOrder$ENSG),]
#gene.list.ncbiid <- bitr(resOrder$ENSG, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")


#ref <- read.table("/BioII/lulab_b/baopengfei/shared_reference/hg38/gene.gtf",sep = "\t",header = T)
# library(gprofiler2)
# tmp <- gprofiler2::gconvert(query = ref$ensg, organism = "hsapiens", target = "ENTREZGENE_ACC" )
# tmp <- tmp[,c("input","target")]
# colnames(tmp) <- c("ensg","ENTREZID")
# ref$ENTREZID <- tmp$ENTREZID[match(ref$ensg,tmp$ensg)]
# ref <- ref[!duplicated(ref$ensg),]
# table(duplicated(ref$ENTREZID))
# write.table(ref,"/BioII/lulab_b/baopengfei/shared_reference/geneset/encodev27_gene.txt",quote = F,sep = "\t",col.names = T,row.names = F)
gene.list.ncbiid <- read.table("/BioII/lulab_b/baopengfei/shared_reference/geneset/encodev27_gene.txt",sep = "\t",header = T)
gene.list.ncbiid <- gene.list.ncbiid[!is.na(gene.list.ncbiid$ENTREZID),]
gene.list.ncbiid <- gene.list.ncbiid[,c("ensg","ENTREZID")]

resOrder <- merge(x = resOrder, y = gene.list.ncbiid, by.x="ENSG",by.y="ensg")
resOrder <- resOrder[!duplicated(resOrder$ENTREZID),]
resOrder <- resOrder[order(resOrder$log2FoldChange,-resOrder$pvalue,decreasing = T),]  # order rows by logFC and pVal
resOrderPval <- resOrder[order(-resOrder$pvalue,resOrder$log2FoldChange,decreasing = T),]  # order rows by pVal then logFC

# get all geneList (rank)
gene.list.ncbiid <- resOrder$log2FoldChange
names(gene.list.ncbiid) <- resOrder$ENTREZID
gene.list.ensgid <- resOrder$log2FoldChange
names(gene.list.ensgid) <- resOrder$ENSG

if(!adjustPval){
  # get top geneList
  gene.top.ncbiid <- resOrderPval$ENTREZID[resOrderPval$pvalue<=cutoff]   # pvalue
  # get top geneList (up-regulated)
  gene.top.ncbiid.up <- resOrderPval$ENTREZID[resOrderPval$pvalue<=cutoff & resOrder$log2FoldChange>0]
  # get top geneList (down-regulated)
  gene.top.ncbiid.down <- resOrderPval$ENTREZID[resOrderPval$pvalue<=cutoff & resOrder$log2FoldChange<0]
} else {
 gene.top.ncbiid <- resOrderPval$ENTREZID[resOrderPval$padj<=cutoff]   # padj
 gene.top.ncbiid.up <- resOrderPval$ENTREZID[resOrderPval$padj<=cutoff & resOrder$log2FoldChange>0]
 gene.top.ncbiid.down <- resOrderPval$ENTREZID[resOrderPval$padj<=cutoff & resOrder$log2FoldChange<0]
}




# 1.Over-Representation Analysis
message("start ORA enrich...")
# 1.1 ORA all gene
enrich.GO<- enrichGO(
  gene = gene.top.ncbiid, 
  #universe = names(gene.list.ncbiid),
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
  #universe = names(gene.list.ncbiid),
  pvalueCutoff = pValue,
  pAdjustMethod = pAdjust,
  qvalueCutoff = qValue
)
#enrich.KEGG@result

enrich.DO <- enrichDO(
  gene  = gene.top.ncbiid, 
  ont = "DO",  
  pvalueCutoff = pValue,
  pAdjustMethod = pAdjust,
  qvalueCutoff = qValue,
  #universe = names(gene.list.ncbiid),
  minGSSize = 5,
  maxGSSize = 500,
  readable = FALSE  
  )

# 1.2 ORA up gene
enrich.GO.up<- enrichGO(
  gene = gene.top.ncbiid.up, 
  #universe = names(gene.list.ncbiid),
  OrgDb = org.Hs.eg.db,
  ont = "ALL", 
  keyType = "ENTREZID",  
  pvalueCutoff = pValue,
  pAdjustMethod = pAdjust, 
  qvalueCutoff = qValue
)

enrich.KEGG.up <- enrichKEGG(
  gene = gene.top.ncbiid.up, 
  keyType = "ncbi-geneid", 
  #universe = names(gene.list.ncbiid),
  pvalueCutoff = pValue,
  pAdjustMethod = pAdjust,
  qvalueCutoff = qValue
)

enrich.DO.up <- enrichDO(
  gene  = gene.top.ncbiid.up, 
  ont = "DO",  
  pvalueCutoff = pValue,
  pAdjustMethod = pAdjust,
  qvalueCutoff = qValue,
  #universe = names(gene.list.ncbiid),
  minGSSize = 5,
  maxGSSize = 500,
  readable = FALSE  
  )

# 1.3 ORA down gene
enrich.GO.down<- enrichGO(
  gene = gene.top.ncbiid.down, 
  #universe = names(gene.list.ncbiid),
  OrgDb = org.Hs.eg.db,
  ont = "ALL", 
  keyType = "ENTREZID",  
  pvalueCutoff = pValue,
  pAdjustMethod = pAdjust, 
  qvalueCutoff = qValue
)

enrich.KEGG.down <- enrichKEGG(
  gene = gene.top.ncbiid.down, 
  keyType = "ncbi-geneid", 
  #universe = names(gene.list.ncbiid),
  pvalueCutoff = pValue,
  pAdjustMethod = pAdjust,
  qvalueCutoff = qValue
)
#enrich.KEGG@result

enrich.DO.down <- enrichDO(
  gene  = gene.top.ncbiid.down, 
  ont = "DO",  
  pvalueCutoff = pValue,
  pAdjustMethod = pAdjust,
  qvalueCutoff = qValue,
  #universe = names(gene.list.ncbiid),
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
#gse.GO@result

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
#gse.KEGG@result

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
types <- list(enrich.GO,enrich.KEGG,enrich.DO,
  enrich.GO.up,enrich.KEGG.up,enrich.DO.up,
  enrich.GO.down,enrich.KEGG.down,enrich.DO.down,
  gse.GO,gse.KEGG,gse.DO)
names(types) <- c("enrich.GO","enrich.KEGG","enrich.DO",
  "enrich.GO.up","enrich.KEGG.up","enrich.DO.up",
  "enrich.GO.down","enrich.KEGG.down","enrich.DO.down",
  "gse.GO","gse.KEGG","gse.DO")

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
# try(
# for (i in 1:length(types)){
#   pdf(paste0(outdir,"/plot/emapplot-",names(types)[i],".pdf"))
#   print(enrichplot::emapplot(types[[i]],showCategory = 20))
#   dev.off()
# }
# )
## 4.4.cnetplot top gene and gene set enrich? (genes centric relations)
# try(
# for (i in 1:length(types)){
#   pdf(paste0(outdir,"/plot/cnetplot-",names(types)[i],".pdf"))
#   print(enrichplot::cnetplot(types[[i]], showCategory = 5))
#   dev.off()
# }
# )

}
