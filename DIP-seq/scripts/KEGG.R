library(clusterProfiler)
library(org.Hs.eg.db)
args<-commandArgs(TRUE)
infile <- args[1]
outfile <- args[2]
#bg.ids.path <- args[3]
gene.list <- as.character(read.table(infile)[,1])
#bg.entrez.ids <- as.character(read.table(bg.ids.path)[,1])
#kegg <- enrichKEGG(gene.list,organism  = 'human', keyType="kegg",universe=bg.entrez.ids,pAdjustMethod = 'BH',pvalueCutoff = 1, qvalueCutoff = 1)
kegg <- enrichKEGG(gene.list,organism  = 'human', keyType="kegg",pAdjustMethod = 'BH',pvalueCutoff = 1, qvalueCutoff = 1)
write.table(kegg,file=outfile,sep="\t",quote=FALSE)
