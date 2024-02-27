# get matrix and label for roc
getwd()
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq")
tpm <- read.table("./output/lulab/matrix/TPM_matrix_promoter.txt",stringsAsFactors = F,header = T, sep = "\t",check.names = F)
tpm[1:4,1:4]

tpm$gene_id <- as.character(lapply(strsplit(x = tpm$gene_id,"\\."), function(x) x[1]))
tpm$gene_id <- as.character(lapply(strsplit(x = tpm$gene_id,"\\_"), function(x) x[2]))

dim(tpm)
length(unique(tpm$gene_id))
tpm <- tpm[!duplicated(tpm$gene_id),]

pd <- read.table("./meta/lulab/sample_table.txt",stringsAsFactors = F,check.names = F,sep = "\t",header = T)
pd <- pd[pd$overall_QC=="Y",]

tpm.filtr <- tpm[,colnames(tpm) %in% pd$data_ID]
tpm.filtr$ID <- tpm$gene_id
tpm.filtr <- tpm.filtr[,c(ncol(tpm.filtr),1:(ncol(tpm.filtr)-1))]
dim(tpm.filtr)
tpm.filtr[1:4,1:4]
table(pd$group)

write.table(tpm.filtr,"./output/lulab/dmr/tpm-matrix-promoter-filter-ensg.txt",quote = F,sep = "\t",row.names = F,col.names = T)
write.table(pd$group,"./output/lulab/dmr/71-label.txt",quote = F,sep = "\t",row.names = F,col.names = F)
