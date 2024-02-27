
rm(list = ls())
options(stringsAsFactors = F)
setwd("/BioII/lulab_b/baopengfei/projects/exOmics")

meta <- read.table("../multi-omics-explore/meta/lulab/sample_table_pair.txt",header = T,check.names = F,sep = "\t")
#meta[[is.na(meta)]]

datatype <- "DNA" #DIP
region <- "promoter150TSS50"  #promoter150TSS50  gene
pos.grp <- "CRC"
neg.grp <- "NC"
feature.num <- 200
criteria <- "passQCstringent"

# for (criteria in c("passQCstringent","passQCloose")){
# for(pos.grp in c("CRC","STAD")){
medip.criteria <- paste0("medip_",criteria)
wgs.criteria <- paste0("wgs_",criteria)


medip.id <- meta$medip_data_id[meta[[medip.criteria]]=="Y"]  
medip.id <- medip.id[c(grep(pos.grp,medip.id),grep(neg.grp,medip.id))]
wgs.id <- meta$wgs_data_id[meta[[wgs.criteria]]=="Y"]  
wgs.id <- wgs.id[c(grep(pos.grp,wgs.id),grep(neg.grp,wgs.id))]

cnv.tpm <- read.table(paste0("./",datatype,"-seq/output/lulab/matrix/CPM_matrix_",region,".txt"),header = T,row.names = 1,check.names = F,sep = "\t")
#cnv.tpm[1:3,1:3]
cnv.tpm <- log2(cnv.tpm+1)

cnv.diff <- read.table(paste0("./",datatype,"-seq/output/lulab/diff/",region,"/",pos.grp,"vs",neg.grp,"-",criteria,"-deseq2_wald/",pos.grp,"vs",neg.grp,"-",criteria,"-deseq2_wald.txt"),header = T,row.names = 1,check.names = F,sep = "\t")
cnv.diff <- cnv.diff[order(-cnv.diff$pvalue,cnv.diff$log2FoldChange,decreasing = T),]
#table(is.na(cnv.diff$pvalue))
top200 <- head(rownames(cnv.diff),feature.num)

cnv.tpm <- cnv.tpm[top200,wgs.id]
cnv.tpm <- cnv.tpm[,wgs.id]

suppressPackageStartupMessages(library(pheatmap))
annotation_col <- as.data.frame(colnames(cnv.tpm))
#annotation_col$sampletype <- neg.grp
annotation_col$sampletype[grep(pattern = neg.grp,x = annotation_col[,1])] <- neg.grp

annotation_col$sampletype[grep(pattern = pos.grp,x = annotation_col[,1])] <- pos.grp
#annotation_col$sampletype[grep(pattern = rm.grp,x = annotation_col[,1])] <- rm.grp
annotation_col <- annotation_col[annotation_col$sampletype==pos.grp | annotation_col$sampletype==neg.grp,]
  
rownames(annotation_col) <- annotation_col[,1]
annotation_col <- annotation_col[-1]

pheatmap(mat = cnv.tpm,
         annotation_col = annotation_col, #data frame that specifies the annotations shown on left side of the heatmap.
         scale = "row",
         cluster_cols = T,annotation_colors = list(sampletype = c( NC = "steelblue", CRC = "orange", STAD = "firebrick")),
         cluster_rows = T,
         show_colnames=F, show_rownames=F, 
         colorRampPalette(c("steelblue","white","red"))(100),
         fontsize_row = 5,
         filename = paste0("./",datatype,"-",region,"-",pos.grp,"vs",neg.grp,"-",criteria,"-",feature.num,"-heatmap.pdf"))
# }
# }
