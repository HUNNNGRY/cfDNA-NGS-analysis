# sum diff table for candidates selection
# last 210115 by bpf


# sum medip promoter CPM-TMM (passQCstringent) -----------------------------------------------------
cmp <- c("STADvsNC","CRCvsNC","CRC_STADvsNC")
for (i in cmp){
  #i <- "STADvsNC"
  print(i)
  tmp <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/output/lulab/diff-medip-TMM/promoter/",i,"-passQCstringent-edger_exact/",i,"-passQCstringent-edger_exact.txt"),sep = "\t",header = T)
  tmp <- tmp[order(-tmp$pvalue,abs(tmp$log2FoldChange),decreasing = T),]
  #tmp$ENSG <- unlist(lapply(strsplit(tmp$id,"_"),function(x) x[2]))
  tmp$ENSG <- unlist(lapply(strsplit(rownames(tmp),"|",fixed = T),function(x) x[1]))
  tmp$ENSG <- gsub("promoter_","",tmp$ENSG)
  colnames(tmp)[3] <- "log2CPM_Mean"
  
  ## only keep top500
  #tmp <- tmp[1:500,]
  
  ref <- read.table("/BioII/lulab_b/baopengfei/shared_reference/hg38/promoter.gtf",sep = "\t",header = F)
  #ref[1:3,]
  colnames(ref) <- c("chr","source","region","start","end","score","strand","phase","ENSG")
  
  ref[["promoter_region"]] <- paste0(ref$chr,":",ref$start,"-",ref$end)
  ref$ENSG <- gsub("gene_id promoter_","",ref$ENSG)
  ref <- ref[,c("chr","strand","ENSG","promoter_region")]
  ref <- ref[!duplicated(ref$ENSG),]
  
  tmp <- dplyr::left_join(tmp,ref)         
  str(tmp)
  
  ## expr mat 
  mat <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/output/lulab/matrix/CPM-TMM_matrix_promoter.txt",sep = "\t",header = T,row.names = 1, check.names = F)
  samp <- read.table("/BioII/lulab_b/baopengfei/projects/multi-omics-explore/meta/lulab/multiomic-medip-passQCstringent.txt",sep = "\t",header = F)$V1
  mat <- mat[,samp]
  #mat[1:4,1:4]
  mat <- log2(mat+1)
  mat$ENSG <- rownames(mat)
  mat$ENSG <- unlist(lapply(strsplit(rownames(mat),"|",fixed = T),function(x) x[1]))
  mat$ENSG <- gsub("promoter_","",mat$ENSG)
  
  ## join with expr mat
  l <- colnames(tmp)
  tmp[["log2CPM-TMM-promoter"]] <- ""
  tmp <- dplyr::left_join(tmp,mat)         
  
  ## add gini (only interpretable for non-negative quantities)
  #avoid minus value
  #tmp.exp <- 2^tmp[,grep("*-PKU-",colnames(tmp))]     # only for minus value
  tmp.exp <- tmp[,grep("*-PKU-",colnames(tmp))]
  
  tmp[["CRC_gini"]] <- apply(X = tmp.exp[,grep("CRC-PKU-",colnames(tmp.exp))], MARGIN = 1, FUN = function(x) edgeR::gini(as.numeric(x)))
  tmp[["NC_gini"]] <- apply(X = tmp.exp[,grep("NC-PKU-",colnames(tmp.exp))], MARGIN = 1, FUN = function(x) edgeR::gini(as.numeric(x)))
  tmp[["STAD_gini"]] <- apply(X = tmp.exp[,grep("STAD-PKU-",colnames(tmp.exp))], MARGIN = 1, FUN = function(x) edgeR::gini(as.numeric(x)))
  #hist(tmp$NC_gini)
  tmp <- tmp[,c(l,"CRC_gini","NC_gini","STAD_gini","log2CPM-TMM-promoter",samp)]
  
  #colnames(tmp)[1] <- "id"
  write.table(tmp,paste0("/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/output/lulab/diff-medip-TMM/promoter/",i,"-passQCstringent-edger_exact/",i,"-passQCstringent-edger_exact.all.txt"),quote = F,sep = "\t",col.names = T,row.names = F)
  rio::export(tmp,paste0("/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/output/lulab/diff-medip-TMM/promoter/",i,"-passQCstringent-edger_exact/",i,"-passQCstringent-edger_exact.all.xlsx"))   # top500Pvalue
}




# 
# mat <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/count_matrix_promoter150TSS50.txt",sep = "\t",header = T,row.names = 1,check.names = F)
# ratio <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/count_matrix_promoter150TSS50.txt.summary",sep = "\t",header = T,row.names = 1,check.names = F)
# colnames(ratio) <- gsub("output/lulab/bam-sorted-deduped/","",colnames(ratio),fixed = T)
# colnames(ratio) <- gsub(".bam","",colnames(ratio),fixed = T)
# ratio <- ratio[,colnames(mat)]
# 
# cor.test(as.numeric(ratio[1,]),as.numeric(colSums(mat)),method = "pearson")  # method = "spearman"



# sum medip promoter CPM-TMM (passQCstringent,filter QIAGEN kit, 2201117) -----------------------------------------------------
cmp <- c("STADvsNC","CRC_STADvsNC")  # CRCvsNC
for (i in cmp){
  #i <- "STADvsNC"
  print(i)
  tmp <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/output/lulab/diff-medip-TMM/promoter/",i,"-passQCstringent-QIAGEN-edger_exact/",i,"-passQCstringent-QIAGEN-edger_exact.txt"),sep = "\t",header = T)  # ,row.names = 1
  tmp <- tmp[order(-tmp$pvalue,abs(tmp$log2FoldChange),decreasing = T),]
  #tmp$ENSG <- unlist(lapply(strsplit(tmp$id,"_"),function(x) x[2]))
  tmp$ENSG <- unlist(lapply(strsplit(tmp$id,"|",fixed = T),function(x) x[1]))
  tmp$ENSG <- gsub("promoter_","",tmp$ENSG)
  rownames(tmp) <- tmp$ENSG
  colnames(tmp)[3] <- "log2CPM_Mean"
  
  ## only keep top500
  #tmp <- tmp[1:500,]
  
  ref <- read.table("/BioII/lulab_b/baopengfei/shared_reference/hg38/promoter.gtf",sep = "\t",header = F)
  #ref[1:3,]
  colnames(ref) <- c("chr","source","region","start","end","score","strand","phase","ENSG")
  
  ref[["promoter_region"]] <- paste0(ref$chr,":",ref$start,"-",ref$end)
  ref$ENSG <- gsub("gene_id promoter_","",ref$ENSG)
  ref <- ref[,c("chr","strand","ENSG","promoter_region")]
  ref <- ref[!duplicated(ref$ENSG),]
  
  tmp <- dplyr::left_join(tmp,ref)         
  str(tmp)
  
  ## expr mat 
  mat <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/output/lulab/matrix/CPM-TMM_matrix_promoter.txt",sep = "\t",header = T,row.names = 1, check.names = F)
  samp <- read.table("/BioII/lulab_b/baopengfei/projects/multi-omics-explore/meta/lulab/multiomic-medip-passQCstringent-QIAGEN.txt",sep = "\t",header = F)$V1
  mat <- mat[,samp]
  #mat[1:4,1:4]
  mat <- log2(mat+1)
  #mat <- as.data.frame(t(scale(t(mat))))  # z-score
  mat$ENSG <- rownames(mat)
  mat$ENSG <- unlist(lapply(strsplit(rownames(mat),"|",fixed = T),function(x) x[1]))
  mat$ENSG <- gsub("promoter_","",mat$ENSG)
  
  ## join with expr mat
  l <- colnames(tmp)
  tmp[["log2CPM-TMM-promoter-passQCstringent-QIAGEN_matrix"]] <- ""
  tmp <- dplyr::left_join(tmp,mat)         
  
  ## add gini (only interpretable for non-negative quantities)
  #avoid minus value
  #tmp.exp <- 2^tmp[,grep("*-PKU-",colnames(tmp))]     # only for minus value
  tmp.exp <- tmp[,grep("*-PKU-",colnames(tmp))]
  
  tmp[["CRC_gini"]] <- apply(X = tmp.exp[,grep("CRC-PKU-",colnames(tmp.exp))], MARGIN = 1, FUN = function(x) edgeR::gini(as.numeric(x)))
  tmp[["NC_gini"]] <- apply(X = tmp.exp[,grep("NC-PKU-",colnames(tmp.exp))], MARGIN = 1, FUN = function(x) edgeR::gini(as.numeric(x)))
  tmp[["STAD_gini"]] <- apply(X = tmp.exp[,grep("STAD-PKU-",colnames(tmp.exp))], MARGIN = 1, FUN = function(x) edgeR::gini(as.numeric(x)))
  #hist(tmp$NC_gini)
  tmp <- tmp[,c(l,"CRC_gini","NC_gini","STAD_gini","log2CPM-TMM-promoter-passQCstringent-QIAGEN_matrix",samp)]
  
  #colnames(tmp)[1] <- "id"
  write.table(tmp,paste0("/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/output/lulab/diff-medip-TMM/promoter/",i,"-passQCstringent-QIAGEN-edger_exact/",i,"-passQCstringent-QIAGEN-edger_exact.all.txt"),quote = F,sep = "\t",col.names = T,row.names = F)
  rio::export(tmp,paste0("/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/output/lulab/diff-medip-TMM/promoter/",i,"-passQCstringent-QIAGEN-edger_exact/",i,"-passQCstringent-QIAGEN-edger_exact.all.xlsx"))   # top500Pvalue
}




# 
# mat <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/count_matrix_promoter150TSS50.txt",sep = "\t",header = T,row.names = 1,check.names = F)
# ratio <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/count_matrix_promoter150TSS50.txt.summary",sep = "\t",header = T,row.names = 1,check.names = F)
# colnames(ratio) <- gsub("output/lulab/bam-sorted-deduped/","",colnames(ratio),fixed = T)
# colnames(ratio) <- gsub(".bam","",colnames(ratio),fixed = T)
# ratio <- ratio[,colnames(mat)]
# 
# cor.test(as.numeric(ratio[1,]),as.numeric(colSums(mat)),method = "pearson")  # method = "spearman"
