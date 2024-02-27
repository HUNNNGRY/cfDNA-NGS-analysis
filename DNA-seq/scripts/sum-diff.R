# sum diff table for candidates selection
# last 210103 by bpf


# sum NOV-countSumCPM -----------------------------------------------------
cmp <- c("STADvsNC","CRCvsNC","CRC_STADvsNC")
for (i in cmp){
tmp <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/diff-wgs-TMM/NOVsum/",i,"-passQCstringent-edger_exact/",i,"-passQCstringent-edger_exact.txt"),sep = "\t",header = T)
tmp <- tmp[order(-tmp$pvalue,abs(tmp$log2FoldChange),decreasing = T),]
tmp$ENSG <- unlist(lapply(strsplit(tmp$id,"_"),function(x) x[2]))
tmp$ENSG <- unlist(lapply(strsplit(tmp$ENSG,"|",fixed = T),function(x) x[1]))

## only keep top500
#tmp <- tmp[1:500,]

ref <- read.table("/BioII/lulab_b/baopengfei/shared_reference/hg38/tss-exon1.gtf",sep = "\t",header = T)
#ref[1:3,]
ref$TSS.start <- ifelse(ref$strand=="+",ref$tss-150,ref$tss-50)
ref$TSS.end <- ifelse(ref$strand=="+",ref$tss+50,ref$tss+150)
ref$TSS.nucleosome.region <- paste0(ref$chr,":",ref$TSS.start,"-",ref$TSS.end)

ref$Exon1end.start <- ifelse(ref$strand=="+",ref$exon1end-300,ref$exon1end+100)
ref$Exon1end.end <- ifelse(ref$strand=="+",ref$exon1end-100,ref$exon1end+300)
ref$Exon1end.nucleosome.region <- paste0(ref$chr,":",ref$Exon1end.start,"-",ref$Exon1end.end)
#table(ref$strand)

colnames(ref)
ref <- ref[,c("chr","strand","ensg","TSS.nucleosome.region","Exon1end.nucleosome.region")]
colnames(ref)[3] <- "ENSG"

tmp <- dplyr::left_join(tmp,ref)         


## expr mat 
mat <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/CPM_matrix_NOVsum.txt",sep = "\t",header = T,check.names = F)
samp <- read.table("/BioII/lulab_b/baopengfei/projects/multi-omics-explore/meta/lulab/multiomic-wgs-passQCstringent.txt",sep = "\t",header = F)$V1
mat <- mat[,samp]
#mat[1:4,1:4]
mat <- log2(mat+1)
#colnames(tmp)[3] <- "log2CPM.Mean"
mat$id <- rownames(mat)


## join with expr mat
l <- colnames(tmp)
tmp$log2CPM.mat <- ""
tmp <- dplyr::left_join(tmp,mat)         

## add gini (only interpretable for non-negative quantities)
tmp$CRC.gini <- apply(X = tmp[,grep("CRC-PKU-",colnames(tmp))], MARGIN = 1, FUN = function(x) edgeR::gini(as.numeric(x)))
tmp$NC.gini <- apply(X = tmp[,grep("NC-PKU-",colnames(tmp))], MARGIN = 1, FUN = function(x) edgeR::gini(as.numeric(x)))
tmp$STAD.gini <- apply(X = tmp[,grep("STAD-PKU-",colnames(tmp))], MARGIN = 1, FUN = function(x) edgeR::gini(as.numeric(x)))
#hist(tmp$NC.gini)
tmp <- tmp[,c(l,"CRC.gini","NC.gini","STAD.gini","log2CPM.mat",samp)]

write.table(tmp,paste0("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/diff-wgs-TMM/NOVsum/",i,"-passQCstringent-edger_exact/",i,"-passQCstringent-edger_exact.top500Pvalue.txt"),quote = F,sep = "\t",col.names = T,row.names = F)
rio::export(tmp,paste0("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/diff-wgs-TMM/NOVsum/",i,"-passQCstringent-edger_exact/",i,"-passQCstringent-edger_exact.top500Pvalue.xlsx"))
}




# sum NOV-relDepRatio (old)  -----------------------------------------------------
cmp <- c("STADvsNC","CRCvsNC","CRC_STADvsNC")
for (i in cmp){
#i <- "STADvsNC"
tmp <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/diff-wgs-TMM/NucOccuVarRelDepRatioMean/",i,"-passQCstringent-wilcox/",i,"-passQCstringent-wilcox.txt"),sep = "\t",header = T)
tmp <- tmp[order(-tmp$pvalue,abs(tmp$log2FoldChange),decreasing = T),]
#tmp$ENSG <- unlist(lapply(strsplit(tmp$id,"_"),function(x) x[2]))
#tmp$ENSG <- unlist(lapply(strsplit(tmp$ENSG,"|",fixed = T),function(x) x[1]))
colnames(tmp)[1] <- "ENSG"

## only keep top500
#tmp <- tmp[1:500,]

ref <- read.table("/BioII/lulab_b/baopengfei/shared_reference/hg38/tss-exon1.gtf",sep = "\t",header = T)
#ref[1:3,]
ref$TSS.start <- ifelse(ref$strand=="+",ref$tss-150,ref$tss-50)
ref$TSS.end <- ifelse(ref$strand=="+",ref$tss+50,ref$tss+150)
ref$TSS.nucleosome.region <- paste0(ref$chr,":",ref$TSS.start,"-",ref$TSS.end)

ref$Exon1end.start <- ifelse(ref$strand=="+",ref$exon1end-300,ref$exon1end+100)
ref$Exon1end.end <- ifelse(ref$strand=="+",ref$exon1end-100,ref$exon1end+300)
ref$Exon1end.nucleosome.region <- paste0(ref$chr,":",ref$Exon1end.start,"-",ref$Exon1end.end)
#table(ref$strand)

colnames(ref)
ref <- ref[,c("chr","strand","ensg","TSS.nucleosome.region","Exon1end.nucleosome.region")]
colnames(ref)[3] <- "ENSG"

tmp <- dplyr::left_join(tmp,ref)         


## expr mat 
mat <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/NucOccuVarRelDepRatio_matrix_aver.txt",sep = "\t",header = T,row.names = 1, check.names = F)
samp <- read.table("/BioII/lulab_b/baopengfei/projects/multi-omics-explore/meta/lulab/multiomic-wgs-passQCstringent.txt",sep = "\t",header = F)$V1
mat <- mat[,samp]
#mat[1:4,1:4]
#mat <- log2(mat+1)
#colnames(tmp)[3] <- "log2CPM.Mean"
mat$ENSG <- rownames(mat)


## join with expr mat
l <- colnames(tmp)
tmp$log2CPM.mat <- ""
tmp <- dplyr::left_join(tmp,mat)         

## add gini (only interpretable for non-negative quantities)
tmp$CRC.gini <- apply(X = tmp[,grep("CRC-PKU-",colnames(tmp))], MARGIN = 1, FUN = function(x) edgeR::gini(as.numeric(x)))
tmp$NC.gini <- apply(X = tmp[,grep("NC-PKU-",colnames(tmp))], MARGIN = 1, FUN = function(x) edgeR::gini(as.numeric(x)))
tmp$STAD.gini <- apply(X = tmp[,grep("STAD-PKU-",colnames(tmp))], MARGIN = 1, FUN = function(x) edgeR::gini(as.numeric(x)))
#hist(tmp$NC.gini)
tmp <- tmp[,c(l,"CRC.gini","NC.gini","STAD.gini","log2CPM.mat",samp)]

colnames(tmp)[1] <- "id"
write.table(tmp,paste0("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/diff-wgs-TMM/NucOccuVarRelDepRatioMean/",i,"-passQCstringent-wilcox/",i,"-passQCstringent-wilcox.top500Pvalue.txt"),quote = F,sep = "\t",col.names = T,row.names = F)
rio::export(tmp,paste0("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/diff-wgs-TMM/NucOccuVarRelDepRatioMean/",i,"-passQCstringent-wilcox/",i,"-passQCstringent-wilcox.top500Pvalue.xlsx"))

}




# sum NOV-WPS -----------------------------------------------------
cmp <- c("STADvsNC","CRCvsNC","CRC_STADvsNC")
for (i in cmp){
#i <- "STADvsNC"
tmp <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/diff-wgs-TMM/NucOccuVarWPS/",i,"-passQCstringent-wilcox/",i,"-passQCstringent-wilcox.txt"),sep = "\t",header = T)
tmp <- tmp[order(-tmp$pvalue,abs(tmp$log2FoldChange),decreasing = T),]
#tmp$ENSG <- unlist(lapply(strsplit(tmp$id,"_"),function(x) x[2]))
#tmp$ENSG <- unlist(lapply(strsplit(tmp$ENSG,"|",fixed = T),function(x) x[1]))
colnames(tmp)[1] <- "ENSG"

## only keep top500
#tmp <- tmp[1:500,]

ref <- read.table("/BioII/lulab_b/baopengfei/shared_reference/hg38/tss-exon1.gtf",sep = "\t",header = T)
#ref[1:3,]
ref$TSS.start <- ifelse(ref$strand=="+",ref$tss-150,ref$tss-50)
ref$TSS.end <- ifelse(ref$strand=="+",ref$tss+50,ref$tss+150)
ref$TSS.nucleosome.region <- paste0(ref$chr,":",ref$TSS.start,"-",ref$TSS.end)

ref$Exon1end.start <- ifelse(ref$strand=="+",ref$exon1end-300,ref$exon1end+100)
ref$Exon1end.end <- ifelse(ref$strand=="+",ref$exon1end-100,ref$exon1end+300)
ref$Exon1end.nucleosome.region <- paste0(ref$chr,":",ref$Exon1end.start,"-",ref$Exon1end.end)
#table(ref$strand)

colnames(ref)
ref <- ref[,c("chr","strand","ensg","TSS.nucleosome.region","Exon1end.nucleosome.region")]
colnames(ref)[3] <- "ENSG"
#ref$ENSG <- unlist(lapply(strsplit(ref$ENSG,".",fixed=T),function(x) x[1]))
ref <- ref[!duplicated(ref$ENSG),]

tmp <- dplyr::left_join(tmp,ref)         


## expr mat 
mat <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/WPS_matrix_gene-aver.txt",sep = "\t",header = T,row.names = 1, check.names = F)
samp <- read.table("/BioII/lulab_b/baopengfei/projects/multi-omics-explore/meta/lulab/multiomic-wgs-passQCstringent.txt",sep = "\t",header = F)$V1
mat <- mat[,samp]
#mat[1:4,1:4]
#mat <- log2(mat+1)
#colnames(tmp)[3] <- "log2CPM.Mean"
mat$ENSG <- rownames(mat)


## join with expr mat
l <- colnames(tmp)
tmp$log2CPM.mat <- ""
tmp <- dplyr::left_join(tmp,mat)         

## add gini (only interpretable for non-negative quantities)
#avoid minus value
tmp.exp <- 2^tmp[,grep("*-PKU-",colnames(tmp))]

tmp$CRC.gini <- apply(X = tmp.exp[,grep("CRC-PKU-",colnames(tmp.exp))], MARGIN = 1, FUN = function(x) edgeR::gini(as.numeric(x)))
tmp$NC.gini <- apply(X = tmp.exp[,grep("NC-PKU-",colnames(tmp.exp))], MARGIN = 1, FUN = function(x) edgeR::gini(as.numeric(x)))
tmp$STAD.gini <- apply(X = tmp.exp[,grep("STAD-PKU-",colnames(tmp.exp))], MARGIN = 1, FUN = function(x) edgeR::gini(as.numeric(x)))
#hist(tmp$NC.gini)
tmp <- tmp[,c(l,"CRC.gini","NC.gini","STAD.gini","log2CPM.mat",samp)]

colnames(tmp)[1] <- "id"
write.table(tmp,paste0("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/diff-wgs-TMM/NucOccuVarWPS/",i,"-passQCstringent-wilcox/",i,"-passQCstringent-wilcox.top500Pvalue.txt"),quote = F,sep = "\t",col.names = T,row.names = F)
rio::export(tmp,paste0("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/diff-wgs-TMM/NucOccuVarWPS/",i,"-passQCstringent-wilcox/",i,"-passQCstringent-wilcox.top500Pvalue.xlsx"))

}







# sum NOV-TSS CPM (2201115) -----------------------------------------------------
cmp <- c("STADvsNC","CRCvsNC","CRC_STADvsNC")
for (r in c("promoter150TSS50","promoter300100exon1end")){
for (i in cmp){
  #i <- "STADvsNC"
  #r <- "promoter300100exon1end"
  tmp <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/diff-wgs-TMM/",r,"/",i,"-passQCstringent-edger_exact/",i,"-passQCstringent-edger_exact.txt"),sep = "\t",header = T)
  tmp <- tmp[order(-tmp$pvalue,abs(tmp$log2FoldChange),decreasing = T),]
  tmp$ENSG <- unlist(lapply(strsplit(tmp$id,"_"),function(x) x[2]))
  tmp$ENSG <- unlist(lapply(strsplit(tmp$ENSG,"|",fixed = T),function(x) x[1]))
  colnames(tmp)[3] <- "log2CPM_Mean"
  
  ## only keep top500
  #tmp <- tmp[1:500,]
  
  ref <- read.table("/BioII/lulab_b/baopengfei/shared_reference/hg38/tss-exon1.gtf",sep = "\t",header = T)
  #ref[1:3,]
  ref$TSS.start <- ifelse(ref$strand=="+",ref$tss-150,ref$tss-50)
  ref$TSS.end <- ifelse(ref$strand=="+",ref$tss+50,ref$tss+150)
  ref$TSS.nucleosome.region <- paste0(ref$chr,":",ref$TSS.start,"-",ref$TSS.end)
  
  ref$Exon1end.start <- ifelse(ref$strand=="+",ref$exon1end-300,ref$exon1end+100)
  ref$Exon1end.end <- ifelse(ref$strand=="+",ref$exon1end-100,ref$exon1end+300)
  ref$Exon1end.nucleosome.region <- paste0(ref$chr,":",ref$Exon1end.start,"-",ref$Exon1end.end)
  #table(ref$strand)
  
  colnames(ref)
  ref <- ref[,c("chr","strand","ensg","TSS.nucleosome.region","Exon1end.nucleosome.region")]
  colnames(ref)[3] <- "ENSG"
  
  tmp <- dplyr::left_join(tmp,ref)         
  
  
  ## expr mat 
  mat <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/CPM-TMM_matrix_",r,".txt"),sep = "\t",header = T,check.names = F,row.names = 1)
  samp <- read.table("/BioII/lulab_b/baopengfei/projects/multi-omics-explore/meta/lulab/multiomic-wgs-passQCstringent.txt",sep = "\t",header = F)$V1
  mat <- mat[,samp]
  #mat[1:4,1:4]
  mat <- log2(mat+1)
  mat$ENSG <- rownames(mat)
  mat$ENSG <- gsub("promoter150TSS50_|promoter300100exon1end_","",mat$ENSG)
  mat$ENSG <- gsub("|201","",mat$ENSG,fixed = T)
  
  
  ## join with expr mat
  l <- colnames(tmp)
  tmp[["log2CPM-TMM-passQCstringent_matrix"]] <- ""
  tmp <- dplyr::left_join(tmp,mat)         
  
  ## add gini (only interpretable for non-negative quantities)
  tmp[["CRC_gini"]] <- apply(X = tmp[,grep("CRC-PKU-",colnames(tmp))], MARGIN = 1, FUN = function(x) edgeR::gini(as.numeric(x)))
  tmp[["NC_gini"]] <- apply(X = tmp[,grep("NC-PKU-",colnames(tmp))], MARGIN = 1, FUN = function(x) edgeR::gini(as.numeric(x)))
  tmp[["STAD_gini"]] <- apply(X = tmp[,grep("STAD-PKU-",colnames(tmp))], MARGIN = 1, FUN = function(x) edgeR::gini(as.numeric(x)))
  #hist(tmp$NC.gini)
  tmp <- tmp[,c(l,"CRC_gini","NC_gini","STAD_gini","log2CPM-TMM-passQCstringent_matrix",samp)]
  
  write.table(tmp,paste0("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/diff-wgs-TMM/",r,"/",i,"-passQCstringent-edger_exact/",i,"-passQCstringent-edger_exact.all.txt"),quote = F,sep = "\t",col.names = T,row.names = F)
  rio::export(tmp,paste0("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/diff-wgs-TMM/",r,"/",i,"-passQCstringent-edger_exact/",i,"-passQCstringent-edger_exact.all.xlsx"))
}
}





# sum NOV-TSS CPM (filter QIAGEN kit, 2201117) -----------------------------------------------------
cmp <- c("STADvsNC","CRC_STADvsNC")  # CRCvsNC
for (r in c("promoter150TSS50","promoter300100exon1end")){
  for (i in cmp){
    #i <- "STADvsNC"
    #r <- "promoter300100exon1end"
    tmp <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/diff-wgs-TMM/",r,"/",i,"-passQCstringent-QIAGEN-edger_exact/",i,"-passQCstringent-QIAGEN-edger_exact.txt"),sep = "\t",header = T)
    tmp <- tmp[order(-tmp$pvalue,abs(tmp$log2FoldChange),decreasing = T),]
    tmp$ENSG <- unlist(lapply(strsplit(tmp$id,"_"),function(x) x[2]))
    tmp$ENSG <- unlist(lapply(strsplit(tmp$ENSG,"|",fixed = T),function(x) x[1]))
    colnames(tmp)[3] <- "log2CPM_Mean"
    
    ## only keep top500
    #tmp <- tmp[1:500,]
    
    ref <- read.table("/BioII/lulab_b/baopengfei/shared_reference/hg38/tss-exon1.gtf",sep = "\t",header = T)
    #ref[1:3,]
    ref$TSS.start <- ifelse(ref$strand=="+",ref$tss-150,ref$tss-50)
    ref$TSS.end <- ifelse(ref$strand=="+",ref$tss+50,ref$tss+150)
    ref$TSS.nucleosome.region <- paste0(ref$chr,":",ref$TSS.start,"-",ref$TSS.end)
    
    ref$Exon1end.start <- ifelse(ref$strand=="+",ref$exon1end-300,ref$exon1end+100)
    ref$Exon1end.end <- ifelse(ref$strand=="+",ref$exon1end-100,ref$exon1end+300)
    ref$Exon1end.nucleosome.region <- paste0(ref$chr,":",ref$Exon1end.start,"-",ref$Exon1end.end)
    #table(ref$strand)
    
    colnames(ref)
    ref <- ref[,c("chr","strand","ensg","TSS.nucleosome.region","Exon1end.nucleosome.region")]
    colnames(ref)[3] <- "ENSG"
    
    tmp <- dplyr::left_join(tmp,ref)         
    
    
    ## expr mat 
    mat <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/CPM-TMM_matrix_",r,".txt"),sep = "\t",header = T,check.names = F,row.names = 1)
    samp <- read.table("/BioII/lulab_b/baopengfei/projects/multi-omics-explore/meta/lulab/multiomic-wgs-passQCstringent-QIAGEN.txt",sep = "\t",header = F)$V1
    mat <- mat[,samp]
    #mat[1:4,1:4]
    mat <- log2(mat+1)
    #mat <- as.data.frame(t(scale(t(mat))))  # z-score
    mat$ENSG <- rownames(mat)
    mat$ENSG <- gsub("promoter150TSS50_|promoter300100exon1end_","",mat$ENSG)
    mat$ENSG <- gsub("|201","",mat$ENSG,fixed = T)
    
    
    ## join with expr mat
    l <- colnames(tmp)
    tmp[["log2CPM-TMM-passQCstringent-QIAGEN_matrix"]] <- ""
    tmp <- dplyr::left_join(tmp,mat)         
    
    ## add gini (only interpretable for non-negative quantities)
    tmp[["CRC_gini"]] <- apply(X = tmp[,grep("CRC-PKU-",colnames(tmp))], MARGIN = 1, FUN = function(x) edgeR::gini(as.numeric(x)))
    tmp[["NC_gini"]] <- apply(X = tmp[,grep("NC-PKU-",colnames(tmp))], MARGIN = 1, FUN = function(x) edgeR::gini(as.numeric(x)))
    tmp[["STAD_gini"]] <- apply(X = tmp[,grep("STAD-PKU-",colnames(tmp))], MARGIN = 1, FUN = function(x) edgeR::gini(as.numeric(x)))
    #hist(tmp$NC.gini)
    tmp <- tmp[,c(l,"CRC_gini","NC_gini","STAD_gini","log2CPM-TMM-passQCstringent-QIAGEN_matrix",samp)]
    
    write.table(tmp,paste0("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/diff-wgs-TMM/",r,"/",i,"-passQCstringent-QIAGEN-edger_exact/",i,"-passQCstringent-QIAGEN-edger_exact.all.txt"),quote = F,sep = "\t",col.names = T,row.names = F)
    rio::export(tmp,paste0("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/diff-wgs-TMM/",r,"/",i,"-passQCstringent-QIAGEN-edger_exact/",i,"-passQCstringent-QIAGEN-edger_exact.all.xlsx"))
  }
}


# sum NOV-relDepRatio (filter QIAGEN kit, 2201117) -----------------------------------------------------
cmp <- c("STADvsNC","CRC_STADvsNC")  # CRCvsNC
for (r in c("promoter150TSS50","promoter300100exon1end")){
  for (i in cmp){
    #i <- "STADvsNC"
    #r <- "promoter300100exon1end"
    tmp <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/diff-wgs-NucOccuVarRelDepRatio/",r,"/",i,"-passQCstringent-QIAGEN-wilcox/",i,"-passQCstringent-QIAGEN-wilcox.txt"),sep = "\t",header = T)
    tmp <- tmp[order(-tmp$pvalue,abs(tmp$log2FoldChange),decreasing = T),]
    tmp$ENSG <- unlist(lapply(strsplit(tmp$id,"_"),function(x) x[2]))
    tmp$ENSG <- unlist(lapply(strsplit(tmp$ENSG,"|",fixed = T),function(x) x[1]))
    #colnames(tmp)[3] <- "RelativeDepthRatio"  # no such col for wilcox
    
    ## only keep top500
    #tmp <- tmp[1:500,]
    
    ref <- read.table("/BioII/lulab_b/baopengfei/shared_reference/hg38/tss-exon1.gtf",sep = "\t",header = T)
    #ref[1:3,]
    ref$TSS.start <- ifelse(ref$strand=="+",ref$tss-150,ref$tss-50)
    ref$TSS.end <- ifelse(ref$strand=="+",ref$tss+50,ref$tss+150)
    ref$TSS.nucleosome.region <- paste0(ref$chr,":",ref$TSS.start,"-",ref$TSS.end)
    
    ref$Exon1end.start <- ifelse(ref$strand=="+",ref$exon1end-300,ref$exon1end+100)
    ref$Exon1end.end <- ifelse(ref$strand=="+",ref$exon1end-100,ref$exon1end+300)
    ref$Exon1end.nucleosome.region <- paste0(ref$chr,":",ref$Exon1end.start,"-",ref$Exon1end.end)
    #table(ref$strand)
    
    colnames(ref)
    ref <- ref[,c("chr","strand","ensg","TSS.nucleosome.region","Exon1end.nucleosome.region")]
    colnames(ref)[3] <- "ENSG"
    
    tmp <- dplyr::left_join(tmp,ref)         
    
    
    ## expr mat 
    mat <- read.table(paste0("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/NucOccuVarRelDepRatio_matrix_",r,".txt"),sep = "\t",header = T,check.names = F,row.names = 1)
    samp <- read.table("/BioII/lulab_b/baopengfei/projects/multi-omics-explore/meta/lulab/multiomic-wgs-passQCstringent-QIAGEN.txt",sep = "\t",header = F)$V1
    mat <- mat[,samp]
    #mat[1:4,1:4]
    #mat <- log2(mat+1)
    #mat <- as.data.frame(t(scale(t(mat))))  # z-score
    mat$ENSG <- rownames(mat)
    mat$ENSG <- gsub("promoter150TSS50_|promoter300100exon1end_","",mat$ENSG)
    mat$ENSG <- gsub("|201","",mat$ENSG,fixed = T)
    
    
    ## join with expr mat
    l <- colnames(tmp)
    tmp[["NucOccuVarRelDepRatio-passQCstringent-QIAGEN_matrix"]] <- ""
    tmp <- dplyr::left_join(tmp,mat)         
    
    ## add gini (only interpretable for non-negative quantities)
    tmp[["CRC_gini"]] <- apply(X = tmp[,grep("CRC-PKU-",colnames(tmp))], MARGIN = 1, FUN = function(x) edgeR::gini(as.numeric(x)))
    tmp[["NC_gini"]] <- apply(X = tmp[,grep("NC-PKU-",colnames(tmp))], MARGIN = 1, FUN = function(x) edgeR::gini(as.numeric(x)))
    tmp[["STAD_gini"]] <- apply(X = tmp[,grep("STAD-PKU-",colnames(tmp))], MARGIN = 1, FUN = function(x) edgeR::gini(as.numeric(x)))
    #hist(tmp$NC.gini)
    tmp <- tmp[,c(l,"CRC_gini","NC_gini","STAD_gini","NucOccuVarRelDepRatio-passQCstringent-QIAGEN_matrix",samp)]
    
    write.table(tmp,paste0("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/diff-wgs-NucOccuVarRelDepRatio/",r,"/",i,"-passQCstringent-QIAGEN-wilcox/",i,"-passQCstringent-QIAGEN-wilcox.all.txt"),quote = F,sep = "\t",col.names = T,row.names = F)
    rio::export(tmp,paste0("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/diff-wgs-NucOccuVarRelDepRatio/",r,"/",i,"-passQCstringent-QIAGEN-wilcox/",i,"-passQCstringent-QIAGEN-wilcox.all.xlsx"))
  }
}

