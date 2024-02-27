# prepare cfMeDIP mat
# last 220107 by bpf


# correct promoter&genebody CPM correctWGS -------------------------------------------
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/")
options(stringsAsFactors = F)

for (region in c("promoter","gene")){
print(region)
pro <- read.table(paste0("output/lulab/matrix/CPM-TMM_matrix_",region,".txt"),sep = "\t",header = T,check.names = F,row.names = 1)
wgs <- read.table(paste0("../DNA-seq/output/lulab/matrix/CPM-TMM_matrix_",region,".txt"),sep = "\t",header = T,check.names = F,row.names = 1)
colnames(wgs) <- gsub("wgs","me",colnames(wgs))
pro <- pro[, colnames(pro) %in% colnames(wgs)]
table(colnames(pro) %in% colnames(wgs))
wgs <- wgs[,colnames(pro)]
table(colnames(wgs) %in% colnames(pro))  # should be all true
table(rownames(wgs) %in% rownames(pro)) # should be all true

pro.correctWGS <- pro/(wgs+1)  # fold enrichment socre, like MACS2 FE mode ?

hist(as.vector(as.matrix(pro.correctWGS)),breaks = 9000,xlim = c(0,10))
table(as.matrix(pro.correctWGS)==0)
pro.correctWGS.out <- data.frame("gene_id"=rownames(pro.correctWGS),pro.correctWGS,check.names = F)
write.table(pro.correctWGS.out,paste0("output/lulab/matrix/CPM-TMM_matrix_",region,".correctWGS.txt"),sep = "\t",quote = F,col.names = T,row.names = F)
# FALSE    TRUE 
# 4510874  443606 
#10% zeros
}




# correct promoter&genebody CPM correctGC&correctWGS -------------------------------------------
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/")
options(stringsAsFactors = F)

for (region in c("promoter","gene")){
  print(region)
  pro <- read.table(paste0("output/lulab/matrix/CPM-TMM_matrix_",region,".correctGC.txt"),sep = "\t",header = T,check.names = F,row.names = 1)
  wgs <- read.table(paste0("../DNA-seq/output/lulab/matrix/CPM-TMM_matrix_",region,".correctGC.txt"),sep = "\t",header = T,check.names = F,row.names = 1)
  colnames(wgs) <- gsub("wgs","me",colnames(wgs))
  pro <- pro[, colnames(pro) %in% colnames(wgs)]
  table(colnames(pro) %in% colnames(wgs))
  wgs <- wgs[,colnames(pro)]
  table(colnames(wgs) %in% colnames(pro))  # should be all true
  table(rownames(wgs) %in% rownames(pro)) # should be all true
  
  pro.correctWGS <- pro/(wgs+1)  # fold enrichment socre, like MACS2 FE mode ?
  
  hist(as.vector(as.matrix(pro.correctWGS)),breaks = 9000,xlim = c(0,10))
  table(as.matrix(pro.correctWGS)==0)
  pro.correctWGS.out <- data.frame("gene_id"=rownames(pro.correctWGS),pro.correctWGS,check.names = F)
  write.table(pro.correctWGS.out,paste0("output/lulab/matrix/CPM-TMM_matrix_",region,".correctGC.correctWGS.txt"),sep = "\t",quote = F,col.names = T,row.names = F)
  # FALSE    TRUE 
  # 4510874  443606 
  #10% zeros
}
