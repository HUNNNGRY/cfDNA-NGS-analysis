#adapted from txRBP-EDA.R
#setwd("/BioII/lulab_b/baopengfei/gitsoft/cfDNA-2016/workflow/scripts/WPS")
options(scipen = 4,stringsAsFactors = F,digits = 5)
options(bedtools.path = "/BioII/lulab_b/baopengfei/biosoft")
library(bedtoolsr)
#

# test filter gn small domains by bp/length (AGO2,G4,RBPhotspot,Alu * WPS, lulab_WSQ-JYF,2401) ----------------------------------------
## get ref
options(stringsAsFactors = F,digits = 5) # scipen = 4,
#setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))

#ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F,header = F)
colnames(ref) <- c("transcript_id","tx.length")
head(ref)
#table(ref$transcript_type)


dst <- "lulab_WSQ-JYF" #"GSE71008_NCpool"
#dedup <- "all"
configDir <- paste0("/BioII/lulab_b/baopengfei/gitsoft/cfDNA-2016/config/",dst)
#smp <- "NCpool"
# dst <- "GSE71008"
# dedup <- "all"

peak.list <- list()
regions.df <- read.table(paste0(configDir,"/regions_original.tsv"),quote = "\"",header = T)
for(i in 1:nrow(regions.df)){
  peak.list[[regions.df$target[i]]] <- regions.df$path[i]
}
peak.list2 <- list()
for( method in names(peak.list)){
  # method <- "Alu_intron"
  print(method)
  m2 <- method # tolower(method)
  
  ## get peak 
  # peak <- Sys.glob(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_",dedup,"/",m2,"_by_sample/b5*_p*/intersect/",smp,".bed6")) # clam raw peak contain "," !!
  #peak <- Sys.glob(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_",dedup,"/",m2,"/b5*_p*_11RNA.bed")) # use consensus peak to aviod multi-smp-comb
  peak <- peak.list[[method]]
  
  print(peak)
  SRSF1 <- data.table::fread(peak,data.table = F,sep = "\t",header = F,stringsAsFactors = F)
  
  colnames(SRSF1) <- c("chr","start","end","name","score","strand")
  SRSF1$len <- SRSF1$end-SRSF1$start
  # hist(SRSF1$len,breaks = 1000,xlim = c(0,500))
  
  # SRSF1 <- SRSF1[SRSF1$len>=10 & SRSF1$len<=60,] # small !!!
  #SRSF1 <- SRSF1[SRSF1$len>=20 & SRSF1$len<=400,] # FTC_long !!!
  SRSF1 <- SRSF1[SRSF1$len>=10 & SRSF1$len<=350,] # Alu [0,350] !!!
  
  # # small: recenter & expand: 10 nt long from both side of central (total 20 nt); 50 nt each side 
  # peakLen <- 20 
  # flankLen <- 50 
  
  # Alu: recenter & expand: 10 nt long from both side of central (total 20 nt); 200 nt each side
  peakLen <- 20
  flankLen <- 200
  
  mid <- as.integer((SRSF1$start+SRSF1$end)*0.5)
  SRSF1$start <- mid - (peakLen+flankLen*2)*0.5 # fixed peak length: 20; fixed flank length: 50
  SRSF1$end <- SRSF1$start + (peakLen+flankLen*2)
  
  
  # # small: filter boundary peaks (wind protection: 16 nt)
  # protection <- 15
  
  # Alu: filter boundary peaks (wind protection: 20 nt)
  protection <- 20
  
  SRSF1$txLen <- ref$tx.length[match(SRSF1$chr,ref$transcript_id)]
  #SRSF1$RNA <- ref$transcript_type[match(SRSF1$chr,ref$transcript_id)]
  SRSF1 <- SRSF1[SRSF1$start>=(protection/2) & SRSF1$end+(protection/2)<=SRSF1$txLen,] 
  
  SRSF1 <- SRSF1[order(SRSF1$score,decreasing = T),]
  #SRSF1[1:3,]
  SRSF1$name <- SRSF1$name
  table(duplicated(SRSF1$name))
  SRSF1 <- SRSF1[!duplicated(SRSF1$name),]
  SRSF1$score <- 1
  SRSF1 <- SRSF1[, c("chr","start","end","name","score","strand")]
  
  #optional/TODO: rm blackRegion-overlapped regions
  #bedtools intersect -v -a {input.Region_file} -b blackList.bed
  
  rbp.df <- SRSF1
  #rbp.df2 <- rbp.df[!grepl("chr",rbp.df$chr,fixed = T),] # rm gn chr
  newPath <- paste0(configDir,"/",method,".bed")
  peak.list2[[method]] <- newPath
  data.table::fwrite(rbp.df,newPath,quote = F,sep = "\t",col.names = F,row.names = F)
}
df <- as.data.frame(t(as.data.frame(peak.list2)))
colnames(df) <- "path"
df$target <- rownames(df)
df$path <- paste0("\"",df$path,"\"")
data.table::fwrite(df[,c("target","path")],paste0(configDir,"/regions.tsv"),quote = F,sep = "\t",col.names = T,row.names = F)
#





# test filter gn small domains by folds (AGO2,G4,RBPhotspot,Alu * WPS, lulab_WSQ-JYF,2401, newest) ----------------------------------------
## get ref
options(stringsAsFactors = F,digits = 5) # scipen = 4,
#setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
setwd("/BioII/lulab_b/baopengfei/gitsoft/cfDNA-2016")


#ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F)
ref <- data.table::fread("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID",data.table = F,sep = '\t',check.names = F,stringsAsFactors = F,header = F)
colnames(ref) <- c("transcript_id","tx.length")
head(ref)
#table(ref$transcript_type)


#dst <- "lulab_WSQ-JYF" #"GSE71008_NCpool"
#dedup <- "all"
#configDir <- paste0("/BioII/lulab_b/baopengfei/gitsoft/cfDNA-2016/config/",dst)
configDir <- paste0("/BioII/lulab_b/baopengfei/gitsoft/cfDNA-2016/","regions")
configDirOut <- paste0("/BioII/lulab_b/baopengfei/gitsoft/cfDNA-2016/","regions_foldExt")
dir.create(configDirOut)
#smp <- "NCpool"
# dst <- "GSE71008"
# dedup <- "all"

peak.list <- list()
regions.df <- read.table(paste0(configDir,"/regions_original.tsv"),quote = "\"",header = T)
for(i in 1:nrow(regions.df)){
  peak.list[[regions.df$target[i]]] <- regions.df$path[i]
}
peak.list2 <- list()
for( method in names(peak.list)){
  # method <- "Alu_intron"
  print(method)
  m2 <- method # tolower(method)
  
  ## get peak 
  # peak <- Sys.glob(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_",dedup,"/",m2,"_by_sample/b5*_p*/intersect/",smp,".bed6")) # clam raw peak contain "," !!
  #peak <- Sys.glob(paste0("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/",dst,"/call_peak_",dedup,"/",m2,"/b5*_p*_11RNA.bed")) # use consensus peak to aviod multi-smp-comb
  peak <- peak.list[[method]]
  
  print(peak)
  SRSF1 <- data.table::fread(peak,data.table = F,sep = "\t",header = F,stringsAsFactors = F)
  
  colnames(SRSF1) <- c("chr","start","end","name","score","strand")
  SRSF1$len <- SRSF1$end-SRSF1$start
  # hist(SRSF1$len,breaks = 1000,xlim = c(0,500))
  
  # SRSF1 <- SRSF1[SRSF1$len>=10 & SRSF1$len<=60,] # small !!!
  #SRSF1 <- SRSF1[SRSF1$len>=20 & SRSF1$len<=400,] # FTC_long !!!
  SRSF1 <- SRSF1[SRSF1$len>=10 & SRSF1$len<=350,] # Alu [0,350] !!!
  
  # # small: recenter & expand: 10 nt long from both side of central (total 20 nt); 50 nt each side 
  # peakLen <- 20 
  # flankLen <- 50 
  
  # Alu: recenter & expand: 10 nt long from both side of central (total 20 nt); 200 nt each side
  # peakLen <- 20
  # flankLen <- 200
  foldExt <- 1 # extend to fold of region distance from each boundary 
  # mid <- as.integer((SRSF1$start+SRSF1$end)*0.5)
  # SRSF1$start <- mid - (peakLen+flankLen*2)*0.5 # fixed peak length: 20; fixed flank length: 50
  # SRSF1$end <- SRSF1$start + (peakLen+flankLen*2)
  # 
  # 
  # # # small: filter boundary peaks (wind protection: 16 nt)
  # # protection <- 15
  # 
  # # Alu: filter boundary peaks (wind protection: 20 nt)
  # protection <- 20
  # 
  # SRSF1$txLen <- ref$tx.length[match(SRSF1$chr,ref$transcript_id)]
  # #SRSF1$RNA <- ref$transcript_type[match(SRSF1$chr,ref$transcript_id)]
  # SRSF1 <- SRSF1[SRSF1$start>=(protection/2) & SRSF1$end+(protection/2)<=SRSF1$txLen,] 
  #SRSF1[1:3,]
  SRSF1.wid <- SRSF1$end - SRSF1$start
  SRSF1 <- bedtoolsr::bt.slop(i = SRSF1, s = T, g = ref, l = foldExt, r = foldExt, pct = T)
  colnames(SRSF1) <- c("chr","start","end","name","score","strand")
  
  SRSF1.wid2 <- SRSF1$end - SRSF1$start
  table(SRSF1.wid2==SRSF1.wid*(foldExt*2+1))
  SRSF1 <- SRSF1[SRSF1.wid2==SRSF1.wid*(foldExt*2+1),]
  
  SRSF1 <- SRSF1[order(SRSF1$score,decreasing = T),]
  #SRSF1[1:3,]
  SRSF1$name <- SRSF1$name
  table(duplicated(SRSF1$name))
  SRSF1 <- SRSF1[!duplicated(SRSF1$name),]
  SRSF1$score <- 1
  SRSF1 <- SRSF1[, c("chr","start","end","name","score","strand")]
  
  #optional/TODO: rm blackRegion-overlapped regions
  #bedtools intersect -v -a {input.Region_file} -b blackList.bed
  
  rbp.df <- SRSF1
  #rbp.df2 <- rbp.df[!grepl("chr",rbp.df$chr,fixed = T),] # rm gn chr
  newPath <- paste0(configDirOut,"/",method,".bed")
  peak.list2[[method]] <- newPath
  data.table::fwrite(rbp.df,newPath,quote = F,sep = "\t",col.names = F,row.names = F)
}
df <- as.data.frame(t(as.data.frame(peak.list2)))
colnames(df) <- "path"
df$target <- rownames(df)
df$path <- paste0("\"",df$path,"\"")
data.table::fwrite(df[,c("target","path")],paste0(configDirOut,"/regions.tsv"),quote = F,sep = "\t",col.names = T,row.names = F)
#
