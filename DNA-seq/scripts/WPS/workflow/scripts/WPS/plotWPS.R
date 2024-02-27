



setwd("/BioII/lulab_b/baopengfei/gitsoft/cfDNA-2016/")
read.shuffle.tx <- function(x,flank_len=2,region_len){
  #x <- "./results/intermediate/FTC_small_tx_eg/table/GRCh38/normalize/G4_eg--csFTA-10_1_STARTS_normalizedFlank.tsv.gz"
  #x <- files[1]
  data <- data.table::fread(x, data.table = F, stringsAsFactors = F, sep = "\t", header=FALSE,check.names = F)
  #data <- data.table::fread(x, data.table = F, stringsAsFactors = F, sep = ",", header=FALSE,check.names = F) # for not normalized value
  
  sample <- unlist(sapply(strsplit(x,"--",fixed = T),"[",2))
  sample <- unlist(sapply(strsplit(sample,"_",fixed = T),"[",1))
  
  #data[[sample]] <- matrixStats::rowMeans2(as.matrix(data[,flank_len:(flank_len+region_len*5-1)]),na.rm = T) # flank_len:(flank_len+region_len-1)
  #data[[sample]] <- matrixStats::rowMeans2(as.matrix(data[,2:(region_len+1)])) # flank_len:(flank_len+region_len-1)
  
  #data <- data[,c("V1",sample),drop=FALSE]
  data[["sample"]] <- sample
  return(data)
}


dst <- "FTC_small_tx_eg" # TCGA-LIHC_small, FTC_small_tx, FTC_long_tx
regions <- "AGO2_eg" # c("AGO2_eg","RBPhotspot_eg","G4iM_eg")
region_len <- 20 # long: 100, small: 20
classes <- c("WPS","WPS_v2","WPS_v3","COV","STARTS")
samples.table <- read.table(paste0("./config/",dst,"/samples.tsv"),sep = "\t",header = T,stringsAsFactors = F)
samples <- samples.table$sample

res2 <- list()
for (region in regions){
  #region <- "RBPhotspot"
  print(region)
  for (t in classes){
    # t <- "WPS_v2"
    print(t)
    files <- Sys.glob(paste0(paste0("./results/intermediate/",dst,"/table/GRCh38/normalize/",region,"--*_",t,"_normalizedFlank.tsv.gz")))  # CRC-PKU-10-wgs     _promoter150TSS50.bed  _20001000TSS.bed  _promoter300100exon1end.bed  _exon1end10002000.bed
    #files <- Sys.glob(paste0(paste0("./results/intermediate/",dst,"/table/GRCh38/target/",region,"--*_",t,".csv.gz")))  # CRC-PKU-10-wgs     _promoter150TSS50.bed  _20001000TSS.bed  _promoter300100exon1end.bed  _exon1end10002000.bed
    
    # run the anonymous function defined above
    wps.tss <- lapply(files, function(x) read.shuffle.tx(x,flank_len = 2, region_len = region_len))
    wps.tss <- do.call("rbind", wps.tss)
    wps.tss <- wps.tss[,!duplicated(colnames(wps.tss))]
    #wps.tss <- cbind(rownames(wps.tss),wps.tss)
    
    #wps.tss$V1 <- gsub("promoter150TSS50_|promoter300100exon1end_","",wps.tss$V1)
    #wps.tss$V1 <- unlist(lapply(strsplit(wps.tss$V1,".",fixed=T),function(x) x[1]))
    # wps.tss <- wps.tss[!duplicated(wps.tss$V1),]
    wps.tss <- wps.tss[!duplicated(wps.tss$sample),]
    #rownames(wps.tss) <- wps.tss$sample
    #colnames(wps.tss)[1] <- "gene_id"
    wps.tss$V1 <- NULL
    # wps <- cbind(gene_id=rownames(wps.tss),wps.tss)
    wps <- wps.tss
    wps <- na.omit(wps)
    #print(wps[1:3,1:3])
    wps$type <- t
    res2[[t]] <- wps
  }
}
res2.df <- do.call(rbind,res2)
#mat <- res2.df[,1:(ncol(res2.df)-2)]
# mat.scale <- t(scale(t(mat)))
# res2.df[,1:(ncol(res2.df)-2)] <- mat.scale
res2.df[1:3,1:3]

df <- reshape2::melt(data = res2.df, id.var=c("sample","type"))
df$variable <- as.integer(gsub("V","",df$variable))

df <- as_tibble(df) %>% 
  mutate(group=substr(sample,1,2)) %>% 
  dplyr::group_by(group,variable,type) %>%  
  dplyr::summarise(mean.value=mean(value,trim=0.1,na.rm=T),median.value=median(value,trim=0.1,na.rm=T))
df$group[df$group=="cs"] <- "cf_small"
df$group[df$group=="FT"] <- "EV_small"
ggplot(df,aes(x=variable,y=mean.value))+
  geom_line()+
  theme_bw() + 
  theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 16,color ="black",face="bold"), 
        axis.text = element_text(size= 20,color = "black"),
        #panel.grid=element_blank(),
        # panel.grid.major.x=element_blank(),
        # panel.grid.minor.x = element_blank(),
        # panel.grid.major.y = element_line(color = "black"), #size= 1,
        #panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(), # angle = 45, hjust = 1 
        legend.position = "right",
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16),
        strip.text = element_text(size= 16) )+
  facet_grid(type~group,scales = "free_y")

