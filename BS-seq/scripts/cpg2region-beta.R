#!/usr/bin/R

setwd("/BioII/lulab_b/baopengfei/projects/exOmics/BS-seq")
dataset <- "GSE149438"
gtf <- read.table("genome/gtf/gencodev27_onlypromoter.gtf",stringsAsFactors=F) 
sample <- read.table(paste0("metadata/",dataset,"/sample_ids_test.txt"),stringsAsFactors=F)
region <- matrix(data=NA,nrow=length(gtf[,10]),ncol=length(sample$V1),dimnames=list(c(gtf[,10]),c(sample$V1)))


for (id in sample$V1)
{
print(paste0("start:",id))
cpg <- read.table(paste0("output/",dataset,"/",id,".deduplicated.bismark.cov.gz"),stringsAsFactors=F)
	#mapply(region,cpg2region)
	for (i in 1:nrow(region))
	{
		print(i)
		list <- rownames(cpg)[cpg[,2] >= gtf[i,4] & cpg[,2] <= gtf[i,5]]
		region[i,id] <- mean(cpg[list,4])
	}
}