#!/usr/bin/R
#suppressPackageStartupMessages(require(parallel))
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/BS-seq")
suppressPackageStartupMessages(require(foreach))
suppressPackageStartupMessages(require(doParallel))
registerDoParallel(cores=12)

#getOption("cores")
#option(cores=8)
dataset <- "GSE149438"
gtf <- read.table("genome/gtf/gencodev27_onlypromoter.gtf",stringsAsFactors=F) 
sample <- read.table(paste0("metadata/",dataset,"/sample_ids_test.txt"),stringsAsFactors=F)
region <- matrix(data=NA,nrow=length(gtf[,10]),ncol=length(sample$V1),dimnames=list(c(gtf[,10]),c(sample$V1)))


for (id in sample$V1)
{
print(paste0("start:",id))
cpg <- read.table(paste0("output/",dataset,"/",id,".deduplicated.bismark.cov.gz"),stringsAsFactors=F)
	#mapply(region,cpg2region)
	foreach (i=1:nrow(region)) %dopar% 
	{
		print(i)
		list <- rownames(cpg)[cpg[,2] >= gtf[i,4] & cpg[,2] <= gtf[i,5]]
		region[i,id] <- mean(cpg[list,4])
	}
}


## single thread:
# for (id in sample$V1)
# {
# print(paste0("start:",id))
# cpg <- read.table(paste0("output/",dataset,"/",id,".deduplicated.bismark.cov.gz"),stringsAsFactors=F)
# 	for (i in 1:nrow(region))
# 	{
# 		print(i)
# 		list <- rownames(cpg)[cpg[,2] >= gtf[i,4] & cpg[,2] <= gtf[i,5]]
# 		region[i,id] <- mean(cpg[list,4])
# 	}
# }

## func:
# cpg2region <- function(region,cpg,gtf,i) {
# 	list <- rownames(cpg)[cpg[,2] >= gtf[i,4] & cpg[,2] <= gtf[i,5]]
# 	region[i,id] <- mean(cpg[list,4])	
# 	return(region)
# }