#sva test
#last update 210328
setwd("/BioII/lulab_b/baopengfei/projects/methylation")
#### fast start ####
{
library(sva)
cdata <- as.matrix(cdata)
ph <- read.table("sif.txt", header = T, sep = "\t", row.names = 1)
modcombat = model.matrix(~1, data = ph)
combat_edata = ComBat(dat=cdata, batch=ph$batch, mod=modcombat, par.prior=TRUE, prior.plots=TRUE)
#此步骤会产生一个和输入文件相同的数据矩阵
}

#### manual of sva ####
# 1.setup data 
library(sva)
library(bladderbatch)
data(bladderdata)
library(pamr)
library(limma)
#sva package assumes there are two types of variables that are being considered: 
#(1) adjustment variables and (2) variables of interest.
pheno = pData(bladderEset)
edata = exprs(bladderEset)
mod = model.matrix(~as.factor(cancer), data=pheno)  # the full model matrix - including both the adjustment variables and the variable of interest
mod0 = model.matrix(~1,data=pheno)    # The null model - contains only the adjustment variables.

# 2.apply sva func.
n.sv = num.sv(edata,mod,method="leek")
n.sv

svobj = sva(edata,mod,mod0,n.sv=n.sv)

# 3.adjust surrogate var using f.pvalue
pValues = f.pvalue(edata,mod,mod0)
qValues = p.adjust(pValues,method="BH")
table(qValues <= 0.05)

modSv = cbind(mod,svobj$sv)
mod0Sv = cbind(mod0,svobj$sv)
pValuesSv = f.pvalue(edata,modSv,mod0Sv)
qValuesSv = p.adjust(pValuesSv,method="BH")
table(qValuesSv <= 0.05)

# 4.adjust for surrogate var using limma
fit = lmFit(edata,modSv)

contrast.matrix <- cbind("C1"=c(-1,1,0,rep(0,svobj$n.sv)),"C2"=c(0,-1,1,rep(0,svobj$n.sv)))  # miss info
fitContrasts = contrasts.fit(fit,contrast.matrix)

eb = eBayes(fitContrasts)
topTableF(eb, adjust="BH")

# 5.apply combat func. to adjst for known batches
batch = pheno$batch
modcombat = model.matrix(~1, data=pheno)
combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE,prior.plots = T)
pValuesComBat = f.pvalue(combat_edata,mod,mod0)
qValuesComBat = p.adjust(pValuesComBat,method="BH")

# 6.count data: ComBat-Seq  (sva > 3.36)
{#ComBat-Seq is an improved model based on the ComBat framework, which speciﬁcally targets RNA-Seq count data. 
#It uses a negative binomial regression to model the count matrix
count_matrix <- matrix(rnbinom(400, size=10, prob=0.1),nrow=50, ncol=8)
batch <- c(rep(1, 4), rep(2, 4))
adjusted <- ComBat_seq(count_matrix, batch=batch, group=NULL)
}

#### test in-house cfmeidp data ####
# 1.read tpm matrix
tpm.pro <- read.delim(paste0("./data/lulab/07counts/TPM_matrix_promoter.txt"),row.names=1,check.names = F)
vsd.pro <- read.delim(paste0("./output/lulab/all-CRCvsNC-promoter-DESeq2-designBygroup-vsd-blindTrue.txt"),row.names=1,check.names = F)


# 2.read meta table
{meta1 <- read.table("./meta/lulab/sample_table.txt",sep = "\t",header = T,stringsAsFactors = F,fill = T)
colnames(meta1)[colnames(meta1)=="Group"] <- "group"
colnames(meta1)[colnames(meta1)=="Sample_id"] <- "id"
rownames(meta1) <- meta1$id
meta1$QC <- "failQC"
meta1$QC[meta1$relH.3 >= 2.0 & meta1$GoGe.1.7 >= 1.5 & meta1$saturation > 0.9] <- "passQC"
table(meta1$QC)
meta1 <- meta1[order(meta1$group,meta1$id),]

meta <- meta1[meta1$QC=="passQC",]
filter.group <- meta$group %in% c("STAD")
filter.sample <- meta$id %in% c("NC-PKU-mix0-me","CRC-PKU-mix0-me","NC-PKU-mix1-me","CRC-PKU-mix1-me")
filter.colum <- 
meta <- meta[!filter.group & !filter.sample,]
table(meta$group)
str(meta)
meta$group <- factor(x = meta$group,levels=c("NC","CRC"))
meta$Lib.batch <- as.factor(meta$Lib.batch)
meta <- meta[,colSums(is.na(meta)) == 0]   # remove colums that contains at least one NA
}
x <- vsd.pro
x <- x[,match(rownames(meta),colnames(x))]
all(colnames(x) %in% rownames(meta))
all(colnames(x) == rownames(meta))
#dim(x)
dim(meta)

# 3.run sva
pheno <- meta
#colnames(meta)
mod = model.matrix(~as.factor(group), data=pheno)  # the full model matrix: adjustment variables + variable of interest
mod0 = model.matrix(~1,data=pheno)    # The null model: only the adjustment variables.

# 2.apply sva func.
#class(x)
edata <- as.matrix(x)

n.sv = num.sv(edata,mod,method="leek")
n.sv    #in-house TPM data has 35 sr ???, vsd has only 1
#table(colSums(is.nan(edata)) == 0)

svobj = sva(dat = edata,mod = mod,mod0 = mod0,n.sv=n.sv)  # TPM error: 'x' contains missing values
#svaseq(edata,mod,mod0,n.sv=n.sv)  #estimating surrogate variables for count based RNA-seq data.

# 3.adjust surrogate var using f.pvalue
pValues = f.pvalue(edata,mod,mod0)
qValues = p.adjust(pValues,method="BH")
table(qValues <= 0.05)

modSv = cbind(mod,svobj$sv)
mod0Sv = cbind(mod0,svobj$sv)
pValuesSv = f.pvalue(edata,modSv,mod0Sv)
qValuesSv = p.adjust(pValuesSv,method="BH")
table(qValuesSv <= 0.05)  # 略微增加

# 4.adjust for surrogate var using limma
{
fit = limma::lmFit(edata,modSv)

contrast.matrix <- cbind("C1"=c(-1,1,0,rep(0,svobj$n.sv)),"C2"=c(0,-1,1,rep(0,svobj$n.sv)))  # miss info
fitContrasts = limma::contrasts.fit(fit,contrast.matrix)   # error...

eb = eBayes(fitContrasts)
topTableF(eb, adjust="BH")
}
# 5.apply combat func. to adjst for known batches
batch = as.factor(pheno$Operator)
modcombat = model.matrix(~ batch + group, data=pheno)
combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE,prior.plots = T)  #Error: 'x' contains missing values
pValuesComBat = f.pvalue(combat_edata,mod,mod0)
qValuesComBat = p.adjust(pValuesComBat,method="BH")

cdata <- vsd.pro
ph <- meta
modcombat = model.matrix(~as.factor(ph$group), data = ph)
combat_edata = ComBat(dat=cdata, batch=ph$Sample.batch)
#Error in tcrossprod(t(design), as.matrix(dat)) : 
#  non-conformable arguments

