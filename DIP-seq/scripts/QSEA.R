
## seems can only give an estimated beta value on a window/Bin

options(stringsAsFactors = F)
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/")
library("BSgenome")
library("qsea")
library("MEDIPS")
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)

#BiocManager::install("MEDIPSData")


# test
library("BiocParallel")
register(MulticoreParam(workers=12))

data(samplesNSCLC, package="MEDIPSData")
knitr::kable(samples_NSCLC)
path=system.file("extdata", package="MEDIPSData")
samples_NSCLC$file_name=paste0(path,"/",samples_NSCLC$file_name )
#or
samples_NSCLC <- read.table("./meta/lulab/testQSEA.txt",header = T)

## 1.All relevant information of the enrichment experiment, including sample information, the genomic read coverage, CpG density are stored in a “qseaSet” objec
qseaSet=createQseaSet(sampleTable=samples_NSCLC, 
                      BSgenome="BSgenome.Hsapiens.UCSC.hg38", 
                      chr.select=paste0("chr", 20:22), 
                      window_size=500)
qseaSet

## 2.read the alignment files and compute the MeDIP coverage for each window
qseaSet=addCoverage(qseaSet, uniquePos=TRUE, paired=TRUE)

## 3.normalize CNV
qseaSet=addCNV(qseaSet, file_name="file_name",window_size=2e6, 
               paired=TRUE, parallel=FALSE, MeDIP=TRUE)

## 4.normalize lib (TMM)
qseaSet=addLibraryFactors(qseaSet)


## 5.Estimating model parameters for transformation to absolute methylation values
qseaSet=addPatternDensity(qseaSet, "CG", name="CpG")

## 6.From the regions without CpGs we can estimate the coverage offset from background reads.
qseaSet = addOffset(qseaSet, enrichmentPattern = "CpG")

## 7.add estimate the relative enrichment
data(tcga_luad_lusc_450kmeth, package="MEDIPSData")

wd=findOverlaps(tcga_luad_lusc_450kmeth, getRegions(qseaSet), select="first")
signal=as.matrix(mcols(tcga_luad_lusc_450kmeth)[,rep(1:2,3)])
qseaSet=addEnrichmentParameters(qseaSet, enrichmentPattern="CpG", 
                                windowIdx=wd, signal=signal)
#or
wd=which(getRegions(qseaSet)$CpG_density>1 &
           getRegions(qseaSet)$CpG_density<15)
signal=(15-getRegions(qseaSet)$CpG_density[wd])*.55/15+.25
qseaSet_blind=addEnrichmentParameters(qseaSet, enrichmentPattern="CpG", 
                                      windowIdx=wd, signal=signal)

## 8.QC
### At first, we check the estimated fraction of background reads:
getOffset(qseaSet, scale="fraction")
plotEPmatrix(qseaSet)

## 9.##Exploratory Analysis
plotCNV(qseaSet)

data(annotation, package="MEDIPSData")
CGIprom=intersect(ROIs[["CpG Island"]], ROIs[["TSS"]],ignore.strand=TRUE)
pca_cgi=getPCA(qseaSet, norm_method="beta", ROIs=CGIprom)
col=rep(c("red", "green"), 3)
plotPCA(pca_cgi, bg=col, main="PCA plot based on CpG Island Promoters")
