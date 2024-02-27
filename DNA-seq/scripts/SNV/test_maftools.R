# oncoplot
# last 230919 by bpf (adapted from tyh)
# used for cfDNA SNV analysis visualization tutorial
library(maftools)
library(BSgenome.Hsapiens.UCSC.hg38)

setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq")
annovar_merge <- $1 # output/test_stepByStep/mutect2-vcf-filtered/annotation/annovar_merge.vcf
annovar_merge <- "/BioII/lulab_b/baopengfei/projects/multi-omics-explore/exOmics/DNA-seq/output/lulab/vcf-filtered/annotation/all_sample.csv" # /BioII/lulab_b/baopengfei/projects/multi-omics-explore/exOmics/DNA-seq/output/lulab/vcf-filtered/annotation/all_sample.csv
onc_driver <- "/BioII/lulab_b/baopengfei/shared_reference/geneset/IntOGen-Drivers-20200213/Compendium_Cancer_Genes.tsv"
outdir <- "output/test_stepByStep/mutect2-vcf-filtered/annotation"
cancer_type <- "COREAD"

options(stringsAsFactors = F)
m <-  data.table::fread(annovar_merge,check.names = F, sep = "\t", header = T)
m <- as.data.frame(m)
# head(m[,1:8],40)

d <-  data.table::fread(onc_driver,check.names = F, sep = "\t",header = T)
d <- as.data.frame(d)
table(d$CANCER_TYPE)
d <- d[d$CANCER_TYPE==cancer_type,] # HC: HCC; HNSC: Head neck cancer; COREAD: CRC
#write.table(unique(d$SYMBOL),"./ref/CRC_driverGene.txt",quote = F,col.names = F,row.names = F,sep = "\t")

#m.crc <- m[grep(unique(d$SYMBOL),m$Gene.refGene),]   
m.crc <- m[m$Gene.refGene %in% unique(d$SYMBOL),] # ,collapse="|"
unique(m.crc$Gene.refGene)
table(m.crc$Tumor_Sample_Barcode)

write.table(m.crc,paste0(outdir,"/",cancer_type,"_driverGene_table.txt"),quote = T,col.names = T,row.names = F,sep = ",")
write.table(m,paste0(outdir,"/all_Gene_table.txt"),quote = T,col.names = T,row.names = F,sep = ",")


###read group meta
pos_smps_path <- "data/test/meta_data/sample_ids_CRC.txt"
neg_smps_path <- "data/test/meta_data/sample_ids_NC.txt"

ct <- read.delim(pos_smps_path,sep = "\t",header = F)$V1 # pos sample_ids
pt <- read.delim(neg_smps_path,sep = "\t",header = F)$V1 # neg sample_ids


## waterfall plot by maftools
library(maftools)
## Summary
# annovar.laml.all <- annovarToMaf( annovar = "output/wes-mutect2-FilterMutectCalls/annotation/all_Gene_table.txt",   # ./scripts/CRC_sample_CRC_COSMIC_cancergenes_filtered.txt
#                                   Center = 'NA',
#                                   refBuild = 'hg38', 
#                                   tsbCol = 'Tumor_Sample_Barcode', 
#                                   #table = 'refGene',
#                                   sep = ",", MAFobj = F)

annovar.laml <- annovarToMaf( annovar = paste0(outdir,"/all_Gene_table.txt"),   # ./scripts/CRC_sample_CRC_COSMIC_cancergenes_filtered.txt
                              Center = 'NA',
                              refBuild = 'hg38', 
                              tsbCol = 'Tumor_Sample_Barcode', 
                              #table = 'refGene',
                              sep = ",", MAFobj = T)
laml=annovar.laml
# unique(annovar.laml@data$Tumor_Sample_Barcode)
# getSampleSummary(laml)
# getGeneSummary(laml)
# getFields(laml)


var.annovar.maf = annovarToMaf( annovar = paste0(outdir,"/",cancer_type,"_driverGene_table.txt"),   # ./scripts/CRC_sample_CRC_COSMIC_cancergenes_filtered.txt
                                Center = 'NA',
                                refBuild = 'hg38', 
                                tsbCol = 'Tumor_Sample_Barcode', 
                                #table = 'refGene',
                                sep = ",")


## op2:ct+pt
pos_lab <- "CRC"
neg_lab <- "NC"
var.annovar.maf.1 <- var.annovar.maf[var.annovar.maf$Tumor_Sample_Barcode %in% c(ct,pt),]
write.table(var.annovar.maf.1,file=paste0(outdir,"/","var_annovar_",pos_lab,"vs",neg_lab,".maf"),quote= T,sep="\t",row.names=F)
var_maf = read.maf(maf = paste0(outdir,"/","var_annovar_",pos_lab,"vs",neg_lab,".maf")
                         # useAll = T # use all mutaions or Only Somatic
                         )


## op2:ct+pt
var_maf@clinical.data$Type <- c(rep("CT",9),rep("PT",8))
var_maf@clinical.data


pdf(paste0(outdir,"/plotmafSummary.pdf"),width = 11,height = 9)
plotmafSummary(maf = var_maf, textSize = 0.6, fs = 1, titleSize = c(1,0.5),showBarcodes = T, rmOutlier = TRUE, addStat = 'median') # addStat="mean",
dev.off()


## op2:ct+pt
pdf(paste0(outdir,"/oncoplot.ct.pt.pdf"),width = 10,height = 12)
oncoplot(maf = var_maf, top = 40 ,showTumorSampleBarcodes = T,clinicalFeatures = "Type",colbar_pathway = T,
         sortByAnnotation = T, # sortByMutation = T,
         draw_titv = TRUE)  #colbar_pathway=F
dev.off()


#keepGeneOrder=T, GeneOrderSort=T,
#sampleOrder
laml.titv = titv(maf = var_maf, plot = FALSE, useSyn = TRUE)
plotTiTv(res = laml.titv)

pdf(paste0(outdir,"/somaticInteractions.pdf",width = 6,height = 7)
somaticInteractions(maf = var_maf, top = 10, pvalue = c(0.05, 0.1))
dev.off()

# #lollipop plot for APC / TP53 
# lollipopPlot(
#   maf = var_maf,
#   gene = 'TP53',
#   AACol = 'aaChange',
#   showMutationRate = TRUE,
#   #labelPos = "all",
#   #refSeqID = "NM_000546",
#   printCount = TRUE,
#   showDomainLabel = TRUE,
#   legendTxtSize = 1.2
# )









### TODO:
# ## mutational signature (SNP)
# rm(list=ls())
# options(stringsAsFactors=FALSE)
# ###切换镜像
# options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
# options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# ###安装R包
# ##install.packages('deconstructSigs')
# ###BiocManager::install('BSgenome')
# ###BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')
# library(deconstructSigs)
# library(BSgenome.Hsapiens.UCSC.hg38)
# ###https://github.com/raerose01/deconstructSigs
# 
# ###读入数据
# var.annovar.maf = annovarToMaf( annovar = paste0(outdir,"/all_Gene_table.txt"),   # ./scripts/CRC_sample_CRC_COSMIC_cancergenes_filtered.txt
#                                 Center = 'NA',
#                                 refBuild = 'hg38', 
#                                 tsbCol = 'Tumor_Sample_Barcode', 
#                                 #table = 'refGene',
#                                 sep = ",")
# write.table(var.annovar.maf,file=paste0(outdir,"/all_Gene_table.maf"),quote= T,sep="\t",row.names=F)
# 
# maf = read.maf(maf = paste0(outdir,"/all_Gene_table.maf"))
# maf <- as.data.frame(maf@data)
# #maf=read.table('var_annovar.maf',header = T,sep = '\t',quote = "")
# 
# table(maf$Variant_Type) 
# maf <- maf[maf$Variant_Type=="SNP",]
# maf[1:3,]
# colnames(maf)
# 
# ###构建肿瘤突变数据框，需要5列信息: sample.ID,chr,pos,ref,alt 
# sample.mut.ref <- data.frame(Sample=maf[,6], 
#                              chr = maf[,1],
#                              pos = maf[,2],
#                              ref = maf[,4],
#                              alt = maf[,5])
# #?mut.to.sigs.input
# #str(sample.mut.ref)
# #sample.mut.ref$Sample <- as.character(sample.mut.ref$Sample)
# 
# table(sample.mut.ref$Sample %in% ct) # 3 types + all
# sample.mut.ref.ct <- sample.mut.ref[sample.mut.ref$Sample %in% ct,]
# table(sample.mut.ref.ct$Sample)
# sample.mut.ref.ct$Sample <- "CT"
#   
# ###run
# sigs.input <- deconstructSigs::mut.to.sigs.input(mut.ref = sample.mut.ref.ct,
#                                 chr = "chr", 
#                                 pos = "pos", 
#                                 ref = "ref", 
#                                 alt = "alt",
#                                 bsg = BSgenome.Hsapiens.UCSC.hg38)
# 
# 
# ###read signatures.cosmic
# signatures.cosmic <- read.table("/BioII/lulab_b/baopengfei/shared_reference/geneset/COSMIC/COSMIC_v3.2_SBS_GRCh38.txt",stringsAsFactors = F,sep = "\t",header = T,row.names = 1)
# signatures.cosmic <- as.data.frame(t(signatures.cosmic))
# #signatures.cosmic need tobe a df, where rows are signatures (SBS1), columns are trinucleotide contexts (A[C>G]T)
# 
# ###plot
# sigs.output <- deconstructSigs::whichSignatures(tumor.ref = sigs.input,
#                                signatures.ref = signatures.cosmic, 
#                                sample.id="CT",#sample.id = 'RT', 'PT', 'CT'  (treat one group as one sample)
#                                contexts.needed = TRUE)
# 
# 
# # Plot output
# pdf("./CT.deconstructSigs.pdf",width = 7,height = 7)
# deconstructSigs::plotSignatures(sigs.output)
# dev.off()
# 
# #sample freq
# #ref freq
# #delta freq (sample-ref)
