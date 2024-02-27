#!/usr/bin/env Rscript
#.libPaths(new=c("/BioII/lulab_b/baopengfei/R/x86_64-pc-linux-gnu-library/3.5", "/BioII/lulab_b/baopengfei/anaconda3/lib/R/library", "/opt/microsoft/ropen/3.5.1/lib64/R/library"))
#.libPaths()

message("Load required packages ...")
suppressPackageStartupMessages(library(MEDIPS))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressPackageStartupMessages(library(argparse))
message("all packages loaded \n")

parser <- ArgumentParser(description='WGS quality control step')
parser$add_argument('-i', '--input', type='character', required=TRUE,help='Input bam file')
parser$add_argument('-w','--window-size',type='integer',default=300,help='Window size for calculation')
parser$add_argument('-sr', '--saturation', type='character',required=TRUE,help='saturation-out-path')
parser$add_argument('-er', '--enrich', type='character',required=TRUE,help='enrich-out-path')
#parser$add_argument('-cr', '--coverage', type='character',required=TRUE,help='coverage-out-path')
args <- parser$parse_args()

BSgenome <- 'BSgenome.Hsapiens.UCSC.hg38'
uniq <- 1
chr.select <- c(paste0('chr',1:22)) 

#1.enrich
message("start enrich")
er = MEDIPS.CpGenrich(file = args$input, BSgenome = BSgenome, chr.select = chr.select,uniq = uniq,paired=TRUE)
er <- t(data.frame(er))
write.table(er,file=args$enrich,sep="\t",col.names=FALSE,quote=FALSE)
message("end enrich: \n")
message(paste0(args$enrich,"\t"))

#2.saturation
for (ws in c("300","3000","30000")) 
{
message(paste0("start ",ws," saturation"))
sr = MEDIPS.saturation(file=args$input, BSgenome = BSgenome, chr.select = chr.select, uniq = uniq, window_size = as.numeric(ws), paired=TRUE)

saturation.table <- sr$distinctSets
colnames(saturation.table) <- c('subset','correlation')
saturation.table <- transform(saturation.table,data=rep("observed",nrow(saturation.table)))

estimate.table <- sr$estimation
colnames(estimate.table) <- c('subset','correlation')
estimate.table <- transform(estimate.table,data=rep("estimated",nrow(estimate.table)))

saturation.table <-rbind(saturation.table,estimate.table)
write.table(saturation.table, file=paste(args$saturation,ws,sep="."), sep='\t',row.names = FALSE,quote=FALSE)
message(paste0("end ",ws," saturation:\n"))
message(paste(args$saturation,ws,sep="."))
}
