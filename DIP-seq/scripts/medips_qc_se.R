#.libPaths(new=c("/BioII/lulab_b/baopengfei/R/x86_64-pc-linux-gnu-library/3.5", "/BioII/lulab_b/baopengfei/anaconda3/lib/R/library", "/opt/microsoft/ropen/3.5.1/lib64/R/library"))

.libPaths()

message("Load required packages ...")
suppressPackageStartupMessages(library(MEDIPS))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressPackageStartupMessages(library(argparse))
message("all packages loaded")

parser <- ArgumentParser(description='MEDIPS quality control step')
parser$add_argument('-i', '--input', type='character', required=TRUE,help='Input bam file')
parser$add_argument('-w','--window-size',type='integer',default=300,help='Window size for calculation')
parser$add_argument('-sr', '--saturation', type='character',required=TRUE,help='saturation-out-path')
parser$add_argument('-er', '--enrich', type='character',required=TRUE,help='enrich-out-path')
parser$add_argument('-cr', '--coverage', type='character',required=TRUE,help='coverage-out-path')
args <- parser$parse_args()

BSgenome <- 'BSgenome.Hsapiens.UCSC.hg38'
uniq <- 1
chr.select <- c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12',
            'chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22')

#1.enrich
message("start enrich")
er = MEDIPS.CpGenrich(file = args$input, BSgenome = BSgenome, chr.select = chr.select,uniq = uniq,paired=FALSE)
er <- t(data.frame(er))
write.table(er,file=args$enrich,sep="\t",col.names=FALSE,quote=FALSE)
message("end enrich")

#2.saturation
message("start saturation")
sr = MEDIPS.saturation(file=args$input, BSgenome = BSgenome, chr.select = chr.select, uniq = uniq, window_size = args$window_size,paired=FALSE)

saturation.table <- sr$distinctSets
colnames(saturation.table) <- c('subset','correlation')
saturation.table <- transform(saturation.table,data=rep("observed",nrow(saturation.table)))

estimate.table <- sr$estimation
colnames(estimate.table) <- c('subset','correlation')
estimate.table <- transform(estimate.table,data=rep("estimated",nrow(estimate.table)))

saturation.table <-rbind(saturation.table,estimate.table)
write.table(saturation.table, file=args$saturation, sep='\t',row.names = FALSE,quote=FALSE)
message("end saturation")

#3.coverage
message("start coverage")
cr<-MEDIPS.seqCoverage(file=args$input,pattern="CG",BSgenome=BSgenome,uniq=uniq,paired=FALSE)
png(file=args$coverage,type="cairo")
MEDIPS.plotSeqCoverage(cr,type="hist",t = 15,main="Sequence pattern coverage histogram")
dev.off()
message("end coverage")