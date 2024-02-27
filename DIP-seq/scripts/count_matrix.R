
.libPaths()

# Combine the featureCounts output per sample into one matrix.
#
# Usage: Rscript combine-featurecounts.R [file, file, ...] > output
#
# Each input file is the featureCounts results for one sample. The gene names
# are extracted from column 1 and the counts from column 7.
#
# The results are exported as tab-separated columns to standard out.

suppressPackageStartupMessages(library("edgeR"))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("stringr"))

# Input files ------------------------------------------------------------------
parser <- ArgumentParser(description='MEDIPS quality control step')
parser$add_argument('-i', '--input', type='character', required=TRUE,nargs='+',help='Input single sample count files')
parser$add_argument('-o','--output',type='character',help='filename for matrix output')
parser$add_argument('-s','--surfix',type='character',default="gene",help='gene/promoter/cgi for DMR')
args <- parser$parse_args()
f_counts <- args$input
surfix <- args$surfix
outfile <- args$output

print(f_counts)
print(surfix)
print(outfile)
# For testing:
# f_counts <- Sys.glob("/scratch/midway2/jdblischak/cardioqtl-counts/*txt")
raw <- readDGE(f_counts, columns = c(1, 7), comment.char = "#")
counts <- as.matrix(raw)
# Change column names from filenames to sample names ---------------------------
colnames(counts) <- colnames(counts) %>%
  str_split("/", simplify = TRUE) %>%
  `[`(, ncol(.)) %>%
  str_replace(paste("_",surfix,sep=""), "")

# Sort the rownames and colnames -----------------------------------------------
#counts[,"rown"] <- rownames(counts)
#counts <- counts[,c("rown",colnames(counts))]
counts <- counts[order(rownames(counts)), order(as.numeric(colnames(counts)))]
# Write to standard out --------------------------------------------------------
# head(counts)
write.table(counts, file = outfile, quote = FALSE, sep = "\t",col.names = T, row.names = T)
