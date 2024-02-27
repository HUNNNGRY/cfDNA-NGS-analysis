# prepare wgs mat
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq")

# NOV ---------------------------------------------------------------------
## combine 2 regions
nuleosome1 <- read.table("output/lulab/matrix/count_matrix_promoter150TSS50.txt",sep = "\t",header = T,row.names = 1,check.names = F)
nuleosome2 <- read.table("output/lulab/matrix/count_matrix_promoter300100exon1end.txt",sep = "\t",header = T,row.names = 1,check.names = F)

## Keep ensg only
rownames(nuleosome1) <- as.character(lapply(strsplit(rownames(nuleosome1),"\\|"), function(x) x[1]))
rownames(nuleosome1) <- sub("promoter150TSS50_","",rownames(nuleosome1) )

rownames(nuleosome2) <- as.character(lapply(strsplit(rownames(nuleosome2),"\\|"), function(x) x[1]))
rownames(nuleosome2) <- sub("promoter300100exon1end_","",rownames(nuleosome2) )

nuleosome1 <- nuleosome1[order(rownames(nuleosome1)),]
nuleosome2 <- nuleosome2[order(rownames(nuleosome2)),]

all(colnames(nuleosome1) == colnames(nuleosome2))
all(rownames(nuleosome1) == rownames(nuleosome2))
nuleosome2 <- nuleosome2[rownames(nuleosome1),]


## sum
nuleosome.combine <- nuleosome1+nuleosome2
rownames(nuleosome.combine) <- paste0("NOVsum_",rownames(nuleosome.combine),"|602")
write.table(nuleosome.combine,"output/lulab/matrix/count_matrix_NOVsum.txt",sep = "\t",row.names = T,col.names = T,quote = F)

## max
#str(nuleosome1)
nuleosome.max <- pmax(as.matrix(nuleosome1),as.matrix(nuleosome2))
rownames(nuleosome.max) <- paste0("NOVmax_",rownames(nuleosome.max),"|301")
write.table(nuleosome.max,"output/lulab/matrix/count_matrix_NOVmax.txt",sep = "\t",row.names = T,col.names = T,quote = F)
#nuleosome2[1:3,1:3]

## min
nuleosome.min <- pmin(as.matrix(nuleosome1),as.matrix(nuleosome2))
rownames(nuleosome.min) <- paste0("NOVmin_",rownames(nuleosome.min),"|301")
write.table(nuleosome.min,"output/lulab/matrix/count_matrix_NOVmin.txt",sep = "\t",row.names = T,col.names = T,quote = F)


# get TMM CPM
$ cnode
$/usr/bin/Rscript scripts/run-NormCountMat.R  \
  -i output/lulab/matrix/count_matrix_NOVsum.txt \
  -o output/lulab/matrix/CPMtmm_matrix_NOVsum.txt.txt
$/usr/bin/Rscript scripts/run-NormCountMat.R  \
  -i output/lulab/matrix/count_matrix_NOVmax.txt \
  -o output/lulab/matrix/CPMtmm_matrix_NOVmax.txt.txt
$/usr/bin/Rscript scripts/run-NormCountMat.R  \
  -i output/lulab/matrix/count_matrix_NOVmin.txt \
  -o output/lulab/matrix/CPMtmm_matrix_NOVmin.txt.txt

