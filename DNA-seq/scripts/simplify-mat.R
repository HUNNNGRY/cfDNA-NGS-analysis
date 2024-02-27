# simply matrix:  
# 1. gene name from promoter_ENSG0029382.1|2001 to ENSG0029382
# 2. filter genes
# 3. filter samples
# last 211027 by bpf
library(argparse)
parser <- ArgumentParser(description='simply matrix: convert row gene name to simpler version, filter rows/genes and cols/samples')
parser$add_argument('-i', '--matrix', type='character', required=TRUE,
                    help='input matrix. usually Rows are genes. Columns are samples.')
parser$add_argument('-s', '--sample', type='character', default="all",
                    help='cols/samples (eg. NC-PKU-12-wgs) to keep, need in txt format (default: all)')
parser$add_argument('-g', '--gene', type='character', default="all",
                    help='rows/genes (eg. ENSG00000238) to keep, need in txt format (default: all)')
parser$add_argument('-s.s', '--sample.surfix', type='character', default="",
                    help='surfix to added after cols/samples id (default: "")')
parser$add_argument('-g.s', '--gene.surfix', type='character', default="",
                    help='surfix to added after rows/genes id (default: "")')
parser$add_argument('-s.d', '--sample.delete', type='character', default="none",
                    help='cols/samples (eg. NC-PKU-12-wgs) to delete/remove, need in txt format (default: none)')
parser$add_argument('-g.d', '--gene.delete', type='character', default="none",
                    help='rows/genes (eg. ENSG00000238) to delete/remove, need in txt format (default: none)')
parser$add_argument('-o', '--outfile', type='character', default="./simplif_matrix.txt",
                    help='output file path, ENSG may be simplified to ^ENSG00009999, default:./simplif_matrix.txt')
parser$add_argument('--format', type='character', default="tsv",
                    help='output file format, tsv or csv, default:tsv')
args <- parser$parse_args()

ma <- args$matrix
sa <- args$sample
ge <- args$gene
sa.d <- args$sample.delete
ge.d <- args$gene.delete
sa.s <- args$sample.surfix
ge.s <- args$gene.surfix
ou <- args$outfile
format <- args$format

# # test
# ma <- "/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/count_matrix_promoter300100exon1end.txt"
# sa <- "/BioII/lulab_b/baopengfei/sample.txt"
# ge <- "/BioII/lulab_b/baopengfei/shared_reference/geneset/driver_ensg.txt"
# ou <- "/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/outfile.txt"


if (sa!="all") {
message(sa)
sam <- read.delim(sa,header = F,stringsAsFactors = F)$V1
}
if (ge!="all") {
message(ge)
gen <- read.delim(ge,header = F,stringsAsFactors = F)$V1
}
if (sa.d!="none") {
message(sa.d)
sam.d <- read.delim(sa.d,header = F,stringsAsFactors = F)$V1
}
if (ge.d!="none") {
message(ge.d)
gen.d <- read.delim(ge.d,header = F,stringsAsFactors = F)$V1
}


message(paste0("read in ",ma))
mat <- read.csv(ma,sep = "\t",check.names = F,header = TRUE, row.names = 1)
#mat[1:4,1:4]

# strip possible prefix
if (length(grep("^ENSG",rownames(mat)))>=2) {
  message("input gene id is ^ENSG*")
  mat$ENSG <- rownames(mat)
} else {
  message("input gene id is not ^ENSG*")
  #mat$ENSG <- unlist(lapply(strsplit(rownames(mat),"_",fixed = TRUE),function(x) x[2]))
  prefix <- paste("gene","promoter","promoter5k","promoter150TSS50","promoter300100exon1end","gene_exon1","CDS","UTR3","UTR5","CpG_island","LINE","SINE","retroposon",sep = "_|")
  prefix <- paste0(prefix,"_")
  mat$ENSG <- sub(prefix,"",rownames(mat) )
}

# strip possible surfix (ensg length and ensg edition)
mat$ENSG <- unlist(lapply(strsplit(mat$ENSG,".",fixed = TRUE),function(x) x[1]))
mat <- mat[!duplicated(mat$ENSG),]
rownames(mat) <- mat$ENSG
mat <- mat[,colnames(mat)!="ENSG",drop=F] # prevent only one sample left

# keep gene and sample
if (ge=="all") {
  message("keep all genes")
} else {
  message("keep ",length(gen)," genes")
  mat <- mat[gen,]
}

if (sa=="all") {
  message("keep all samples")
} else {
  message("keep ",length(sam)," samples")
  mat <- mat[,sam]
}


# delete gene and sample
if (ge.d=="none") {
  message("keep all genes")
} else {
  message("delete ",length(gen.d)," genes")
  mat <- mat[!(rownames(mat) %in% gen.d),]
}

if (sa.d=="none") {
  message("keep all samples")
} else {
  message("delete ",length(sam.d)," samples")
  mat <- mat[,!(colnames(mat) %in% sam.d)]
}

rownames(mat) <- paste0(rownames(mat),ge.s)
colnames(mat) <- paste0(colnames(mat),sa.s)

message("final rows: ",nrow(mat))
message("final cols: ",ncol(mat))
message(paste0("write to ",ou))
if(format=="csv"){
  write.table(mat,ou,sep = ",",quote = F,col.names = T,row.names = T)
}else{
  write.table(mat,ou,sep = "\t",quote = F,col.names = T,row.names = T)
}

# get cancer driver ensg
#hallmark <- read.csv("/BioII/lulab_b/baopengfei/shared_reference/geneset/cancer_hallmark/PATH_ID_NAME_hallmark.tsv",sep = "\t",check.names = F,header = TRUE)
# driver <- read.csv("/BioII/lulab_b/baopengfei/shared_reference/geneset/IntOGen-Drivers-20200213/Compendium_Cancer_Genes.tsv",sep = "\t",check.names = F,header = TRUE)
# gene.list.ncbiid <- read.table("/BioII/lulab_b/baopengfei/shared_reference/geneset/encodev27_gene_v2.txt",sep = "\t",header = T)
# #gene.list.ncbiid[1:3,]
# colnames(gene.list.ncbiid)[12] <- "SYMBOL"
# gene.list.ncbiid <- gene.list.ncbiid[,c("SYMBOL","ensg"),drop=F]
# #length(unique(driver$SYMBOL))
# #table(unique(driver$SYMBOL) %in% gene.list.ncbiid$name)
# #table(duplicated(gene.list.ncbiid$name))
# driver <- driver[!duplicated(driver$SYMBOL),"SYMBOL",drop=F]
# driver <- dplyr::left_join(driver,gene.list.ncbiid)
# 
# write.table(driver$ensg,"/BioII/lulab_b/baopengfei/shared_reference/geneset/driver_ensg.txt",sep = "\t",quote = F,col.names = F,row.names = F)







