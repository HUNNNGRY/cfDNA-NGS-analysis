# correct WGS GC (run)
# last 220111 by bpf
#correct GC-bias in count mat not CPM mat 

options(stringsAsFactors = F)

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library((BSgenome.Hsapiens.UCSC.hg38)))
parser <- ArgumentParser(description='GO/KEGG/DO geneset enrichment by ORA/GSEA using diff table (network independent version)')
parser$add_argument('-i', '--matrix', type='character', required=TRUE,
                    help='input count matrix. Rows are genes (^ENSG*, *_ENSG*). Columns are samples')
parser$add_argument('-r', '--bed', type='character', required=TRUE,
                    help='reference bed file path')
parser$add_argument('-o', '--outfile', type='character', required=TRUE,
                    help='output GC-content corrected count matrix. Rows are genes (^ENSG*, *_ENSG*). Columns are samples')
args <- parser$parse_args()
matrix <- args$matrix
bed <- args$bed
outfile <- args$outfile


#get gene&promoter count GC adjusted
correct <- function(sample.id){
  # print(sample.id)
  
  pos <- hg38$pos[match(g,hg38$ENSG)]
  # chr <- hg38$chr[match(g,hg38$ENSG)]
  # start <- hg38$start[match(g,hg38$ENSG)]
  # end <- hg38$end[match(g,hg38$ENSG)]
  
  test <- cbind(pos,gene[[sample.id]])
  #table(pos==rownames(g))
  test <- as.data.frame(test)
  colnames(test) <- c("gene","bc")  # bin100kb_chr1:29555235-34534399|4979164
  #test$gene <- gsub("bin100kb_","",test$gene)
  test$chr <- unlist(lapply(strsplit(test$gene,":",fixed = T), function(x) x[1]))
  test$start <- as.numeric(unlist(lapply(strsplit(test$gene,":|-",perl = T), function(x) x[2])))
  test$end <- as.numeric(unlist(lapply(strsplit(test$gene,":|-|\\|",perl = T), function(x) x[3])))
  test <- test[,c("chr","start","end","bc")]
  test$bc <- as.numeric(test$bc)
  #str(test)
  
  #correct.gc <- PopSV::correct.GC(bc.f = test, gc.df = hg38)  # , outfile.prefix = "./pop"
  gc.df <- hg38
  bc.f <- test
  
  gc.class = cut(gc.df$GCcontent, breaks = seq(0, 1, 0.02), include.lowest = TRUE)
  samp.ii = unlist(tapply(1:length(gc.class), gc.class, function(e) e[sample(1:length(e),
                                                                             min(c(length(e), 500)))]))
  if(!is.data.frame(bc.f) & length(bc.f)==1 & is.character(bc.f)){
    bc.df = utils::read.table(bc.f, as.is = TRUE, header = TRUE)
  } else {
    bc.df = bc.f
  }
  
  if(any(bc.df$chr != gc.df$chr) | any(bc.df$start != gc.df$start)){
    bc.df = merge(bc.df, gc.df[,c('chr','start','end','GCcontent')])
    bc.df = bc.df[order(as.character(bc.df$chr), bc.df$start),]
  } else {
    bc.df$GCcontent = gc.df$GCcontent
  }
  
  lo = stats::loess(bc ~ GCcontent, data = bc.df[samp.ii, ])
  bc.df$bc = mean(bc.df$bc, na.rm = TRUE) * bc.df$bc/stats::predict(lo, newdata = bc.df)
  bc.df$bc = round(bc.df$bc, digits = 2)
  if (any(bc.df$bc < 0, na.rm = TRUE))
    bc.df$bc[bc.df$bc < 0] = 0
  bc.df$GCcontent = NULL
  colnames(bc.df)[4] <- sample.id
  return(bc.df)
}

#from PopSV
getGC <- function(bins.df, genome=BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19) {
  if(!all(c("chr","start","end") %in% colnames(bins.df))){
    stop("Missing column in 'bin.df'. 'chr', 'start' and 'end' are required.")
  }
  bins.df$chunk = rep(1:ceiling(nrow(bins.df)/1000), each = 1000)[1:nrow(bins.df)]
  addGC <- function(df) {
    if (!grepl("chr", df$chr[1])) {
      chrs = paste("chr", df$chr, sep = "")
    } else {
      chrs = df$chr
    }
    seq.l = Biostrings::getSeq(genome, chrs, df$start, df$end)
    lf = Biostrings::letterFrequency(seq.l, letters = c("G", "C"))
    df$GCcontent = rowSums(lf)/(df$end - df$start + 1)
    df[which(!is.na(df$GCcontent)), ]
  }
  ## bins.df = dplyr::do(dplyr::group_by(bins.df,chunk),addGC(.))
  chunk = . = NULL  ## Uglily appease R checks
  bins.df = dplyr::do(dplyr::group_by(bins.df, chunk), addGC(.))
  bins.df$chunk = NULL
  return(as.data.frame(bins.df))
}
  
message("read in ",bed)
#ref <- read.table("/BioII/lulab_b/baopengfei/shared_reference/hg38/gene.bed",header = F)
ref <- read.table(bed,header = F)
ref <- ref[,1:4]
colnames(ref) <- c("chr","start","end","ENSG")
# str(ref)

hg38 <- getGC(bins.df = ref, genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
hg38$pos <- paste0(hg38$chr,":",hg38$start,"-",hg38$end)
# table(duplicated(hg38$pos))  # 13 TURE
hg38$ID <- paste0(hg38$pos,":",hg38$ENSG)
# table(duplicated(hg38$ID))  # 0 TURE

message("read in ",matrix)
#gene <- read.table("./output/lulab/matrix/count_matrix_promoter.txt",header = T, row.names = 1, check.names = F)
gene <- read.table(matrix,header = T, row.names = 1, check.names = F)
g <- unlist(lapply(strsplit(rownames(gene),"|",fixed = T), function(x) x[1]))

dat <- lapply(colnames(gene),correct)
dat1 <- do.call("cbind",dat)
# table(duplicated(colnames(dat1)))
dat1 <- dat1[,!duplicated(colnames(dat1))]

dat1 <- dat1[,!(colnames(dat1) %in% c("chr","start","end"))]
dat1 <- round(dat1)   # keep interger
row.names(dat1) <- g
dat.l <- cbind(rownames(gene),dat1)
colnames(dat.l)[1] <- "gene_id"

#options ---
message("write to ",outfile)
#write.table(dat.l,"./output/lulab/matrix/count_matrix_promoter.correctGC.txt",sep = "\t",quote = F,row.names = F,col.names = T)
write.table(dat.l,outfile,sep = "\t",quote = F,row.names = F,col.names = T)

