# run GSVA
# last 211028 by pengfei
# KEGG geneset online or offline mode 
# gene ID conversion refernece offline

options(stringsAsFactors = F)  # , digits = 6
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(GSVA))
suppressPackageStartupMessages(library(progress))
suppressPackageStartupMessages(library(clusterProfiler))

parser <- ArgumentParser(description='GSVA: get pathway/genesets matrix from TPM matrix')
parser$add_argument('-i', '--matrix', type='character', required=TRUE,
                    help='input matrix. Rows are non-duplicated genes(*ENSG* or hgnc_symbol). Columns are samples')
parser$add_argument('-d', '--distrubution', type='character', required=TRUE,
                    help='input matrix distrubution.  TPM (log transformed or not) using Gaussian ; counts using Poisson')
parser$add_argument('-g', '--geneset', type='character', required=TRUE,
                    help='geneset tab file : contain gene and geneset name (pathway/geneset)')
parser$add_argument('-m', '--method', type='character', default="ssgsea",
                    help='method=c("gsva", "ssgsea", "zscore", "plage"; "mean"), default:ssgsea, gsva is not that powerful for small samplesize')
parser$add_argument('-o', '--outfile', type='character', default="./run-gsva-Pathway.txt",
                    help='outfile name')
parser$add_argument('-p', '--cores', type='integer', default=0,
                    help='parallel cores, default: use all available')
#parser$add_argument('--online', dest="online", action='store_TRUE',
#                    help='whether to use online real-time updated KEGG, offline will use pre-saved PATH_ID_NAME_KEGGplusHallmark.txt file and -g is not necessary (default=FALSE)')
#parser.set_defaults(online=True)

args <- parser$parse_args()
matrix <- args$matrix
d <- args$distrubution
g <- args$geneset
method <- args$method
outfile <- args$outfile
cores <- args$cores
#online <- args$online
online <- FALSE
message(paste0("online: ",online))

# # test
# setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/")
# matrix <- "/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/CPM_matrix_gene.txt"
# g <- "/BioII/lulab_b/baopengfei/projects/methylation/ref/ensembl_entrez_symbol_pathway.txt"
# method <- "gsva"
# outfile <- "./gsva.txt"
# cores <-  12
# d <- "Gaussian"

## read input mat
message("reading input matrix: ",matrix)
exp <- read.table(matrix,sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)

#see if max < 50 in matrix file !
exp <- as.matrix(exp)
if(d=="Gaussian" & (max(exp) < 50)){
  message("input mat is log transformed, using Gaussian.")
} else if (d=="Gaussian" & (max(exp) >= 50)){
  message("input mat is not log transformed, converting to log and using Gaussian.")
  exp <- log2(exp+1)
} else if (d=="Poisson"){
  message("input mat is count, using Poisson.")
} else {
  message("unknown input mat")
}

#message("test2")
if(length(grep("ENSG",rownames(exp)))>=2)
{
message("input gene id is ENSG based, will transform to ENTREZID first")
## 1. convert rowname from *ENSG* to ENSG****
prefix <- paste("gene","promoter","promoter5k","promoter150TSS50","promoter300100exon1end","gene_exon1","CDS","UTR3","UTR5","CpG_island","LINE","SINE","retroposon",sep = "_|")
prefix <- paste0(prefix,"_")
rownames(exp) <- as.character(lapply(strsplit(rownames(exp),"\\|"), function(x) x[1]))
rownames(exp) <- sub(prefix,"",rownames(exp) )
dup <- duplicated(as.character(lapply(strsplit(rownames(exp),"\\."), function(x) x[1])))
exp <- exp[!dup,]
rownames(exp) <- as.character(lapply(strsplit(rownames(exp),"\\."), function(x) x[1]))

## 2. read in ref id convert file (ENSG**** and hgnc_symbol) and geneset file (contain pathway or geneset)
reference <-read.table("/BioII/lulab_b/baopengfei/shared_reference/geneset/encodev27_gene_v2.txt",sep = "\t",header = T,check.names = F)
#table(duplicated(reference$name)) # 1679 TRUE
#table(duplicated(reference$name[match(rownames(exp), reference$ensg)])) # 45 TRUE
reference <- reference[!duplicated(reference$ensg) & !duplicated(reference$ENTREZID),]

#convert rowname from ENSG**** to ENTREZID
exp <- exp[rownames(exp) %in% reference$ensg,]
rownames(exp) <- reference$ENTREZID[match(rownames(exp), reference$ensg)]  # many NA names but does not matter
} else 
{
  message("input gene id is considered ENTREZID based")
}


# 3. convert geneset/pathway df to list
if(online){message("online mode: KEGG2ENTREZ")
tmp <- download_KEGG(species="hsa")
tmp1 <- tmp$KEGGPATHID2EXTID  # online use entrez id in priority !!!
tmp2 <- tmp$KEGGPATHID2NAME
#ENTREZID        KEGGID  DESCRPTION
colnames(tmp1) <- c("KEGGID","ENTREZID")
colnames(tmp2) <- c("KEGGID","DESCRPTION")
geneset <- dplyr::left_join(tmp1,tmp2)
geneset <- geneset[,c("DESCRPTION","ENTREZID")]
} else
{message("offline mode: KEGG2ENTREZ")
geneset <- read.table(g,sep = "\t",header = T,check.names = F)
#a <- ifelse("short_DESCRPTION" %in% colnames(geneset),"short_DESCRPTION","DESCRPTION")  # USE short_DESCRPTION in priority 
a <- "DESCRPTION"                               # USE DESCRPTION in priority 
#geneset <- geneset[,c(a,"hgnc_symbol")]         # offline use hgnc_symbol 
geneset <- geneset[,c(a,"ENTREZID")]         # offline use ENTREZID 
}

colnames(geneset) <- c("geneset","gene")
geneset <- geneset[!duplicated(geneset),]

geneset <- geneset[geneset$gene!="" & geneset$geneset!="",]
geneset.l <-  base::split(geneset$gene,geneset$geneset)  

## 4. run gsva or mean
if(method!="mean"){
message("run gsva...")
re <- GSVA::gsva(exp, 
           geneset.l,
           method=method,  # method=c("gsva", "ssgsea", "zscore", "plage"), gsva is not that powerful for small samplesize 
           verbose=T,
           kcdf=d, # Gaussian: microarray fluorescent units in logarithmic scale, RNA-seq log-CPMs, log-RPKMs or log-TPMs
           parallel.sz=cores)  # mx.diff=F ?
#re <- t(re)
} else if(method=="mean"){
## mean func.
message("run mean...")
getaMeanPathwayMatFromTPMMat <- function(exp,PATH_ID_NAME){
  { 
    ##change sample name (contail group info)
    #rownames(TPM) <- sub(pattern = "Promoter_",replacement = "",x = rownames(TPM))
    log2TPM <- exp
    # colData <- colData[rownames(colData) %in% colnames(log2TPM),]
    # #groups <- unique(colData$group)  # Group or group !!!  without group selection
    # #groups <- groups
    # log2TPM <- log2TPM[,match(rownames(colData),colnames(log2TPM))]
    # #all(colnames(log2TPM) == rownames(colData))
    # 
    # colData$group_sample_id <- NA
    # colData$group_sample_id <- rownames(colData)   # already has group in data id
    # colData <- colData[rownames(colData) %in% colnames(log2TPM),]
    # colnames(log2TPM) <- colData$group_sample_id   
  }
  
  # set progress
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent eta: :eta",
    total = length(unique(PATH_ID_NAME$geneset)), clear = FALSE, width= 60) 
  
  #get all genes in each pathway in KEGG
  n=1
  pathway_count={}
  while(n <= length(unique(PATH_ID_NAME$geneset))){
    pathway_name <- as.character(unique(PATH_ID_NAME$geneset)[n])
    pathway <- subset(PATH_ID_NAME,geneset==pathway_name)
    candidate <- unique(pathway$gene)  
  
    #summary one pathway total count in each sample and make 1 row dataframe
    j=1
    pathway_gene_count={}
    while(j<=length(candidate)){
      target <- candidate[j]
      if(length(grep(target,rownames(log2TPM)))==0) {
        #print(paste0("No ",target," in this dataset."))
        j=j+1
      } else {
        temp <- log2TPM[grep(target,rownames(log2TPM),fixed=TRUE),]
        pathway_gene_count <- rbind(pathway_gene_count,temp)
        j=j+1
      }
    }
    if(is.null(pathway_gene_count)){
      n=n+1
      next;
    } else {
      one_pathway <- as.data.frame(t(colMeans(pathway_gene_count)))   ## colSums: sum up all log2TPM ( colMeans ...)
      rownames(one_pathway) <- pathway_name
      pathway_count <- rbind(pathway_count,one_pathway)
      n=n+1
    }
    pb$tick()
    Sys.sleep(1 / 100)
  }
  return(pathway_count)
}  
re <- getaMeanPathwayMatFromTPMMat(exp,geneset)
} else {
message("enrich/integration method not supported, please use (gsva,ssgsea,zscore, plage; or mean)")
}

## 4. write output
message("writing outfile: ",outfile)
write.table(re,outfile,sep = "\t",quote = F,col.names = T,row.names = T)
message("done!")
