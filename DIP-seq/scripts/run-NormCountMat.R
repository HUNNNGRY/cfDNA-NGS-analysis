# norm count file
# last 220218 by pengfei
# b.p.f@qq.com
# 2023.02.08: update externalLibSize and geneLength
# todo: standard parse params

suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(argparse))
#suppressPackageStartupMessages(library(EDASeq))

parser <- ArgumentParser(description='normalize count matrix using edgeR func. (lib size + TMM factor)')
parser$add_argument('-i', '--matrix', type='character', required=TRUE,
                    help='input count matrix. Rows are genes. Columns are samples.')
parser$add_argument('-s', '--sample', type='character', default="all",
                    help='tsv/txt file: pass filter colnames/samples, default: all (keep all colnames/samples), optional: txt path of valid samples, will change mat to the same order')
parser$add_argument('-g', '--group', type='character', default="FALSE",
                    help='group info of colnames/samples, default: FALSE, (consider as one group), optional: txt path of group, in same order of valid samples')
parser$add_argument('-m', '--method', type='character', default="TMM",
                    help='normalize method: "TMM","TMMwsp","RLE","upperquartile","none", default: TMM')
parser$add_argument('-f', '--filterByExpr', type='character', default="FALSE",
                    help='remove low expr records/genes, default: FALSE, optional: TRUE')
parser$add_argument('-p', '--priCount', type='integer', default=2,
                    help='priCount: int var for cpm(), default: 2')
parser$add_argument('--geneLength', type='character', default="FALSE",
                    help='if provide geneLength txt path, then cal rpkm mat as output instead of cpm, default: FALSE, optional: txt path of geneLength, in same order of valid genes')
parser$add_argument('-e', '--externalLibSize', type='character', default="FALSE",
                    help='to build DEGlist replace colSums, default: FALSE, optional: txt path of libsize, in same order of valid samples')
parser$add_argument('-l', '--log', type='character', default="FALSE",
                    help='whether to log cpm result, default: FALSE, optional: TRUE')
parser$add_argument('-o', '--outfile', type='character', required=TRUE,
                    help='output file path')
args <- parser$parse_args()

for(i in 1:length(args)){
  v.name <- names(args)[i]
  assign(v.name,args[[i]])
  message(paste0(v.name,": ",get(v.name)))
}

ma <- args$matrix
sa <- (args$sample)
gr <- (args$group)
fi <- toupper(args$filterByExpr)
me <- args$method
pr <- args$priCount
gl <- (args$geneLength)
ex <- (args$externalLibSize)
lo <- args$log
ou <- args$outfile


# # test multi-omics (old)
# {
# setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq")
# ma <- "output/lulab/matrix/archive/count_matrix_gene.txt"
# sa <- "all"
# gr <- "FALSE"
# fi <- FALSE
# me <- "TMM"
# lo <- FALSE
# ou <- "output/lulab/matrix/archive/CPM_matrix_gene.txt1"
# }

## test exPeak (new)
# {
# setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/")
# ma <- "output/GSE71008_NCpool/call_peak_all/count_matrix/expeak_b5_d50_p1.txt"
# sa <- "all"
# gr <- "FALSE"
# fi <- FALSE
# me <- "TMM"
# pr <- 1
# gl <- "output/GSE71008_NCpool/call_peak_all/count_matrix/expeak_b5_d50_p1.peakLength"
# ex <- "output/GSE71008_NCpool/call_peak_all/count_matrix/libSize"
# lo <- "FALSE"
# ou <- "output/GSE71008_NCpool/call_peak_all/count_matrix/expeak_b5_d50_p1_RPKM.txt"
# }

## test piranha (new)
# {
#   setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/")
#   ma <- "output/GSE71008/call_peak_all/count_matrix/piranha_b5_p01.txt"
#   sa <- "all"
#   gr <- "tmp/tmpGroup" # "FALSE"
#   fi <- FALSE
#   me <- "TMM"
#   pr <- 1
#   gl <- "output/GSE71008/call_peak_all/count_matrix/piranha_b5_p01.peakLength"
#   ex <- "FALSE" #"output/GSE71008/call_peak_all/count_matrix/libSize"
#   lo <- "FALSE"
#   ou <- "output/GSE71008/call_peak_all/count_matrix/piranha_b5_p01_RPKM.txt"
# }


message(paste0("read input mat from ",ma))
count.matrix <- read.csv(ma,sep = "\t",header = TRUE, check.names = F,row.names = 1)
#count.matrix[1:5,1:5]

if(sa=="all"){
  message(paste0("keep all samples from input mat"))
  s <- colnames(count.matrix)
} else {
  message(paste0("only keep samples from path: ",sa))
  s <- read.delim(sa,header = F,stringsAsFactors = F)$V1
}
#count.matrix <- count.matrix[,colnames(count.matrix) %in% s,,drop=F] # do not change col order
count.matrix <- count.matrix[,s,drop=F] # change col order

if(gr=="FALSE"){
  # message("using group info from sample prefix")
  # g <- unlist(lapply(strsplit(s,"-|_|\\."),function(x) x[1]))
  message("not using group info (all as 1)")
  g <- rep(1,length(s))
} else {
  message(paste0("using group info from path: ",gr))
  g <- read.delim(gr,header = F,stringsAsFactors = F)$V1
}
tissue.types <- factor(g,level=unique(g)) 

# put normal at reference level ????
if (ex=="FALSE"){
  message("use rowSums as lib.size in DGEList")
  y <- edgeR::DGEList(counts=count.matrix,group=tissue.types)  # def: lib.size=colSums()
}else{
  message("use provided path as lib.size in DGEList")
  libSize <- read.delim(ex,header = F,stringsAsFactors = F)$V1
  y <- edgeR::DGEList(counts=count.matrix,group=tissue.types,lib.size=libSize)
}
print(paste0("dim: ",dim(y$counts)))
print(paste0("grp: ",table(tissue.types)))

# remove low expressed gene
if(fi=="TRUE"){
  message("filter low expr by group info")
  keep <- edgeR::filterByExpr(y,group = tissue.types)
  if (ex=="FALSE"){
    message("libSize not kept !")
    y <- y[keep, , keep.lib.sizes=FALSE]   # gene num decrease (may set TRUE； https://support.bioconductor.org/p/9143866/)
  } else{
    message("libSize is kept !")
    y <- y[keep, , keep.lib.sizes=TRUE]   # keep.lib.sizes=FALSE, the lib.size for each sample (cf. the y$samples data.frame) will be recalculated to be the sum of the counts left
  }
}else if(fi=="FALSE"){
  message("do not filter low expr")
}
#EDASeq::plotRLE(edgeR::cpm(y))


# calculate scaling factor for library size
message(paste0("cal scaling factor using ",me))
y <- edgeR::calcNormFactors(y, method=me) # method="RLE" etc. is also applicable
if(gl=='FALSE'){
  message(paste0("gene length not provided, will output cpm mat"))
  y.cpm <- as.data.frame(edgeR::cpm(y,normalized.lib.sizes = TRUE,
                                    log = lo, 
                                    prior.count = pr) 
                                    # group=y$sample$group
                         )
  # gene.length=null
} else {
  message(paste0("gene length provided, will output rpkm mat"))
  geneLen <- read.delim(gl,header = F,stringsAsFactors = F)$V1
  y.cpm <- as.data.frame(edgeR::rpkm(y,normalized.lib.sizes = TRUE,
                                     log = lo, 
                                     prior.count = pr, 
                                     gene.length = geneLen)
                                     # group=y$sample$group
                         )
}
# y$counts[1:3,1:3]
# head(colSums(y$counts))
# head(y$samples)
# head(colSums(y.cpm))


#note: In the edgeR package, RPKM values are usually calculated using the cpm function. This is because the cpm function (counts per million) is more flexible and can calculate different types of expression standardization, while the rpkm function (reads per kilobase of transcript per million reads) can only calculate RPKM values.
#normalized.lib.sizes represents the normalized sequencing depth, while lib.sizes represents the raw sequencing depth.
#cpmByGroup/rpkmByGroup will use group param and produce merged sample by group, but in cpm/rpkm group param seem not needed 

#EDASeq::plotRLE(edgeR::cpm(y))
#limma::plotMDS(y.cpm,labels=y$samples$group)
#limma::plotMDS(y.cpm,labels=y$sample$group,color=y$sample$group)

y.cpm <- cbind(rownames(y.cpm),y.cpm)
colnames(y.cpm)[1] <- "gene_id"

# output
message(paste0("write file to ",ou))
write.table(y.cpm, file = ou, sep = "\t",col.names = T,row.names = F,quote = F)
