import argparse,sys,os,pysam,pyBigWig 
import numpy as np 
import pandas as pd
from collections import defaultdict
from bx.intervals.intersection import Intersecter, Interval

parser = argparse.ArgumentParser()
parser.add_argument("-i","--inputBam", dest="inputBam", help="Input bam (align file) path.")
# parser.add_argument("-l","--lengthSR", dest="lengthSR", help="Length of SE full reads (RNA sm default 20)",default=20,type=int)
# parser.add_argument("-m","--merged", dest="merged", help="Assume reads are merged (default Off)",default=False,action="store_true")
parser.add_argument("-o","--outputBigWig", dest="outputBigWig", help="Output end depth bigwig path",default='readEnd.bigWig') # reserve atleast 6 digits 
parser.add_argument("--minInsert", dest="minInsSize", help="Minimum read length threshold to consider (def None)",default=-1,type=int) # template_length
parser.add_argument("--maxInsert", dest="maxInsSize", help="Maxmum read length threshold to consider (def None)",default=-1,type=int) # template_length
# parser.add_argument("--max_length", dest="max_length", help="Assumed maximum insert size (default 1000)",default=1000,type=int) # not used
# parser.add_argument("--downsample", dest="downsample", help="Ratio to down sample reads (default 1)",default=1,type=float)
# parser.add_argument("-v","--verbose", dest="verbose", help="Turn debug output on",default=False,action="store_true")
# parser.add_argument("--libType", dest="libType", help="RNA library stranded type, forward: read1 is forward strand of original RNA template; reverse: read1 is reverse strand of original RNA template. SE only has 'forward' libtype (def: reverse)",default="reverse")
options = parser.parse_args()

minInsSize,maxInsSize = None,None
if options.minInsSize > 0 and options.maxInsSize > 0 and options.minInsSize < options.maxInsSize:
    minInsSize = options.minInsSize
    maxInsSize = options.maxInsSize
    if options.verbose: sys.stderr.write("Using min/max length cutoffs: %d/%d\n"%(minInsSize,maxInsSize))

def aln_length(cigarlist):
    tlength = 0
    for operation,length in cigarlist:
        if operation == 0 or operation == 2 or operation == 3 or operation >= 6: tlength += length
    return tlength


# 打开bam文件和bigwig文件
bamfile=options.inputBam
bwfile=options.outputBigWig

if os.path.exists(bamfile) and (os.path.exists(bamfile.replace(".bam",".bai")) or os.path.exists(bamfile+".bai")):
    # bam_file = pysam.Samfile( bamfile, "rb" )
    bam_file = pysam.AlignmentFile(bamfile, "rb") 
    bw_file = pyBigWig.open(bwfile, "w")

    # 设置bigwig文件的元数据
    chroms = bam_file.references 
    chrom_lengths = bam_file.lengths 
    chrom_sizes = [(c, l) for c, l in zip(chroms, chrom_lengths)] 
    bw_file.addHeader(list(chrom_sizes))
    
    for chrom in chroms: 
        prefix = ""
        if chrom.startswith("chr"):
            prefix = "chr"
            break
        if chrom == "MT":
            chrom="M"
        chrom_length = chrom_lengths[bam_file.get_tid(chrom)]
        chrom_positions = np.array([pos for pos in range(chrom_length)], dtype=np.uint32)
        chrom_values = np.zeros(chrom_length, dtype=np.uint32)  
        
        
        for read in bam_file.fetch(prefix+chrom,0,chrom_length): # fetch returns a AlignedSegment object which represents a single read along with its fields and optional tags
            rstart = read.pos+1 # 1-based
            lseq = aln_length(read.cigar)
            # lseq_ext = max(lengthSR,lseq)
            lseq_ext = lseq
            rend = rstart+lseq_ext-1 # end included
            if minInsSize != None and ((lseq_ext < minInsSize) or (lseq_ext > maxInsSize)): continue
            # original code below only consider SE 5' end exact (3' end not exact)
            #consider reserve or forward stranded RNA bam, SE only has "forward" options.libType (2204,bpf)
            if (rstart >= 0 and rstart <= chrom_length):
                chrom_values[rstart-1]+=1 # [0]
            if (rend >= 0 and rend <= chrom_length):
                chrom_values[rend-1]+=1 # [0]
        
        bw_file.addEntries(chrom, chrom_positions, values=chrom_values)
        
    # 关闭文件
    bam_file.close() 
    bw_file.close()
else:
    sys.stderr.write("not valid input bam\n")


# Traceback (most recent call last):
#   File "/BioII/lulab_b/baopengfei/gitsoft/cfDNA-2016/workflow/scripts/WPS/bam2readEndCov.py", line 76, in <module>
#     bw_file.addEntries(chrom, chrom_positions, values=chrom_values)
# RuntimeError: You must provide a valid set of entries. These can be comprised of any of the following


# V1
# import pysam
# import pyBigWig

# # 打开bam文件和bigwig文件
# bam_file = pysam.AlignmentFile("your_bam_file.bam", "rb")
# bw_file = pyBigWig.open("your_bigwig_file.bw", "w")

# # 设置bigwig文件的元数据
# chroms = bam_file.references
# chrom_lengths = bam_file.lengths
# chrom_sizes = [(c, l) for c, l in zip(chroms, chrom_lengths)]
# bw_file.addHeader(list(chrom_sizes))

# # 遍历bam文件中的每个read，提取其起始和终止位点计数
# counts = {}
# for read in bam_file:
#     if read.is_unmapped:
#         continue
#     start, end = read.reference_start, read.reference_end
#     if start not in counts:
#         counts[start] = 0
#     if end not in counts:
#         counts[end] = 0
#     counts[start] += 1
#     counts[end] += 1

# # 将起始和终止位点计数写入bigwig文件
# for chrom in chroms:
#     chrom_counts = [(pos, count) for pos, count in counts.items() if bam_file.get_reference_name(pos) == chrom]
#     chrom_counts.sort()
#     chrom_values = [count for pos, count in chrom_counts]
#     chrom_positions = [pos for pos, count in chrom_counts]
#     bw_file.addEntries(chrom, chrom_positions, values=chrom_values)

# # 关闭文件
# bam_file.close()
# bw_file.close()