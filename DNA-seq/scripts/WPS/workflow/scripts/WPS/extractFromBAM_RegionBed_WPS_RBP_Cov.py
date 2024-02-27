#!/usr/bin/env python

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *03.06.2014
"""

import sys, os
import argparse
import gzip
import pysam
import random
import math
import numpy as np
import pandas as pd

from collections import defaultdict
from bx.intervals.intersection import Intersecter, Interval

def isSoftClipped(cigar):
  #Op BAM Description
  #M 0 alignment match (can be a sequence match or mismatch)
  #I 1 insertion to the reference
  #D 2 deletion from the reference
  #N 3 skipped region from the reference
  #S 4 soft clipping (clipped sequences present in SEQ)
  #H 5 hard clipping (clipped sequences NOT present in SEQ)
  #P 6 padding (silent deletion from padded reference)  ???
  #= 7 sequence match
  #X 8 sequence mismatch
  for (op,count) in cigar:
    if op in [4,5,6]: return True
  return False

def aln_length(cigarlist):
  tlength = 0
  for operation,length in cigarlist:
    if operation == 0 or operation == 2 or operation == 3 or operation >= 6: tlength += length
  return tlength

def regionFileParser(infile):
  for line in infile:
    ########
    # implement proper bedfile reading
    ########
    chrom,start,end,cid,score,strand = line.split() # positions should be 0-based and end non-inclusive
    if chrom.startswith("chr"):
      chrom = chrom.replace("chr","")
    yield chrom, start, end, cid, score, strand
  return 

def parseRegion(regionstr):
  try:
    chrom=regionstr.split(":")[0]
    start,end=(regionstr.split(":")[-1]).split("-")
    start=str(int(start)-1)
    if chrom.startswith("chr"):
      chrom = chrom.replace("chr","")
    yield chrom, start, end, regionstr, "0", "+"
    return
  except:
    sys.stderr.write("Failed parsing coordinate string %s\n"%(regionstr))
  return


parser = argparse.ArgumentParser()
parser.add_argument('files', nargs='+', help='BAM files with aligned reads.')
parser.add_argument("-i","--input", dest="input", help="Use regions transcript file (def transcriptAnno.tsv)",default="transcriptAnno.tsv")
#parser.add_argument("-l","--lengthSR", dest="lengthSR", help="Length of full reads (default 76)",default=76,type=int)
#parser.add_argument("-l","--lengthSR", dest="lengthSR", help="Length of full reads (default 150)",default=150,type=int)
# parser.add_argument("-l","--lengthSR", dest="lengthSR", help="Length of SE full reads (default 36)",default=36,type=int)
parser.add_argument("-l","--lengthSR", dest="lengthSR", help="Length of SE full reads (RNA sm default 20)",default=20,type=int)
parser.add_argument("-m","--merged", dest="merged", help="Assume reads are merged (default Off)",default=False,action="store_true")
#parser.add_argument("-t","--trimmed", dest="trimmed", help="Assume reads are trimmed (default Off)",default=False,action="store_true")
# parser.add_argument("-w","--protection", dest="protection", help="Base pair protection window size (default 120)",default=120,type=int)
parser.add_argument("-w","--protection", dest="protection", help="Base pair protection window size (RNA sm default 15), too long will leads to signal shift (towards right)",default=15,type=int)
#parser.add_argument("-c","--method", dest="method", help="Type of protection score to calculate. The default is: %(default)s and your choices are: %(choices)s.",default="WPSv1",type=str, choices=("WPSv1","WPSv2","WPSv3","WPSv4","WPSv5"))
parser.add_argument("-o","--outfile", dest="outfile", help="Outfile prefix (def 'block_%%s.tsv.gz')",default='block_%s.tsv.gz') # reserve atleast 6 digits 
parser.add_argument("-e","--empty", dest="empty", help="Keep files of empty blocks (def Off)",default=False,action="store_true")
parser.add_argument("--minInsert", dest="minInsSize", help="Minimum read length threshold to consider (def None)",default=-1,type=int) # template_length
parser.add_argument("--maxInsert", dest="maxInsSize", help="Maxmum read length threshold to consider (def None)",default=-1,type=int) # template_length
parser.add_argument("--max_length", dest="max_length", help="Assumed maximum insert size (default 1000)",default=1000,type=int) # not used
parser.add_argument("--downsample", dest="downsample", help="Ratio to down sample reads (default 0.5)",default=0.5,type=float)
parser.add_argument("--onefile", dest="onefile", help="Print as single (WPS only ?) output to stdout (default OFF)",default=False,action="store_true")
parser.add_argument("-v","--verbose", dest="verbose", help="Turn debug output on",default=False,action="store_true")
parser.add_argument("--libType", dest="libType", help="RNA library stranded type, forward: read1 is forward strand of original RNA template; reverse: read1 is reverse strand of original RNA template; unstranded: TCGA-mRNA-seq?. SE only has 'forward' libtype (def: reverse)",default="reverse")


options = parser.parse_args()

minInsSize,maxInsSize = None,None
if options.minInsSize > 0 and options.maxInsSize > 0 and options.minInsSize < options.maxInsSize:
  minInsSize = options.minInsSize
  maxInsSize = options.maxInsSize
  if options.verbose: sys.stderr.write("Using min/max length cutoffs: %d/%d\n"%(minInsSize,maxInsSize))
  
options.outfile = options.outfile.strip("""\'""")
outfiles = {}
if not options.onefile:
  outfiles = {'WPS':gzip.open(options.outfile%"WPS",'wt'), 'WPS_v2':gzip.open(options.outfile%"WPS_v2",'wt'), 'COV':gzip.open(options.outfile%"COV",'wt'), 'STARTS':gzip.open(options.outfile%"STARTS",'wt'), 'length':gzip.open(options.outfile%"length",'wt')}

if options.verbose:
    sys.stderr.write("library type: %s\n"%(options.libType))
    sys.stderr.write("Downsample Ratio: %d\n"%(options.downsample))
protection = options.protection//2
libType = options.libType
downsample = options.downsample
merged = options.merged
lengthSR = options.lengthSR
#validChroms = set(map(str,range(1,23)+["X","Y"]))
#validChroms = set(map(str,range(1,23)+["X","Y"]))
#validChroms = [str(i) for i in range(1, 23)] + ["X","Y"]
#validChroms = [str(i) for i in range(1, 23)] + ["X","Y","MT"]  # need remove for tx map

if os.path.exists(options.input):
  if ".gz" in options.input:
    infile = gzip.open(options.input,"r")
  else:
    infile = open(options.input)
  regionIterator = regionFileParser(infile)
else:
  regionIterator = parseRegion(options.input)

# add insert length
regions = pd.read_csv(options.input, header=None)
region_num = regions.shape[0]
regions_length = pd.DataFrame(data=np.zeros((maxInsSize,region_num)),index=range(maxInsSize),dtype=np.int)  # ,columns=regions.iloc[:,0]
current_regionNum = 0

for chrom,start,end,cid,score,strand in regionIterator:
  # if chrom not in validChroms: continue # need remove for tx map
  #regionStart,regionEnd = int(start)-300,int(end)+300
  regionStart,regionEnd = int(start),int(end)
  #if regionStart < 1: continue
  l_bound=regionStart-protection-1
  r_bound=regionEnd+protection+1
  #outchrom = chrom.replace("chr","")
  #if outchrom.startswith('gi|'):
    #NCfield = outchrom.split("|")[-2]
    #if NCfield.startswith("NC"):
      #outchrom = "%d"%(int(NCfield.split(".")[0].split("_")[-1]))
  if options.verbose:
    sys.stderr.write("Processing region: %s:%s-%s %s %s %s\n"%(chrom,start,end,cid,score,strand))
  
  posRange = defaultdict(lambda:[0,0])
  filteredReads = Intersecter()

  for bamfile in options.files: # treat all bam files as one sample, i.e. merge all reads
    if options.verbose: sys.stderr.write("Reading %s\n"%bamfile)
    bamfile = bamfile.strip("""\'""")
    if os.path.exists(bamfile) and (os.path.exists(bamfile.replace(".bam",".bai")) or os.path.exists(bamfile+".bai")):
      input_file = pysam.Samfile( bamfile, "rb" )
      prefix = ""
      for tchrom in input_file.references:
        if tchrom.startswith("chr"):
          prefix = "chr"
          break
      if options.verbose: sys.stderr.write("Retrieving reads...\n")
      if chrom == "MT":
        chrom="M"
      for read in input_file.fetch(prefix+chrom,max(regionStart-protection-1,0),regionEnd+protection+1): # fetch returns a AlignedSegment object which represents a single read along with its fields and optional tags
        #default filter duplicates reads
        # if read.is_duplicate or read.is_qcfail or read.is_unmapped: continue
        #default filter softclip reads
        # if isSoftClipped(read.cigar): continue
        if read.is_paired:
          #default only keep proper paired reads
          # if read.mate_is_unmapped: continue
          # if read.rnext != read.tid: continue          
          
          #only add specific stranded reads, if not, then skip this read. (2204,bpf)
          # if (strand=="-"): 
          #   if(libType == "reverse"): 
          #     if (read.is_reverse):
          #       if (read.is_read1): continue
          #     else:
          #       if (read.is_read2): continue 
          #   if(libType == "forward"): 
          #     if (read.is_reverse):
          #       if (read.is_read2): continue
          #     else:
          #       if (read.is_read1): continue
          #if (strand=="+"): 
          if(libType == "reverse"):
            if (read.is_reverse and read.is_read2):continue
            elif ((not read.is_reverse) and read.is_read1): continue 
          elif(libType == "forward"):
            if (read.is_reverse and read.is_read1): continue
            elif ((not read.is_reverse) and read.is_read2): continue  

          # if (strand=="-" and libType == "reverse" and read.is_read1 and read.is_reverse): continue
          # if (strand=="-" and libType == "reverse" and read.is_read2 and (not read.is_reverse)): continue
          # if (strand=="-" and libType == "forward" and read.is_read1 and (not read.is_reverse)): continue
          # if (strand=="-" and libType == "forward" and read.is_read2 and read.is_reverse): continue
          # if (strand=="+" and libType == "reverse" and read.is_read1 and (not read.is_reverse)): continue
          # if (strand=="+" and libType == "reverse" and read.is_read2 and read.is_reverse): continue
          # if (strand=="+" and libType == "forward" and read.is_read1 and read.is_reverse): continue
          # if (strand=="+" and libType == "forward" and read.is_read2 and (not read.is_reverse)): continue
          #if (strand=="*" or strand=="."):
          
          if read.is_read1 or (read.is_read2 and ((read.pnext+read.qlen < l_bound) or (read.pnext > r_bound)) ): # approximate mate1,2 have similar qlen for faster calculation
          # if read.is_read1 or (read.is_read2 and ((read.pnext+aln_length(input_file.mate(read).cigar) < regionStart-protection-1) or (read.pnext > regionEnd+protection+1)) ): # avoid count PE twice: only consider R1, if R1 not in region then consider R2. Note: aln_length(read.mate().cigar) better than read.mate().qlen, latter one may has insertion/deletion...  (2204,bpf)
            if read.isize == 0: continue
            #default downsample (need after PE each frag/insert selected)
            if downsample != None and random.random() >= downsample: continue #  seems will not decrease run time
            rstart = min(read.pos,read.pnext)+1 # 1-based
            lseq = abs(read.isize)
            rend = rstart+lseq-1 # end included
            if minInsSize != None and ((lseq < minInsSize) or (lseq > maxInsSize)): continue
            
            #if options.verbose: sys.stderr.write("Adding interval for read...\n")
            filteredReads.add_interval(Interval(rstart,rend))
            
            #add lenth freq count (only count original region overlapped reads)
            if (rstart <= regionEnd) and (rend >= regionStart):
              regions_length.iloc[lseq-1, current_regionNum] = regions_length.iloc[lseq-1, current_regionNum] + 1
            
            #print read.qname,rstart,rend,rend-rstart,abs(read.isize)
            for i in range(rstart,rend+1):
              if i >= regionStart and i <= regionEnd:
                posRange[i][0]+=1
            if rstart >= regionStart and rstart <= regionEnd:
              posRange[rstart][1]+=1
            if rend >= regionStart and rend <= regionEnd:
              posRange[rend][1]+=1
        else: # for SE bam input 
          #default downsample (need after PE selected)
          if downsample != None and random.random() >= downsample: continue
          # if libType == "reverse":
          #   # if options.verbose: sys.stderr.write("Warning: consider SE only has 'forward' options.libType.\n")
          #   sys.stderr.write("Warning: consider SE only has 'forward' options.libType.\n")
          #   libType == "forward"
          rstart = read.pos+1 # 1-based
          lseq = aln_length(read.cigar)
          lseq_ext = max(lengthSR,lseq)
          rend = rstart+lseq_ext-1 # end included
          if minInsSize != None and ((lseq_ext < minInsSize) or (lseq_ext > maxInsSize)): continue
          
          filteredReads.add_interval(Interval(rstart,rend))
          #print read.qname,rstart,rend,rend-rstart,aln_length(read.cigar)
          for i in range(rstart,rend+1):
            if i >= regionStart and i <= regionEnd:
              posRange[i][0]+=1
          
          #add lenth freq count (only count original region overlapped reads)
          if (rstart <= regionEnd) and (rend >= regionStart):
            regions_length.iloc[lseq-1, current_regionNum] += 1 # use original read length
            
          #if ((options.merged or read.qname.startswith('M_')) or ((options.trimmed or read.qname.startswith('T_')) and read.qlen <= options.lengthSR-10)):

          # original code below only consider SE 5' end exact (3' end not exact)
          #consider reserve or forward stranded RNA bam, SE only has "forward" options.libType (2204,bpf)
          if merged:
            # if (strand=="+" and read.is_reverse): continue #  or (strand=="-" and (not read.is_reverse))
            if (rstart >= regionStart and rstart <= regionEnd):
              posRange[rstart][1]+=1
            if (rend >= regionStart and rend <= regionEnd):
              posRange[rend][1]+=1
          # elif (read.is_reverse):
          #   #only add specific stranded reads, if not, then skip this read. (2204,bpf)
          #   if (strand=="+"): continue
          #   if rend >= regionStart and rend <= regionEnd:
          #     posRange[rend][1]+=1
          elif ((not read.is_reverse)):
            # if (strand=="-"): continue
            if (rstart >= regionStart and rstart <= regionEnd):
              posRange[rstart][1]+=1
    else:
      sys.stderr.write("File without BAM index skipping: %s\n"%(bamfile))
      # if options.verbose: sys.stderr.write("Warning: without input reads, output files will be populated with zero/empty values.\n")
      sys.stderr.write("Warning: without input reads, output files will be populated with zero/empty values.\n")


  if options.verbose: sys.stderr.write("Evaluating posRange vector...\n")
  #filename = options.outfile%cid
  #outfile = gzip.open(filename,'w')
  cov_sites = 0
  outLines = []
  wps_list = []
  wpsv2_list = []
  cov_list = []
  starts_list = []

  for pos in range(regionStart,regionEnd+1):
    rstart,rend = pos-protection,pos+protection
    gcount,bcount,ecount = 0.0,0.0,0.0
    for read in filteredReads.find(rstart,rend):
      if (read.start > rstart):
        bcount += 1.0
      #elif below means only count start as ends  
      elif (read.end < rend): 
        ecount += 1.0
      else: 
        gcount += 1.0
    #startCount seems contains both start and end read counts
    covCount,startCount = posRange[pos]
    
    #cov_sites not write to outfile ?
    cov_sites += covCount
    # wpsValue = gcount-(bcount+ecount)
    wpsValue = gcount - (bcount + ecount) # 
    # wpsValue = gcount/(bcount+ecount)
    wpsValue2 = (gcount + 1) / (bcount + ecount + 1) # pseudo count : 0.1
    #if (options.method != "WPSv1") and (2*gcount+bcount+ecount > 1):
    #  if options.method == "WPSv2":
    #    wpsValue = 2.0*gcount/(2.0*gcount+bcount+ecount)
    #  elif options.method == "WPSv3":
    #    wpsValue = 2.0*gcount/(2.0*gcount+bcount+ecount)
    #    wpsValue = wpsValue * (1.0 - math.sqrt( (wpsValue*(1.0-wpsValue)) / (2.0*gcount+bcount+ecount-1.0) ))
    #  elif options.method == "WPSv4":
    #    wpsValue = max(2.0*gcount-abs(bcount-ecount),0)/(2.0*gcount+bcount+ecount)
    #  elif options.method == "WPSv5":
    #    wpsValue = max(2.0*gcount-abs(bcount-ecount),0)/(2.0*gcount+bcount+ecount)
    #    wpsValue = wpsValue * (1.0 - math.sqrt( (wpsValue*(1.0-wpsValue)) / (2.0*gcount+bcount+ecount-1.0) ))
    #elif (options.method != "WPSv1"): 
    #  wpsValue = 0.0
    if options.onefile:
      outLines.append("%s\t%d\t%.4f\t%.4f\t%.4f\n"%(chrom,pos,covCount,startCount,wpsValue))
    else: 
      wps_list.append(wpsValue)
      wpsv2_list.append(wpsValue2)
      cov_list.append(covCount)
      starts_list.append(startCount)

  #alway keep direction in 5'-->3' (left-->right) in final matrix
  if options.onefile:
    if strand == "-": outLines = outLines[::-1]
    for line in outLines: sys.stdout.write(line)
  else:
    if strand == "-": wps_list = wps_list[::-1]
    if strand == "-": wpsv2_list = wpsv2_list[::-1]
    if strand == "-": cov_list = cov_list[::-1]
    if strand == "-": starts_list = starts_list[::-1]
    #span - (start + end)
    outfiles['WPS'].write(cid+","+",".join(map(lambda x:str(round(x,5)).replace(".0",""),wps_list))+"\n")
    #span / (start + end)
    outfiles['WPS_v2'].write(cid+","+",".join(map(lambda x:str(round(x,5)).replace(".0",""),wpsv2_list))+"\n")
    #total depth (span + start + end), pre-computed
    outfiles['COV'].write(cid+","+",".join(map(lambda x:str(round(x,5)).replace(".0",""),cov_list))+"\n")
    #depth (start + end), pre-computed
    outfiles['STARTS'].write(cid+","+",".join(map(lambda x:str(round(x,5)).replace(".0",""),starts_list))+"\n")
  current_regionNum += 1

#write insert length (one bam as input)     
regions_length.to_csv(outfiles['length'])

if not options.onefile:
  for name,filestream in outfiles.items():
    filestream.close()

