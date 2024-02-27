#!/bin/bash
## get gene cnv from wisecondorx (ref:inputfile;gene.bed;outfile; outlog; sample_id; log2ratio or zscore )

if [ "$6" == "log2ratio" ] ; then
	cat $1 | grep -v chr | awk '{print "chr" $1 "\t" $2 "\t" $3 "\t" $1":"$2"-"$3  "\t"  $4}' | sort -k1,1 -k2,2n |  bedtools map -a $2 -b - -null NA -o mean  -c 5 | awk -v i=$5 'BEGIN {print "gene_id" "\t" i } {print $4 "\t" $7}' > $3 2> $4
elif [ "$6" == "zscore" ] ; then
	cat $1 | grep -v chr | awk '{print "chr" $1 "\t" $2 "\t" $3 "\t" $1":"$2"-"$3  "\t"  $5}' | sort -k1,1 -k2,2n |  bedtools map -a $2 -b - -null NA -o mean  -c 5 | awk -v i=$5 'BEGIN {print "gene_id" "\t" i } {print $4 "\t" $7}' > $3 2> $4
else
	echo "please choose log2ratio or zscore"
fi
