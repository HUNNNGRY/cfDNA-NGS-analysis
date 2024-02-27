# in-house adapted NDR
## generate flank reference region bed file
### promoter150TSS50 --> 20001000TSS
cd /BioII/lulab_b/baopengfei/shared_reference/hg38
bedtools flank -g /BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/genome/chrom.size -s -l 2000 -r 0 -i  promoter150TSS50.bed > tmp1
bedtools flank -g /BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/genome/chrom.size -s -l 1000 -r 0 -i  promoter150TSS50.bed > tmp2
bedtools subtract  -g /BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/genome/chrom.size  -a  tmp1  -b tmp2  > tmp3
mv tmp3 /BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/ref/NDR/20001000TSS.bed
sort -k1,1 -k2,2n /BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/ref/NDR/20001000TSS.bed > /BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/ref/NDR/20001000TSS_sort.bed

### 2000TSS2000
bedtools slop -g /BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/genome/chrom.size -s -l 1851 -r 1951 -i  promoter150TSS50.bed >  promoter2000TSS2000.bed

### exon1end10002000 --> exon1end10002000
cd /BioII/lulab_b/baopengfei/shared_reference/hg38
bedtools flank -g /BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/genome/chrom.size -s -l 0 -r 2000 -i  promoter300100exon1end.bed > tmp1
bedtools flank -g /BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/genome/chrom.size -s -l 0 -r 1000 -i  promoter300100exon1end.bed > tmp2
bedtools subtract  -g /BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/genome/chrom.size  -a  tmp1  -b tmp2  > tmp3
mv tmp3 /BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/ref/NDR/exon1end10002000.bed
sort -k1,1 -k2,2n /BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/ref/NDR/exon1end10002000.bed > /BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/ref/NDR/exon1end10002000_sort.bed

### 2000exon1end2000
bedtools slop -g /BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/genome/chrom.size -s -l 1701 -r 2101 -i  promoter300100exon1end.bed >  promoter2000exon1end2000.bed

## lulab: cal dataset depth: TSS,exon1,ref_TSS,ref_exon1
cd  /BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq
cat  data/lulab/meta_data/sample_ids.txt | parallel -k -I % -j 10 "bedtools coverage -g genome/chrom.size -sorted -a ref/gtf/promoter150TSS50.bed -b output/lulab/bam-sorted-deduped/%.bam -mean >   output/lulab/NDR/depth/%_promoter150TSS50.bed ;
bedtools coverage -g genome/chrom.size -sorted -a ref/gtf/promoter300100exon1end.bed -b output/lulab/bam-sorted-deduped/%.bam -mean > output/lulab/NDR/depth/%_promoter300100exon1end.bed;
bedtools coverage -g genome/chrom.size -sorted -a  ref/NDR/exon1end10002000_sort.bed -b output/lulab/bam-sorted-deduped/%.bam -mean > output/lulab/NDR/depth/%_exon1end10002000.bed;
bedtools coverage -g genome/chrom.size -sorted -a ref/NDR/20001000TSS_sort.bed -b output/lulab/bam-sorted-deduped/%.bam -mean > output/lulab/NDR/depth/%_20001000TSS.bed"

## get matrix (draft)
sum-NOV_Rel_Dep.R
