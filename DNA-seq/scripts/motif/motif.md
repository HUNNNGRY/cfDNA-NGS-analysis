## get 5 end motif
cd  /BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq
mkdir output/lulab/motif5
for i in `cat meta/lulab/sample_ids.txt`
do 
echo "start $i at `date`"
python /BioII/lulab_b/baopengfei/gitsoft/cfdna-main/cfdna_practical.py motif-analysis -b output/lulab/04bam_dedup/${i}.bam -s  False -m 4 | head -n1 | cut -d "{" -f 4 | cut -d "}" -f 1 | tr "," "\n" | tr -d " " > output/lulab/motif5/${i}_motif5.txt
done

## join outfile to matrix 
sum-motif.Rmd 
