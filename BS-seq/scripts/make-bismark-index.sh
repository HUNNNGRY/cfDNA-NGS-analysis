#!/bin/bash
bsub -q Z-LU -n 8 "bismark_genome_preparation --path_to_aligner /BioII/lulab_b/jinyunfan/anaconda3/envs/bioinfo_py36/bin --bowtie2 --parallel 8 --genomic_composition --verbose genome/index/bismark > bismark-index.log 2>&1"
