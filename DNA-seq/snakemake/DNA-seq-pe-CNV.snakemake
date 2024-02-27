shell.prefix("set -x; set -e;")
include: 'DNA-seq-common.snakemake'


################################################################
# CNV (option1: mask repeats depth by jinyunfan)
################################################################
#Read counts of each bins across genome with repeats excluded
rule getSegmentCoverage:
    input:
        bam=bam_dir+"/{sample_id}.bam",
        gappedbins=gn_dir+"/hg38.bins."+cnv_binsize+".norepeats.bed",
        bins=gn_dir+"/hg38.bins."+cnv_binsize+".bed",
        chromsize=gn_dir+"/chrom.size"
    output:
        coverage=outdir+"/segment-coverage/{sample_id}.bed"
    threads: 5
    log: outdir+"/log/{sample_id}/getSegmentCoverage.log"
    conda:
        "envs/CNV_WisecondorX.yml"
    shell:
        """
        bedtools coverage -g {input.chromsize} -sorted -a {input.gappedbins} -b {input.bam} -mean > {output.coverage}.gap ;\
        bedtools map -a {input.bins} -b {output.coverage}.gap -c 4 -o mean > {output.coverage} ;\
        rm {output.coverage}.gap > {log} 2>&1
        """


################################################################
# CNV (option2: get 2 matrices from wisecondorx)
################################################################
## get npz files from bam files (wisecondorx auto correct GC)
rule wisecondorxConvert:
    input:
        bam=bam_dir+"/{sample_id}.bam",
    output:
        npz=outdir+"/wisecondorx/convert/{sample_id}.npz"
    conda:
        "envs/CNV_WisecondorX.yml"
    params:
        bs=int(cnv_binsize)*0.05
    log: 
        outdir+"/log/{sample_id}/wisecondorx-convert.log"
    shell:
        """
        WisecondorX convert --normdup --binsize {params.bs} {input.bam} {output.npz} > {log} 2>&1
        """

## get reference npz fils from normal sample npz files (at least 10 samples)
rule wisecondorxGetRef:
    input:
        npz=expand(outdir+"/wisecondorx/convert/{sample_id}.npz",sample_id=sample_ids)
    output:
        ref=outdir+"/wisecondorx/ref/reference_nipt_"+cnv_binsize+".npz"
    conda:
        "envs/CNV_WisecondorX.yml"
    params:
        nc=outdir+"/wisecondorx/convert/{"+",".join(ref_sample_ids)+"}.npz",
        bs=int(cnv_binsize)
    threads: 
        16
    log: 
        outdir+"/log/wisecondorx-getRef.log"
    shell:
        """
        WisecondorX newref {params.nc} {output.ref} --binsize {params.bs} --cpus 4 --nipt > {log} 2>&1
        """

## predict CNV use disease group and reference npz files
rule wisecondorxPredCNV:
    input:
        npz=outdir+"/wisecondorx/convert/{sample_id}.npz",
        ref=outdir+"/wisecondorx/ref/reference_nipt_"+cnv_binsize+".npz",
        gn_blacklist_path=gn_dir+"/hg38.blacklist.bed"
    output:
        cnv=outdir+"/wisecondorx/CNV/"+cnv_binsize+"/{sample_id}_segments.bed"
    conda:
        "envs/CNV_WisecondorX.yml"
    params:
        cnv_path=outdir+"/wisecondorx/CNV/"+cnv_binsize+"/{sample_id}",
        # bl="ref/gtf/hg38_blacklist.bed"
    log: 
        outdir+"/log/{sample_id}/wisecondorx-predCNV.log"
    shell:
        """
        WisecondorX predict {input.npz} {input.ref} {params.cnv_path} --blacklist {input.gn_blacklist_path} --plot --bed > {log} 2>&1
        """

## map segment level CNV to gene level
### prepare single sample CNV 
rule wisecondorxGeneCNV:
    input:
        cnv=outdir+"/wisecondorx/CNV/"+cnv_binsize+"/{sample_id}_segments.bed"
    output:
        log2ratio=outdir+"/wisecondorx/CNV/"+cnv_binsize+"/{sample_id}_gene_log2ratio.bed",
        zscore=outdir+"/wisecondorx/CNV/"+cnv_binsize+"/{sample_id}_gene_zscore.bed",
    conda:
        "envs/CNV_WisecondorX.yml"
    params:
        gene=gtf_dir+"/gene.bed",
        #rm {params.gene}
    log:
        log1=outdir+"/wisecondorx/CNV/"+cnv_binsize+"/log/{sample_id}_gene_log2ratio.log",
        log2=outdir+"/wisecondorx/CNV/"+cnv_binsize+"/log/{sample_id}_gene_zscore.log"
    shell:
        """
        bash scripts/getWisecondorxGeneCNV.sh {input.cnv} {params.gene} {output.log2ratio} {log.log1} {wildcards.sample_id} log2ratio
        bash scripts/getWisecondorxGeneCNV.sh {input.cnv} {params.gene} {output.zscore} {log.log2} {wildcards.sample_id} zscore
        """

### join single CNV file to CNV matrix
rule wisecondorxGeneCNVmat:
    input:
        log2ratio=expand(outdir+"/wisecondorx/CNV/"+cnv_binsize+"/{sample_id}_gene_log2ratio.bed",sample_id=sample_ids),
        zscore=expand(outdir+"/wisecondorx/CNV/"+cnv_binsize+"/{sample_id}_gene_zscore.bed",sample_id=sample_ids)
    output:
        log2ratio=outdir+"/matrix/CNVlog2ratio_matrix_gene.txt",
        zscore=outdir+"/matrix/CNVzscore_matrix_gene.txt"
    conda:
        "envs/CNV_WisecondorX.yml"
    params:
        log2ratio=outdir+"/wisecondorx/CNV/"+cnv_binsize+"/*_gene_log2ratio.bed",
        zscore=outdir+"/wisecondorx/CNV/"+cnv_binsize+"/*_gene_zscore.bed",
        gene=outdir+"/wisecondorx/CNV/"+cnv_binsize+"/*_gene_*.bed",
        
        # log2ratio=outdir+"/wisecondorx/CNV/"+cnv_binsize+"/"+sample_ids2+"_gene_log2ratio.bed",
        # zscore=outdir+"/wisecondorx/CNV/"+cnv_binsize+"/"+sample_ids2+"_gene_zscore.bed",
        # gene1=outdir+"/wisecondorx/CNV/"+cnv_binsize+"/{"+",".join(sample_ids)+"}_gene_{zscore"+","+"log2ratio}.bed", # zscore error 
        # tmp1=outdir+"/wisecondorx/CNV/"+cnv_binsize+"/{"+",".join(sample_ids)+"}_gene_zscore.bed",
        # tmp2=outdir+"/wisecondorx/CNV/"+cnv_binsize+"/{"+",".join(sample_ids)+"}_gene_log2ratio.bed",

        tmp=outdir+"/wisecondorxGeneCNVmat_tmp.txt"
    log:
        log1=outdir+"/matrix/log/CNVlog2ratio_matrix_gene.txt",
        log2=outdir+"/matrix/log/CNVzscore_matrix_gene.txt"
    shell:
        """
        bash scripts/multijoin.sh {output.log2ratio} {params.log2ratio} > {log.log1} 2>&1
        bash scripts/multijoin.sh {output.zscore} {params.zscore} > {log.log2} 2>&1

        cat {output.log2ratio} | tr " " "\t" > {params.tmp} && mv {params.tmp} {output.log2ratio}
        cat {output.zscore} | tr " " "\t" > {params.tmp} && mv {params.tmp} {output.zscore}
        rm {params.gene}
        """



################################################################
# CNV (option3: get matrix from CNVkit)
################################################################
## in one command

# rule CNVkit_batch:
#     input:
#         bam=expand(bam_dir+"/{sample_id}.bam",sample_id=treat_sample_ids),
#     output:
#         ref=outdir+"/cnvkit/ref/my_reference.cnn",
#         cnr=expand(outdir+"/cnvkit/CNV/{sample_id}.cnr",sample_id=treat_sample_ids),
#         # cns=expand(outdir+"/cnvkit/CNV/{sample_id}.cns",sample_id=treat_sample_ids),
#     conda:
#         "envs/DNA.yml"
#     params:
#         bam=bam_dir+"/{"+",".join(treat_sample_ids)+"}.bam",
#         nc=bam_dir+"/{"+",".join(ref_sample_ids)+"}.bam",
#         out=outdir+"/cnvkit/CNV",
#         # bs=int(cnv_binsize)
#     threads: 
#         8
#     log: 
#         outdir+"/log/cnvkit-batch.log"
#     shell:
#         """
#         cnvkit.py batch \
#             {params.bam} \
#             --normal {params.nc} \
#             --method wgs \
#             -p 4 \
#             --fasta {gn_fa_path} \
#             --output-reference {output.ref} \
#             --output-dir {params.out} \
#             --scatter > {log} 2>&1
#         """

## CNVkit step-by-step: build ctrl samples reference
rule CNVkitGetRef:
    input:
        bam=expand(bam_dir+"/{sample_id}.bam",sample_id=ref_sample_ids),
        binGn_normsk_noblacklist=gn_dir+"/hg38.bins."+cnv_binsize+".norepeats.noblacklist.bed",
        gene=gtf_dir+"/gene.bed",
        gn_fa_path=gn_dir+"/genome.fa"
    output:
        ref=outdir+"/CNVkit/ref/"+cnv_binsize+"_my_reference.cnn",
    conda:
        "envs/CNV_CNVkit.yml"
    # params:
    #     gene=gtf_dir+"/gene.bed",
    threads: 
        12
    log: 
        outdir+"/log/CNVkit-getRef.log"
    shell:
        """
        cnvkit.py batch \
            --normal {input.bam} \
            --processes 6 \
            --fasta {input.gn_fa_path} \
            --targets {input.binGn_normsk_noblacklist}\
            --target-avg-size {cnv_binsize} \
            --annotate {input.gene} \
            --output-reference {output} > {log} 2>&1 
        """

## CNVkit step-by-step: predict CNV use disease (and ctrl) samples and reference
rule CNVkitPredCNV:
    input:
        bam=expand(bam_dir+"/{sample_id}.bam",sample_id=sample_ids), # treat_sample_ids
        ref=outdir+"/CNVkit/ref/"+cnv_binsize+"_my_reference.cnn",
    output:
        cnr=outdir+"/CNVkit/CNV/"+cnv_binsize+"/{sample_id}.cnr",
        cns=outdir+"/CNVkit/CNV/"+cnv_binsize+"/{sample_id}.cns"
    conda:
        "envs/CNV_CNVkit.yml"
    params:
        cnv_out_dir=outdir+"/CNVkit/CNV/"+cnv_binsize,
        # bl="ref/gtf/hg38_blacklist.bed"
    threads: 
        12
    log: 
        outdir+"/log/{sample_id}/CNVkit-predCNV.log"
    shell:
        """
        cnvkit.py batch \
            {input.bam} \
            --seq-method wgs \
            --reference {input.ref} \
            --processes 8 \
            --output-dir {params.cnv_out_dir} \
            --scatter --diagram
        """