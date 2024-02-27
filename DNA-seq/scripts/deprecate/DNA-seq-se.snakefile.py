shell.prefix('set -x; set -e;')
## not adapted for dataset not using dedup bam , like pe
#config = {"dataset":"lulab","BQSR":False}

sample_ids = open("metadata/{}/sample_ids.txt".format(config["dataset"])).read().strip().split("\n")
indir = "data/{}".format(config["dataset"])
outdir = "output/{}".format(config["dataset"])
dataset = "{}".format(config["dataset"])

rule all:
    input:
       qc0 = expand("{outdir}/qc0/{sample_id}/{sample_id}_fastqc.html",outdir=[outdir],sample_id=sample_ids),
       qc1 = expand("{outdir}/qc1/{sample_id}/{sample_id}_fastqc.html",outdir=[outdir],sample_id=sample_ids),
       vcf = expand("{outdir}/vcf-filtered/{sample_id}.vcf.gz",outdir=[outdir],sample_id=sample_ids),
       #bwa = expand("{outdir}/bwa_unmapped/{sample_id}.bam",outdir=[outdir],sample_id=sample_ids),
       stat = expand("{outdir}/stats/{sample_id}.txt",outdir=[outdir],sample_id=sample_ids),
       flagstat  =expand("{outdir}/flagstat/{sample_id}.txt",outdir=[outdir],sample_id=sample_ids),
       microbe = expand("{outdir}/microbe/report/{sample_id}.txt",outdir=outdir,sample_id=sample_ids),
       counts = expand("{outdir}/bin-counts/{sample_id}.txt",outdir=[outdir],sample_id=sample_ids)
       #CNV = expand("{outdir}/CNV/CNV.txt",outdir=outdir)
       


rule qc0:
    input:
        fastq1=indir+'/{sample_id}.fastq.gz'
        #fastq2=indir+'/{sample_id}_2.fastq.gz'
    output:
        report1='{outdir}/qc0/{sample_id}/{sample_id}_fastqc.html'
        #report2='{outdir}/qc0/{sample_id}/{sample_id}_2_fastqc.html'
    params:
        outdir='{outdir}/qc0/{sample_id}'
    shell:
        """
        fastqc -o {params.outdir} {input.fastq1} 
        """

rule trimming:
    input:
        fastq1=indir+'/{sample_id}.fastq.gz'
        #fastq2=indir+'/{sample_id}_2.fastq.gz'
    output:
        fastq1 = '{outdir}/cutadapt/{sample_id}.fastq.gz',
        #fastq2 = '{outdir}/cutadapt/{sample_id}_2.fastq.gz',
        report1 = '{outdir}/log/{sample_id}/trimming_statistics.txt'
        #report2 = '{outdir}/log/{sample_id}/trimming_statistics_2.txt',
    params:
        outdir='{outdir}/cutadapt',
        quality = 30 
    threads:
        4
    log:
        log = '{outdir}/log/{sample_id}/trimming.txt'
    shell:
        """
        trim_galore --phred33 --cores {threads} --quality {params.quality} -o {params.outdir} --basename {wildcards.sample_id} {input.fastq1} > {log.log} 2>&1
        mv {params.outdir}/{wildcards.sample_id}_trimmed.fq.gz {output.fastq1} 
        mv {params.outdir}/{wildcards.sample_id}.fastq.gz_trimming_report.txt {output.report1}
        """

rule qc1:
    input:
        fastq1='{outdir}/cutadapt/{sample_id}.fastq.gz'
        #fastq2='{outdir}/cutadapt/{sample_id}_2.fastq.gz'
    output:
        report1='{outdir}/qc1/{sample_id}/{sample_id}_fastqc.html'
        #report2='{outdir}/qc1/{sample_id}/{sample_id}_2_fastqc.html'
    params:
        outdir='{outdir}/qc1/{sample_id}'
    shell:
        """
        fastqc -o {params.outdir} {input.fastq1}
        """

rule bwa_alignment:
    input:
        fastq1 = '{outdir}/cutadapt/{sample_id}.fastq.gz',
        #fastq2 = '{outdir}/cutadapt/{sample_id}_2.fastq.gz',
        index = 'genome/bwa-mem2-index2/genome.0123' 
    output:
        bam = "{outdir}/bwa_bam/{sample_id}.bam",
    params:
        index = "genome/bwa-mem2-index2/genome"
    threads:
        12
    log:
        '{outdir}/log/{sample_id}/bwa-alignment.txt'
    shell:
        '''
        /BioII/lulab_b/baopengfei/biosoft/bwa-mem2/bwa-mem2 mem -T 0 -t {threads} {params.index} {input.fastq1} | samtools view -b > {output.bam} 2> {log}
        '''

rule getUnaligned:
    input:
        bam = "{outdir}/bwa_bam/{sample_id}.bam"
    output:
        unmapped_1 = "{outdir}/bwa_unmapped/{sample_id}.fastq.gz"
        #unmapped_2 = "{outdir}/bwa_unmapped/{sample_id}_2.fastq.gz",
    log:
        '{outdir}/log/{sample_id}/get-unaligned.txt'
    shell:
        '''
        samtools fastq -s {output.unmapped_1} -f 13 {input.bam} 2> {log}        #not fully tested
        '''

rule bowtie2_alignment:
    input:
        fastq1 = '{outdir}/cutadapt/{sample_id}.fastq.gz',
        #fastq2 = '{outdir}/cutadapt/{sample_id}_2.fastq.gz',
        index='genome/bowtie2-index/genome.1.bt2'
    log:
        '{outdir}/log/{sample_id}/bowtie2.log'
    params:
        index='genome/bowtie2-index/genome'
        #unmapped_path='{outdir}/unmapped/{sample_id}.fastq.gz'
    output:
        bam="{outdir}/bam/{sample_id}.bam",
        unmapped1='{outdir}/unmapped/{sample_id}.fastq.gz'
        #unmapped2='{outdir}/unmapped/{sample_id}_2.fastq.gz'
    threads:
        12 
    shell:
        """
        bowtie2 -p {threads} -U {input.fastq1} \
         --no-unal --un-gz {output.unmapped1} -x {params.index} 2> {log} | samtools view -b > {output.bam}
        """



rule flagstat:
    input:
        bam = "{outdir}/bam/{sample_id}.bam"
    output:
        flagstat = "{outdir}/flagstat/{sample_id}.txt"
    shell:
        """
        samtools flagstat {input.bam} > {output.flagstat}
        """ 


rule sort:
    input:
        bam = "{outdir}/bam/{sample_id}.bam"
    output:
        bam = "{outdir}/bam-sorted/{sample_id}.bam",
        bai = "{outdir}/bam-sorted/{sample_id}.bam.bai"
    threads: 2
    shell:
        """
        samtools view -h -F 0x4  {input.bam} | samtools sort  -@ {threads}  > {output.bam}
        samtools index {output.bam}
        """

rule dedup:
    input:
        bam = "{outdir}/bam-sorted/{sample_id}.bam"
    output:
        bam = "{outdir}/bam-sorted-deduped/{sample_id}.bam",
        bai = "{outdir}/bam-sorted-deduped/{sample_id}.bam.bai",
        metrics = "{outdir}/log/{sample_id}/dedup-metrics.txt"
    log:
        "{outdir}/log/{sample_id}/MarkDuplicates.log"
    shell:
        """
         gatk MarkDuplicates --REMOVE_DUPLICATES true \
            --ASSUME_SORT_ORDER coordinate \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics} > {log} 2>&1
        samtools index {output.bam}
        """


rule getBamStatistics:
    input:
        bam = "{outdir}/bam-sorted-deduped/{sample_id}.bam"
    output:
        summary = "{outdir}/stats/{sample_id}.txt"
    shell:
        """
        samtools stats {input.bam} > {output.summary}
        """


rule addReadsGroup:
    input:
        bam = "{outdir}/bam-sorted-deduped/{sample_id}.bam"
    output:
        bam = "{outdir}/bam-sorted-deduped-RG/{sample_id}.bam",
        bai = "{outdir}/bam-sorted-deduped-RG/{sample_id}.bam.bai"
    log:
        '{outdir}/log/{sample_id}/addReadsGroup.log'
    shell:
        """
        gatk AddOrReplaceReadGroups --java-options -Xmx4G \
        --INPUT {input.bam} --OUTPUT {output.bam} -SO coordinate \
        --RGLB library --RGPL illumina --RGPU HiSeq2000 --RGSM {wildcards.sample_id} > {log} 2>&1
        samtools index {output.bam}
        """

rule baseQualityRecalibration:
    input:
        bam='{outdir}/bam-sorted-deduped-RG/{sample_id}.bam',
        dbSNP='/Share2/home/lulab/jinyunfan/data/SNP/dbSNP/All_20180418.vcf.gz',
        reference="genome/fasta/hg38.fa"
    output: 
        bam='{outdir}/bam-sorted-deduped-RG-BQSR/{sample_id}.bam',
        bai = '{outdir}/bam-sorted-deduped-RG-BQSR/{sample_id}.bam.bai',
        grpFile='{outdir}/BQSR-grp/{sample_id}.grp'
    log:
        prepareGrp = '{outdir}/log/{sample_id}/BaseRecalibrator.log',
        BQSR = '{outdir}/log/{sample_id}/ApplyBQSR.log'
    params:
        tmp = "tmp"
    shell:
        """
        gatk BaseRecalibrator --java-options -Xmx4G --input {input.bam} \
        --output {output.grpFile} \
        --tmp-dir {params.tmp} \
        --known-sites {input.dbSNP} --reference {input.reference} > {log.prepareGrp} 2>&1

        gatk ApplyBQSR -R {input.reference} -I {input.bam}  \
        --bqsr-recal-file {output.grpFile} --tmp-dir {params.tmp} -O {output.bam} > {log.BQSR} 2>&1
        samtool index {output.bam}
        """

rule HaplotypeCaller:
    input:
        bam='{outdir}/bam-sorted-deduped-RG-BQSR/{sample_id}.bam' if config["BQSR"] else '{outdir}/bam-sorted-deduped-RG/{sample_id}.bam',
        reference="genome/fasta/hg38.fa"
    output:
        vcf="{outdir}/vcf/{sample_id}.vcf.gz"
    log:
        '{outdir}/log/{sample_id}/HaplotypeCaller.log'
    params:
        tmp = "tmp"
    shell:
        """
        gatk HaplotypeCaller --java-options -Xmx4G -R {input.reference} \
        -I {input.bam} -O {output.vcf} \
        --tmp-dir {params.tmp}  > {log} 2>&1
        """

rule filterVcf:
    input:
        vcf="{outdir}/vcf/{sample_id}.vcf.gz",
        reference="genome/fasta/hg38.fa",
    output:
        vcf="{outdir}/vcf-filtered/{sample_id}.vcf.gz",
    log:
        '{outdir}/log/{sample_id}/variant-filtering.log'
    shell:
        """
        gatk VariantFiltration --java-options -Xmx4G \
        -R {input.reference} -V {input.vcf} \
        -window 35 -cluster 3 \
        --filter-name FS20 -filter "FS > 20.0" \
        --filter-name QD2 -filter "QD < 2.0" \
        --filter-name DP10 -filter "DP < 10.0" \
        --filter-name QUAL20 -filter "QUAL < 20.0" -O {output.vcf} > {log} 2>&1
        """ 
    
rule countMicrobe:
    input:
        unmapped_1 = "{outdir}/unmapped/{sample_id}.fastq.gz",
        #unmapped_2 = "{outdir}/unmapped/{sample_id}_2.fastq.gz",
        database = "/Share2/home/lulab/jinyunfan/data/kraken2db/standard-db"
    output:
        report = "{outdir}/microbe/report/{sample_id}.txt",
        unclassified_1 = "{outdir}/microbe/unclassified/{sample_id}.fastq.gz",
        #unclassified_2 = "{outdir}/microbe/unclassified/{sample_id}_2.fastq.gz",
        assignment = "{outdir}/microbe/assignment/{sample_id}.txt.gz"
    params:
        unclassified="{outdir}/microbe/unclassified/{sample_id}.fastq",
    threads:
        20
    log:
        '{outdir}/log/{sample_id}/kraken2-classification.log'
    shell:
        """
        ## se not fully tested

        LANG=C perl -e exit
        kraken2 --db {input.database} --threads {threads} --unclassified-out {params.unclassified} --report {output.report}  --use-names  {input.unmapped_1} > {output.assignment} 2> {log}
        gzip  {output.assignment}
        gzip  {wildcards.outdir}/microbe/unclassified/{wildcards.sample_id}.fastq 
        """

rule collectReadCounts:
    input:
        bam='{outdir}/bam-sorted-deduped-RG-BQSR/{sample_id}.bam' if config["BQSR"] else '{outdir}/bam-sorted-deduped-RG/{sample_id}.bam',
        reference = "genome/fasta/hg38.fa",
        intervals = "genome/interval-list/hg38.interval_list" 
    output:
        counts = "{outdir}/bin-counts/{sample_id}.txt"
    log:
        '{outdir}/log/{sample_id}/collect-count-by-bins.log'
    shell:
        """
         gatk CollectReadCounts  -L {input.intervals}  -R {input.reference} -imr OVERLAPPING_ONLY  \
         -I {input.bam} --format TSV -O {output.counts} > {log} 2>&1
        """

rule callCNV:
    input:
        paths = expand("{outdir}/bin-counts/{sample_id}.txt",sample_id=sample_ids,outdir=[outdir]),
        intervals = "genome/interval-list/hg38.interval_list"
    output:
        CNV = "{outdir}/CNV/CNV.txt"
    run:
        import subprocess
        prefix = output.CNV
        cmd = ["gatk","DetermineGermlineContigPloidy","-L",input.intervals,"--interval-merging-rule","OVERLAPPING_ONLY","--output",".","--output-prefix",prefix,"--verbosity","DEBUG"]
        for path in input.paths:
            cmd += ["-I",path]
        print(" ".join(cmd))
        #subprocess.run(cmd)
