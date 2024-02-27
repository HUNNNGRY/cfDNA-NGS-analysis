shell.prefix('set -x; set -e;')

sample_ids = open("data/{}/meta_data/sample_ids.txt".format(config["dataset"])).read().strip().split("\n")
indir = "data/{}/fastq".format(config["dataset"])
outdir = "output/{}".format(config["dataset"])
dataset = "{}".format(config["dataset"])
bam_dir = outdir+"/bam" if not config["remove_duplications"] else outdir+"/bam-sorted-deduped"
RG_dir = outdir+"/bam-sorted-RG" if not config["remove_duplications"] else outdir+'/bam-sorted-deduped-RG'
BQSR_dir = outdir+"/bam-sorted-RG-BQSR" if not config["remove_duplications"] else outdir+'/bam-sorted-deduped-RG-BQSR'

def request(config,outdir,sample_ids):
    output = dict()
    output["qc0"] = expand('{outdir}/qc0/{sample_id}/{sample_id}_1_fastqc.html',outdir=outdir,sample_id=sample_ids)
    output["qc1"] = expand('{outdir}/qc1/{sample_id}/{sample_id}_1_fastqc.html',outdir=outdir,sample_id=sample_ids)
    if config["remove_duplications"]: 
        output["bam"] = expand("{outdir}/bam-sorted-deduped/{sample_id}.bam",outdir=outdir,sample_id=sample_ids)
    else:
        output["bam"] = expand("{outdir}/bam/{sample_id}.bam",outdir=outdir,sample_id=sample_ids)
    output["bigwig"] = expand("{outdir}/wig/{sample_id}.bigwig",outdir=[outdir],sample_id=sample_ids),
    if config["SNV"]: 
        output["SNV"] = expand("{outdir}/vcf-filtered/{sample_id}.vcf.gz",outdir=[outdir],sample_id=sample_ids),
        if config["BQSR"]:
            output["BQSR"] = expand(BQSR_dir+"/{sample_id}.bam",sample_id=sample_ids), 
    if config["fragsize"]: 
        output["fragsize"] = expand("{outdir}/fragment-length/histogram.txt",outdir=[outdir]),
    if config["CNV"]: 
        output["CNV"] = expand("{outdir}/segment-coverage/{sample_id}.bed",outdir=[outdir],sample_id=sample_ids),
    output["wisecondorx_cnv"] = expand(outdir+"/wisecondorx/CNV/100000/{sample_id}_segments.bed",outdir=[outdir],sample_id=sample_ids),
    output["wisecondorx_gene_log2R"] = outdir+"/matrix/CNVlog2ratio_matrix_gene.txt",
    output["wisecondorx_gene"] = outdir+"/matrix/CNVzscore_matrix_gene.txt", 
    output["gene_matrix"]= outdir+"/matrix/count_matrix_gene.txt",
    output["gene_matrix_log"]= outdir+"/matrix/count_matrix_gene.txt.summary",
    output["promoter150TSS50_matrix"]= outdir+"/matrix/count_matrix_promoter150TSS50.txt",
    output["promoter150TSS50_matrix_log"]= outdir+"/matrix/count_matrix_promoter150TSS50.txt.summary"
    output["promoter300exon1end100_matrix"]= outdir+"/matrix/count_matrix_promoter300exon1end100.txt",
    output["promoter300exon1end100_matrix_log"]= outdir+"/matrix/count_matrix_promoter300exon1end100.txt.summary"
    output["stat"] = expand("{outdir}/stats/{sample_id}.txt",outdir=[outdir],sample_id=sample_ids),
    output["flagstat"]  =expand("{outdir}/flagstat/{sample_id}.txt",outdir=[outdir],sample_id=sample_ids),
#    output["depth"] = expand("{outdir}/stats/coverage_depth.txt",outdir=[outdir]),
    output["wgs_qc"]=expand("{outdir}/wgs_qc/quality-control.txt",outdir=[outdir]),
    if config["microbe"]: 
        output["microbe"] = expand("{outdir}/microbe/report/{sample_id}.txt",outdir=outdir,sample_id=sample_ids),
    if config["structure_variations"]: 
        output["structure_variations"] = expand("{outdir}/struatural-variation/{sample_id}/results/variants/candidateSV.vcf.gz",outdir=outdir,sample_id=sample_ids)
    return list(output.values())

rule all:
    input:request(config,outdir,sample_ids)

rule qc0:
    input:
        fastq1=indir+'/{sample_id}_1.fastq.gz',
        fastq2=indir+'/{sample_id}_2.fastq.gz'
    output:
        report1='{outdir}/qc0/{sample_id}/{sample_id}_1_fastqc.html',
        report2='{outdir}/qc0/{sample_id}/{sample_id}_2_fastqc.html'
    params:
        outdir='{outdir}/qc0/{sample_id}'
    shell:
        """
        fastqc -o {params.outdir} {input.fastq1} 
        fastqc -o {params.outdir} {input.fastq2}
        """

rule trimming:
    input:
        fastq1=indir+'/{sample_id}_1.fastq.gz',
        fastq2=indir+'/{sample_id}_2.fastq.gz'
    output:
        fastq1 = '{outdir}/cutadapt/{sample_id}_1.fastq.gz',
        fastq2 = '{outdir}/cutadapt/{sample_id}_2.fastq.gz',
        report1 = '{outdir}/log/{sample_id}/trimming_statistics_1.txt',
        report2 = '{outdir}/log/{sample_id}/trimming_statistics_2.txt'
    params:
        outdir='{outdir}/cutadapt',
        quality = 30 
    threads:
        4
    log:
        log = '{outdir}/log/{sample_id}/trimming.txt'
    shell:
        """
        trim_galore --phred33 --paired  --cores {threads} --quality {params.quality} -o {params.outdir} --basename {wildcards.sample_id} {input.fastq1} {input.fastq2} > {log.log} 2>&1
        mv {params.outdir}/{wildcards.sample_id}_val_1.fq.gz {output.fastq1}
        mv {params.outdir}/{wildcards.sample_id}_val_2.fq.gz {output.fastq2}
        mv {params.outdir}/{wildcards.sample_id}_1.fastq.gz_trimming_report.txt {output.report1}
        mv {params.outdir}/{wildcards.sample_id}_2.fastq.gz_trimming_report.txt {output.report2}
        """


# Quality control after adapter trimming
rule qc1:
    input:
        fastq1='{outdir}/cutadapt/{sample_id}_1.fastq.gz',
        fastq2='{outdir}/cutadapt/{sample_id}_2.fastq.gz'
    output:
        report1='{outdir}/qc1/{sample_id}/{sample_id}_1_fastqc.html',
        report2='{outdir}/qc1/{sample_id}/{sample_id}_2_fastqc.html'
    params:
        outdir='{outdir}/qc1/{sample_id}'
    shell:
        """
        fastqc -o {params.outdir} {input.fastq1}
        fastqc -o {params.outdir} {input.fastq2}
        """


# Mapping with bwa
rule bwa_alignment:
    input:
        fastq1 = '{outdir}/cutadapt/{sample_id}_1.fastq.gz',
        fastq2 = '{outdir}/cutadapt/{sample_id}_2.fastq.gz',
        index = 'genome/bwa-mem2-index/genome.0123' 
    output:
        bam = "{outdir}/bam/{sample_id}.bam",
        sam = temp("{outdir}/bam/{sample_id}.sam")
    params:
        index = "genome/bwa-mem2-index/genome",
    threads:
        6
    log:
        '{outdir}/log/{sample_id}/bwa-alignment.txt'
    shell:
        '''
        /BioII/lulab_b/baopengfei/biosoft/bwa-mem2/bwa-mem2 mem -T 0 -t {threads} {params.index} {input.fastq1} {input.fastq2} -o {output.sam} > {log} 2>&1 
        samtools view -b {output.sam} > {output.bam}
        '''

# Get unmapped reads from bwa's bam file
rule getUnaligned:
    input:
        bam = "{outdir}/bam/{sample_id}.bam"
    output:
        unmapped_1 = "{outdir}/unmapped/{sample_id}_1.fastq.gz",
        unmapped_2 = "{outdir}/unmapped/{sample_id}_2.fastq.gz"
    threads: 4    
    log:
        '{outdir}/log/{sample_id}/get-unaligned.txt'
    shell:
        '''
        samtools fastq -@ {threads} -1 {output.unmapped_1} -2 {output.unmapped_2} -0 /dev/null -s /dev/null -f 13 {input.bam} 2> {log}
        '''

# Generate flag statistics in bam file
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
    params:
        keep_proper_pair="-f 2" if config["onlykeep_properpair"] else "",
    threads: 4
    shell:
        """
        samtools view -h {params.keep_proper_pair} -F 0x4  {input.bam} | samtools sort  -@ {threads}  > {output.bam}
        samtools index -@ {threads} {output.bam}
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
    threads: 20
    shell:
        """
         gatk MarkDuplicates --REMOVE_DUPLICATES true \
            --ASSUME_SORT_ORDER coordinate \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics} > {log} 2>&1
        samtools index -@ {threads} {output.bam}
        """

rule wig:
    input:
        bam = "{outdir}/bam/{sample_id}.bam" if not config["remove_duplications"] else "{outdir}/bam-sorted-deduped/{sample_id}.bam"
    output:
        bigwig = "{outdir}/wig/{sample_id}.bigwig"
    log:
	    log = "{outdir}/wig/log/{sample_id}.log"
    threads:
        10 
    params:
        bin = 50
    shell:
        """
        bamCoverage --binSize {params.bin} --numberOfProcessors {threads} --extendReads --normalizeUsing CPM -b {input.bam} -o {output.bigwig} > {log.log} 2>&1 
        """

rule wgs_qc:
    input:
        bam = "{outdir}/bam/{sample_id}.bam" if not config["remove_duplications"] else "{outdir}/bam-sorted-deduped/{sample_id}.bam"
    output:
        saturation_qc = "{outdir}/wgs_qc/saturation/{sample_id}.txt.300",
        saturations_30k = "{outdir}/wgs_qc/saturation/{sample_id}.txt.3000",
        saturation_3M = "{outdir}/wgs_qc/saturation/{sample_id}.txt.30000",
        enrich_qc = "{outdir}/wgs_qc/enrich/{sample_id}.txt",
    conda:
        "envs/r35.yml"
    log:
        '{outdir}/wgs_qc/log/{sample_id}_wgsqc.log'
    params:
        '{outdir}/wgs_qc/saturation/{sample_id}.txt'
    shell:
        """
        Rscript scripts/wgs_qc.R -i {input.bam} -sr {params} -er {output.enrich_qc} \
            > {log} 2>&1
        """

rule summary_wgs_qc:
    input:
        enrichments = expand("{outdir}/wgs_qc/enrich/{sample_id}.txt",outdir=[outdir],sample_id=sample_ids),
        saturations = expand("{outdir}/wgs_qc/saturation/{sample_id}.txt.300",outdir=[outdir],sample_id=sample_ids),
        saturations_30k = expand("{outdir}/wgs_qc/saturation/{sample_id}.txt.3000",outdir=[outdir],sample_id=sample_ids),
        saturation_3M = expand("{outdir}/wgs_qc/saturation/{sample_id}.txt.30000",outdir=[outdir],sample_id=sample_ids),
    output:
        summary = "{outdir}/wgs_qc/quality-control.txt"
    run:
        import pandas as pd
        sample_ids = [path.split("/")[-1].split(".")[0] for path in input.enrichments]
        records = []
        for sample_id in sample_ids:
            enrichment_path = wildcards.outdir + "/wgs_qc/enrich/{}.txt".format(sample_id)
            with open(enrichment_path) as f:
                for line in f:
                    key,value = line.strip().split("\t")
                    if key == "enrichment.score.relH":
                        relH = value
                    elif key == "enrichment.score.GoGe":
                        GoGe = value
            saturation_path = wildcards.outdir + "/wgs_qc/saturation/{}.txt.300".format(sample_id)
            sat_df = pd.read_csv(saturation_path,sep="\t")
            es_sat_df = sat_df[sat_df["data"]=="estimated"]
            estimated_saturation = es_sat_df.sort_values(by="subset")["correlation"].iloc[-1]
            ob_sat_df = sat_df[sat_df["data"]=="observed"]
            observed_saturation = ob_sat_df.sort_values(by="subset")["correlation"].iloc[-1]

            saturation_path1 = wildcards.outdir + "/wgs_qc/saturation/{}.txt.3000".format(sample_id)
            sat_df1 = pd.read_csv(saturation_path1,sep="\t")
            es_sat_df1 = sat_df1[sat_df1["data"]=="estimated"]
            estimated_saturation1 = es_sat_df1.sort_values(by="subset")["correlation"].iloc[-1]
            ob_sat_df1 = sat_df1[sat_df1["data"]=="observed"]
            observed_saturation1 = ob_sat_df1.sort_values(by="subset")["correlation"].iloc[-1]

            saturation_path2 = wildcards.outdir + "/wgs_qc/saturation/{}.txt.30000".format(sample_id)
            sat_df2 = pd.read_csv(saturation_path2,sep="\t")
            es_sat_df2 = sat_df2[sat_df2["data"]=="estimated"]
            estimated_saturation2 = es_sat_df2.sort_values(by="subset")["correlation"].iloc[-1]
            ob_sat_df2 = sat_df2[sat_df2["data"]=="observed"]
            observed_saturation2 = ob_sat_df2.sort_values(by="subset")["correlation"].iloc[-1]

            records.append((sample_id,relH,GoGe,estimated_saturation,observed_saturation,estimated_saturation1,observed_saturation1,estimated_saturation2,observed_saturation2))
        table = pd.DataFrame.from_records(records)
        table.columns = ["sample_id","enrichment.score.relH","enrichment.score.GoGe","estimated.max.saturation","observed.max.saturation","estimated.max.saturation.3k","observed.max.saturation.3k","estimated.max.saturation.30k","observed.max.saturation.30k"]
        table.to_csv(output.summary,sep="\t",index=False)

rule getBamStatistics:
    input:
        bam = bam_dir+"/{sample_id}.bam"
    output:
        summary = outdir+"/stats/{sample_id}.txt"
    shell:
        """
        samtools stats {input.bam} > {output.summary}
        """


rule count_matrix_gene:
    input:
        bam = expand(bam_dir+"/{sample_id}.bam",outdir=outdir,sample_id=sample_ids)
    output:
        gene_matrix = outdir+"/matrix/count_matrix_gene.txt",
        gene_TPM = outdir+"/matrix/TPM_matrix_gene.txt",        
        gene_sum = outdir+"/matrix/count_matrix_gene.txt.summary",
    log:
        log2=outdir+'/matrix/log/count_matrix_gene.log',
    threads:
        20
    params:
        bam = bam_dir+"/*.bam",
        tmp = outdir+"/matrix/tmp",
        gtf2='ref/gtf/gene.gtf',
    shell:
        """
        # gene
        featureCounts -T {threads} -O -t gene -g gene_id -M -p  \
            -a  {params.gtf2} \
            -o {output.gene_matrix} {params.bam} \
            > {log.log2} 2>&1
        Rscript scripts/multiFeatureCounts2countMat.R   \
            -i {output.gene_matrix}  \
            -o {params.tmp} ;
        mv {params.tmp} {output.gene_matrix};
        Rscript  scripts/TPM.R \
            -i {output.gene_matrix} \
            -o  {output.gene_TPM} 
        """


rule count_matrix_promoter150TSS50:
    input:
        bam = expand(bam_dir+"/{sample_id}.bam",outdir=outdir,sample_id=sample_ids)
    output:      
        promoter150TSS50_matrix = outdir+"/matrix/count_matrix_promoter150TSS50.txt",
        promoter150TSS50_TPM = outdir+"/matrix/TPM_matrix_promoter150TSS50.txt",
        promoter150TSS50_sum = outdir+"/matrix/count_matrix_promoter150TSS50.txt.summary",
    log:
        log4=outdir+'/matrix/log/count_matrix_promoter150TSS50.log',
    threads:
        20
    params:
        bam = bam_dir+"/*.bam",
        tmp = outdir+"/matrix/tmp",
        gtf4 = 'ref/gtf/promoter150TSS50.gtf',
    shell:
        """
        # promoter150TSS50
        featureCounts -T {threads} -O -t promoter150TSS50 -g gene_id -M -p  \
            -a  {params.gtf4} \
            -o {output.promoter150TSS50_matrix} {params.bam} \
            > {log.log4} 2>&1
        Rscript scripts/multiFeatureCounts2countMat.R   \
            -i {output.promoter150TSS50_matrix}  \
            -o {params.tmp} ;
        mv {params.tmp} {output.promoter150TSS50_matrix};
        Rscript  scripts/TPM.R \
            -i {output.promoter150TSS50_matrix} \
            -o  {output.promoter150TSS50_TPM} 
        """


rule count_matrix_promoter300exon1end100:
    input:
        bam = expand(bam_dir+"/{sample_id}.bam",outdir=outdir,sample_id=sample_ids)
    output:
        promoter300exon1end100_matrix = outdir+"/matrix/count_matrix_promoter300exon1end100.txt",
        promoter300exon1end100_TPM = outdir+"/matrix/TPM_matrix_promoter300exon1end100.txt",
        promoter300exon1end100_sum = outdir+"/matrix/count_matrix_promoter300exon1end100.txt.summary"
    log:
        log5=outdir+'/matrix/log/count_matrix_promoter300exon1end100.log',
    threads:
        20
    params:
        bam = bam_dir+"/*.bam",
        tmp = outdir+"/matrix/tmp",
        gtf5='ref/gtf/promoter300exon1end100.gtf',
    shell:
        """
        # promoter300exon1end100
        featureCounts -T {threads} -O -t promoter300exon1end100 -g gene_id -M -p  \
            -a  {params.gtf5} \
            -o {output.promoter300exon1end100_matrix} {params.bam} \
            > {log.log5} 2>&1
        Rscript scripts/multiFeatureCounts2countMat.R   \
            -i {output.promoter300exon1end100_matrix}  \
            -o {params.tmp} ;
        mv {params.tmp} {output.promoter300exon1end100_matrix};
        Rscript  scripts/TPM.R \
            -i {output.promoter300exon1end100_matrix} \
            -o  {output.promoter300exon1end100_TPM} 
        """



# Summarize histogram of fragment length from samtools stats output
rule getFragmentSize:
    input:
        stats = expand("{outdir}/stats/{sample_id}.txt",outdir=outdir,sample_id=sample_ids)
    output:
        hist = "{outdir}/fragment-length/histogram.txt"
    run:
        import re
        import pandas as pd
        records = []
        for path in input.stats:
            sample_id = re.sub(".txt$","",path.split("/")[-1])
            with open(path) as f:
                for line in f:
                    if not line.startswith("IS"):
                        continue
                    #IS, insert size, pairs total, inward oriented pairs, outward oriented pairs, other pairs
                    fields = line.strip().split("\t")
                    IS, N = fields[1],fields[3]
                    records.append((sample_id,int(IS),int(N)))
        fsHistTable = pd.DataFrame.from_records(records)
        fsHistTable.columns = ["sample_id","insertion-size","count"]
        fsHist = fsHistTable.pivot(index="insertion-size",columns="sample_id",values="count")
        fsHist.to_csv(output.hist,sep="\t")
                                                                    

# Read counts of each bins across genome with repeats excluded
rule getSegmentCoverage:
    input:
        bam=bam_dir+'/{sample_id}.bam',
        gappedbins = "genome/hg38.bins."+config["binsize"]+".norepeats.bed",
        bins = "genome/hg38.bins."+config["binsize"]+".bed",
        chromsize = "genome/chrom.size"
    output:
        coverage = "{outdir}/segment-coverage/{sample_id}.bed"
    threads: 20
    log: '{outdir}/log/{sample_id}/getSegmentCoverage.log'
    shell:
        """
        bedtools coverage -g {input.chromsize} -sorted -a {input.gappedbins} -b {input.bam} -counts > {output.coverage}.tmp ;\
        bedtools map -a {input.bins} -b {output.coverage}.tmp  -c 4 -o sum > {output.coverage} ;\
        rm {output.coverage}.tmp > {log} 2>&1
        """

# CNV of each bins across genome using wisecondorx
rule wisecondorxConvert:
    input:
        bam = bam_dir+'/{sample_id}.bam',
    output:
        npz = outdir+"/wisecondorx/convert/{sample_id}.npz"
    log: 
        outdir+'/log/{sample_id}/wisecondorx-convert.log'
    shell:
        """
        WisecondorX convert {input.bam} {output.npz} > {log} 2>&1
        """

rule wisecondorxGetRef:
    input:
        npz = expand(outdir+"/wisecondorx/convert/{sample_id}.npz",sample_id=sample_ids)
    output:
        ref = temp(outdir+"/wisecondorx/ref/reference_100000.npz")
    params:
        nc = outdir+"/wisecondorx/convert/NC*.npz",
        bs = 100000
    threads: 
        20
    log: 
        outdir+'/log/wisecondorx-getRef.log'
    shell:
        """
        WisecondorX newref {params.nc} {output.ref} --binsize {params.bs} --cpus {threads} > {log} 2>&1
        """

rule wisecondorxPredCNV:
    input:
        npz = outdir+"/wisecondorx/convert/{sample_id}.npz",
        ref = outdir+"/wisecondorx/ref/reference_100000.npz"
    output:
        cnv = outdir+"/wisecondorx/CNV/100000/{sample_id}_segments.bed"
    params:
        cnv_path = outdir+"/wisecondorx/CNV/100000/{sample_id}",
	bl = "ref/gtf/hg38_blacklist.bed"
    log: 
        outdir+'/log/{sample_id}/wisecondorx-predCNV.log'
    shell:
        """
        WisecondorX predict {input.npz} {input.ref} {params.cnv_path} --blacklist {params.bl} --plot --bed > {log} 2>&1
        """

rule wisecondorxGeneCNV:
    input:
        cnv = outdir+"/wisecondorx/CNV/100000/{sample_id}_segments.bed"
    output:
        log2ratio = outdir+"/wisecondorx/CNV/100000/{sample_id}_gene_log2ratio.bed",
        zscore = outdir+"/wisecondorx/CNV/100000/{sample_id}_gene_zscore.bed",
    params:
        gene = "ref/gtf/gene.bed",
        #rm {params.gene}
    log:
        log1 = outdir+"/wisecondorx/CNV/100000/log/{sample_id}_gene_log2ratio.log",
        log2 = outdir+"/wisecondorx/CNV/100000/log/{sample_id}_gene_zscore.log"
    shell:
        """
        bash scripts/getWisecondorxGeneCNV.sh {input.cnv} {params.gene} {output.log2ratio} {log.log1} {wildcards.sample_id}
        bash scripts/getWisecondorxGeneCNV.sh {input.cnv} {params.gene} {output.zscore} {log.log2} {wildcards.sample_id}
        """

rule wisecondorxGeneCNVmat:
    input:
        log2ratio = expand(outdir+"/wisecondorx/CNV/100000/{sample_id}_gene_log2ratio.bed",sample_id=sample_ids),
        zscore = expand(outdir+"/wisecondorx/CNV/100000/{sample_id}_gene_zscore.bed",sample_id=sample_ids)
    output:
        log2ratio = outdir+"/matrix/CNVlog2ratio_matrix_gene.txt",
        zscore = outdir+"/matrix/CNVzscore_matrix_gene.txt"
    params:
        log2ratio = outdir+"/wisecondorx/CNV/100000/*_gene_log2ratio.bed",
        zscore = outdir+"/wisecondorx/CNV/100000/*_gene_zscore.bed",
        gene = outdir+"/wisecondorx/CNV/100000/*_gene_*.bed"
    log:
        log1 = outdir+"/matrix/log/CNVlog2ratio_matrix_gene.txt",
        log2 = outdir+"/matrix/log/CNVzscore_matrix_gene.txt"
    shell:
        """
        bash scripts/multijoin.sh {output.log2ratio} {params.log2ratio} > {log.log1} 2>&1
        bash scripts/multijoin.sh {output.zscore} {params.zscore} > {log.log2} 2>&1
        rm {params.gene}
        """


# Add a dummy read group as gatk's caller requires this
rule addReadsGroup:
    input:
        bam=outdir+"/bam/{sample_id}.bam" if not config["remove_duplications"] else outdir+"/bam-sorted-deduped/{sample_id}.bam"
    output:
        bam=RG_dir+"/{sample_id}.bam",
        bai=RG_dir+"/{sample_id}.bam.bai"
    params:
        id="{sample_id}"
    threads: 10
    log:
        log=outdir+"/log/{sample_id}/addReadsGroup.log"
    shell:
        """
        gatk AddOrReplaceReadGroups --java-options -Xmx15G --TMP_DIR temp\
        --INPUT {input.bam} --OUTPUT {output.bam} -SO coordinate \
        --RGLB library --RGPL illumina --RGPU HiSeq2000 --RGSM {params.id} > {log.log} 2>&1
        samtools index -@ {threads} {output.bam}
        """

rule baseQualityRecalibration:
    input:
        bam=RG_dir+"/{sample_id}.bam",
        dbSNP='ref/SNV_ref/All_20180418.vcf.gz',
        reference="genome/fasta/hg38.fa"
    output: 
        bam=BQSR_dir+'/{sample_id}.bam',
        bai=BQSR_dir+'/{sample_id}.bam.bai',
        grpFile=BQSR_dir+'/BQSR-grp/{sample_id}.grp'
    log:
        prepareGrp=outdir+'/log/{sample_id}/BaseRecalibrator.log',
        BQSR=outdir+'/log/{sample_id}/ApplyBQSR.log'
    threads: 20
    params:
        tmp = "tmp"
    shell:
        """
	# --java-options -Xmx15G
        gatk BaseRecalibrator --input {input.bam} \
        --output {output.grpFile} \
        --tmp-dir {params.tmp} \
        --known-sites {input.dbSNP} --reference {input.reference} > {log.prepareGrp} 2>&1

        gatk ApplyBQSR -R {input.reference} -I {input.bam}  \
        --bqsr-recal-file {output.grpFile} --tmp-dir {params.tmp} -O {output.bam} > {log.BQSR} 2>&1
        samtool index -@ {threads} {output.bam}
        """
        
# Call SNV
rule HaplotypeCaller:
    input:
        bam=BQSR_dir+"/{sample_id}.bam" if config["BQSR"] else RG_dir+'/{sample_id}.bam',
        reference="genome/fasta/hg38.fa"
    output:
        vcf = "{outdir}/vcf/{sample_id}.vcf.gz"
    log:
        '{outdir}/log/{sample_id}/HaplotypeCaller.log'
    params:
        tmp = "tmp"
    threads: 40
    shell:
        """
	# omit the --java-options -Xmxn15G option from the Java command line then a default value will be used.
        gatk HaplotypeCaller   -R {input.reference} \
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
    threads: 20
    shell:
        """
        gatk VariantFiltration --java-options -Xmx15G \
        -R {input.reference} -V {input.vcf} \
        -window 35 -cluster 3 \
        --filter-name FS20 -filter "FS > 20.0" \
        --filter-name QD2 -filter "QD < 2.0" \
        --filter-name DP10 -filter "DP < 10.0" \
        --filter-name QUAL20 -filter "QUAL < 20.0" -O {output.vcf} > {log} 2>&1
        """ 


# Classifiy microbe with kraken2    
rule countMicrobe:
    input:
        unmapped_1 = "{outdir}/unmapped/{sample_id}_1.fastq.gz",
        unmapped_2 = "{outdir}/unmapped/{sample_id}_2.fastq.gz",
        database = "/Share2/home/lulab/jinyunfan/data/kraken2db/standard-db"
    output:
        report = "{outdir}/microbe/report/{sample_id}.txt",
        unclassified_1 = "{outdir}/microbe/unclassified/{sample_id}_1.fastq.gz",
        unclassified_2 = "{outdir}/microbe/unclassified/{sample_id}_2.fastq.gz",
        assignment = "{outdir}/microbe/assignment/{sample_id}.txt.gz"
    params:
        unclassified="{outdir}/microbe/unclassified/{sample_id}#.fastq",
    threads:
        40
    log:
        '{outdir}/log/{sample_id}/kraken2-classification.log'
    shell:
        """
        LANG=C perl -e exit
        kraken2 --db {input.database} --paired --threads 15 --unclassified-out {params.unclassified} --report {output.report}  --use-names  {input.unmapped_1} {input.unmapped_2}  >  {output.assignment} 2> {log}
        gzip  {output.assignment}
        gzip  {wildcards.outdir}/microbe/unclassified/{wildcards.sample_id}_1.fastq 
        gzip  {wildcards.outdir}/microbe/unclassified/{wildcards.sample_id}_2.fastq 
        """


# Config manta structural variation caller
rule prepareMantaConfig:
    input:
        bam=BQSR_dir+'/{sample_id}.bam' if config["BQSR"] else "{outdir}/bam-sorted-deduped/{sample_id}.bam", 
        reference = "genome/fasta/hg38.fa"
    output:
        config = "{outdir}/struatural-variation/{sample_id}/runWorkflow.py.config.pickle",
        script = "{outdir}/struatural-variation/{sample_id}/runWorkflow.py" 
    log: '{outdir}/log/{sample_id}/prepareMantaConfig.log'
    shell:
        """
        configManta.py --bam {input.bam} --referenceFasta {input.reference} --runDir "{wildcards.outdir}/struatural-variation/{wildcards.sample_id}"  > {log} 2>&1
        """

# Call structure variation
rule runMantan:
    input:
        bam=BQSR_dir+'/{sample_id}.bam' if config["BQSR"] else "{outdir}/bam-sorted-deduped/{sample_id}.bam", 
        reference = "genome/fasta/hg38.fa",
        config = "{outdir}/struatural-variation/{sample_id}/runWorkflow.py.config.pickle",
        script = "{outdir}/struatural-variation/{sample_id}/runWorkflow.py"
    output:
        sv = "{outdir}/struatural-variation/{sample_id}/results/variants/candidateSV.vcf.gz"
    threads: 40
    log: '{outdir}/log/{sample_id}/runMantan.log'
    shell:
        """
        {input.script} -j 15 > {log} 2>&1
        """


