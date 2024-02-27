shell.prefix('set -x; set -e;')

sample_ids = open("data/{}/meta_data/sample_ids.txt".format(config["dataset"])).read().strip().split("\n")
indir = "data/{}/fastq".format(config["dataset"])
outdir = "output/{}".format(config["dataset"])
dataset = "{}".format(config["dataset"])

rule all:
    input:
        #multiqc_fastp = expand("{outdir}/01fastp/log/multiqc_report.html",outdir=[outdir]),
        bigwig = expand("{outdir}/05wig/{sample_id}.bigwig",outdir=[outdir],sample_id=sample_ids),
        summary = expand("{outdir}/06medips_qc/quality-control.txt",outdir=[outdir]),
        gene_matrix= expand("{outdir}/07counts/count_matrix_gene.txt",outdir=[outdir]),
        promoter_matrix= expand("{outdir}/07counts/count_matrix_promoter.txt",outdir=[outdir]),
        CGI_matrix= expand("{outdir}/07counts/count_matrix_CGI.txt",outdir=[outdir])

rule fastp:
    input:
        mate1=indir+'/{sample_id}.fastq.gz',
        #mate2=indir+'/{sample_id}_2.fastq.gz'
    output:
        fastq1='{outdir}/01fastp/{sample_id}.fastq.gz',
        #astq2='{outdir}/01fastp/{sample_id}_2.fastq.gz'
    params: 
        html='{outdir}/01fastp/log/{sample_id}_fastp.html',
        json='{outdir}/01fastp/log/{sample_id}_fastp.json'
    log:
        '{outdir}/01fastp/log/{sample_id}_fastp.log'
    threads:
        config["bowtie2_threads"]
    shell:
        """
        fastp --thread {threads} -i {input.mate1} \
        -o {output.fastq1} \
        -h {params.html} -j {params.json} \
        > {log} 2>&1
        """

# rule multiqc_fastp:
#     input:  
#         html = outdir + "/01fastp/log/{}_fastp.html".format(sample_id),
#         json = outdir + "/01fastp/log/{}_fastp.json".format(sample_id)
#     output:
#         multiqc_fastp = "{outdir}/01fastp/log/multiqc_report.html"
#     params:
#         "{outdir}/01fastp/log"
#     # wrapper:
#     #     "0.70.0/bio/multiqc" 
#     shell:
#         """
#         cd {params}
#         multiqc .
#         cd -
#         """

rule map_lambda:
    input:
        fastq1=outdir+'/01fastp/{sample_id}.fastq.gz',
        #fastq2=outdir+'/01fastp/{sample_id}_2.fastq.gz',
        index='ref/lambda_spikein/lambda_spikein.1.bt2'
    log:
        '{outdir}/02bam/log/{sample_id}_nolambda_bowtie2.log'
    params:
        unmapped_path='{outdir}/unmapped/{sample_id}_nolambda.fastq.gz',
        index='ref/lambda_spikein/lambda_spikein'
    output:
        bam='{outdir}/unmapped/{sample_id}_lambda.bam',
        unmapped1=temp('{outdir}/unmapped/{sample_id}_nolambda.fastq.gz'),
        #unmapped2=temp('{outdir}/unmapped/{sample_id}_nolambda_2.fastq.gz')
    threads:
        config["bowtie2_threads"]
    shell:
        """
        bowtie2 -p {threads} -U {input.fastq1} \
         --no-unal --un-gz {params.unmapped_path} -x {params.index} 2> {log} | \
         samtools view -b | samtools sort -m 256M -@ {threads} > {output.bam}
        """

rule map_hg38:
    input:
        fastq1=outdir+'/unmapped/{sample_id}_nolambda.fastq.gz',
        #fastq2=outdir+'/unmapped/{sample_id}_nolambda_2.fastq.gz',
        index='ref/hg38/genome.1.bt2'
    log:
        '{outdir}/02bam/log/{sample_id}_bowtie2.log'
    params:
        unmapped_path='{outdir}/unmapped/{sample_id}.fastq.gz',
        index='ref/hg38/genome'
    output:
        bam='{outdir}/02bam/{sample_id}.bam',
        unmapped1='{outdir}/unmapped/{sample_id}.fastq.gz',
        #unmapped2='{outdir}/unmapped/{sample_id}_2.fastq.gz'
    threads:
        config["bowtie2_threads"]
    shell:
        """
        #align&sort
        bowtie2 -p {threads} -U {input.fastq1} \
         --no-unal --un-gz {params.unmapped_path} -x {params.index} 2> {log} | \
         samtools view -b | samtools sort -m 256M -@ {threads} > {output.bam}
        """

#rule bam_filter:
    # input:
    #     bam="{outdir}/02bam/{sample_id}.bam"
    # output:
    #     bam=temp("{outdir}/03bam_filter/{sample_id}.bam")
    # params:
    #     min_size=config["insertion_min"],
    #     max_size=config["insertion_max"]
    # shell:
    #     """
    #     # PE with " -f 2 "
    #     (samtools view -H {input.bam}; samtools view -F 4 {input.bam} | \
    #     awk -v min_size={params.min_size} -v max_size={params.max_size} \
    #     '{{if($9>0){{size=$9}}else{{size=-$9}};if(size>=min_size&&size<=max_size){{print $0}}}}') | \
    #     samtools view -b > {output.bam}
    #     """

rule bam_dedup:
    input:
        bam = "{outdir}/02bam/{sample_id}.bam"
    output:
        bam = "{outdir}/04bam_dedup/{sample_id}.bam"
    params:
        bai = "{outdir}/04bam_dedup/{sample_id}.bam.bai",
        metrics = "{outdir}/04bam_dedup/log/{sample_id}_dedup.txt"
    log:
        '{outdir}/04bam_dedup/log/{sample_id}_dedup.log'
    threads:
        config["bowtie2_threads"]
    shell:
        """
        java -jar /data/baopengfei/biosoft/picard/picard.jar \
            MarkDuplicates REMOVE_DUPLICATES=true \
            ASSUME_SORT_ORDER=coordinate \
            I={input.bam} \
            O={output.bam} \
            M={params.metrics} \
            > {log} 2>&1
            
        #samtools index -@ {threads} {output.bam}
        """

rule wig:
    input:
        bam = "{outdir}/02bam/{sample_id}.bam" if not config["dedup"] else "{outdir}/04bam_dedup/{sample_id}.bam",
        #bai = "{outdir}/02bam/{sample_id}.bam.bai" if not config["dedup"] else "{outdir}/04bam_dedup/{sample_id}.bam.bai"
    output:
        bigwig = "{outdir}/05wig/{sample_id}.bigwig"
    threads:
        config["bowtie2_threads"]
    params:
        binsize = config["binsize"]
    shell:
        """
        bamCoverage --numberOfProcessors {threads} --extendReads 200 --normalizeUsing CPM -b {input.bam} -o {output.bigwig}
        """
        
rule medips_qc:
    input:
        bam = "{outdir}/02bam/{sample_id}.bam" if not config["dedup"] else "{outdir}/04bam_dedup/{sample_id}.bam"
    output:
        saturation_qc = "{outdir}/06medips_qc/saturation/{sample_id}.txt",
        enrich_qc = "{outdir}/06medips_qc/enrich/{sample_id}.txt",
        coverage = "{outdir}/06medips_qc/coverage/{sample_id}.png"
    conda:
        "envs/r35.yml"
    log:
        '{outdir}/06medips_qc/log/{sample_id}_medipsqc.log'
    params:
        #qcdir = "{outdir}/06medips_qc"
        binsize = config["binsize"]
    shell:
        """
        echo $PATH
        which -a R
        Rscript scripts/medips_qc_se.R -i {input.bam} -w {params.binsize} -sr {output.saturation_qc} -er {output.enrich_qc} -cr {output.coverage} \
	    > {log} 2>&1
        """

rule summary_qc:
    input:
        enrichments = expand("{outdir}/06medips_qc/enrich/{sample_id}.txt",outdir=[outdir],sample_id=sample_ids),
        saturations = expand("{outdir}/06medips_qc/saturation/{sample_id}.txt",outdir=[outdir],sample_id=sample_ids),
        coverage = expand("{outdir}/06medips_qc/coverage/{sample_id}.png",outdir=[outdir],sample_id=sample_ids)
    output:
        summary = "{outdir}/06medips_qc/quality-control.txt"
    run:
        import pandas as pd
        sample_ids = [path.split("/")[-1].split(".")[0] for path in input.enrichments]
        records = []
        for sample_id in sample_ids:
            enrichment_path = wildcards.outdir + "/06medips_qc/enrich/{}.txt".format(sample_id)
            with open(enrichment_path) as f:
                for line in f:
                    key,value = line.strip().split("\t")
                    if key == "enrichment.score.relH":
                        relH = value
                    elif key == "enrichment.score.GoGe":
                        GoGe = value
            saturation_path = wildcards.outdir + "/06medips_qc/saturation/{}.txt".format(sample_id)
            sat_df = pd.read_csv(saturation_path,sep="\t")
            es_sat_df = sat_df[sat_df["data"]=="estimated"]
            estimated_saturation = es_sat_df.sort_values(by="subset")["correlation"].iloc[-1]
            ob_sat_df = sat_df[sat_df["data"]=="observed"]
            observed_saturation = ob_sat_df.sort_values(by="subset")["correlation"].iloc[-1]
            records.append((sample_id,relH,GoGe,estimated_saturation,observed_saturation))
        table = pd.DataFrame.from_records(records)
        table.columns = ["sample_id","enrichment.score.relH","enrichment.score.GoGe","estimated.max.saturation","observed.max.saturation"]
        table.to_csv(output.summary,sep="\t",index=False)

rule count_gene:
    input:
        bam = "{outdir}/02bam/{sample_id}.bam" if not config["dedup"] else "{outdir}/04bam_dedup/{sample_id}.bam"
    output:
        gene1 = "{outdir}/07counts/{sample_id}_gene",
        gene2 = "{outdir}/07counts/{sample_id}_gene.summary"
    log:
        '{outdir}/07counts/log/{sample_id}_gene.log'
    params:
        gtf='ref/gtf/gencodev27_onlygene.gtf'
    threads:
        config["bowtie2_threads"]
    shell:
        """
        featureCounts -T {threads} -O -t gene -g gene_id -M -p  \
        -a {params.gtf} \
        -o {output.gene1} {input.bam} \
        > {log} 2>&1
        """

rule count_promoter:
    input:
        bam = "{outdir}/02bam/{sample_id}.bam" if not config["dedup"] else "{outdir}/04bam_dedup/{sample_id}.bam"
    output:
        promoter1 = "{outdir}/07counts/{sample_id}_promoter",
        promoter2 = "{outdir}/07counts/{sample_id}_promoter.summary"
    log:
        '{outdir}/07counts/log/{sample_id}_promoter.log'
    params:
        gtf='ref/gtf/gencodev27_onlypromoter.gtf'
    threads:
        config["bowtie2_threads"]
    shell:
        """
        featureCounts -T {threads} -O -t promoter -g gene_id -M -p  \
        -a {params.gtf} \
        -o {output.promoter1} {input.bam} \
        > {log} 2>&1
        """

rule count_CGI:
    input:
        bam = "{outdir}/02bam/{sample_id}.bam" if not config["dedup"] else "{outdir}/04bam_dedup/{sample_id}.bam"
    output:
        CGI1 = "{outdir}/07counts/{sample_id}_CGI",
        CGI2 = "{outdir}/07counts/{sample_id}_CGI.summary"
    log:
        '{outdir}/07counts/log/{sample_id}_CGI.log'
    params:
        gtf='ref/gtf/CGI.gtf'
    threads:
        config["bowtie2_threads"]
    shell:
        """
        featureCounts -T {threads} -O -t CGI -g gene_id -M -p  \
        -a {params.gtf} \
        -o {output.CGI1} {input.bam} \
        > {log} 2>&1
        """

rule combine_matrix:
    input:
        CGI_count = expand("{outdir}/07counts/{sample_id}_CGI",outdir=[outdir],sample_id=sample_ids),
        promoter_count = expand("{outdir}/07counts/{sample_id}_promoter",outdir=[outdir],sample_id=sample_ids),
        gene_count = expand("{outdir}/07counts/{sample_id}_gene",outdir=[outdir],sample_id=sample_ids)
    output:
        #gene = "{outdir}/07counts/gene.txt",
        #gene_clean = expand("{outdir}/07counts/{sample_id}_gene_clean",outdir=[outdir],sample_id=sample_ids),
        gene_matrix = "{outdir}/07counts/count_matrix_gene.txt",
        promoter_matrix = "{outdir}/07counts/count_matrix_promoter.txt",
        CGI_matrix = "{outdir}/07counts/count_matrix_CGI.txt"
    log:
        '{outdir}/07counts/log/count_matrix.log'
        #gene_log = '{outdir}/07counts/log/count_matrix_gene.log',
        #promoter_log = '{outdir}/07counts/log/count_matrix_promoter.log',
        #CGI_log = '{outdir}/07counts/log/count_matrix_CGI.log'
    params:
        outdir = '{outdir}'
    conda:
        "envs/py27.yml"
    shell:
        """
        echo "outdir is {params.outdir}"
        python scripts/count_matrix.py {params.outdir} > {log} 2>&1
        """