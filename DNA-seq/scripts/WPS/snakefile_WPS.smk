#
# 174:         awk 'BEGIN{{ FS="\\t"; OFS="\\t" }}{{ print $1,$2-{params.flank},$3+{params.flank},$4,$5,$6 }}' |
# test add shell.prefix due to plot_overlay.py error
#shell.prefix("set -x; set -e;")

from snakemake.utils import validate
import pandas as pd
from os.path import exists


#configfile: "config/config.yml"


validate(config, schema="workflow/schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample_id"]
validate(samples, schema="workflow/schemas/samples_DNA.schema.yaml")
regions = pd.read_csv(config["regions"], sep="\t").set_index("target", drop=False)
regions.index.names = ["region_id"]
validate(regions, schema="workflow/schemas/regions.schema.yaml")


def get_WPS_ref(sample):
    ref_samples = samples["ref_samples"][sample].split(",")
    genomes = samples["genome_build"][ref_samples].values.tolist()
    return expand(
        "results/intermediate/{{ID}}/table/{GENOME}/target/{{target_region}}--{ref_SAMPLE}_WPS.csv.gz",
        zip,
        ref_SAMPLE=ref_samples,
        GENOME=genomes,
    )


def get_COV_ref(sample):
    ref_samples = samples["ref_samples"][sample].split(",")
    genomes = samples["genome_build"][ref_samples].values.tolist()
    return expand(
        "results/intermediate/{{ID}}/table/{GENOME}/target/{{target_region}}--{ref_SAMPLE}_COV.csv.gz",
        zip,
        ref_SAMPLE=ref_samples,
        GENOME=genomes,
    )


def get_STARTS_ref(sample):
    ref_samples = samples["ref_samples"][sample].split(",")
    genomes = samples["genome_build"][ref_samples].values.tolist()
    return expand(
        "results/intermediate/{{ID}}/table/{GENOME}/target/{{target_region}}--{ref_SAMPLE}_STARTS.csv.gz",
        zip,
        ref_SAMPLE=ref_samples,
        GENOME=genomes,
    )


def get_WPS_background_ref(sample):
    ref_samples = samples["ref_samples"][sample].split(",")
    genomes = samples["genome_build"][ref_samples].values.tolist()
    return expand(
        "results/intermediate/{{ID}}/table/{GENOME}/background/{{target_region}}--{ref_SAMPLE}_WPS.background.csv.gz",
        zip,
        ref_SAMPLE=ref_samples,
        GENOME=genomes,
    )


def get_COV_background_ref(sample):
    ref_samples = samples["ref_samples"][sample].split(",")
    genomes = samples["genome_build"][ref_samples].values.tolist()
    return expand(
        "results/intermediate/{{ID}}/table/{GENOME}/background/{{target_region}}--{ref_SAMPLE}_COV.background.csv.gz",
        zip,
        ref_SAMPLE=ref_samples,
        GENOME=genomes,
    )


def get_STARTS_background_ref(sample):
    ref_samples = samples["ref_samples"][sample].split(",")
    genomes = samples["genome_build"][ref_samples].values.tolist()
    return expand(
        "results/intermediate/{{ID}}/background/{{GENOME}}/table/{{target_region}}--{ref_SAMPLE}_STARTS.background.csv.gz",
        zip,
        ref_SAMPLE=ref_samples,
        GENOME=genomes,
    )

#only use first row lenth, treat all other regions has the same length
def get_length(input):
    if exists(input):
        df = pd.read_csv(input, sep="\t", header=None)
        length = df[2] - df[1]
        return length[0]
    else:
        length=2000
        #print(f"{input} does not exist: using length of {length}")
        return length
	# length: may used for enough calculating extended region ;only return 1st value ?


rule all:
    input:
        expand(
            expand(
                "results/intermediate/{ID}/regions/{GENOME}/target_region/{target_region}_blacklist-excluded.bed.gz",
                zip,
                ID=samples["ID"],
                GENOME=samples["genome_build"],
                allow_missing=True,
            ),
            target_region=regions["target"],
        ),
        expand(
            expand(
                "results/intermediate/{ID}/regions/{GENOME}/background/{target_region}_background_regions.bed.gz",
                zip,
                ID=samples["ID"],
                GENOME=samples["genome_build"],
                allow_missing=True,
            ),
            target_region=regions["target"],
        ),
        expand(expand(
            "results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_WPS.csv.gz",
            zip,
            SAMPLE=samples["sample"],
            ID=samples["ID"],
            GENOME=samples["genome_build"],
            allow_missing=True,
        ),target_region=regions["target"],
        ),
        expand(expand(
            "results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_COV.csv.gz",
            zip,
            SAMPLE=samples["sample"],
            ID=samples["ID"],
            GENOME=samples["genome_build"],
            allow_missing=True,
        ),target_region=regions["target"],
        ),
        expand(expand(
            "results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_STARTS.csv.gz",
            zip,
            SAMPLE=samples["sample"],
            ID=samples["ID"],
            GENOME=samples["genome_build"],
            allow_missing=True,
        ),target_region=regions["target"],
        ),

        ## add normalize
        expand(expand(
            "results/intermediate/{ID}/table/{GENOME}/normalize/{target_region}--{SAMPLE}_{type}_{normalize}.tsv.gz",
            zip,
            SAMPLE=samples["sample"],
            ID=samples["ID"],
            GENOME=samples["genome_build"],
            allow_missing=True,
        ),target_region=regions["target"],type=["WPS","COV","STARTS"],normalize=["normalized"] # "normalized",
        ),


        #expand(
        #    expand(
        #        "results/plots/overlays/{ID}/{target_region}--{SAMPLE}_overlays.pdf",
        #        zip,
        #        SAMPLE=samples["sample"],
        #        ID=samples["ID"],
        #        allow_missing=True,
        #    ),
        #    target_region=regions["target"],
        #),


rule add_flanks:
    input:
        Region_file=lambda wildcards: regions["path"][wildcards.target_region],
    output:
        "results/intermediate/{ID}/regions/{GENOME}/target_region/{target_region}_flanking.bed.gz", 
    params:
        flank=1000, # orignail fixed version
    shell:
        """
        (fname={input.Region_file};
        if [[ $fname == *.gz ]];
        then zcat $fname;
        else cat $fname;
        fi;) | 
        awk 'BEGIN{{ FS="\\t"; OFS="\\t" }} $2-{params.flank} >= 0 {{ print $1,$2-{params.flank},$3+{params.flank},$4,$5,$6 }} $2-{params.flank} < 0 {{ print $1,0,$3+{params.flank},$4,$5,$6 }}' |
        gzip -c > {output}
        """ #"""
        #zcat {input.Region_file} | awk 'BEGIN{{ FS="\\t"; OFS="\\t" }}{{ print $1,$2-{params.flank},$3+{params.flank},$4,$5,$6 }}' > {output}
        #"""

rule exclude_blacklist:
    input:
        Region_file="results/intermediate/{ID}/regions/{GENOME}/target_region/{target_region}_flanking.bed.gz",
        blacklist=lambda wildcards: config[wildcards.GENOME]["universal_blacklist"],
    output:
        "results/intermediate/{ID}/regions/{GENOME}/target_region/{target_region}_blacklist-excluded.bed.gz",
    conda:
        "workflow/envs/background.yml"
    shell:
        """
        bedtools intersect -v -a {input.Region_file} -b {input.blacklist} |gzip -c > {output}
        """


#shuffle 1000 bg regions
rule generate_random_background:
    input:
        region="results/intermediate/{ID}/regions/{GENOME}/target_region/{target_region}_blacklist-excluded.bed.gz",
        genome=lambda wildcards: config[wildcards.GENOME]["genome_autosomes"],
        gap=lambda wildcards: config[wildcards.GENOME]["universal_blacklist"],
    output:
        "results/intermediate/{ID}/regions/{GENOME}/background/{target_region}_background_regions.bed.gz",
    params:
        length=lambda wildcards, input: get_length(input.region),
    conda:
        "workflow/envs/background.yml"
    shell:
        """
        bedtools random -n 1000 -l {params.length} -g {input.genome} | \
        bedtools shuffle -i stdin -g {input.genome} -excl {input.gap} -noOverlapping | \
        gzip -c > {output}
        """

rule extract_counts:
    input:
        target="results/intermediate/{ID}/regions/{GENOME}/target_region/{target_region}_blacklist-excluded.bed.gz",
        BAMFILE=lambda wildcards: samples["path"][wildcards.SAMPLE],
    output:
        WPS="results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_WPS.csv.gz",
        COV="results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_COV.csv.gz",
        STARTS="results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_STARTS.csv.gz",
        length="results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_length.csv.gz",
    params:
        minRL=config["minRL"], # usually target len + protection len (1bp TSS + 120bp = 121bp)
        maxRL=config["maxRL"],
        bpProtection=config["bpProtection"],
        out_pre="results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_%s.csv.gz",
    conda:
        "workflow/envs/cfDNA.yml"
    shell:
        """
        workflow/scripts/WPS/extractFromBAM_RegionBed_WPS_Cov.py \
        --minInsert={params.minRL} \
        --maxInsert={params.maxRL} \
        --protection={params.bpProtection} \
        -i {input.target} \
        -o {params.out_pre} {input.BAMFILE}
        """


rule extract_counts_background:
    input:
        background="results/intermediate/{ID}/regions/{GENOME}/background/{target_region}_background_regions.bed.gz",
        BAMFILE=lambda wildcards: samples["path"][wildcards.SAMPLE],
    output:
        WPS="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_WPS.background.csv.gz",
        COV="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_COV.background.csv.gz",
        STARTS="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_STARTS.background.csv.gz",
        length="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_length.background.csv.gz",
    params:
        minRL=config["minRL"],
        maxRL=config["maxRL"],
        out_pre="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_%s.background.csv.gz",
    conda:
        "workflow/envs/cfDNA.yml"
    shell:
        """
        workflow/scripts/WPS/extractFromBAM_RegionBed_WPS_Cov.py \
        --minInsert={params.minRL} \
        --maxInsert={params.maxRL} \
        -i {input.background} \
        -o {params.out_pre} {input.BAMFILE}
        """

"""
rule plot_overlays:
    input:
        WPS=lambda wc: "results/intermediate/{{ID}}/table/{GENOME}/target/{{target_region}}--{{SAMPLE}}_WPS.csv.gz".format(GENOME=samples["genome_build"].loc[samples["sample"] == wc.SAMPLE].values[0]),
        WPS_ref=lambda wildcards: get_WPS_ref(wildcards.SAMPLE),
        COV=lambda wc: "results/intermediate/{{ID}}/table/{GENOME}/target/{{target_region}}--{{SAMPLE}}_COV.csv.gz".format(GENOME=samples["genome_build"].loc[samples["sample"] == wc.SAMPLE].values[0]),
        COV_ref=lambda wildcards: get_COV_ref(wildcards.SAMPLE),
        #STARTS=lambda wc: "results/intermediate/{{ID}}/table/{GENOME}/target/{{target_region}}--{{SAMPLE}}_STARTS.csv.gz".format(GENOME=samples["genome_build"].loc[samples["sample"] == wc.SAMPLE].values[0]),
        #STARTS_ref=lambda wildcards: get_STARTS_ref(wildcards.SAMPLE),
        
        WPS_back=lambda wc: "results/intermediate/{{ID}}/table/{GENOME}/background/{{target_region}}--{{SAMPLE}}_WPS.background.csv.gz".format(GENOME=samples["genome_build"].loc[samples["sample"] == wc.SAMPLE].values[0]),
        WPS_back_ref=lambda wildcards: get_WPS_background_ref(wildcards.SAMPLE),
        COV_back=lambda wc: "results/intermediate/{{ID}}/table/{GENOME}/background/{{target_region}}--{{SAMPLE}}_COV.background.csv.gz".format(GENOME=samples["genome_build"].loc[samples["sample"] == wc.SAMPLE].values[0]),
        COV_back_ref=lambda wildcards: get_COV_background_ref(wildcards.SAMPLE),
        #STARTS_back=lambda wc: "results/intermediate/{{ID}}/table/{GENOME}/background/{{target_region}}--{{SAMPLE}}_STARTS.background.csv.gz".format(GENOME=samples["genome_build"].loc[samples["sample"] == wc.SAMPLE].values[0]),
        #STARTS_back_ref=lambda wildcards: get_STARTS_background_ref(wildcards.SAMPLE),

    output:
        "results/plots/overlays/{ID}/{target_region}--{SAMPLE}_overlays.pdf",
    params:
        target="{target_region}",
        sample="{SAMPLE}",
        ref_IDs=lambda wildcards: samples["ref_samples"][wildcards.SAMPLE].split(","),
        overlay_mode = config["plotting"]["overlay_mode"],
        smoothing = config["plotting"]["smoothing"],
        rolling = config["plotting"]["rolling"],
        background_norm = config["plotting"]["background_norm"]
    conda:
        "workflow/envs/overlays.yml"
    script:
        "workflow/scripts/WPS/overlays.py"
"""


## add normalize by shuffle 1000 bg regions
#results/intermediate/lulab/table/GRCh38/background/Exon1end--STAD-PKU-9-wgs_COV.background.csv.gz
rule normalize_WPS_by_shuffle1000regions:
    input:  
        target_WPS="results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_WPS.csv.gz",
        background_WPS="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_WPS.background.csv.gz",
        target_COV="results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_COV.csv.gz",
        background_COV="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_COV.background.csv.gz",
        target_STARTS="results/intermediate/{ID}/table/{GENOME}/target/{target_region}--{SAMPLE}_STARTS.csv.gz",
        background_STARTS="results/intermediate/{ID}/table/{GENOME}/background/{target_region}--{SAMPLE}_STARTS.background.csv.gz",
    output:
        output_WPS="results/intermediate/{ID}/table/{GENOME}/normalize/{target_region}--{SAMPLE}_WPS_normalized.tsv.gz",
        output_COV="results/intermediate/{ID}/table/{GENOME}/normalize/{target_region}--{SAMPLE}_COV_normalized.tsv.gz",
        output_STARTS="results/intermediate/{ID}/table/{GENOME}/normalize/{target_region}--{SAMPLE}_STARTS_normalized.tsv.gz",
    conda:
        "workflow/envs/overlays.yml"
    script:
        """workflow/scripts/expression_analysis/normalize.py"""

