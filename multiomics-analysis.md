---

# cfDNA Methylation/cfMeDIP (DIP-seq)

## pre-process

### required files

- smk: snakemake/DIP-seq-pe.snakemake
- cfg: config/test.yaml
- env:
  - snakemake/envs/DIP.yml
  - snakemake/envs/py27.yml
  - snakemake/envs/r35.yml

### test run

dst="test"
snakemake --rerun-incomplete --keep-going --printshellcmds --reason   --use-conda --nolock --latency-wait 200 --restart-times 1 --jobs 16   --snakefile snakemake/DIP-seq-pe.snakemake   --configfile config/${dst}.yaml \

> run-${dst}_methylation.log 2>&1

## methylation level

included within snakemake/DIP-seq-pe.snakemake (rule: count_matrix_gene)

---

# cfDNA WGS (DNA-seq)

## pre-process

### required files

- smk:
  - snakemake/DNA-seq-common.snakemake
  - snakemake/DNA-seq-pe.snakemake
- cfg: config/test.yaml
- env:
  - snakemake/envs/DNA.yml
    - DIP-seq/snakemake/envs/py27.yml
    - DIP-seq/snakemake/envs/r35.yml

### test run

dst="test"
snakemake --rerun-incomplete --keep-going --printshellcmds --reason   --use-conda --nolock --latency-wait 200 --restart-times 1 --jobs 16   --snakefile snakemake/DNA-seq-pe.snakemake   --configfile config/${dst}.yaml \

> run-${dst}_WGS.log 2>&1

## cfDNA WGS CNV

### required files

- smk: snakemake/DNA-seq-pe-CNV.snakemake
- cfg: config/test.yaml
- env:
  - snakemake/envs/DNA_CNV_WisecondorX.yml
  - snakemake/envs/DNA_CNV_CNVkit.yml

### test run

included within pre-process (snakemake/DNA-seq-pe-CNV.snakemake)

## cfDNA WGS FragSize

### required files

- smk: snakemake/DNA-seq-pe-FragSize.snakemake
- cfg: config/test.yaml
- env:
  - snakemake/envs/DNA_FragSize.yml

### test run

included within pre-process (snakemake/DNA-seq-pe-FragSize.snakemake)

## cfDNA WGS WPS

(modified from https://github.com/kircherlab/cfDNA)

### required files

- smk: snakemake/WPS/snakefile_WPS.smk
- cfg: snakemake/WPS/config/config.yml
- env: snakemake/WPS/workflow/envs/cfDNA.yml

### test run

```bash
snakemake --rerun-incomplete --keep-going --printshellcmds --reason \
  --use-conda --nolock --latency-wait 200 --restart-times 1 --jobs 16 \
  --snakefile snakemake/WPS/snakefile_WPS.smk \
  --configfile snakemake/WPS/config/config.yml \
  > run_WPS.log 2>&1
```
