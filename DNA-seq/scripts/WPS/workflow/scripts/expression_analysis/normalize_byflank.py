import pandas as pd
import numpy as np
from scipy import stats
import gzip


target_WPS = snakemake.input["target_WPS"]
target_WPS_v2 = snakemake.input["target_WPS_v2"]
target_COV = snakemake.input["target_COV"]
target_STARTS = snakemake.input["target_STARTS"]


background_WPS0 = snakemake.input["background_WPS0"]
background_WPS0_v2 = snakemake.input["background_WPS0_v2"]
background_COV0 = snakemake.input["background_COV0"]
background_STARTS0 = snakemake.input["background_STARTS0"]

background_WPS1 = snakemake.input["background_WPS1"]
background_WPS1_v2 = snakemake.input["background_WPS1_v2"]
background_COV1 = snakemake.input["background_COV1"]
background_STARTS1 = snakemake.input["background_STARTS1"]

# background_WPS2 = snakemake.input["background_WPS2"]
# background_WPS2_v2 = snakemake.input["background_WPS2_v2"]
# background_COV2 = snakemake.input["background_COV2"]
# background_STARTS2 = snakemake.input["background_STARTS2"]

# background_WPS3 = snakemake.input["background_WPS3"]
# background_WPS3_v2 = snakemake.input["background_WPS3_v2"]
# background_COV3 = snakemake.input["background_COV3"]
# background_STARTS3 = snakemake.input["background_STARTS3"]

# background_WPS4 = snakemake.input["background_WPS4"]
# background_WPS4_v2 = snakemake.input["background_WPS4_v2"]
# background_COV4 = snakemake.input["background_COV4"]
# background_STARTS4 = snakemake.input["background_STARTS4"]


output_WPS = snakemake.output["output_WPS"]
output_WPS_v2 = snakemake.output["output_WPS_v2"]
output_COV = snakemake.output["output_COV"]
output_STARTS = snakemake.output["output_STARTS"]

index_choice=None

# add if WPS then delta not FC 
def normalize_sample(path_a: str, path_b0: str, path_b1: str, method: str="FoldChange"): # path_b2: str, path_b3: str, path_b4: str,
    """Reads .csv file, calculates mean over all rows and divides by trimmed mean.

    Args:
        path (str): Path to a .csv file

    Returns:
        [type]: [description]
    """
    #with open(path_a, 'rb') as fd:
    #    gzip_fd = gzip.GzipFile(fileobj=fd)
    #path_a = "/BioII/lulab_b/baopengfei/gitsoft/cfDNA-2016/results/intermediate/lulab_WSQ-JYF/table/GRCh38/target/PUM2_mRNA_UTR3--NC_ChQ-10_WPS.csv.gz"
    #path_b = "/BioII/lulab_b/baopengfei/gitsoft/cfDNA-2016/results/intermediate/lulab_WSQ-JYF/table/GRCh38/background/PUM2_mRNA_UTR3--NC_ChQ-10_WPS.backgroundFlank.csv.gz"
    sample_a = pd.read_csv(path_a,compression='gzip', encoding='ISO-8859–1', header=None) 
    if pd.api.types.is_string_dtype(sample_a[0]):
        sample_a = sample_a.set_index(0)
        sample_a.index.name="ID"
    sample_b0 = pd.read_csv(path_b0,compression='gzip', encoding='ISO-8859–1', header=None).mean(axis=1) # mean of each row single base
    sample_b1 = pd.read_csv(path_b1,compression='gzip', encoding='ISO-8859–1', header=None).mean(axis=1) # mean of each row single base
    # sample_b2 = pd.read_csv(path_b2,compression='gzip', encoding='ISO-8859–1', header=None).mean(axis=1) # mean of each row single base
    # sample_b3 = pd.read_csv(path_b3,compression='gzip', encoding='ISO-8859–1', header=None).mean(axis=1) # mean of each row single base
    # sample_b4 = pd.read_csv(path_b4,compression='gzip', encoding='ISO-8859–1', header=None).mean(axis=1) # mean of each row single base

    sample = sample_a
    if method == "FoldChange":
        #loop each column/base/position
        for i in range(0,sample_a.shape[1]-1):
            #sys.stderr.write("using FoldChange normalize...\n")
            sample.iloc[:,i] = sample_a.iloc[:,i].values / (0.5*(sample_b0.values+sample_b1.values)+0.1) # pseudo count: 0.1
    else:
        #print("WPS using Delta")  
        for i in range(0,sample_a.shape[1]-1):
            #sys.stderr.write("using Delta normalize...\n")
            sample.iloc[:,i] = sample_a.iloc[:,i].values - 0.5*(sample_b0.values+sample_b1.values)
    return sample.round(4)
#        sample = sample_a / stats.trim_mean(sample_b, 0.1) # + 0.1 # pseudo count

normalized_WPS = normalize_sample(path_a=target_WPS, path_b0=background_WPS0, path_b1=background_WPS1, method="Delta")
normalized_WPS_v2 = normalize_sample(path_a=target_WPS_v2, path_b0=background_WPS0_v2, path_b1=background_WPS1_v2, method="FoldChange")
normalized_COV = normalize_sample(path_a=target_COV, path_b0=background_COV0, path_b1=background_COV1, method="FoldChange")
normalized_STARTS = normalize_sample(path_a=target_STARTS, path_b0=background_STARTS0, path_b1=background_STARTS1, method="FoldChange")

## add gzip 
if normalized_WPS.index.is_object():
    normalized_WPS.to_csv(output_WPS, sep="\t",header=None,compression='gzip')
else:
    normalized_WPS.to_csv(output_WPS, sep="\t",header=None, index=None,compression='gzip')
if normalized_WPS_v2.index.is_object():
    normalized_WPS_v2.to_csv(output_WPS_v2, sep="\t",header=None,compression='gzip')
else:
    normalized_WPS_v2.to_csv(output_WPS_v2, sep="\t",header=None, index=None,compression='gzip')

if normalized_COV.index.is_object():
    normalized_COV.to_csv(output_COV, sep="\t",header=None,compression='gzip')
else:
    normalized_COV.to_csv(output_COV, sep="\t",header=None, index=None,compression='gzip')

if normalized_STARTS.index.is_object():
    normalized_STARTS.to_csv(output_STARTS, sep="\t",header=None,compression='gzip')
else:
    normalized_STARTS.to_csv(output_STARTS, sep="\t",header=None, index=None,compression='gzip')




