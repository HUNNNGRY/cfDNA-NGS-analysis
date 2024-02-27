import pandas as pd
import numpy as np
from scipy import stats
import gzip


target_WPS = snakemake.input["target_WPS"]
background_WPS = snakemake.input["background_WPS"]

target_COV = snakemake.input["target_COV"]
background_COV = snakemake.input["background_COV"]

target_STARTS = snakemake.input["target_STARTS"]
background_STARTS = snakemake.input["background_STARTS"]

output_WPS = snakemake.output["output_WPS"]
output_COV = snakemake.output["output_COV"]
output_STARTS = snakemake.output["output_STARTS"]

index_choice=None

# add if WPS then delta not FC 
def normalize_sample(path_a: str, path_b: str, method: str="FoldChange"):
    """Reads .csv file, calculates mean over all rows and divides by trimmed mean.

    Args:
        path (str): Path to a .csv file

    Returns:
        [type]: [description]
    """
    #with open(path_a, 'rb') as fd:
    #    gzip_fd = gzip.GzipFile(fileobj=fd)
    sample_a = pd.read_csv(path_a,compression='gzip', encoding='ISO-8859–1', header=None) 
    if pd.api.types.is_string_dtype(sample_a[0]):
        sample_a = sample_a.set_index(0)
        sample_a.index.name="ID"
    sample_b = pd.read_csv(path_b,compression='gzip', encoding='ISO-8859–1', header=None).mean(axis=1)
    if method == "FoldChange":
        sample = sample_a / stats.trim_mean(sample_b, 0.1) # + 0.1 # pseudo count
    else:
	#print("WPS using Delta")  
        #sys.stderr.write("using Delta normalize...\n")
        sample = sample_a - stats.trim_mean(sample_b, 0.1)   
    return sample.round(4)

#normalized_WPS = normalize_sample(path_a=target_WPS, path_b=background_WPS, method="FoldChange")
normalized_WPS = normalize_sample(path_a=target_WPS, path_b=background_WPS, method="Delta")
normalized_COV = normalize_sample(path_a=target_COV, path_b=background_COV, method="FoldChange")
normalized_STARTS = normalize_sample(path_a=target_STARTS, path_b=background_STARTS, method="FoldChange")

## add gzip 
if normalized_WPS.index.is_object():
    normalized_WPS.to_csv(output_WPS, sep="\t",header=None,compression='gzip')
else:
    normalized_WPS.to_csv(output_WPS, sep="\t",header=None, index=None,compression='gzip')

if normalized_COV.index.is_object():
    normalized_COV.to_csv(output_COV, sep="\t",header=None,compression='gzip')
else:
    normalized_COV.to_csv(output_COV, sep="\t",header=None, index=None,compression='gzip')

if normalized_STARTS.index.is_object():
    normalized_STARTS.to_csv(output_STARTS, sep="\t",header=None,compression='gzip')
else:
    normalized_STARTS.to_csv(output_STARTS, sep="\t",header=None, index=None,compression='gzip')




