#!/usr/bin/env python
import argparse
import os
from tqdm import tqdm
import pandas as pd

parser = argparse.ArgumentParser(description='Summarize featureCount output')
parser.add_argument('--indir', '-i', type=str, required=True,help='dir contains input featureCount result')
parser.add_argument('--output','-o',type=str, required=True,help='output count matrix')
parser.add_argument('--surfix','-s',type=str, default='.txt',required=False,help='sample featurecounts file surfix that will be stripped in sampleID (.txt, _promoter, ...) ' )
args = parser.parse_args()

records = []

files = [ each for each in os.listdir(args.indir) if not (each.endswith(".summary") | each.endswith("log"))]

print(f"Find {len(files)} samples .")

for file in tqdm(files):
    input = args.indir + "/" + file
    with open(input) as f:
        _ = next(f)
        _ = next(f)
        for line in f:
            line = line.strip()
            data = line.split("\t")
            gene_id,length,count = data[0],data[5],data[6]
            gene_id = gene_id + "|" + length
            sample_id = file.replace(args.surfix,"").strip()
            records.append((sample_id,gene_id,count))

counts = pd.DataFrame.from_records(records)
counts.columns = ["sample_id","gene_id","counts"]
matrix = counts.pivot(index="gene_id",columns="sample_id",values="counts")
matrix.to_csv(args.output,sep="\t")

