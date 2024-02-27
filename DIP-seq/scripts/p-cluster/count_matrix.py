#!/usr/bin/env python
# -*- coding: utf-8 -*-
# older DIP-seq version

#参考https://divingintogeneticsandgenomics.rbind.io/post/merge-featurecount-table-from-rnaseq/
import os
import csv
import glob
import sys

regions = ["gene","promoter","CGI"]     #实际上这三种可以选择性计算。后续可添加到config文件中进行自由选择
for region in regions:
    in_path = ""
    outdir = sys.argv[1]
    input_list = [outdir,"/07counts/*",region]
    files = glob.glob(in_path.join(input_list))     #找到多少count文件就用多少，而不一定必须找到对应每个样本（expand）
    list_column = []
    n = 1
    for file in files:
        #print n
        column_data = []
        with open(file, 'r') as f:
            reader = csv.reader(f, delimiter = "\t")
                # skip the comment line
            comment = next(reader)
            if n <= 1:
                for row in reader:
                    # for the first file, keep the gene column as well
                    m = 1
                    if m <= 1:
                        sample_name = row[6].split("/")[-1].split(".")[0]
                        print (sample_name)
                        column_data.append(row[0] + '\t' + sample_name)
                    else:
                        column_data.append(row[0] + '\t' + row[6])
                    m = m + 1
            else:
                for row in reader:
                    m = 1
                    if m <= 1:
                        sample_name = row[6].split("/")[-1].split(".")[0]
                        print (sample_name)
                        column_data.append(sample_name)
                    else:
                        column_data.append(row[6])
                    m = m + 1
            n = n + 1
        list_column.append(column_data)
        
    # This creates a list of row lists from the list of column lists
    # If any of the column lists are too short they will be padded with None
    # map function is a gem :)
    rows = map(None, *list_column)
    out_path = ""
    output_list = [outdir,"/07counts/count_matrix_",region,".txt"]

    with open(out_path.join(output_list),'w') as f_out:
        for row in rows:
            f_out.write('\t'.join(row))
            f_out.write('\n')
