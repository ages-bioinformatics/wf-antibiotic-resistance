#!/usr/bin/env python

import sys
import pandas as pd
from Bio import SeqIO
from collections import defaultdict

assembly_file = sys.argv[1]
xfinder_result_files = sys.argv[2:]

ANNOTATED_OFFSET = 50000

dfs = []
for xfinder_result_file in xfinder_result_files:    
    df = pd.read_csv(xfinder_result_file, sep="\t")
    # transform resfinder output to same format as armfinder table output
    if "Contig" in df.columns:
        df[["Start","Stop"]] = df["Position in contig"].str.extract(r"(\d+)\.\.(\d+)", expand=True)

        df["Contig id"] = df["Contig"].astype(str).str.split(" ", expand=True)[0]
    df["Contig id"] = df["Contig id"].astype(str)
    dfs.append(df.copy())

df = pd.concat(dfs)

if df.empty:
    sys.stderr.write("No resistance genes found, exiting\n")
    exit(0)

df["lower_limit"] = df["Start"].astype(int).apply(lambda k: max(0, k - ANNOTATED_OFFSET))
df["upper_limit"] = df["Stop"].astype(int) + ANNOTATED_OFFSET

keep = defaultdict(list)
for key, subset in df.groupby("Contig id"):
    subset = subset.sort_values("lower_limit")
    start = subset["lower_limit"].iloc[0]
    stop = subset["upper_limit"].iloc[0]
    for i, row in subset.iterrows():
        if row["lower_limit"] < stop:
            stop = row["upper_limit"]
        else:
            keep[key].append((start, stop))
            start = row["lower_limit"]
            stop = row["upper_limit"]
    keep[key].append((start, stop))

kept_records = []
for record in SeqIO.parse(assembly_file, "fasta"):
    if record.id in keep:
        for start, stop in keep[record.id]:
            stop = min(stop, len(record))
            new_record = record[start:stop]
            new_record.id = record.id + f":{start}-{stop}"
            kept_records.append(new_record)

SeqIO.write(kept_records, sys.stdout, "fasta")
