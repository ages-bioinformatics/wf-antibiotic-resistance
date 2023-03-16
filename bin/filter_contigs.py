#!/usr/bin/env python

import sys
import pandas as pd
from Bio import SeqIO

resfinder_result = sys.argv[1]
assembly_file = sys.argv[2]

df = pd.read_csv(resfinder_result, sep="\t")
if df.empty:
    sys.stderr.write("No resistance genes found, exiting\n")
    exit(0)
df["Contig"] = df["Contig"].str.split(" ", expand=True)[0]
contigs_with_resistance_genes = set(df["Contig"].unique())

kept_records = []
for record in SeqIO.parse(assembly_file, "fasta"):
    if record.id in contigs_with_resistance_genes:
        kept_records.append(record)


SeqIO.write(kept_records, sys.stdout, "fasta")
