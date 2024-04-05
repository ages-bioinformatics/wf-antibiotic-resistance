#!/usr/bin/env python
import argparse
import difflib
import pandas as pd
import sys
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--translation_table', help = "path/to/translation_table.tsv", required=True)
parser.add_argument('--tool', help = "AR tool [resfinder|amrfinder]", required=True)
parser.add_argument('--name', help = "MLST scheme name to be translated to species", required=True)

def main():
    args = parser.parse_args()
    df = pd.read_csv(args.translation_table, sep="\t")
    tool = difflib.get_close_matches(args.tool, df.columns, n=1)
    if len(tool) < 1:
        sys.stdout.write("")
        exit(0)
    if tool[0] != args.tool:
        sys.stderr.write(f"INFO: using column {tool}, requested tool: {args.tool}\n")
    df = df.replace([np.nan],"")
    translator = df.set_index("mlst_scheme_name")[tool[0]].to_dict()
    print(translator.get(args.name, ""))


if __name__=='__main__':
    main()
