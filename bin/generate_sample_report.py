#!/usr/bin/env python3

import sys
import pandas as pd

def main():
    if len(sys.argv) < 4:
        sys.exit("Usage: generate_sample_report.py <sample_id> <output_csv> <input_csvs...>")

    sample_id = sys.argv[1]
    out_csv = sys.argv[2]
    input_files = sys.argv[3:]

    dfs = []
    for f in input_files:
        df = pd.read_csv(f)
        dfs.append(df)

    # concatenate columns side-by-side
    merged = pd.concat(dfs, axis=1)

    # remove duplicate columns (happens when headers repeat)
    merged = merged.loc[:, ~merged.columns.duplicated()]

    # add Sample_ID column at front
    merged.insert(0, "Sample_ID", sample_id)

    merged.to_csv(out_csv, index=False)

if __name__ == "__main__":
    main()