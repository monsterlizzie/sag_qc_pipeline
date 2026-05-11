#!/usr/bin/env python3

import sys
import pandas as pd

QC_COLUMNS = [
    "Sample_ID",
    "Read_QC",
    "Assembly_QC",
    "Taxonomy_QC",
    "Mapping_QC",
    "Overall_QC",
    "Bases",
    "Contigs#",
    "Assembly_Length",
    "Seq_Depth",
    "Ref_Cov_%",
    "Het-SNP#",
    "SAG_Species",
    "SAG_Species_%",
    "Top_Non_SAG_Species",
    "Top_Non_SAG_Species_%"
]

def main():
    if len(sys.argv) < 3:
        sys.exit("Usage: generate_overall_report.py <output_csv> <sample_qc_csv...>")

    out_csv = sys.argv[1]
    qc_files = sys.argv[2:]

    dfs = []
    for f in qc_files:
        df = pd.read_csv(f, dtype=str)
        dfs.append(df)

    merged = pd.concat(dfs, ignore_index=True)

    for col in QC_COLUMNS:
        if col not in merged.columns:
            merged[col] = "NA"

    merged = merged[QC_COLUMNS]
    merged = merged.fillna("NA")
    merged = merged.sort_values("Sample_ID")

    merged.to_csv(out_csv, index=False)
    print(f"Wrote {out_csv} with {len(merged)} samples")

if __name__ == "__main__":
    main()