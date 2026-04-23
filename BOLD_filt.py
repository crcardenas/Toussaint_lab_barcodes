# if necessary
# conda create -n py3.10 anaconda::python=3.10 conda-forge::pandas

import os
import argparse
import pandas as pd

def get_best_species(tsv,top_n):
    sorted_tsv = tsv.sort_values(by=["Query ID", "ID%"], ascending=[True, False])

    records = sorted_tsv.groupby('Query ID', group_keys=False).head(top_n)
    return records
#def write_output_files():

def main():
    parser = argparse.ArgumentParser(description="Filter top match for samples BOLD hit")
    parser.add_argument("-i", "--input", required=True, help="Input TSV file")
    parser.add_argument("-o", "--output", default="output.tsv", help="Output file name")
    parser.add_argument( "--top", type=int, default=1, help="Number of top hits per Query ID")
    
    args = parser.parse_args()

    tsv = pd.read_csv(args.input, sep="\t")
    records = get_best_species(tsv,top_n=args.top)
    records.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()