import os
import glob
import argparse

#def concatenate_fasta(input_dir, output_file):
#    fasta_files = glob.glob(os.path.join(input_dir, "*_barcode.fasta"))

#    with open(output_file, "w") as outfile:
#        for fasta in fasta_files:
#            sample_id = os.path.basename(fasta).split("_")[0]  # e.g., "CBX0339"
#            with open(fasta, "r") as infile:
#                lines = infile.readlines()
#                if lines:
#                    header = lines[0].strip()  # first line (FASTA header)
#                    sequence = "".join(line.strip() for line in lines[1:])  # join sequence lines
#                    # Write in desired format
#                    outfile.write(f">{sample_id} | {header[1:]}\n")
#                    outfile.write(f"{sequence}\n")

def iter_fasta_records(path):
    header = None
    seq_parts = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        yield header, "".join(seq_parts)

def concatenate_fasta(input_dir, output_file):
    fasta_files = sorted(glob.glob(os.path.join(input_dir, "*_barcode.fasta")))

    with open(output_file, "w") as out:
        for fasta in fasta_files:
            sample_id = os.path.basename(fasta).rsplit("_")[0]
            sample_source =  os.path.basename(fasta).rsplit("_", 3)[1]
            for in_header, seq in iter_fasta_records(fasta):
                if not seq:
                    continue
                out.write(f">{sample_id} | {sample_source} {in_header}\n")
                out.write(f"{seq}\n")




def main():
    parser = argparse.ArgumentParser(
        description="Concatenate *_barcode.fasta files into one fasta."
    )
    parser.add_argument(
        "-i", "--input", required=True, 
        help="Directory containing *_barcode.fasta files."
    )
    parser.add_argument(
        "-p", "--prefix", default="barcode_loci.fasta",
        help="Name of the output fasta file (default: barcode_loci.fasta)."
    )
    parser.add_argument(
        "-o", "--output", default=".",
        help="Directory where the output file will be saved (default: current directory)."
    )

    args = parser.parse_args()

    # Ensure output directory exists
    os.makedirs(args.output, exist_ok=True)

    output_file = os.path.join(args.output, args.prefix)
    concatenate_fasta(args.input, output_file)

    print(f"Concatenated FASTA written to {output_file}")

if __name__ == "__main__":
    main()