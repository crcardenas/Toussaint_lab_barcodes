import os
import glob
import argparse

def iter_fasta_records(path):
    header = None
    seq_parts = []
    with open(path) as file:
        for line in file:
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
            software_source = os.path.basename(fasta).rsplit("_", 3)[1]
            for in_header, seq in iter_fasta_records(fasta):
                if not seq:
                    continue
                out.write(f">{sample_id}_{software_source}_barcode\n")
                out.write(f"{seq}\n")

def sample_id_from_concat_header(header):
    return header.split("|", 1)[0].strip()

def filter_and_keep_longest_per_sample(input_fasta, output_fasta, min_len=329):
    longest = {}  # sample_id -> (seq_len, header, seq)

    for header, seq in iter_fasta_records(input_fasta):
        if len(seq) <= min_len:
            continue

        sample_id = sample_id_from_concat_header(header)
        prev = longest.get(sample_id)
        if prev is None or len(seq) > prev[0]:
            longest[sample_id] = (len(seq), header, seq)

    with open(output_fasta, "w") as out:
        for sample_id in sorted(longest):
            _, header, seq = longest[sample_id]
            out.write(f">{header}\n")
            out.write(f"{seq}\n")

def main():
    parser = argparse.ArgumentParser(
        description="Concatenate *_barcode.fasta files, optionally filter and keep longest per sample."
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Directory containing *_barcode.fasta files."
    )
    parser.add_argument(
        "-p", "--prefix", required=True,
        help="Output FASTA filename for concatenated output (required)."
    )
    parser.add_argument(
        "-o", "--output", default=".",
        help="Directory where outputs will be saved (default: current directory)."
    )

    parser.add_argument(
        "--longest_only", action="store_true",
        help="After concatenation, remove reads less than the --min-len and keep the longest read per sample."
    )
    parser.add_argument(
        "--min_len", type=int, default=329,
        help="Minimum length to retain when using --longest-only (default: 329)."
    )
    parser.add_argument(
        "--filtered_prefix", required=True,
        help="Output FASTA filename for filtered output (default: barcode_loci_longest.fasta)."
    )

    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    output_file = os.path.join(args.output, args.prefix)
    concatenate_fasta(args.input, output_file)
    print(f"Concatenated FASTA written to {output_file}")

    if args.longest_only:
        filtered_file = os.path.join(args.output, args.filtered_prefix)
        filter_and_keep_longest_per_sample(output_file, filtered_file, min_len=args.min_len)
        print(f"Filtered FASTA written to {filtered_file}")

if __name__ == "__main__":
    main()
