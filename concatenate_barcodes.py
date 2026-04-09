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

def iter_concatenated_records(input_dir):
    fasta_files = sorted(glob.glob(os.path.join(input_dir, "*_barcode.fasta")))

    for fasta in fasta_files:
        sample_id = os.path.basename(fasta).rsplit("_")[0]
        software_source = os.path.basename(fasta).rsplit("_", 3)[1]
        for in_header, seq in iter_fasta_records(fasta):
            if not seq:
                continue
            yield f"{sample_id}_{software_source}_barcode", seq

def sample_id_from_concat_header(header):
    return header.split("_")[0].strip()

def filter_and_keep_longest_per_sample(records, min_len=329):
    longest = {}

    for header, seq in records:
        if len(seq) <= min_len:
            continue

        sample_id = sample_id_from_concat_header(header)
        prev = longest.get(sample_id)
        if prev is None or len(seq) > prev[0]:
            longest[sample_id] = (len(seq), header, seq)

    for sample_id in sorted(longest):
        _, header, seq = longest[sample_id]
        yield header, seq

def write_output_files(records, output_path, max_seqs=100):
    records = list(records)
    if not records:
        return

    base, ext = os.path.splitext(output_path)
    if not ext:
        ext = ".fasta"

    # if <= max_seqs → write single file exactly as given
    if len(records) <= max_seqs:
        with open(output_path, "w") as out:
            for header, seq in records:
                out.write(f">{header}\n{seq}\n")
        return

    # otherwise split into multiple files
    for i in range(0, len(records), max_seqs):
        file_number = (i // max_seqs) + 1
        out_file = f"{base}_{file_number}{ext}"
        with open(out_file, "w") as out:
            for header, seq in records[i:i+max_seqs]:
                out.write(f">{header}\n{seq}\n")

def main():
    parser = argparse.ArgumentParser(
        description="Concatenate *_barcode.fasta files and write FASTA (split if >100 sequences)."
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Directory containing *_barcode.fasta files."
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output FASTA file; which can include the path to a directory (e.g., ./example/test/concatenated_barcodes.fasta). You must include a file extension"
    )
    parser.add_argument(
        "-s","--subset", type=int, default=100,
        help="Max sequences per output file (default: 100). Always used by default"
    )
    parser.add_argument(
        "--longest_only", action="store_true",
        help="Keep only longest read per sample."
    )
    parser.add_argument(
        "--min_len", type=int, default=329,
        help="Can be set when using --longest_only flag is in use to set a minimum length to retain (default: 329)."
    )

    args = parser.parse_args()

    # ensure directory exists
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)

    records = iter_concatenated_records(args.input)

    if args.longest_only:
        records = filter_and_keep_longest_per_sample(records, min_len=args.min_len)

    write_output_files(records, args.output, max_seqs=args.subset)

    print(f"FASTA written based on {args.output}")

if __name__ == "__main__":
    main()
