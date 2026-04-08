import os
import glob
import argparse
import difflib
from Bio.Seq import Seq


def parse_fasta(fasta_file):
    sequences = {}
    header = None
    seq_chunks = []

    with open(fasta_file, "r") as record:
        for line in record:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    sequences[header] = "".join(seq_chunks)
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)

        if header is not None:
            sequences[header] = "".join(seq_chunks)

    return sequences

def parse_header(header):
    parts = header.split("_")
    sample_id = parts[0]
    seq_type = parts[1]
    return sample_id, seq_type

def calculate_similarity(seq1, seq2):
    seq3 = str(Seq(seq2).reverse_complement())
    similarity = difflib.SequenceMatcher(None, seq1, seq2).ratio()
    similarity_complement = difflib.SequenceMatcher(None, seq1, seq3).ratio()
    return similarity, similarity_complement

def compare_sequences(sequences):
    samples = {}

    for header, sequence in sequences.items():
        sample_id, seq_type = parse_header(header)

        if sample_id not in samples:
            samples[sample_id] = {"GO": None, "MF": None}

        samples[sample_id][seq_type] = sequence

    results = []

    for sample_id, seqs in samples.items():
        go_seq = seqs["GO"]
        mf_seq = seqs["MF"]

        if go_seq is not None and mf_seq is not None:
            similarity, similarity_complement = calculate_similarity(go_seq, mf_seq)
            results.append({
                "sample_id": sample_id,
                "similarity": similarity,
                "similarity_complement": similarity_complement,
                "GO_length": len(go_seq),
                "MF_length": len(mf_seq),
            })

    return results

def write_results(results, output_file):
    with open(output_file, "w") as out:
        out.write("sample_id\tsimilarity\tsimilarity_complement\tGO_length\tMF_length\n")
        for row in results:
            out.write(
                f"{row['sample_id']}\t"
                f"{row['similarity']:.4f}\t"
                f"{row['similarity_complement']:.4f}\t"
                f"{row['GO_length']}\t"
                f"{row['MF_length']}\n"
            )

def main():
    parser = argparse.ArgumentParser(
        description="Compare GO and MF sequences for each sample in a multi-FASTA file."
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Concatenated FASTA file."
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output file name"
    )

    args = parser.parse_args()

    sequences = parse_fasta(args.input)
    results = compare_sequences(sequences)

    output_file = f"{args.output}.tsv"
    write_results(results, output_file)

    print(f"Results written to {output_file}")


if __name__ == "__main__":
    main()
