import argparse
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


def extract_kmers(sequence, kmer_size):
    kmers = []
    n_kmers = len(sequence) - kmer_size + 1

    for i in range(n_kmers):
        kmer = sequence[i:i + kmer_size]
        kmers.append(kmer)

    return kmers


def jaccard_similarity(seq_a, seq_b):
    set_a = set(seq_a)
    set_b = set(seq_b)

    intersection = len(set_a.intersection(set_b))
    union = len(set_a.union(set_b))

    if union == 0:
        return 0.0

    return intersection / union


def compare_sequences(sequences, kmer_size):
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
            mf_seq_revcomp = str(Seq(seqs["MF"]).reverse_complement())

            go_kmers = extract_kmers(go_seq, kmer_size)
            mf_kmers = extract_kmers(mf_seq, kmer_size)
            mf_revcomp_kmers = extract_kmers(mf_seq_revcomp, kmer_size)

            similarity = jaccard_similarity(go_kmers, mf_kmers)
            similairty_revcomp = jaccard_similarity(go_kmers, mf_revcomp_kmers)

            results.append({
                "sample_id": sample_id,
                "kmer_jaccard_similarity": similarity,
                "kmer_jaccard_sim_revcomp": similairty_revcomp,
                "GO_length": len(go_seq),
                "MF_length": len(mf_seq),
            })

    return results


def write_results(results, output_file):
    with open(output_file, "w") as out:
        out.write("sample_id\tkmer_jaccard_similarity\tkmer_jaccard_sim_revcom\tGO_length\tMF_length\n")
        for row in results:
            out.write(
                f"{row['sample_id']}\t"
                f"{row['kmer_jaccard_similarity']:.4f}\t"
                f"{row['kmer_jaccard_sim_revcomp']:.4f}\t"
                f"{row['GO_length']}\t"
                f"{row['MF_length']}\n"
            )


def main():
    parser = argparse.ArgumentParser(
        description="Compare GO and MF sequences for each sample using kmer Jaccard similarity."
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Concatenated FASTA file."
    )
    parser.add_argument(
        "-k", "--kmer_size", type=int, default=21,
        help="kmer size, typical sizes used: 21, 33 and 55"
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output file name"
    )

    args = parser.parse_args()

    sequences = parse_fasta(args.input)
    results = compare_sequences(sequences, args.kmer_size)
    write_results(results, args.output)


if __name__ == "__main__":
    main()
