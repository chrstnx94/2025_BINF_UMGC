# BISF 617 Final Project
# Team Members: Christian Figueroa-Perez, Simbiat Yusuf

import os
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt

def find_orfs(sequence, min_length=100):
    # Team Member Name: Simbiat Yusuf, Christian Figueroa-Perez

    # Find ORFs in all 6 reading frames of a sequence
    orfs = []
    seq_len = len(sequence)

    for strand, nuc in [(+1, sequence), (-1, sequence.reverse_complement())]:
        for frame in range(3):
            # Slice the sequence from the current frame
            sub_seq = nuc[frame:]
            # Trim to full codons
            codon_len = len(sub_seq) - (len(sub_seq) % 3)
            sub_seq = sub_seq[:codon_len]
            # Translate
            trans = str(sub_seq.translate(to_stop=False))
            trans_len = len(trans)

            aa_start = 0
            while aa_start < trans_len:
                aa_start = trans.find("M", aa_start)
                if aa_start == -1:
                    break
                aa_end = trans.find("*", aa_start)
                if aa_end == -1:
                    break
                orf_len = aa_end - aa_start
                if orf_len * 3 >= min_length:
                    start = frame + aa_start * 3
                    end = frame + aa_end * 3 + 3
                    if strand == -1:
                        start, end = seq_len - end, seq_len - start
                    aa_seq = trans[aa_start:aa_end]
                    orfs.append((start, end, strand, aa_seq))
                aa_start = aa_end + 1

    return orfs

def process_fasta(input_file, output_file, min_length=100):
    # Team Member Name: Simbiat Yusuf, Christian Figueroa-Perez

    # Parse input file (FASTA), find all ORFs ≥ min_length,
    # write them as a FASTA to output file, and print to terminal.
    # Returns a summary dict for visualization.

    summary = {}

    with open(output_file, 'w') as out_f:
        for record in SeqIO.parse(input_file, "fasta"):
            orfs = find_orfs(record.seq, min_length)
            fwd = [o for o in orfs if o[2] == +1]
            rev = [o for o in orfs if o[2] == -1]
            summary[record.id] = {'forward': fwd, 'reverse': rev}

            print(f"\n> {record.id}")
            print(f"{len(fwd)} ORFs in forward, {len(rev)} in reverse\n")

            for i, (start, end, strand, aa_seq) in enumerate(orfs, 1):
                strand_sym = '+' if strand == +1 else '-'
                header = f"{record.id}_orf{i} | {start}-{end} | strand:{strand_sym}"
                formatted = f">{header}\n{aa_seq}"
                print(formatted + "\n")         # Print to terminal
                out_f.write(formatted + "\n")   # Write to file

    return summary


def create_visualization(orf_summary):
    # Team Member Name: Christian Figueroa-Perez

    # Given a dict {seq_id: {'forward': [...], 'reverse': [...]}, ...},
    # plot counts of forward vs. reverse ORFs per sequence.

    vis_dir = os.path.join("output", "visualization")
    os.makedirs(vis_dir, exist_ok=True)
    output_path = os.path.join(vis_dir, "orf_visualization.png")

    labels = []
    fwd_counts = []
    rev_counts = []

    for seq_id, data in orf_summary.items():
        labels.append(seq_id)
        fwd_counts.append(len(data['forward']))
        rev_counts.append(len(data['reverse']))

    if not labels:
        print("No ORFs found; skipping visualization.")
        return

    x = range(len(labels))
    plt.figure(figsize=(10, 6))
    plt.bar(x, fwd_counts, width=0.4, label="Forward ORFs", color="skyblue", align="center")
    plt.bar(x, rev_counts, width=0.4, label="Reverse ORFs", color="salmon",  align="edge")
    plt.xticks(x, labels, rotation=45, ha="right")
    plt.ylabel("ORF count")
    plt.title("Six‐frame ORF Counts per Sequence")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path)
    print(f"Saved visualization to {output_path}")


def main():
    # Team Member Name: Christian Figueroa-Perez
    fasta_path = input("Path to input FASTA: ").strip()
    default_len = 100
    prompt = f"Minimum ORF length in nucleotides (default={default_len}): "
    ml = input(prompt).strip()
    min_length = int(ml) if ml.isdigit() else default_len

    # Ensure output dirs exist
    orf_dir = os.path.join("output", "orfs")
    os.makedirs(orf_dir, exist_ok=True)
    orf_fasta = os.path.join(orf_dir, "orf_output.fasta")

    # Find, write ORFs, and get summary
    summary = process_fasta(fasta_path, orf_fasta, min_length)
    print(f"Wrote ORFs to {orf_fasta}")

    # Visualize counts
    create_visualization(summary)


if __name__ == "__main__":
    main()
