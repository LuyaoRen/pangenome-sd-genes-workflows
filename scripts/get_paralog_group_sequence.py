from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Get sequences of paralogs based on results from phylogeny tree", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-f", "--fasta", type=str, help="merged fasta file")
    parser.add_argument("-p", "--paralog", type=str, help="paralogs assigned file")
    parser.add_argument("-n", "--cluster", type=str, help="clsuter id")
    parser.add_argument("-o", "--output_director", type=str, help="output director")
    parser.add_argument("-c", "--file_check", type=str, help="output file check")
    args = parser.parse_args()

    # Load fasta
    record_dict = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))

    # Parse paralog assignment
    paralogs = {}
    with open(args.paralog, "r") as f:
        for line in f:
            strings = line.strip().split('\t')
            clade = strings[0]
            seq_id = strings[1][::-1].replace("_", ":", 1)[::-1].replace("_1_", "#1#").replace("_2_", "#2#")
            if strings[8] in ('0', 'NA'):
                if clade not in paralogs:
                    paralogs[clade] = {}
                if seq_id in record_dict:
                    paralogs[clade][seq_id] = record_dict[seq_id].seq
                else:
                    print(f"Warning: {seq_id} not found in fasta!")

    os.makedirs(args.output_director, exist_ok=True)

    for clade, paralog_seqs in paralogs.items():
        records = []
        chm13_seqs = {seq_id: sequence for seq_id, sequence in paralog_seqs.items() if seq_id.startswith("CHM13")}
        grch38_seqs = {seq_id: sequence for seq_id, sequence in paralog_seqs.items() if seq_id.startswith("GRCh38")}
        other_seqs = {seq_id: sequence for seq_id, sequence in paralog_seqs.items() if not seq_id.startswith("CHM13")}
        for seq_id, sequence in {**chm13_seqs, **grch38_seqs, **other_seqs}.items():
            record = SeqRecord(sequence, id=seq_id, description="")
            records.append(record)
        output_path = os.path.join(args.output_director, f"{args.cluster}_{clade}.fa")
        SeqIO.write(records, output_path, "fasta")

    list_path = os.path.join(args.output_director, f"{args.cluster}.fa.list")
    with open(list_path, "w") as fout:
        for clade in paralogs.keys():
            fa_path = os.path.join(args.output_director, f"{args.cluster}_{clade}.fa")
            fout.write(f"{fa_path}\n")

