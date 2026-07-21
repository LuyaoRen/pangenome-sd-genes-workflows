from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import argparse
import pandas as pd

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Extract shared regions of paralogs using reference genome (CHM13 and Chimpanzee), then based on their shared regions to extract shared MSA from samples", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-m", "--msa", type=str, help="msa fasta file")
	parser.add_argument("-o", "--output", help="output file name", type=str)
	parser.add_argument("-l", "--fragment_length", help="minimum contigous fragment length", type=int)
	parser.add_argument("-g", "--gene_name", help="gene used to locate intron region", type=str)
	parser.add_argument("-gl", "--gene_loc", help="location of gene used to locate intron region", type=str)


	args = parser.parse_args()

	sequences = [record.seq for record in SeqIO.parse(args.msa, "fasta") if "chr" in record.id or "NC_" in record.id]

	record_dict = SeqIO.to_dict(SeqIO.parse(args.msa, "fasta"))
	ref_seq = record_dict[args.gene_loc]
	chromosome = args.gene_loc.split(":")[0].split('_')[-1]
	ref_start_pos = int(args.gene_loc.split(":")[1].split("-")[0])
	pair = {}
	for index, letter in enumerate(ref_seq):
		if letter != "-":
			pair[ref_start_pos] = index+1
			ref_start_pos=ref_start_pos+1
		else:
			ref_start_pos = ref_start_pos

	intron = pd.read_table('/net/eichler/vol28/projects/sd_variants/nobackups/gene_annotation/chm13v2.0_RefSeq_Liftoff_v5.2_gene_model_intron.txt',header=None)
	gene_intron = intron[intron[4] == args.gene_name]

	sequences = [record.seq for record in SeqIO.parse(args.msa, "fasta") if "chr" in record.id]
	common_seq = ''
	for i in range(len(sequences[0])):
		chars_at_pos = [seq[i] for seq in sequences]
		if '-' in chars_at_pos:
			common_seq += '-'
		else:
			first_char = chars_at_pos[0]
			if all(char == first_char for char in chars_at_pos):
				common_seq += first_char
			else:
				common_seq += '-'
#common_seq = ''
#n_seq = len(sequences)
#threshold = int(round(n_seq * ((n_seq-1)/n_seq)))
#for i in range(len(sequences[0])):
#    chars_at_pos = [seq[i] for seq in sequences]
#    chars_no_gap = [c for c in chars_at_pos if c != '-']
#    if len(chars_no_gap) < threshold:
#        common_seq += '-'
#        continue
#    count = Counter(chars_no_gap)
#    most_common_base, count_most_common = count.most_common(1)[0]
#    if count_most_common >= threshold:
#        common_seq += most_common_base
#    else:
#        common_seq += '-'


	segments = []

	for index, row in gene_intron.iterrows():
		intron_segment = common_seq[pair[row[1]]-1:pair[row[2]]]
		segments.append(intron_segment)

	intron_segments = []

	for seg in segments:
		merged_segments = []
		current_segments = ''
		segs = re.split(r'(-+)', seg)
		for ss in segs:
			if ss.count('-') >= 50 and all(char == '-' for char in ss):
				merged_segments.append(current_segments)
				current_segments = ''
			else:
				current_segments += ss
		merged_segments.append(current_segments)
		intron_segments.append(merged_segments)


	intron_segments = [item for sublist in intron_segments for item in sublist]

	filtered_segments = [item for item in intron_segments if len(item) >= args.fragment_length]
	filtered_segments = [item for item in filtered_segments if item.count('-')/len(item) < 0.3]


	segments_pos = []
	for i in filtered_segments:
		start_pos = common_seq.find(i) + 1
		end_pos = start_pos + len(i) - 1
		seg_pos = [start_pos,end_pos]
		segments_pos.append(seg_pos)


	shared_region_fa = []
	for i in range(len(segments_pos)):
		sequence_str = sequences[0][segments_pos[i][0]-1:segments_pos[i][1]].replace('-','')
		sequence_record = SeqRecord(Seq(sequence_str), id=f"chm13.{i}", description="shared_region")
		shared_region_fa.append(sequence_record)

	SeqIO.write(shared_region_fa, args.output, "fasta")

