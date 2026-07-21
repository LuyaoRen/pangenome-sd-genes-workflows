from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Extract shared regions of paralogs using the CHM13 and priamtes reference genomes. Based on these shread regions, extract the corresponding shared regions from the samples", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-f", "--fasta", type=str, help="fasta file")
	parser.add_argument("-b", "--bed", type=str, help="bed file")
	parser.add_argument("-o", "--output", help="output file name", type=str)

	args = parser.parse_args()

#    shared_region_seq = ''
#    for record in SeqIO.parse(args.shared_fasta,"fasta"):
#        shared_region_seq += str(record.seq)

	sample_pos = {}
	for line in open(args.bed,"r"):
		strings = line.strip().split('\t')
		if strings[0] not in sample_pos:
			sample_pos[strings[0]] = {}
			sample_pos[strings[0]]['sequence'] = []
			sample_pos[strings[0]]['strand'] = []
		sample_pos[strings[0]]['sequence'].append([int(strings[1]),int(strings[2])])
		sample_pos[strings[0]]['strand'].append(strings[4])

	sample_shared_region_fa = []
	for record in SeqIO.parse(args.fasta, "fasta"):
		if record.id in sample_pos:
			pos_list = sample_pos[record.id]['sequence']
			strand_list = sample_pos[record.id]['strand']
			seq = ''
			for pos, strand in zip(pos_list, strand_list):
				subseq = record.seq[pos[0]:pos[1]]
				if strand == '-':
					subseq = subseq.reverse_complement()
				seq += subseq
			sequence_record = SeqRecord(seq, id=record.id)
			sample_shared_region_fa.append(sequence_record)
		else:
			pass

	SeqIO.write(sample_shared_region_fa,args.output, "fasta")
