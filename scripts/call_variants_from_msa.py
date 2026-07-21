from Bio import SeqIO
import pandas as pd
import numpy as np
import re
import argparse


#### functions
def remove_common_hyphens(A,B):
	new_A = ''
	new_B = ''
	for a,b in zip(A,B):
		if a != '-' or b != '-':
			new_A += a
			new_B += b
	return new_A, new_B


# call variants
def call_snv(ref,query):
	variants = []
	for i in range(len(ref)):
		#variants
		if i == 0 and (ref[i] == '-' or query[i] == '-'):
			pass
		# mismatch or SNV 
		if ref[i] != query [i] and ref[i] != '-' and query[i] != '-':
			variants.append({'CHROM': chromosome,
				'START': i + 1,
				'END': i + 1,
				'GT': ref[i].upper() + "/" + query[i].upper()})
	return pd.DataFrame(variants)	

def call_insertion(ref, query):
	variants = []
	insertion_start = -1
	insertion_detected = False
	for i in range(len(ref)):
		if ref[i] == '-' and query[i] != '-':
			if not insertion_detected:
				insertion_start = i
				insertion_detected = True
		elif insertion_detected and ref[i] != '-':
			alt_allele = query[insertion_start:i]
			variants.append({'CHROM': chromosome,
					'START': i + 1,
					'END': insertion_start,
					'GT': "-/" + alt_allele.upper()})
			insertion_detected = False
	if insertion_detected:
		alt_allele = query[insertion_start:len(query)]
		variants.append({'CHROM': chromosome,
					'START': i + 1,
					'END': insertion_start,
					'GT': "-/" + alt_allele.upper()})
	return pd.DataFrame(variants)


def call_deletion(ref, query):
	deletions = []
	deletion_start = -1
	deletion_detected = False
	for i in range(len(ref)):
		if query[i] == '-' and ref[i] != '-':
			if not deletion_detected:
				deletion_start = i
				deletion_detected = True
		elif deletion_detected and query[i] != '-':
			ref_allele = ref[deletion_start:i]
			deletions.append({'CHROM': chromosome,
					'START': deletion_start + 1,
					'END': i,
					'GT': ref_allele.upper() + '/-'})
			deletion_detected = False
	if deletion_detected:
		ref_allele = ref[deletion_start:len(query)]
		deletions.append({'CHROM': chromosome,
					'START': deletion_start + 1,
					'END': i,
					'GT': ref_allele.upper() + '/-'})
	return pd.DataFrame(deletions)


# convert coordinates from MSA to reference genome
def convert_cooredinates(variants,ref,ref_start_pos):
	pair = {}
	for index, letter in enumerate(ref):
		if letter != "-":
			pair[index+1] = ref_start_pos
			ref_start_pos=ref_start_pos+1
		else:
			ref_start_pos = ref_start_pos
	#the last deletion
	last_item = list(pair.items())[-1]
	pair[last_item[0]+1] = last_item[1] + 1
	#the first insertion
	pair[0] = list(pair.items())[0][1] - 1
	#
	variants['START'] = variants['START'].replace(pair)
	variants['END'] = variants['END'].replace(pair)
	# the last insertion
	if variants['START'].iloc[-1] == len(ref):
		variants['START'].iloc[-1] = variants['END'].iloc[-1] + 1
	return variants

# keep variants in exon regions
def is_in_any_region(start, end, exon_region):
	new_start = min(start, end)
	new_end = max(start, end)
	return ((new_start >= exon_region[1]) & (new_end <= exon_region[2])).any()

def filter_by_length(row):
	parts = row['GT'].split('/')
	return len(parts[0]) <= 50 and len(parts[1]) <= 50

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Call variants from MSA to default VEP input, it will use CHM13 sequence in the MSA as reference genome", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-m", "--msa", type=str, help="msa fasta file")
	parser.add_argument("-o", "--output", help="output file name", type=str)
	#parser.add_argument("-e", "--exonbed", help="exon_bed_region", type=str)

	args = parser.parse_args()

	## check if sequence has been reversed 
	record_ref_id = [record.id for record in SeqIO.parse(args.msa, "fasta") if "chr" in record.id]
	if "_R_" in record_ref_id[0]:
		msa_sequences = {}
		for record in SeqIO.parse(args.msa, "fasta"):
			record.seq = record.seq.reverse_complement()
			msa_sequences[record.id] = record.seq
			reverse = True
	else:
		msa_sequences = {}
		for record in SeqIO.parse(args.msa, "fasta"):
			msa_sequences[record.id] = record.seq
			reverse = False

	## call variants and convert coordinates
	ref_seq = msa_sequences[record_ref_id[0]]
	chromosome = record_ref_id[0].split(":")[0].split('_')[-1]
	ref_start_pos = int(record_ref_id[0].split(":")[1].split("-")[0])
	merged_vcf = pd.DataFrame(columns = ['CHROM', 'START', 'END','GT'])

	for seq_id,seq in msa_sequences.items():
		if seq_id == record_ref_id[0]:
			pass
		else:
			print(seq_id)
			ref, query = remove_common_hyphens(ref_seq,seq)
			variants = pd.concat([call_snv(ref,query), call_insertion(ref,query), call_deletion(ref,query)], axis=0, ignore_index=True)
			if variants.empty:
				pass
			else:
				variants_sorted = variants.sort_values(by='START')
				variants_sorted = convert_cooredinates(variants_sorted,ref,ref_start_pos)
				variants_sorted[seq_id.replace('_R_','')] = '1'
				merged_vcf = pd.merge(merged_vcf, variants_sorted, on=['CHROM', 'START', 'END','GT'], how='outer')
	merged_vcf = merged_vcf[merged_vcf.apply(filter_by_length, axis=1)]
	merged_vcf.fillna('0', inplace=True)
	merged_vcf.to_csv(args.output,index=False,sep="\t")
	# filter variants not in exon regions
	#exon_region = pd.read_table(args.exonbed,header=None)
	#filtered_variants_df = merged_vcf[merged_vcf.apply(lambda row: is_in_any_region(row['START'], row['END'], exon_region), axis=1)]
	#filtered_variants_df.fillna('0', inplace=True)
	#filtered_variants_df.to_csv(args.output,index=False,sep="\t")




