import pandas as pd
from io import StringIO
import argparse

def pick_consequence(val):
    val = val.replace('splice_region_variant,', '').replace(',splice_region_variant', '').replace('splice_region_variant', '')
    parts = val.split(',')
    if 'missense_variant' in parts and len(parts) > 1:
        return [p for p in parts if p != 'missense_variant'][0]
    return parts[0]

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Extract shared regions of paralogs using reference genome (CHM13 and Chimpanzee), then based on their shared regions to extract shared MSA from samples", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-f", "--variants", type=str, help="*.variants.txt")
	parser.add_argument("-v", "--vep_output", type=str, help="*.variants.vep_output.txt")
	parser.add_argument("-p", "--paralog_group", type=str, help="clade_information.txt")
	parser.add_argument("-s", "--snp_genotype", type=str, help="*.snp_genotype.txt")
	parser.add_argument("-l", "--lollipop", type=str, help="*.lollipop.txt")
	parser.add_argument("-d", "--detailed_lollipop", type=str, help="*.detailed.lollipop.txt")

	args = parser.parse_args()
	clade = args.variants.split('/')[-1].split('_')[2].split('.')[0]
	clade_df = pd.read_table(args.paralog_group,header=None)
	#subset_df = clade_df[(clade_df[0] == clade) & (clade_df[1].str.contains("CHM13"))]
	#gene_name = subset_df[5].iloc[0]
	subset_df = clade_df[(clade_df[0] == clade) & (clade_df[1].str.contains("CHM13", na=False))]
	if len(subset_df) > 0:
		gene_name = subset_df[5].iloc[0]
	else:
		subset_df = clade_df[(clade_df[0] == clade) & (clade_df[1].str.contains("GRCh38", na=False))]
		if len(subset_df) > 0:
			gene_name = subset_df[1].iloc[0].replace("GRCh38_", "", 1)
		else:
			gene_name = clade_df.loc[clade_df[0] == clade, 5].iloc[0]

	freq = pd.read_table(args.variants)
#	selected_columns = list(freq.columns[:4]) + [col for col in freq.columns if "s1" in col or "s2" in col]
#	freq = freq[selected_columns]

	num = []
	for index2,gt in freq.iloc[:, 4:].iterrows():
		count_1 = sum(gt == 1)
		num.append(count_1)

	variant_id = freq['CHROM'].astype(str) + "_" + freq['START'].astype(str) + "_" + freq['GT'].astype(str)

	freq_sum = pd.DataFrame({
	    'variant_id': variant_id,
	    'num': num
	})

	freq_sum['hap'] = freq.apply(lambda row: ", ".join([col for col in freq.columns if row[col] == 1]), axis=1)

	with open(args.vep_output, 'r') as f:
	    lines = f.readlines()

	# Filter out lines that start with "##"
	filtered_lines = [line for line in lines if not line.startswith("##")]

	# Convert the filtered lines into a pandas DataFrame
	vep = pd.read_csv(StringIO(''.join(filtered_lines)), delimiter='\t')

	freq['variant_id'] = freq['CHROM'].astype(str) + "_" + freq['START'].astype(str) + "_" + freq['GT'].astype(str)
	gt = freq.merge(vep, left_on="variant_id", right_on="#Uploaded_variation", how="right")

	gt = gt[~gt["GT"].str.contains("-", na=False)]
	gt = gt[gt['CDS_position'] != '-']

	rows = []
	for index, row in gt.iterrows():
	    if '/' not in row['GT']:
	        continue
	    ref, alt = row['GT'].split('/')
	    row['GT'] = ref
	    for i in range(4, freq.shape[1]):
	        if row.iloc[i] == 0:
	            row.iloc[i] = ref
	        elif row.iloc[i] == 1:
	            row.iloc[i] = alt
	    rows.append(row.copy())

	new_df = pd.DataFrame(rows)
	if new_df.empty:
		print("Warning: new_df is empty. Writing header only.")
		empty_df = pd.DataFrame(columns=["##CHROM", "POS", "REF"])
		empty_df.to_csv(args.snp_genotype, index=None, sep="\t")
	else:
		columns_to_keep = [0, 1] + list(range(3, freq.shape[1] - 1))
		new_df = new_df.iloc[:, columns_to_keep]
		new_df.rename(columns={new_df.columns[0]: "##CHROM",
		                       new_df.columns[1]: "POS",
		                       new_df.columns[2]: "REF"}, inplace=True)
		new_df.to_csv(args.snp_genotype, index=None, sep="\t")


	df = freq_sum.merge(vep, left_on="variant_id", right_on="#Uploaded_variation", how="right")
	df['tag'] = '-'
	df.loc[df['Extra'].str.contains('likely_pathogenic|probably_damaging|possibly_damaging', na=False), 'tag'] = 'predicted_deleterious'
	df.loc[df['Extra'].str.contains('benign', na=False), 'tag'] = '-'

	sub_df = df[df['Gene'] == gene_name]
	sub_df = sub_df[
	(sub_df['CDS_position'] != '-') | 
	(sub_df['Extra'].str.contains('IMPACT=HIGH', na=False))
	]
	sub_df['CHROM'] = sub_df['Location'].str.split(':').str[0]
	sub_df['POS'] = sub_df['Location'].str.split(':').str[1]
	sub_df = sub_df[['CHROM','POS','Consequence','tag','Protein_position','Amino_acids','num','hap']]
	sub_df = sub_df[sub_df["num"] != 0]
	sub_df['Consequence'] = sub_df['Consequence'].apply(pick_consequence)
	sub_df.loc[sub_df['Consequence'].str.contains('stop_gained|start_lost|stop_lost|frameshift_variant|splice_donor_variant|splice_acceptor_variant', na=False), 'tag'] = 'deleterious'
	sub_df['POS'] = sub_df['POS'].apply(lambda x: str(x).split('-')[0] if isinstance(x, str) else x)
	long_df = (
	    sub_df.assign(hap=sub_df["hap"].str.split(","))
	    .explode("hap")
	    .reset_index(drop=True)
	)
	long_df['hap'] = long_df['hap'].replace({'#': '_', ':': '_'}, regex=True)
	sub_df.to_csv(args.lollipop,sep="\t",index=None)
	long_df.to_csv(args.detailed_lollipop,sep="\t",index=None)
