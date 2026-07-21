import pandas as pd

configfile: "config.yaml"

manifest_path = config['tbl']
human_gene_dict = config['human_gene']
human_ref = config['human_ref']

cluster_gene_map = {
    (cluster, gene): pos
    for cluster, genes in human_gene_dict.items()
    for gene, pos in genes.items()
}

manifest_df = pd.read_csv(manifest_path, sep="\t", index_col='sample')
samples = manifest_df.index.tolist()

clusters = list(human_gene_dict.keys())
genes_per_cluster = {c: list(g.keys()) for c, g in human_gene_dict.items()}

valid_clusters = clusters
valid_genes = sorted({gene for genes in human_gene_dict.values() for gene in genes})
valid_samples = manifest_df.index.tolist()

wildcard_constraints:
    cluster = "|".join(valid_clusters),
    gene = "|".join(valid_genes),
    sample = "|".join(valid_samples)

def find_asm(wildcards):
    return manifest_df.at[wildcards.sample, 'asm']

def find_flagger(wildcards):
    return manifest_df.at[wildcards.sample, 'flagger']

def find_human_search_pos(wildcards):
    gene_pos = cluster_gene_map[(wildcards.cluster, wildcards.gene)]
    chrom, coords = gene_pos.split(":")
    start, end = coords.split("-")
    return f"{chrom}:{start}-{end}"
        
rule all:
    input:
        gene_fas=[
            f"results/{cluster}/gene_sequence/{gene}.chm13.fa"
            for cluster in clusters
            for gene in genes_per_cluster[cluster]
        ],
        merged_pafs=[
            f"results/{cluster}/paralogs_bed/{sample}.paf"
            for cluster in clusters
            for sample in samples
        ],
        trees=[
            f"results/{cluster}/msa/msa.shared_region.trimal.fa.treefile"
            for cluster in clusters
        ],
        clades=[
            f"results/{cluster}/tree/{cluster}_clade_information.txt"
            for cluster in clusters
        ],

rule get_human_gene_seq:
    input:
        ref=human_ref,
    output:
        fa="results/{cluster}/gene_sequence/{gene}.chm13.fa",
    params:
        gene_pos_opts=find_human_search_pos,
    threads: 1
    resources:
        hrs = 12,
        mem = 4,
    shell:
        """
        samtools faidx {input.ref} {params.gene_pos_opts} > {output.fa}
        """

rule merge_gene_seq:
    input:
        fa=lambda wildcards: [
            f"results/{wildcards.cluster}/gene_sequence/{gene}.chm13.fa"
            for gene in genes_per_cluster[wildcards.cluster]
        ],
    output:
        merged_fa="results/{cluster}/gene_sequence/all_genes_merged.chm13.fa",
    threads: 1
    resources:
        hrs = 8,
        mem = 2,
    shell:
        """
        cat {input.fa} > {output.merged_fa}
        """

rule shared_region:
    input:
        chm13_fa="results/{cluster}/gene_sequence/all_genes_merged.chm13.fa",
    output:
        msa="results/{cluster}/gene_sequence/chm13.msa.fa",
        shared_region_fa="results/{cluster}/gene_sequence/chm13.shared_region.fa",
        chm13_paf="results/{cluster}/gene_sequence/chm13.shared_region.paf",
        bed="results/{cluster}/gene_sequence/chm13.shared_region.bed",
    params:
        fragment_length=config['fragment_length'],
        gene_name=lambda wildcards: list(human_gene_dict[wildcards.cluster].keys())[0],
        gene_loc=lambda wildcards: human_gene_dict[wildcards.cluster][list(human_gene_dict[wildcards.cluster].keys())[0]],
    threads: 12
    resources:
        hrs = 24,
        mem = 16,
    shell:
        """
        mafft --adjustdirection --thread 12 --auto --reorder {input.chm13_fa} > {output.msa}
        python scripts/big_gene_shared_region_only_intron.py -m {output.msa} -o {output.shared_region_fa} -l {params.fragment_length} -g {params.gene_name} -gl {params.gene_loc}
        python scripts/get_shared_paf.py -s {output.shared_region_fa} -t {input.chm13_fa} -p {output.chm13_paf}
        #cat {output.chm13_paf} | awk '{{print $6"\t"$8"\t"$9"\t"$1"\t"$5}}' > {output.bed}
        awk 'BEGIN{{OFS="\t"}} {{print $6, $8, $9, $1, $5}}' {output.chm13_paf} | \
        awk 'BEGIN{{FS=OFS="\t"}} {{
            key = $1 FS $4
            count[key]++
            lines[NR] = $0
            keys[NR] = key
        }} END {{
            for (i = 1; i <= NR; i++) {{
                if (count[keys[i]] == 1)
                    print lines[i]
            }}
        }}' > {output.bed}
        """

rule extract_sequence:
    input:
        asm=find_asm,
        query="results/{cluster}/gene_sequence/all_genes_merged.chm13.fa",
    output:
        paf="results/{cluster}/paralogs_bed/{sample}.paf",
    threads: 4
    resources:
        hrs = 24,
        mem = 8,
    shell:
        """
        minimap2 -x asm20 -c --secondary=yes -p 0.3 -N 10000 --eqx -t 4 -r 500 -K 500M {input.asm} {input.query} > {output.paf}
        """

rule paralogs_loc:
    input:
        paf="results/{cluster}/paralogs_bed/{sample}.paf",
        flagger_bed=find_flagger,
    output:
        bed="results/{cluster}/paralogs_bed/{sample}.bed",
        filtered_bed="results/{cluster}/paralogs_bed/{sample}.filtered.bed",
        flagger_annotated_bed="results/{cluster}/paralogs_bed/{sample}.filtered.flagger.bed",
    threads: 4
    resources:
        hrs = 24,
        mem = 8,
    shell:
        """
        if [ -s {input.paf} ]; then
            rb stat --paf {input.paf} > {output.bed}
            python scripts/best_mapping_human_primates_rustybam.py -b {output.bed} -f {output.filtered_bed}
            if [ -s {input.flagger_bed} ]; then
                if [ -s {output.filtered_bed} ]; then
                    bedtools annotate -i {output.filtered_bed} -files {input.flagger_bed} > {output.flagger_annotated_bed}
                else
                    touch {output.flagger_annotated_bed}
                fi
            else
                touch {output.flagger_annotated_bed}
            fi
        else
            touch {output.bed}; touch {output.filtered_bed}; touch {output.flagger_annotated_bed}
        fi
        """

rule get_sequence:
    input:
        asm=find_asm,
        bed=rules.paralogs_loc.output.filtered_bed,
        shared_fa="results/{cluster}/gene_sequence/chm13.shared_region.fa",
    output:
        fa="results/{cluster}/get_sequence/{sample}.fa",
        paf="results/{cluster}/get_sequence/{sample}.shared.paf",
        bed="results/{cluster}/get_sequence/{sample}.shared.bed",
    threads: 4
    resources:
        hrs = 24,
        mem = 8,
    shell:
        """
        if [ -s {input.bed} ]; then
            samtools faidx -r <(awk '{{print $1":"$2"-"$3}}' {input.bed}) {input.asm} | awk -v s={wildcards.sample} '/^>/{{$0=">"s"_"substr($0,2)}}1' > {output.fa}
            python scripts/get_shared_paf.py -s {input.shared_fa} -t {output.fa} -p {output.paf}
            #cat {output.paf} | awk '{{print $6"\t"$8"\t"$9"\t"$1"\t"$5}}' > {output.bed}
            awk 'BEGIN{{OFS="\t"}} {{print $6, $8, $9, $1, $5}}' {output.paf} | \
            awk 'BEGIN{{FS=OFS="\t"}} {{
                key = $1 FS $4
                count[key]++
                lines[NR] = $0
                keys[NR] = key
            }} END {{
                for (i = 1; i <= NR; i++) {{
                    if (count[keys[i]] == 1)
                        print lines[i]
                }}
            }}' > {output.bed}
        else
            touch {output.fa}; touch {output.paf}; touch {output.bed}
        fi
        """

rule mafft_msa:
    input:
        sample_bed=[
            f"results/{{cluster}}/get_sequence/{sample}.shared.bed" for sample in samples
        ],
        sample_fa=[
            f"results/{{cluster}}/get_sequence/{sample}.fa" for sample in samples
        ],
        flagger_bed=[
            f"results/{{cluster}}/paralogs_bed/{sample}.filtered.flagger.bed" for sample in samples
        ],
        filtered_bed=[
            f"results/{{cluster}}/paralogs_bed/{sample}.filtered.bed" for sample in samples
        ],
    output:
        merged_bed="results/{cluster}/msa/sample.chm13.chimp.bed",
        merged_fa="results/{cluster}/msa/sample.chm13.chimp.merged.fa",
        shared_region_merged="results/{cluster}/msa/sample.chm13.chimp.shared_region.merged.fa",
        msa="results/{cluster}/msa/msa.shared_region.fa",
        trimmed_msa="results/{cluster}/msa/msa.shared_region.trimal.fa",
        flagger_merged_bed="results/{cluster}/msa/flagger.merged.bed",
        filtered_merged_bed="results/{cluster}/msa/filtered.merged.bed",
    threads: 12
    resources:
        hrs = 24,
        mem = 12,
    shell:
        """
        cat {input.filtered_bed} > {output.filtered_merged_bed}
        cat {input.flagger_bed} > {output.flagger_merged_bed}
        cat {input.sample_fa} > {output.merged_fa}
        cat {input.sample_bed} > {output.merged_bed}
        python scripts/merge_shared_region_fasta.py -f {output.merged_fa} -b {output.merged_bed} -o {output.shared_region_merged}
        mafft --adjustdirection --thread 12 --auto --reorder {output.shared_region_merged} > {output.msa}
        trimal -in {output.msa} -out {output.trimmed_msa} -gappyout -keepheader
        """

rule phylo_tree:
    input:
        msa="results/{cluster}/msa/msa.shared_region.trimal.fa",
    output:
        model="results/{cluster}/msa/msa.shared_region.trimal.fa.model.gz",
        mldist="results/{cluster}/msa/msa.shared_region.trimal.fa.mldist",
        bionj="results/{cluster}/msa/msa.shared_region.trimal.fa.bionj",
        splits_nex="results/{cluster}/msa/msa.shared_region.trimal.fa.splits.nex",
        contree="results/{cluster}/msa/msa.shared_region.trimal.fa.contree",
        iqtree="results/{cluster}/msa/msa.shared_region.trimal.fa.iqtree",
        treefile="results/{cluster}/msa/msa.shared_region.trimal.fa.treefile",
        ckp="results/{cluster}/msa/msa.shared_region.trimal.fa.ckp.gz",
        log="results/{cluster}/msa/msa.shared_region.trimal.fa.log",
    threads: 12
    resources:
        hrs = 24,
        mem = 12,
    shell:
        """
        iqtree2 -T AUTO -s {input.msa} -m MFP -B 1000 --threads-max 12
        """

rule phylo_plot:
    input:
        tree_file="results/{cluster}/msa/msa.shared_region.trimal.fa.treefile",
        paralog_file="results/{cluster}/msa/filtered.merged.bed",
        gene_bed=config['gene_bed'],
        flagger_file="results/{cluster}/msa/flagger.merged.bed",
        region_file="results/{cluster}/gene_sequence/chm13.shared_region.bed",
    output:
        p1="results/{cluster}/tree/{cluster}_evolutionary_tree.pdf",
        p2="results/{cluster}/tree/{cluster}_hprc_hgsvc_tree.pdf",
        clade_info="results/{cluster}/tree/{cluster}_clade_information.txt",
        tree="results/{cluster}/tree/{cluster}_rooted.tree",
        cn_count="results/{cluster}/tree/{cluster}_cn_counts.txt",
        all_paralog="results/{cluster}/tree/{cluster}_all_paralogs.txt",
    threads: 1
    resources:
        hrs = 8,
        mem = 4,
    shell:
        """
        Rscript scripts/tree_plot.R --tree {input.tree_file} \
        --paralog {input.paralog_file} \
        --gene_bed {input.gene_bed} \
        --flagger {input.flagger_file} \
        --region {input.region_file} \
        --p1_output {output.p1} \
        --p2_output {output.p2} \
        --summary {output.clade_info} \
        --rooted_tree {output.tree} \
        --cn_count {output.cn_count} \
        --all_paralog {output.all_paralog}
        """