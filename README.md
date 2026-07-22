# pangenome-sd-genes-workflows

Snakemake workflows for constructing phylogenetic trees from the pangenome (`phylo_tree`) and identifying and classifying variants for each paralog (`paralog_variants`).

## Overview

This repository contains two linked Snakemake pipelines for characterizing paralogous SD genes located within segmental duplications (SDs) across a panel of phased genome assemblies (HPRC/HGSVC):

1. **`phylo_tree`** — extracts each paralog's sequence from every phased haplotype assembly, builds a multiple sequence alignment (MSA) of the shared intronic region, infers a maximum-likelihood phylogenetic tree (IQ-TREE2), and assigns every haplotype copy to a clade/paralog group.
2. **`paralog_variants`** — takes the clade assignments and tree from `phylo_tree`, regroups sequences per clade, calls variants against the clade consensus, annotates them with Ensembl VEP, and summarizes copy-number and fixation status per gene family.

**`paralog_variants` consumes the outputs of `phylo_tree`** (merged fasta, clade table, tree file) — always run `phylo_tree` first for a given SD gene family before running `paralog_variants` on it.

## Repository layout

```
SD_gene_scripts/
├── phylo_tree/                          # Pipeline 1
│   ├── tree_method_only_intron_sge.smk  # Snakefile
│   ├── config.yaml
│   ├── env.cfg
│   └── manifest.txt                   # example sample manifest
├── paralog_variants/                    # Pipeline 2
│   ├── Snakefile
│   ├── config.yaml
│   ├── env.cfg
│   └── manifest.txt                     # example cluster manifest
└── scripts/                             # helper Python/R scripts shared by both pipelines
```

## Environment requirements

Both pipelines were built for an HPC environment using environment modules (see each `env.cfg`). Equivalent tools should be installed via conda/mamba if running elsewhere.

| Tool | phylo_tree | paralog_variants |
|---|---|---|
| Snakemake + Python 3 (pandas) | ✓ | ✓ |
| samtools | ✓ | |
| minimap2 | ✓ | |
| mafft | ✓ | ✓ |
| bedtools | ✓ | |
| trimAl | ✓ (used but **not listed** in `env.cfg` — make sure it's on `PATH`) | |
| IQ-TREE2 | ✓ | |
| rustybam (`rb`) | ✓ | |
| R (≥4.4, optparse + tree/plotting packages e.g. ape/ggtree) | ✓ | ✓ |
| Ensembl VEP (111) + AlphaMissense & PolyPhen_SIFT plugins | | ✓ |

```bash
# phylo_tree
module load miniconda/24.7.1 samtools/1.21 minimap2/2.28 mafft/7.525 bedtools/2.31.1 iqtree/2.1.2 R/4.4.1 rustybam/0.1.33

# paralog_variants
module load miniconda/24.7.1 mafft/7.525 R/4.4.1 ensembl-vep/111.0
```

## What you need to prepare

### 1. `phylo_tree`

Edit `phylo_tree/config.yaml`:

| Field | Description |
|---|---|
| `human_ref` | Reference FASTA (e.g. `T2T-CHM13v2.fasta`) |
| `tbl` | Path to your sample manifest (see below) |
| `gene_bed` | CHM13 gene-model BED used only to annotate the tree plots (format below) |
| `fragment_length` | Minimum contiguous fragment length (bp) used when defining the shared intronic region for the MSA |
| `human_gene` | Dict of `{cluster_id: {gene_symbol: "chrom:start-end"}}` — reference coordinates of every paralog/gene in each gene cluster you want to analyze |

Sample manifest (`tbl`, tab-separated, one row per haplotype):

| Column | Description |
|---|---|
| `sample` | Sample/haplotype ID (e.g. `HG00235_h1`) |
| `asm` | Path to that haplotype's phased assembly FASTA |
| `project` | Source project label (e.g. `HPRC-release2`) |
| `flagger` | Path to a Flagger BED file flagging unreliable/collapsed regions in that assembly |

`gene_bed` (e.g. `chm13v2.0_RefSeq_Liftoff_v5.2.gene.bed`) — tab-separated, **no header row**, one row per gene:

| Column # | Description |
|---|---|
| 1 | Chromosome |
| 2 | Gene start |
| 3 | Gene end |
| 4 | Gene symbol |

Example:
```
chr15	28419321	28443019	ARHGAP11B
chr15	30411628	30436433	ARHGAP11A
```

### 2. `paralog_variants`

Edit `paralog_variants/config.yaml`:

| Field | Description |
|---|---|
| `tbl` | Path to your cluster manifest (see below) |
| `sex` | `TRUE`/`FALSE` — whether to apply sex-chromosome-specific logic in the fixation step |
| `gene_model` | Path to the gene-model file used by `fixation.R` (format below) |
| `sample_info` | Path to the sample metadata file used by `fixation.R` (format below) |
| `gff_chm13` | GFF3 annotation passed to `vep --gff` for CHM13-referenced clusters |
| `gff_GRCh38` | GFF3 annotation passed to `vep --gff` for GRCh38-referenced clusters |
| `ref_chm13` | Reference FASTA passed to `vep --fasta` for CHM13-referenced clusters |
| `ref_GRCh38` | Reference FASTA passed to `vep --fasta` for GRCh38-referenced clusters |
| `alpha_missense_chm13` | AlphaMissense score file for the VEP `AlphaMissense` plugin, CHM13 branch |
| `alpha_missense_hg38` | AlphaMissense score file for the VEP `AlphaMissense` plugin, GRCh38 branch |
| `polyphen_chm13` | PolyPhen/SIFT score database for the VEP `PolyPhen_SIFT` plugin, CHM13 branch |

Cluster manifest (`tbl`, tab-separated, one row per gene cluster). Except for `cluster`, these are normally the **direct outputs of `phylo_tree`** for the same cluster:

| Column | Description | Typically comes from |
|---|---|---|
| `cluster` | Cluster ID, must match the `phylo_tree` cluster | your config |
| `paralog_table` | Clade/paralog assignment table | `phylo_tree`'s `results/{cluster}/tree/{cluster}_clade_information.txt` |
| `merged_fasta_file` | Merged multi-sample paralog fasta | `phylo_tree`'s `results/{cluster}/msa/sample.chm13.chimp.merged.fa` |
| `tree` | Phylogenetic tree file | `phylo_tree`'s `results/{cluster}/msa/msa.shared_region.trimal.fa.treefile` |

File formats for `gene_model` and `sample_info`:

- **`sample_info`** (e.g. `all_sample_info.txt`) — one row per sample, tab-separated, with a header row:

  | Column | Description |
  |---|---|
  | `sample` | Sample ID, must match the `sample` values in the `phylo_tree` manifest |
  | `sex` | `1` = Male, `2` = Female |

  Example:
  ```
  sample	sex
  HG00096	1
  HG00097	2
  HG00099	2
  ```

  `fixation.R` only classifies samples by **sex**; it does not use population/superpopulation information, so this file does not need those columns.

- **`gene_model`** (e.g. `chm13v2.0_RefSeq_Liftoff_v5.2_gene_model_exons_cds_utr.full.txt`) — used by `fixation.R` to draw the lollipop plots' gene/transcript track. Tab-separated, **no header row**, one row per exon/CDS/UTR feature:

  | Column # | Description |
  |---|---|
  | 1 | Chromosome |
  | 2 | Feature start (0-based) |
  | 3 | Feature end |
  | 4 | Feature type: `exon`, `CDS`, `5UTR`, or `3UTR` |
  | 5 | Gene symbol |
  | 6 | Transcript ID (RefSeq accession) |
  | 7 | Strand (`+`/`-`) |

  Example:
  ```
  chr15	28419321	28420142	exon	ARHGAP11B	NM_001039841.3	+
  chr15	28420014	28420142	CDS	ARHGAP11B	NM_001039841.3	+
  chr15	28419321	28420013	5UTR	ARHGAP11B	NM_001039841.3	+
  chr15	28429418	28429420	3UTR	ARHGAP11B	NM_001039841.3	+
  ```
  `fixation.R` filters this table to the gene matching each clade's `main_paralog` and to non-`exon` rows (i.e., `CDS`/`5UTR`/`3UTR`) to build the transcript model shown under the lollipop plot.

## How to run

```bash

# 1) phylo_tree
module load miniconda/24.7.1 samtools/1.21 minimap2/2.28 mafft/7.525 bedtools/2.31.1 iqtree/2.1.2 R/4.4.1 rustybam/0.1.33
snakemake -s phylo_tree/tree_method_only_intron_sge.smk \
          --configfile phylo_tree/config.yaml \
          --cores <N> -p

# 2) paralog_variants (after updating its manifest.txt with the phylo_tree outputs above)
module load miniconda/4.12.0 mafft/7.525 R/4.4.1 ensembl-vep/111.0
snakemake -s paralog_variants/Snakefile \
          --configfile paralog_variants/config.yaml \
          --cores <N> -p
```

For grid/cluster execution (SGE, etc.), submit via Snakemake's `--cluster`/`--profile` options — every rule already declares `resources: hrs=, mem=` for a cluster profile to consume.

## Outputs

### `phylo_tree` → `results/{cluster}/...`

| Directory/file | Contents |
|---|---|
| `gene_sequence/` | Extracted reference gene sequences, merged fasta, shared (intronic) region fasta/PAF/BED |
| `paralogs_bed/`, `get_sequence/` | Per-sample paralog coordinates and sequences, filtered by best reciprocal mapping and Flagger quality |
| `msa/` | Merged multi-sample paralog fasta, MSA (`msa.shared_region.fa`), trimmed MSA (`msa.shared_region.trimal.fa`), and the ML tree (IQ-TREE2 outputs: `.treefile`, `.contree`, `.iqtree`, `.mldist`, model/bootstrap files) |
| `tree/{cluster}_evolutionary_tree.pdf`, `_hprc_hgsvc_tree.pdf` | Annotated tree plots |
| `tree/{cluster}_clade_information.txt` | Per-haplotype clade/paralog assignment — **feeds into `paralog_variants`** |
| `tree/{cluster}_rooted.tree`, `_cn_counts.txt`, `_all_paralogs.txt` | Rooted tree, per-species copy-number counts, full paralog listing |

### `paralog_variants` → `results/{cluster}/variants/...`

| Directory/file | Contents |
|---|---|
| `{cluster}_{clade}.fa` / `.msa` | Per-clade paralog sequences and alignment |
| `{cluster}_{clade}.variants.txt`, `.vep_input.txt`, `.vep_output.txt` | Called variants and VEP functional annotation (incl. AlphaMissense / PolyPhen-SIFT scores) |
| `{cluster}_{clade}.snp_genotype.txt`, `.lollipop.txt`, `.detailed.lollipop.txt` | Per-sample genotype and variant lollipop summary tables |
| `{cluster}.cnv_plot.pdf`, `{cluster}.bar_plot.pdf` | Copy-number heatmap and variant summary bar plot |
| `{cluster}_final_variants_summary.txt` | Final fixation/variant summary table for monoclade groups in the gene family |


