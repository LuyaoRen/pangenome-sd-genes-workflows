# pangenome-sd-genes-workflows

Snakemake workflows for constructing phylogenetic trees from the pangenome (`phylo_tree`) and identifying and classifying variants for each paralog (`paralog_variants`).

## Overview

This repository contains two linked Snakemake pipelines for characterizing paralogous, protein-coding genes located within segmental duplications (SDs) across a panel of phased genome assemblies (HPRC/HGSVC):

1. **`phylo_tree`** — extracts each paralog's sequence from every phased haplotype assembly, builds a multiple sequence alignment (MSA) of the shared region, infers a maximum-likelihood phylogenetic tree (IQ-TREE2), and assigns every haplotype copy to a clade/paralog group.
2. **`paralog_variants`** — takes the clade assignments and tree from `phylo_tree`, regroups sequences per clade, calls variants against the clade consensus, annotates them with Ensembl VEP, and summarizes copy-number and fixation status per gene family.

**`paralog_variants` consumes the outputs of `phylo_tree`** (merged fasta, clade table, tree file) — always run `phylo_tree` first for a given gene cluster before running `paralog_variants` on it.

## Repository layout

```
SD_gene_scripts/
├── phylo_tree/                          # Pipeline 1
│   ├── tree_method_only_intron_sge.smk  # Snakefile
│   ├── config.yaml
│   ├── env.cfg
│   └── manifest.*.txt                   # example sample manifest
├── paralog_variants/                    # Pipeline 2
│   ├── Snakefile
│   ├── config.yaml
│   ├── env.cfg
│   └── manifest.txt                     # example cluster manifest
└── scripts/                             # helper Python/R scripts shared by both pipelines
```

Both Snakefiles call helper scripts using the relative path `scripts/...`, so **always launch `snakemake` from the `SD_gene_scripts/` repository root**, pointing at the Snakefile with `-s`.

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
module load miniconda/4.12.0 mafft/7.525 R/4.4.1 ensembl-vep/111.0
```

## What you need to prepare

### 1. `phylo_tree`

Edit `phylo_tree/config.yaml`:

| Field | Description |
|---|---|
| `human_ref` | Reference FASTA (e.g. `T2T-CHM13v2.fasta`), samtools-indexed |
| `tbl` | Path to your sample manifest (see below) |
| `fragment_length` | Minimum contiguous fragment length (bp) used when defining the shared intronic region for the MSA |
| `human_gene` | Dict of `{cluster_id: {gene_symbol: "chrom:start-end"}}` — reference coordinates of every paralog/gene in each gene cluster you want to analyze |

Sample manifest (`tbl`, tab-separated, one row per haplotype):

| Column | Description |
|---|---|
| `sample` | Sample/haplotype ID (e.g. `HG00235_h1`) |
| `asm` | Path to that haplotype's phased assembly FASTA (samtools/minimap2-indexable) |
| `project` | Source project label (e.g. `HPRC-release2`), for bookkeeping only |
| `flagger` | Path to a Flagger BED file flagging unreliable/collapsed regions in that assembly |

Also required (currently hardcoded in the Snakefile's `phylo_plot` rule — update the path for your system):
- `chm13v2.0_RefSeq_Liftoff_v5.2.gene.bed` — CHM13 gene-model BED, used to annotate the tree plots

### 2. `paralog_variants`

Edit `paralog_variants/config.yaml`:

| Field | Description |
|---|---|
| `tbl` | Path to your cluster manifest (see below) |
| `sex` | `TRUE`/`FALSE` — whether to apply sex-chromosome-specific logic in the fixation step |

Cluster manifest (`tbl`, tab-separated, one row per gene cluster). Except for `cluster`, these are normally the **direct outputs of `phylo_tree`** for the same cluster:

| Column | Description | Typically comes from |
|---|---|---|
| `cluster` | Cluster ID, must match the `phylo_tree` cluster | your config |
| `paralog_table` | Clade/paralog assignment table | `phylo_tree`'s `results/{cluster}/tree/{cluster}_clade_information.txt` |
| `merged_fasta_file` | Merged multi-sample paralog fasta | `phylo_tree`'s `results/{cluster}/msa/sample.chm13.chimp.merged.fa` |
| `tree` | Phylogenetic tree file | `phylo_tree`'s `results/{cluster}/msa/msa.shared_region.trimal.fa.treefile` |

Also required (currently hardcoded in the Snakefile and must be supplied/adjusted for your system):
- `all_sample_info.txt` — sample/haplotype metadata table (species, sex, population, etc.) used by `fixation.R`
- `chm13v2.0_RefSeq_Liftoff_v5.2_gene_model_exons_cds_utr.full.txt` — gene model used by `fixation.R`
- VEP reference/annotation set, auto-selected per cluster based on the first sequence ID in the MSA:
  - If CHM13-referenced: `chm13v2.0_RefSeq_Liftoff_v5.2.mane_select.gff3.gz`, `T2T-CHM13v2.fasta`, `AlphaMissense_chm13Lifted.header.tsv.gz`, `homo_sapiens_pangenome_PolyPhen_SIFT_20240502.db`
  - If GRCh38-referenced: `gencode.v49.annotation.gff3.mane_select.gz`, `hg38.no_alt.fa`, `AlphaMissense_hg38.tsv.gz`

## How to run

```bash
# from the SD_gene_scripts/ repository root

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
| `{cluster}_{clade}.snp_genotype.txt`, `.lollipop.txt`, `.detailed.lollipop.txt` | Per-sample genotype and variant "lollipop" tables |
| `{cluster}.cnv_plot.pdf`, `{cluster}.bar_plot.pdf` | Copy-number heatmap and variant summary bar plot |
| **`{cluster}_final_variants_summary.txt`** | Final per-gene-cluster fixation/variant summary table — top-level pipeline output |
