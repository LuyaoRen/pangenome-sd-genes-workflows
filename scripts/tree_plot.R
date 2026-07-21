#!/usr/bin/env Rscript
library(RColorBrewer)
library(optparse)
library(ggtree)
library(ape)
library(stringr)
library(dplyr)
library(readr)
library(tidyr)
library(GenomicRanges)
library(ggplot2)
library(patchwork)

# ==== Parse command-line arguments ====
option_list = list(
  make_option(c("-t", "--tree"), type="character", help="Input tree file", metavar="character"),
  make_option(c("-p", "--paralog"), type="character", help="Paralog file", metavar="character"),
  make_option(c("-g", "--gene_bed"), type="character", help="Gene annotation BED file", metavar="character"),
  make_option(c("-f", "--flagger"), type="character", help="Flagger BED file", metavar="character"),
  make_option(c("-r", "--region"), type="character", help="Shared region file", metavar="character"),
  make_option("--p1_output", type="character", help="Output p1 pdf", metavar="character"),
  make_option("--p2_output", type="character", help="Output p2 pdf", metavar="character"),
  make_option("--summary", type="character", help="Output clade summary file", metavar="character"),
  make_option("--rooted_tree", type="character", help="Output rooted tree file", metavar="character"),
  make_option("--cn_count", type="character", help="Output copy number counts in primates", metavar="character"),
  make_option("--all_paralog", type="character", help="Output all paralogs in the tree", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# ==== Input files ====
tree_file <- opt$tree
paralog_file <- opt$paralog
gene_bed_file <- opt$gene_bed
flagger_file <- opt$flagger
region_file <- opt$region
p1_output_file <- opt$p1_output
p2_output_file <- opt$p2_output
clade_summary_file <- opt$summary
rooted_tree_file <- opt$rooted_tree
cn_count_file <- opt$cn_count
all_paralog_file <- opt$all_paralog

tree_file <- 'msa.shared_region.fa.treefile'
paralog_file <- 'filtered.merged.bed'
gene_bed_file <- 'chm13v2.0_RefSeq_Liftoff_v5.2.gene.bed'
flagger_file <- 'flagger.merged.bed'
region_file <- 'chm13.shared_region.bed'

# ==== Define function to get best overlap ====
#get_best_overlap <- function(query_gr, subject_gr, subject_gene_names) {
#  hits <- findOverlaps(query_gr, subject_gr)
#  overlap_df <- data.frame(
#    query = queryHits(hits),
#    subject = subjectHits(hits),
#    width = width(pintersect(query_gr[queryHits(hits)], subject_gr[subjectHits(hits)]))
#  )
#  best_hits <- overlap_df %>%
#    group_by(query) %>%
#    slice_max(order_by = width, n = 1, with_ties = FALSE) %>%
#    ungroup()
#  
#  result <- rep(NA_character_, length(query_gr))
#  result[best_hits$query] <- subject_gene_names[best_hits$subject]
#  return(result)
#}

get_best_overlap <- function(query_gr, subject_gr, subject_gene_names) {
  hits <- findOverlaps(query_gr, subject_gr)
  
  overlap_width <- width(pintersect(query_gr[queryHits(hits)], subject_gr[subjectHits(hits)]))
  subject_width <- width(subject_gr[subjectHits(hits)])
  overlap_ratio <- overlap_width / subject_width
  
  overlap_df <- data.frame(
    query = queryHits(hits),
    subject = subjectHits(hits),
    ratio = overlap_ratio
  )
  
  best_hits <- overlap_df %>%
    group_by(query) %>%
    slice_max(order_by = ratio, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  result <- rep(NA_character_, length(query_gr))
  result[best_hits$query] <- subject_gene_names[best_hits$subject]
  return(result)
}

parse_tip_labels <- function(tip_labels) {
  is_grch38 <- grepl("^GRCh38_h1_", tip_labels)
  is_mfa <- grepl("^MFA_h1_", tip_labels)
  
  tip_labels_non_grch38_mfa <- tip_labels[!(is_grch38 | is_mfa)]
  tip_info_non_grch38 <- str_match(tip_labels_non_grch38_mfa, "([^_]+)_([0-9]+)-([0-9]+)$")
  
  tip_df_non_grch38 <- data.frame(
    label = tip_labels_non_grch38_mfa,
    contig = tip_info_non_grch38[,2],
    start = as.integer(tip_info_non_grch38[,3]),
    end = as.integer(tip_info_non_grch38[,4]),
    stringsAsFactors = FALSE
  )
  
  tip_labels_grch38 <- tip_labels[is_grch38]
  tip_labels_grch38_clean <- sub("^GRCh38_h1_", "", tip_labels_grch38)
  tip_info_grch38 <- str_match(tip_labels_grch38_clean, "(.+)_([0-9]+)-([0-9]+)$")
  
  tip_df_grch38 <- data.frame(
    label = tip_labels_grch38,
    contig = tip_info_grch38[,2],
    start = as.integer(tip_info_grch38[,3]),
    end = as.integer(tip_info_grch38[,4]),
    stringsAsFactors = FALSE
  )
  
  tip_labels_mfa <- tip_labels[is_mfa]
  tip_labels_mfa_clean <- sub("^MFA_h1_", "", tip_labels_mfa)
  tip_info_mfa <- str_match(tip_labels_mfa_clean, "(.+)_([0-9]+)-([0-9]+)$")
  
  tip_df_mfa <- data.frame(
    label = tip_labels_mfa,
    contig = tip_info_mfa[,2],
    start = as.integer(tip_info_mfa[,3]),
    end = as.integer(tip_info_mfa[,4]),
    stringsAsFactors = FALSE
  )
  
  tip_df <- bind_rows(tip_df_non_grch38, tip_df_grch38, tip_df_mfa)
  return(tip_df)
}

#split_clades <- function(tree, node, ref_genes) {
#  tips_in_clade <- extract.clade(tree, node)$tip.label
#  ref_in_clade <- intersect(tips_in_clade, ref_genes)
#  
#  if (length(ref_in_clade) <= 1) {
#    return(list(node))
#  }
#  
#  ref_mrca <- MRCA(tree, ref_in_clade)
#  if (ref_mrca != node) {
#    return(list(node))
#  }
#  
#  child_nodes <- tree$edge[tree$edge[, 1] == node, 2]
#  result <- list()
#  for (child in child_nodes) {
#    result <- c(result, split_clades(tree, child, ref_genes))
#  }
#  return(result)
#}

split_clades <- function(tree, node, ref_genes) {

  if (node <= Ntip(tree)) return(list(node))  
  
  tips_in_clade <- extract.clade(tree, node)$tip.label
  ref_in_clade  <- intersect(tips_in_clade, ref_genes)
  
  if (length(ref_in_clade) <= 1) return(list(node))
  
  ref_mrca <- getMRCA(tree, ref_in_clade)    # 用 ape::getMRCA
  if (is.null(ref_mrca) || ref_mrca != node) return(list(node))
  
  child_nodes <- tree$edge[tree$edge[,1] == node, 2]
  
  result <- list()
  
  for (child in child_nodes[child_nodes > Ntip(tree)]) {
    result <- c(result, split_clades(tree, child, ref_genes))
  }
  return(result)
}


get_gene_colors <- function(genes, palette = "Set1") {
  n <- length(genes)
  max_n <- brewer.pal.info[palette, "maxcolors"]
  if (n <= max_n) {
    cols <- brewer.pal(n, palette)
  } else {
    cols <- colorRampPalette(brewer.pal(max_n, palette))(n)
  }
  setNames(cols, genes)
}


# ==== Load tree ====
tree <- read.tree(tree_file)
tip_labels <- tree$tip.label


# ==== Extract contig/start/end from tip labels ====
tip_df <- parse_tip_labels(tip_labels)
colnames(tip_df) <- c("label", "contig", "start", "end")
#tip_df$label <- tip_labels
tip_df <- tip_df %>% mutate(start = as.integer(start), end = as.integer(end))

# ==== Load predicted paralogs file ====
paralog_df <- read_tsv(paralog_file, col_names = c("chrom", "start", "end", "matched_region", "ref_len", "identity")) %>%
  mutate(start = as.integer(start), end = as.integer(end), chrom = sub(".*#", "", chrom))

merged_df <- left_join(tip_df, paralog_df, by = c("contig" = "chrom", "start", "end"))

# ==== Load gene annotations and replace '_' with '-' ====
gene_bed <- read_tsv(gene_bed_file, col_names = c("chr", "start", "end", "gene")) %>%
  mutate(gene = gsub("_", "-", gene))

# ==== matched_region -> predicted gene name ====
matched_df <- as.data.frame(str_match(merged_df$matched_region, "(chr[0-9XYM]+):([0-9]+)-([0-9]+)"))
colnames(matched_df) <- c("full", "chr", "start", "end")
matched_df$predicted_gene <- NA
valid_rows <- !is.na(matched_df$chr)
matched_gr <- GRanges(seqnames = matched_df$chr[valid_rows], ranges = IRanges(start = as.numeric(matched_df$start[valid_rows]), end = as.numeric(matched_df$end[valid_rows])))
gene_gr <- GRanges(seqnames = gene_bed$chr, ranges = IRanges(start = gene_bed$start, end = gene_bed$end), gene = gene_bed$gene)
matched_df$predicted_gene <- get_best_overlap(matched_gr, gene_gr, gene_gr$gene)
merged_df$predicted_gene <- matched_df$predicted_gene

# ==== Find actual self-mapped gene for CHM13 tips ====
tip_gr <- GRanges(seqnames = merged_df$contig, ranges = IRanges(start = merged_df$start, end = merged_df$end))
merged_df$self_gene <- get_best_overlap(tip_gr, gene_gr, gene_gr$gene)

# ==== Update tip labels ====
merged_df$final_gene <- ifelse(
  grepl("CHM13_h1_", merged_df$label) & !is.na(merged_df$predicted_gene) & !is.na(merged_df$self_gene),
  paste0(merged_df$predicted_gene, "(best_match)_", merged_df$self_gene, "(refseq)"),
  ifelse(!is.na(merged_df$predicted_gene), merged_df$predicted_gene, merged_df$matched_region)
)

merged_df$new_label <- ifelse(!is.na(merged_df$final_gene), paste0(merged_df$label, "_", merged_df$final_gene), merged_df$label)

label_map <- setNames(merged_df$new_label, merged_df$label)
tree$tip.label <- ifelse(tree$tip.label %in% names(label_map), label_map[tree$tip.label], tree$tip.label)
leaf_names <- tree$tip.label

# ==== Evolutionary tree (primate-only) ====
to_drop <- leaf_names[
  sapply(strsplit(leaf_names, "_"), function(x) {
    prefix <- x[1]
    any(grepl("^NA|^HG|^GRCh38", prefix))
  })
]
evolutionary_tree <- drop.tip(tree, to_drop)
evolutionary_tree_leaf_names <- evolutionary_tree$tip.label

outgroup <- NA
for (prefix in c("MFA", "SSY", "PAB", "PPY", "GGO", "PTR", "PPA")) {
  match <- evolutionary_tree_leaf_names[grepl(prefix, evolutionary_tree_leaf_names)]
  if (length(match) > 0) {
    outgroup <- match[1]
    break
  }
}

evolutionary_tree_rooted <- root(evolutionary_tree, which(evolutionary_tree$tip.label == outgroup))
evolutionary_tree_tip_labels <- evolutionary_tree_rooted$tip.label
species <- str_extract(evolutionary_tree_tip_labels, "^[^_]+")
species_counts <- data.frame(species = species) %>%
  count(species, name = "number") %>%
  arrange(desc(number))
write.table(species_counts,cn_count_file,sep="\t",col.names = F,row.names = F,quote = F)

p1 <- ggtree(evolutionary_tree_rooted) + geom_treescale() + geom_tiplab()
max_x <- max(p1$data$x, na.rm = TRUE)
p1 <- p1 + xlim(0, max_x * 2)

# ==== Sample tree (non-primate + CHM13 only) ====
ref_genes <- leaf_names[grepl("CHM13", leaf_names)]
to_keep <- leaf_names[
  sapply(strsplit(leaf_names, "_"), function(x) {
    prefix <- x[1]
    any(grepl("^NA|^HG|^GRCh38|^CHM13", prefix))
  })
]
to_keep <- c(to_keep, outgroup)

sample_tree <- keep.tip(tree, to_keep)
sample_tree_rooted <- root(sample_tree, which(sample_tree$tip.label == outgroup))
sample_tree_rooted <- drop.tip(sample_tree_rooted, outgroup)

write.tree(sample_tree_rooted, file = rooted_tree_file)
# ==== Shared region summary ====

region_df <- read_tsv(region_file, col_names = FALSE)

colnames(region_df)[1:3] <- c("matched_region", "start", "end")

region_length_summary <- region_df %>%
  mutate(region_length = end - start) %>%   # 每一行小区间长度
  group_by(matched_region) %>%
  summarise(total_region_length = sum(region_length, na.rm = TRUE)) %>%
  ungroup()

region_length <- round(region_length_summary$total_region_length[1] / 1000, 1)

sample_gene_df <- merged_df %>% filter(new_label %in% sample_tree_rooted$tip.label)
p2 <- ggtree(sample_tree_rooted)
p2$data <- left_join(p2$data, sample_gene_df[, c("new_label", "predicted_gene")], by = c("label" = "new_label"))


# ==== MRCA and clade coloring ====
mrca_root <- MRCA(sample_tree_rooted, ref_genes)
clade_nodes <- split_clades(sample_tree_rooted, mrca_root, ref_genes)
children_nodes <- unlist(clade_nodes)
clade_list <- list()
plot_layers <- list()
bootstrap_layers <- list()
colors <- RColorBrewer::brewer.pal(length(clade_nodes), "Paired")

for (i in seq_along(clade_nodes)) {
  node <- clade_nodes[[i]]
  tips <- extract.clade(sample_tree_rooted, node)$tip.label
  clade_list[[paste0("clade", i)]] <- tips
  plot_layers[[i]] <- geom_hilight(node=node, fill=colors[i], alpha=0.4)
  
  bs_val <- sample_tree_rooted$node.label[node - length(sample_tree_rooted$tip.label)]
  if (!is.na(bs_val)) {
    x <- max(p2$data$x[p2$data$node == node], na.rm = TRUE)
    y <- mean(p2$data$y[p2$data$node == node], na.rm = TRUE)
    bootstrap_layers[[i]] <- annotate("label", x = x, y = y, label = bs_val, size = 3, fill = "white")
  }
}

# ==== Prepare data for coloring tip points ====
for (layer in plot_layers) p2 <- p2 + layer

# ==== Assign gene color and annotate highlights ====
gene_levels <- sort(unique(na.omit(p2$data$predicted_gene)))
gene_colors <- get_gene_colors(gene_levels, "Set1")

highlight_tips <- leaf_names[grepl("CHM13_h1_", leaf_names) | grepl("GRCh38_h1_", leaf_names)]
p2 <- p2 +
  geom_tippoint(aes(color = predicted_gene), size = 2.5, na.rm = TRUE) +
  geom_tiplab(aes(subset = label %in% highlight_tips), size = 2.5, fontface = "bold", color = "red") +
  scale_color_manual(values = gene_colors, na.translate = FALSE, name = "Best Match") +
  ggtitle(paste0("Outgroup: ", outgroup, " | Regions for tree: ", region_length, " kb")) +
  theme_tree2()

# ==== Output clade and quality info ====
# Annotate clade summaries
flagger_df <- read_tsv(flagger_file, col_names = FALSE)
colnames(flagger_df)[1:7] <- c("raw_contig", "start", "end", "target_region", "ref_len", "identity", "score")
flagger_df <- flagger_df %>% mutate(contig = sub(".*#", "", raw_contig))

clade_output <- bind_rows(lapply(seq_along(clade_nodes), function(i) {
  tips <- extract.clade(sample_tree_rooted, children_nodes[i])$tip.label
  df <- merged_df %>% filter(new_label %in% tips) %>% select(label, contig, start, end, predicted_gene, ref_len, identity)
  mutate(df, clade = paste0("clade", i)) %>% select(clade, everything())
}))


qc_merged <- left_join(clade_output, flagger_df %>% select(contig, start, end, score), by = c("contig", "start", "end"))

all_paralog <- merged_df[
  -grep("MFA|SSY|PAB|PPY|GGO|PTR|PPA", merged_df$label),
  c('label','contig','start','end','predicted_gene','ref_len','identity')
]
all_paralog <- left_join(all_paralog, flagger_df %>% select(contig, start, end, score), by = c("contig", "start", "end"))

write_tsv(all_paralog, all_paralog_file, col_names = FALSE)

qc_summary <- qc_merged %>%
  group_by(clade) %>%
  summarise(total = n(), failed = sum(score != 0, na.rm = TRUE), missing = sum(is.na(score)), passed = total - failed - missing) %>%
  mutate(text_label = paste0("total: ", total, "\npass_flagger: ", passed, "\nfailed_flagger: ", failed, "\nnot_available: ", missing))

for (i in seq_along(children_nodes)) {
  node_id <- children_nodes[i]
  tips <- extract.clade(sample_tree_rooted, node_id)$tip.label
  x_pos <- mean(p2$data$x[p2$data$label %in% tips], na.rm = TRUE)
  y_pos <- mean(p2$data$y[p2$data$label %in% tips], na.rm = TRUE)
  label <- qc_summary$text_label[qc_summary$clade == paste0("clade", i)]
  
  # ==== Add clade label ====
  p2 <- p2 + annotate("text", x = x_pos + 0.01, y = y_pos, label = label, hjust = 0, size = 3)
  
}

# ==== Show bootstrap on children nodes ====
bootstrap_layers <- list()

for (child in children_nodes) {
  parent_id <- p2$data$parent[p2$data$node == child] 
  bs_val <- p2$data$label[p2$data$node == parent_id]
  
  if (!is.na(bs_val) && bs_val != "") {
    bs_score <- round(as.numeric(bs_val), 1)
    x <- p2$data$x[p2$data$node == parent_id]
    y <- p2$data$y[p2$data$node == parent_id]
    
    bootstrap_layers[[as.character(parent_id)]] <- 
      annotate("text", x = x, y = y, 
               label = bs_score, 
               size = 5, color = "black", hjust = -0.2)
  }
}

for (layer in bootstrap_layers) {
  p2 <- p2 + layer
}


failed_labels <- qc_merged %>% filter(score > 0) %>% pull(label)

matched_failed_data <- p2$data %>%
  filter(sapply(label, function(x) any(startsWith(x, failed_labels))))

p2 <- p2 +
  geom_point(
    data = matched_failed_data,
    aes(x = x, y = y),
    shape = 21,       
    fill = NA,       
    color = "black",  
    size = 2.5, 
    stroke = 0.8,
    show.legend = FALSE
  )


#check length
merged_df <- merged_df %>%
  separate(matched_region, into = c("chr", "coords"), sep = ":") %>%
  separate(coords, into = c("best_match_start", "best_match_end"), sep = "-")

merged_df$best_match_start <- as.numeric(merged_df$best_match_start)
merged_df$best_match_end <- as.numeric(merged_df$best_match_end)

merged_df$length_check <- ifelse(
  abs(merged_df$ref_len - (merged_df$best_match_end - merged_df$best_match_start)) /
    (merged_df$best_match_end - merged_df$best_match_start) < 0.4, 1, 0 )

failed_labels <- merged_df %>% filter(length_check == 0) %>% pull(label)

matched_failed_data <- p2$data %>%
  filter(sapply(label, function(x) any(startsWith(x, failed_labels))))

p2 <- p2 +
  geom_point(
    data = matched_failed_data,
    aes(x = x, y = y),
    shape = 4,       
    fill = NA,       
    color = "black",  
    size = 2, 
    stroke = 0.8,
    alpha = 0.5,
    show.legend = FALSE
  )

qc_merged$score[qc_merged$label %in% failed_labels] <- 2
write_tsv(qc_merged, clade_summary_file, col_names = FALSE)

# ==== Save plots ====
pdf(p1_output_file, height=15, width=15)
plot(p1)
dev.off()

pdf(p2_output_file, height=15, width=10)
plot(p2)
dev.off()




