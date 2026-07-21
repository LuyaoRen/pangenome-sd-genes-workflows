#!/usr/bin/env Rscript
library(dplyr)
library(ape)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(purrr)
library(GenomicRanges)
library(trackViewer)
library(grid)
library(optparse)

option_list <- list(
  make_option(c("--calde_info"), type = "character", help = "Path to clade_information.txt"),
  make_option(c("--anno"), type = "character", help = "Path to all_sample_info.txt"),
  make_option(c("--detailed_lollipop"), type = "character", help = "Comma-separated lollipop detailed files"),
  make_option(c("--lollipop"), type = "character", help = "Comma-separated lollipop summary files"),
  make_option(c("--tree"), type = "character", help = "Path to phylogeny tree"),
  make_option(c("--gene_model"), type = "character", help = "Path to gene model annotation file"),
  make_option(c("--cnv_plot"), type = "character", help = "Path to output PDF for CNV heatmap"),
  make_option(c("--bar_plot"), type = "character", help = "Path to output PDF for variant barplot"),
  make_option(c("--fixation_table"), type = "character", help = "Path to output CSV fixation summary"),
  make_option(c("--sex"), action = "store_true", default = FALSE,help = "Whether to include sex-specific logic [default %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# cnvplot
df <- read.table(opt$calde_info,header=FALSE)
#df <- read.table('~/cluster_1_clade_information.txt.m',header=FALSE)
colnames(df) <- c('clade','tip_label','contig','start','end','paralog','length','sequence_identity','flagger_overlap')
df <- df %>%
  separate(tip_label, into = c("sample", "hap"), sep = "_", extra = "drop", remove = FALSE)
df$flagger <- ifelse(df$flagger_overlap > 0, 'fail', 'pass')
df$flagger[is.na(df$flagger)] <- "pass"
df <- df[!grepl('chr', df$contig), ]
#anno <- read.table(opt$anno,header=TRUE)
anno <- read.table('all_sample_info.txt',header=TRUE)
df <- merge(df, anno, by='sample', all=TRUE)

clade_main_paralog <- df %>%
  filter(!(is.na(clade) & is.na(paralog))) %>%
  group_by(clade, paralog) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(clade) %>%
  slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
  dplyr::select(clade, main_paralog = paralog)

df <- df %>%
  left_join(clade_main_paralog, by = "clade") %>%
  mutate(paralog = main_paralog) %>%
  dplyr::select(-main_paralog)

# Create full combination of sample x paralog x hap
all_combos <- expand_grid(
  sample = unique(na.omit(df$sample)),
  paralog = unique(na.omit(df$paralog)),
  hap = unique(na.omit(df$hap))
)

# Count presence per sample/paralog/hap
count_df <- df %>%
  filter(!is.na(sample) & !is.na(paralog) & !is.na(hap)) %>%
  group_by(sample, paralog, hap) %>%
  summarise(num = n(), .groups = "drop")

# Prepare main plot dataframe
plot_df <- all_combos %>%
  left_join(count_df, by = c("sample", "paralog", "hap")) %>%
  mutate(num = replace_na(num, 0)) %>%
  left_join(distinct(df[, c("sample", "sex")]), by = "sample")



# Add sex label and sorting order
plot_df$sex_label <- factor(ifelse(plot_df$sex == 1, "Male", "Female"), levels = c("Male", "Female"))
plot_df$order <- factor(plot_df$sample, levels = plot_df %>%
                          distinct(sample, sex_label) %>%
                          arrange(sex_label, sample) %>%
                          pull(sample))


if (opt$sex) {
  plot_df <- plot_df %>%
    group_by(paralog, sample) %>%
    mutate(num = if (any(sex == 1 & hap == "h1" & num > 0 & any(hap == "h2" & num == 0))) {
      ifelse(hap == "h1", num[hap == "h2"], ifelse(hap == "h2", num[hap == "h1"], num))
    } else {
      num
    }) %>%
    ungroup()
}

paralog_order <- names(sort(table(plot_df$paralog[plot_df$num > 0]), decreasing = TRUE))
plot_df$paralog <- factor(plot_df$paralog, levels = rev(paralog_order))
plot_df$num <- factor(plot_df$num)

# Prepare separate bar dfs
bar_df_sex <- plot_df %>% dplyr::select(order, sex_label)
bar_df_sex$group <- "Sex"
colnames(bar_df_sex)[2] <- "value"

bar_df <- bar_df_sex
bar_df$group <- factor(bar_df$group, levels = c("Sex"))

# Define colors
sex_colors <- c("Male" = "#4ABDAC", "Female" = "#FC4A1A")
bar_colors <- sex_colors

# h1: CN heatmap
h1 <- ggplot(plot_df, aes(order, paralog, fill = num)) + 
  geom_tile() +
  scale_fill_manual(name = "Copy Number", values = c("0"='white', "1" = "black", "2" = "#9E0142", "3" = "#D53E4F", "4" = "#FDAE61","5" = "#FFFFBF","6" = "#ABDDA4","7" = "#66C2A5","8" = "#3288BD","9" = "#5E4FA2")) +
  facet_grid(hap ~ ., scales = "free_y") +
  theme_minimal(base_size = 10) +
  theme(
    strip.text.y = element_text(size = 15, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    panel.spacing.x = unit(1, "mm"),
    legend.position = "none",
    plot.margin = margin(t = 0, r = 0, b = 15, l = 0)
  )

# h2: bar plots
h2 <- ggplot(bar_df, aes(x = order, fill = value)) +
  geom_bar(stat = "count") +
  facet_wrap(~ group, ncol = 1, strip.position = "left") +
  scale_fill_manual(name = "", values = bar_colors) +
  theme_void() +
  theme(strip.text = element_blank(), legend.position = "none")

# h3: CN count bar
cn_num <- plot_df %>%
  group_by(paralog, num, hap) %>%
  summarise(freq = n(), .groups = "drop")
cn_num$paralog <- factor(cn_num$paralog, levels = rev(paralog_order))

h3 <- ggplot(cn_num, aes(x = freq, y = paralog, fill = num)) +
  geom_bar(stat = "identity") +
  facet_grid(hap ~ .) +
  scale_x_reverse() +
  theme_half_open() +
  scale_fill_manual(
    name = "Copy Number",
    values = c("0"='white', "1" = "black", "2" = "#9E0142", "3" = "#D53E4F", "4" = "#FDAE61","5" = "#FFFFBF","6" = "#ABDDA4","7" = "#66C2A5","8" = "#3288BD","9" = "#5E4FA2"),
    guide = guide_legend(direction = "vertical")
  ) +
  theme(
    axis.line.y = element_blank(),
    strip.text = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 15, r = 0, b = 0, l = 0)
  )

legend_copy <- get_legend(
  ggplot(plot_df, aes(x = order, fill = num)) +
    geom_bar() +
    scale_fill_manual(name = "Copy Number", values = c("0"='white', "1" = "black", "2" = "#9E0142", "3" = "#D53E4F", "4" = "#FDAE61","5" = "#FFFFBF","6" = "#ABDDA4","7" = "#66C2A5","8" = "#3288BD","9" = "#5E4FA2")) +
    theme(legend.position = "right")
)

legend_sex <- get_legend(
  ggplot(bar_df_sex, aes(x = order, fill = value)) +
    geom_bar() +
    scale_fill_manual(name = "Sex", values = sex_colors) +
    theme(legend.position = "right")
)

plot_main <- plot_grid(h2, h1, ncol = 1, align = 'v', rel_heights = c(1, 10), axis = 'lr')
plot_combined <- plot_grid(h3, plot_main, ncol = 2, rel_widths = c(3, 10))
legend_all <- plot_grid(legend_sex, legend_copy, ncol = 1, align = 'v')
plot_final <- plot_grid(plot_combined, legend_all, ncol = 2, rel_widths = c(10, 2))

pdf(opt$cnv_plot,height = 4*length(table(cn_num$paralog))/2,width = 18)
#pdf('cnv_plot.pdf',height = 4*length(table(cn_num$paralog))/2,width = 18)
plot(plot_final)
dev.off()

# barplot
file_list <- unlist(strsplit(opt$detailed_lollipop, ","))
file_list <- c('cluster_245_clade1.detailed.lollipop.txt','cluster_245_clade2.detailed.lollipop.txt',
               'cluster_245_clade3.detailed.lollipop.txt','cluster_245_clade4.detailed.lollipop.txt',
               'cluster_245_clade5.detailed.lollipop.txt','cluster_245_clade6.detailed.lollipop.txt')

consequence_table <- do.call(rbind, lapply(file_list, function(file) {
  read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "",comment.char = "",fill = TRUE)
}))

consequence_summary <- consequence_table[-grep('-',consequence_table$tag),]
clean_string <- function(x) {
  x <- trimws(x)                      
  x <- gsub("[[:space:]]+", "", x)    
  x <- iconv(x, to = "UTF-8") 
  return(x)
}

df$tip_label <- clean_string(df$tip_label)
consequence_summary$hap <- clean_string(consequence_summary$hap)

consequence_summary <- consequence_summary %>%
  group_by(hap) %>%
  summarise(annotation = paste(tag, collapse = ","), .groups = "drop")


idx <- match(df$tip_label, consequence_summary$hap)
matched <- which(!is.na(idx))
df$tag <- 'no_mutation'
df$tag[matched] <- consequence_summary$annotation[idx[matched]]

df <- df %>%
  mutate(summary_tag = case_when(
    flagger == "fail" ~ "fail",
    flagger == "pass" & grepl("\\bdeleterious\\b", tag) ~ "deleterious",
    flagger == "pass" &
      grepl("\\bpredicted_deleterious\\b", tag) &
      !grepl("\\bdeleterious\\b", tag) ~ "predicted_deleterious",
    flagger == "pass" & tag == "no_mutation" ~ "pass",
    TRUE ~ NA_character_
  ))

copy_number <- as.data.frame(table(df[,c('paralog','sample','hap','summary_tag')]))
copy_number$summary_tag <- factor(copy_number$summary_tag, levels = c("fail", "pass","predicted_deleterious","deleterious"))
copy_number$sex <- anno$sex[match(copy_number$sample,anno$sample)]
copy_number$sex_label <- factor(ifelse(copy_number$sex == 1, "Male", "Female"), levels = c("Male", "Female"))
if (opt$sex) {
  copy_number <- copy_number %>%
    group_by(paralog, sample) %>%
    group_modify(~ {
      df_group <- .x
      has_male_h1_pos <- any(df_group$sex == 1 & df_group$hap == "h1" & df_group$Freq > 0)
      has_h2_zero <- any(df_group$hap == "h2" & df_group$Freq == 0)
      if (has_male_h1_pos && has_h2_zero && all(c("h1", "h2") %in% df_group$hap)) {
        idx_h1 <- which(df_group$hap == "h1")
        idx_h2 <- which(df_group$hap == "h2")
        temp_freq <- df_group$Freq[idx_h1]
        df_group$Freq[idx_h1] <- df_group$Freq[idx_h2]
        df_group$Freq[idx_h2] <- temp_freq
        temp_tag <- df_group$summary_tag[idx_h1]
        df_group$summary_tag[idx_h1] <- df_group$summary_tag[idx_h2]
        df_group$summary_tag[idx_h2] <- temp_tag
      }
      df_group
    }) %>%
    ungroup()
}
mycolor <- c(
  fail = "#9C9B97",
  pass = "#F1C787",
  predicted_deleterious = "#1098F7",
  deleterious = "#B21309"
)
p1 <- ggplot(data=copy_number, aes(x=sample, y=Freq,fill=summary_tag)) +
  geom_bar(stat="identity")+
  coord_flip()+
  scale_y_continuous(breaks = seq(min(copy_number$Freq), max(copy_number$Freq), by = 1)) + 
  theme_light()+
  theme(legend.position = "none")+
  scale_fill_manual(values=mycolor)+
  ylab('Copy Number')+
  xlab('')+
  facet_grid(sex_label ~ paralog + hap, space = "free", scales = "free")+
  theme(
    strip.text.x = element_text(size = 15, face = "italic", color = "black"), 
    strip.text.y = element_text(size = 15, face = "plain", color = "black"), 
    strip.background = element_rect(fill = "gray90"),
    axis.text.y = element_text(size = 5),   
    axis.text.x = element_text(size = 15), 
    axis.title.x = element_text(size = 15)
  )

pdf(opt$bar_plot,height = 13,width = 8*length(table(cn_num$paralog))/2)
plot(p1)
dev.off()


######### fixation table#########
##### copy number
homo_del <- df %>%
  filter(!is.na(clade)) %>%
  group_by(clade) %>%
  summarise('homo_del' = 298 - n_distinct(sample), .groups = "drop")

df <- df %>%
  group_by(sample, clade) %>%
  filter(!any(flagger == "fail")) %>%
  ungroup()

# 0. clade to gene

# 1. total sample n
n_sample_per_clade <- df %>%
  filter(!is.na(clade)) %>%
  group_by(clade) %>%
  summarise(n_sample = n_distinct(sample), .groups = "drop")

# 2. cn
df_wide <- df %>%
  filter(!is.na(sample) & !is.na(clade) & !is.na(hap) & !is.na(summary_tag)) %>%
  group_by(sample, clade, hap) %>%
  summarise(
    tag_num = n(),
    tag_info = paste(unique(summary_tag), collapse = ","),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = hap,
    values_from = c(tag_num, tag_info),
    names_glue = "{hap}_{.value}"
  ) %>%
  dplyr::rename(
    h1_num = h1_tag_num,
    h2_num = h2_tag_num,
    h1_info = h1_tag_info,
    h2_info = h2_tag_info
  )
df_wide$sex <- anno$sex[match(df_wide$sample,anno$sample)]
df_wide$sex_label <- factor(ifelse(df_wide$sex == 1, "Male", "Female"), levels = c("Male", "Female"))

if (opt$sex == TRUE) {
  df_wide <- df_wide %>%
    mutate(
      h1_info = as.character(h1_info),
      h2_info = as.character(h2_info)
    )
  swap_flag <- with(df_wide, sex == 1 & !is.na(h1_num) & is.na(h2_num))
  h1_num_old  <- df_wide$h1_num
  h2_num_old  <- df_wide$h2_num
  h1_info_old <- as.character(df_wide$h1_info)
  h2_info_old <- as.character(df_wide$h2_info)
  df_wide <- df_wide %>%
    mutate(
      h1_num  = ifelse(swap_flag, h2_num_old, h1_num_old),
      h2_num  = ifelse(swap_flag, h1_num_old, h2_num_old),
      h1_info = ifelse(swap_flag, h2_info_old, h1_info_old),
      h2_info = ifelse(swap_flag, h1_info_old, h2_info_old)
    )

  all_categories <- c(
    "Female_1+1", "Female_1+2", "Female_1+n", "Female_n+n", 
    "Female_0+1", "Female_0+2", "Female_0+n",
    "Male_0+na", "Male_1+na", "Male_2+na", "Male_n+na"
  )
  cn_summary <- df_wide %>%
    mutate(
      h1_num = replace_na(h1_num, 0),
      h2_num = replace_na(h2_num, 0),
      sex_label = ifelse(sex == 1, "Male", "Female")
    ) %>%
    mutate(
      cn_category = case_when(
        sex == 2 ~ case_when(
          h1_num == 1 & h2_num == 1 ~ "Female_1+1",
          h1_num == 1 & h2_num == 2 ~ "Female_1+2",
          h1_num == 2 & h2_num == 1 ~ "Female_1+2",
          h1_num == 1 & h2_num > 2  ~ "Female_1+n",
          h1_num > 2  & h2_num == 1 ~ "Female_1+n",
          h1_num > 1 & h2_num > 1   ~ "Female_n+n",
          h1_num == 0 & h2_num == 1 ~ "Female_0+1",
          h1_num == 1 & h2_num == 0 ~ "Female_0+1",
          h1_num == 0 & h2_num == 2 ~ "Female_0+2",
          h1_num == 2 & h2_num == 0 ~ "Female_0+2",
          h1_num == 0 & h2_num > 2  ~ "Female_0+n",
          h1_num > 2 & h2_num == 0  ~ "Female_0+n",
          TRUE ~ NA_character_
        ),
        sex == 1 ~ case_when( 
          h2_num == 0 ~ "Male_0+na",
          h2_num == 1 ~ "Male_1+na",
          h2_num == 2 ~ "Male_2+na",
          h2_num > 2  ~ "Male_n+na",
          TRUE ~ NA_character_
        ),
        TRUE ~ NA_character_
      ),
      cn_category = factor(cn_category, levels = all_categories)
    ) %>%
    group_by(clade, cn_category) %>%
    summarise(n = n(), .groups = "drop") %>%
    complete(clade, cn_category = all_categories, fill = list(n = 0)) %>%
    pivot_wider(names_from = cn_category, values_from = n)
}else{
  all_categories <- c("1+1", "1+2", "1+n", "n+n", "0+1", "0+2", "0+n")
  cn_summary <- df_wide %>%
  mutate(
    h1_num = replace_na(h1_num, 0),
    h2_num = replace_na(h2_num, 0)
  ) %>%
  mutate(
    cn_category = case_when(
      h1_num == 1 & h2_num == 1 ~ "1+1",
      h1_num == 1 & h2_num == 2 ~ "1+2",
      h1_num == 2 & h2_num == 1 ~ "1+2",
      h1_num == 1 & h2_num > 2  ~ "1+n",
      h1_num > 2  & h2_num == 1  ~ "1+n",
      h1_num > 1 & h2_num > 1   ~ "n+n",
      h1_num == 0 & h2_num == 1 ~ "0+1",
      h1_num == 1 & h2_num == 0 ~ "0+1",
      h1_num == 0 & h2_num == 2 ~ "0+2",
      h1_num == 2 & h2_num == 0 ~ "0+2",
      h1_num == 0 & h2_num > 2  ~ "0+n",
      h1_num > 2 & h2_num == 0  ~ "0+n",
      TRUE ~ NA_character_
    ),
    cn_category = factor(cn_category, levels = all_categories)
  ) %>%
  group_by(clade, cn_category) %>%
  summarise(n = n(), .groups = "drop") %>%
  complete(clade, cn_category = all_categories, fill = list(n = 0)) %>%
  pivot_wider(names_from = cn_category, values_from = n)
}

# 3. function
df_wide <- df_wide %>%
  mutate(
    cn = case_when(
      is.na(h1_num) & is.na(h2_num) ~ "homo_del",
      is.na(h1_num) | is.na(h2_num) ~ "hap",
      TRUE ~ "dip"
    )
  )

if(opt$sex){
  all_classes <- c(
    "Female_dip_func", "Female_dip_non_func", "Female_hap_func", "Female_hap_non_func",
    "Female_lof_dip_func", "Female_lof_dip_non_func", "Female_lof_hap_func", "Female_lof_hap_non_func",
    "Male_hap_func", "Male_hap_non_func", "Male_lof_hap_func", "Male_lof_hap_non_func"
  )
  
  func_summary <- df_wide %>%
    mutate(
      h1_info = as.character(h1_info),
      h2_info = as.character(h2_info),
      
      Female_dip_func = sex == 2 & cn == "dip" & (grepl("pass", h1_info) | grepl("pass", h2_info)),
      Female_dip_non_func = sex == 2 & cn == "dip" & !(grepl("pass", h1_info) | grepl("pass", h2_info)),
      Female_hap_func = sex == 2 & cn == "hap" & (grepl("pass", h1_info) | grepl("pass", h2_info)),
      Female_hap_non_func = sex == 2 & cn == "hap" & !(grepl("pass", h1_info) | grepl("pass", h2_info)),
      
      Female_lof_dip_func = sex == 2 & cn == "dip" &
        (!grepl("\\bdeleterious\\b", h1_info) | !grepl("\\bdeleterious\\b", h2_info)),
      
      Female_lof_dip_non_func = sex == 2 & cn == "dip" &
        (grepl("\\bdeleterious\\b", h1_info) & grepl("\\bdeleterious\\b", h2_info)),
      
      Female_lof_hap_func = sex == 2 & cn == "hap" & {
        info <- ifelse(is.na(h1_info), h2_info, h1_info)
        !grepl("\\bdeleterious\\b", info)
      },
      
      Female_lof_hap_non_func = sex == 2 & cn == "hap" & {
        info <- ifelse(is.na(h1_info), h2_info, h1_info)
        grepl("\\bdeleterious\\b", info)
      },
      
      Male_hap_func = sex == 1 & cn == "hap" & (grepl("pass", h1_info) | grepl("pass", h2_info)),
      Male_hap_non_func = sex == 1 & cn == "hap" & !(grepl("pass", h1_info) | grepl("pass", h2_info)),
      
      Male_lof_hap_func = sex == 1 & cn == "hap" & {
        info <- ifelse(is.na(h1_info), h2_info, h1_info)
        !grepl("\\bdeleterious\\b", info)
      },
      
      Male_lof_hap_non_func = sex == 1 & cn == "hap" & {
        info <- ifelse(is.na(h1_info), h2_info, h1_info)
        grepl("\\bdeleterious\\b", info)
      }
    ) %>%
    select(clade, all_of(all_classes)) %>%
    pivot_longer(
      cols = all_of(all_classes),
      names_to = "class",
      values_to = "match"
    ) %>%
    filter(match) %>%
    count(clade, class) %>%
    complete(clade, class = all_classes, fill = list(n = 0)) %>%
    pivot_wider(
      names_from = class,
      values_from = n
    )
}else{
  all_classes <- c(
    "dip_func", "dip_non_func", "hap_func", "hap_non_func",
    "lof_dip_func", "lof_dip_non_func", "lof_hap_func", "lof_hap_non_func"
  )
  
  func_summary <- df_wide %>%
    mutate(
      dip_func = cn == "dip" & (grepl("pass", h1_info) | grepl("pass", h2_info)),
      dip_non_func = cn == "dip" & !(grepl("pass", h1_info) | grepl("pass", h2_info)),
      hap_func = cn == "hap" & (grepl("pass", h1_info) | grepl("pass", h2_info)),
      hap_non_func = cn == "hap" & !(grepl("pass", h1_info) | grepl("pass", h2_info)),
      
      lof_dip_func = cn == "dip" &
        (!grepl("\\bdeleterious\\b", h1_info) | !grepl("\\bdeleterious\\b", h2_info)),
      
      lof_dip_non_func = cn == "dip" &
        (grepl("\\bdeleterious\\b", h1_info) & grepl("\\bdeleterious\\b", h2_info)),
      
      lof_hap_func = cn == "hap" & {
        info <- ifelse(is.na(h1_info), h2_info, h1_info)
        !grepl("\\bdeleterious\\b", info)
      },
      
      lof_hap_non_func = cn == "hap" & {
        info <- ifelse(is.na(h1_info), h2_info, h1_info)
        grepl("\\bdeleterious\\b", info)
      }
    ) %>%
    select(clade, all_of(all_classes)) %>%
    pivot_longer(
      cols = all_of(all_classes),
      names_to = "class",
      values_to = "match"
    ) %>%
    filter(match) %>%
    count(clade, class) %>%
    complete(clade, class = all_classes, fill = list(n = 0)) %>%
    pivot_wider(
      names_from = class,
      values_from = n
    )
}

#variant summary
lollipop_list <- unlist(strsplit(opt$lollipop, ",")) 
lollipop_list <- c('cluster_43_clade1.lollipop.txt','cluster_43_clade2.lollipop.txt')

variant_list <- list()  

for (i in lollipop_list) {
  lollipop_df <- read.delim(i, header = TRUE)
  clade <- sub(".*_(clade\\d+)\\.lollipop\\.txt$", "\\1", basename(i))
  total_n <- nrow(lollipop_df)
  mis_n <- sum(lollipop_df$Consequence == 'missense_variant')
  lof_n <- sum(lollipop_df$tag == 'deleterious')  # 注意拼写
  syn_n <- total_n - lof_n - mis_n
  mis_deleterious_n <- sum(lollipop_df$Consequence == 'missense_variant' & lollipop_df$tag == 'predicted_deleterious')
  mis_begin_n <- sum(lollipop_df$Consequence == 'missense_variant' & lollipop_df$tag == '-')
  variant_list[[clade]] <- data.frame(
    clade = clade,
    total = total_n,
    syn = syn_n,
    mis = mis_n,
    mis_deleterious = mis_deleterious_n,
    mis_begin = mis_begin_n,
    lof = lof_n
  )
}

variant_summary <- bind_rows(variant_list)


hap_summary <- df_wide %>%
  select(clade, h1_info, h2_info) %>%
  pivot_longer(cols = c(h1_info, h2_info), names_to = "hap", values_to = "info") %>%
  mutate(
    is_mis_lof = info %in% c("predicted_deleterious", "deleterious"),
    is_lof = info == "deleterious"
  ) %>%
  group_by(clade) %>%
  summarise(hap_total = sum(!is.na(info)),
    non_func_hap_mis_lof = sum(is_mis_lof, na.rm = TRUE),
    non_func_hap_lof = sum(is_lof, na.rm = TRUE),
    .groups = "drop"
  )

#popgenome
#cluster.fasta <- readData(opt$fasta_folder)
#cluster.fasta <- readData('test_fa/')
#cluster.fasta = diversity.stats(cluster.fasta,pi = TRUE, keep.site.info = FALSE)
#cluster.fasta = neutrality.stats(cluster.fasta)
#pi <- cluster.fasta@nuc.diversity.within/cluster.fasta@n.sites
#tajima <- cluster.fasta@Tajima.D
#n_sites <- cluster.fasta@n.sites

#pop_genetic_summary <- data.frame(
#  clade = sapply(strsplit(cluster.fasta@region.names, "_"), `[`, 3),
#  pi = as.numeric(pi),
#  tajima_D = as.numeric(tajima),
#  CDS_length = as.numeric(n_sites)
#)

# genetic distance
#tree <- read.tree(opt$tree) 
tree <- read.tree('msa.shared_region.fa.treefile') 
#df <- read.table(opt$calde_info,header=F)
df <- read.table('cluster_43_clade_information.txt',header=F)
df <- df[df$V9 <= 0 | is.na(df$V9), ]
colnames(df)[1:2] <- c("clade", "tip")
dist_matrix <- cophenetic.phylo(tree)
pairwise_df <- df %>%
  select(clade, tip) %>%
  distinct() %>%
  group_by(clade) %>%
  summarise(pairs = list(combn(tip, 2, simplify = FALSE)), .groups = "drop") %>%
  unnest(pairs) %>%
  mutate(
    tip1 = sapply(pairs, `[`, 1),
    tip2 = sapply(pairs, `[`, 2),
    dist = mapply(function(a, b) dist_matrix[a, b], tip1, tip2)
  ) %>%
  select(clade, tip1, tip2, dist)

clade_stats <- pairwise_df %>%
  group_by(clade) %>%
  summarise(
    genetic_distance_mean = mean(dist, na.rm = TRUE),
    genetic_distance_sd = sd(dist, na.rm = TRUE)
  )



final_df <- purrr::reduce(
  list(clade_main_paralog,n_sample_per_clade,homo_del,cn_summary, func_summary, variant_summary, hap_summary,clade_stats),
  ~ merge(.x, .y, by = "clade", all = TRUE)
)

write.table(final_df, file = opt$fixation_table, row.names = FALSE,sep="\t",quote = F)

# lollipop
lollipop_list <- unlist(strsplit(opt$lollipop, ",")) 

for(i in lollipop_list){
  df <- read.delim(i, header = TRUE)[, 1:7]
  if (nrow(df) == 0) next  
  plot_name <- gsub("txt","pdf",i)
  clade <- sub(".*_(clade\\d+)\\.lollipop\\.txt$", "\\1", basename(i))
  gene_name <- clade_main_paralog %>%
    filter(clade == !!clade) %>%
    pull(main_paralog) %>%
    as.character()
  
  df$color <- '#FBD178'  
  df$color[df$Consequence == 'missense_variant'] <- '#79ADD6'
  df$color[df$Consequence %in% c('stop_gained', 'start_lost', 'stop_lost',
                                 'frameshift_variant', 'splice_donor_variant',
                                 'splice_acceptor_variant')] <- '#EB716B'
  
  df$ref <- sapply(strsplit(df$Amino_acids, "/"), `[`, 1)
  df$alter <- sapply(strsplit(df$Amino_acids, "/"), `[`, 2)
  df$alter[is.na(df$alter)] <- ''
  df$anno <- paste0(df$ref, df$Protein_position, df$alter)
  df$border <- "black"
  df$anno <- ifelse(df$tag == "predicted_deleterious", paste0(df$anno, "(LD)"), df$anno)
  # gene model
  gene_model_raw <- read.table('202503/chm13v2.0_RefSeq_Liftoff_v5.2_gene_model_exons_cds_utr.full.txt')
  #gene_model_raw <- read.table(opt$gene_model)
  gene_model <- gene_model_raw[gene_model_raw$V5 == gene_name & gene_model_raw$V4 != "exon", ]
  chromo <- unique(gene_model$V1)
  text_anno <- paste0(unique(gene_model$V5)[1], " ", unique(gene_model$V6)[1])
  gene_model$height <- ifelse(gene_model$V4 == "CDS", 0.03, 0.015)
  gene_model <- gene_model[, c(2, 3, 8)]
  colnames(gene_model) <- c("start", "end", "height")
  gene_model$color <- ifelse(gene_model$height == 0.03, "#C2DFEA", "#a2a9af")
  
  sample.gr <- GRanges(chromo, IRanges(df$POS, width = 1, names = df$anno))
  sample.gr$score <- as.numeric(df$num)
  sample.gr$color <- df$color
  sample.gr$border <- df$border
  sample.gr$label.parameter.rot <- 45
  sample.gr$label.parameter.label <- names(sample.gr)
  sample.gr$label.parameter.label[which(df$color=='#FBD178')] <- NA
  
  features <- GRanges(chromo, IRanges(gene_model$start, width = gene_model$end - gene_model$start))
  features$fill <- gene_model$color
  features$height <- gene_model$height
  
  legend <- list(labels = c('Synonymous', 'Missense', 'Nonsense'),
                 fill = c('#FBD178', '#79ADD6', '#EB716B'))
  
  pdf(plot_name, height = 8, width = 12)
  lolliplot(sample.gr, features, lollipop_style_switch_limit = 1,ylab.gp=gpar(cex=1.2),
            legend = legend, ylab = "#n Haplotype\n", xaxis = FALSE)
  grid.text(text_anno, x = 0.5, y = 0.01, just = "bottom",  gp=gpar(cex=1.1))
  dev.off()
}


