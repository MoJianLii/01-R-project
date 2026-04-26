# =========================================================
# neuron_interaction_topjournal_R_output_02.R
# 目的：围绕4条核心结论输出可直接用于论文的高质量可视化
# 路径：E:/zaw/2603
# 输出：E:/zaw/2603/R_output_02
# =========================================================

required_pkgs <- c(
  "data.table", "dplyr", "tidyr", "stringr", "purrr", "readr", "tibble",
  "forcats", "ggplot2", "ggrepel", "scales", "sandwich", "lmtest"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(readr)
  library(tibble)
  library(forcats)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(sandwich)
  library(lmtest)
})

# ---------------------------
# Paths
# ---------------------------
base_dir <- "E:/zaw/2603"
out_dir <- file.path(base_dir, "R_output_02")
fig_dir <- file.path(out_dir, "figures")
tbl_dir <- file.path(out_dir, "tables")
meta_dir <- file.path(out_dir, "metadata")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tbl_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(meta_dir, showWarnings = FALSE, recursive = TRUE)

find_input_csv <- function(base_dir) {
  preferred <- c(
    file.path(base_dir, "neuron-neuron-partner.csv"),
    file.path(base_dir, "neuron_neuron_partner.csv"),
    file.path(base_dir, "partner.csv")
  )
  preferred <- preferred[file.exists(preferred)]
  if (length(preferred) > 0) return(preferred[1])
  
  hits <- list.files(
    base_dir,
    pattern = "(partner|neuron).*(csv)$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(hits) > 0) return(hits[1])
  
  stop("未在E:/zaw/2603中找到输入CSV。请把原始数据放到该目录下，并命名为neuron-neuron-partner.csv或类似名称。")
}

input_csv <- find_input_csv(base_dir)

# ---------------------------
# Global styles
# ---------------------------
base_family <- "sans"

pair_cols <- c(
  "E-E" = "#3C78D8",
  "E-I" = "#D1495B",
  "I-I" = "#2A9D8F"
)

pref_cols <- c(
  "same-type preferring" = "#2A9D8F",
  "cross-type preferring" = "#D1495B"
)

bg_cols <- c(
  "Higher than block E-I background" = "#C44E52",
  "Lower than block E-I background" = "#4C72B0"
)

pub_theme <- function(base_size = 12) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0, margin = margin(b = 8)),
      plot.subtitle = element_text(size = base_size, hjust = 0, color = "#444444", margin = margin(b = 10)),
      plot.caption = element_text(size = base_size - 1, color = "#555555", hjust = 0),
      axis.title = element_text(face = "bold", size = base_size + 0.5),
      axis.text = element_text(color = "black", size = base_size),
      axis.line = element_line(linewidth = 0.55, color = "black"),
      axis.ticks = element_line(linewidth = 0.5, color = "black"),
      axis.ticks.length = grid::unit(-1.8, "mm"),
      axis.text.x = element_text(margin = margin(t = 7)),
      axis.text.y = element_text(margin = margin(r = 7)),
      strip.background = element_rect(fill = "#F2F2F2", color = NA),
      strip.text = element_text(face = "bold", size = base_size),
      panel.grid = element_blank(),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size),
      legend.key.width = grid::unit(1.2, "cm"),
      plot.margin = margin(10, 14, 10, 10)
    )
}

theme_set(pub_theme(12))

save_pub <- function(p, filename, width, height) {
  png(file.path(fig_dir, paste0(filename, ".png")), width = width, height = height, units = "in", res = 600, bg = "white")
  print(p)
  dev.off()
  
  pdf(file.path(fig_dir, paste0(filename, ".pdf")), width = width, height = height, bg = "white", useDingbats = FALSE)
  print(p)
  dev.off()
}

fmt_num <- function(x, digits = 4) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

fmt_p <- function(x) {
  ifelse(
    is.na(x), "NA",
    ifelse(x < 1e-3, formatC(x, format = "e", digits = 2), formatC(x, format = "f", digits = 4))
  )
}

canonical_pair <- function(a, b) {
  ifelse(a < b, paste(a, b, sep = " || "), paste(b, a, sep = " || "))
}

clean_subclass <- function(x) {
  str_replace(x, "^\\d+\\s+", "")
}

rank_biserial_paired <- function(x, y) {
  d <- x - y
  d <- d[d != 0]
  if (length(d) == 0) return(NA_real_)
  pos <- sum(d > 0)
  neg <- sum(d < 0)
  (pos - neg) / length(d)
}

safe_paired_wilcox <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 3) {
    return(list(p_value = NA_real_, n = length(x), mean_diff = NA_real_, median_diff = NA_real_, rank_biserial = NA_real_))
  }
  wt <- suppressWarnings(wilcox.test(x, y, paired = TRUE, exact = FALSE))
  list(
    p_value = wt$p.value,
    n = length(x),
    mean_diff = mean(x - y),
    median_diff = median(x - y),
    rank_biserial = rank_biserial_paired(x, y)
  )
}

paired_wilcox_df <- function(df, a, b, metric_name) {
  tmp <- df %>% select(all_of(c(a, b))) %>% drop_na()
  res <- safe_paired_wilcox(tmp[[a]], tmp[[b]])
  tibble(
    metric = metric_name,
    comparison = paste(a, "vs", b),
    n = res$n,
    mean_diff = res$mean_diff,
    median_diff = res$median_diff,
    rank_biserial = res$rank_biserial,
    p_value = res$p_value
  )
}

pick_focus_row <- function(df, tokens) {
  if (nrow(df) == 0) return(df[0, ])
  keep <- rep(TRUE, nrow(df))
  for (tok in tokens) {
    keep <- keep & str_detect(df$subpair_name_raw, fixed(tok, ignore_case = TRUE))
  }
  out <- df[keep, , drop = FALSE]
  if (nrow(out) == 0) return(df[0, ])
  out %>% slice(1)
}

make_empty_plot <- function(title_text) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = title_text, size = 5, fontface = "bold") +
    xlim(0, 1) + ylim(0, 1) +
    theme_void(base_family = base_family)
}

# ---------------------------
# Read data
# ---------------------------
df <- fread(input_csv)
df[, pair_key := canonical_pair(cluster.1_label, cluster.2_label)]

# ---------------------------
# Unique clusters
# ---------------------------
c1 <- copy(df[, .(
  label = cluster.1_label,
  slide = cluster.1_slide,
  layer = cluster.1_layer,
  region = cluster.1_region,
  cell_Neruon_type = cluster.1_cell_Neruon_type,
  subclass = cluster.1_subclass,
  total_cell_num = cluster.1_total_cell_num,
  Glut_Neruon_cell_ids_num = cluster.1_Glut_Neruon_cell_ids_num,
  GABA_Neruon_cell_ids_num = cluster.1_GABA_Neruon_cell_ids_num,
  cauchy_combination_p = cluster.1_cauchy_combination_p,
  E_I_Ratio = cluster.1_E_I_Ratio
)])

c2 <- copy(df[, .(
  label = cluster.2_label,
  slide = cluster.2_slide,
  layer = cluster.2_layer,
  region = cluster.2_region,
  cell_Neruon_type = cluster.2_cell_Neruon_type,
  subclass = cluster.2_subclass,
  total_cell_num = cluster.2_total_cell_num,
  Glut_Neruon_cell_ids_num = cluster.2_Glut_Neruon_cell_ids_num,
  GABA_Neruon_cell_ids_num = cluster.2_GABA_Neruon_cell_ids_num,
  cauchy_combination_p = cluster.2_cauchy_combination_p,
  E_I_Ratio = cluster.2_E_I_Ratio
)])

clusters <- unique(rbindlist(list(c1, c2), use.names = TRUE))
clusters[, block := paste(slide, layer, sep = "|")]
clusters[, subclass_clean := clean_subclass(subclass)]

# ---------------------------
# Unique undirected pairs
# ---------------------------
pair_list <- lapply(split(df, df$pair_key), function(sub) {
  stopifnot(nrow(sub) == 2)
  labs <- sort(unique(c(sub$cluster.1_label, sub$cluster.2_label)))
  a <- labs[1]
  b <- labs[2]
  ra <- sub[sub$cluster.1_label == a, ][1, ]
  type_a <- ifelse(ra$cluster.1_cell_Neruon_type == "Glut", "E", "I")
  type_b <- ifelse(ra$cluster.2_cell_Neruon_type == "Glut", "E", "I")
  
  sub_raw_sorted <- sort(c(ra$cluster.1_subclass, ra$cluster.2_subclass))
  sub_clean_sorted <- sort(c(clean_subclass(ra$cluster.1_subclass), clean_subclass(ra$cluster.2_subclass)))
  
  tibble(
    label_a = a,
    label_b = b,
    type_a_raw = ra$cluster.1_cell_Neruon_type,
    type_b_raw = ra$cluster.2_cell_Neruon_type,
    type_a = type_a,
    type_b = type_b,
    pair_type = paste(sort(c(type_a, type_b)), collapse = "-"),
    subclass_a = ra$cluster.1_subclass,
    subclass_b = ra$cluster.2_subclass,
    subclass_a_clean = clean_subclass(ra$cluster.1_subclass),
    subclass_b_clean = clean_subclass(ra$cluster.2_subclass),
    subpair_name_raw = paste(sub_raw_sorted, collapse = " × "),
    subpair_name_clean = paste(sub_clean_sorted, collapse = " × "),
    layer = ra$cluster.1_layer,
    slide = ra$cluster.1_slide,
    region_a = ra$cluster.1_region,
    region_b = ra$cluster.2_region,
    size_a = ra$cluster.1_total_cell_num,
    size_b = ra$cluster.2_total_cell_num,
    EI_ratio_a = ra$cluster.1_E_I_Ratio,
    EI_ratio_b = ra$cluster.2_E_I_Ratio,
    p_a = ra$cluster.1_cauchy_combination_p,
    p_b = ra$cluster.2_cauchy_combination_p,
    overlap = ra$overlap_cell,
    union = ra$union_cell,
    jaccard = ra$jaccard,
    overlap_prop_a = ra$cluster.1.overlap.percent,
    overlap_prop_b = ra$cluster.2.overlap.percent
  )
})

pairs_u <- bind_rows(pair_list) %>%
  mutate(
    block = paste(slide, layer, sep = "|"),
    size_small = pmin(size_a, size_b),
    size_large = pmax(size_a, size_b),
    size_ratio = size_large / size_small,
    overlap_prop_small = overlap / size_small,
    overlap_prop_large = overlap / size_large,
    nesting_gap = overlap_prop_small - overlap_prop_large,
    nested_90 = overlap_prop_small >= 0.90,
    pair_type = factor(pair_type, levels = c("E-E", "E-I", "I-I")),
    layer = factor(layer, levels = c("layer 1", "layer 2/3", "layer 4", "layer 5", "layer 6a", "layer 6b"))
  )

# ---------------------------
# Long edge table for cluster-centric preference
# ---------------------------
edges_long <- bind_rows(
  pairs_u %>% transmute(
    label = label_a,
    type_raw = type_a_raw,
    subclass = subclass_a,
    subclass_clean = subclass_a_clean,
    layer = layer,
    partner_label = label_b,
    partner_type = type_b_raw,
    jaccard = jaccard
  ),
  pairs_u %>% transmute(
    label = label_b,
    type_raw = type_b_raw,
    subclass = subclass_b,
    subclass_clean = subclass_b_clean,
    layer = layer,
    partner_label = label_a,
    partner_type = type_a_raw,
    jaccard = jaccard
  )
)

# ---------------------------
# Core summaries
# ---------------------------
block_metrics <- pairs_u %>%
  group_by(block, slide, layer, pair_type) %>%
  summarise(
    n_pairs = n(),
    mean_jaccard = mean(jaccard, na.rm = TRUE),
    mean_overlap_prop_small = mean(overlap_prop_small, na.rm = TRUE),
    mean_overlap_prop_large = mean(overlap_prop_large, na.rm = TRUE),
    mean_nesting_gap = mean(nesting_gap, na.rm = TRUE),
    prop_nested90 = mean(nested_90, na.rm = TRUE),
    .groups = "drop"
  )

blocks_complete <- block_metrics %>% count(block) %>% filter(n == 3) %>% pull(block)
block_metrics_complete <- block_metrics %>% filter(block %in% blocks_complete)

block_small_wide <- block_metrics_complete %>%
  select(block, slide, layer, pair_type, mean_overlap_prop_small) %>%
  pivot_wider(names_from = pair_type, values_from = mean_overlap_prop_small) %>%
  drop_na()

block_gap_wide <- block_metrics_complete %>%
  select(block, slide, layer, pair_type, mean_nesting_gap) %>%
  pivot_wider(names_from = pair_type, values_from = mean_nesting_gap) %>%
  drop_na()

block_j_wide <- block_metrics_complete %>%
  select(block, slide, layer, pair_type, mean_jaccard) %>%
  pivot_wider(names_from = pair_type, values_from = mean_jaccard) %>%
  drop_na()

block_nested_wide <- block_metrics_complete %>%
  select(block, slide, layer, pair_type, prop_nested90) %>%
  pivot_wider(names_from = pair_type, values_from = prop_nested90) %>%
  drop_na()

# ---------------------------
# Tests: geometry family
# 与原思路一致：9个检验一起做BH，便于和聊天中的FDR口径一致
# ---------------------------
overall_tests <- bind_rows(
  paired_wilcox_df(block_j_wide, "E-E", "E-I", "block_mean_jaccard"),
  paired_wilcox_df(block_j_wide, "E-E", "I-I", "block_mean_jaccard"),
  paired_wilcox_df(block_j_wide, "E-I", "I-I", "block_mean_jaccard"),
  paired_wilcox_df(block_small_wide, "E-E", "E-I", "block_mean_overlap_prop_small"),
  paired_wilcox_df(block_small_wide, "E-E", "I-I", "block_mean_overlap_prop_small"),
  paired_wilcox_df(block_small_wide, "E-I", "I-I", "block_mean_overlap_prop_small"),
  paired_wilcox_df(block_gap_wide, "E-E", "E-I", "block_mean_nesting_gap"),
  paired_wilcox_df(block_gap_wide, "E-E", "I-I", "block_mean_nesting_gap"),
  paired_wilcox_df(block_gap_wide, "E-I", "I-I", "block_mean_nesting_gap")
) %>% mutate(fdr_bh = p.adjust(p_value, method = "BH"))

# ---------------------------
# Tests: layer-wise near-complete nesting family
# ---------------------------
nested_layer_tests <- block_metrics_complete %>%
  group_by(layer) %>%
  group_split() %>%
  map_dfr(function(dat) {
    layer_name <- unique(dat$layer)
    piv <- dat %>%
      select(block, slide, pair_type, prop_nested90) %>%
      pivot_wider(names_from = pair_type, values_from = prop_nested90) %>%
      drop_na()
    
    bind_rows(
      paired_wilcox_df(piv, "E-E", "E-I", "prop_nested90"),
      paired_wilcox_df(piv, "E-E", "I-I", "prop_nested90"),
      paired_wilcox_df(piv, "E-I", "I-I", "prop_nested90")
    ) %>% mutate(layer = as.character(layer_name))
  }) %>%
  mutate(fdr_bh = p.adjust(p_value, method = "BH"))

# ---------------------------
# Cluster-centric preference
# ---------------------------
cluster_pref <- edges_long %>%
  mutate(same_type_partner = type_raw == partner_type) %>%
  group_by(label, type_raw, subclass, subclass_clean, layer, same_type_partner) %>%
  summarise(mean_jaccard = mean(jaccard, na.rm = TRUE), n = n(), .groups = "drop") %>%
  pivot_wider(names_from = same_type_partner, values_from = mean_jaccard, names_prefix = "same_") %>%
  rename(cross_type_mean_jaccard = same_FALSE, same_type_mean_jaccard = same_TRUE) %>%
  drop_na(cross_type_mean_jaccard, same_type_mean_jaccard) %>%
  mutate(same_minus_cross = same_type_mean_jaccard - cross_type_mean_jaccard)

exc_pref <- cluster_pref %>% filter(type_raw == "Glut")
gaba_pref <- cluster_pref %>% filter(type_raw == "Gaba")

safe_pref_test <- function(dat) {
  if (nrow(dat) < 3) {
    return(tibble(
      n_clusters = nrow(dat),
      mean_same_minus_cross = mean(dat$same_minus_cross, na.rm = TRUE),
      median_same_minus_cross = median(dat$same_minus_cross, na.rm = TRUE),
      p_value = NA_real_
    ))
  }
  tibble(
    n_clusters = nrow(dat),
    mean_same_minus_cross = mean(dat$same_minus_cross, na.rm = TRUE),
    median_same_minus_cross = median(dat$same_minus_cross, na.rm = TRUE),
    p_value = suppressWarnings(wilcox.test(dat$same_type_mean_jaccard, dat$cross_type_mean_jaccard, paired = TRUE, exact = FALSE)$p.value)
  )
}

exc_sub_tests <- exc_pref %>%
  group_by(subclass, subclass_clean) %>%
  group_modify(~ safe_pref_test(.x)) %>%
  ungroup() %>%
  filter(n_clusters >= 5) %>%
  mutate(fdr_bh = p.adjust(p_value, method = "BH"))

gaba_sub_tests <- gaba_pref %>%
  group_by(subclass, subclass_clean) %>%
  group_modify(~ safe_pref_test(.x)) %>%
  ungroup() %>%
  filter(n_clusters >= 5) %>%
  mutate(fdr_bh = p.adjust(p_value, method = "BH"))

exc_sub_layer_tests <- exc_pref %>%
  group_by(subclass, subclass_clean, layer) %>%
  group_modify(~ safe_pref_test(.x)) %>%
  ungroup() %>%
  filter(n_clusters >= 3) %>%
  mutate(fdr_bh = p.adjust(p_value, method = "BH"))

# ---------------------------
# E-I subclass pair vs within-block E-I background
# ---------------------------
subpair_tests <- pairs_u %>%
  group_by(subpair_name_raw, subpair_name_clean) %>%
  group_split() %>%
  map_dfr(function(sub) {
    sp <- unique(sub$subpair_name_raw)
    sp_clean <- unique(sub$subpair_name_clean)
    if (nrow(sub) < 5 || n_distinct(sub$block) < 3) return(NULL)
    
    sp_block <- sub %>%
      group_by(block, pair_type) %>%
      summarise(sp_mean = mean(jaccard, na.rm = TRUE), .groups = "drop")
    
    other_bg <- pairs_u %>%
      filter(subpair_name_raw != sp) %>%
      group_by(block, pair_type) %>%
      summarise(other_mean = mean(jaccard, na.rm = TRUE), .groups = "drop")
    
    merged <- inner_join(sp_block, other_bg, by = c("block", "pair_type"))
    if (nrow(merged) < 3) return(NULL)
    
    wt <- safe_paired_wilcox(merged$sp_mean, merged$other_mean)
    
    tibble(
      subpair_name_raw = sp,
      subpair_name_clean = sp_clean,
      pair_type = unique(as.character(merged$pair_type))[1],
      n_edges = nrow(sub),
      n_blocks = nrow(merged),
      mean_jaccard = mean(sub$jaccard, na.rm = TRUE),
      median_jaccard = median(sub$jaccard, na.rm = TRUE),
      mean_diff_vs_bg = wt$mean_diff,
      median_diff_vs_bg = wt$median_diff,
      p_value = wt$p_value
    )
  }) %>%
  mutate(fdr_bh = p.adjust(p_value, method = "BH"))

ei_subpairs_sig <- subpair_tests %>%
  filter(pair_type == "E-I", is.finite(fdr_bh), fdr_bh < 0.05) %>%
  mutate(
    direction = ifelse(mean_diff_vs_bg > 0, "Higher than block E-I background", "Lower than block E-I background")
  ) %>%
  arrange(desc(abs(mean_diff_vs_bg)))

# ---------------------------
# Key result extraction for the 4 chat conclusions
# ---------------------------
get_test_row <- function(df, metric_name, comp_name) {
  out <- df %>% filter(metric == metric_name, comparison == comp_name)
  if (nrow(out) == 0) {
    return(tibble(metric = metric_name, comparison = comp_name, n = NA_real_, mean_diff = NA_real_, median_diff = NA_real_, rank_biserial = NA_real_, p_value = NA_real_, fdr_bh = NA_real_))
  }
  out %>% slice(1)
}

res_overlap_EE_EI <- get_test_row(overall_tests, "block_mean_overlap_prop_small", "E-E vs E-I")
res_overlap_EI_II <- get_test_row(overall_tests, "block_mean_overlap_prop_small", "E-I vs I-I")
res_gap_EE_EI <- get_test_row(overall_tests, "block_mean_nesting_gap", "E-E vs E-I")
res_gap_EI_II <- get_test_row(overall_tests, "block_mean_nesting_gap", "E-I vs I-I")

get_layer_nested_row <- function(layer_name, comp_name = "E-E vs E-I") {
  out <- nested_layer_tests %>% filter(as.character(layer) == layer_name, comparison == comp_name)
  if (nrow(out) == 0) {
    return(tibble(metric = "prop_nested90", comparison = comp_name, n = NA_real_, mean_diff = NA_real_, median_diff = NA_real_, rank_biserial = NA_real_, p_value = NA_real_, layer = layer_name, fdr_bh = NA_real_))
  }
  out %>% slice(1)
}

res_nested_l23 <- get_layer_nested_row("layer 2/3")
res_nested_l4 <- get_layer_nested_row("layer 4")

focus_006 <- exc_sub_tests %>% filter(str_detect(subclass, "^006\\s")) %>% slice(1)
focus_007_l23 <- exc_sub_layer_tests %>% filter(str_detect(subclass, "^007\\s"), as.character(layer) == "layer 2/3") %>% slice(1)
focus_007_l4 <- exc_sub_layer_tests %>% filter(str_detect(subclass, "^007\\s"), as.character(layer) == "layer 4") %>% slice(1)

focus_pair_1 <- pick_focus_row(ei_subpairs_sig, c("007", "L2/3 IT CTX", "Vip"))
focus_pair_2 <- pick_focus_row(ei_subpairs_sig, c("032", "L5 NP CTX", "Pvalb"))
focus_pair_3 <- pick_focus_row(ei_subpairs_sig, c("032", "L5 NP CTX", "Sst"))
focus_pair_4 <- pick_focus_row(ei_subpairs_sig, c("002", "IT EP-CLA", "Sst"))
focus_pair_5 <- pick_focus_row(ei_subpairs_sig, c("021", "L4 RSP-ACA", "Pvalb"))
focus_pair_6 <- pick_focus_row(ei_subpairs_sig, c("007", "L2/3 IT CTX", "Lamp5"))
focus_pair_7 <- pick_focus_row(ei_subpairs_sig, c("004", "L6 IT CTX", "Vip"))

key_results <- bind_rows(
  tibble(
    section = "1_geometry",
    item = c(
      "E-I smaller-cluster coverage > E-E",
      "E-I smaller-cluster coverage > I-I",
      "E-I nesting gap > E-E",
      "E-I nesting gap > I-I"
    ),
    estimate = c(
      -res_overlap_EE_EI$mean_diff,
      res_overlap_EI_II$mean_diff,
      -res_gap_EE_EI$mean_diff,
      res_gap_EI_II$mean_diff
    ),
    p_value = c(res_overlap_EE_EI$p_value, res_overlap_EI_II$p_value, res_gap_EE_EI$p_value, res_gap_EI_II$p_value),
    fdr_bh = c(res_overlap_EE_EI$fdr_bh, res_overlap_EI_II$fdr_bh, res_gap_EE_EI$fdr_bh, res_gap_EI_II$fdr_bh),
    note = c(
      "block-level mean smaller overlap proportion",
      "block-level mean smaller overlap proportion",
      "block-level mean nesting gap",
      "block-level mean nesting gap"
    )
  ),
  tibble(
    section = "2_layer_near_complete_nesting",
    item = c(
      "layer 2/3: E-I near-complete nesting > E-E",
      "layer 4: E-I near-complete nesting > E-E"
    ),
    estimate = c(-res_nested_l23$mean_diff, -res_nested_l4$mean_diff),
    p_value = c(res_nested_l23$p_value, res_nested_l4$p_value),
    fdr_bh = c(res_nested_l23$fdr_bh, res_nested_l4$fdr_bh),
    note = c("prop overlap_small >= 0.90", "prop overlap_small >= 0.90")
  ),
  tibble(
    section = "3_exc_subclass_preference",
    item = c(
      if (nrow(focus_006) > 0) paste0(focus_006$subclass, ": same-type minus cross-type") else "006 subclass not found",
      if (nrow(focus_007_l23) > 0) paste0(focus_007_l23$subclass, " in layer 2/3") else "007 layer 2/3 not found",
      if (nrow(focus_007_l4) > 0) paste0(focus_007_l4$subclass, " in layer 4") else "007 layer 4 not found"
    ),
    estimate = c(
      if (nrow(focus_006) > 0) focus_006$mean_same_minus_cross else NA_real_,
      if (nrow(focus_007_l23) > 0) focus_007_l23$mean_same_minus_cross else NA_real_,
      if (nrow(focus_007_l4) > 0) focus_007_l4$mean_same_minus_cross else NA_real_
    ),
    p_value = c(
      if (nrow(focus_006) > 0) focus_006$p_value else NA_real_,
      if (nrow(focus_007_l23) > 0) focus_007_l23$p_value else NA_real_,
      if (nrow(focus_007_l4) > 0) focus_007_l4$p_value else NA_real_
    ),
    fdr_bh = c(
      if (nrow(focus_006) > 0) focus_006$fdr_bh else NA_real_,
      if (nrow(focus_007_l23) > 0) focus_007_l23$fdr_bh else NA_real_,
      if (nrow(focus_007_l4) > 0) focus_007_l4$fdr_bh else NA_real_
    ),
    note = c("overall excitatory subclass test", "layer-specific raw p is most informative", "layer-specific raw p is most informative")
  ),
  bind_rows(focus_pair_1, focus_pair_2, focus_pair_3, focus_pair_4, focus_pair_5, focus_pair_6, focus_pair_7) %>%
    transmute(
      section = "4_specific_EI_pairs",
      item = subpair_name_raw,
      estimate = mean_diff_vs_bg,
      p_value = p_value,
      fdr_bh = fdr_bh,
      note = direction
    )
)

# ---------------------------
# Export tables
# ---------------------------
dataset_overview <- tibble(
  metric = c(
    "input_csv", "raw directional rows", "unique undirected pairs", "unique clusters",
    "Glut clusters", "Gaba clusters", "slides", "layers", "complete slide-layer blocks"
  ),
  value = c(
    basename(input_csv), nrow(df), nrow(pairs_u), nrow(clusters),
    sum(clusters$cell_Neruon_type == "Glut"),
    sum(clusters$cell_Neruon_type == "Gaba"),
    n_distinct(clusters$slide),
    n_distinct(clusters$layer),
    length(blocks_complete)
  )
)

pairtype_desc <- pairs_u %>%
  group_by(pair_type) %>%
  summarise(
    n_pairs = n(),
    mean_jaccard = mean(jaccard, na.rm = TRUE),
    median_jaccard = median(jaccard, na.rm = TRUE),
    mean_overlap_prop_small = mean(overlap_prop_small, na.rm = TRUE),
    median_overlap_prop_small = median(overlap_prop_small, na.rm = TRUE),
    mean_overlap_prop_large = mean(overlap_prop_large, na.rm = TRUE),
    mean_nesting_gap = mean(nesting_gap, na.rm = TRUE),
    prop_nested90 = mean(nested_90, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(dataset_overview, file.path(tbl_dir, "table_00_dataset_overview.csv"))
write_csv(pairtype_desc, file.path(tbl_dir, "table_01_pairtype_descriptives.csv"))
write_csv(overall_tests, file.path(tbl_dir, "table_02_overall_geometry_tests.csv"))
write_csv(nested_layer_tests, file.path(tbl_dir, "table_03_layer_nested90_tests.csv"))
write_csv(exc_sub_tests, file.path(tbl_dir, "table_04_exc_subclass_preference.csv"))
write_csv(exc_sub_layer_tests, file.path(tbl_dir, "table_05_exc_subclass_layerwise_preference.csv"))
write_csv(gaba_sub_tests, file.path(tbl_dir, "table_06_gaba_subclass_preference.csv"))
write_csv(ei_subpairs_sig, file.path(tbl_dir, "table_07_significant_EI_subpairs_vs_background.csv"))
write_csv(key_results, file.path(tbl_dir, "table_08_key_results_from_chat_conclusions.csv"))

# ---------------------------
# Figure 1: spatial geometry is different, not simply stronger
# ---------------------------
geom_small_long <- block_metrics_complete %>%
  select(block, pair_type, value = mean_overlap_prop_small)

geom_gap_long <- block_metrics_complete %>%
  select(block, pair_type, value = mean_nesting_gap)

geom_scatter <- block_metrics_complete %>%
  select(block, pair_type, mean_overlap_prop_small, mean_overlap_prop_large)

ann_small <- c(
  paste0("E-I > E-E: FDR = ", fmt_p(res_overlap_EE_EI$fdr_bh)),
  paste0("E-I > I-I: FDR = ", fmt_p(res_overlap_EI_II$fdr_bh))
)
ann_gap <- c(
  paste0("E-I > E-E: FDR = ", fmt_p(res_gap_EE_EI$fdr_bh)),
  paste0("E-I > I-I: FDR = ", fmt_p(res_gap_EI_II$fdr_bh))
)

p1a <- ggplot(geom_small_long, aes(x = pair_type, y = value, color = pair_type)) +
  geom_line(aes(group = block), color = "grey78", alpha = 0.22, linewidth = 0.35) +
  geom_boxplot(aes(fill = pair_type), width = 0.34, alpha = 0.16, outlier.shape = NA, color = "black", linewidth = 0.45) +
  geom_point(position = position_jitter(width = 0.06, height = 0), size = 1.1, alpha = 0.55) +
  stat_summary(aes(group = 1), fun = median, geom = "line", linewidth = 1.15, color = "black") +
  stat_summary(fun = median, geom = "point", size = 3.1, color = "black") +
  scale_color_manual(values = pair_cols) +
  scale_fill_manual(values = pair_cols) +
  scale_y_continuous(labels = label_number(accuracy = 0.01), expand = expansion(mult = c(0.05, 0.18))) +
  labs(
    title = "Smaller-cluster coverage is highest for E-I",
    subtitle = "Block-level paired comparison across slide×layer units",
    x = NULL,
    y = "Mean smaller-cluster overlap proportion"
  ) +
  annotate("text", x = 1.05, y = max(geom_small_long$value, na.rm = TRUE) * 1.12, label = ann_small[1], hjust = 0, size = 3.7) +
  annotate("text", x = 1.05, y = max(geom_small_long$value, na.rm = TRUE) * 1.04, label = ann_small[2], hjust = 0, size = 3.7)

p1b <- ggplot(geom_gap_long, aes(x = pair_type, y = value, color = pair_type)) +
  geom_line(aes(group = block), color = "grey78", alpha = 0.22, linewidth = 0.35) +
  geom_boxplot(aes(fill = pair_type), width = 0.34, alpha = 0.16, outlier.shape = NA, color = "black", linewidth = 0.45) +
  geom_point(position = position_jitter(width = 0.06, height = 0), size = 1.1, alpha = 0.55) +
  stat_summary(aes(group = 1), fun = median, geom = "line", linewidth = 1.15, color = "black") +
  stat_summary(fun = median, geom = "point", size = 3.1, color = "black") +
  scale_color_manual(values = pair_cols) +
  scale_fill_manual(values = pair_cols) +
  scale_y_continuous(labels = label_number(accuracy = 0.01), expand = expansion(mult = c(0.05, 0.18))) +
  labs(
    title = "E-I shows the largest nesting gap",
    subtitle = "A larger gap means the smaller cluster is nested inside a much larger one",
    x = NULL,
    y = "Mean nesting gap"
  ) +
  annotate("text", x = 1.05, y = max(geom_gap_long$value, na.rm = TRUE) * 1.12, label = ann_gap[1], hjust = 0, size = 3.7) +
  annotate("text", x = 1.05, y = max(geom_gap_long$value, na.rm = TRUE) * 1.04, label = ann_gap[2], hjust = 0, size = 3.7)

p1c <- ggplot(geom_scatter, aes(x = mean_overlap_prop_large, y = mean_overlap_prop_small, color = pair_type, fill = pair_type)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.55, color = "grey50") +
  stat_ellipse(geom = "polygon", alpha = 0.10, linewidth = 0.55, level = 0.70, show.legend = FALSE) +
  geom_point(size = 2.0, alpha = 0.78, stroke = 0) +
  stat_summary(aes(group = pair_type), fun = mean, geom = "point", size = 4.2, shape = 21, color = "black", stroke = 0.55, show.legend = FALSE) +
  scale_color_manual(values = pair_cols) +
  scale_fill_manual(values = pair_cols) +
  scale_x_continuous(labels = label_number(accuracy = 0.01), expand = expansion(mult = c(0.02, 0.03))) +
  scale_y_continuous(labels = label_number(accuracy = 0.01), expand = expansion(mult = c(0.02, 0.06))) +
  coord_equal() +
  labs(
    title = "The geometry is hierarchical rather than symmetric",
    x = "Mean overlap proportion of the larger cluster",
    y = "Mean overlap proportion of the smaller cluster"
  )

save_pub(p1a + theme(legend.position = "none"), "fig01a_smaller_cluster_coverage", 5.6, 5.2)
save_pub(p1b + theme(legend.position = "none"), "fig01b_nesting_gap", 8, 5.2)
save_pub(p1c + theme(legend.position = "none"), "fig01c_hierarchical_geometry_scatter", 8, 5.2)

# ---------------------------
# Figure 2: near-complete nesting in layer 2/3 and layer 4
# ---------------------------
layer_heat <- block_metrics_complete %>%
  group_by(layer, pair_type) %>%
  summarise(value = mean(prop_nested90, na.rm = TRUE), .groups = "drop")

p2a <- ggplot(layer_heat, aes(x = pair_type, y = layer, fill = value)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", value)), size = 4.2, fontface = "bold") +
  scale_fill_gradientn(colours = c("#F4F6FB", "#C9D9F2", "#8FB3E3", "#4F84C4", "#1F4E79")) +
  labs(
    title = "Near-complete nesting is strongest in shallow-to-middle layers",
    subtitle = "Definition: smaller-cluster overlap ≥ 90%",
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  theme(legend.position = "right")

make_layer_nested_plot <- function(layer_name, test_row) {
  dat <- block_metrics_complete %>% filter(as.character(layer) == layer_name)
  if (nrow(dat) == 0) return(make_empty_plot(paste0(layer_name, " not found")))
  
  ggplot(dat, aes(x = pair_type, y = prop_nested90, color = pair_type)) +
    geom_line(aes(group = block), color = "grey78", alpha = 0.25, linewidth = 0.35) +
    geom_boxplot(aes(fill = pair_type), width = 0.34, alpha = 0.16, outlier.shape = NA, color = "black", linewidth = 0.45) +
    geom_point(position = position_jitter(width = 0.06, height = 0), size = 1.2, alpha = 0.62) +
    stat_summary(aes(group = 1), fun = median, geom = "line", linewidth = 1.12, color = "black") +
    stat_summary(fun = median, geom = "point", size = 3.0, color = "black") +
    scale_color_manual(values = pair_cols) +
    scale_fill_manual(values = pair_cols) +
    scale_y_continuous(labels = label_percent(accuracy = 1), expand = expansion(mult = c(0.05, 0.18))) +
    labs(
      title = paste0(layer_name, ": E-I near-complete nesting"),
      subtitle = paste0("E-I > E-E: FDR = ", fmt_p(test_row$fdr_bh)),
      x = NULL,
      y = "Fraction of pairs with overlap_small ≥ 90%"
    )
}

p2b <- make_layer_nested_plot("layer 2/3", res_nested_l23)
p2c <- make_layer_nested_plot("layer 4", res_nested_l4)

save_pub(p2a, "fig02a_near_complete_nesting_heatmap", 8, 5.3)
save_pub(p2b + theme(legend.position = "none"), "fig02b_layer23_near_complete_nesting", 5.6, 5.3)
save_pub(p2c + theme(legend.position = "none"), "fig02c_layer4_near_complete_nesting", 5.6, 5.3)

# ---------------------------
# Figure 3: excitatory subclass partner preference
# ---------------------------
exc_plot_df <- exc_sub_tests %>%
  mutate(
    preference = ifelse(mean_same_minus_cross < 0, "cross-type preferring", "same-type preferring"),
    label_flag = case_when(
      str_detect(subclass, "^006\\s") ~ TRUE,
      str_detect(subclass, "^007\\s") ~ TRUE,
      fdr_bh < 0.05 & mean_same_minus_cross < 0 ~ TRUE,
      TRUE ~ FALSE
    ),
    subclass_plot = fct_reorder(subclass, mean_same_minus_cross)
  )

p3a <- ggplot(exc_plot_df, aes(x = mean_same_minus_cross, y = subclass_plot, color = preference)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.55, color = "grey45") +
  geom_segment(aes(x = 0, xend = mean_same_minus_cross, y = subclass_plot, yend = subclass_plot), linewidth = 0.85, alpha = 0.85) +
  geom_point(aes(size = ifelse(fdr_bh < 0.05, 3.3, 2.3)), stroke = 0.5, fill = "white", shape = 21) +
  geom_text_repel(
    data = exc_plot_df %>% filter(label_flag),
    aes(label = paste0(subclass, "\nΔ=", fmt_num(mean_same_minus_cross, 4), "; FDR=", fmt_p(fdr_bh))),
    size = 3.3,
    min.segment.length = 0,
    segment.alpha = 0.6,
    box.padding = 0.28,
    point.padding = 0.15,
    max.overlaps = Inf,
    seed = 123
  ) +
  scale_color_manual(values = pref_cols) +
  scale_size_identity() +
  labs(
    title = "Some common excitatory subclasses prefer inhibitory partners",
    subtitle = "Negative values indicate stronger cross-type than same-type spatial coupling",
    x = "same-type mean Jaccard − cross-type mean Jaccard",
    y = NULL
  )

focus_007_dat <- exc_pref %>%
  filter(str_detect(subclass, "^007\\s"), as.character(layer) %in% c("layer 2/3", "layer 4")) %>%
  mutate(layer = factor(as.character(layer), levels = c("layer 2/3", "layer 4")))

focus_007_ann <- tibble(
  layer = factor(c("layer 2/3", "layer 4"), levels = c("layer 2/3", "layer 4")),
  p_lab = c(
    if (nrow(focus_007_l23) > 0) paste0("p = ", fmt_p(focus_007_l23$p_value)) else "p = NA",
    if (nrow(focus_007_l4) > 0) paste0("p = ", fmt_p(focus_007_l4$p_value)) else "p = NA"
  )
)

if (nrow(focus_007_dat) > 0) {
  p3b <- ggplot(focus_007_dat, aes(x = layer, y = same_minus_cross)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.55, color = "grey45") +
    geom_boxplot(width = 0.40, outlier.shape = NA, fill = "#E8EEF7", color = "black", linewidth = 0.5) +
    geom_jitter(aes(color = layer), width = 0.08, height = 0, size = 2.0, alpha = 0.75, show.legend = FALSE) +
    stat_summary(fun = mean, geom = "point", size = 3.6, shape = 23, fill = "#D1495B", color = "black") +
    geom_text(
      data = focus_007_ann,
      aes(x = layer, y = max(focus_007_dat$same_minus_cross, na.rm = TRUE) + 0.012, label = p_lab),
      inherit.aes = FALSE,
      size = 3.7,
      fontface = "bold"
    ) +
    scale_color_manual(values = c("layer 2/3" = "#3C78D8", "layer 4" = "#D1495B")) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.18))) +
    labs(
      title = "007 L2/3 IT CTX Glut keeps the cross-type tendency in layer 2/3 and layer 4",
      subtitle = "Cluster-level paired test within the focal excitatory subclass",
      x = NULL,
      y = "same-type mean Jaccard − cross-type mean Jaccard"
    )
} else {
  p3b <- make_empty_plot("Target subclass 007 was not found in the current dataset")
}

save_pub(p3a + theme(legend.position = "none"), "fig03a_exc_subclass_partner_preference", 8.4, 6.6)
save_pub(p3b + theme(legend.position = "none"), "fig03b_subclass_007_layerwise_preference", 8, 6.6)

# ---------------------------
# Figure 4: specific E-I pairs, not generic inhibitory roles
# 优先显示聊天中点名的配对；如果缺失，则回退为按绝对效应排序
# ---------------------------
focus_pairs_from_chat <- bind_rows(
  focus_pair_1, focus_pair_2, focus_pair_3, focus_pair_4,
  focus_pair_5, focus_pair_6, focus_pair_7
) %>% distinct(subpair_name_raw, .keep_all = TRUE)

if (nrow(focus_pairs_from_chat) < 4) {
  pos_top <- ei_subpairs_sig %>% filter(mean_diff_vs_bg > 0) %>% arrange(desc(mean_diff_vs_bg)) %>% slice_head(n = 4)
  neg_top <- ei_subpairs_sig %>% filter(mean_diff_vs_bg < 0) %>% arrange(mean_diff_vs_bg) %>% slice_head(n = 4)
  fig4_df <- bind_rows(pos_top, neg_top)
} else {
  fig4_df <- focus_pairs_from_chat
}

fig4_df <- fig4_df %>%
  mutate(
    pair_label = str_wrap(subpair_name_raw, width = 34),
    pair_label = fct_reorder(pair_label, mean_diff_vs_bg),
    text_hjust = ifelse(mean_diff_vs_bg >= 0, 0, 1),
    x_text = mean_diff_vs_bg + ifelse(mean_diff_vs_bg >= 0, 0.006, -0.006)
  )

p4 <- ggplot(fig4_df, aes(x = mean_diff_vs_bg, y = pair_label, color = direction)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.55, color = "grey45") +
  geom_segment(aes(x = 0, xend = mean_diff_vs_bg, y = pair_label, yend = pair_label), linewidth = 1.0, alpha = 0.9) +
  geom_point(size = 3.2) +
  geom_text(
    aes(x = x_text, hjust = text_hjust, label = paste0("ΔJ = ", fmt_num(mean_diff_vs_bg, 4), "\nFDR = ", fmt_p(fdr_bh))),
    size = 3.45,
    lineheight = 0.95,
    show.legend = FALSE
  ) +
  scale_color_manual(values = bg_cols) +
  scale_x_continuous(expand = expansion(mult = c(0.16, 0.24))) +
  labs(
    title = "Specific E-I pairs can be above or below the block-level E-I background",
    subtitle = "This argues against treating inhibitory subclasses as a single uniform spatial role",
    x = "Δ mean Jaccard vs within-block E-I background",
    y = NULL
  )

save_pub(p4, "fig04_specific_EI_pairs_above_or_below_background", 12.0, 7.2)

# ---------------------------
# Figure manifest
# ---------------------------
figure_manifest <- tibble(
  file_stem = c(
    "fig01_EI_geometry_not_symmetric_but_nested",
    "fig02_near_complete_nesting_layer23_layer4",
    "fig03_exc_subclass_prefers_inhibitory_partners",
    "fig04_specific_EI_pairs_above_or_below_background"
  ),
  content = c(
    "主图1：证明三类互作不是简单强弱差异，而是空间几何类型不同；重点突出E-I在small-cluster coverage和nesting gap上最高。",
    "主图2：展示near-complete nesting在不同layer中的分布，并重点强调layer 2/3和layer 4中E-I > E-E。",
    "主图3：展示兴奋性亚类的伙伴偏好，突出006总体更偏向抑制伙伴，以及007在layer 2/3和layer 4中的同样趋势。",
    "主图4：展示具体E-I配对相对block内E-I背景的偏离，强调空间耦合由具体配对决定，而不是由抑制亚类统一决定。"
  )
)

write_csv(figure_manifest, file.path(meta_dir, "figure_manifest.csv"))
write_lines(
  c(
    paste0("Input CSV: ", normalizePath(input_csv, winslash = "/", mustWork = FALSE)),
    paste0("Output dir: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE)),
    "",
    "Generated figures:",
    paste0("- ", figure_manifest$file_stem, ": ", figure_manifest$content)
  ),
  file.path(meta_dir, "run_summary.txt")
)

# ---------------------------
# Optional cluster-robust regression table retained
# ---------------------------
mod_df <- pairs_u %>%
  mutate(
    jaccard_clip = pmin(pmax(jaccard, 1e-4), 1 - 1e-4),
    logit_jaccard = log(jaccard_clip / (1 - jaccard_clip)),
    log_size_small = log(size_small),
    log_size_ratio = log(size_ratio),
    layer_chr = as.character(layer),
    pair_type_chr = as.character(pair_type)
  )

fit <- lm(logit_jaccard ~ pair_type_chr + layer_chr + log_size_small + log_size_ratio, data = mod_df)
robust_vcov <- sandwich::vcovCL(fit, cluster = mod_df$block)
ct <- lmtest::coeftest(fit, vcov. = robust_vcov)

reg_table <- data.frame(
  term      = rownames(ct),
  estimate  = as.numeric(ct[, 1]),
  std_error = as.numeric(ct[, 2]),
  statistic = as.numeric(ct[, 3]),
  p_value   = as.numeric(ct[, 4]),
  row.names = NULL,
  check.names = FALSE
)

readr::write_csv(reg_table, file.path(tbl_dir, "table_09_cluster_robust_regression.csv"))

message("All outputs are ready under: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
