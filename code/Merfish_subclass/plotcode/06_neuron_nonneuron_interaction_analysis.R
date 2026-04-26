# 清空环境
rm(list = ls())
gc()
# =========================================================
# neuron_nonneuron_R_output_06_topjournal_fixed_v2.R
# Full reproducible workflow for neuron–nonneuron partner data
# Input : E:/zaw/2603/neuron-nonneuron-partner.csv
# Output: E:/zaw/2603/R_output_06
# =========================================================

required_pkgs <- c(
  "dplyr", "tidyr", "stringr", "purrr", "readr", "tibble",
  "ggplot2", "forcats", "broom", "sandwich", "lmtest", "scales"
)
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(forcats)
  library(broom)
  library(sandwich)
  library(lmtest)
  library(scales)
})

options(stringsAsFactors = FALSE, scipen = 999)

# Guard against function masking from other packages
select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter
mutate <- dplyr::mutate
summarise <- dplyr::summarise
arrange <- dplyr::arrange
distinct <- dplyr::distinct
count <- dplyr::count
left_join <- dplyr::left_join
inner_join <- dplyr::inner_join
transmute <- dplyr::transmute
pull <- dplyr::pull
first <- dplyr::first
n_distinct <- dplyr::n_distinct

# ---------------------------
# Paths
# ---------------------------
root_dir <- "E:/zaw/2603"
input_csv <- file.path(root_dir, "neuron-nonneuron-partner.csv")
out_dir <- file.path(root_dir, "R_output_06")
fig_dir <- file.path(out_dir, "figures")
tbl_dir <- file.path(out_dir, "tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tbl_dir, showWarnings = FALSE, recursive = TRUE)
if (!file.exists(input_csv)) stop("Input file not found: ", input_csv)

# ---------------------------
# Helpers
# ---------------------------
canonical_pair <- function(a, b) {
  ifelse(a < b, paste(a, b, sep = " || "), paste(b, a, sep = " || "))
}

clean_subclass <- function(x) {
  stringr::str_replace(as.character(x), "^\\d+\\s+", "")
}

rank_biserial_paired <- function(x, y) {
  d <- x - y
  d <- d[!is.na(d) & d != 0]
  if (length(d) == 0) return(NA_real_)
  (sum(d > 0) - sum(d < 0)) / length(d)
}

paired_wilcox_tbl <- function(dat, a, b) {
  tmp <- dat %>% dplyr::select(all_of(c(a, b))) %>% tidyr::drop_na()
  if (nrow(tmp) == 0) return(NULL)
  wt <- suppressWarnings(wilcox.test(tmp[[a]], tmp[[b]], paired = TRUE, exact = FALSE))
  tibble(
    n = nrow(tmp),
    mean_diff = mean(tmp[[a]] - tmp[[b]]),
    median_diff = median(tmp[[a]] - tmp[[b]]),
    rank_biserial = rank_biserial_paired(tmp[[a]], tmp[[b]]),
    p_value = wt$p.value
  )
}

tidy_cluster_model <- function(fit, cluster_vec, logistic = TRUE) {
  vc <- sandwich::vcovCL(fit, cluster = cluster_vec)
  ct <- lmtest::coeftest(fit, vcov. = vc)
  out <- tibble(
    term = rownames(ct),
    estimate = as.numeric(ct[, 1]),
    std_error = as.numeric(ct[, 2]),
    stat = as.numeric(ct[, 3]),
    p_value = as.numeric(ct[, 4])
  )
  out <- out %>%
    dplyr::mutate(
      conf_low = estimate - 1.96 * std_error,
      conf_high = estimate + 1.96 * std_error
    )
  if (logistic) {
    out <- out %>%
      dplyr::mutate(
        OR = exp(estimate),
        OR_low = exp(conf_low),
        OR_high = exp(conf_high)
      )
  }
  out
}

fmt_p_short <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-300) return("<1e-300")
  if (p < 1e-4) return(format(p, scientific = TRUE, digits = 2))
  sprintf("%.4f", p)
}

wrap_plot_title <- function(x, width = 54) stringr::str_wrap(x, width = width)
wrap_plot_subtitle <- function(x, width = 92) stringr::str_wrap(x, width = width)

save_pub <- function(p, filename, width, height, dpi = 600) {
  ggsave(
    filename = file.path(fig_dir, paste0(filename, ".png")),
    plot = p, width = width, height = height, dpi = dpi,
    bg = "white", limitsize = FALSE
  )
  ggsave(
    filename = file.path(fig_dir, paste0(filename, ".pdf")),
    plot = p, width = width, height = height,
    bg = "white", limitsize = FALSE
  )
}

pick_pairs_by_pattern <- function(tab, patterns, sign = c("both", "pos", "neg")) {
  sign <- match.arg(sign)
  out <- purrr::map_dfr(patterns, function(pt) {
    x <- tab %>% dplyr::filter(stringr::str_detect(subpair_name_clean, stringr::regex(pt, ignore_case = TRUE)))
    if (nrow(x) == 0) return(tibble())
    if (sign == "pos") x <- x %>% dplyr::filter(mean_diff_vs_bg > 0)
    if (sign == "neg") x <- x %>% dplyr::filter(mean_diff_vs_bg < 0)
    if (nrow(x) == 0) return(tibble())
    x %>% dplyr::arrange(p_value, dplyr::desc(abs(mean_diff_vs_bg))) %>% dplyr::slice(1)
  })
  out %>% dplyr::distinct(subpair_name_clean, .keep_all = TRUE)
}

plot_theme <- ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    plot.title.position = "plot",
    plot.margin = ggplot2::margin(18, 28, 14, 18),
    plot.title = ggplot2::element_text(face = "bold", size = 16.5, lineheight = 1.06, margin = ggplot2::margin(b = 8)),
    plot.subtitle = ggplot2::element_text(size = 11.8, colour = "#4D4D4D", lineheight = 1.10, margin = ggplot2::margin(b = 12)),
    axis.title = ggplot2::element_text(size = 13.5, colour = "black"),
    axis.text = ggplot2::element_text(size = 11.5, colour = "black"),
    legend.position = "top",
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 11.2)
  )

theme_set(plot_theme)

col_neuron <- c("Glut" = "#C44E52", "Gaba" = "#4E79A7")
col_guest <- c("Not exact-contained" = "#BDBDBD", "Exact-contained" = "#C44E52")
col_posneg <- c("Higher than background" = "#C44E52", "Lower than background" = "#4E79A7")
col_shared <- c("Shared by both" = "#7A4FB3", "One neuron type only" = "#BDBDBD")

# ---------------------------
# Read data and build pair table
# ---------------------------
df <- readr::read_csv(input_csv, show_col_types = FALSE) %>%
  dplyr::mutate(pair_key = canonical_pair(cluster.1_label, cluster.2_label))

pairs <- df %>%
  dplyr::group_by(pair_key) %>%
  dplyr::filter(cluster.1_label == pmin(cluster.1_label, cluster.2_label)) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    label_a = cluster.1_label,
    label_b = cluster.2_label,
    slide = cluster.1_slide,
    layer = cluster.1_layer,
    region_a = cluster.1_region,
    region_b = cluster.2_region,
    type_a_raw = cluster.1_cell_Neruon_type,
    type_b_raw = cluster.2_cell_Neruon_type,
    subclass_a = cluster.1_subclass,
    subclass_b = cluster.2_subclass,
    subclass_a_clean = clean_subclass(cluster.1_subclass),
    subclass_b_clean = clean_subclass(cluster.2_subclass),
    size_a = cluster.1_total_cell_num,
    size_b = cluster.2_total_cell_num,
    overlap = overlap_cell,
    union = union_cell,
    jaccard = jaccard,
    overlap_prop_a = cluster.1.overlap.percent,
    overlap_prop_b = cluster.2.overlap.percent
  ) %>%
  dplyr::mutate(
    type_pair = purrr::map2_chr(type_a_raw, type_b_raw, ~ paste(sort(c(.x, .y)), collapse = " × ")),
    block = paste(slide, layer, sep = "|"),
    size_small = pmin(size_a, size_b),
    size_large = pmax(size_a, size_b),
    size_ratio = size_large / size_small,
    overlap_prop_small = overlap / size_small,
    overlap_prop_large = overlap / size_large,
    nesting_gap = overlap_prop_small - overlap_prop_large,
    exact_small_contained = overlap_prop_small == 1,
    near95_small_contained = overlap_prop_small >= 0.95,
    near90_small_contained = overlap_prop_small >= 0.90
  )

clusters <- dplyr::bind_rows(
  pairs %>% dplyr::transmute(label = label_a, slide, layer, type_raw = type_a_raw, subclass = subclass_a, subclass_clean = subclass_a_clean, size = size_a),
  pairs %>% dplyr::transmute(label = label_b, slide, layer, type_raw = type_b_raw, subclass = subclass_b, subclass_clean = subclass_b_clean, size = size_b)
) %>% dplyr::distinct() %>% dplyr::mutate(block = paste(slide, layer, sep = "|"))

mixed <- pairs %>%
  dplyr::filter((type_a_raw == "NonNeuron" & type_b_raw %in% c("Glut", "Gaba")) |
                  (type_b_raw == "NonNeuron" & type_a_raw %in% c("Glut", "Gaba"))) %>%
  dplyr::mutate(
    nn_on_a = type_a_raw == "NonNeuron",
    neuron_type = ifelse(nn_on_a, type_b_raw, type_a_raw),
    nn_subclass_full = ifelse(nn_on_a, subclass_a, subclass_b),
    nn_subclass = ifelse(nn_on_a, subclass_a_clean, subclass_b_clean),
    neuron_subclass_full = ifelse(nn_on_a, subclass_b, subclass_a),
    neuron_subclass = ifelse(nn_on_a, subclass_b_clean, subclass_a_clean),
    size_nn = ifelse(nn_on_a, size_a, size_b),
    size_neuron = ifelse(nn_on_a, size_b, size_a),
    overlap_prop_nn = ifelse(nn_on_a, overlap_prop_a, overlap_prop_b),
    overlap_prop_neuron = ifelse(nn_on_a, overlap_prop_b, overlap_prop_a),
    label_nn = ifelse(nn_on_a, label_a, label_b),
    label_neuron = ifelse(nn_on_a, label_b, label_a),
    nn_in_neuron_exact = overlap_prop_nn == 1,
    neuron_in_nn_exact = overlap_prop_neuron == 1,
    exact_any = nn_in_neuron_exact | neuron_in_nn_exact,
    both_exact = nn_in_neuron_exact & neuron_in_nn_exact,
    smaller_identity = dplyr::case_when(
      size_nn < size_neuron ~ "NN",
      size_neuron < size_nn ~ "Neuron",
      TRUE ~ "Equal"
    ),
    log_size_neuron = log(size_neuron),
    log_size_nn = log(size_nn),
    log_size_ratio = log(size_ratio),
    jaccard_clip = pmin(pmax(jaccard, 1e-4), 1 - 1e-4),
    logit_jaccard = log(jaccard_clip / (1 - jaccard_clip)),
    subpair_name = paste(neuron_subclass_full, nn_subclass_full, sep = " || "),
    subpair_name_clean = paste(neuron_subclass, nn_subclass, sep = " × ")
  )

major_nn <- mixed %>% dplyr::count(nn_subclass, sort = TRUE) %>% dplyr::filter(n >= 100) %>% dplyr::pull(nn_subclass)

# ---------------------------
# Pair-specific background tests
# ---------------------------
subpair_bg_jaccard <- mixed %>%
  dplyr::group_by(neuron_type, subpair_name) %>%
  dplyr::group_split() %>%
  purrr::map_dfr(function(sub) {
    nt <- unique(sub$neuron_type)
    sp_name <- unique(sub$subpair_name)
    sp_clean <- unique(sub$subpair_name_clean)
    if (nrow(sub) < 20 || dplyr::n_distinct(sub$block) < 10) return(NULL)
    sp <- sub %>% dplyr::group_by(block) %>% dplyr::summarise(sp_mean = mean(jaccard), .groups = "drop")
    bg <- mixed %>% dplyr::filter(neuron_type == nt, subpair_name != sp_name) %>%
      dplyr::group_by(block) %>% dplyr::summarise(bg_mean = mean(jaccard), .groups = "drop")
    merged <- dplyr::inner_join(sp, bg, by = "block")
    if (nrow(merged) < 10) return(NULL)
    wt <- suppressWarnings(wilcox.test(merged$sp_mean, merged$bg_mean, paired = TRUE, exact = FALSE))
    tibble(
      neuron_type = nt,
      subpair_name = sp_name,
      subpair_name_clean = sp_clean,
      n_edges = nrow(sub),
      n_blocks = nrow(merged),
      mean_jaccard = mean(sub$jaccard),
      mean_diff_vs_bg = mean(merged$sp_mean - merged$bg_mean),
      median_diff_vs_bg = median(merged$sp_mean - merged$bg_mean),
      p_value = wt$p.value
    )
  }) %>% dplyr::mutate(fdr_bh = p.adjust(p_value, method = "BH"))

subpair_bg_exact <- mixed %>%
  dplyr::group_by(neuron_type, subpair_name) %>%
  dplyr::group_split() %>%
  purrr::map_dfr(function(sub) {
    nt <- unique(sub$neuron_type)
    sp_name <- unique(sub$subpair_name)
    sp_clean <- unique(sub$subpair_name_clean)
    if (nrow(sub) < 20 || dplyr::n_distinct(sub$block) < 10) return(NULL)
    sp <- sub %>% dplyr::group_by(block) %>% dplyr::summarise(sp_rate = mean(nn_in_neuron_exact), .groups = "drop")
    bg <- mixed %>% dplyr::filter(neuron_type == nt, subpair_name != sp_name) %>%
      dplyr::group_by(block) %>% dplyr::summarise(bg_rate = mean(nn_in_neuron_exact), .groups = "drop")
    merged <- dplyr::inner_join(sp, bg, by = "block")
    if (nrow(merged) < 10) return(NULL)
    wt <- suppressWarnings(wilcox.test(merged$sp_rate, merged$bg_rate, paired = TRUE, exact = FALSE))
    tibble(
      neuron_type = nt,
      subpair_name = sp_name,
      subpair_name_clean = sp_clean,
      n_edges = nrow(sub),
      n_blocks = nrow(merged),
      exact_rate = mean(sub$nn_in_neuron_exact),
      mean_diff_vs_bg = mean(merged$sp_rate - merged$bg_rate),
      median_diff_vs_bg = median(merged$sp_rate - merged$bg_rate),
      p_value = wt$p.value
    )
  }) %>% dplyr::mutate(fdr_bh = p.adjust(p_value, method = "BH"))

# ---------------------------
# Cluster-level host / guest analyses
# ---------------------------
neuron_cluster_stats <- mixed %>%
  dplyr::group_by(label_neuron, neuron_type, neuron_subclass_full, neuron_subclass, slide, layer) %>%
  dplyr::summarise(
    neuron_size = dplyr::first(size_neuron),
    n_mixed_edges = dplyr::n(),
    mean_jaccard = mean(jaccard),
    mean_nesting_gap = mean(nesting_gap),
    n_exact_nn_guests = sum(nn_in_neuron_exact),
    prop_exact_nn_guests = mean(nn_in_neuron_exact),
    n_exact_neuron_in_nn = sum(neuron_in_nn_exact),
    n_unique_nn_subclasses = dplyr::n_distinct(nn_subclass),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    host_exact_any = n_exact_nn_guests > 0,
    multi_exact_nn_guests = n_exact_nn_guests >= 2,
    log_neuron_size = log(neuron_size),
    block = paste(slide, layer, sep = "|")
  )

nnnn_edges <- pairs %>% dplyr::filter(type_pair == "NonNeuron × NonNeuron")

nnnn_long <- dplyr::bind_rows(
  nnnn_edges %>% dplyr::transmute(label_nn = label_a, jaccard = jaccard, overlap_prop_self = overlap_prop_a),
  nnnn_edges %>% dplyr::transmute(label_nn = label_b, jaccard = jaccard, overlap_prop_self = overlap_prop_b)
)

nnnn_cluster <- nnnn_long %>%
  dplyr::group_by(label_nn) %>%
  dplyr::summarise(
    n_nnn_edges = dplyr::n(),
    mean_nnn_jaccard = mean(jaccard, na.rm = TRUE),
    mean_nnn_overlap_self = mean(overlap_prop_self, na.rm = TRUE),
    .groups = "drop"
  )

nn_cluster_stats <- mixed %>%
  dplyr::group_by(label_nn, nn_subclass_full, nn_subclass, slide, layer) %>%
  dplyr::summarise(
    nn_size = dplyr::first(size_nn),
    n_mixed_edges = dplyr::n(),
    mean_jaccard = mean(jaccard),
    mean_nesting_gap = mean(nesting_gap),
    n_exact_neuron_hosts = sum(nn_in_neuron_exact),
    prop_exact_neuron_hosts = mean(nn_in_neuron_exact),
    n_exact_nn_hosts = sum(neuron_in_nn_exact),
    n_unique_neuron_subclasses = dplyr::n_distinct(neuron_subclass_full),
    n_unique_neuron_types = dplyr::n_distinct(neuron_type),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    guest_exact_any = n_exact_neuron_hosts > 0,
    multi_host_exact = n_exact_neuron_hosts >= 2,
    log_nn_size = log(nn_size),
    block = paste(slide, layer, sep = "|")
  ) %>%
  dplyr::left_join(nnnn_cluster, by = "label_nn") %>%
  dplyr::mutate(
    n_nnn_edges = ifelse(is.na(n_nnn_edges), 0L, n_nnn_edges),
    mean_nnn_jaccard = ifelse(is.na(mean_nnn_jaccard), 0, mean_nnn_jaccard),
    mean_nnn_overlap_self = ifelse(is.na(mean_nnn_overlap_self), 0, mean_nnn_overlap_self),
    nn_subclass_model = ifelse(nn_subclass %in% major_nn, nn_subclass, "OtherRareNN")
  )

fit_host <- glm(
  host_exact_any ~ neuron_type + log_neuron_size + n_mixed_edges + layer,
  data = neuron_cluster_stats, family = binomial()
)
reg_host <- tidy_cluster_model(fit_host, neuron_cluster_stats$block, logistic = TRUE)

fit_guest <- glm(
  guest_exact_any ~ log_nn_size + n_mixed_edges + mean_nnn_jaccard + layer + nn_subclass_model,
  data = nn_cluster_stats, family = binomial()
)
reg_guest <- tidy_cluster_model(fit_guest, nn_cluster_stats$block, logistic = TRUE)

fit_multihost <- glm(
  multi_host_exact ~ log_nn_size + n_mixed_edges + mean_nnn_jaccard + layer + nn_subclass_model,
  data = nn_cluster_stats, family = binomial()
)
reg_multihost <- tidy_cluster_model(fit_multihost, nn_cluster_stats$block, logistic = TRUE)

exact_guest_both_types <- mixed %>%
  dplyr::filter(nn_in_neuron_exact) %>%
  dplyr::group_by(label_nn) %>%
  dplyr::summarise(
    n_exact_hosts = dplyr::n_distinct(label_neuron),
    n_glut_exact_hosts = sum(neuron_type == "Glut"),
    n_gaba_exact_hosts = sum(neuron_type == "Gaba"),
    has_glut_exact = any(neuron_type == "Glut"),
    has_gaba_exact = any(neuron_type == "Gaba"),
    nn_subclass = dplyr::first(nn_subclass),
    nn_size = dplyr::first(size_nn),
    block = dplyr::first(block),
    .groups = "drop"
  ) %>%
  dplyr::mutate(both_neuron_types_exact = has_glut_exact & has_gaba_exact)

exact_guest_both_types_summary <- exact_guest_both_types %>%
  dplyr::group_by(nn_subclass) %>%
  dplyr::summarise(
    n_exact_guest_clusters = dplyr::n(),
    both_neuron_types_rate = mean(both_neuron_types_exact),
    mean_n_exact_hosts = mean(n_exact_hosts),
    median_n_exact_hosts = median(n_exact_hosts),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(both_neuron_types_rate))

overall_shared_summary <- tibble(
  n_exact_guest_clusters = nrow(exact_guest_both_types),
  n_shared_by_both_types = sum(exact_guest_both_types$both_neuron_types_exact),
  shared_rate = mean(exact_guest_both_types$both_neuron_types_exact)
)

# ---------------------------
# Selected exemplars
# ---------------------------
jaccard_interest_patterns <- c(
  "L2/3 IT RSP.*Oligo",
  "L4 RSP-ACA.*Oligo",
  "L2/3 IT CTX.*Endo",
  "L4/5 IT CTX.*Astro-TE",
  "Lamp5.*VLMC",
  "Lamp5.*SMC"
)
exact_interest_patterns <- c(
  "Lamp5.*OPC",
  "Lamp5.*Endo",
  "Lamp5.*VLMC",
  "Lamp5.*SMC"
)

sel_jaccard_pairs <- pick_pairs_by_pattern(subpair_bg_jaccard, jaccard_interest_patterns, sign = "both") %>%
  dplyr::mutate(direction = ifelse(mean_diff_vs_bg > 0, "Higher than background", "Lower than background"))
if (nrow(sel_jaccard_pairs) < 8) {
  extra <- subpair_bg_jaccard %>% dplyr::filter(fdr_bh < 0.05) %>% dplyr::arrange(dplyr::desc(abs(mean_diff_vs_bg))) %>% dplyr::slice_head(n = 8)
  sel_jaccard_pairs <- dplyr::bind_rows(sel_jaccard_pairs, extra) %>% dplyr::distinct(subpair_name_clean, .keep_all = TRUE) %>%
    dplyr::slice_head(n = 10) %>% dplyr::mutate(direction = ifelse(mean_diff_vs_bg > 0, "Higher than background", "Lower than background"))
}

sel_exact_pairs <- pick_pairs_by_pattern(subpair_bg_exact, exact_interest_patterns, sign = "both") %>%
  dplyr::mutate(direction = ifelse(mean_diff_vs_bg > 0, "Higher than background", "Lower than background"))
if (nrow(sel_exact_pairs) < 6) {
  extra <- subpair_bg_exact %>% dplyr::filter(fdr_bh < 0.05) %>% dplyr::arrange(dplyr::desc(abs(mean_diff_vs_bg))) %>% dplyr::slice_head(n = 8)
  sel_exact_pairs <- dplyr::bind_rows(sel_exact_pairs, extra) %>% dplyr::distinct(subpair_name_clean, .keep_all = TRUE) %>%
    dplyr::slice_head(n = 8) %>% dplyr::mutate(direction = ifelse(mean_diff_vs_bg > 0, "Higher than background", "Lower than background"))
}

# ---------------------------
# Save tables
# ---------------------------
write_csv(reg_guest, file.path(tbl_dir, "table_01_reg_guest_exact_any.csv"))
write_csv(reg_host, file.path(tbl_dir, "table_02_reg_host_exact_any.csv"))
write_csv(reg_multihost, file.path(tbl_dir, "table_03_reg_multi_host_exact.csv"))
write_csv(subpair_bg_jaccard, file.path(tbl_dir, "table_04_subpair_bg_jaccard.csv"))
write_csv(subpair_bg_exact, file.path(tbl_dir, "table_05_subpair_bg_exact.csv"))
write_csv(nn_cluster_stats, file.path(tbl_dir, "table_06_nn_cluster_stats.csv"))
write_csv(neuron_cluster_stats, file.path(tbl_dir, "table_07_neuron_cluster_stats.csv"))
write_csv(exact_guest_both_types, file.path(tbl_dir, "table_08_exact_guest_both_types.csv"))
write_csv(exact_guest_both_types_summary, file.path(tbl_dir, "table_09_exact_guest_both_types_summary.csv"))
write_csv(overall_shared_summary, file.path(tbl_dir, "table_10_exact_guest_both_types_overall.csv"))
write_csv(sel_jaccard_pairs, file.path(tbl_dir, "table_11_selected_pair_jaccard_examples.csv"))
write_csv(sel_exact_pairs, file.path(tbl_dir, "table_12_selected_pair_exact_examples.csv"))

# ---------------------------
# Figures
# ---------------------------
forest_guest <- reg_guest %>%
  dplyr::filter(term %in% c("log_nn_size", "mean_nnn_jaccard", "n_mixed_edges")) %>%
  dplyr::mutate(
    term_label = dplyr::recode(
      term,
      `log_nn_size` = "log(non-neuron cluster size)",
      `mean_nnn_jaccard` = "Mean NN-NN Jaccard",
      `n_mixed_edges` = "Number of neuron-NN partners"
    ),
    term_label = factor(term_label, levels = rev(term_label))
  )

fig1 <- ggplot(forest_guest, aes(x = OR, y = term_label)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.55, colour = "#7F7F7F") +
  geom_errorbarh(aes(xmin = OR_low, xmax = OR_high), height = 0.16, linewidth = 0.7, colour = "#4D4D4D") +
  geom_point(size = 2.9, colour = "#C44E52") +
  scale_x_log10(labels = label_number(accuracy = 0.01)) +
  labs(
    title = wrap_plot_title("Exact-contained non-neuron clusters are smaller and less self-cohesive within the non-neuronal network"),
    subtitle = wrap_plot_subtitle(
      paste0(
        "NN guest model: OR(log_nn_size) = ", sprintf("%.3f", reg_guest$OR[reg_guest$term == "log_nn_size"]),
        ", p = ", fmt_p_short(reg_guest$p_value[reg_guest$term == "log_nn_size"]),
        "; OR(mean_nnn_jaccard) = ", sprintf("%.4f", reg_guest$OR[reg_guest$term == "mean_nnn_jaccard"]),
        ", p = ", fmt_p_short(reg_guest$p_value[reg_guest$term == "mean_nnn_jaccard"])
      )
    ),
    x = "Adjusted odds ratio (log scale)",
    y = NULL
  ) +
  plot_theme
save_pub(fig1, "fig01_nn_guest_forest_topjournal", 8.8, 5.8)

nn_guest_plot <- nn_cluster_stats %>%
  dplyr::mutate(group = ifelse(guest_exact_any, "Exact-contained", "Not exact-contained"))

fig2 <- ggplot(nn_guest_plot, aes(group, nn_size, fill = group)) +
  geom_violin(width = 0.86, alpha = 0.18, colour = NA, trim = FALSE) +
  geom_boxplot(width = 0.28, outlier.shape = NA, linewidth = 0.7, alpha = 0.9) +
  scale_y_log10(labels = label_number(accuracy = 1)) +
  scale_fill_manual(values = col_guest) +
  labs(
    title = wrap_plot_title("Exact-contained non-neuron clusters are not random samples: they are markedly smaller"),
    subtitle = wrap_plot_subtitle(
      paste0(
        "Cluster-level model effect: OR(log_nn_size) = ",
        sprintf("%.3f", reg_guest$OR[reg_guest$term == "log_nn_size"]),
        "; p = ", fmt_p_short(reg_guest$p_value[reg_guest$term == "log_nn_size"])
      )
    ),
    x = NULL,
    y = "Non-neuron cluster size (log scale)"
  ) +
  plot_theme +
  theme(legend.position = "none")
save_pub(fig2, "fig02_nn_guest_size_topjournal", 8.2, 5.8)

fig3 <- ggplot(nn_guest_plot, aes(group, mean_nnn_jaccard, fill = group)) +
  geom_violin(width = 0.86, alpha = 0.18, colour = NA, trim = FALSE) +
  geom_boxplot(width = 0.28, outlier.shape = NA, linewidth = 0.7, alpha = 0.9) +
  scale_fill_manual(values = col_guest) +
  labs(
    title = wrap_plot_title("Exact-contained non-neuron clusters are less cohesive within the NN-NN interaction scaffold"),
    subtitle = wrap_plot_subtitle(
      paste0(
        "Cluster-level model effect: OR(mean_nnn_jaccard) = ",
        sprintf("%.4f", reg_guest$OR[reg_guest$term == "mean_nnn_jaccard"]),
        "; p = ", fmt_p_short(reg_guest$p_value[reg_guest$term == "mean_nnn_jaccard"])
      )
    ),
    x = NULL,
    y = "Mean NN-NN Jaccard"
  ) +
  plot_theme +
  theme(legend.position = "none")
save_pub(fig3, "fig03_nn_guest_nnn_cohesion_topjournal", 8.2, 5.8)

forest_host <- reg_host %>%
  dplyr::filter(term %in% c("neuron_typeGlut", "log_neuron_size", "n_mixed_edges")) %>%
  dplyr::mutate(
    term_label = dplyr::recode(
      term,
      `neuron_typeGlut` = "Glut vs Gaba",
      `log_neuron_size` = "log(neuron cluster size)",
      `n_mixed_edges` = "Mixed degree"
    ),
    term_label = factor(term_label, levels = rev(term_label))
  )

fig4 <- ggplot(forest_host, aes(x = OR, y = term_label)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.55, colour = "#7F7F7F") +
  geom_errorbarh(aes(xmin = OR_low, xmax = OR_high), height = 0.16, linewidth = 0.7, colour = "#4D4D4D") +
  geom_point(size = 2.9, colour = "#4E79A7") +
  scale_x_log10(labels = label_number(accuracy = 0.01)) +
  labs(
    title = wrap_plot_title("Glut clusters are more likely to become exact hosts of non-neuron islands"),
    subtitle = wrap_plot_subtitle(
      paste0(
        "After controlling for neuron size and mixed degree, OR(Glut host_exact_any) = ",
        sprintf("%.2f", reg_host$OR[reg_host$term == "neuron_typeGlut"]),
        "; p = ", fmt_p_short(reg_host$p_value[reg_host$term == "neuron_typeGlut"])
      )
    ),
    x = "Adjusted odds ratio (log scale)",
    y = NULL
  ) +
  plot_theme
save_pub(fig4, "fig04_glut_host_forest_topjournal", 8.6, 5.6)

shared_df <- tibble(
  group = c("Shared by both", "One neuron type only"),
  rate = c(
    mean(exact_guest_both_types$both_neuron_types_exact),
    1 - mean(exact_guest_both_types$both_neuron_types_exact)
  )
)

fig5 <- ggplot(shared_df, aes(group, rate, fill = group)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = percent(rate, accuracy = 0.1)), vjust = -0.25, size = 4.8) +
  scale_fill_manual(values = col_shared) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = wrap_plot_title("A subset of exact-contained non-neuron niches is shared by both Glut and Gaba hosts"),
    subtitle = wrap_plot_subtitle(
      paste0(
        "Among all exact guest NN clusters, ",
        percent(mean(exact_guest_both_types$both_neuron_types_exact), accuracy = 0.1),
        " have both Glut and Gaba exact hosts."
      )
    ),
    x = NULL,
    y = "Fraction of exact guest non-neuron clusters"
  ) +
  plot_theme +
  theme(legend.position = "none")
save_pub(fig5, "fig05_shared_exact_guest_rate_topjournal", 7.6, 5.6)

sel_jaccard_pairs <- sel_jaccard_pairs %>% dplyr::arrange(mean_diff_vs_bg)
fig6 <- ggplot(sel_jaccard_pairs, aes(x = mean_diff_vs_bg, y = forcats::fct_reorder(subpair_name_clean, mean_diff_vs_bg), fill = direction)) +
  geom_col() +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  scale_fill_manual(values = col_posneg) +
  labs(
    title = wrap_plot_title("What matters is not the coarse label, but the specific neuron × non-neuron pair"),
    subtitle = wrap_plot_subtitle("Representative pairs with Jaccard significantly above or below the within-neuron-type background."),
    x = "Δ mean Jaccard vs background",
    y = NULL
  ) +
  plot_theme
save_pub(fig6, "fig06_selected_pair_jaccard_examples_topjournal", 10.2, 6.4)

sel_exact_pairs <- sel_exact_pairs %>% dplyr::arrange(mean_diff_vs_bg)
fig7 <- ggplot(sel_exact_pairs, aes(x = mean_diff_vs_bg, y = forcats::fct_reorder(subpair_name_clean, mean_diff_vs_bg), fill = direction)) +
  geom_col() +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  scale_fill_manual(values = col_posneg) +
  labs(
    title = wrap_plot_title("The same non-neuron subclass can switch to a different spatial role under different neuron subclasses"),
    subtitle = wrap_plot_subtitle("Representative pairs with exact containment significantly above or below the within-neuron-type background."),
    x = "Δ exact containment rate vs background",
    y = NULL
  ) +
  plot_theme
save_pub(fig7, "fig07_selected_pair_exact_examples_topjournal", 10.2, 6.0)

shared_subclass_plot <- exact_guest_both_types_summary %>%
  dplyr::filter(n_exact_guest_clusters >= 10) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::arrange(both_neuron_types_rate)

fig8 <- ggplot(shared_subclass_plot, aes(x = both_neuron_types_rate, y = forcats::fct_reorder(nn_subclass, both_neuron_types_rate))) +
  geom_col(fill = "#7A4FB3") +
  labs(
    title = wrap_plot_title("Shared exact guest niches are concentrated in a subset of non-neuron subclasses"),
    subtitle = wrap_plot_subtitle("Top non-neuron subclasses ranked by the proportion of exact guest clusters shared by both Glut and Gaba hosts."),
    x = "Fraction shared by both neuron types",
    y = NULL
  ) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  plot_theme
save_pub(fig8, "fig08_shared_exact_guest_subclasses_topjournal", 8.6, 6.0)

figure_manifest <- tibble(
  filename = c(
    "fig01_nn_guest_forest_topjournal",
    "fig02_nn_guest_size_topjournal",
    "fig03_nn_guest_nnn_cohesion_topjournal",
    "fig04_glut_host_forest_topjournal",
    "fig05_shared_exact_guest_rate_topjournal",
    "fig06_selected_pair_jaccard_examples_topjournal",
    "fig07_selected_pair_exact_examples_topjournal",
    "fig08_shared_exact_guest_subclasses_topjournal"
  ),
  summary = c(
    "Cluster-level model for exact-contained non-neuron guests",
    "Exact-contained non-neuron guests are smaller",
    "Exact-contained non-neuron guests are less cohesive in NN-NN space",
    "Glut clusters are more likely to become exact hosts",
    "Overall fraction of exact guest NN clusters shared by both Glut and Gaba",
    "Selected pair-level Jaccard deviations from background",
    "Selected pair-level exact containment deviations from background",
    "Non-neuron subclasses enriched for cross-neuron-type exact sharing"
  )
)
write_csv(figure_manifest, file.path(out_dir, "figure_manifest.csv"))

message("All outputs written to: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
