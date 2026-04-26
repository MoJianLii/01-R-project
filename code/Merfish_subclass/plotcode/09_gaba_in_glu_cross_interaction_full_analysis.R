# 清空环境
rm(list = ls())
gc()
# ============================================================
# GABA fully contained in GLU versus other spatial interactions
# Cross-dataset comparative analysis for:
#   1) neuron-neuron interactions
#   2) neuron-nonneuron interactions
#
# Input files expected in the working directory:
#   neuron-neuron-partner.csv
#   neuron-nonneuron-partner.csv
#
# Main outputs:
#   gaba_in_glu_cross_results/
#   plus CSV tables and PNG figures
#
# Author: OpenAI GPT-5.4 Pro (equivalent R reproduction script)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(ggplot2)
})

options(stringsAsFactors = FALSE, scipen = 999)

# dplyr verb aliases to avoid conflicts with AnnotationDbi / other packages
select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter
mutate <- dplyr::mutate
summarise <- dplyr::summarise
arrange <- dplyr::arrange
distinct <- dplyr::distinct
count <- dplyr::count
left_join <- dplyr::left_join
right_join <- dplyr::right_join
full_join <- dplyr::full_join
inner_join <- dplyr::inner_join
pull <- dplyr::pull
transmute <- dplyr::transmute
first <- dplyr::first


# -----------------------------
# Paths
# -----------------------------
base_dir <- "E:/zaw/2603"

nn_path   <- file.path(base_dir, "neuron-neuron-partner.csv")
nonn_path <- file.path(base_dir, "neuron-nonneuron-partner.csv")

outdir <- file.path(base_dir, "R_output_09")
figdir <- file.path(outdir, "figures")
tabdir <- file.path(outdir, "tables")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(nn_path)) {
  stop("File not found: ", nn_path)
}
if (!file.exists(nonn_path)) {
  stop("File not found: ", nonn_path)
}

# -----------------------------
# Helpers
# -----------------------------
type_code <- c("Glut" = "E", "Gaba" = "I", "NonNeuron" = "N")

make_pair_id <- function(a, b) {
  ifelse(a <= b, paste(a, b, sep = "|||"), paste(b, a, sep = "|||"))
}

neglog10_safe <- function(x) {
  -log10(pmax(x, 1e-300))
}

fisher_compare <- function(a_full, a_total, b_full, b_total) {
  tab <- matrix(c(a_full, a_total - a_full, b_full, b_total - b_full),
                nrow = 2, byrow = TRUE)
  ft <- fisher.test(tab)
  tibble(
    OR = unname(ft$estimate),
    p = ft$p.value
  )
}

mwu_compare <- function(x, y) {
  wt <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
  tibble(
    x_median = median(x, na.rm = TRUE),
    y_median = median(y, na.rm = TRUE),
    p = wt$p.value
  )
}

signed_rank_compare <- function(a, b) {
  ok <- !(is.na(a) | is.na(b))
  a <- a[ok]
  b <- b[ok]
  if (length(a) == 0) {
    return(tibble(
      n = 0, A_median = NA_real_, B_median = NA_real_,
      median_diff = NA_real_, A_gt_B = NA_integer_,
      A_lt_B = NA_integer_, A_eq_B = NA_integer_, p = NA_real_
    ))
  }
  wt <- suppressWarnings(wilcox.test(a, b, paired = TRUE, exact = FALSE))
  tibble(
    n = length(a),
    A_median = median(a, na.rm = TRUE),
    B_median = median(b, na.rm = TRUE),
    median_diff = median(a - b, na.rm = TRUE),
    A_gt_B = sum(a > b, na.rm = TRUE),
    A_lt_B = sum(a < b, na.rm = TRUE),
    A_eq_B = sum(a == b, na.rm = TRUE),
    p = wt$p.value
  )
}

safe_mcnemar <- function(x, y) {
  tab <- table(factor(x, levels = c(FALSE, TRUE)),
               factor(y, levels = c(FALSE, TRUE)))
  mt <- mcnemar.test(tab, correct = FALSE)
  tibble(
    false_false = unname(tab[1, 1]),
    false_true  = unname(tab[1, 2]),
    true_false  = unname(tab[2, 1]),
    true_true   = unname(tab[2, 2]),
    p = mt$p.value
  )
}

prop_ge <- function(x, k) {
  tibble(
    yes = sum(x >= k, na.rm = TRUE),
    total = length(x),
    prop = mean(x >= k, na.rm = TRUE)
  )
}

pairwise_fullrate_fisher <- function(df, motif_a, motif_b) {
  a <- df %>% filter(motif == motif_a) %>%
    summarise(full_n = sum(full_containment), total_n = n())
  b <- df %>% filter(motif == motif_b) %>%
    summarise(full_n = sum(full_containment), total_n = n())
  tibble(
    motif_a = motif_a,
    motif_b = motif_b,
    a_rate = a$full_n / a$total_n,
    b_rate = b$full_n / b$total_n
  ) %>%
    bind_cols(fisher_compare(a$full_n, a$total_n, b$full_n, b$total_n))
}

bh <- function(p) p.adjust(p, method = "BH")

fmt_pct <- function(x, digits = 1) sprintf(paste0("%.", digits, "f%%"), x * 100)

# -----------------------------
# Load and prepare
# -----------------------------
nn   <- read_csv(nn_path, show_col_types = FALSE)
nonn <- read_csv(nonn_path, show_col_types = FALSE)

prepare_directed <- function(df, source_name) {
  tibble(
    source = source_name,
    guest_label = df$`cluster.1_label`,
    host_label = df$`cluster.2_label`,
    slide = df$`cluster.1_slide`,
    guest_layer = df$`cluster.1_layer`,
    host_layer = df$`cluster.2_layer`,
    guest_region = df$`cluster.1_region`,
    host_region = df$`cluster.2_region`,
    guest_type = df$`cluster.1_cell_Neruon_type`,
    host_type = df$`cluster.2_cell_Neruon_type`,
    guest_subclass = df$`cluster.1_subclass`,
    host_subclass = df$`cluster.2_subclass`,
    guest_size = df$`cluster.1_total_cell_num`,
    host_size = df$`cluster.2_total_cell_num`,
    guest_p = df$`cluster.1_cauchy_combination_p`,
    host_p = df$`cluster.2_cauchy_combination_p`,
    guest_ei = df$`cluster.1_E_I_Ratio`,
    host_ei = df$`cluster.2_E_I_Ratio`,
    overlap_cell = df$overlap_cell,
    union_cell = df$union_cell,
    jaccard = df$jaccard,
    guest_overlap_pct = df$`cluster.1.overlap.percent`,
    host_overlap_pct = df$`cluster.2.overlap.percent`
  ) %>%
    mutate(
      full_containment = guest_overlap_pct == 1,
      partial = guest_overlap_pct < 1 & host_overlap_pct < 1,
      guest_fill = guest_overlap_pct,
      host_fill = host_overlap_pct,
      size_ratio_host_over_guest = host_size / guest_size,
      log_size_ratio = log2(size_ratio_host_over_guest),
      guest_logp = neglog10_safe(guest_p),
      host_logp = neglog10_safe(host_p),
      logp_gap = host_logp - guest_logp,
      overlap_asym = guest_overlap_pct - host_overlap_pct,
      motif = paste0(type_code[guest_type], "\u2282", type_code[host_type]),
      pair_id = make_pair_id(guest_label, host_label)
    )
}

combined <- bind_rows(
  prepare_directed(nn, "neuron-neuron"),
  prepare_directed(nonn, "neuron-nonneuron")
)

full    <- combined %>% filter(full_containment)
partial <- combined %>% filter(partial)
eligible <- combined %>% filter(host_size >= guest_size)

motif_order <- c("I\u2282E", "E\u2282I", "I\u2282N", "N\u2282I",
                 "E\u2282N", "N\u2282E", "I\u2282I", "E\u2282E", "N\u2282N")

clusters_unique <- bind_rows(
  combined %>% transmute(label = guest_label, type = guest_type, subclass = guest_subclass, layer = guest_layer, slide = slide),
  combined %>% transmute(label = host_label, type = host_type, subclass = host_subclass, layer = host_layer, slide = slide)
) %>% distinct(label, .keep_all = TRUE)

# ============================================================
# Table 01: full containment rates by motif
# ============================================================
table01 <- combined %>%
  group_by(motif) %>%
  summarise(full_n = sum(full_containment), total_n = n(), full_rate = full_n / total_n, .groups = "drop") %>%
  mutate(motif = factor(motif, levels = motif_order)) %>%
  arrange(motif)

write_csv(table01, file.path(tabdir, "table01_full_containment_rates_by_motif.csv"))

# ============================================================
# Table 02: size-eligible rates by motif
# ============================================================
table02 <- eligible %>%
  group_by(motif) %>%
  summarise(full_n = sum(full_containment), total_n = n(),
            full_rate_if_host_ge_guest = full_n / total_n, .groups = "drop") %>%
  mutate(motif = factor(motif, levels = motif_order)) %>%
  arrange(motif)

write_csv(table02, file.path(tabdir, "table02_size_eligible_full_rates_by_motif.csv"))

# ============================================================
# Table 03: I⊂E versus selected motif rate tests (raw + size-eligible)
# ============================================================
selected_compare <- c("E\u2282I", "I\u2282N", "N\u2282I", "E\u2282N", "N\u2282E", "I\u2282I", "E\u2282E", "N\u2282N")

get_rate_row <- function(df_all, df_eligible, motif_b) {
  a <- df_all %>% filter(motif == "I\u2282E") %>% summarise(full_n = sum(full_containment), total_n = n())
  b <- df_all %>% filter(motif == motif_b)      %>% summarise(full_n = sum(full_containment), total_n = n())
  ae <- df_eligible %>% filter(motif == "I\u2282E") %>% summarise(full_n = sum(full_containment), total_n = n())
  be <- df_eligible %>% filter(motif == motif_b)      %>% summarise(full_n = sum(full_containment), total_n = n())
  
  raw <- fisher_compare(a$full_n, a$total_n, b$full_n, b$total_n)
  elig <- fisher_compare(ae$full_n, ae$total_n, be$full_n, be$total_n)
  
  tibble(
    compare_to = motif_b,
    I_in_E_rate_raw = a$full_n / a$total_n,
    compare_rate_raw = b$full_n / b$total_n,
    OR_raw = raw$OR,
    p_raw = raw$p,
    I_in_E_rate_if_host_ge_guest = ae$full_n / ae$total_n,
    compare_rate_if_host_ge_guest = be$full_n / be$total_n,
    OR_if_host_ge_guest = elig$OR,
    p_if_host_ge_guest = elig$p
  )
}

table03 <- map_dfr(selected_compare, ~get_rate_row(combined, eligible, .x)) %>%
  mutate(
    fdr_raw = bh(p_raw),
    fdr_if_host_ge_guest = bh(p_if_host_ge_guest)
  )

write_csv(table03, file.path(tabdir, "table03_i_in_e_vs_selected_rate_tests.csv"))

# ============================================================
# Table 04: geometry summary for all full-containment motifs
# ============================================================
table04 <- full %>%
  group_by(motif) %>%
  summarise(
    full_n = n(),
    guest_size_median = median(guest_size),
    host_size_median = median(host_size),
    host_to_guest_size_ratio_median = median(size_ratio_host_over_guest),
    overlap_cell_median = median(overlap_cell),
    jaccard_median = median(jaccard),
    host_fill_median = median(host_fill),
    host_log10p_median = median(host_logp),
    guest_log10p_median = median(guest_logp),
    log10p_gap_median = median(logp_gap),
    .groups = "drop"
  ) %>%
  mutate(motif = factor(motif, levels = motif_order)) %>%
  arrange(motif)

write_csv(table04, file.path(tabdir, "table04_full_containment_geometry_summary.csv"))

# ============================================================
# Table 05: I⊂E versus selected motifs on key continuous metrics
# ============================================================
metrics_for_test <- c("guest_size", "host_size", "size_ratio_host_over_guest",
                      "overlap_cell", "jaccard", "host_fill",
                      "host_logp", "guest_logp", "logp_gap")

table05 <- purrr::map_dfr(selected_compare, function(compare_to_now) {
  purrr::map_dfr(metrics_for_test, function(metric_now) {
    x <- full %>%
      dplyr::filter(motif == "I\u2282E") %>%
      dplyr::pull(!!rlang::sym(metric_now))
    
    y <- full %>%
      dplyr::filter(motif == compare_to_now) %>%
      dplyr::pull(!!rlang::sym(metric_now))
    
    tibble(
      compare_to = compare_to_now,
      metric = metric_now,
      I_in_E_median = median(x, na.rm = TRUE),
      compare_median = median(y, na.rm = TRUE),
      p = suppressWarnings(wilcox.test(x, y, exact = FALSE)$p.value)
    )
  })
}) %>%
  dplyr::mutate(fdr = bh(p)) %>%
  dplyr::select(compare_to, metric, I_in_E_median, compare_median, p, fdr)

write_csv(table05, file.path(tabdir, "table05_i_in_e_vs_selected_metric_tests.csv"))

# ============================================================
# Table 06/07: partial asymmetry
# ============================================================
table06 <- partial %>%
  group_by(motif) %>%
  summarise(
    partial_n = n(),
    asym_median = median(overlap_asym),
    jaccard_median = median(jaccard),
    .groups = "drop"
  ) %>%
  mutate(motif = factor(motif, levels = motif_order)) %>%
  arrange(motif)

write_csv(table06, file.path(tabdir, "table06_partial_asymmetry_by_motif.csv"))

partial_compare <- c("I\u2282N", "I\u2282I", "N\u2282E", "E\u2282I", "N\u2282I", "E\u2282N")
table07 <- map_dfr(partial_compare, function(m) {
  x <- partial %>% filter(motif == "I\u2282E") %>% pull(overlap_asym)
  y <- partial %>% filter(motif == m) %>% pull(overlap_asym)
  tibble(
    compare_to = m,
    I_in_E_asym_median = median(x, na.rm = TRUE),
    compare_asym_median = median(y, na.rm = TRUE),
    p = suppressWarnings(wilcox.test(x, y, exact = FALSE)$p.value)
  )
}) %>%
  mutate(fdr = bh(p))

write_csv(table07, file.path(tabdir, "table07_i_in_e_partial_asymmetry_tests.csv"))

# ============================================================
# Table 08: slide-level rate comparisons
# ============================================================
slide_rates <- combined %>%
  group_by(slide, motif) %>%
  summarise(full_rate = mean(full_containment), .groups = "drop") %>%
  pivot_wider(names_from = motif, values_from = full_rate)

slide_comp_pairs <- tribble(
  ~A, ~B,
  "I\u2282E", "E\u2282I",
  "I\u2282E", "I\u2282N",
  "I\u2282E", "N\u2282I",
  "I\u2282E", "N\u2282E",
  "I\u2282E", "E\u2282E"
)

table08 <- pmap_dfr(slide_comp_pairs, function(A, B) {
  sub <- slide_rates %>%
    select(all_of(c(A, B))) %>%
    drop_na()
  tibble(
    A = A,
    B = B,
    n_slides = nrow(sub),
    A_rate_median = median(sub[[A]]),
    B_rate_median = median(sub[[B]]),
    A_gt_B_slides = sum(sub[[A]] > sub[[B]]),
    A_lt_B_slides = sum(sub[[A]] < sub[[B]]),
    p = suppressWarnings(wilcox.test(sub[[A]], sub[[B]], paired = TRUE, exact = FALSE)$p.value)
  )
})

write_csv(table08, file.path(tabdir, "table08_slide_level_rate_comparisons.csv"))

# ============================================================
# Same-GABA guest matched analyses
# ============================================================
gaba_full <- full %>% filter(guest_type == "Gaba")

guest_host_agg <- gaba_full %>%
  group_by(guest_label, host_type) %>%
  summarise(
    n_hosts = n_distinct(host_label),
    host_size_med = median(host_size),
    host_logp_med = median(host_logp),
    host_fill_med = median(host_fill),
    jaccard_med = median(jaccard),
    size_ratio_med = median(size_ratio_host_over_guest),
    overlap_med = median(overlap_cell),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = host_type,
              values_from = c(n_hosts, host_size_med, host_logp_med, host_fill_med, jaccard_med, size_ratio_med, overlap_med))

table09_rows <- list(
  tibble(comparison = "same_GABA_guest_Glut_host_vs_NonNeuron_host", metric = "host_size") %>%
    bind_cols(signed_rank_compare(guest_host_agg$host_size_med_Glut, guest_host_agg$host_size_med_NonNeuron)),
  tibble(comparison = "same_GABA_guest_Glut_host_vs_NonNeuron_host", metric = "host_log10p") %>%
    bind_cols(signed_rank_compare(guest_host_agg$host_logp_med_Glut, guest_host_agg$host_logp_med_NonNeuron)),
  tibble(comparison = "same_GABA_guest_Glut_host_vs_NonNeuron_host", metric = "host_fill") %>%
    bind_cols(signed_rank_compare(guest_host_agg$host_fill_med_Glut, guest_host_agg$host_fill_med_NonNeuron)),
  tibble(comparison = "same_GABA_guest_Glut_host_vs_NonNeuron_host", metric = "host_to_guest_size_ratio") %>%
    bind_cols(signed_rank_compare(guest_host_agg$size_ratio_med_Glut, guest_host_agg$size_ratio_med_NonNeuron)),
  tibble(comparison = "same_GABA_guest_Glut_host_vs_NonNeuron_host", metric = "n_hosts") %>%
    bind_cols(signed_rank_compare(guest_host_agg$n_hosts_Glut, guest_host_agg$n_hosts_NonNeuron)),
  
  tibble(comparison = "same_GABA_guest_Glut_host_vs_Gaba_host", metric = "host_size") %>%
    bind_cols(signed_rank_compare(guest_host_agg$host_size_med_Glut, guest_host_agg$host_size_med_Gaba)),
  tibble(comparison = "same_GABA_guest_Glut_host_vs_Gaba_host", metric = "host_log10p") %>%
    bind_cols(signed_rank_compare(guest_host_agg$host_logp_med_Glut, guest_host_agg$host_logp_med_Gaba)),
  tibble(comparison = "same_GABA_guest_Glut_host_vs_Gaba_host", metric = "host_fill") %>%
    bind_cols(signed_rank_compare(guest_host_agg$host_fill_med_Glut, guest_host_agg$host_fill_med_Gaba)),
  tibble(comparison = "same_GABA_guest_Glut_host_vs_Gaba_host", metric = "host_to_guest_size_ratio") %>%
    bind_cols(signed_rank_compare(guest_host_agg$size_ratio_med_Glut, guest_host_agg$size_ratio_med_Gaba)),
  tibble(comparison = "same_GABA_guest_Glut_host_vs_Gaba_host", metric = "n_hosts") %>%
    bind_cols(signed_rank_compare(guest_host_agg$n_hosts_Glut, guest_host_agg$n_hosts_Gaba))
)

table09 <- bind_rows(table09_rows)
write_csv(table09, file.path(tabdir, "table09_matched_same_gaba_guest_comparisons.csv"))

# ============================================================
# Same-GLU host matched analyses
# ============================================================
glu_host_full <- full %>% filter(host_type == "Glut")

host_guest_agg <- glu_host_full %>%
  group_by(host_label, guest_type) %>%
  summarise(
    n_guests = n_distinct(guest_label),
    guest_size_med = median(guest_size),
    guest_logp_med = median(guest_logp),
    host_fill_med = median(host_fill),
    jaccard_med = median(jaccard),
    size_ratio_med = median(size_ratio_host_over_guest),
    overlap_med = median(overlap_cell),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = guest_type,
              values_from = c(n_guests, guest_size_med, guest_logp_med, host_fill_med, jaccard_med, size_ratio_med, overlap_med))

table10_rows <- list(
  tibble(comparison = "same_Glut_host_Gaba_guest_vs_NonNeuron_guest", metric = "guest_size") %>%
    bind_cols(signed_rank_compare(host_guest_agg$guest_size_med_Gaba, host_guest_agg$guest_size_med_NonNeuron)),
  tibble(comparison = "same_Glut_host_Gaba_guest_vs_NonNeuron_guest", metric = "guest_log10p") %>%
    bind_cols(signed_rank_compare(host_guest_agg$guest_logp_med_Gaba, host_guest_agg$guest_logp_med_NonNeuron)),
  tibble(comparison = "same_Glut_host_Gaba_guest_vs_NonNeuron_guest", metric = "host_fill") %>%
    bind_cols(signed_rank_compare(host_guest_agg$host_fill_med_Gaba, host_guest_agg$host_fill_med_NonNeuron)),
  tibble(comparison = "same_Glut_host_Gaba_guest_vs_NonNeuron_guest", metric = "host_to_guest_size_ratio") %>%
    bind_cols(signed_rank_compare(host_guest_agg$size_ratio_med_Gaba, host_guest_agg$size_ratio_med_NonNeuron)),
  tibble(comparison = "same_Glut_host_Gaba_guest_vs_NonNeuron_guest", metric = "n_guests") %>%
    bind_cols(signed_rank_compare(host_guest_agg$n_guests_Gaba, host_guest_agg$n_guests_NonNeuron)),
  
  tibble(comparison = "same_Glut_host_Gaba_guest_vs_Glut_guest", metric = "guest_size") %>%
    bind_cols(signed_rank_compare(host_guest_agg$guest_size_med_Gaba, host_guest_agg$guest_size_med_Glut)),
  tibble(comparison = "same_Glut_host_Gaba_guest_vs_Glut_guest", metric = "guest_log10p") %>%
    bind_cols(signed_rank_compare(host_guest_agg$guest_logp_med_Gaba, host_guest_agg$guest_logp_med_Glut)),
  tibble(comparison = "same_Glut_host_Gaba_guest_vs_Glut_guest", metric = "host_fill") %>%
    bind_cols(signed_rank_compare(host_guest_agg$host_fill_med_Gaba, host_guest_agg$host_fill_med_Glut)),
  tibble(comparison = "same_Glut_host_Gaba_guest_vs_Glut_guest", metric = "host_to_guest_size_ratio") %>%
    bind_cols(signed_rank_compare(host_guest_agg$size_ratio_med_Gaba, host_guest_agg$size_ratio_med_Glut)),
  tibble(comparison = "same_Glut_host_Gaba_guest_vs_Glut_guest", metric = "n_guests") %>%
    bind_cols(signed_rank_compare(host_guest_agg$n_guests_Gaba, host_guest_agg$n_guests_Glut))
)

table10 <- bind_rows(table10_rows)
write_csv(table10, file.path(tabdir, "table10_matched_same_glu_host_comparisons.csv"))

# ============================================================
# Multiplicity
# ============================================================
guest_mult <- full %>%
  group_by(motif, guest_label) %>%
  summarise(n_hosts = n_distinct(host_label), .groups = "drop") %>%
  group_by(motif) %>%
  summarise(
    guest_n = n(),
    median_hosts_per_guest = median(n_hosts),
    mean_hosts_per_guest = mean(n_hosts),
    pct_ge2_hosts = mean(n_hosts >= 2),
    pct_ge3_hosts = mean(n_hosts >= 3),
    .groups = "drop"
  )

host_mult <- full %>%
  group_by(motif, host_label) %>%
  summarise(n_guests = n_distinct(guest_label), .groups = "drop") %>%
  group_by(motif) %>%
  summarise(
    host_n = n(),
    median_guests_per_host = median(n_guests),
    mean_guests_per_host = mean(n_guests),
    pct_ge2_guests = mean(n_guests >= 2),
    pct_ge3_guests = mean(n_guests >= 3),
    .groups = "drop"
  )

table11 <- guest_mult %>%
  left_join(host_mult, by = "motif") %>%
  mutate(motif = factor(motif, levels = motif_order)) %>%
  arrange(motif)

write_csv(table11, file.path(tabdir, "table11_multiplicity_summary_by_motif.csv"))

get_multiplicity_test <- function(kind = c("hosts_with_>=2_guests", "guests_with_>=2_hosts"),
                                  motif_b) {
  kind <- match.arg(kind)
  if (kind == "hosts_with_>=2_guests") {
    a_raw <- full %>% filter(motif == "I\u2282E") %>% group_by(host_label) %>% summarise(v = n_distinct(guest_label), .groups = "drop") %>% pull(v)
    b_raw <- full %>% filter(motif == motif_b)     %>% group_by(host_label) %>% summarise(v = n_distinct(guest_label), .groups = "drop") %>% pull(v)
  } else {
    a_raw <- full %>% filter(motif == "I\u2282E") %>% group_by(guest_label) %>% summarise(v = n_distinct(host_label), .groups = "drop") %>% pull(v)
    b_raw <- full %>% filter(motif == motif_b)     %>% group_by(guest_label) %>% summarise(v = n_distinct(host_label), .groups = "drop") %>% pull(v)
  }
  
  a <- prop_ge(a_raw, 2)
  b <- prop_ge(b_raw, 2)
  ft <- fisher_compare(a$yes, a$total, b$yes, b$total)
  
  tibble(
    comparison = paste(kind, "I\u2282E vs", motif_b),
    I_in_E_pct = a$prop,
    compare_pct = b$prop,
    OR = ft$OR,
    p = ft$p
  )
}

mult_compare <- c("I\u2282N", "I\u2282I", "E\u2282I", "E\u2282E", "N\u2282E")
table12 <- bind_rows(
  map_dfr(mult_compare, ~get_multiplicity_test("hosts_with_>=2_guests", .x)),
  map_dfr(mult_compare, ~get_multiplicity_test("guests_with_>=2_hosts", .x))
) %>% mutate(fdr = bh(p))

write_csv(table12, file.path(tabdir, "table12_multiplicity_tests.csv"))

# ============================================================
# Host-host overlap around multi-host guests
# ============================================================
nn_lookup <- combined %>%
  filter(source == "neuron-neuron") %>%
  group_by(pair_id) %>%
  summarise(jaccard = max(jaccard), .groups = "drop")

nonn_lookup <- combined %>%
  filter(source == "neuron-nonneuron") %>%
  group_by(pair_id) %>%
  summarise(jaccard = max(jaccard), .groups = "drop")

gg_background <- combined %>%
  filter(guest_type == "Glut", host_type == "Glut") %>%
  distinct(pair_id, jaccard) %>%
  pull(jaccard)

ii_background <- combined %>%
  filter(guest_type == "Gaba", host_type == "Gaba") %>%
  distinct(pair_id, jaccard) %>%
  pull(jaccard)

nn_background <- combined %>%
  filter(guest_type == "NonNeuron", host_type == "NonNeuron") %>%
  distinct(pair_id, jaccard) %>%
  pull(jaccard)

get_hosthost_df <- function(motif_name) {
  d <- full %>% filter(motif == motif_name)
  if (nrow(d) == 0) return(tibble(guest_label = character(), pair_id = character(), jaccard = numeric()))
  host_type_here <- unique(d$host_type)
  lookup <- if (host_type_here == "NonNeuron") nonn_lookup else nn_lookup
  
  out <- d %>%
    group_by(guest_label) %>%
    summarise(hosts = list(unique(host_label)), .groups = "drop") %>%
    mutate(
      pair_id = map(hosts, function(h) {
        if (length(h) < 2) return(character())
        combn(h, 2, FUN = function(z) make_pair_id(z[1], z[2]), simplify = TRUE)
      })
    ) %>%
    dplyr::select(guest_label, pair_id) %>%
    tidyr::unnest(pair_id)
  
  out %>% left_join(lookup, by = "pair_id")
}

hosthost_IE <- get_hosthost_df("I\u2282E")
hosthost_IN <- get_hosthost_df("I\u2282N")
hosthost_II <- get_hosthost_df("I\u2282I")
hosthost_NE <- get_hosthost_df("N\u2282E")
hosthost_EE <- get_hosthost_df("E\u2282E")

compare_hosthost <- function(x, bg, motif_name, bg_name) {
  wt <- suppressWarnings(wilcox.test(x$jaccard, bg, exact = FALSE))
  tibble(
    motif = motif_name,
    hosthost_pair_n = sum(!is.na(x$jaccard)),
    hosthost_jaccard_median = median(x$jaccard, na.rm = TRUE),
    background = bg_name,
    background_n = length(bg),
    background_jaccard_median = median(bg, na.rm = TRUE),
    p_vs_background = wt$p.value
  )
}

table13 <- bind_rows(
  compare_hosthost(hosthost_IE, gg_background, "I\u2282E", "GG_background"),
  compare_hosthost(hosthost_IN, nn_background, "I\u2282N", "NN_background"),
  compare_hosthost(hosthost_II, ii_background, "I\u2282I", "II_background"),
  compare_hosthost(hosthost_NE, gg_background, "N\u2282E", "GG_background"),
  compare_hosthost(hosthost_EE, gg_background, "E\u2282E", "GG_background"),
  
  tibble(
    motif = "I\u2282E_vs_I\u2282N",
    hosthost_pair_n = sum(!is.na(hosthost_IE$jaccard)),
    hosthost_jaccard_median = median(hosthost_IE$jaccard, na.rm = TRUE),
    background = "I\u2282N",
    background_n = sum(!is.na(hosthost_IN$jaccard)),
    background_jaccard_median = median(hosthost_IN$jaccard, na.rm = TRUE),
    p_vs_background = suppressWarnings(wilcox.test(hosthost_IE$jaccard, hosthost_IN$jaccard, exact = FALSE)$p.value)
  ),
  tibble(
    motif = "I\u2282E_vs_N\u2282E",
    hosthost_pair_n = sum(!is.na(hosthost_IE$jaccard)),
    hosthost_jaccard_median = median(hosthost_IE$jaccard, na.rm = TRUE),
    background = "N\u2282E",
    background_n = sum(!is.na(hosthost_NE$jaccard)),
    background_jaccard_median = median(hosthost_NE$jaccard, na.rm = TRUE),
    p_vs_background = suppressWarnings(wilcox.test(hosthost_IE$jaccard, hosthost_NE$jaccard, exact = FALSE)$p.value)
  ),
  tibble(
    motif = "I\u2282E_vs_E\u2282E",
    hosthost_pair_n = sum(!is.na(hosthost_IE$jaccard)),
    hosthost_jaccard_median = median(hosthost_IE$jaccard, na.rm = TRUE),
    background = "E\u2282E",
    background_n = sum(!is.na(hosthost_EE$jaccard)),
    background_jaccard_median = median(hosthost_EE$jaccard, na.rm = TRUE),
    p_vs_background = suppressWarnings(wilcox.test(hosthost_IE$jaccard, hosthost_EE$jaccard, exact = FALSE)$p.value)
  )
)

write_csv(table13, file.path(tabdir, "table13_hosthost_overlap_tests.csv"))

# ============================================================
# GABA cluster-level preference for GLU versus NonNeuron hosts
# ============================================================
gaba_clusters <- clusters_unique %>% filter(type == "Gaba")

gaba_guest_full <- full %>%
  filter(guest_type == "Gaba") %>%
  count(guest_label, host_type, name = "n") %>%
  pivot_wider(names_from = host_type, values_from = n, values_fill = 0)

gaba_status <- gaba_clusters %>%
  transmute(guest_label = label, guest_subclass = subclass, guest_layer = layer, slide = slide) %>%
  left_join(gaba_guest_full, by = "guest_label") %>%
  mutate(
    Gaba = coalesce(Gaba, 0L),
    Glut = coalesce(Glut, 0L),
    NonNeuron = coalesce(NonNeuron, 0L),
    any_Glut_host = Glut > 0,
    any_NonNeuron_host = NonNeuron > 0,
    any_Gaba_host = Gaba > 0
  )

mcn_gaba <- safe_mcnemar(gaba_status$any_Glut_host, gaba_status$any_NonNeuron_host)

table14 <- tibble(
  n_gaba_clusters = nrow(gaba_status),
  any_Glut_host_n = sum(gaba_status$any_Glut_host),
  any_NonNeuron_host_n = sum(gaba_status$any_NonNeuron_host),
  both_n = sum(gaba_status$any_Glut_host & gaba_status$any_NonNeuron_host),
  Glut_only_n = sum(gaba_status$any_Glut_host & !gaba_status$any_NonNeuron_host),
  NonNeuron_only_n = sum(!gaba_status$any_Glut_host & gaba_status$any_NonNeuron_host),
  McNemar_p = mcn_gaba$p
)

write_csv(table14, file.path(tabdir, "table14_gaba_cluster_glut_vs_nonneuron_host_preference.csv"))

table15 <- gaba_status %>%
  group_by(guest_layer) %>%
  summarise(
    n_gaba_clusters = n(),
    pct_any_Glut_host = mean(any_Glut_host),
    pct_any_NN_host = mean(any_NonNeuron_host),
    glut_only = sum(any_Glut_host & !any_NonNeuron_host),
    nn_only = sum(!any_Glut_host & any_NonNeuron_host),
    mcnemar_p = safe_mcnemar(any_Glut_host, any_NonNeuron_host)$p,
    .groups = "drop"
  ) %>%
  mutate(fdr = bh(mcnemar_p))

write_csv(table15, file.path(tabdir, "table15_gaba_layer_glut_vs_nonneuron_host_preference.csv"))

table16 <- gaba_status %>%
  group_by(guest_subclass) %>%
  summarise(
    n_clusters = n(),
    pct_any_Glut = mean(any_Glut_host),
    pct_any_NonNeuron = mean(any_NonNeuron_host),
    delta = pct_any_Glut - pct_any_NonNeuron,
    both = mean(any_Glut_host & any_NonNeuron_host),
    glut_only = sum(any_Glut_host & !any_NonNeuron_host),
    nn_only = sum(!any_Glut_host & any_NonNeuron_host),
    mcnemar_p = safe_mcnemar(any_Glut_host, any_NonNeuron_host)$p,
    .groups = "drop"
  ) %>%
  mutate(fdr = bh(mcnemar_p)) %>%
  arrange(desc(delta))

write_csv(table16, file.path(tabdir, "table16_gaba_subclass_glut_vs_nonneuron_host_preference.csv"))

# ============================================================
# GLU host-level preference for GABA versus NonNeuron guests
# ============================================================
glu_clusters <- clusters_unique %>% filter(type == "Glut")

glu_host_full_counts <- full %>%
  filter(host_type == "Glut") %>%
  count(host_label, guest_type, name = "n") %>%
  pivot_wider(names_from = guest_type, values_from = n, values_fill = 0)

glu_status <- glu_clusters %>%
  transmute(host_label = label, host_subclass = subclass, host_layer = layer, slide = slide) %>%
  left_join(glu_host_full_counts, by = "host_label") %>%
  mutate(
    Gaba = coalesce(Gaba, 0L),
    Glut = coalesce(Glut, 0L),
    NonNeuron = coalesce(NonNeuron, 0L),
    any_Gaba_guest = Gaba > 0,
    any_NonNeuron_guest = NonNeuron > 0,
    any_Glut_guest = Glut > 0
  )

table17 <- tibble(
  n_glu_clusters = nrow(glu_status),
  any_Gaba_guest_n = sum(glu_status$any_Gaba_guest),
  any_NonNeuron_guest_n = sum(glu_status$any_NonNeuron_guest),
  both_n = sum(glu_status$any_Gaba_guest & glu_status$any_NonNeuron_guest),
  Gaba_only_n = sum(glu_status$any_Gaba_guest & !glu_status$any_NonNeuron_guest),
  NonNeuron_only_n = sum(!glu_status$any_Gaba_guest & glu_status$any_NonNeuron_guest),
  McNemar_p = safe_mcnemar(glu_status$any_Gaba_guest, glu_status$any_NonNeuron_guest)$p
)

write_csv(table17, file.path(tabdir, "table17_glu_cluster_gaba_vs_nonneuron_guest_preference.csv"))

table18 <- glu_status %>%
  group_by(host_layer) %>%
  summarise(
    n_glu_clusters = n(),
    pct_any_Gaba_guest = mean(any_Gaba_guest),
    pct_any_NonNeuron_guest = mean(any_NonNeuron_guest),
    gaba_only = sum(any_Gaba_guest & !any_NonNeuron_guest),
    nn_only = sum(!any_Gaba_guest & any_NonNeuron_guest),
    mcnemar_p = safe_mcnemar(any_Gaba_guest, any_NonNeuron_guest)$p,
    .groups = "drop"
  ) %>%
  mutate(fdr = bh(mcnemar_p))

write_csv(table18, file.path(tabdir, "table18_glu_layer_gaba_vs_nonneuron_guest_preference.csv"))

table19 <- glu_status %>%
  group_by(host_subclass) %>%
  summarise(
    n_clusters = n(),
    pct_any_Gaba_guest = mean(any_Gaba_guest),
    pct_any_NonNeuron_guest = mean(any_NonNeuron_guest),
    delta = pct_any_Gaba_guest - pct_any_NonNeuron_guest,
    both = mean(any_Gaba_guest & any_NonNeuron_guest),
    gaba_only = sum(any_Gaba_guest & !any_NonNeuron_guest),
    nn_only = sum(!any_Gaba_guest & any_NonNeuron_guest),
    mcnemar_p = safe_mcnemar(any_Gaba_guest, any_NonNeuron_guest)$p,
    .groups = "drop"
  ) %>%
  mutate(fdr = bh(mcnemar_p)) %>%
  arrange(desc(delta))

write_csv(table19, file.path(tabdir, "table19_glu_subclass_gaba_vs_nonneuron_guest_preference.csv"))

# ============================================================
# Enrichment analyses
# ============================================================
# Table 20: GABA subclass enrichment for any GLU host
table20 <- gaba_status %>%
  group_by(guest_subclass) %>%
  group_modify(~{
    if (nrow(.x) < 10) return(tibble())
    a <- sum(.x$any_Glut_host)
    b <- nrow(.x) - a
    rest <- gaba_status %>% filter(guest_subclass != .y$guest_subclass)
    c <- sum(rest$any_Glut_host)
    d <- nrow(rest) - c
    ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))
    tibble(
      n_clusters = nrow(.x),
      pct_any_Glut_host = a / nrow(.x),
      odds_ratio = unname(ft$estimate),
      p = ft$p.value
    )
  }) %>%
  ungroup() %>%
  mutate(fdr = bh(p)) %>%
  arrange(desc(odds_ratio))

write_csv(table20, file.path(tabdir, "table20_gaba_subclass_enrichment_for_any_glut_host.csv"))

# Table 21: GLU host subclass enrichment for any GABA guest
table21 <- glu_status %>%
  group_by(host_subclass) %>%
  group_modify(~{
    if (nrow(.x) < 10) return(tibble())
    a <- sum(.x$any_Gaba_guest)
    b <- nrow(.x) - a
    rest <- glu_status %>% filter(host_subclass != .y$host_subclass)
    c <- sum(rest$any_Gaba_guest)
    d <- nrow(rest) - c
    ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))
    tibble(
      n_clusters = nrow(.x),
      pct_any_Gaba_guest = a / nrow(.x),
      odds_ratio = unname(ft$estimate),
      p = ft$p.value
    )
  }) %>%
  ungroup() %>%
  mutate(fdr = bh(p)) %>%
  arrange(desc(odds_ratio))

write_csv(table21, file.path(tabdir, "table21_glu_host_subclass_enrichment_for_any_gaba_guest.csv"))

# Table 22: pair-specific full enrichment within I-E overlaps
ie <- combined %>% filter(guest_type == "Gaba", host_type == "Glut") %>%
  mutate(is_full = full_containment)

table22 <- ie %>%
  count(guest_subclass, host_subclass, wt = is_full, name = "full_n") %>%
  left_join(ie %>% count(guest_subclass, host_subclass, name = "total_n"),
            by = c("guest_subclass", "host_subclass")) %>%
  filter(total_n >= 20, full_n >= 5) %>%
  rowwise() %>%
  mutate(
    rest_full = sum(ie$is_full) - full_n,
    rest_total = nrow(ie) - total_n,
    ft = list(fisher.test(matrix(c(full_n, total_n - full_n,
                                   rest_full, rest_total - rest_full),
                                 nrow = 2, byrow = TRUE))),
    full_rate = full_n / total_n,
    odds_ratio = unname(ft$estimate),
    p = ft$p.value
  ) %>%
  ungroup() %>%
  select(guest_subclass, host_subclass, full_n, total_n, full_rate, odds_ratio, p) %>%
  mutate(fdr = bh(p)) %>%
  arrange(fdr, desc(odds_ratio))

write_csv(table22, file.path(tabdir, "table22_i_in_e_pair_specific_full_enrichment.csv"))

# Table 23: host subclass enrichment within I-N overlaps
in_df <- combined %>% filter(guest_type == "Gaba", host_type == "NonNeuron") %>%
  mutate(is_full = full_containment)

table23 <- in_df %>%
  group_by(host_subclass) %>%
  summarise(
    full_n = sum(is_full),
    total_n = n(),
    .groups = "drop"
  ) %>%
  filter(total_n >= 20, full_n >= 5) %>%
  rowwise() %>%
  mutate(
    rest_full = sum(in_df$is_full) - full_n,
    rest_total = nrow(in_df) - total_n,
    ft = list(fisher.test(matrix(c(full_n, total_n - full_n,
                                   rest_full, rest_total - rest_full),
                                 nrow = 2, byrow = TRUE))),
    full_rate = full_n / total_n,
    odds_ratio = unname(ft$estimate),
    p = ft$p.value
  ) %>%
  ungroup() %>%
  select(host_subclass, full_n, total_n, full_rate, odds_ratio, p) %>%
  mutate(fdr = bh(p)) %>%
  arrange(fdr, desc(odds_ratio))

write_csv(table23, file.path(tabdir, "table23_i_in_n_host_subclass_full_enrichment.csv"))

# ============================================================
# Stratified (size-ratio) analyses with Mantel-Haenszel
# ============================================================
# GABA guest: GLU host vs NonNeuron host
gaba_gn <- combined %>%
  filter(guest_type == "Gaba", host_type %in% c("Glut", "NonNeuron")) %>%
  filter(host_size >= guest_size) %>%
  mutate(
    is_Glut_host = host_type == "Glut",
    ratio_bin = ntile(size_ratio_host_over_guest, 5)
  )

table24 <- gaba_gn %>%
  group_by(ratio_bin, is_Glut_host) %>%
  summarise(
    full_n = sum(full_containment),
    total_n = n(),
    .groups = "drop"
  ) %>%
  mutate(rate = full_n / total_n) %>%
  select(ratio_bin, is_Glut_host, rate, total_n) %>%
  pivot_wider(names_from = is_Glut_host, values_from = c(rate, total_n),
              names_glue = "{ifelse(is_Glut_host, 'Glut_host', 'NonNeuron_host')}_{.value}") %>%
  rename(
    Glut_host_rate = Glut_host_rate,
    NonNeuron_host_rate = NonNeuron_host_rate,
    n_pairs_Glut = Glut_host_total_n,
    n_pairs_NonNeuron = NonNeuron_host_total_n
  ) %>%
  mutate(n_pairs = n_pairs_Glut + n_pairs_NonNeuron)

write_csv(table24, file.path(tabdir, "table24_size_ratio_stratified_gaba_glut_vs_nonneuron.csv"))

# Mantel-Haenszel
make_2x2xk <- function(df, row_binary, full_binary, strata) {
  lev <- sort(unique(df[[strata]]))
  arr <- array(0, dim = c(2, 2, length(lev)),
               dimnames = list(
                 row = c("TRUE", "FALSE"),
                 col = c("TRUE", "FALSE"),
                 strata = as.character(lev)
               ))
  for (s in seq_along(lev)) {
    d <- df %>% filter(.data[[strata]] == lev[s])
    tab <- table(factor(d[[row_binary]], levels = c(TRUE, FALSE)),
                 factor(d[[full_binary]], levels = c(TRUE, FALSE)))
    arr[, , s] <- tab
  }
  arr
}

mh_gn_array <- make_2x2xk(gaba_gn, "is_Glut_host", "full_containment", "ratio_bin")
mh_gn <- mantelhaen.test(mh_gn_array)

table25 <- tibble(
  comparison = "GABA guest, Glut host vs NonNeuron host, stratified by size-ratio quintiles",
  pooled_OR = unname(mh_gn$estimate),
  ci_low = mh_gn$conf.int[1],
  ci_high = mh_gn$conf.int[2],
  CMH_statistic = unname(mh_gn$statistic),
  CMH_p = mh_gn$p.value
)

write_csv(table25, file.path(tabdir, "table25_cmh_gaba_glut_vs_nonneuron_size_stratified.csv"))

# I⊂E vs E⊂I
ie_dir <- combined %>%
  filter(guest_type %in% c("Gaba", "Glut"),
         host_type %in% c("Gaba", "Glut"),
         guest_type != host_type) %>%
  filter(host_size >= guest_size) %>%
  mutate(
    guest_is_Gaba = guest_type == "Gaba",
    ratio_bin = ntile(size_ratio_host_over_guest, 5)
  )

table26 <- ie_dir %>%
  group_by(ratio_bin, guest_is_Gaba) %>%
  summarise(
    full_n = sum(full_containment),
    total_n = n(),
    .groups = "drop"
  ) %>%
  mutate(rate = full_n / total_n) %>%
  select(ratio_bin, guest_is_Gaba, rate, total_n) %>%
  pivot_wider(names_from = guest_is_Gaba, values_from = c(rate, total_n),
              names_glue = "{ifelse(guest_is_Gaba, 'I_in_E', 'E_in_I')}_{.value}") %>%
  rename(
    I_in_E_rate = I_in_E_rate,
    E_in_I_rate = E_in_I_rate,
    n_pairs_I = I_in_E_total_n,
    n_pairs_E = E_in_I_total_n
  ) %>%
  mutate(n_pairs = n_pairs_I + n_pairs_E)

write_csv(table26, file.path(tabdir, "table26_size_ratio_stratified_i_in_e_vs_e_in_i.csv"))

mh_ie_array <- make_2x2xk(ie_dir, "guest_is_Gaba", "full_containment", "ratio_bin")
mh_ie <- mantelhaen.test(mh_ie_array)

table27 <- tibble(
  comparison = "I⊂E vs E⊂I, stratified by size-ratio quintiles",
  pooled_OR = unname(mh_ie$estimate),
  ci_low = mh_ie$conf.int[1],
  ci_high = mh_ie$conf.int[2],
  CMH_statistic = unname(mh_ie$statistic),
  CMH_p = mh_ie$p.value
)

write_csv(table27, file.path(tabdir, "table27_cmh_i_in_e_vs_e_in_i_size_stratified.csv"))

# ============================================================
# Directional binomial tests
# ============================================================
binom_row <- function(label, a, b) {
  bt <- binom.test(a, a + b, p = 0.5)
  tibble(
    comparison = label,
    A_n = a,
    B_n = b,
    A_fraction = a / (a + b),
    binom_p_two_sided = bt$p.value
  )
}

table28 <- bind_rows(
  binom_row("I⊂E_vs_E⊂I", 1857, 287),
  binom_row("I⊂N_vs_N⊂I", 1109, 1352),
  binom_row("E⊂N_vs_N⊂E", 1150, 6062)
)

write_csv(table28, file.path(tabdir, "table28_directional_binomial_tests.csv"))

# ============================================================
# A compact key-findings table
# ============================================================
table00 <- tribble(
  ~Finding, ~I_in_E, ~Comparator, ~Statistic,
  "1. I⊂E full rate", "1857 / 8457 = 21.96%", "vs E⊂I 287 / 8457 = 3.39%", "Fisher OR=8.01; p=4.25e-317; binomial p=8.29e-281",
  "2. Size-opportunity adjusted", "28.46% when host_size>=guest_size", "vs I⊂N 15.50%; E⊂E 20.59%", "Fisher p(I⊂E vs I⊂N)=1.09e-75; p(I⊂E vs E⊂E)=3.18e-26",
  "3. Host significance", "host_log10p median 7.14", "vs I⊂N 2.23; E⊂I 1.90", "MWU p=8.37e-260 and 9.76e-126",
  "4. Host-guest significance gap", "median 5.43", "vs I⊂N 0.57; E⊂E 2.62", "MWU p=2.54e-257 and 2.03e-52",
  "5. Same GABA guest: GLU host size", "527", "vs matched NonNeuron host 336", "paired Wilcoxon p=3.39e-50 (n=489)",
  "6. Same GABA guest: GLU host fill", "0.199", "vs matched NonNeuron host 0.335", "paired Wilcoxon p=2.40e-39 (n=489)",
  "7. Same GLU host: GABA guest size", "135", "vs matched NonNeuron guest 94", "paired Wilcoxon p=1.38e-37 (n=958)",
  "8. Same GLU host: n guests", "median GABA guests=1", "vs matched NonNeuron guests=3", "paired Wilcoxon p=5.23e-110 (n=958)",
  "9. Partial asymmetry", "median +0.146", "vs I⊂N +0.0009", "MWU p=4.79e-298; slide-level p=3.02e-16",
  "10. GABA cluster-level preference", "1359 with any GLU host", "vs 820 with any NonNeuron host", "McNemar p=2.94e-56"
)
write_csv(table00, file.path(tabdir, "table00_key_findings_summary.csv"))


# ============================================================
# Figures
# ============================================================
fmt_p <- function(p) {
  if (is.na(p)) return('NA')
  if (p < 1e-300) return('<1e-300')
  if (p < 1e-4) return(formatC(p, format = 'e', digits = 2))
  sprintf('%.4f', p)
}

fmt_pct_lab <- function(x, digits = 1) sprintf(paste0('%.', digits, 'f%%'), 100 * x)
wrap_title <- function(x, width = 58) stringr::str_wrap(x, width = width)
wrap_subtitle <- function(x, width = 92) stringr::str_wrap(x, width = width)

pal_motif <- c(
  "I⊂E" = "#B22222",
  "E⊂I" = "#D98C8C",
  "I⊂N" = "#4E79A7",
  "N⊂I" = "#9BBCE0",
  "E⊂N" = "#59A14F",
  "N⊂E" = "#A8D39E",
  "I⊂I" = "#7F3C8D",
  "E⊂E" = "#F28E2B",
  "N⊂N" = "#8C8C8C"
)

col_compare <- c(
  "Glut host" = "#B22222",
  "NonNeuron host" = "#4E79A7",
  "GABA guest" = "#B22222",
  "NonNeuron guest" = "#4E79A7"
)

figure_theme <- ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    plot.title.position = 'plot',
    plot.caption.position = 'plot',
    plot.title = ggplot2::element_text(face = 'bold', size = 16.5, colour = 'black', margin = ggplot2::margin(b = 8)),
    plot.subtitle = ggplot2::element_text(size = 11.5, colour = '#4D4D4D', lineheight = 1.08, margin = ggplot2::margin(b = 10)),
    axis.title = ggplot2::element_text(size = 13.5, colour = 'black'),
    axis.text = ggplot2::element_text(size = 11.2, colour = 'black'),
    axis.line = ggplot2::element_line(linewidth = 0.55, colour = 'black'),
    axis.ticks = ggplot2::element_line(linewidth = 0.5, colour = 'black'),
    legend.position = 'top',
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 11),
    plot.margin = ggplot2::margin(18, 22, 14, 18),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank()
  )

theme_set(figure_theme)

save_pub <- function(plot_obj, filename, width, height, dpi = 600) {
  ggplot2::ggsave(filename = file.path(figdir, filename), plot = plot_obj,
                  width = width, height = height, dpi = dpi, bg = 'white', limitsize = FALSE)
  pdf_name <- sub('\\.png$', '.pdf', filename)
  ggplot2::ggsave(filename = file.path(figdir, pdf_name), plot = plot_obj,
                  width = width, height = height, bg = 'white', device = cairo_pdf, limitsize = FALSE)
}

# Figure 1 ---------------------------------------------------
fig01_stats_a <- table03 %>% dplyr::filter(compare_to == 'E⊂I')
fig01_stats_b <- table03 %>% dplyr::filter(compare_to == 'I⊂N')
fig01_df <- table01 %>%
  dplyr::mutate(
    motif = factor(motif, levels = motif_order),
    lab = paste0(fmt_pct_lab(full_rate), '\n', 'n=', scales::comma(full_n))
  )

fig01 <- ggplot(fig01_df, aes(x = motif, y = full_rate, fill = motif)) +
  geom_col(width = 0.72, colour = 'white', linewidth = 0.35) +
  geom_text(aes(label = lab), vjust = -0.28, size = 3.65, lineheight = 0.95) +
  scale_fill_manual(values = pal_motif, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.18))) +
  labs(
    x = NULL,
    y = 'Full containment rate',
    title = wrap_title('I⊂E is the dominant full-containment motif across directed interaction classes'),
    subtitle = wrap_subtitle(
      paste0(
        'I⊂E = ', fmt_pct_lab(fig01_df$full_rate[fig01_df$motif == 'I⊂E'], 2),
        '; vs E⊂I OR = ', sprintf('%.2f', fig01_stats_a$OR_raw), ', FDR = ', fmt_p(fig01_stats_a$fdr_raw),
        '; vs I⊂N OR = ', sprintf('%.2f', fig01_stats_b$OR_raw), ', FDR = ', fmt_p(fig01_stats_b$fdr_raw)
      )
    )
  ) +
  guides(fill = 'none')
save_pub(fig01, 'figure01_full_containment_rate_by_motif.png', 9.6, 5.6)

# Figure 2 ---------------------------------------------------
fig02_df <- table04 %>%
  dplyr::filter(motif %in% c('N⊂E', 'I⊂E', 'E⊂E', 'I⊂N', 'I⊂I')) %>%
  dplyr::mutate(
    motif = factor(motif, levels = c('N⊂E', 'I⊂E', 'E⊂E', 'I⊂N', 'I⊂I')),
    lab = sprintf('%.3f', host_fill_median)
  )
fig02_p1 <- table05 %>% dplyr::filter(compare_to == 'I⊂N', metric == 'host_fill')
fig02_p2 <- table05 %>% dplyr::filter(compare_to == 'E⊂E', metric == 'host_fill')

fig02 <- ggplot(fig02_df, aes(x = motif, y = host_fill_median, fill = motif)) +
  geom_col(width = 0.72, colour = 'white', linewidth = 0.35) +
  geom_text(aes(label = lab), vjust = -0.25, size = 3.7) +
  scale_fill_manual(values = pal_motif) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(
    x = NULL,
    y = 'Median host fill (= Jaccard under full containment)',
    title = wrap_title('I⊂E sits in a distinct containment-tightness regime'),
    subtitle = wrap_subtitle(
      paste0('Higher fill than I⊂N (FDR = ', fmt_p(fig02_p1$fdr), ') and E⊂E (FDR = ', fmt_p(fig02_p2$fdr), ').')
    )
  ) +
  guides(fill = 'none')
save_pub(fig02, 'figure02_host_fill_spectrum_selected_motifs.png', 8.2, 5.4)

# Figure 3 ---------------------------------------------------
fig03_df <- table04 %>%
  dplyr::filter(motif %in% c('I⊂E', 'N⊂E', 'E⊂E', 'I⊂N', 'N⊂N')) %>%
  dplyr::mutate(
    motif = factor(motif, levels = c('I⊂E', 'N⊂E', 'E⊂E', 'I⊂N', 'N⊂N')),
    lab = sprintf('%.2f', log10p_gap_median)
  )
fig03_p1 <- table05 %>% dplyr::filter(compare_to == 'I⊂N', metric == 'logp_gap')
fig03_p2 <- table05 %>% dplyr::filter(compare_to == 'E⊂E', metric == 'logp_gap')

fig03 <- ggplot(fig03_df, aes(x = motif, y = log10p_gap_median, fill = motif)) +
  geom_col(width = 0.72, colour = 'white', linewidth = 0.35) +
  geom_text(aes(label = lab), vjust = -0.25, size = 3.7) +
  scale_fill_manual(values = pal_motif) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(
    x = NULL,
    y = 'Median host_log10p - guest_log10p',
    title = wrap_title('I⊂E has the largest host–guest significance gap'),
    subtitle = wrap_subtitle(
      paste0('Compared with I⊂N: FDR = ', fmt_p(fig03_p1$fdr), '; compared with E⊂E: FDR = ', fmt_p(fig03_p2$fdr), '.')
    )
  ) +
  guides(fill = 'none')
save_pub(fig03, 'figure03_host_guest_significance_gap_selected_motifs.png', 8.6, 5.4)

# Figure 4 ---------------------------------------------------
fig04_df <- guest_host_agg %>% dplyr::select(host_size_med_Glut, host_size_med_NonNeuron) %>% tidyr::drop_na()
lim04 <- max(fig04_df$host_size_med_Glut, fig04_df$host_size_med_NonNeuron, na.rm = TRUE)
fig04_stat <- table09 %>% dplyr::filter(comparison == 'same_GABA_guest_Glut_host_vs_NonNeuron_host', metric == 'host_size')
fig04_cor <- suppressWarnings(cor.test(fig04_df$host_size_med_NonNeuron, fig04_df$host_size_med_Glut, method = 'spearman'))
label04 <- paste0(
  'n = ', nrow(fig04_df),
  '\npaired p = ', fmt_p(fig04_stat$p),
  '\nmedian(GLU) = ', sprintf('%.0f', fig04_stat$A_median),
  '\nmedian(NN) = ', sprintf('%.0f', fig04_stat$B_median),
  '\nSpearman ρ = ', sprintf('%.2f', unname(fig04_cor$estimate)), ', p = ', fmt_p(fig04_cor$p.value)
)

fig04 <- ggplot(fig04_df, aes(x = host_size_med_NonNeuron, y = host_size_med_Glut)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, linewidth = 0.6, colour = '#7F7F7F') +
  geom_point(size = 1.8, alpha = 0.52, colour = '#B22222') +
  annotate('label', x = lim04 * 0.10, y = lim04 * 0.92, hjust = 0, vjust = 1,
           label = label04, size = 3.5, label.size = 0.25, fill = 'white') +
  coord_equal(xlim = c(0, lim04 * 1.02), ylim = c(0, lim04 * 1.02), clip = 'off') +
  labs(
    x = 'Matched NonNeuron host median size',
    y = 'Matched Glut host median size',
    title = wrap_title('For the same GABA guest, matched Glut hosts are systematically larger')
  )
save_pub(fig04, 'figure04_matched_same_gaba_guest_host_size_glut_vs_nonneuron.png', 6.3, 6.1)

# Figure 5 ---------------------------------------------------
fig05_df <- guest_host_agg %>% dplyr::select(host_fill_med_Glut, host_fill_med_NonNeuron) %>% tidyr::drop_na()
lim05 <- max(fig05_df$host_fill_med_Glut, fig05_df$host_fill_med_NonNeuron, na.rm = TRUE)
fig05_stat <- table09 %>% dplyr::filter(comparison == 'same_GABA_guest_Glut_host_vs_NonNeuron_host', metric == 'host_fill')
fig05_cor <- suppressWarnings(cor.test(fig05_df$host_fill_med_NonNeuron, fig05_df$host_fill_med_Glut, method = 'spearman'))
label05 <- paste0(
  'n = ', nrow(fig05_df),
  '\npaired p = ', fmt_p(fig05_stat$p),
  '\nmedian(GLU) = ', sprintf('%.3f', fig05_stat$A_median),
  '\nmedian(NN) = ', sprintf('%.3f', fig05_stat$B_median),
  '\nSpearman ρ = ', sprintf('%.2f', unname(fig05_cor$estimate)), ', p = ', fmt_p(fig05_cor$p.value)
)

fig05 <- ggplot(fig05_df, aes(x = host_fill_med_NonNeuron, y = host_fill_med_Glut)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, linewidth = 0.6, colour = '#7F7F7F') +
  geom_point(size = 1.8, alpha = 0.52, colour = '#B22222') +
  annotate('label', x = lim05 * 0.08, y = lim05 * 0.94, hjust = 0, vjust = 1,
           label = label05, size = 3.5, label.size = 0.25, fill = 'white') +
  coord_equal(xlim = c(0, lim05 * 1.02), ylim = c(0, lim05 * 1.02), clip = 'off') +
  labs(
    x = 'Matched NonNeuron host median fill',
    y = 'Matched Glut host median fill',
    title = wrap_title('For the same GABA guest, Glut hosts form looser shells than NonNeuron hosts')
  )
save_pub(fig05, 'figure05_matched_same_gaba_guest_host_fill_glut_vs_nonneuron.png', 6.3, 6.1)

# Figure 6 ---------------------------------------------------
fig06_df <- host_guest_agg %>% dplyr::select(guest_size_med_Gaba, guest_size_med_NonNeuron) %>% tidyr::drop_na()
lim06 <- max(fig06_df$guest_size_med_Gaba, fig06_df$guest_size_med_NonNeuron, na.rm = TRUE)
fig06_stat <- table10 %>% dplyr::filter(comparison == 'same_Glut_host_Gaba_guest_vs_NonNeuron_guest', metric == 'guest_size')
fig06_cor <- suppressWarnings(cor.test(fig06_df$guest_size_med_NonNeuron, fig06_df$guest_size_med_Gaba, method = 'spearman'))
label06 <- paste0(
  'n = ', nrow(fig06_df),
  '\npaired p = ', fmt_p(fig06_stat$p),
  '\nmedian(GABA) = ', sprintf('%.0f', fig06_stat$A_median),
  '\nmedian(NN) = ', sprintf('%.0f', fig06_stat$B_median),
  '\nSpearman ρ = ', sprintf('%.2f', unname(fig06_cor$estimate)), ', p = ', fmt_p(fig06_cor$p.value)
)

fig06 <- ggplot(fig06_df, aes(x = guest_size_med_NonNeuron, y = guest_size_med_Gaba)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, linewidth = 0.6, colour = '#7F7F7F') +
  geom_point(size = 1.8, alpha = 0.52, colour = '#B22222') +
  annotate('label', x = lim06 * 0.10, y = lim06 * 0.92, hjust = 0, vjust = 1,
           label = label06, size = 3.5, label.size = 0.25, fill = 'white') +
  coord_equal(xlim = c(0, lim06 * 1.02), ylim = c(0, lim06 * 1.02), clip = 'off') +
  labs(
    x = 'Matched NonNeuron guest median size',
    y = 'Matched GABA guest median size',
    title = wrap_title('For the same Glut host, GABA guests are larger than matched NonNeuron guests')
  )
save_pub(fig06, 'figure06_matched_same_glu_host_guest_size_gaba_vs_nonneuron.png', 6.3, 6.1)

# Figure 7 ---------------------------------------------------
fig07_slide <- partial %>%
  dplyr::filter(motif == 'I⊂E') %>%
  dplyr::group_by(slide) %>%
  dplyr::summarise(slide_med = median(overlap_asym, na.rm = TRUE), .groups = 'drop')

fig07_slide_stat <- tibble::tibble(
  median_signed_asym = median(fig07_slide$slide_med, na.rm = TRUE),
  slide_level_p_vs_zero = {
    x <- fig07_slide$slide_med
    x <- x[!is.na(x)]
    x_nz <- x[x != 0]
    if (length(x_nz) == 0) 1 else suppressWarnings(wilcox.test(x_nz, mu = 0, exact = FALSE)$p.value)
  }
)

fig07_ref <- table06 %>%
  dplyr::filter(motif == 'I⊂E') %>%
  dplyr::transmute(motif, asym_median, fdr = NA_real_)

fig07_comp <- table06 %>%
  dplyr::filter(motif %in% c('I⊂N', 'I⊂I', 'N⊂E', 'E⊂I', 'E⊂N')) %>%
  dplyr::left_join(
    table07 %>% dplyr::rename(motif = compare_to),
    by = 'motif'
  ) %>%
  dplyr::select(motif, asym_median, fdr)

fig07_df <- dplyr::bind_rows(fig07_ref, fig07_comp) %>%
  dplyr::mutate(
    motif = factor(motif, levels = c('I⊂E', 'I⊂N', 'I⊂I', 'N⊂E', 'E⊂I', 'E⊂N')),
    stat_lab = dplyr::case_when(
      as.character(motif) == 'I⊂E' ~ paste0('median=', sprintf('%.3f', asym_median), '\nref'),
      TRUE ~ paste0('median=', sprintf('%.3f', asym_median), '\nFDR=', purrr::map_chr(fdr, fmt_p))
    ),
    vjust_lab = ifelse(asym_median >= 0, -0.24, 1.14)
  )

fig07 <- ggplot(fig07_df, aes(x = motif, y = asym_median, fill = motif)) +
  geom_col(width = 0.72, colour = 'white', linewidth = 0.35) +
  geom_hline(yintercept = 0, linewidth = 0.6, colour = '#7F7F7F') +
  geom_text(aes(label = stat_lab, vjust = vjust_lab), size = 3.35, lineheight = 0.95, show.legend = FALSE) +
  scale_fill_manual(values = pal_motif) +
  labs(
    x = NULL,
    y = 'Median signed overlap asymmetry\n(guest overlap − host overlap)',
    title = wrap_title('Partial-overlap directionality is strongest and most stable for I⊂E'),
    subtitle = wrap_subtitle(
      paste0('I⊂E slide-level p vs 0 = ', fmt_p(fig07_slide_stat$slide_level_p_vs_zero),
             '; median slide asymmetry = ', sprintf('%.3f', fig07_slide_stat$median_signed_asym), '.')
    )
  ) +
  guides(fill = 'none')
save_pub(fig07, 'figure07_partial_overlap_asymmetry_selected_motifs.png', 8.8, 5.6)

# Figure 8 ---------------------------------------------------
fig08_df <- table24 %>%
  tidyr::pivot_longer(cols = c(Glut_host_rate, NonNeuron_host_rate), names_to = 'host_kind', values_to = 'rate') %>%
  dplyr::mutate(
    host_kind = dplyr::recode(host_kind,
                              Glut_host_rate = 'Glut host',
                              NonNeuron_host_rate = 'NonNeuron host'),
    label = fmt_pct_lab(rate)
  )
fig08_stat <- table25

fig08 <- ggplot(fig08_df, aes(x = factor(ratio_bin), y = rate, colour = host_kind, group = host_kind)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.5) +
  geom_text(aes(label = label), size = 3.2, vjust = -0.75, show.legend = FALSE) +
  scale_colour_manual(values = col_compare[c('Glut host','NonNeuron host')]) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0.02, 0.15))) +
  labs(
    x = 'Host/guest size-ratio quintile',
    y = 'Full containment rate',
    colour = NULL,
    title = wrap_title('Across size-ratio strata, Glut hosts remain more likely than NonNeuron hosts to fully contain GABA guests'),
    subtitle = wrap_subtitle(
      paste0('CMH pooled OR = ', sprintf('%.2f', fig08_stat$pooled_OR),
             ' (95% CI ', sprintf('%.2f', fig08_stat$ci_low), '–', sprintf('%.2f', fig08_stat$ci_high),
             '), p = ', fmt_p(fig08_stat$CMH_p), '.')
    )
  )
save_pub(fig08, 'figure08_size_ratio_stratified_gaba_glut_vs_nonneuron.png', 8.4, 5.4)

# Figure 9 ---------------------------------------------------
fig09_df <- table16 %>%
  dplyr::filter(n_clusters >= 50) %>%
  dplyr::mutate(guest_subclass_short = stringr::str_replace(guest_subclass, ' Gaba$', '')) %>%
  dplyr::slice_max(order_by = pct_any_Glut, n = 10) %>%
  dplyr::arrange(delta) %>%
  dplyr::mutate(guest_subclass_short = factor(guest_subclass_short, levels = guest_subclass_short))

xmax09 <- max(fig09_df$pct_any_Glut, fig09_df$pct_any_NonNeuron, na.rm = TRUE)
fig09 <- ggplot(fig09_df, aes(y = guest_subclass_short)) +
  geom_segment(aes(x = pct_any_NonNeuron, xend = pct_any_Glut, yend = guest_subclass_short),
               linewidth = 0.8, colour = '#BDBDBD') +
  geom_point(aes(x = pct_any_NonNeuron, colour = 'NonNeuron host'), size = 2.7) +
  geom_point(aes(x = pct_any_Glut, colour = 'Glut host'), size = 2.9) +
  geom_text(aes(
    x = pmax(pct_any_Glut, pct_any_NonNeuron) + 0.035,
    label = paste0('FDR=', purrr::map_chr(fdr, fmt_p))
  ),
  hjust = 0, size = 3.2, colour = '#444444') +
  scale_colour_manual(values = col_compare[c('Glut host','NonNeuron host')]) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, min(1, xmax09 + 0.16)),
                     expand = expansion(mult = c(0.01, 0.01))) +
  coord_cartesian(clip = 'off') +
  labs(
    x = 'Fraction of GABA clusters with any full host',
    y = NULL,
    colour = NULL,
    title = wrap_title('Main GABA subclasses show a systematic preference for full Glut hosts over full NonNeuron hosts'),
    subtitle = wrap_subtitle('Each row shows matched subclass-level probabilities; right-side labels report subclass-level FDR.')
  )
save_pub(fig09, 'figure09_main_gaba_subclass_host_preference.png', 9.4, 5.8)

# Figure 10 --------------------------------------------------
fig10_df <- table19 %>%
  dplyr::filter(n_clusters >= 50) %>%
  dplyr::mutate(host_subclass_short = stringr::str_replace(host_subclass, ' Glut$', '')) %>%
  dplyr::slice_max(order_by = pct_any_Gaba_guest, n = 8) %>%
  dplyr::arrange(delta) %>%
  dplyr::mutate(host_subclass_short = factor(host_subclass_short, levels = host_subclass_short))

xmax10 <- max(fig10_df$pct_any_Gaba_guest, fig10_df$pct_any_NonNeuron_guest, na.rm = TRUE)
fig10 <- ggplot(fig10_df, aes(y = host_subclass_short)) +
  geom_segment(aes(x = pct_any_NonNeuron_guest, xend = pct_any_Gaba_guest, yend = host_subclass_short),
               linewidth = 0.8, colour = '#BDBDBD') +
  geom_point(aes(x = pct_any_NonNeuron_guest, colour = 'NonNeuron guest'), size = 2.7) +
  geom_point(aes(x = pct_any_Gaba_guest, colour = 'GABA guest'), size = 2.9) +
  geom_text(aes(x = pmax(pct_any_Gaba_guest, pct_any_NonNeuron_guest) + 0.035,
                label = paste0('FDR=', purrr::map_chr(fdr, fmt_p))),
            hjust = 0, size = 3.2, colour = '#444444') +
  scale_colour_manual(values = col_compare[c('GABA guest','NonNeuron guest')]) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, min(1, xmax10 + 0.16)),
                     expand = expansion(mult = c(0.01, 0.01))) +
  coord_cartesian(clip = 'off') +
  labs(
    x = 'Fraction of Glut host clusters with any full guest',
    y = NULL,
    colour = NULL,
    title = wrap_title('Selected Glut subclasses are more likely to host full GABA guests than full NonNeuron guests'),
    subtitle = wrap_subtitle('Each row shows subclass-level host propensity; right-side labels report subclass-level FDR.')
  )
save_pub(fig10, 'figure10_selected_glu_subclass_guest_preference.png', 9.4, 5.6)

message('Analysis complete. Results written to: ', normalizePath(outdir))
