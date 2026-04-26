# 清空环境
rm(list = ls())
gc()
options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(scales)
  library(forcats)
  library(tibble)
  library(pROC)
})

# ============================================================
# R_output_11
# Causal-priority model for containment motifs across:
#   1) neuron-neuron
#   2) neuron-nonneuron
#
# Main biological questions addressed:
#   - Why Gaba⊂Glut is the highest-priority causal motif
#   - Why smaller⊂larger Glut is the second-tier candidate
#   - Why Gaba⊂NonNeuron is more consistent with a geometric/byproduct regime
#
# Output:
#   base_dir/R_output_11/
#     ├─ tables/
#     └─ figures/
# ============================================================

# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------
base_dir <- "E:/zaw/2603"

nn_candidates <- c(
  file.path(base_dir, "neuron-neuron-partner.csv"),
  file.path(base_dir, "neuron-neuron-partner.txt")
)
nonn_candidates <- c(
  file.path(base_dir, "neuron-nonneuron-partner.csv"),
  file.path(base_dir, "neuron-nonneuron-partner.txt")
)

nn_file <- nn_candidates[file.exists(nn_candidates)][1]
nonn_file <- nonn_candidates[file.exists(nonn_candidates)][1]

if (is.na(nn_file) || !nzchar(nn_file)) {
  stop("Cannot find neuron-neuron partner file under: ", base_dir)
}
if (is.na(nonn_file) || !nzchar(nonn_file)) {
  stop("Cannot find neuron-nonneuron partner file under: ", base_dir)
}

out_dir <- file.path(base_dir, "R_output_11")
fig_dir <- file.path(out_dir, "figures")
table_dir <- file.path(out_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------
# Formatting helpers
# ------------------------------------------------------------------
fmt_p <- function(p) {
  if (length(p) != 1) stop("fmt_p expects a scalar p-value")
  if (is.na(p)) return("NA")
  if (p < 1e-300) return("<1e-300")
  if (p < 1e-4) return(format(p, scientific = TRUE, digits = 2))
  sprintf("%.4f", p)
}
fmt_pct <- function(x, digits = 1) sprintf(paste0("%.", digits, "f%%"), 100 * x)
fmt_num <- function(x, digits = 3) sprintf(paste0("%.", digits, "f"), x)
wrap_title <- function(x, width = 68) stringr::str_wrap(x, width = width)
wrap_subtitle <- function(x, width = 102) stringr::str_wrap(x, width = width)

col_motif <- c(
  "Gaba⊂Glut" = "#C44E52",
  "smaller⊂Glut (GG)" = "#4E79A7",
  "Gaba⊂NonNeuron" = "#59A14F",
  "NonNeuron⊂Glut" = "#9C755F",
  "smaller⊂Gaba (II)" = "#B07AA1"
)
col_binary <- c("Higher" = "#C44E52", "Lower" = "#4E79A7")
col_partner <- c("Glut host" = "#C44E52", "NonNeuron host" = "#59A14F")

base_theme <- function() {
  theme_classic(base_size = 12) +
    theme(
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.margin = margin(18, 26, 16, 16),
      plot.title = element_text(face = "bold", size = 16, lineheight = 1.02, margin = margin(b = 8)),
      plot.subtitle = element_text(size = 11.4, colour = "#4D4D4D", lineheight = 1.10, margin = margin(b = 10)),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11, colour = "black"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 10.8),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 12)
    )
}

save_pub <- function(p, filename, width = 8.2, height = 5.4, dpi = 600) {
  ggsave(file.path(fig_dir, filename), p, width = width, height = height, dpi = dpi,
         bg = "white", limitsize = FALSE)
  pdf_name <- sub("\\.png$", ".pdf", filename)
  ggsave(file.path(fig_dir, pdf_name), p, width = width, height = height,
         bg = "white", device = cairo_pdf, limitsize = FALSE)
}

# ------------------------------------------------------------------
# Data helpers
# ------------------------------------------------------------------
region_breadth <- function(x) {
  if (is.na(x) || x == "") return(NA_integer_)
  length(str_split(as.character(x), ",", simplify = FALSE)[[1]])
}

safe_auc <- function(y, score) {
  y <- as.integer(y)
  ok <- is.finite(score) & !is.na(y)
  y <- y[ok]
  score <- score[ok]
  if (length(unique(y)) < 2) return(NA_real_)
  out <- tryCatch(as.numeric(pROC::auc(response = y, predictor = score)), error = function(e) NA_real_)
  out
}

paired_wilcox_p <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  suppressWarnings(wilcox.test(x, y, paired = TRUE, alternative = "two.sided", exact = FALSE)$p.value)
}

one_sample_wilcox_p <- function(x, mu = 0) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  suppressWarnings(wilcox.test(x, mu = mu, exact = FALSE)$p.value)
}

cramers_v_binary_by_group <- function(group_code, y) {
  tab <- table(group_code, y)
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NA_real_)
  chi <- suppressWarnings(chisq.test(tab, correct = FALSE)$statistic)
  n <- sum(tab)
  as.numeric(sqrt(chi / n))
}

conditional_perm_specificity <- function(df, outcome_col, pair_cols, strata_cols, n_perm = 1000, seed = 1) {
  set.seed(seed)
  dat <- df %>%
    dplyr::select(dplyr::all_of(c(outcome_col, pair_cols, strata_cols))) %>%
    tidyr::drop_na()
  if (nrow(dat) == 0) {
    return(list(obs_v = NA_real_, perm_mean_v = NA_real_, perm_sd_v = NA_real_, z = NA_real_, p = NA_real_, n_pairs = 0, n_strata = 0))
  }
  
  pair_id <- interaction(dat[, pair_cols], drop = TRUE, lex.order = TRUE)
  strata <- interaction(dat[, strata_cols], drop = TRUE, lex.order = TRUE)
  y <- as.integer(dat[[outcome_col]])
  
  obs_v <- cramers_v_binary_by_group(pair_id, y)
  strata_idx <- split(seq_len(nrow(dat)), strata)
  perm_v <- numeric(n_perm)
  
  for (b in seq_len(n_perm)) {
    y_perm <- y
    for (idx in strata_idx) {
      if (length(idx) > 1) y_perm[idx] <- sample(y_perm[idx], replace = FALSE)
    }
    perm_v[b] <- cramers_v_binary_by_group(pair_id, y_perm)
  }
  
  list(
    obs_v = obs_v,
    perm_mean_v = mean(perm_v, na.rm = TRUE),
    perm_sd_v = sd(perm_v, na.rm = TRUE),
    z = (obs_v - mean(perm_v, na.rm = TRUE)) / sd(perm_v, na.rm = TRUE),
    p = (sum(perm_v >= obs_v, na.rm = TRUE) + 1) / (sum(is.finite(perm_v)) + 1),
    n_pairs = nlevels(pair_id),
    n_strata = nlevels(strata)
  )
}

partial_asym_stats <- function(df, fill_a, fill_b, label_text) {
  partial <- df %>% dplyr::filter(.data[[fill_a]] < 1, .data[[fill_b]] < 1)
  partial <- partial %>% dplyr::mutate(asym = .data[[fill_a]] - .data[[fill_b]])
  slide_median <- partial %>%
    dplyr::group_by(slide) %>%
    dplyr::summarise(med = median(asym, na.rm = TRUE), .groups = "drop")
  
  tibble(
    name = label_text,
    n_partial = nrow(partial),
    pair_median = median(partial$asym, na.rm = TRUE),
    pair_mean = mean(partial$asym, na.rm = TRUE),
    pair_wilcoxon_p = one_sample_wilcox_p(partial$asym, 0),
    slide_n = nrow(slide_median),
    slide_positive_n = sum(slide_median$med > 0, na.rm = TRUE),
    slide_median_of_medians = median(slide_median$med, na.rm = TRUE),
    slide_wilcoxon_p = one_sample_wilcox_p(slide_median$med, 0)
  )
}

mcnemar_exact <- function(x, y) {
  b <- sum(x == 1 & y == 0, na.rm = TRUE)
  c <- sum(x == 0 & y == 1, na.rm = TRUE)
  if ((b + c) == 0) return(list(b_10 = b, c_01 = c, p = NA_real_))
  p <- binom.test(min(b, c), b + c, p = 0.5, alternative = "two.sided")$p.value
  list(b_10 = b, c_01 = c, p = p)
}

dedupe_unordered <- function(df) {
  key <- purrr::pmap_chr(
    df[, c("cluster.1_label", "cluster.2_label")],
    function(cluster.1_label, cluster.2_label) paste(sort(c(cluster.1_label, cluster.2_label)), collapse = " || ")
  )
  df[!duplicated(key), , drop = FALSE]
}

prep_df <- function(df, source_name) {
  num_cols <- c(
    "cluster.1_total_cell_num", "cluster.2_total_cell_num",
    "cluster.1_E_I_Ratio", "cluster.2_E_I_Ratio",
    "overlap_cell", "union_cell", "jaccard",
    "cluster.1.overlap.percent", "cluster.2.overlap.percent",
    "cluster.1_cauchy_combination_p", "cluster.2_cauchy_combination_p"
  )
  for (cc in num_cols) df[[cc]] <- suppressWarnings(as.numeric(df[[cc]]))
  df$cluster.1_cauchy_combination_p[df$cluster.1_cauchy_combination_p == 0] <- .Machine$double.xmin
  df$cluster.2_cauchy_combination_p[df$cluster.2_cauchy_combination_p == 0] <- .Machine$double.xmin
  df$logp1 <- -log10(df$cluster.1_cauchy_combination_p)
  df$logp2 <- -log10(df$cluster.2_cauchy_combination_p)
  df$region_breadth1 <- vapply(df$cluster.1_region, region_breadth, integer(1))
  df$region_breadth2 <- vapply(df$cluster.2_region, region_breadth, integer(1))
  df$source <- source_name
  df
}

side_df <- function(df, idx) {
  tibble::tibble(
    label = df[[paste0("cluster.", idx, "_label")]],
    slide = df[[paste0("cluster.", idx, "_slide")]],
    layer = df[[paste0("cluster.", idx, "_layer")]],
    region = df[[paste0("cluster.", idx, "_region")]],
    type = df[[paste0("cluster.", idx, "_cell_Neruon_type")]],
    subclass = df[[paste0("cluster.", idx, "_subclass")]],
    n = df[[paste0("cluster.", idx, "_total_cell_num")]],
    logp = df[[paste0("logp", idx)]],
    fill = df[[paste0("cluster.", idx, ".overlap.percent")]],
    breadth = df[[paste0("region_breadth", idx)]]
  )
}

orient_cross_pair <- function(df, type_a, type_b, prefix_a, prefix_b) {
  keep <- (df$cluster.1_cell_Neruon_type == type_a & df$cluster.2_cell_Neruon_type == type_b) |
    (df$cluster.1_cell_Neruon_type == type_b & df$cluster.2_cell_Neruon_type == type_a)
  x <- df[keep, , drop = FALSE]
  out_list <- vector("list", nrow(x))
  
  for (i in seq_len(nrow(x))) {
    row <- x[i, , drop = FALSE]
    t1 <- row$cluster.1_cell_Neruon_type
    if (t1 == type_a) {
      a <- side_df(row, 1)
      b <- side_df(row, 2)
    } else {
      a <- side_df(row, 2)
      b <- side_df(row, 1)
    }
    out_list[[i]] <- tibble::tibble(
      slide = a$slide,
      overlap_cell = row$overlap_cell,
      jaccard = row$jaccard,
      !!paste0(prefix_a, "_label") := a$label,
      !!paste0(prefix_a, "_subclass") := a$subclass,
      !!paste0(prefix_a, "_layer") := a$layer,
      !!paste0(prefix_a, "_n") := a$n,
      !!paste0(prefix_a, "_logp") := a$logp,
      !!paste0(prefix_a, "_fill") := a$fill,
      !!paste0(prefix_a, "_region_breadth") := a$breadth,
      !!paste0(prefix_b, "_label") := b$label,
      !!paste0(prefix_b, "_subclass") := b$subclass,
      !!paste0(prefix_b, "_layer") := b$layer,
      !!paste0(prefix_b, "_n") := b$n,
      !!paste0(prefix_b, "_logp") := b$logp,
      !!paste0(prefix_b, "_fill") := b$fill,
      !!paste0(prefix_b, "_region_breadth") := b$breadth
    )
  }
  dplyr::bind_rows(out_list)
}

orient_same_type_by_size <- function(df, type_name) {
  keep <- df$cluster.1_cell_Neruon_type == type_name & df$cluster.2_cell_Neruon_type == type_name
  x <- df[keep, , drop = FALSE]
  out_list <- vector("list", nrow(x))
  
  for (i in seq_len(nrow(x))) {
    row <- x[i, , drop = FALSE]
    s1 <- side_df(row, 1)
    s2 <- side_df(row, 2)
    if (s1$n <= s2$n) {
      small <- s1; large <- s2
    } else {
      small <- s2; large <- s1
    }
    out_list[[i]] <- tibble::tibble(
      slide = small$slide,
      overlap_cell = row$overlap_cell,
      jaccard = row$jaccard,
      small_label = small$label,
      small_subclass = small$subclass,
      small_layer = small$layer,
      small_n = small$n,
      small_logp = small$logp,
      small_fill = small$fill,
      small_region_breadth = small$breadth,
      large_label = large$label,
      large_subclass = large$subclass,
      large_layer = large$layer,
      large_n = large$n,
      large_logp = large$logp,
      large_fill = large$fill,
      large_region_breadth = large$breadth
    )
  }
  dplyr::bind_rows(out_list)
}

# ------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------
nn <- readr::read_csv(nn_file, show_col_types = FALSE)
nonn <- readr::read_csv(nonn_file, show_col_types = FALSE)

nn_u <- prep_df(dedupe_unordered(nn), "neuron-neuron")
nonn_u <- prep_df(dedupe_unordered(nonn), "neuron-nonneuron")

nn_u$pair_type <- apply(
  nn_u[, c("cluster.1_cell_Neruon_type", "cluster.2_cell_Neruon_type")],
  1, function(z) paste(sort(z), collapse = "-")
)
nonn_u$pair_type <- apply(
  nonn_u[, c("cluster.1_cell_Neruon_type", "cluster.2_cell_Neruon_type")],
  1, function(z) paste(sort(z), collapse = "-")
)

GI <- nn_u[nn_u$pair_type == "Gaba-Glut", , drop = FALSE]
GG <- nn_u[nn_u$pair_type == "Glut-Glut", , drop = FALSE]
II <- nn_u[nn_u$pair_type == "Gaba-Gaba", , drop = FALSE]
GN <- nonn_u[nonn_u$pair_type == "Glut-NonNeuron", , drop = FALSE]
IN <- nonn_u[nonn_u$pair_type == "Gaba-NonNeuron", , drop = FALSE]
NN <- nonn_u[nonn_u$pair_type == "NonNeuron-NonNeuron", , drop = FALSE]

GIo <- orient_cross_pair(nn_u, "Gaba", "Glut", "G", "E")
GNo <- orient_cross_pair(nonn_u, "Glut", "NonNeuron", "E", "N")
INo <- orient_cross_pair(nonn_u, "Gaba", "NonNeuron", "G", "N")
GGs <- orient_same_type_by_size(GG, "Glut")
IIs <- orient_same_type_by_size(II, "Gaba")

GIo <- GIo %>% mutate(y_G_in_E = as.integer(G_fill == 1), log_size_ratio = log2((E_n + 1) / (G_n + 1)), size_bin = dplyr::ntile(log_size_ratio, 10))
GNo <- GNo %>% mutate(y_N_in_E = as.integer(N_fill == 1), log_size_ratio = log2((E_n + 1) / (N_n + 1)), size_bin = dplyr::ntile(log_size_ratio, 10))
INo <- INo %>% mutate(y_G_in_N = as.integer(G_fill == 1), log_size_ratio = log2((N_n + 1) / (G_n + 1)), size_bin = dplyr::ntile(log_size_ratio, 10))
GGs <- GGs %>% mutate(y_small_in_large = as.integer(small_fill == 1), log_size_ratio = log2((large_n + 1) / (small_n + 1)), size_bin = dplyr::ntile(log_size_ratio, 10))
IIs <- IIs %>% mutate(y_small_in_large = as.integer(small_fill == 1), log_size_ratio = log2((large_n + 1) / (small_n + 1)), size_bin = dplyr::ntile(log_size_ratio, 10))

# ------------------------------------------------------------------
# Residual identity and precursor analyses
# ------------------------------------------------------------------
res_GI <- conditional_perm_specificity(GIo, "y_G_in_E", c("G_subclass", "E_subclass"), c("slide", "size_bin", "G_layer", "E_layer"), n_perm = 1000, seed = 11)
res_GN <- conditional_perm_specificity(GNo, "y_N_in_E", c("N_subclass", "E_subclass"), c("slide", "size_bin", "N_layer", "E_layer"), n_perm = 1000, seed = 12)
res_IN <- conditional_perm_specificity(INo, "y_G_in_N", c("G_subclass", "N_subclass"), c("slide", "size_bin", "G_layer", "N_layer"), n_perm = 1000, seed = 13)
res_GG <- conditional_perm_specificity(GGs, "y_small_in_large", c("small_subclass", "large_subclass"), c("slide", "size_bin", "small_layer", "large_layer"), n_perm = 1000, seed = 14)
res_II <- conditional_perm_specificity(IIs, "y_small_in_large", c("small_subclass", "large_subclass"), c("slide", "size_bin", "small_layer", "large_layer"), n_perm = 1000, seed = 15)

as_GI <- partial_asym_stats(GIo, "G_fill", "E_fill", "GI")
as_GN <- partial_asym_stats(GNo, "N_fill", "E_fill", "GN")
as_IN <- partial_asym_stats(INo, "G_fill", "N_fill", "IN")
as_GG <- partial_asym_stats(GGs, "small_fill", "large_fill", "GG")
as_II <- partial_asym_stats(IIs, "small_fill", "large_fill", "II")

motif_summary <- tibble::tibble(
  motif = c("Gaba⊂Glut", "NonNeuron⊂Glut", "Gaba⊂NonNeuron", "smaller⊂Glut (GG)", "smaller⊂Gaba (II)"),
  pair_class = c("GI", "GN", "IN", "GG", "II"),
  raw_rate = c(mean(GIo$G_fill == 1), mean(GNo$N_fill == 1), mean(INo$G_fill == 1), mean(GGs$small_fill == 1), mean(IIs$small_fill == 1)),
  eligible_rate = c(
    mean(GIo$G_fill[GIo$E_n >= GIo$G_n] == 1, na.rm = TRUE),
    mean(GNo$N_fill[GNo$E_n >= GNo$N_n] == 1, na.rm = TRUE),
    mean(INo$G_fill[INo$N_n >= INo$G_n] == 1, na.rm = TRUE),
    mean(GGs$small_fill == 1, na.rm = TRUE),
    mean(IIs$small_fill == 1, na.rm = TRUE)
  ),
  size_auc = c(
    safe_auc(GIo$y_G_in_E, GIo$log_size_ratio),
    safe_auc(GNo$y_N_in_E, GNo$log_size_ratio),
    safe_auc(INo$y_G_in_N, INo$log_size_ratio),
    safe_auc(GGs$y_small_in_large, GGs$log_size_ratio),
    safe_auc(IIs$y_small_in_large, IIs$log_size_ratio)
  ),
  identity_v = c(res_GI$obs_v, res_GN$obs_v, res_IN$obs_v, res_GG$obs_v, res_II$obs_v),
  identity_perm_z = c(res_GI$z, res_GN$z, res_IN$z, res_GG$z, res_II$z),
  identity_perm_p = c(res_GI$p, res_GN$p, res_IN$p, res_GG$p, res_II$p),
  partial_asym_median = c(as_GI$pair_median, as_GN$pair_median, as_IN$pair_median, as_GG$pair_median, as_II$pair_median),
  partial_slide_p = c(as_GI$slide_wilcoxon_p, as_GN$slide_wilcoxon_p, as_IN$slide_wilcoxon_p, as_GG$slide_wilcoxon_p, as_II$slide_wilcoxon_p),
  n_pairs = c(nrow(GIo), nrow(GNo), nrow(INo), nrow(GGs), nrow(IIs))
)
write_csv(motif_summary, file.path(table_dir, "table01_motif_summary.csv"))
write_csv(bind_rows(as_GI, as_GN, as_IN, as_GG, as_II), file.path(table_dir, "table02_partial_asymmetry_summary.csv"))

# ------------------------------------------------------------------
# Size-matched GABA partner preference: Glut host vs NonNeuron host
# ------------------------------------------------------------------
gi_elig <- GIo %>%
  dplyr::filter(E_n >= G_n) %>%
  dplyr::transmute(G_label, G_subclass, G_layer, partner = "Glut", full = as.integer(G_fill == 1), size_ratio = E_n / G_n)

in_elig <- INo %>%
  dplyr::filter(N_n >= G_n) %>%
  dplyr::transmute(G_label, G_subclass, G_layer, partner = "NonNeuron", full = as.integer(G_fill == 1), size_ratio = N_n / G_n)

gcomb <- dplyr::bind_rows(gi_elig, in_elig) %>% dplyr::mutate(size_bin = dplyr::ntile(size_ratio, 10))

# subclass/layer/size-bin stratified MH pooled OR
strat_g <- gcomb %>% dplyr::mutate(stratum = interaction(G_subclass, G_layer, size_bin, drop = TRUE))
mh_tables <- split(strat_g, strat_g$stratum)
mh_tables <- mh_tables[vapply(mh_tables, function(x) length(unique(x$partner)) > 1, logical(1))]
mh_array <- array(0, dim = c(2, 2, length(mh_tables)), dimnames = list(c("Glut", "NonNeuron"), c("0", "1"), NULL))
kk <- 1
for (nm in names(mh_tables)) {
  x <- mh_tables[[nm]]
  ct <- table(x$partner, x$full)
  ct2 <- matrix(0, nrow = 2, ncol = 2, dimnames = list(c("Glut", "NonNeuron"), c("0", "1")))
  ct2[rownames(ct), colnames(ct)] <- ct
  mh_array[, , kk] <- ct2
  kk <- kk + 1
}
mh_gaba <- stats::mantelhaen.test(mh_array)

cmh_gaba_summary <- tibble::tibble(
  comparison = "Comparable GABA pairs: Glut host vs NonNeuron host",
  pooled_OR = unname(mh_gaba$estimate),
  ci_low = if (!is.null(mh_gaba$conf.int)) mh_gaba$conf.int[1] else NA_real_,
  ci_high = if (!is.null(mh_gaba$conf.int)) mh_gaba$conf.int[2] else NA_real_,
  p = mh_gaba$p.value,
  n_strata = length(mh_tables)
)
write_csv(cmh_gaba_summary, file.path(table_dir, "table03_cmh_gaba_host_preference.csv"))

# ------------------------------------------------------------------
# Same GABA cluster: Glut host vs NonNeuron host
# ------------------------------------------------------------------
rate_by_guest <- function(df, guest_col, full_vec) {
  tibble::tibble(guest = df[[guest_col]], full = full_vec) %>%
    dplyr::group_by(guest) %>%
    dplyr::summarise(rate = mean(full), n_pairs = dplyr::n(), n_full = sum(full), .groups = "drop")
}

g_glu_rate_elig <- rate_by_guest(GIo[GIo$E_n >= GIo$G_n, ], "G_label", GIo$G_fill[GIo$E_n >= GIo$G_n] == 1)
g_nn_rate_elig  <- rate_by_guest(INo[INo$N_n >= INo$G_n, ], "G_label", INo$G_fill[INo$N_n >= INo$G_n] == 1)
shared_g_elig <- intersect(g_glu_rate_elig$guest, g_nn_rate_elig$guest)

paired_rate_elig <- dplyr::inner_join(
  dplyr::rename(g_glu_rate_elig[g_glu_rate_elig$guest %in% shared_g_elig, ], rate_glu = rate, n_pairs_glu = n_pairs, n_full_glu = n_full),
  dplyr::rename(g_nn_rate_elig[g_nn_rate_elig$guest %in% shared_g_elig, ], rate_nn = rate, n_pairs_nn = n_pairs, n_full_nn = n_full),
  by = "guest"
)

paired_rate_summary <- tibble::tibble(
  n_shared_gaba = nrow(paired_rate_elig),
  median_rate_glu = median(paired_rate_elig$rate_glu, na.rm = TRUE),
  median_rate_nn = median(paired_rate_elig$rate_nn, na.rm = TRUE),
  p_paired = paired_wilcox_p(paired_rate_elig$rate_glu, paired_rate_elig$rate_nn)
)
write_csv(paired_rate_elig, file.path(table_dir, "table04_shared_gaba_eligible_rate_pairs.csv"))
write_csv(paired_rate_summary, file.path(table_dir, "table05_shared_gaba_eligible_rate_summary.csv"))

# ------------------------------------------------------------------
# Exact directed containment counts: GABA→GLU vs GLU→GABA
# ------------------------------------------------------------------
gi_dir_counts <- tibble::tibble(
  direction = c("GABA in Glut", "Glut in GABA"),
  exact_n = c(sum(GIo$G_fill == 1, na.rm = TRUE), sum(GIo$E_fill == 1, na.rm = TRUE)),
  total_pairs = nrow(GIo)
) %>%
  dplyr::mutate(rate = exact_n / total_pairs)

binom_gi <- binom.test(gi_dir_counts$exact_n[1], sum(gi_dir_counts$exact_n), p = 0.5, alternative = "greater")
write_csv(gi_dir_counts, file.path(table_dir, "table06_gi_directional_exact_counts.csv"))
write_csv(tibble::tibble(p = binom_gi$p.value), file.path(table_dir, "table07_gi_direction_binom.csv"))

# ------------------------------------------------------------------
# Priority score table
# ------------------------------------------------------------------
priority_tbl <- motif_summary %>%
  dplyr::filter(motif %in% c("Gaba⊂Glut", "smaller⊂Glut (GG)", "Gaba⊂NonNeuron")) %>%
  dplyr::mutate(
    priority_tier = c("Tier 1", "Tier 2", "Tier 3"),
    interpretation = c(
      "Best causal candidate",
      "Field-level hierarchical candidate",
      "Mostly geometric/byproduct regime"
    )
  )
write_csv(priority_tbl, file.path(table_dir, "table08_priority_tiers.csv"))

# ------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------
# Fig 1: priority map
fig1_df <- priority_tbl %>%
  dplyr::mutate(label = paste0(motif, "\nAUC=", fmt_num(size_auc, 3), "; z=", fmt_num(identity_perm_z, 2)))

p1 <- ggplot(fig1_df, aes(x = size_auc, y = identity_perm_z, size = abs(partial_asym_median), colour = motif)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.6, colour = "#888888") +
  geom_vline(xintercept = 0.75, linetype = 2, linewidth = 0.6, colour = "#BBBBBB") +
  geom_point(alpha = 0.95) +
  geom_text(aes(label = label), nudge_x = 0.015, hjust = 0, size = 3.4, show.legend = FALSE) +
  scale_colour_manual(values = col_motif[names(col_motif) %in% fig1_df$motif]) +
  scale_size_continuous(range = c(4.5, 10)) +
  labs(
    title = wrap_title("Containment motifs separate into causal-priority tiers"),
    subtitle = wrap_subtitle("x-axis: size-only discriminability (AUC). y-axis: residual identity signal after conditioning on slide/layer/size bin. Point size scales with precursor partial asymmetry."),
    x = "Size-only AUC for exact containment",
    y = "Residual identity permutation z-score",
    colour = NULL,
    size = "|Partial asym|"
  ) +
  base_theme()
save_pub(p1, "fig1_priority_map.png", width = 8.8, height = 5.8)

# Fig 2: size-controlled OR forest
fig2_df <- cmh_gaba_summary %>%
  dplyr::mutate(
    comparison = factor(comparison, levels = comparison),
    label = paste0("OR=", fmt_num(pooled_OR, 2), " [", fmt_num(ci_low, 2), ", ", fmt_num(ci_high, 2), "]\np=", fmt_p(p))
  )

p2 <- ggplot(fig2_df, aes(x = pooled_OR, y = comparison)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.6, colour = "#888888") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.18, linewidth = 0.9, colour = col_motif[["Gaba⊂Glut"]]) +
  geom_point(size = 3.8, colour = col_motif[["Gaba⊂Glut"]]) +
  geom_text(aes(x = ci_high * 1.05, label = label), hjust = 0, size = 3.6, colour = "#444444") +
  scale_x_log10(labels = scales::label_number(accuracy = 0.01), expand = expansion(mult = c(0.03, 0.34))) +
  coord_cartesian(clip = "off") +
  labs(
    title = wrap_title("Gaba⊂Glut remains enriched after size matching"),
    subtitle = wrap_subtitle("Only comparable pairs are kept, and the pooled odds ratio is stratified by GABA subclass, layer, and host/guest size bin."),
    x = "Mantel–Haenszel pooled odds ratio (log scale)",
    y = NULL
  ) +
  base_theme()
save_pub(p2, "fig2_size_control_cmh_forest.png", width = 8.4, height = 4.8)

# Fig 3: residual identity permutation comparison
fig3_df <- priority_tbl %>%
  dplyr::mutate(
    motif = factor(motif, levels = c("Gaba⊂Glut", "smaller⊂Glut (GG)", "Gaba⊂NonNeuron")),
    stat_lab = paste0("z=", fmt_num(identity_perm_z, 2), "\np=", purrr::map_chr(identity_perm_p, fmt_p))
  )

p3 <- ggplot(fig3_df, aes(x = motif, y = identity_perm_z, fill = motif)) +
  geom_col(width = 0.66, alpha = 0.96) +
  geom_text(aes(label = stat_lab), vjust = -0.28, size = 3.4) +
  scale_fill_manual(values = col_motif[names(col_motif) %in% levels(fig3_df$motif)]) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = wrap_title("Residual identity signal distinguishes causal motifs from geometric byproducts"),
    subtitle = wrap_subtitle("Permutation is conditioned on slide, layer, and size bin. Gaba⊂Glut and smaller⊂Glut retain significant subclass-specific structure, whereas Gaba⊂NonNeuron does not."),
    x = NULL,
    y = "Permutation z-score for residual identity specificity"
  ) +
  base_theme() +
  theme(legend.position = "none")
save_pub(p3, "fig3_residual_identity_barplot.png", width = 8.4, height = 5.0)

# Fig 4: precursor partial asymmetry
fig4_df <- bind_rows(as_GI, as_GG, as_IN) %>%
  dplyr::mutate(
    motif = factor(c("Gaba⊂Glut", "smaller⊂Glut (GG)", "Gaba⊂NonNeuron"), levels = c("Gaba⊂Glut", "smaller⊂Glut (GG)", "Gaba⊂NonNeuron")),
    stat_lab = paste0("median=", fmt_num(slide_median_of_medians, 3), "\n", slide_positive_n, "/", slide_n, " slides > 0\np=", purrr::map_chr(slide_wilcoxon_p, fmt_p))
  )

p4 <- ggplot(fig4_df, aes(x = motif, y = slide_median_of_medians, fill = motif)) +
  geom_col(width = 0.66, alpha = 0.96) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.6, colour = "#888888") +
  geom_text(aes(label = stat_lab), vjust = ifelse(fig4_df$slide_median_of_medians >= 0, -0.25, 1.15), size = 3.3) +
  scale_fill_manual(values = col_motif[names(col_motif) %in% levels(fig4_df$motif)]) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.20))) +
  labs(
    title = wrap_title("Partial-overlap asymmetry appears early for Gaba⊂Glut and smaller⊂Glut, but not for Gaba⊂NonNeuron"),
    subtitle = wrap_subtitle("Slide-level medians are tested against zero. Positive values indicate that the candidate guest already overlaps more asymmetrically at the partial stage."),
    x = NULL,
    y = "Slide-level median signed asymmetry"
  ) +
  base_theme() +
  theme(legend.position = "none")
save_pub(p4, "fig4_partial_precursor_barplot.png", width = 9.0, height = 6.0)

# Fig 5: paired shared-GABA host preference
fig5_long <- paired_rate_elig %>%
  dplyr::transmute(
    guest,
    `Glut host` = rate_glu,
    `NonNeuron host` = rate_nn
  ) %>%
  tidyr::pivot_longer(-guest, names_to = "host_class", values_to = "rate") %>%
  dplyr::mutate(host_class = factor(host_class, levels = c("NonNeuron host", "Glut host")))

fig5_stat <- tibble::tibble(
  x = 1.5,
  y = max(fig5_long$rate, na.rm = TRUE) * 1.07,
  label = paste0(
    "n=", scales::comma(nrow(paired_rate_elig)),
    "; median Glut=", fmt_num(median(paired_rate_elig$rate_glu, na.rm = TRUE), 2),
    "; median NonNeuron=", fmt_num(median(paired_rate_elig$rate_nn, na.rm = TRUE), 2),
    "\np=", fmt_p(paired_rate_summary$p_paired)
  )
)

p5 <- ggplot(fig5_long, aes(x = host_class, y = rate, group = guest)) +
  geom_line(colour = "#CFCFCF", alpha = 0.28, linewidth = 0.35) +
  geom_point(aes(colour = host_class), alpha = 0.80, size = 1.8) +
  geom_boxplot(aes(fill = host_class), width = 0.34, outlier.shape = NA, alpha = 0.90, linewidth = 0.8) +
  geom_text(data = fig5_stat, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 3.3) +
  scale_colour_manual(values = col_partner) +
  scale_fill_manual(values = col_partner) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0.02, 0.16))) +
  labs(
    title = wrap_title("Within the same GABA cluster, Glut behaves more like the true exact host than NonNeuron"),
    subtitle = wrap_subtitle("Only GABA clusters with comparable Glut and NonNeuron partners are kept. Each line is one shared GABA cluster."),
    x = NULL,
    y = "Exact containment rate for the GABA cluster",
    colour = NULL,
    fill = NULL
  ) +
  base_theme()
save_pub(p5, "fig5_shared_gaba_paired_rates.png", width = 8.6, height = 5.6)

# Fig 6: directional exact count within GI
fig6_df <- gi_dir_counts %>%
  dplyr::mutate(direction = factor(direction, levels = c("GABA in Glut", "Glut in GABA")))
sub6 <- paste0("Among ", scales::comma(nrow(GIo)), " unique GABA–Glut pairs: ",
               "directional OR=", fmt_num((fig6_df$exact_n[1] / (nrow(GIo) - fig6_df$exact_n[1])) / (fig6_df$exact_n[2] / (nrow(GIo) - fig6_df$exact_n[2])), 2),
               "; binomial p=", fmt_p(binom_gi$p.value), ".")

p6 <- ggplot(fig6_df, aes(x = direction, y = exact_n, fill = direction)) +
  geom_col(width = 0.62, alpha = 0.96) +
  geom_text(aes(label = paste0(scales::comma(exact_n), "\n(", fmt_pct(rate, 1), ")")), vjust = -0.15, size = 3.8) +
  scale_fill_manual(values = c("GABA in Glut" = col_motif[["Gaba⊂Glut"]], "Glut in GABA" = "#4E79A7")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(
    title = wrap_title("Directionality inside GABA–Glut overlap is strongly biased toward Gaba⊂Glut"),
    subtitle = wrap_subtitle(sub6),
    x = NULL,
    y = "Number of exact directional edges"
  ) +
  base_theme() +
  theme(legend.position = "none")
save_pub(p6, "fig6_gi_directional_exact_counts.png", width = 8.0, height = 5.2)

# Fig 7: motif-specific comparison dashboard
fig7_df <- priority_tbl %>%
  dplyr::transmute(
    motif,
    `Residual identity z` = identity_perm_z,
    `Partial precursor` = partial_asym_median,
    `Size-only AUC` = size_auc
  ) %>%
  tidyr::pivot_longer(-motif, names_to = "metric", values_to = "value") %>%
  dplyr::group_by(metric) %>%
  dplyr::mutate(zv = as.numeric(scale(value))) %>%
  dplyr::ungroup()

p7 <- ggplot(fig7_df, aes(x = metric, y = fct_reorder(motif, zv, .fun = mean), fill = zv)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f", value)), size = 3.4) +
  scale_fill_gradient2(low = "#4E79A7", mid = "#F2F2F2", high = "#C44E52", midpoint = 0) +
  labs(
    title = wrap_title("The three priority tiers separate consistently across orthogonal evidence axes"),
    subtitle = wrap_subtitle("Values are shown in their native scales; tile colour reflects z-scored ranking within each metric."),
    x = NULL,
    y = NULL,
    fill = "Within-metric\nrank z"
  ) +
  base_theme()
save_pub(p7, "fig7_priority_dashboard_heatmap.png", width = 8.6, height = 4.6)

# Fig 8: concise manuscript-style schematic summary
fig8_df <- tibble::tibble(
  tier = factor(c("Tier 1", "Tier 2", "Tier 3"), levels = c("Tier 1", "Tier 2", "Tier 3")),
  motif = c("Gaba⊂Glut", "smaller⊂Glut (GG)", "Gaba⊂NonNeuron"),
  y = c(3, 2, 1),
  key_result = c(
    paste0("Size-matched OR=", fmt_num(cmh_gaba_summary$pooled_OR, 2), "; residual z=", fmt_num(res_GI$z, 2), "; partial med=", fmt_num(as_GI$slide_median_of_medians, 3)),
    paste0("Residual z=", fmt_num(res_GG$z, 2), "; partial med=", fmt_num(as_GG$slide_median_of_medians, 3)),
    paste0("AUC=", fmt_num(motif_summary$size_auc[motif_summary$motif == "Gaba⊂NonNeuron"], 3), "; residual p=", fmt_p(res_IN$p), "; partial p=", fmt_p(as_IN$slide_wilcoxon_p))
  ),
  interpretation = c(
    "Best causal motif",
    "Hierarchical excitatory field rule",
    "Mostly geometric/support-shell byproduct"
  )
)

p8 <- ggplot(fig8_df, aes(x = tier, y = y, colour = motif)) +
  geom_point(size = 6) +
  geom_text(aes(label = motif), nudge_y = 0.28, size = 3.8, show.legend = FALSE) +
  geom_text(aes(label = key_result), nudge_y = -0.18, size = 3.1, show.legend = FALSE, colour = "#444444") +
  geom_text(aes(label = interpretation), nudge_y = -0.47, size = 3.2, show.legend = FALSE, colour = "#666666") +
  scale_colour_manual(values = col_motif[names(col_motif) %in% fig8_df$motif]) +
  scale_y_continuous(limits = c(0.4, 3.5), breaks = NULL) +
  labs(
    title = wrap_title("Causal-priority interpretation of containment motifs"),
    subtitle = wrap_subtitle("Tier 1 is the best candidate for direct causal experiments, Tier 2 is a field-level hierarchy rule, and Tier 3 is best interpreted as a geometric/support-system outcome."),
    x = NULL,
    y = NULL,
    colour = NULL
  ) +
  base_theme()
save_pub(p8, "fig8_priority_schematic_summary.png", width = 9.5, height = 4.9)

# ------------------------------------------------------------------
# Text summary for quick manuscript/PPT reuse
# ------------------------------------------------------------------
summary_txt <- c(
  paste0("Tier 1 candidate: Gaba⊂Glut. Size-matched pooled OR = ", fmt_num(cmh_gaba_summary$pooled_OR, 2),
         ", 95% CI [", fmt_num(cmh_gaba_summary$ci_low, 2), ", ", fmt_num(cmh_gaba_summary$ci_high, 2),
         "], p = ", fmt_p(cmh_gaba_summary$p), ". Residual identity permutation z = ", fmt_num(res_GI$z, 2),
         ", p = ", fmt_p(res_GI$p), ". Partial asymmetry median = ", fmt_num(as_GI$slide_median_of_medians, 3),
         ", slide-level p = ", fmt_p(as_GI$slide_wilcoxon_p), "."),
  paste0("Tier 2 candidate: smaller⊂larger Glut. Residual identity permutation z = ", fmt_num(res_GG$z, 2),
         ", p = ", fmt_p(res_GG$p), ". Partial asymmetry median = ", fmt_num(as_GG$slide_median_of_medians, 3),
         ", slide-level p = ", fmt_p(as_GG$slide_wilcoxon_p), "."),
  paste0("Tier 3 regime: Gaba⊂NonNeuron. Size-only AUC = ", fmt_num(motif_summary$size_auc[motif_summary$motif == "Gaba⊂NonNeuron"], 3),
         ". Residual identity permutation z = ", fmt_num(res_IN$z, 2), ", p = ", fmt_p(res_IN$p),
         ". Partial asymmetry median = ", fmt_num(as_IN$slide_median_of_medians, 3),
         ", slide-level p = ", fmt_p(as_IN$slide_wilcoxon_p), ".")
)
writeLines(summary_txt, con = file.path(out_dir, "README_R_output_11_summary.txt"))

message("R_output_11 finished: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
