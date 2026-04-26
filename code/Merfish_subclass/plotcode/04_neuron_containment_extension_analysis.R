# 清空环境
rm(list = ls())
gc()
# =========================================================
# neuron_containment_extension_analysis_topjournal.R
# Top-journal visualization version
# Input : E:/zaw/2603/neuron-neuron-partner.csv
# Output: E:/zaw/2603/R_output_04
# =========================================================

required_pkgs <- c(
  "data.table", "dplyr", "tidyr", "purrr", "stringr",
  "ggplot2", "readr", "tibble", "sandwich", "lmtest", "scales"
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
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(sandwich)
  library(lmtest)
  library(scales)
})

options(scipen = 999)


# ---------------------------
# Protect dplyr verbs from method-dispatch conflicts
# ---------------------------
select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter
mutate <- dplyr::mutate
summarise <- dplyr::summarise
arrange <- dplyr::arrange
count <- dplyr::count
distinct <- dplyr::distinct
left_join <- dplyr::left_join
inner_join <- dplyr::inner_join
pull <- dplyr::pull
transmute <- dplyr::transmute
first <- dplyr::first


# ---------------------------
# Paths
# ---------------------------
base_dir <- "E:/zaw/2603"
input_csv <- file.path(base_dir, "neuron-neuron-partner.csv")
out_dir <- file.path(base_dir, "R_output_04")
fig_dir <- file.path(out_dir, "figures")
tbl_dir <- file.path(out_dir, "tables")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tbl_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------
# Helpers
# ---------------------------
canonical_pair <- function(a, b) {
  ifelse(a < b, paste(a, b, sep = " || "), paste(b, a, sep = " || "))
}

clean_subclass <- function(x) {
  stringr::str_replace(x, "^\\d+\\s+", "")
}

rank_biserial_paired <- function(x, y) {
  d <- x - y
  d <- d[d != 0]
  if (length(d) == 0) return(NA_real_)
  (sum(d > 0) - sum(d < 0)) / length(d)
}

paired_wilcox_from_wide <- function(dat, a, b) {
  # 修复处：显式调用 dplyr::select 防止被 AnnotationDbi/MASS 覆盖
  sub <- dat %>% dplyr::select(all_of(c(a, b))) %>% tidyr::drop_na()
  if (nrow(sub) == 0) return(NULL)
  wt <- suppressWarnings(wilcox.test(sub[[a]], sub[[b]], paired = TRUE, exact = FALSE))
  tibble(
    n = nrow(sub),
    mean_diff = mean(sub[[a]] - sub[[b]]),
    median_diff = median(sub[[a]] - sub[[b]]),
    rank_biserial = rank_biserial_paired(sub[[a]], sub[[b]]),
    p_value = wt$p.value
  )
}

robust_coef_table <- function(fit, cluster, exponentiate = FALSE) {
  V <- sandwich::vcovCL(fit, cluster = cluster)
  ct <- lmtest::coeftest(fit, vcov. = V)
  out <- data.frame(
    term = rownames(ct),
    estimate = as.numeric(ct[, 1]),
    std.error = as.numeric(ct[, 2]),
    statistic = as.numeric(ct[, 3]),
    p.value = as.numeric(ct[, 4]),
    row.names = NULL,
    check.names = FALSE
  )
  out <- out %>%
    mutate(
      conf.low = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error
    )
  if (exponentiate) {
    out <- out %>%
      mutate(
        OR = exp(estimate),
        OR_low = exp(conf.low),
        OR_high = exp(conf.high)
      )
  }
  out
}

fisher_enrichment <- function(category_vec, outcome_vec, min_n = 20) {
  # 修复：强制将输入转换为原子向量，防止 list-column 报错
  category_vec <- as.character(unlist(category_vec))
  outcome_vec <- as.logical(unlist(outcome_vec))
  
  stopifnot(length(category_vec) == length(outcome_vec))
  
  dat <- tibble(category = category_vec, outcome = outcome_vec)
  total_pos <- sum(dat$outcome)
  total_neg <- sum(!dat$outcome)
  
  dat %>%
    dplyr::count(category, name = "n_edges") %>%
    filter(n_edges >= min_n) %>%
    mutate(
      n_exact = map_int(category, ~sum(dat$outcome[dat$category == .x])),
      rate_exact = n_exact / n_edges,
      odds_ratio = map_dbl(category, ~{
        pos <- sum(dat$outcome[dat$category == .x])
        neg <- sum(!dat$outcome[dat$category == .x])
        other_pos <- total_pos - pos
        other_neg <- total_neg - neg
        unname(fisher.test(matrix(c(pos, neg, other_pos, other_neg), nrow = 2, byrow = TRUE))$estimate)
      }),
      p_value = map_dbl(category, ~{
        pos <- sum(dat$outcome[dat$category == .x])
        neg <- sum(!dat$outcome[dat$category == .x])
        other_pos <- total_pos - pos
        other_neg <- total_neg - neg
        fisher.test(matrix(c(pos, neg, other_pos, other_neg), nrow = 2, byrow = TRUE))$p.value
      })
    ) %>%
    mutate(
      fdr_bh = p.adjust(p_value, method = "BH"),
      direction = ifelse(odds_ratio > 1, "enriched", "depleted")
    ) %>%
    arrange(fdr_bh, desc(odds_ratio))
}

mw_test <- function(df, group_col, metric) {
  g1 <- df %>% filter(.data[[group_col]]) %>% pull(.data[[metric]]) %>% na.omit()
  g0 <- df %>% filter(!.data[[group_col]]) %>% pull(.data[[metric]]) %>% na.omit()
  if (length(g1) == 0 || length(g0) == 0) return(NULL)
  wt <- wilcox.test(g1, g0, exact = FALSE)
  tibble(
    metric = metric,
    n_group1 = length(g1),
    n_group0 = length(g0),
    mean_group1 = mean(g1),
    mean_group0 = mean(g0),
    median_group1 = median(g1),
    median_group0 = median(g0),
    mean_diff = mean(g1) - mean(g0),
    p_value = wt$p.value
  )
}

fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 1e-300, "<1e-300",
                ifelse(p < 1e-4, format(p, scientific = TRUE, digits = 3),
                       sprintf("%.6f", p))))
}

fmt_p_short <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-300) return("<1e-300")
  if (p < 1e-4) return(format(p, scientific = TRUE, digits = 2))
  sprintf("%.4f", p)
}

save_pub <- function(p, filename, width, height, dpi = 600) {
  ggsave(file.path(fig_dir, paste0(filename, ".png")), p,
         width = width, height = height, dpi = dpi, bg = "white")
  ggsave(file.path(fig_dir, paste0(filename, ".pdf")), p,
         width = width, height = height, bg = "white")
}

add_sig_bracket <- function(df, x1, x2, y, label) {
  bind_rows(df,
            tibble(x = x1, xend = x2, y = y, yend = y, label = label, xlab = (x1 + x2) / 2, ylab = y),
            tibble(x = x1, xend = x1, y = y, yend = y * 0.93, label = NA_character_, xlab = NA_real_, ylab = NA_real_),
            tibble(x = x2, xend = x2, y = y, yend = y * 0.93, label = NA_character_, xlab = NA_real_, ylab = NA_real_))
}

# ---------------------------
# Plot style
# ---------------------------
base_family <- "sans"
col_pair <- c("E-E" = "#4E79A7", "E-I" = "#A06AB4", "I-I" = "#59A14F")
col_dir <- c("GABA in Glu" = "#C44E52", "Glu in GABA" = "#4E79A7")
col_small <- c("Smaller cluster = GABA" = "#C44E52", "Smaller cluster = Glu" = "#4E79A7")
col_binary <- c("Exact GABA-in-Glu" = "#C44E52", "Other E-I" = "#BDBDBD")
col_layer <- c("#E8D9EE", "#D3B8E0", "#BF95D1", "#AA73C3", "#954FB4", "#7F2EA5")
col_multi <- c("GABA guest in Glu" = "#C44E52", "Glu guest in GABA" = "#4E79A7")

pub_theme <- function() {
  theme_classic(base_size = 12, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0),
      plot.subtitle = element_text(size = 11.5, colour = "#4D4D4D", margin = margin(b = 8)),
      axis.title = element_text(size = 14, colour = "black"),
      axis.text = element_text(size = 12, colour = "black"),
      axis.line = element_line(linewidth = 0.7, colour = "black"),
      axis.ticks = element_line(linewidth = 0.7, colour = "black"),
      axis.ticks.length = unit(0.18, "cm"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 11),
      panel.grid = element_blank(),
      plot.margin = margin(12, 16, 12, 12)
    )
}

theme_set(pub_theme())

# ---------------------------
# Read data
# ---------------------------
df <- fread(input_csv)
df[, pair_key := canonical_pair(cluster.1_label, cluster.2_label)]

# Unique clusters
c1 <- copy(df[, .(
  label = cluster.1_label,
  slide = cluster.1_slide,
  layer = cluster.1_layer,
  region = cluster.1_region,
  cell_Neruon_type = cluster.1_cell_Neruon_type,
  subclass = cluster.1_subclass,
  total_cell_num = cluster.1_total_cell_num,
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
  cauchy_combination_p = cluster.2_cauchy_combination_p,
  E_I_Ratio = cluster.2_E_I_Ratio
)])

clusters <- unique(rbindlist(list(c1, c2), use.names = TRUE))
clusters[, block := paste(slide, layer, sep = "|")]
clusters[, subclass_clean := clean_subclass(subclass)]

# Unique undirected pairs
pair_list <- lapply(split(df, df$pair_key), function(sub) {
  stopifnot(nrow(sub) == 2)
  labs <- sort(unique(c(sub$cluster.1_label, sub$cluster.2_label)))
  a <- labs[1]
  b <- labs[2]
  ra <- sub[sub$cluster.1_label == a, ][1, ]
  tibble(
    label_a = a,
    label_b = b,
    slide = ra$cluster.1_slide,
    layer = ra$cluster.1_layer,
    region_a = ra$cluster.1_region,
    region_b = ra$cluster.2_region,
    type_a_raw = ra$cluster.1_cell_Neruon_type,
    type_b_raw = ra$cluster.2_cell_Neruon_type,
    subclass_a = ra$cluster.1_subclass,
    subclass_b = ra$cluster.2_subclass,
    subclass_a_clean = clean_subclass(ra$cluster.1_subclass),
    subclass_b_clean = clean_subclass(ra$cluster.2_subclass),
    size_a = ra$cluster.1_total_cell_num,
    size_b = ra$cluster.2_total_cell_num,
    ov_a = ra$cluster.1.overlap.percent,
    ov_b = ra$cluster.2.overlap.percent,
    overlap = ra$overlap_cell,
    union = ra$union_cell,
    jaccard = ra$jaccard,
    p_a = ra$cluster.1_cauchy_combination_p,
    p_b = ra$cluster.2_cauchy_combination_p,
    ei_ratio_a = ra$cluster.1_E_I_Ratio,
    ei_ratio_b = ra$cluster.2_E_I_Ratio
  )
})

pairs <- bind_rows(pair_list) %>%
  mutate(
    pair_type = map2_chr(type_a_raw, type_b_raw, ~paste(sort(c(ifelse(.x == "Glut", "E", "I"), ifelse(.y == "Glut", "E", "I"))), collapse = "-")),
    block = paste(slide, layer, sep = "|"),
    size_small = pmin(size_a, size_b),
    size_large = pmax(size_a, size_b),
    overlap_prop_small = overlap / size_small,
    overlap_prop_large = overlap / size_large,
    exact_smaller_contained = abs(overlap_prop_small - 1) < 1e-12,
    near95_smaller_contained = overlap_prop_small >= 0.95
  )

# Oriented E-I table
ei <- pairs %>%
  filter(pair_type == "E-I") %>%
  mutate(
    glu_label = ifelse(type_a_raw == "Glut", label_a, label_b),
    gaba_label = ifelse(type_a_raw == "Glut", label_b, label_a),
    glu_subclass = ifelse(type_a_raw == "Glut", subclass_a_clean, subclass_b_clean),
    gaba_subclass = ifelse(type_a_raw == "Glut", subclass_b_clean, subclass_a_clean),
    glu_size = ifelse(type_a_raw == "Glut", size_a, size_b),
    gaba_size = ifelse(type_a_raw == "Glut", size_b, size_a),
    glu_overlap_prop = ifelse(type_a_raw == "Glut", ov_a, ov_b),
    gaba_overlap_prop = ifelse(type_a_raw == "Glut", ov_b, ov_a),
    glu_region = ifelse(type_a_raw == "Glut", region_a, region_b),
    gaba_region = ifelse(type_a_raw == "Glut", region_b, region_a),
    glu_p = ifelse(type_a_raw == "Glut", p_a, p_b),
    gaba_p = ifelse(type_a_raw == "Glut", p_b, p_a),
    glu_ei_ratio = ifelse(type_a_raw == "Glut", ei_ratio_a, ei_ratio_b),
    gaba_ei_ratio = ifelse(type_a_raw == "Glut", ei_ratio_b, ei_ratio_a),
    same_region = glu_region == gaba_region,
    gaba_in_glu_exact = abs(gaba_overlap_prop - 1) < 1e-12,
    glu_in_gaba_exact = abs(glu_overlap_prop - 1) < 1e-12,
    gaba_in_glu_99 = gaba_overlap_prop >= 0.99,
    gaba_in_glu_95 = gaba_overlap_prop >= 0.95,
    glu_in_gaba_99 = glu_overlap_prop >= 0.99,
    glu_in_gaba_95 = glu_overlap_prop >= 0.95,
    size_ratio_glu_to_gaba = glu_size / gaba_size,
    log_size_ratio_glu_to_gaba = log(size_ratio_glu_to_gaba),
    gaba_is_smaller = gaba_size < glu_size,
    glu_is_smaller = glu_size < gaba_size,
    same_size = glu_size == gaba_size,
    smaller_is_gaba = as.integer(gaba_is_smaller),
    smaller_contained_exact = case_when(
      gaba_is_smaller ~ gaba_in_glu_exact,
      glu_is_smaller ~ glu_in_gaba_exact,
      TRUE ~ FALSE
    ),
    smaller_contained_95 = case_when(
      gaba_is_smaller ~ gaba_in_glu_95,
      glu_is_smaller ~ glu_in_gaba_95,
      TRUE ~ FALSE
    ),
    smaller_type = case_when(
      gaba_is_smaller ~ "GABA_smaller",
      glu_is_smaller ~ "GLU_smaller",
      TRUE ~ "same_size"
    ),
    subpair_clean = paste(glu_subclass, gaba_subclass, sep = " × "),
    jaccard_clip = pmin(pmax(jaccard, 1e-4), 1 - 1e-4),
    logit_jaccard = log(jaccard_clip / (1 - jaccard_clip))
  )

# ---------------------------
# Tables
# ---------------------------
overview <- tibble(
  metric = c(
    "Raw directional rows",
    "Unique undirected pairs",
    "E-I undirected pairs",
    "Exact GABA in Glu",
    "Exact Glu in GABA",
    "Near-exact GABA in Glu (>=99%)",
    "Near-exact GABA in Glu (>=95%)",
    "E-I with smaller cluster = GABA",
    "E-I with smaller cluster = Glu",
    "Exact smaller-cluster containment rate when smaller = GABA",
    "Exact smaller-cluster containment rate when smaller = Glu",
    ">=95% smaller-cluster containment rate when smaller = GABA",
    ">=95% smaller-cluster containment rate when smaller = Glu"
  ),
  value = c(
    nrow(df),
    nrow(pairs),
    nrow(ei),
    sum(ei$gaba_in_glu_exact),
    sum(ei$glu_in_gaba_exact),
    sum(ei$gaba_in_glu_99),
    sum(ei$gaba_in_glu_95),
    sum(ei$gaba_is_smaller),
    sum(ei$glu_is_smaller),
    mean(ei$gaba_in_glu_exact[ei$gaba_is_smaller]),
    mean(ei$glu_in_gaba_exact[ei$glu_is_smaller]),
    mean(ei$gaba_in_glu_95[ei$gaba_is_smaller]),
    mean(ei$glu_in_gaba_95[ei$glu_is_smaller])
  )
)

pairtype_containment <- pairs %>%
  group_by(pair_type) %>%
  summarise(
    n_pairs = n(),
    exact_smaller_contained_n = sum(exact_smaller_contained),
    exact_smaller_contained_rate = mean(exact_smaller_contained),
    near95_smaller_contained_n = sum(near95_smaller_contained),
    near95_smaller_contained_rate = mean(near95_smaller_contained),
    mean_jaccard = mean(jaccard),
    median_jaccard = median(jaccard),
    .groups = "drop"
  )

n_exact_dir <- sum(ei$gaba_in_glu_exact) + sum(ei$glu_in_gaba_exact)
k_gaba_dir <- sum(ei$gaba_in_glu_exact)
binom_exact <- binom.test(k_gaba_dir, n_exact_dir, p = 0.5, alternative = "greater")

tab_smaller_exact <- matrix(
  c(
    sum(ei$gaba_in_glu_exact[ei$gaba_is_smaller]),
    sum(!ei$gaba_in_glu_exact[ei$gaba_is_smaller]),
    sum(ei$glu_in_gaba_exact[ei$glu_is_smaller]),
    sum(!ei$glu_in_gaba_exact[ei$glu_is_smaller])
  ),
  nrow = 2, byrow = TRUE
)
fisher_smaller_exact <- fisher.test(tab_smaller_exact)

tab_smaller_95 <- matrix(
  c(
    sum(ei$gaba_in_glu_95[ei$gaba_is_smaller]),
    sum(!ei$gaba_in_glu_95[ei$gaba_is_smaller]),
    sum(ei$glu_in_gaba_95[ei$glu_is_smaller]),
    sum(!ei$glu_in_gaba_95[ei$glu_is_smaller])
  ),
  nrow = 2, byrow = TRUE
)
fisher_smaller_95 <- fisher.test(tab_smaller_95)

direction_tests <- tibble(
  test = c(
    "Binomial among exact-direction E-I edges: GABA_in_Glu vs Glu_in_GABA",
    "Fisher: exact containment of smaller cluster, smaller GABA vs smaller Glu",
    "Fisher: >=95% containment of smaller cluster, smaller GABA vs smaller Glu"
  ),
  n_edges = c(n_exact_dir, sum(tab_smaller_exact), sum(tab_smaller_95)),
  k_gaba_in_glu = c(k_gaba_dir, sum(ei$gaba_in_glu_exact[ei$gaba_is_smaller]), sum(ei$gaba_in_glu_95[ei$gaba_is_smaller])),
  proportion_gaba_in_glu = c(k_gaba_dir / n_exact_dir, mean(ei$gaba_in_glu_exact[ei$gaba_is_smaller]), mean(ei$gaba_in_glu_95[ei$gaba_is_smaller])),
  odds_ratio = c(NA_real_, unname(fisher_smaller_exact$estimate), unname(fisher_smaller_95$estimate)),
  p_value = c(binom_exact$p.value, fisher_smaller_exact$p.value, fisher_smaller_95$p.value)
)

tmp <- ei %>%
  filter(gaba_is_smaller | glu_is_smaller) %>%
  mutate(
    smaller_contained_exact_int = as.integer(smaller_contained_exact),
    smaller_contained_95_int = as.integer(smaller_contained_95),
    abs_size_ratio = ifelse(gaba_is_smaller, glu_size / gaba_size, gaba_size / glu_size),
    log_abs_size_ratio = log(abs_size_ratio)
  )

fit_smaller_exact <- glm(smaller_contained_exact_int ~ smaller_is_gaba + log_abs_size_ratio + layer,
                         data = tmp, family = binomial())
fit_smaller_95 <- glm(smaller_contained_95_int ~ smaller_is_gaba + log_abs_size_ratio + layer,
                      data = tmp, family = binomial())

smaller_exact_logit <- robust_coef_table(fit_smaller_exact, cluster = tmp$block, exponentiate = TRUE)
smaller_95_logit <- robust_coef_table(fit_smaller_95, cluster = tmp$block, exponentiate = TRUE)

exact_vs_other_desc <- ei %>%
  mutate(group = ifelse(gaba_in_glu_exact, "GABA_in_Glu_exact", "Other_EI")) %>%
  group_by(group) %>%
  summarise(
    n_edges = n(),
    mean_jaccard = mean(jaccard),
    median_jaccard = median(jaccard),
    mean_size_ratio_glu_to_gaba = mean(size_ratio_glu_to_gaba),
    median_size_ratio_glu_to_gaba = median(size_ratio_glu_to_gaba),
    mean_glu_size = mean(glu_size),
    median_glu_size = median(glu_size),
    mean_gaba_size = mean(gaba_size),
    median_gaba_size = median(gaba_size),
    mean_overlap = mean(overlap),
    median_overlap = median(overlap),
    same_region_rate = mean(same_region),
    mean_glu_overlap_prop = mean(glu_overlap_prop),
    mean_gaba_overlap_prop = mean(gaba_overlap_prop),
    .groups = "drop"
  )

metric_list <- c("jaccard", "size_ratio_glu_to_gaba", "glu_size", "gaba_size",
                 "overlap", "same_region", "glu_overlap_prop", "gaba_overlap_prop",
                 "glu_ei_ratio", "gaba_ei_ratio")

block_metric_tests <- purrr::map_dfr(metric_list, function(metric) {
  tmp_metric <- ei %>%
    mutate(group = ifelse(gaba_in_glu_exact, "exact", "other")) %>%
    group_by(block, group) %>%
    summarise(value = mean(.data[[metric]], na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = group, values_from = value) %>%
    drop_na(exact, other)
  if (nrow(tmp_metric) < 10) return(NULL)
  wt <- suppressWarnings(wilcox.test(tmp_metric$exact, tmp_metric$other, paired = TRUE, exact = FALSE))
  tibble(
    metric = metric,
    n_blocks = nrow(tmp_metric),
    mean_diff_exact_minus_other = mean(tmp_metric$exact - tmp_metric$other),
    median_diff_exact_minus_other = median(tmp_metric$exact - tmp_metric$other),
    rank_biserial = rank_biserial_paired(tmp_metric$exact, tmp_metric$other),
    p_value = wt$p.value
  )
}) %>%
  mutate(fdr_bh = p.adjust(p_value, method = "BH"))

fit_jaccard <- lm(logit_jaccard ~ as.integer(gaba_in_glu_exact) + log_size_ratio_glu_to_gaba + layer, data = ei)
jaccard_regression <- robust_coef_table(fit_jaccard, cluster = ei$block, exponentiate = FALSE)

layer_rate <- ei %>%
  group_by(layer) %>%
  summarise(
    n_edges = n(),
    exact_n = sum(gaba_in_glu_exact),
    exact_rate = mean(gaba_in_glu_exact),
    near95_n = sum(gaba_in_glu_95),
    near95_rate = mean(gaba_in_glu_95),
    mean_size_ratio = mean(size_ratio_glu_to_gaba),
    .groups = "drop"
  )

layer_tests <- bind_rows(lapply(sort(unique(ei$layer)), function(lay) {
  in_layer <- ei$layer == lay
  out <- list()
  for (outcome in c("gaba_in_glu_exact", "gaba_in_glu_95")) {
    a <- sum(ei[[outcome]][in_layer])
    b <- sum(!ei[[outcome]][in_layer])
    c <- sum(ei[[outcome]][!in_layer])
    d <- sum(!ei[[outcome]][!in_layer])
    ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))
    out[[outcome]] <- tibble(
      layer = lay,
      outcome = ifelse(outcome == "gaba_in_glu_exact", "exact", "near95"),
      n_in_layer = sum(in_layer),
      rate_in_layer = a / sum(in_layer),
      odds_ratio = unname(ft$estimate),
      p_value = ft$p.value
    )
  }
  bind_rows(out)
})) %>%
  mutate(fdr_bh = p.adjust(p_value, method = "BH"))

fit_layer_exact <- glm(as.integer(gaba_in_glu_exact) ~ log_size_ratio_glu_to_gaba + layer,
                       data = ei, family = binomial())
layer_exact_logit <- robust_coef_table(fit_layer_exact, cluster = ei$block, exponentiate = TRUE)

block_exact <- pairs %>%
  group_by(block, pair_type) %>%
  summarise(exact_rate = mean(exact_smaller_contained), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = pair_type, values_from = exact_rate) %>%
  drop_na()

pairtype_exact_tests <- bind_rows(
  paired_wilcox_from_wide(block_exact, "E-E", "E-I") %>% dplyr::mutate(comparison = "E-E vs E-I"),
  paired_wilcox_from_wide(block_exact, "E-E", "I-I") %>% dplyr::mutate(comparison = "E-E vs I-I"),
  paired_wilcox_from_wide(block_exact, "E-I", "I-I") %>% dplyr::mutate(comparison = "E-I vs I-I")
)
pairtype_exact_tests <- pairtype_exact_tests %>%
  dplyr::mutate(fdr_bh = p.adjust(p_value, method = "BH"))
pairtype_exact_tests <- pairtype_exact_tests[, c("comparison", setdiff(names(pairtype_exact_tests), "comparison"))]

layer_pairtype_exact_tests <- bind_rows(lapply(sort(unique(pairs$layer)), function(lay) {
  piv <- pairs %>%
    filter(layer == lay) %>%
    group_by(slide, pair_type) %>%
    summarise(exact_rate = mean(exact_smaller_contained), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = pair_type, values_from = exact_rate)
  out <- list()
  for (comp in list(c("E-E", "E-I"), c("E-E", "I-I"), c("E-I", "I-I"))) {
    a <- comp[1]; b <- comp[2]
    if (all(c(a, b) %in% names(piv))) {
      # 修复处：显式调用 dplyr::select
      sub <- piv %>% dplyr::select(all_of(c(a, b))) %>% tidyr::drop_na()
      if (nrow(sub) >= 5) {
        wt <- suppressWarnings(wilcox.test(sub[[a]], sub[[b]], paired = TRUE, exact = FALSE))
        out[[paste(a, b)]] <- tibble(
          layer = lay,
          comparison = paste(a, "vs", b),
          n_slides = nrow(sub),
          mean_diff = mean(sub[[a]] - sub[[b]]),
          median_diff = median(sub[[a]] - sub[[b]]),
          rank_biserial = rank_biserial_paired(sub[[a]], sub[[b]]),
          p_value = wt$p.value
        )
      }
    }
  }
  bind_rows(out)
})) %>%
  mutate(fdr_bh = p.adjust(p_value, method = "BH"))

host_glu_subclass_enrichment <- fisher_enrichment(ei$glu_subclass, ei$gaba_in_glu_exact, min_n = 50)
guest_gaba_subclass_enrichment <- fisher_enrichment(ei$gaba_subclass, ei$gaba_in_glu_exact, min_n = 50)
subclass_pair_fisher_enrichment <- fisher_enrichment(ei$subpair_clean, ei$gaba_in_glu_exact, min_n = 25)

subclass_pair_block_background <- bind_rows(lapply(split(ei, ei$subpair_clean), function(sub) {
  sp <- unique(sub$subpair_clean)
  if (nrow(sub) < 20 || n_distinct(sub$block) < 10) return(NULL)
  sp_block <- sub %>%
    group_by(block) %>%
    summarise(sp_rate = mean(gaba_in_glu_exact), .groups = "drop")
  other_bg <- ei %>%
    filter(subpair_clean != sp) %>%
    group_by(block) %>%
    summarise(other_rate = mean(gaba_in_glu_exact), .groups = "drop")
  merged <- inner_join(sp_block, other_bg, by = "block")
  if (nrow(merged) < 10) return(NULL)
  wt <- suppressWarnings(wilcox.test(merged$sp_rate, merged$other_rate, paired = TRUE, exact = FALSE))
  tibble(
    subpair_clean = sp,
    n_edges = nrow(sub),
    n_blocks = nrow(merged),
    pair_exact_rate = mean(sub$gaba_in_glu_exact),
    mean_diff_vs_block_EI_bg = mean(merged$sp_rate - merged$other_rate),
    median_diff_vs_block_EI_bg = median(merged$sp_rate - merged$other_rate),
    rank_biserial = rank_biserial_paired(merged$sp_rate, merged$other_rate),
    p_value = wt$p.value
  )
})) %>%
  mutate(
    fdr_bh = p.adjust(p_value, method = "BH"),
    direction = ifelse(mean_diff_vs_block_EI_bg > 0, "higher_containment_than_bg", "lower_containment_than_bg")
  ) %>%
  arrange(fdr_bh)

ee <- pairs %>% filter(pair_type == "E-E")
ii <- pairs %>% filter(pair_type == "I-I")

ee_long <- bind_rows(
  ee %>% transmute(label = label_a, subclass = subclass_a_clean, layer = layer, size = size_a, partner = label_b, jaccard = jaccard),
  ee %>% transmute(label = label_b, subclass = subclass_b_clean, layer = layer, size = size_b, partner = label_a, jaccard = jaccard)
)

ii_long <- bind_rows(
  ii %>% transmute(label = label_a, subclass = subclass_a_clean, layer = layer, size = size_a, partner = label_b, jaccard = jaccard),
  ii %>% transmute(label = label_b, subclass = subclass_b_clean, layer = layer, size = size_b, partner = label_a, jaccard = jaccard)
)

glu_cluster_stats <- ei %>%
  group_by(glu_label, glu_subclass, layer) %>%
  summarise(
    glu_size = dplyr::first(glu_size),
    n_gaba_partners = n_distinct(gaba_label),
    n_ei_edges = n(),
    exact_host_count = sum(gaba_in_glu_exact),
    exact_host_rate = mean(gaba_in_glu_exact),
    near95_host_count = sum(gaba_in_glu_95),
    mean_ei_jaccard = mean(jaccard),
    median_ei_jaccard = median(jaccard),
    mean_size_ratio = mean(size_ratio_glu_to_gaba),
    mean_overlap = mean(overlap),
    .groups = "drop"
  ) %>%
  mutate(ever_exact_host = exact_host_count > 0) %>%
  left_join(ee_long %>% group_by(label) %>% summarise(n_ee_partners = n_distinct(partner), mean_ee_jaccard = mean(jaccard), .groups = "drop"),
            by = c("glu_label" = "label")) %>%
  left_join(clusters %>% dplyr::select(label, slide) %>% distinct(), by = c("glu_label" = "label")) %>%
  mutate(log_glu_size = log(glu_size))

gaba_cluster_stats <- ei %>%
  group_by(gaba_label, gaba_subclass, layer) %>%
  summarise(
    gaba_size = dplyr::first(gaba_size),
    n_glu_partners = n_distinct(glu_label),
    n_ei_edges = n(),
    exact_contained_count = sum(gaba_in_glu_exact),
    exact_contained_rate = mean(gaba_in_glu_exact),
    near95_contained_count = sum(gaba_in_glu_95),
    mean_ei_jaccard = mean(jaccard),
    median_ei_jaccard = median(jaccard),
    mean_size_ratio = mean(size_ratio_glu_to_gaba),
    mean_overlap = mean(overlap),
    .groups = "drop"
  ) %>%
  mutate(ever_exact_contained = exact_contained_count > 0) %>%
  left_join(ii_long %>% group_by(label) %>% summarise(n_ii_partners = n_distinct(partner), mean_ii_jaccard = mean(jaccard), .groups = "drop"),
            by = c("gaba_label" = "label")) %>%
  left_join(clusters %>% dplyr::select(label, slide) %>% distinct(), by = c("gaba_label" = "label")) %>%
  mutate(log_gaba_size = log(gaba_size))

glu_cluster_host_tests <- bind_rows(lapply(
  c("glu_size", "n_gaba_partners", "n_ei_edges", "exact_host_rate", "mean_ei_jaccard",
    "mean_size_ratio", "mean_overlap", "n_ee_partners", "mean_ee_jaccard"),
  function(m) mw_test(glu_cluster_stats, "ever_exact_host", m)
)) %>% mutate(fdr_bh = p.adjust(p_value, method = "BH"))

gaba_cluster_guest_tests <- bind_rows(lapply(
  c("gaba_size", "n_glu_partners", "n_ei_edges", "exact_contained_rate", "mean_ei_jaccard",
    "mean_size_ratio", "mean_overlap", "n_ii_partners", "mean_ii_jaccard"),
  function(m) mw_test(gaba_cluster_stats, "ever_exact_contained", m)
)) %>% mutate(fdr_bh = p.adjust(p_value, method = "BH"))

fit_glu_host <- glm(as.integer(ever_exact_host) ~ log_glu_size + mean_ee_jaccard + n_ee_partners + layer,
                    data = glu_cluster_stats %>% drop_na(mean_ee_jaccard, n_ee_partners),
                    family = binomial())
fit_gaba_guest <- glm(as.integer(ever_exact_contained) ~ log_gaba_size + mean_ii_jaccard + n_ii_partners + layer,
                      data = gaba_cluster_stats %>% drop_na(mean_ii_jaccard, n_ii_partners),
                      family = binomial())

glu_cluster_logit <- robust_coef_table(fit_glu_host,
                                       cluster = glu_cluster_stats %>% drop_na(mean_ee_jaccard, n_ee_partners) %>% pull(slide),
                                       exponentiate = TRUE)
gaba_cluster_logit <- robust_coef_table(fit_gaba_guest,
                                        cluster = gaba_cluster_stats %>% drop_na(mean_ii_jaccard, n_ii_partners) %>% pull(slide),
                                        exponentiate = TRUE)

cluster_counts_block <- clusters %>%
  dplyr::count(block, cell_Neruon_type) %>%
  pivot_wider(names_from = cell_Neruon_type, values_from = n, values_fill = 0) %>%
  rename(nE = Glut, nI = Gaba)

edge_counts_block <- pairs %>%
  dplyr::count(block, pair_type) %>%
  pivot_wider(names_from = pair_type, values_from = n, values_fill = 0)

dens <- cluster_counts_block %>%
  left_join(edge_counts_block, by = "block") %>%
  mutate(
    possible_EE = nE * (nE - 1) / 2,
    possible_II = nI * (nI - 1) / 2,
    possible_EI = nE * nI,
    density_EE = `E-E` / possible_EE,
    density_EI = `E-I` / possible_EI,
    density_II = `I-I` / possible_II
  )

block_pair_summary <- pairs %>%
  group_by(block, pair_type) %>%
  summarise(
    mean_jaccard = mean(jaccard),
    exact_rate = mean(exact_smaller_contained),
    near95_rate = mean(near95_smaller_contained),
    .groups = "drop"
  )

block_ei <- ei %>%
  group_by(block) %>%
  summarise(
    ei_n = n(),
    gaba_in_glu_exact_rate = mean(gaba_in_glu_exact),
    gaba_in_glu_95_rate = mean(gaba_in_glu_95),
    mean_jaccard_ei = mean(jaccard),
    mean_size_ratio_ei = mean(size_ratio_glu_to_gaba),
    .groups = "drop"
  )

block_all <- block_ei %>%
  left_join(block_pair_summary %>% dplyr::select(block, pair_type, mean_jaccard) %>% pivot_wider(names_from = pair_type, values_from = mean_jaccard),
            by = "block") %>%
  left_join(block_pair_summary %>% dplyr::select(block, pair_type, exact_rate) %>% pivot_wider(names_from = pair_type, values_from = exact_rate, names_prefix = "exact_rate_"),
            by = "block") %>%
  left_join(block_pair_summary %>% dplyr::select(block, pair_type, near95_rate) %>% pivot_wider(names_from = pair_type, values_from = near95_rate, names_prefix = "near95_rate_"),
            by = "block") %>%
  left_join(dens %>% dplyr::select(block, density_EE, density_EI, density_II), by = "block")

corr_targets <- c("E-E", "E-I", "I-I", "exact_rate_E-E", "exact_rate_E-I", "exact_rate_E-I",
                  "density_EE", "density_EI", "density_II", "mean_size_ratio_ei")

block_correlations <- bind_rows(lapply(corr_targets, function(target) {
  # 修复处：显式调用 dplyr::select
  sub <- block_all %>% dplyr::select(gaba_in_glu_exact_rate, all_of(target)) %>% drop_na()
  if (nrow(sub) < 10) return(NULL)
  ct <- cor.test(sub$gaba_in_glu_exact_rate, sub[[target]], method = "spearman", exact = FALSE)
  tibble(
    target = target,
    n_blocks = nrow(sub),
    spearman_rho = unname(ct$estimate),
    p_value = ct$p.value
  )
})) %>%
  mutate(fdr_bh = p.adjust(p_value, method = "BH")) %>%
  arrange(fdr_bh)

gaba_exact_hosts <- ei %>%
  filter(gaba_in_glu_exact) %>%
  group_by(gaba_label, gaba_subclass, layer) %>%
  summarise(
    n_exact_glu_hosts = n_distinct(glu_label),
    n_host_subclasses = n_distinct(glu_subclass),
    mean_host_size = mean(glu_size),
    guest_size = dplyr::first(gaba_size),
    mean_jaccard = mean(jaccard),
    .groups = "drop"
  ) %>%
  mutate(multi_host = n_exact_glu_hosts >= 2)

glu_exact_hosts_rev <- ei %>%
  filter(glu_in_gaba_exact) %>%
  group_by(glu_label, glu_subclass, layer) %>%
  summarise(
    n_exact_gaba_hosts = n_distinct(gaba_label),
    host_subclasses = n_distinct(gaba_subclass),
    guest_size = dplyr::first(glu_size),
    mean_host_size = mean(gaba_size),
    mean_jaccard = mean(jaccard),
    .groups = "drop"
  ) %>%
  mutate(multi_host = n_exact_gaba_hosts >= 2)

multi_host_ft <- fisher.test(matrix(
  c(sum(gaba_exact_hosts$multi_host), sum(!gaba_exact_hosts$multi_host),
    sum(glu_exact_hosts_rev$multi_host), sum(!glu_exact_hosts_rev$multi_host)),
  nrow = 2, byrow = TRUE
))

multi_host_direction_comparison <- tibble(
  direction = c("GABA guest in Glu", "Glu guest in GABA"),
  n_guest_clusters = c(nrow(gaba_exact_hosts), nrow(glu_exact_hosts_rev)),
  multi_host_n = c(sum(gaba_exact_hosts$multi_host), sum(glu_exact_hosts_rev$multi_host)),
  multi_host_rate = c(mean(gaba_exact_hosts$multi_host), mean(glu_exact_hosts_rev$multi_host))
)

top_exact_gaba_in_glu_pairs <- ei %>%
  filter(gaba_in_glu_exact) %>%
  mutate(pair_name_clean = paste(glu_subclass, gaba_subclass, sep = " × ")) %>%
  arrange(desc(jaccard), desc(overlap), desc(glu_size)) %>%
  # 修复处：显式调用 dplyr::select
  dplyr::select(layer, pair_name_clean, glu_size, gaba_size, overlap, jaccard, size_ratio_glu_to_gaba, glu_label, gaba_label) %>%
  slice_head(n = 20)

top_glu_host_clusters <- ei %>%
  filter(gaba_in_glu_exact) %>%
  group_by(glu_label, glu_subclass, layer) %>%
  summarise(
    host_size = dplyr::first(glu_size),
    n_exact_gaba_guests = n_distinct(gaba_label),
    n_guest_subclasses = n_distinct(gaba_subclass),
    mean_guest_size = mean(gaba_size),
    mean_jaccard = mean(jaccard),
    .groups = "drop"
  ) %>%
  arrange(desc(n_exact_gaba_guests), desc(host_size)) %>%
  slice_head(n = 20)

top_gaba_guest_clusters <- gaba_exact_hosts %>%
  arrange(desc(n_exact_glu_hosts), desc(guest_size)) %>%
  slice_head(n = 20)

dec_dat <- tmp %>%
  mutate(size_ratio_decile = ntile(abs_size_ratio, 10)) %>%
  group_by(size_ratio_decile, smaller_type) %>%
  summarise(
    n = n(),
    exact_rate = mean(smaller_contained_exact),
    median_abs_size_ratio = median(abs_size_ratio),
    .groups = "drop"
  )

size_ratio_decile_rates <- dec_dat

# ---------------------------
# Save tables
# ---------------------------
write_csv(overview, file.path(tbl_dir, "table_01_overview.csv"))
write_csv(pairtype_containment, file.path(tbl_dir, "table_02_pairtype_containment.csv"))
write_csv(direction_tests, file.path(tbl_dir, "table_03_direction_tests.csv"))
write_csv(smaller_exact_logit, file.path(tbl_dir, "table_04_smaller_type_logit_exact.csv"))
write_csv(smaller_95_logit, file.path(tbl_dir, "table_05_smaller_type_logit_95.csv"))
write_csv(exact_vs_other_desc, file.path(tbl_dir, "table_06_exact_vs_other_EI_descriptives.csv"))
write_csv(block_metric_tests, file.path(tbl_dir, "table_07_block_paired_exact_vs_other_metrics.csv"))
write_csv(jaccard_regression, file.path(tbl_dir, "table_08_jaccard_regression_exact_vs_other.csv"))
write_csv(layer_rate, file.path(tbl_dir, "table_09_layer_rates.csv"))
write_csv(layer_tests, file.path(tbl_dir, "table_10_layer_enrichment_tests.csv"))
write_csv(layer_exact_logit, file.path(tbl_dir, "table_11_layer_logit_exact.csv"))
write_csv(pairtype_exact_tests, file.path(tbl_dir, "table_12_pairtype_exact_rate_tests.csv"))
write_csv(layer_pairtype_exact_tests, file.path(tbl_dir, "table_13_layerwise_pairtype_exact_rate_tests.csv"))
write_csv(host_glu_subclass_enrichment, file.path(tbl_dir, "table_14_host_glu_subclass_enrichment.csv"))
write_csv(guest_gaba_subclass_enrichment, file.path(tbl_dir, "table_15_guest_gaba_subclass_enrichment.csv"))
write_csv(subclass_pair_fisher_enrichment, file.path(tbl_dir, "table_16_subclass_pair_fisher_enrichment.csv"))
write_csv(subclass_pair_block_background, file.path(tbl_dir, "table_17_subclass_pair_block_background_tests.csv"))
write_csv(glu_cluster_host_tests, file.path(tbl_dir, "table_18_glu_cluster_host_tests.csv"))
write_csv(gaba_cluster_guest_tests, file.path(tbl_dir, "table_19_gaba_cluster_guest_tests.csv"))
write_csv(glu_cluster_logit, file.path(tbl_dir, "table_20_glu_cluster_logit_ever_host.csv"))
write_csv(gaba_cluster_logit, file.path(tbl_dir, "table_21_gaba_cluster_logit_ever_contained.csv"))
write_csv(block_correlations, file.path(tbl_dir, "table_22_block_correlations.csv"))
write_csv(multi_host_direction_comparison, file.path(tbl_dir, "table_23_multi_host_direction_comparison.csv"))
write_csv(top_exact_gaba_in_glu_pairs, file.path(tbl_dir, "table_24_top_exact_gaba_in_glu_pairs.csv"))
write_csv(top_glu_host_clusters, file.path(tbl_dir, "table_25_top_glu_host_clusters.csv"))
write_csv(top_gaba_guest_clusters, file.path(tbl_dir, "table_26_top_gaba_guest_clusters.csv"))
write_csv(size_ratio_decile_rates, file.path(tbl_dir, "table_27_size_ratio_decile_rates.csv"))

# ---------------------------
# Top-journal figures
# ---------------------------

# ---------- plotting helpers ----------
wrap_plot_title <- function(x, width = 52) stringr::str_wrap(x, width = width)
wrap_plot_subtitle <- function(x, width = 92) stringr::str_wrap(x, width = width)

plot_theme_safe <- ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.margin = ggplot2::margin(20, 28, 14, 18),
    plot.title = ggplot2::element_text(
      face = "bold", size = 16.5, colour = "black",
      lineheight = 1.06, margin = ggplot2::margin(b = 8)
    ),
    plot.subtitle = ggplot2::element_text(
      size = 11.8, colour = "#4D4D4D",
      lineheight = 1.10, margin = ggplot2::margin(b = 12)
    ),
    axis.title = ggplot2::element_text(size = 13.5, colour = "black"),
    axis.text = ggplot2::element_text(size = 11.5, colour = "black"),
    legend.position = "top",
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 11.2)
  )

# Figure 1: Pair-type exact containment rate with paired FDR annotations
block_exact_plot <- pairs %>%
  group_by(block, pair_type) %>%
  summarise(exact_rate = mean(exact_smaller_contained), .groups = "drop") %>%
  mutate(pair_type = factor(pair_type, levels = c("E-E", "E-I", "I-I")))

sig1 <- tibble()
if (nrow(pairtype_exact_tests) > 0) {
  ee_ei <- pairtype_exact_tests %>% filter(comparison == "E-E vs E-I") %>% pull(fdr_bh)
  ei_ii <- pairtype_exact_tests %>% filter(comparison == "E-I vs I-I") %>% pull(fdr_bh)
  ee_ii <- pairtype_exact_tests %>% filter(comparison == "E-E vs I-I") %>% pull(fdr_bh)
  y1 <- max(block_exact_plot$exact_rate, na.rm = TRUE) * 1.12
  y2 <- max(block_exact_plot$exact_rate, na.rm = TRUE) * 1.24
  y3 <- max(block_exact_plot$exact_rate, na.rm = TRUE) * 1.36
  sig1 <- add_sig_bracket(sig1, 1, 2, y1, paste0("FDR=", fmt_p_short(ee_ei)))
  sig1 <- add_sig_bracket(sig1, 2, 3, y2, paste0("FDR=", fmt_p_short(ei_ii)))
  sig1 <- add_sig_bracket(sig1, 1, 3, y3, paste0("FDR=", fmt_p_short(ee_ii)))
}

fig1 <- ggplot(block_exact_plot, aes(pair_type, exact_rate, fill = pair_type)) +
  geom_violin(width = 0.90, alpha = 0.16, color = NA, trim = FALSE) +
  geom_boxplot(width = 0.34, outlier.shape = NA, linewidth = 0.7, alpha = 0.90) +
  geom_jitter(width = 0.12, size = 1.3, alpha = 0.22, color = "#4D4D4D") +
  geom_segment(
    data = sig1,
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE, linewidth = 0.55
  ) +
  geom_text(
    data = sig1 %>% filter(!is.na(label)),
    aes(x = xlab, y = ylab * 1.03, label = label),
    inherit.aes = FALSE, size = 3.8
  ) +
  scale_fill_manual(values = col_pair) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.22))
  ) +
  labs(
    title = wrap_plot_title("E-I has the highest exact-containment rate across the three interaction classes"),
    subtitle = wrap_plot_subtitle(
      paste0(
        "Block-level smaller-cluster exact containment rate. Means: E-I ",
        percent(mean(pairs$exact_smaller_contained[pairs$pair_type == "E-I"]), accuracy = 0.01),
        ", E-E ",
        percent(mean(pairs$exact_smaller_contained[pairs$pair_type == "E-E"]), accuracy = 0.01),
        ", I-I ",
        percent(mean(pairs$exact_smaller_contained[pairs$pair_type == "I-I"]), accuracy = 0.01)
      )
    ),
    x = NULL,
    y = "Exact containment rate of the smaller cluster"
  ) +
  coord_cartesian(clip = "off") +
  plot_theme_safe
save_pub(fig1, "fig01_pairtype_exact_containment_topjournal", 9.6, 6.2)

# Figure 2: Directional exact asymmetry counts
count_df <- tibble(
  direction = factor(c("GABA in Glu", "Glu in GABA"), levels = c("GABA in Glu", "Glu in GABA")),
  n_edges = c(sum(ei$gaba_in_glu_exact), sum(ei$glu_in_gaba_exact)),
  prop = c(sum(ei$gaba_in_glu_exact), sum(ei$glu_in_gaba_exact)) / n_exact_dir
)

fig2 <- ggplot(count_df, aes(direction, n_edges, fill = direction)) +
  geom_col(width = 0.62, color = NA) +
  geom_text(
    aes(label = paste0(comma(n_edges), "\n(", percent(prop, accuracy = 0.1), ")")),
    vjust = -0.28, size = 4.6, lineheight = 0.95
  ) +
  scale_fill_manual(values = col_dir) +
  scale_y_continuous(labels = comma_format(), expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = wrap_plot_title("Exact directional containment is strongly biased toward GABA-in-Glu"),
    subtitle = wrap_plot_subtitle(
      paste0(
        percent(k_gaba_dir / n_exact_dir, accuracy = 0.1),
        " of all exact directional E-I containment edges are GABA in Glu; binomial p = ",
        fmt_p_short(binom_exact$p.value)
      )
    ),
    x = NULL,
    y = "Number of E-I exact-direction edges"
  ) +
  plot_theme_safe +
  theme(legend.position = "none")
save_pub(fig2, "fig02_directional_exact_asymmetry_topjournal", 8.4, 5.8)

# Figure 3: Smaller-cluster comparison
small_df <- tibble(
  smaller_type = factor(
    c("Smaller cluster = GABA", "Smaller cluster = Glu"),
    levels = c("Smaller cluster = GABA", "Smaller cluster = Glu")
  ),
  rate = c(
    mean(ei$gaba_in_glu_exact[ei$gaba_is_smaller]),
    mean(ei$glu_in_gaba_exact[ei$glu_is_smaller])
  ),
  n = c(sum(ei$gaba_is_smaller), sum(ei$glu_is_smaller))
)

fig3 <- ggplot(small_df, aes(smaller_type, rate, fill = smaller_type)) +
  geom_col(width = 0.62) +
  geom_text(
    aes(label = paste0(percent(rate, accuracy = 0.1), "\n(n=", comma(n), ")")),
    vjust = -0.26, size = 4.5, lineheight = 0.95
  ) +
  scale_fill_manual(values = col_small) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.16))) +
  labs(
    title = wrap_plot_title("The directional bias persists even when only the smaller cluster is compared"),
    subtitle = wrap_plot_subtitle(
      paste0(
        "Fisher OR = ",
        sprintf("%.2f", unname(fisher_smaller_exact$estimate)),
        "; p = ",
        fmt_p_short(fisher_smaller_exact$p.value),
        "; adjusted OR = ",
        sprintf("%.2f", smaller_exact_logit$OR[smaller_exact_logit$term == "smaller_is_gaba"]),
        ", p = ",
        fmt_p_short(smaller_exact_logit$p.value[smaller_exact_logit$term == "smaller_is_gaba"])
      )
    ),
    x = NULL,
    y = "Exact containment rate of the smaller cluster"
  ) +
  plot_theme_safe +
  theme(legend.position = "none")
save_pub(fig3, "fig03_smaller_cluster_direction_bias_topjournal", 8.8, 5.8)

# Figure 4: Adjusted OR forest for smaller-cluster exact model
forest_small <- smaller_exact_logit %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term_label = recode(
      term,
      `smaller_is_gaba` = "Smaller cluster = GABA",
      `log_abs_size_ratio` = "log(larger/smaller)",
      `layerlayer 2/3` = "Layer 2/3",
      `layerlayer 4` = "Layer 4",
      `layerlayer 5` = "Layer 5",
      `layerlayer 6a` = "Layer 6a",
      `layerlayer 6b` = "Layer 6b"
    )
  ) %>%
  mutate(term_label = factor(term_label, levels = rev(term_label)))

fig4 <- ggplot(forest_small, aes(x = OR, y = term_label)) +
  geom_vline(xintercept = 1, linetype = 2, color = "#7F7F7F", linewidth = 0.55) +
  geom_errorbarh(aes(xmin = OR_low, xmax = OR_high), height = 0.18, linewidth = 0.7, color = "#4D4D4D") +
  geom_point(size = 2.8, color = "#C44E52") +
  scale_x_log10(labels = label_number(accuracy = 0.01)) +
  labs(
    title = wrap_plot_title("The smaller-GABA effect remains after controlling for size ratio and layer"),
    subtitle = wrap_plot_subtitle(
      paste0(
        "Adjusted OR for smaller cluster = GABA: ",
        sprintf("%.2f", smaller_exact_logit$OR[smaller_exact_logit$term == "smaller_is_gaba"]),
        "; p = ",
        fmt_p_short(smaller_exact_logit$p.value[smaller_exact_logit$term == "smaller_is_gaba"])
      )
    ),
    x = "Adjusted odds ratio (log scale)",
    y = NULL
  ) +
  plot_theme_safe
save_pub(fig4, "fig04_smaller_cluster_logit_forest_topjournal", 9.0, 6.0)

# Figure 5: Jaccard exact vs other E-I, with paired block summary and regression effect
exact_jaccard_plot <- ei %>%
  mutate(group = ifelse(gaba_in_glu_exact, "Exact GABA-in-Glu", "Other E-I")) %>%
  group_by(block, group) %>%
  summarise(value = mean(jaccard), .groups = "drop") %>%
  drop_na() %>%
  mutate(group = factor(group, levels = c("Other E-I", "Exact GABA-in-Glu")))

fig5 <- ggplot(exact_jaccard_plot, aes(group, value, fill = group)) +
  geom_violin(width = 0.86, alpha = 0.16, color = NA, trim = FALSE) +
  geom_boxplot(width = 0.30, outlier.shape = NA, linewidth = 0.7, alpha = 0.9) +
  geom_jitter(width = 0.10, size = 1.2, alpha = 0.20, color = "#4D4D4D") +
  scale_fill_manual(values = col_binary) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    title = wrap_plot_title("Exact GABA-in-Glu edges have higher Jaccard despite larger size ratios"),
    subtitle = wrap_plot_subtitle(
      paste0(
        "Block-level mean Jaccard; adjusted logit(Jaccard) coefficient for exact containment = ",
        sprintf("%.4f", jaccard_regression$estimate[jaccard_regression$term == "as.integer(gaba_in_glu_exact)"]),
        "; p = ",
        fmt_p_short(jaccard_regression$p.value[jaccard_regression$term == "as.integer(gaba_in_glu_exact)"])
      )
    ),
    x = NULL,
    y = "Mean Jaccard per block"
  ) +
  plot_theme_safe +
  theme(legend.position = "none")
save_pub(fig5, "fig05_exact_vs_other_jaccard_topjournal", 9.0, 6.0)

# Figure 6: Layer-wise exact rate
layer_plot <- layer_rate %>%
  mutate(layer = factor(layer, levels = c("layer 1", "layer 2/3", "layer 4", "layer 5", "layer 6a", "layer 6b")))

fig6 <- ggplot(layer_plot, aes(layer, exact_rate, fill = layer)) +
  geom_col(width = 0.70) +
  geom_text(
    aes(label = paste0(percent(exact_rate, accuracy = 0.1), "\n(n=", comma(n_edges), ")")),
    vjust = -0.22, size = 4.1, lineheight = 0.95
  ) +
  scale_fill_manual(values = col_layer) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.16))) +
  labs(
    title = wrap_plot_title("Layer dependence of exact GABA-in-Glu containment"),
    subtitle = wrap_plot_subtitle("Rates among observed E-I pairs; shallow and middle layers show the highest levels"),
    x = NULL,
    y = "Exact GABA-in-Glu rate"
  ) +
  plot_theme_safe +
  theme(legend.position = "none")
save_pub(fig6, "fig06_layer_exact_rate_topjournal", 8.8, 5.8)

# Figure 7: Pair-specific background deviation
sel_pos <- subclass_pair_block_background %>%
  filter(fdr_bh < 0.05, mean_diff_vs_block_EI_bg > 0) %>%
  arrange(desc(mean_diff_vs_block_EI_bg)) %>%
  slice_head(n = 6)

sel_neg <- subclass_pair_block_background %>%
  filter(fdr_bh < 0.05, mean_diff_vs_block_EI_bg < 0) %>%
  arrange(mean_diff_vs_block_EI_bg) %>%
  slice_head(n = 10)

sel_pair_bg <- bind_rows(sel_pos, sel_neg) %>%
  arrange(mean_diff_vs_block_EI_bg) %>%
  mutate(
    subpair_clean = factor(subpair_clean, levels = subpair_clean),
    group = ifelse(
      mean_diff_vs_block_EI_bg > 0,
      "Higher than block E-I background",
      "Lower than block E-I background"
    )
  )

fig7 <- ggplot(sel_pair_bg, aes(mean_diff_vs_block_EI_bg, subpair_clean, fill = group)) +
  geom_col() +
  geom_vline(xintercept = 0, linewidth = 0.55, color = "#7F7F7F") +
  scale_fill_manual(values = c(
    "Higher than block E-I background" = "#C44E52",
    "Lower than block E-I background" = "#4E79A7"
  )) +
  labs(
    title = wrap_plot_title("Containment propensity is strongly pair-specific rather than uniform across inhibitory subclasses"),
    subtitle = wrap_plot_subtitle("Δ exact-containment rate relative to the within-block E-I background"),
    x = "Δ exact-containment rate vs block E-I background",
    y = NULL
  ) +
  plot_theme_safe
save_pub(fig7, "fig07_pair_specific_background_deviation_topjournal", 10.8, 6.4)

# Figure 8: GABA guest-cluster logistic forest
forest_gaba <- gaba_cluster_logit %>%
  filter(term %in% c("log_gaba_size", "mean_ii_jaccard", "n_ii_partners")) %>%
  mutate(
    term_label = recode(
      term,
      `log_gaba_size` = "log(GABA cluster size)",
      `mean_ii_jaccard` = "Mean I-I Jaccard",
      `n_ii_partners` = "Number of I-I partners"
    )
  ) %>%
  mutate(term_label = factor(term_label, levels = rev(term_label)))

fig8 <- ggplot(forest_gaba, aes(x = OR, y = term_label)) +
  geom_vline(xintercept = 1, linetype = 2, color = "#7F7F7F", linewidth = 0.55) +
  geom_errorbarh(aes(xmin = OR_low, xmax = OR_high), height = 0.18, linewidth = 0.7, color = "#4D4D4D") +
  geom_point(size = 2.8, color = "#4E79A7") +
  scale_x_log10(labels = label_number(accuracy = 0.01)) +
  labs(
    title = wrap_plot_title("Exact-contained GABA clusters are smaller and less embedded within inhibitory modules"),
    subtitle = wrap_plot_subtitle(
      paste0(
        "OR(log GABA size) = ",
        sprintf("%.3f", gaba_cluster_logit$OR[gaba_cluster_logit$term == "log_gaba_size"]),
        "; OR(mean I-I Jaccard) = ",
        sprintf("%.3f", gaba_cluster_logit$OR[gaba_cluster_logit$term == "mean_ii_jaccard"])
      )
    ),
    x = "Adjusted odds ratio (log scale)",
    y = NULL
  ) +
  plot_theme_safe
save_pub(fig8, "fig08_gaba_guest_cluster_forest_topjournal", 8.8, 5.6)

# Figure 9: Multi-host direction asymmetry
fig9 <- ggplot(multi_host_direction_comparison, aes(direction, multi_host_rate, fill = direction)) +
  geom_col(width = 0.62) +
  geom_text(
    aes(label = paste0(percent(multi_host_rate, accuracy = 0.1), "\n(n=", comma(n_guest_clusters), ")")),
    vjust = -0.24, size = 4.4, lineheight = 0.95
  ) +
  scale_fill_manual(values = col_multi) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.16))) +
  labs(
    title = wrap_plot_title("Multi-host exact coverage is itself directionally asymmetric"),
    subtitle = wrap_plot_subtitle(
      paste0(
        "Fisher OR = ",
        sprintf("%.2f", unname(multi_host_ft$estimate)),
        "; p = ",
        fmt_p_short(multi_host_ft$p.value)
      )
    ),
    x = NULL,
    y = "Fraction of exact guest clusters with ≥2 exact hosts"
  ) +
  plot_theme_safe +
  theme(legend.position = "none")
save_pub(fig9, "fig09_multi_host_direction_topjournal", 8.2, 5.8)

# ---------------------------
# Figure manifest
# ---------------------------
figure_manifest <- tibble(
  filename = c(
    "fig01_pairtype_exact_containment_topjournal",
    "fig02_directional_exact_asymmetry_topjournal",
    "fig03_smaller_cluster_direction_bias_topjournal",
    "fig04_smaller_cluster_logit_forest_topjournal",
    "fig05_exact_vs_other_jaccard_topjournal",
    "fig06_layer_exact_rate_topjournal",
    "fig07_pair_specific_background_deviation_topjournal",
    "fig08_gaba_guest_cluster_forest_topjournal",
    "fig09_multi_host_direction_topjournal"
  ),
  summary = c(
    "E-I exact containment rate vs E-E and I-I at block level",
    "Exact directional asymmetry: GABA in Glu vs Glu in GABA",
    "Smaller-GABA vs smaller-Glu exact containment rates",
    "Adjusted OR forest for smaller-cluster containment model",
    "Block-level Jaccard comparison for exact GABA-in-Glu vs other E-I edges",
    "Layer-specific exact GABA-in-Glu rates",
    "Pair-specific deviation from block-level E-I background",
    "Cluster-level logistic effects for exact-contained GABA guests",
    "Direction asymmetry of multi-host exact coverage"
  )
)

write_csv(figure_manifest, file.path(out_dir, "figure_manifest.csv"))

message("All outputs saved to: ", normalizePath(out_dir))