# ============================================================
# 兴奋-抑制空间簇进一步统计分析（严格参考代码修订版）
# 主题：聚焦 GABA cluster 完全被 GLU cluster 包含（I_in_E）
# 输入：E:/zaw/2603/neuron-neuron-partner.csv
# 输出：E:/zaw/2603/R_output_03
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
})

options(scipen = 999)

# -------------------------------
# 固定路径
# -------------------------------
base_dir   <- "E:/zaw/2603"
input_file <- file.path(base_dir, "neuron-neuron-partner.csv")
out_dir    <- file.path(base_dir, "R_output_03")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------
# 工具函数
# -------------------------------
fmt_p <- function(p) {
  ifelse(
    is.na(p), NA_character_,
    ifelse(
      p < 1e-300, "<1e-300",
      ifelse(p < 1e-4, format(p, scientific = TRUE, digits = 3),
             sprintf("%.6f", p))
    )
  )
}

fmt_p_plot <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-3) return(paste0("FDR=", sub("e", "×10^", format(p, scientific = TRUE, digits = 2))))
  paste0("FDR=", sprintf("%.3f", p))
}

type_to_ei <- function(x) ifelse(x == "Glut", "E", "I")

pair_type_from_types <- function(x, y) {
  purrr::map2_chr(x, y, ~ paste(sort(c(type_to_ei(.x), type_to_ei(.y))), collapse = "-"))
}

make_pair_key <- function(x, y) paste(pmin(x, y), pmax(x, y), sep = "||")

count_region_proxy <- function(x) stringr::str_count(x, "layer")

build_full_pairs_one_stratum <- function(df_sub) {
  if (nrow(df_sub) < 2) return(tibble())
  cmb <- utils::combn(seq_len(nrow(df_sub)), 2)
  tibble(
    stratum = df_sub$stratum[1],
    slide   = df_sub$slide[1],
    layer   = df_sub$layer[1],
    label_a = df_sub$label[cmb[1, ]],
    label_b = df_sub$label[cmb[2, ]],
    type_a  = df_sub$type[cmb[1, ]],
    type_b  = df_sub$type[cmb[2, ]],
    subclass_a = df_sub$subclass[cmb[1, ]],
    subclass_b = df_sub$subclass[cmb[2, ]],
    n_a = df_sub$n[cmb[1, ]],
    n_b = df_sub$n[cmb[2, ]]
  ) %>%
    mutate(
      pair_key = make_pair_key(label_a, label_b),
      pair_type = pair_type_from_types(type_a, type_b)
    )
}

cmh_binary_compare <- function(df, exposure_col, exposure_yes, outcome_col,
                               outcome_yes = TRUE, strata_col = "stratum",
                               comparison_label = NA_character_) {
  sp <- split(df, df[[strata_col]])
  arr_list <- list()
  strata_keep <- c()
  
  for (s in names(sp)) {
    d <- sp[[s]]
    a <- sum(d[[exposure_col]] == exposure_yes & d[[outcome_col]] == outcome_yes, na.rm = TRUE)
    b <- sum(d[[exposure_col]] == exposure_yes & d[[outcome_col]] != outcome_yes, na.rm = TRUE)
    c <- sum(d[[exposure_col]] != exposure_yes & d[[outcome_col]] == outcome_yes, na.rm = TRUE)
    d0 <- sum(d[[exposure_col]] != exposure_yes & d[[outcome_col]] != outcome_yes, na.rm = TRUE)
    
    if ((a + b) > 0 && (c + d0) > 0 && (a + c) > 0 && (b + d0) > 0) {
      arr_list[[length(arr_list) + 1]] <- matrix(c(a, b, c, d0), nrow = 2, byrow = TRUE)
      strata_keep <- c(strata_keep, s)
    }
  }
  
  if (length(arr_list) == 0) {
    return(tibble(
      comparison = comparison_label,
      n_strata = 0L,
      odds_ratio_MH = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p = NA_real_
    ))
  }
  
  arr <- array(
    0,
    dim = c(2, 2, length(arr_list)),
    dimnames = list(
      exposure = c(as.character(exposure_yes), "other"),
      outcome = c("yes", "no"),
      stratum = strata_keep
    )
  )
  
  for (i in seq_along(arr_list)) arr[, , i] <- arr_list[[i]]
  
  mh <- mantelhaen.test(arr)
  
  tibble(
    comparison = comparison_label,
    n_strata = length(arr_list),
    odds_ratio_MH = unname(mh$estimate),
    ci_low = unname(mh$conf.int[1]),
    ci_high = unname(mh$conf.int[2]),
    p = mh$p.value
  )
}

paired_stratum_wilcox <- function(df, group_col, metric, group1, group0,
                                  strata_cols = c("slide", "layer")) {
  tmp <- df %>%
    transmute(across(all_of(strata_cols)),
              group = .data[[group_col]],
              value = .data[[metric]]) %>%
    group_by(across(all_of(c(strata_cols, "group")))) %>%
    summarise(value = median(value, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = group, values_from = value)
  
  if (!(group1 %in% names(tmp)) || !(group0 %in% names(tmp))) {
    return(tibble(
      metric = metric,
      group1 = group1,
      group0 = group0,
      n_strata = 0L,
      median_diff = NA_real_,
      prop_group1_gt_group0 = NA_real_,
      p = NA_real_
    ))
  }
  
  sub <- tmp %>% select(all_of(strata_cols), all_of(group1), all_of(group0)) %>% drop_na()
  if (nrow(sub) == 0) {
    return(tibble(
      metric = metric,
      group1 = group1,
      group0 = group0,
      n_strata = 0L,
      median_diff = NA_real_,
      prop_group1_gt_group0 = NA_real_,
      p = NA_real_
    ))
  }
  
  wt <- wilcox.test(sub[[group1]], sub[[group0]], paired = TRUE, exact = FALSE)
  
  tibble(
    metric = metric,
    group1 = group1,
    group0 = group0,
    n_strata = nrow(sub),
    median_diff = median(sub[[group1]] - sub[[group0]], na.rm = TRUE),
    prop_group1_gt_group0 = mean(sub[[group1]] > sub[[group0]], na.rm = TRUE),
    p = wt$p.value
  )
}

group_compare_wilcox <- function(df, group_col, metric) {
  g1 <- df %>% filter(.data[[group_col]]) %>% pull(.data[[metric]])
  g0 <- df %>% filter(!.data[[group_col]]) %>% pull(.data[[metric]])
  wt <- wilcox.test(g1, g0, paired = FALSE, exact = FALSE)
  tibble(
    metric = metric,
    n_group1 = length(g1),
    n_group0 = length(g0),
    median_group1 = median(g1, na.rm = TRUE),
    median_group0 = median(g0, na.rm = TRUE),
    p = wt$p.value
  )
}

fisher_feature_enrichment <- function(df, feature_col, event_col = "status",
                                      event_label = "I_in_E",
                                      baseline_label = "partial",
                                      min_count = 10) {
  feat_counts <- table(df[[feature_col]])
  keep <- names(feat_counts)[feat_counts >= min_count]
  
  out <- purrr::map_dfr(keep, function(term) {
    a <- sum(df[[feature_col]] == term & df[[event_col]] == event_label, na.rm = TRUE)
    b <- sum(df[[feature_col]] == term & df[[event_col]] == baseline_label, na.rm = TRUE)
    c <- sum(df[[feature_col]] != term & df[[event_col]] == event_label, na.rm = TRUE)
    d <- sum(df[[feature_col]] != term & df[[event_col]] == baseline_label, na.rm = TRUE)
    ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))
    tibble(
      feature = term,
      count = a + b,
      event = a,
      baseline = b,
      odds_ratio = unname(ft$estimate),
      p = ft$p.value,
      event_rate = a / (a + b)
    )
  }) %>%
    mutate(fdr = p.adjust(p, method = "BH")) %>%
    arrange(fdr, desc(odds_ratio))
  
  out
}

cluster_subclass_enrichment <- function(df, subclass_col, group_col, min_count = 20) {
  sub_counts <- table(df[[subclass_col]])
  keep <- names(sub_counts)[sub_counts >= min_count]
  
  out <- purrr::map_dfr(keep, function(term) {
    a <- sum(df[[subclass_col]] == term & df[[group_col]], na.rm = TRUE)
    b <- sum(df[[subclass_col]] == term & !df[[group_col]], na.rm = TRUE)
    c <- sum(df[[subclass_col]] != term & df[[group_col]], na.rm = TRUE)
    d <- sum(df[[subclass_col]] != term & !df[[group_col]], na.rm = TRUE)
    ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))
    tibble(
      subclass = term,
      count = a + b,
      group_true = a,
      group_false = b,
      odds_ratio = unname(ft$estimate),
      p = ft$p.value,
      rate_true = a / (a + b)
    )
  }) %>%
    mutate(fdr = p.adjust(p, method = "BH")) %>%
    arrange(fdr, desc(odds_ratio))
  
  out
}

# -------------------------------
# 读取与重命名
# -------------------------------
raw <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  rename(
    label1 = `cluster.1_label`,
    slide1 = `cluster.1_slide`,
    layer1 = `cluster.1_layer`,
    region1 = `cluster.1_region`,
    type1 = `cluster.1_cell_Neruon_type`,
    subclass1 = `cluster.1_subclass`,
    n1 = `cluster.1_total_cell_num`,
    glu1 = `cluster.1_Glut_Neruon_cell_ids_num`,
    gaba1 = `cluster.1_GABA_Neruon_cell_ids_num`,
    p1 = `cluster.1_cauchy_combination_p`,
    eir1 = `cluster.1_E_I_Ratio`,
    label2 = `cluster.2_label`,
    slide2 = `cluster.2_slide`,
    layer2 = `cluster.2_layer`,
    region2 = `cluster.2_region`,
    type2 = `cluster.2_cell_Neruon_type`,
    subclass2 = `cluster.2_subclass`,
    n2 = `cluster.2_total_cell_num`,
    glu2 = `cluster.2_Glut_Neruon_cell_ids_num`,
    gaba2 = `cluster.2_GABA_Neruon_cell_ids_num`,
    p2 = `cluster.2_cauchy_combination_p`,
    eir2 = `cluster.2_E_I_Ratio`,
    overlap = `overlap_cell`,
    union = `union_cell`,
    jaccard = `jaccard`,
    ov1 = `cluster.1.overlap.percent`,
    ov2 = `cluster.2.overlap.percent`
  ) %>%
  mutate(pair_key = make_pair_key(label1, label2))

# -------------------------------
# 去重为唯一无向重叠配对
# -------------------------------
pairs <- purrr::imap_dfr(split(raw, raw$pair_key), function(g, key) {
  if (nrow(g) != 2) stop(paste("Pair key does not have exactly 2 rows:", key))
  labs <- sort(c(g$label1[1], g$label2[1]))
  a <- labs[1]
  b <- labs[2]
  r_ab <- g %>% filter(label1 == a, label2 == b)
  r_ba <- g %>% filter(label1 == b, label2 == a)
  if (nrow(r_ab) != 1 || nrow(r_ba) != 1) stop(paste("Cannot find both directions for:", key))
  
  tibble(
    pair_key = key,
    slide = r_ab$slide1,
    layer = r_ab$layer1,
    stratum = paste(r_ab$slide1, r_ab$layer1, sep = "__"),
    label_a = a,
    label_b = b,
    type_a = r_ab$type1,
    type_b = r_ab$type2,
    subclass_a = r_ab$subclass1,
    subclass_b = r_ab$subclass2,
    region_a = r_ab$region1,
    region_b = r_ab$region2,
    n_a = r_ab$n1,
    n_b = r_ab$n2,
    glu_a = r_ab$glu1,
    glu_b = r_ab$glu2,
    gaba_a = r_ab$gaba1,
    gaba_b = r_ab$gaba2,
    p_a = r_ab$p1,
    p_b = r_ab$p2,
    eir_a = r_ab$eir1,
    eir_b = r_ab$eir2,
    overlap = r_ab$overlap,
    union = r_ab$union,
    jaccard = r_ab$jaccard,
    ov_a = r_ab$ov1,
    ov_b = r_ab$ov2
  )
}) %>%
  mutate(
    pair_type = pair_type_from_types(type_a, type_b),
    contain_a_in_b = abs(ov_a - 1) < 1e-12,
    contain_b_in_a = abs(ov_b - 1) < 1e-12,
    any_containment = contain_a_in_b | contain_b_in_a
  )

# -------------------------------
# 唯一 cluster 表
# -------------------------------
clusters <- bind_rows(
  raw %>%
    transmute(label = label1, slide = slide1, layer = layer1, region = region1,
              type = type1, subclass = subclass1, n = n1, glu = glu1, gaba = gaba1,
              p = p1, eir = eir1),
  raw %>%
    transmute(label = label2, slide = slide2, layer = layer2, region = region2,
              type = type2, subclass = subclass2, n = n2, glu = glu2, gaba = gaba2,
              p = p2, eir = eir2)
) %>%
  distinct(label, .keep_all = TRUE) %>%
  mutate(
    stratum = paste(slide, layer, sep = "__"),
    glu_frac = glu / n,
    gaba_frac = gaba / n,
    region_count = count_region_proxy(region)
  )

# -------------------------------
# 重建所有可能配对宇宙
# -------------------------------
all_pairs <- purrr::map_dfr(split(clusters, clusters$stratum), build_full_pairs_one_stratum) %>%
  mutate(
    observed = pair_key %in% pairs$pair_key
  ) %>%
  left_join(pairs %>% select(pair_key, any_containment), by = "pair_key") %>%
  mutate(any_containment = tidyr::replace_na(any_containment, FALSE))

# -------------------------------
# 构建定向 E-I 表（固定 GLU / GABA 方向）
# -------------------------------
ei <- pairs %>%
  filter(pair_type == "E-I") %>%
  mutate(
    glu_label = if_else(type_a == "Glut", label_a, label_b),
    gaba_label = if_else(type_a == "Glut", label_b, label_a),
    glu_subclass = if_else(type_a == "Glut", subclass_a, subclass_b),
    gaba_subclass = if_else(type_a == "Glut", subclass_b, subclass_a),
    glu_region = if_else(type_a == "Glut", region_a, region_b),
    gaba_region = if_else(type_a == "Glut", region_b, region_a),
    glu_n = if_else(type_a == "Glut", n_a, n_b),
    gaba_n = if_else(type_a == "Glut", n_b, n_a),
    glu_glu = if_else(type_a == "Glut", glu_a, glu_b),
    glu_gaba = if_else(type_a == "Glut", gaba_a, gaba_b),
    gaba_glu = if_else(type_a == "Glut", glu_b, glu_a),
    gaba_gaba = if_else(type_a == "Glut", gaba_b, gaba_a),
    glu_p = if_else(type_a == "Glut", p_a, p_b),
    gaba_p = if_else(type_a == "Glut", p_b, p_a),
    glu_eir = if_else(type_a == "Glut", eir_a, eir_b),
    gaba_eir = if_else(type_a == "Glut", eir_b, eir_a),
    glu_ov = if_else(type_a == "Glut", ov_a, ov_b),
    gaba_ov = if_else(type_a == "Glut", ov_b, ov_a)
  ) %>%
  mutate(
    I_in_E = abs(gaba_ov - 1) < 1e-12,
    E_in_I = abs(glu_ov - 1) < 1e-12,
    partial = !(I_in_E | E_in_I),
    status = case_when(
      I_in_E & E_in_I ~ "mutual_full",
      I_in_E ~ "I_in_E",
      E_in_I ~ "E_in_I",
      TRUE ~ "partial"
    ),
    asymmetry_abs = abs(glu_ov - gaba_ov),
    size_ratio_gaba_to_glu = gaba_n / glu_n,
    log_size_ratio = log((gaba_n + 0.5) / (glu_n + 0.5)),
    min_n = pmin(glu_n, gaba_n),
    max_n = pmax(glu_n, gaba_n),
    glu_region_count = count_region_proxy(glu_region),
    gaba_region_count = count_region_proxy(gaba_region),
    region_count_diff_glu_minus_gaba = glu_region_count - gaba_region_count,
    glu_logp = -log10(glu_p),
    gaba_logp = -log10(gaba_p),
    glu_glu_frac = glu_glu / glu_n,
    gaba_gaba_frac = gaba_gaba / gaba_n,
    motif = paste(glu_subclass, gaba_subclass, sep = " || ")
  ) %>%
  filter(status != "mutual_full")

# 所有可能 E-I 配对宇宙
all_ei <- all_pairs %>%
  filter(pair_type == "E-I") %>%
  mutate(
    glu_label = if_else(type_a == "Glut", label_a, label_b),
    gaba_label = if_else(type_a == "Glut", label_b, label_a),
    glu_subclass = if_else(type_a == "Glut", subclass_a, subclass_b),
    gaba_subclass = if_else(type_a == "Glut", subclass_b, subclass_a),
    glu_n = if_else(type_a == "Glut", n_a, n_b),
    gaba_n = if_else(type_a == "Glut", n_b, n_a)
  ) %>%
  left_join(
    ei %>% select(pair_key, overlap, jaccard, glu_ov, gaba_ov, I_in_E, E_in_I),
    by = "pair_key"
  ) %>%
  mutate(
    overlap = replace_na(overlap, 0),
    jaccard = replace_na(jaccard, 0),
    glu_ov = replace_na(glu_ov, 0),
    gaba_ov = replace_na(gaba_ov, 0),
    I_in_E = replace_na(I_in_E, FALSE),
    E_in_I = replace_na(E_in_I, FALSE)
  )

# -------------------------------
# 构建所有定向完全包含事件
# -------------------------------
contain_events <- bind_rows(
  pairs %>%
    filter(abs(ov_a - 1) < 1e-12) %>%
    transmute(
      pair_key, slide, layer, stratum,
      contained_label = label_a,
      container_label = label_b,
      contained_type = type_a,
      container_type = type_b,
      contained_subclass = subclass_a,
      container_subclass = subclass_b,
      contained_n = n_a,
      container_n = n_b,
      overlap, jaccard,
      contained_ov = ov_a,
      container_ov = ov_b
    ),
  pairs %>%
    filter(abs(ov_b - 1) < 1e-12) %>%
    transmute(
      pair_key, slide, layer, stratum,
      contained_label = label_b,
      container_label = label_a,
      contained_type = type_b,
      container_type = type_a,
      contained_subclass = subclass_b,
      container_subclass = subclass_a,
      contained_n = n_b,
      container_n = n_a,
      overlap, jaccard,
      contained_ov = ov_b,
      container_ov = ov_a
    )
) %>%
  mutate(
    event_type = paste0(type_to_ei(contained_type), "_in_", type_to_ei(container_type)),
    size_ratio_contained_to_container = contained_n / container_n,
    asymmetry_abs = 1 - container_ov
  )

# -------------------------------
# 表1：数据总览
# -------------------------------
table_01_data_summary <- tribble(
  ~metric, ~value,
  "原始CSV行数", nrow(raw),
  "去重后唯一重叠配对数", nrow(pairs),
  "唯一cluster数", nrow(clusters),
  "slide×layer分层数", n_distinct(clusters$stratum),
  "所有可能配对数（同slide×layer内）", nrow(all_pairs),
  "重叠配对占全部可能配对比例", mean(all_pairs$observed),
  "观测到E-I重叠配对数", sum(pairs$pair_type == "E-I"),
  "观测到E-E重叠配对数", sum(pairs$pair_type == "E-E"),
  "观测到I-I重叠配对数", sum(pairs$pair_type == "I-I"),
  "存在任一方向100%包含的重叠配对数", sum(pairs$any_containment),
  "重叠配对中包含比例", mean(pairs$any_containment)
)
readr::write_csv(table_01_data_summary, file.path(out_dir, "table_01_data_summary.csv"))

# -------------------------------
# 表2-4：pair type 的包含率 / CMH 检验
# -------------------------------
table_02_pairtype_containment_rates_observed <- pairs %>%
  group_by(pair_type) %>%
  summarise(
    n_pairs = n(),
    n_containment = sum(any_containment),
    containment_rate = mean(any_containment),
    .groups = "drop"
  ) %>%
  arrange(desc(containment_rate))
readr::write_csv(table_02_pairtype_containment_rates_observed,
                 file.path(out_dir, "table_02_pairtype_containment_rates_observed.csv"))

table_03_pairtype_rates_all_possible <- all_pairs %>%
  group_by(pair_type) %>%
  summarise(
    n_possible_pairs = n(),
    n_observed_overlaps = sum(observed),
    overlap_rate = mean(observed),
    n_containment = sum(any_containment),
    containment_rate_all_possible = mean(any_containment),
    .groups = "drop"
  )
readr::write_csv(table_03_pairtype_rates_all_possible,
                 file.path(out_dir, "table_03_pairtype_rates_all_possible.csv"))

table_04_cmh_pairtype_tests <- bind_rows(
  cmh_binary_compare(pairs %>% filter(pair_type %in% c("E-I", "I-I")),
                     "pair_type", "E-I", "any_containment", TRUE,
                     comparison_label = "observed_overlaps_only | containment | E-I vs I-I"),
  cmh_binary_compare(pairs %>% filter(pair_type %in% c("E-E", "I-I")),
                     "pair_type", "E-E", "any_containment", TRUE,
                     comparison_label = "observed_overlaps_only | containment | E-E vs I-I"),
  cmh_binary_compare(pairs %>% filter(pair_type %in% c("E-I", "E-E")),
                     "pair_type", "E-I", "any_containment", TRUE,
                     comparison_label = "observed_overlaps_only | containment | E-I vs E-E"),
  cmh_binary_compare(all_pairs %>% filter(pair_type %in% c("E-I", "I-I")),
                     "pair_type", "E-I", "any_containment", TRUE,
                     comparison_label = "all_possible_pairs | containment | E-I vs I-I"),
  cmh_binary_compare(all_pairs %>% filter(pair_type %in% c("E-E", "I-I")),
                     "pair_type", "E-E", "any_containment", TRUE,
                     comparison_label = "all_possible_pairs | containment | E-E vs I-I"),
  cmh_binary_compare(all_pairs %>% filter(pair_type %in% c("E-I", "E-E")),
                     "pair_type", "E-I", "any_containment", TRUE,
                     comparison_label = "all_possible_pairs | containment | E-I vs E-E"),
  cmh_binary_compare(all_pairs %>% filter(pair_type %in% c("E-I", "I-I")),
                     "pair_type", "E-I", "observed", TRUE,
                     comparison_label = "all_possible_pairs | overlap | E-I vs I-I"),
  cmh_binary_compare(all_pairs %>% filter(pair_type %in% c("E-E", "I-I")),
                     "pair_type", "E-E", "observed", TRUE,
                     comparison_label = "all_possible_pairs | overlap | E-E vs I-I"),
  cmh_binary_compare(all_pairs %>% filter(pair_type %in% c("E-I", "E-E")),
                     "pair_type", "E-I", "observed", TRUE,
                     comparison_label = "all_possible_pairs | overlap | E-I vs E-E")
) %>%
  mutate(p_fmt = fmt_p(p))
readr::write_csv(table_04_cmh_pairtype_tests,
                 file.path(out_dir, "table_04_cmh_pairtype_tests.csv"))

# -------------------------------
# 表5-6：E-I 包含方向性（修正layer分母bug）
# -------------------------------
table_05_ei_directionality_overall <- tibble(
  scope = "observed_EI_overlaps",
  n_EI_pairs = nrow(ei),
  I_in_E = sum(ei$I_in_E),
  E_in_I = sum(ei$E_in_I),
  partial = sum(ei$partial),
  I_in_E_rate_among_EI_overlaps = mean(ei$I_in_E),
  E_in_I_rate_among_EI_overlaps = mean(ei$E_in_I),
  I_in_E_among_containment_events = sum(ei$I_in_E) / (sum(ei$I_in_E) + sum(ei$E_in_I)),
  directionality_exact_binom_p = binom.test(sum(ei$I_in_E),
                                            sum(ei$I_in_E) + sum(ei$E_in_I),
                                            p = 0.5,
                                            alternative = "greater")$p.value,
  all_possible_EI_pairs = nrow(all_ei),
  I_in_E_rate_all_possible_EI = mean(all_ei$I_in_E),
  E_in_I_rate_all_possible_EI = mean(all_ei$E_in_I)
) %>%
  mutate(p_fmt = fmt_p(directionality_exact_binom_p))
readr::write_csv(table_05_ei_directionality_overall,
                 file.path(out_dir, "table_05_ei_directionality_overall.csv"))

table_06_ei_directionality_by_layer <- ei %>%
  group_by(layer) %>%
  summarise(
    n_EI_pairs = n(),
    n_I_in_E = sum(I_in_E),
    n_E_in_I = sum(E_in_I),
    n_partial = sum(partial),
    I_in_E_rate_among_EI_overlaps = sum(I_in_E) / n(),
    E_in_I_rate_among_EI_overlaps = sum(E_in_I) / n(),
    I_in_E_among_containment_events = ifelse(sum(I_in_E) + sum(E_in_I) > 0,
                                             sum(I_in_E) / (sum(I_in_E) + sum(E_in_I)),
                                             NA_real_),
    directionality_p = ifelse(sum(I_in_E) + sum(E_in_I) > 0,
                              binom.test(sum(I_in_E),
                                         sum(I_in_E) + sum(E_in_I),
                                         p = 0.5,
                                         alternative = "greater")$p.value,
                              NA_real_),
    .groups = "drop"
  ) %>%
  arrange(desc(I_in_E_rate_among_EI_overlaps)) %>%
  mutate(p_fmt = fmt_p(directionality_p))
readr::write_csv(table_06_ei_directionality_by_layer,
                 file.path(out_dir, "table_06_ei_directionality_by_layer.csv"))

# -------------------------------
# 表7：I_in_E vs partial 的 stratum 配对 Wilcoxon
# -------------------------------
metrics_iine_partial <- c(
  "jaccard", "overlap", "glu_ov", "gaba_ov", "glu_n", "gaba_n",
  "size_ratio_gaba_to_glu", "asymmetry_abs",
  "glu_region_count", "gaba_region_count", "region_count_diff_glu_minus_gaba",
  "glu_logp", "gaba_logp", "glu_eir", "gaba_eir"
)

table_07_iine_vs_partial_paired_wilcoxon <- purrr::map_dfr(
  metrics_iine_partial,
  ~ paired_stratum_wilcox(
    ei %>% filter(status %in% c("I_in_E", "partial")),
    group_col = "status",
    metric = .x,
    group1 = "I_in_E",
    group0 = "partial"
  )
) %>%
  mutate(p_fmt = fmt_p(p))
readr::write_csv(table_07_iine_vs_partial_paired_wilcoxon,
                 file.path(out_dir, "table_07_iine_vs_partial_paired_wilcoxon.csv"))

# -------------------------------
# 表8：I_in_E 与其他完全包含事件的比较
# -------------------------------
pairwise_event_compare <- function(metric, a, b) {
  paired_stratum_wilcox(
    contain_events %>% filter(event_type %in% c(a, b)),
    group_col = "event_type",
    metric = metric,
    group1 = a,
    group0 = b
  ) %>%
    mutate(comparison = paste(a, "vs", b))
}

table_08_containment_eventtype_pairwise_paired_wilcoxon <- bind_rows(
  purrr::map_dfr(
    c("container_n", "contained_n", "size_ratio_contained_to_container", "container_ov", "jaccard", "overlap"),
    ~ pairwise_event_compare(.x, "I_in_E", "E_in_E")
  ),
  purrr::map_dfr(
    c("container_n", "contained_n", "size_ratio_contained_to_container", "container_ov", "jaccard", "overlap"),
    ~ pairwise_event_compare(.x, "I_in_E", "I_in_I")
  ),
  purrr::map_dfr(
    c("container_n", "contained_n", "size_ratio_contained_to_container", "container_ov", "jaccard", "overlap"),
    ~ pairwise_event_compare(.x, "I_in_E", "E_in_I")
  )
) %>%
  mutate(p_fmt = fmt_p(p))
readr::write_csv(table_08_containment_eventtype_pairwise_paired_wilcoxon,
                 file.path(out_dir, "table_08_containment_eventtype_pairwise_paired_wilcoxon.csv"))

# -------------------------------
# 表9-10：多重嵌套 / multiplicity
# -------------------------------
host_counts_iine <- ei %>% filter(I_in_E) %>% count(glu_label, name = "n_gaba_contained")
guest_counts_iine <- ei %>% filter(I_in_E) %>% count(gaba_label, name = "n_glu_hosts")
host_counts_eini <- ei %>% filter(E_in_I) %>% count(gaba_label, name = "n_glu_contained")
guest_counts_eini <- ei %>% filter(E_in_I) %>% count(glu_label, name = "n_gaba_hosts")

table_09_multiplicity_summary <- bind_rows(
  tibble(
    type = "GLU host (I_in_E)",
    n_objects = nrow(host_counts_iine),
    median_partners = median(host_counts_iine$n_gaba_contained),
    mean_partners = mean(host_counts_iine$n_gaba_contained),
    prop_multi = mean(host_counts_iine$n_gaba_contained >= 2),
    max_partners = max(host_counts_iine$n_gaba_contained)
  ),
  tibble(
    type = "GABA host (E_in_I)",
    n_objects = nrow(host_counts_eini),
    median_partners = median(host_counts_eini$n_glu_contained),
    mean_partners = mean(host_counts_eini$n_glu_contained),
    prop_multi = mean(host_counts_eini$n_glu_contained >= 2),
    max_partners = max(host_counts_eini$n_glu_contained)
  ),
  tibble(
    type = "GABA guest (I_in_E)",
    n_objects = nrow(guest_counts_iine),
    median_partners = median(guest_counts_iine$n_glu_hosts),
    mean_partners = mean(guest_counts_iine$n_glu_hosts),
    prop_multi = mean(guest_counts_iine$n_glu_hosts >= 2),
    max_partners = max(guest_counts_iine$n_glu_hosts)
  ),
  tibble(
    type = "GLU guest (E_in_I)",
    n_objects = nrow(guest_counts_eini),
    median_partners = median(guest_counts_eini$n_gaba_hosts),
    mean_partners = mean(guest_counts_eini$n_gaba_hosts),
    prop_multi = mean(guest_counts_eini$n_gaba_hosts >= 2),
    max_partners = max(guest_counts_eini$n_gaba_hosts)
  )
)
readr::write_csv(table_09_multiplicity_summary,
                 file.path(out_dir, "table_09_multiplicity_summary.csv"))

ft_host <- fisher.test(matrix(c(
  sum(host_counts_iine$n_gaba_contained >= 2), sum(host_counts_iine$n_gaba_contained < 2),
  sum(host_counts_eini$n_glu_contained >= 2), sum(host_counts_eini$n_glu_contained < 2)
), nrow = 2, byrow = TRUE))

ft_guest <- fisher.test(matrix(c(
  sum(guest_counts_iine$n_glu_hosts >= 2), sum(guest_counts_iine$n_glu_hosts < 2),
  sum(guest_counts_eini$n_gaba_hosts >= 2), sum(guest_counts_eini$n_gaba_hosts < 2)
), nrow = 2, byrow = TRUE))

table_10_multiplicity_fisher_tests <- bind_rows(
  tibble(
    comparison = "Host multiplicity: GLU hosts in I_in_E vs GABA hosts in E_in_I",
    odds_ratio = unname(ft_host$estimate),
    p = ft_host$p.value,
    group1_multi_fraction = mean(host_counts_iine$n_gaba_contained >= 2),
    group2_multi_fraction = mean(host_counts_eini$n_glu_contained >= 2)
  ),
  tibble(
    comparison = "Contained multiplicity: GABA guests in I_in_E vs GLU guests in E_in_I",
    odds_ratio = unname(ft_guest$estimate),
    p = ft_guest$p.value,
    group1_multi_fraction = mean(guest_counts_iine$n_glu_hosts >= 2),
    group2_multi_fraction = mean(guest_counts_eini$n_gaba_hosts >= 2)
  )
) %>%
  mutate(p_fmt = fmt_p(p))
readr::write_csv(table_10_multiplicity_fisher_tests,
                 file.path(out_dir, "table_10_multiplicity_fisher_tests.csv"))

# -------------------------------
# 表11-19：motif / subclass 富集
# -------------------------------
ei_iine_partial <- ei %>% filter(status %in% c("I_in_E", "partial"))

table_11_iine_vs_partial_motif_enrichment_full <- fisher_feature_enrichment(
  ei_iine_partial,
  feature_col = "motif",
  event_col = "status",
  event_label = "I_in_E",
  baseline_label = "partial",
  min_count = 10
)
readr::write_csv(table_11_iine_vs_partial_motif_enrichment_full,
                 file.path(out_dir, "table_11_iine_vs_partial_motif_enrichment_full.csv"))

table_12_top_enriched_motifs_iine_vs_partial <- table_11_iine_vs_partial_motif_enrichment_full %>%
  filter(odds_ratio > 1) %>%
  arrange(fdr, desc(odds_ratio)) %>%
  slice_head(n = 20)
readr::write_csv(table_12_top_enriched_motifs_iine_vs_partial,
                 file.path(out_dir, "table_12_top_enriched_motifs_iine_vs_partial.csv"))

table_13_top_depleted_motifs_iine_vs_partial <- table_11_iine_vs_partial_motif_enrichment_full %>%
  filter(odds_ratio < 1) %>%
  arrange(fdr, odds_ratio) %>%
  slice_head(n = 20)
readr::write_csv(table_13_top_depleted_motifs_iine_vs_partial,
                 file.path(out_dir, "table_13_top_depleted_motifs_iine_vs_partial.csv"))

table_14_glu_subclass_enrichment_pairlevel_iine_vs_partial <- fisher_feature_enrichment(
  ei_iine_partial,
  feature_col = "glu_subclass",
  event_col = "status",
  event_label = "I_in_E",
  baseline_label = "partial",
  min_count = 20
)
readr::write_csv(table_14_glu_subclass_enrichment_pairlevel_iine_vs_partial,
                 file.path(out_dir, "table_14_glu_subclass_enrichment_pairlevel_iine_vs_partial.csv"))

table_15_gaba_subclass_enrichment_pairlevel_iine_vs_partial <- fisher_feature_enrichment(
  ei_iine_partial,
  feature_col = "gaba_subclass",
  event_col = "status",
  event_label = "I_in_E",
  baseline_label = "partial",
  min_count = 20
)
readr::write_csv(table_15_gaba_subclass_enrichment_pairlevel_iine_vs_partial,
                 file.path(out_dir, "table_15_gaba_subclass_enrichment_pairlevel_iine_vs_partial.csv"))

# cluster 层面统计
gaba_stats <- ei %>%
  group_by(gaba_label) %>%
  summarise(
    slide = first(slide),
    layer = first(layer),
    gaba_subclass = first(gaba_subclass),
    gaba_n = first(gaba_n),
    gaba_eir = first(gaba_eir),
    ei_partner_count = n(),
    full_container_count = sum(I_in_E),
    any_I_in_E = any(I_in_E),
    any_E_in_I = any(E_in_I),
    median_glu_partner_n = median(glu_n),
    median_jaccard = median(jaccard),
    .groups = "drop"
  ) %>%
  mutate(
    full_container_count_ge2 = full_container_count >= 2
  )

glu_stats <- ei %>%
  group_by(glu_label) %>%
  summarise(
    slide = first(slide),
    layer = first(layer),
    glu_subclass = first(glu_subclass),
    glu_n = first(glu_n),
    glu_eir = first(glu_eir),
    ei_partner_count = n(),
    full_containee_count = sum(I_in_E),
    any_host_I_in_E = any(I_in_E),
    any_contained_E_in_I = any(E_in_I),
    median_gaba_partner_n = median(gaba_n),
    median_jaccard = median(jaccard),
    .groups = "drop"
  ) %>%
  mutate(
    multi_host = full_containee_count >= 2
  )

readr::write_csv(gaba_stats, file.path(out_dir, "table_16_gaba_cluster_stats.csv"))
readr::write_csv(glu_stats, file.path(out_dir, "table_17_glu_cluster_stats.csv"))

table_18_gaba_cluster_subclass_enrichment_ever_contained <- cluster_subclass_enrichment(
  gaba_stats, subclass_col = "gaba_subclass", group_col = "any_I_in_E", min_count = 20
)
readr::write_csv(table_18_gaba_cluster_subclass_enrichment_ever_contained,
                 file.path(out_dir, "table_18_gaba_cluster_subclass_enrichment_ever_contained.csv"))

table_19_glu_cluster_subclass_enrichment_ever_hosting <- cluster_subclass_enrichment(
  glu_stats, subclass_col = "glu_subclass", group_col = "any_host_I_in_E", min_count = 20
)
readr::write_csv(table_19_glu_cluster_subclass_enrichment_ever_hosting,
                 file.path(out_dir, "table_19_glu_cluster_subclass_enrichment_ever_hosting.csv"))

table_20_glu_host_subclass_enrichment_multi_containee <- cluster_subclass_enrichment(
  glu_stats %>% filter(any_host_I_in_E),
  subclass_col = "glu_subclass", group_col = "multi_host", min_count = 10
)
readr::write_csv(table_20_glu_host_subclass_enrichment_multi_containee,
                 file.path(out_dir, "table_20_glu_host_subclass_enrichment_multi_containee.csv"))

table_21_gaba_guest_subclass_enrichment_multi_host <- cluster_subclass_enrichment(
  gaba_stats %>% filter(any_I_in_E),
  subclass_col = "gaba_subclass", group_col = "full_container_count_ge2", min_count = 10
)
readr::write_csv(table_21_gaba_guest_subclass_enrichment_multi_host,
                 file.path(out_dir, "table_21_gaba_guest_subclass_enrichment_multi_host.csv"))

# -------------------------------
# 表22-25：示例性极端案例
# -------------------------------
table_22_top_I_in_E_by_overlap <- ei %>%
  filter(I_in_E) %>%
  arrange(desc(overlap)) %>%
  select(slide, layer, glu_label, gaba_label, glu_subclass, gaba_subclass,
         glu_n, gaba_n, overlap, jaccard, glu_ov) %>%
  slice_head(n = 20)
readr::write_csv(table_22_top_I_in_E_by_overlap,
                 file.path(out_dir, "table_22_top_I_in_E_by_overlap.csv"))

table_23_top_I_in_E_by_extreme_small_ratio <- ei %>%
  filter(I_in_E) %>%
  arrange(size_ratio_gaba_to_glu) %>%
  select(slide, layer, glu_label, gaba_label, glu_subclass, gaba_subclass,
         glu_n, gaba_n, overlap, jaccard, glu_ov, size_ratio_gaba_to_glu) %>%
  slice_head(n = 20)
readr::write_csv(table_23_top_I_in_E_by_extreme_small_ratio,
                 file.path(out_dir, "table_23_top_I_in_E_by_extreme_small_ratio.csv"))

table_24_top_GLU_hosts_by_number_of_contained_GABA <- host_counts_iine %>%
  arrange(desc(n_gaba_contained))
readr::write_csv(table_24_top_GLU_hosts_by_number_of_contained_GABA,
                 file.path(out_dir, "table_24_top_GLU_hosts_by_number_of_contained_GABA.csv"))

table_25_top_GABA_guests_by_number_of_GLU_hosts <- guest_counts_iine %>%
  arrange(desc(n_glu_hosts))
readr::write_csv(table_25_top_GABA_guests_by_number_of_GLU_hosts,
                 file.path(out_dir, "table_25_top_GABA_guests_by_number_of_GLU_hosts.csv"))

# -------------------------------
# 表26-29：cluster 层面组间比较
# -------------------------------
table_26_gaba_cluster_group_comparisons <- bind_rows(
  group_compare_wilcox(gaba_stats, "any_I_in_E", "gaba_n"),
  group_compare_wilcox(gaba_stats, "any_I_in_E", "gaba_eir"),
  group_compare_wilcox(gaba_stats, "any_I_in_E", "ei_partner_count"),
  group_compare_wilcox(gaba_stats, "any_I_in_E", "median_glu_partner_n"),
  group_compare_wilcox(gaba_stats, "any_I_in_E", "median_jaccard")
) %>%
  mutate(group = "GABA ever fully contained by GLU vs never",
         p_fmt = fmt_p(p))
readr::write_csv(table_26_gaba_cluster_group_comparisons,
                 file.path(out_dir, "table_26_gaba_cluster_group_comparisons.csv"))

table_27_glu_cluster_group_comparisons <- bind_rows(
  group_compare_wilcox(glu_stats, "any_host_I_in_E", "glu_n"),
  group_compare_wilcox(glu_stats, "any_host_I_in_E", "glu_eir"),
  group_compare_wilcox(glu_stats, "any_host_I_in_E", "ei_partner_count"),
  group_compare_wilcox(glu_stats, "any_host_I_in_E", "median_gaba_partner_n"),
  group_compare_wilcox(glu_stats, "any_host_I_in_E", "median_jaccard")
) %>%
  mutate(group = "GLU ever fully contains GABA vs never",
         p_fmt = fmt_p(p))
readr::write_csv(table_27_glu_cluster_group_comparisons,
                 file.path(out_dir, "table_27_glu_cluster_group_comparisons.csv"))

# -------------------------------
# 表28-29：描述性统计
# -------------------------------
table_28_status_descriptive_medians <- ei %>%
  group_by(status) %>%
  summarise(
    n = n(),
    jaccard_median = median(jaccard),
    overlap_median = median(overlap),
    glu_n_median = median(glu_n),
    gaba_n_median = median(gaba_n),
    size_ratio_median = median(size_ratio_gaba_to_glu),
    asymmetry_median = median(asymmetry_abs),
    glu_ov_median = median(glu_ov),
    gaba_ov_median = median(gaba_ov),
    glu_region_count_median = median(glu_region_count),
    gaba_region_count_median = median(gaba_region_count),
    glu_logp_median = median(glu_logp),
    gaba_logp_median = median(gaba_logp),
    .groups = "drop"
  )
readr::write_csv(table_28_status_descriptive_medians,
                 file.path(out_dir, "table_28_status_descriptive_medians.csv"))

table_29_eventtype_descriptive_medians <- contain_events %>%
  group_by(event_type) %>%
  summarise(
    n = n(),
    contained_n_median = median(contained_n),
    container_n_median = median(container_n),
    size_ratio_median = median(size_ratio_contained_to_container),
    overlap_median = median(overlap),
    jaccard_median = median(jaccard),
    container_ov_median = median(container_ov),
    .groups = "drop"
  )
readr::write_csv(table_29_eventtype_descriptive_medians,
                 file.path(out_dir, "table_29_eventtype_descriptive_medians.csv"))

# -------------------------------
# 表30-31：logistic 敏感性分析
# -------------------------------
glm_df <- ei %>%
  filter(status %in% c("I_in_E", "partial")) %>%
  mutate(
    y = if_else(status == "I_in_E", 1, 0),
    log_glu_n = log1p(glu_n),
    log_gaba_n = log1p(gaba_n)
  )

fit_layer <- glm(
  y ~ log_glu_n + log_gaba_n + region_count_diff_glu_minus_gaba +
    glu_logp + gaba_logp + factor(layer),
  family = binomial(),
  data = glm_df
)

fit_slide_layer <- glm(
  y ~ log_glu_n + log_gaba_n + region_count_diff_glu_minus_gaba +
    glu_logp + gaba_logp + factor(layer) + factor(slide),
  family = binomial(),
  data = glm_df
)

table_30_logistic_sensitivity_layer_adjusted <- broom::tidy(fit_layer, conf.int = TRUE) %>%
  mutate(
    odds_ratio = exp(estimate),
    or_low = exp(conf.low),
    or_high = exp(conf.high),
    p_fmt = fmt_p(p.value)
  )
readr::write_csv(table_30_logistic_sensitivity_layer_adjusted,
                 file.path(out_dir, "table_30_logistic_sensitivity_layer_adjusted.csv"))

table_31_logistic_sensitivity_slide_layer_adjusted <- broom::tidy(fit_slide_layer, conf.int = TRUE) %>%
  mutate(
    odds_ratio = exp(estimate),
    or_low = exp(conf.low),
    or_high = exp(conf.high),
    p_fmt = fmt_p(p.value)
  )
readr::write_csv(table_31_logistic_sensitivity_slide_layer_adjusted,
                 file.path(out_dir, "table_31_logistic_sensitivity_slide_layer_adjusted.csv"))

# -------------------------------
# 图形主题（尽量贴近参考代码）
# -------------------------------
plot_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "#4D4D4D"),
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black")
  )

pair_cols <- c("E-E" = "#4E79A7", "E-I" = "#9C6BC4", "I-I" = "#59A14F")
status_cols <- c("Partial" = "#BDBDBD", "I in E" = "#CC5A5F", "E in I" = "#5B84B1")

# -------------------------------
# fig01：保留参考代码展示内容，只美化
# -------------------------------
fig01 <- table_02_pairtype_containment_rates_observed %>%
  mutate(pair_type = factor(pair_type, levels = c("E-E", "E-I", "I-I"))) %>%
  ggplot(aes(x = pair_type, y = containment_rate, fill = pair_type)) +
  geom_col(width = 0.70) +
  geom_text(aes(label = scales::percent(containment_rate, accuracy = 0.1)),
            vjust = -0.2, size = 5.0) +
  scale_fill_manual(values = pair_cols) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.10))) +
  labs(
    x = "Pair type",
    y = "Containment rate among observed overlaps",
    title = "Containment is most frequent in E-I overlaps",
    subtitle = "Unique undirected overlap pairs, summarized by interaction class"
  ) +
  plot_theme +
  theme(legend.position = "none")
ggsave(file.path(fig_dir, "fig01_pairtype_containment_rate_observed.png"), fig01,
       width = 7.0, height = 5.2, dpi = 400, bg = "white")

# -------------------------------
# fig02：仍是count bar，只修顶部裁切与格式
# -------------------------------
fig02_df <- ei %>%
  count(status) %>%
  filter(status %in% c("partial", "I_in_E", "E_in_I")) %>%
  mutate(
    status_label = factor(c("Partial", "I in E", "E in I")[match(status, c("partial", "I_in_E", "E_in_I"))],
                          levels = c("Partial", "I in E", "E in I")),
    pct = n / sum(n),
    label = paste0(scales::comma(n), "\n(", scales::percent(pct, accuracy = 0.1), ")")
  )

fig02 <- ggplot(fig02_df, aes(x = status_label, y = n, fill = status_label)) +
  geom_col(width = 0.70) +
  geom_text(aes(label = label), vjust = -0.22, size = 5.0, lineheight = 0.95) +
  scale_fill_manual(values = status_cols) +
  scale_y_continuous(labels = scales::comma_format(),
                     expand = expansion(mult = c(0, 0.16))) +
  labs(
    x = "E-I interaction status",
    y = "Number of E-I overlap pairs",
    title = "Directional asymmetry within E-I overlaps",
    subtitle = "Most E-I overlaps are partial, but one-sided full containment strongly favors I in E"
  ) +
  coord_cartesian(clip = "off") +
  plot_theme +
  theme(legend.position = "none")
ggsave(file.path(fig_dir, "fig02_ei_status_counts.png"), fig02,
       width = 7.2, height = 5.4, dpi = 400, bg = "white")

# -------------------------------
# fig03：保留参考代码展示对象（scatter），只做轻量清理
# -------------------------------
fig03 <- contain_events %>%
  mutate(event_type = factor(event_type, levels = c("I_in_E", "E_in_E", "I_in_I", "E_in_I"))) %>%
  ggplot(aes(x = container_n, y = contained_n, shape = event_type, color = event_type)) +
  geom_point(alpha = 0.35, size = 1.8) +
  scale_x_log10(labels = scales::comma_format()) +
  scale_y_log10(labels = scales::comma_format()) +
  scale_color_manual(values = c("I_in_E" = "#D66A6A",
                                "E_in_E" = "#E4A06F",
                                "I_in_I" = "#8BBF85",
                                "E_in_I" = "#6C92C4")) +
  labs(
    x = "Container cluster size (cells, log scale)",
    y = "Contained cluster size (cells, log scale)",
    title = "Containment events differ in host-domain size hierarchy",
    subtitle = "Each point is one full-containment event"
  ) +
  plot_theme
ggsave(file.path(fig_dir, "fig03_containment_event_size_scatter.png"), fig03,
       width = 7.2, height = 5.8, dpi = 400, bg = "white")

# -------------------------------
# fig04：保留参考代码表达，修正layer rate计算
# -------------------------------
fig04 <- table_06_ei_directionality_by_layer %>%
  filter(layer != "layer 6b") %>%
  ggplot(aes(x = reorder(layer, -I_in_E_rate_among_EI_overlaps),
             y = I_in_E_rate_among_EI_overlaps,
             fill = I_in_E_rate_among_EI_overlaps)) +
  geom_col(width = 0.70) +
  geom_text(aes(label = scales::percent(I_in_E_rate_among_EI_overlaps, accuracy = 0.1)),
            vjust = -0.2, size = 5.0) +
  scale_fill_gradient(low = "#E6C7DA", high = "#E0117F", guide = "none") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.10))) +
  labs(
    x = "Layer",
    y = "I_in_E rate among E-I overlaps",
    title = "I_in_E enrichment is strongest in layer 2/3",
    subtitle = "Rates are computed within each layer among observed E-I overlaps"
  ) +
  plot_theme
ggsave(file.path(fig_dir, "fig04_iine_rate_by_layer.png"), fig04,
       width = 7.2, height = 4.8, dpi = 400, bg = "white")

# -------------------------------
# fig05：严格按参考代码的boxplot主题，不再画灰点
# 加显著性差异 + n与median标注
# -------------------------------
fig05_df <- ei %>%
  filter(status %in% c("partial", "I_in_E", "E_in_I")) %>%
  mutate(
    status = factor(status, levels = c("partial", "I_in_E", "E_in_I")),
    status_label = factor(c("Partial", "I in E", "E in I")[match(status, c("partial", "I_in_E", "E_in_I"))],
                          levels = c("Partial", "I in E", "E in I"))
  )

# pairwise Wilcoxon + BH
p_raw <- c(
  wilcox.test(size_ratio_gaba_to_glu ~ status_label,
              data = fig05_df %>% filter(status_label %in% c("Partial", "I in E")),
              exact = FALSE)$p.value,
  wilcox.test(size_ratio_gaba_to_glu ~ status_label,
              data = fig05_df %>% filter(status_label %in% c("I in E", "E in I")),
              exact = FALSE)$p.value,
  wilcox.test(size_ratio_gaba_to_glu ~ status_label,
              data = fig05_df %>% filter(status_label %in% c("Partial", "E in I")),
              exact = FALSE)$p.value
)
p_adj <- p.adjust(p_raw, method = "BH")

max_y_fig05 <- max(fig05_df$size_ratio_gaba_to_glu, na.rm = TRUE)

sig_df <- tibble(
  x    = c(1, 2, 1),
  xend = c(2, 3, 3),
  y    = c(max_y_fig05 * 1.90, max_y_fig05 * 3.20, max_y_fig05 * 5.40),
  yend = c(max_y_fig05 * 1.70, max_y_fig05 * 2.90, max_y_fig05 * 4.90),
  label = c(fmt_p_plot(p_adj[1]), fmt_p_plot(p_adj[2]), fmt_p_plot(p_adj[3]))
)

box_label_df <- fig05_df %>%
  group_by(status_label) %>%
  summarise(
    n = n(),
    med = median(size_ratio_gaba_to_glu, na.rm = TRUE),
    q3 = quantile(size_ratio_gaba_to_glu, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    x = c(1, 2, 3),
    y = q3 * c(1.58, 1.52, 1.58),
    label = paste0("n=", scales::comma(n), "\nmedian=", sprintf("%.2f", med))
  )

fig05 <- ggplot(fig05_df, aes(x = status_label, y = size_ratio_gaba_to_glu, fill = status_label)) +
  geom_boxplot(outlier.shape = NA, width = 0.56, linewidth = 0.7) +
  geom_text(data = box_label_df,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 3.8, lineheight = 0.92) +
  geom_segment(data = sig_df,
               aes(x = x, xend = xend, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.55) +
  geom_segment(data = sig_df,
               aes(x = x, xend = x, y = y, yend = yend),
               inherit.aes = FALSE, linewidth = 0.55) +
  geom_segment(data = sig_df,
               aes(x = xend, xend = xend, y = y, yend = yend),
               inherit.aes = FALSE, linewidth = 0.55) +
  geom_text(data = sig_df,
            aes(x = (x + xend) / 2, y = y * 1.15, label = label),
            inherit.aes = FALSE, size = 3.7) +
  scale_fill_manual(values = status_cols) +
  scale_y_log10(
    labels = scales::label_number(accuracy = 0.01),
    limits = c(0.01, max_y_fig05 * 8)
  ) +
  labs(
    x = "E-I interaction status",
    y = "GABA / GLU size ratio (log scale)",
    title = "Containment direction tracks guest-to-host size ratio",
    subtitle = "I in E events concentrate at smaller GABA-to-GLU size ratios than partial overlaps"
  ) +
  coord_cartesian(clip = "off") +
  plot_theme +
  theme(
    legend.position = "none",
    plot.margin = margin(15, 20, 12, 12)
  )
ggsave(file.path(fig_dir, "fig05_size_ratio_boxplot.png"), fig05,
       width = 8.0, height = 6.0, dpi = 400, bg = "white")

# -------------------------------
# fig06：仍展示“host multiplicity distribution”
# 但把参考代码的重叠hist改成更适合整数分布的dodged离散柱图
# 展示对象不变：每种host类型在不同n_partner上的host cluster数量
# -------------------------------
host_plot_df <- bind_rows(
  host_counts_iine %>% transmute(kind = "GLU hosts in I_in_E", n_partner = n_gaba_contained),
  host_counts_eini %>% transmute(kind = "GABA hosts in E_in_I", n_partner = n_glu_contained)
)

fig06_df <- host_plot_df %>%
  count(kind, n_partner, name = "n") %>%
  mutate(kind = factor(kind, levels = c("GABA hosts in E_in_I", "GLU hosts in I_in_E")))

fig06 <- ggplot(fig06_df, aes(x = factor(n_partner), y = n, fill = kind)) +
  geom_col(position = position_dodge(width = 0.76), width = 0.68) +
  scale_fill_manual(values = c("GABA hosts in E_in_I" = "#6C92C4",
                               "GLU hosts in I_in_E" = "#CC6B6B")) +
  labs(
    x = "Number of contained partners per host cluster",
    y = "Number of host clusters",
    title = "I_in_E hosts more often contain multiple guests",
    subtitle = paste0(
      "GLU hosts in I_in_E: n=", scales::comma(nrow(host_counts_iine)),
      ", multi-guest rate (≥2)=", scales::percent(mean(host_counts_iine$n_gaba_contained >= 2), accuracy = 0.1),
      "; GABA hosts in E_in_I: n=", scales::comma(nrow(host_counts_eini)),
      ", multi-guest rate (≥2)=", scales::percent(mean(host_counts_eini$n_glu_contained >= 2), accuracy = 0.1),
      "; Fisher ", fmt_p_plot(ft_host$p.value)
    )
  ) +
  plot_theme
ggsave(file.path(fig_dir, "fig06_host_multiplicity_hist.png"), fig06,
       width = 8.2, height = 5.6, dpi = 400, bg = "white")

message("Analysis complete. Outputs written to: ", out_dir)