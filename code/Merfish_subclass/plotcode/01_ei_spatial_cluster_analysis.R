rm(list = ls()); gc()
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
})

# ============================================================
# Top-journal style visualization for neuron-neuron partner analysis
# Working directory: E:/zaw/2603
# ============================================================

# -----------------------------
# 路径
# -----------------------------
base_dir <- "E:/zaw/2603"
path_txt <- file.path(base_dir, "neuron-neuron-partner.txt")
path_csv <- file.path(base_dir, "neuron-neuron-partner.csv")
path_out <- file.path(base_dir, "R_output_01")
dir.create(path_out, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 读入函数
# -----------------------------
read_partner_file <- function(txt_path, csv_path) {
  if (file.exists(txt_path)) {
    message("Reading TXT: ", txt_path)
    read_tsv(txt_path, show_col_types = FALSE)
  } else if (file.exists(csv_path)) {
    message("Reading CSV: ", csv_path)
    read_csv(csv_path, show_col_types = FALSE)
  } else {
    stop("Cannot find input file. Please place neuron-neuron-partner.txt or neuron-neuron-partner.csv in: ", base_dir)
  }
}

# -----------------------------
# 辅助函数
# -----------------------------
subclass_short <- function(x) {
  x %>%
    str_remove("^\\d+\\s+") %>%
    str_replace(" Gaba$", "") %>%
    str_replace(" Glut$", "")
}

fmt_p <- function(p) {
  if (length(p) > 1) {
    return(vapply(p, fmt_p, character(1)))
  }
  if (is.na(p)) return("NA")
  if (p == 0 || p < 1e-300) return("<1e-300")
  if (p < 1e-3) return(formatC(p, format = "e", digits = 2))
  sprintf("%.3f", p)
}

paired_compare <- function(piv, a, b) {
  tmp <- piv %>% dplyr::select(all_of(c(a, b))) %>% na.omit()
  if (nrow(tmp) == 0) {
    return(tibble(n = 0, p = NA_real_, median_diff = NA_real_, mean_diff = NA_real_))
  }
  w <- wilcox.test(tmp[[a]], tmp[[b]], paired = TRUE, exact = FALSE)
  tibble(
    n = nrow(tmp),
    p = w$p.value,
    median_diff = median(tmp[[a]] - tmp[[b]], na.rm = TRUE),
    mean_diff = mean(tmp[[a]] - tmp[[b]], na.rm = TRUE)
  )
}

save_pub <- function(p, filename, width = 7, height = 5) {
  ggsave(
    filename = file.path(path_out, paste0(filename, ".pdf")),
    plot = p, width = width, height = height,
    device = cairo_pdf, dpi = 600, bg = "white"
  )
  ggsave(
    filename = file.path(path_out, paste0(filename, ".png")),
    plot = p, width = width, height = height,
    dpi = 600, bg = "white"
  )
}

theme_pub <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1.5, hjust = 0.5, lineheight = 1.08),
      plot.subtitle = element_text(size = base_size - 0.4, hjust = 0.5, color = "grey20", lineheight = 1.05),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text = element_text(color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      axis.ticks.length = unit(0.18, "cm"),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      panel.border = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      plot.margin = margin(10, 16, 10, 10)
    )
}

pair_cols <- c(
  "Glut-Glut" = "#D95F02",
  "E-I" = "#1B9E77",
  "Gaba-Gaba" = "#7570B3"
)

# -----------------------------
# 读入数据
# -----------------------------
df <- read_partner_file(path_txt, path_csv)

# -----------------------------
# 去重: A-B 与 B-A 视为同一无向 pair
# -----------------------------
df <- df %>%
  mutate(
    pair_key_sorted = if_else(
      `cluster.1_label` < `cluster.2_label`,
      paste(`cluster.1_label`, `cluster.2_label`, sep = "||"),
      paste(`cluster.2_label`, `cluster.1_label`, sep = "||")
    )
  )

dedup <- df %>%
  arrange(pair_key_sorted, `cluster.1_label`, `cluster.2_label`) %>%
  distinct(pair_key_sorted, .keep_all = TRUE) %>%
  mutate(
    pair_type_unordered = case_when(
      `cluster.1_cell_Neruon_type` == "Glut" & `cluster.2_cell_Neruon_type` == "Glut" ~ "Glut-Glut",
      `cluster.1_cell_Neruon_type` == "Gaba" & `cluster.2_cell_Neruon_type` == "Gaba" ~ "Gaba-Gaba",
      TRUE ~ "E-I"
    ),
    max_containment = pmax(`cluster.1.overlap.percent`, `cluster.2.overlap.percent`),
    min_containment = pmin(`cluster.1.overlap.percent`, `cluster.2.overlap.percent`),
    sub1_short = subclass_short(`cluster.1_subclass`),
    sub2_short = subclass_short(`cluster.2_subclass`)
  )

# -----------------------------
# cluster 元数据
# -----------------------------
clusters1 <- dedup %>%
  transmute(
    label       = `cluster.1_label`,
    slide       = `cluster.1_slide`,
    layer       = `cluster.1_layer`,
    region      = `cluster.1_region`,
    cell_type   = `cluster.1_cell_Neruon_type`,
    subclass    = `cluster.1_subclass`,
    total_cell_num = `cluster.1_total_cell_num`,
    Glut_num    = `cluster.1_Glut_Neruon_cell_ids_num`,
    GABA_num    = `cluster.1_GABA_Neruon_cell_ids_num`,
    cauchy_p    = `cluster.1_cauchy_combination_p`,
    E_I_Ratio   = `cluster.1_E_I_Ratio`
  )

clusters2 <- dedup %>%
  transmute(
    label       = `cluster.2_label`,
    slide       = `cluster.2_slide`,
    layer       = `cluster.2_layer`,
    region      = `cluster.2_region`,
    cell_type   = `cluster.2_cell_Neruon_type`,
    subclass    = `cluster.2_subclass`,
    total_cell_num = `cluster.2_total_cell_num`,
    Glut_num    = `cluster.2_Glut_Neruon_cell_ids_num`,
    GABA_num    = `cluster.2_GABA_Neruon_cell_ids_num`,
    cauchy_p    = `cluster.2_cauchy_combination_p`,
    E_I_Ratio   = `cluster.2_E_I_Ratio`
  )

clusters <- bind_rows(clusters1, clusters2) %>%
  distinct(label, .keep_all = TRUE) %>%
  mutate(
    gaba_frac = GABA_num / (GABA_num + Glut_num),
    glut_frac = Glut_num / (GABA_num + Glut_num),
    sub_short = subclass_short(subclass)
  )

# -----------------------------
# 数据总览
# -----------------------------
overview <- tibble(
  directed_pairs    = nrow(df),
  unique_pairs      = nrow(dedup),
  unique_clusters   = nrow(clusters),
  slides            = n_distinct(clusters$slide),
  slide_layer_units = nrow(distinct(clusters, slide, layer)),
  glut_clusters     = sum(clusters$cell_type == "Glut"),
  gaba_clusters     = sum(clusters$cell_type == "Gaba")
)
write_csv(overview, file.path(path_out, "dataset_overview.csv"))

# -----------------------------
# Result block 1: pair type summary
# -----------------------------
pair_summary <- dedup %>%
  group_by(pair_type_unordered) %>%
  summarise(
    n = n(),
    mean_j = mean(jaccard, na.rm = TRUE),
    median_j = median(jaccard, na.rm = TRUE),
    q75_j = quantile(jaccard, 0.75, na.rm = TRUE),
    prop_cont95 = mean(max_containment >= 0.95, na.rm = TRUE),
    prop_cont99 = mean(max_containment >= 0.99, na.rm = TRUE),
    prop_cont100 = mean(max_containment >= 1.00, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(pair_type_unordered = factor(pair_type_unordered, levels = c("Glut-Glut", "E-I", "Gaba-Gaba")))
write_csv(pair_summary, file.path(path_out, "pair_summary.csv"))

agg_pair <- dedup %>%
  group_by(`cluster.1_slide`, `cluster.1_layer`, pair_type_unordered) %>%
  summarise(
    median_j = median(jaccard, na.rm = TRUE),
    mean_j = mean(jaccard, na.rm = TRUE),
    prop_cont95 = mean(max_containment >= 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(block_id = paste(`cluster.1_slide`, `cluster.1_layer`, sep = "|"))

pivot_pair <- agg_pair %>%
  select(block_id, `cluster.1_slide`, `cluster.1_layer`, pair_type_unordered, median_j) %>%
  pivot_wider(names_from = pair_type_unordered, values_from = median_j) %>%
  drop_na()

pivot_cont <- agg_pair %>%
  select(block_id, `cluster.1_slide`, `cluster.1_layer`, pair_type_unordered, prop_cont95) %>%
  pivot_wider(names_from = pair_type_unordered, values_from = prop_cont95) %>%
  drop_na()

wilcox_median_GG_vs_EI <- wilcox.test(pivot_pair$`Glut-Glut`, pivot_pair$`E-I`, paired = TRUE, exact = FALSE)
wilcox_median_II_vs_EI <- wilcox.test(pivot_pair$`Gaba-Gaba`, pivot_pair$`E-I`, paired = TRUE, exact = FALSE)
wilcox_cont_EI_vs_GG   <- wilcox.test(pivot_cont$`E-I`, pivot_cont$`Glut-Glut`, paired = TRUE, exact = FALSE)
wilcox_cont_EI_vs_II   <- wilcox.test(pivot_cont$`E-I`, pivot_cont$`Gaba-Gaba`, paired = TRUE, exact = FALSE)

# -----------------------------
# Result block 2: E-I asymmetry
# -----------------------------
ei <- dedup %>%
  filter(pair_type_unordered == "E-I") %>%
  mutate(
    gaba_sub = if_else(`cluster.1_cell_Neruon_type` == "Gaba", sub1_short, sub2_short),
    glut_sub = if_else(`cluster.1_cell_Neruon_type` == "Glut", sub1_short, sub2_short),
    gaba_overlap_percent = if_else(`cluster.1_cell_Neruon_type` == "Gaba", `cluster.1.overlap.percent`, `cluster.2.overlap.percent`),
    glut_overlap_percent = if_else(`cluster.1_cell_Neruon_type` == "Glut", `cluster.1.overlap.percent`, `cluster.2.overlap.percent`),
    delta_gaba_minus_glut_cont = gaba_overlap_percent - glut_overlap_percent,
    gaba_total_cell_num = if_else(`cluster.1_cell_Neruon_type` == "Gaba", `cluster.1_total_cell_num`, `cluster.2_total_cell_num`),
    glut_total_cell_num = if_else(`cluster.1_cell_Neruon_type` == "Glut", `cluster.1_total_cell_num`, `cluster.2_total_cell_num`),
    log_size_ratio = log(glut_total_cell_num / gaba_total_cell_num),
    max_containment = pmax(gaba_overlap_percent, glut_overlap_percent)
  )

k_pair <- sum(ei$delta_gaba_minus_glut_cont > 0, na.rm = TRUE)
n_pair <- sum(ei$delta_gaba_minus_glut_cont != 0, na.rm = TRUE)
binom_pair <- binom.test(k_pair, n_pair, p = 0.5, alternative = "greater")

ei_slide_layer <- ei %>%
  group_by(`cluster.1_slide`, `cluster.1_layer`) %>%
  summarise(
    median_delta = median(delta_gaba_minus_glut_cont, na.rm = TRUE),
    mean_delta   = mean(delta_gaba_minus_glut_cont, na.rm = TRUE),
    prop_gaba_more = mean(delta_gaba_minus_glut_cont > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(block_id = paste(`cluster.1_slide`, `cluster.1_layer`, sep = "|"))

k_sl <- sum(ei_slide_layer$median_delta > 0, na.rm = TRUE)
n_sl <- nrow(ei_slide_layer)
binom_slide_layer <- binom.test(k_sl, n_sl, p = 0.5, alternative = "greater")

cor_ei <- suppressWarnings(cor.test(
  ei$log_size_ratio,
  ei$delta_gaba_minus_glut_cont,
  method = "spearman",
  exact = FALSE
))

asym_summary <- tibble(
  metric = c(
    "E-I pair-level proportion: GABA containment > Glut containment",
    "E-I slide-layer proportion: median delta > 0",
    "Spearman rho: log(Glut/GABA size ratio) vs asymmetry delta"
  ),
  value = c(k_pair / n_pair, k_sl / n_sl, unname(cor_ei$estimate)),
  p = c(binom_pair$p.value, binom_slide_layer$p.value, cor_ei$p.value),
  n = c(n_pair, n_sl, sum(complete.cases(ei$log_size_ratio, ei$delta_gaba_minus_glut_cont)))
)
write_csv(asym_summary, file.path(path_out, "asymmetry_summary.csv"))
write_csv(ei_slide_layer, file.path(path_out, "ei_slide_layer_summary.csv"))

# -----------------------------
# Result block 3: canonical GABA-GABA motifs
# -----------------------------
gaba_pairs <- dedup %>%
  filter(pair_type_unordered == "Gaba-Gaba") %>%
  mutate(
    subA = if_else(sub1_short < sub2_short, sub1_short, sub2_short),
    subB = if_else(sub1_short < sub2_short, sub2_short, sub1_short),
    pair_subclass = paste(subA, subB, sep = "||"),
    slide_layer = paste(`cluster.1_slide`, `cluster.1_layer`, sep = "|")
  )

canonical_subs <- c("Vip","Lamp5","Sncg","Sst","Pvalb","Pvalb chandelier","Sst Chodl")

canonical_gaba_pair_summary <- gaba_pairs %>%
  filter(subA %in% canonical_subs, subB %in% canonical_subs) %>%
  group_by(subA, subB) %>%
  summarise(
    n = n(),
    mean_j = mean(jaccard, na.rm = TRUE),
    median_j = median(jaccard, na.rm = TRUE),
    p75 = quantile(jaccard, 0.75, na.rm = TRUE),
    p90 = quantile(jaccard, 0.90, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_j), desc(n))
write_csv(canonical_gaba_pair_summary, file.path(path_out, "canonical_gaba_pair_summary.csv"))

canonical_interest <- c("Sst||Vip", "Pvalb||Vip", "Lamp5||Vip", "Pvalb||Sst", "Lamp5||Sst")

canonical_piv <- gaba_pairs %>%
  filter(pair_subclass %in% canonical_interest) %>%
  group_by(slide_layer, pair_subclass) %>%
  summarise(median_j = median(jaccard, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = pair_subclass, values_from = median_j)

canonical_tests <- bind_rows(
  paired_compare(canonical_piv, "Sst||Vip", "Pvalb||Vip") %>% mutate(comparison = "Vip-Sst vs Vip-Pvalb"),
  paired_compare(canonical_piv, "Pvalb||Sst", "Pvalb||Vip") %>% mutate(comparison = "Pvalb-Sst vs Pvalb-Vip"),
  paired_compare(canonical_piv, "Lamp5||Sst", "Lamp5||Vip") %>% mutate(comparison = "Lamp5-Sst vs Lamp5-Vip")
) %>%
  select(comparison, everything())
write_csv(canonical_tests, file.path(path_out, "canonical_pair_tests_selected.csv"))

# -----------------------------
# Result block 4: layer dependence
# -----------------------------
ei_layer_asymmetry <- ei %>%
  group_by(`cluster.1_layer`) %>%
  summarise(
    n = n(),
    prop_gaba_more = mean(delta_gaba_minus_glut_cont > 0, na.rm = TRUE),
    median_delta   = median(delta_gaba_minus_glut_cont, na.rm = TRUE),
    mean_delta     = mean(delta_gaba_minus_glut_cont, na.rm = TRUE),
    prop_cont95    = mean(max_containment >= 0.95, na.rm = TRUE),
    median_j       = median(jaccard, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(ei_layer_asymmetry, file.path(path_out, "ei_layer_asymmetry.csv"))

layer_tests <- tibble(
  pair_type = c("E-I", "E-I", "Glut-Glut", "Glut-Glut", "Gaba-Gaba", "Gaba-Gaba"),
  metric = c("median_j", "prop_cont95", "median_j", "prop_cont95", "median_j", "prop_cont95"),
  H = c(
    kruskal.test(median_j ~ `cluster.1_layer`, data = filter(agg_pair, pair_type_unordered == "E-I"))$statistic,
    kruskal.test(prop_cont95 ~ `cluster.1_layer`, data = filter(agg_pair, pair_type_unordered == "E-I"))$statistic,
    kruskal.test(median_j ~ `cluster.1_layer`, data = filter(agg_pair, pair_type_unordered == "Glut-Glut"))$statistic,
    kruskal.test(prop_cont95 ~ `cluster.1_layer`, data = filter(agg_pair, pair_type_unordered == "Glut-Glut"))$statistic,
    kruskal.test(median_j ~ `cluster.1_layer`, data = filter(agg_pair, pair_type_unordered == "Gaba-Gaba"))$statistic,
    kruskal.test(prop_cont95 ~ `cluster.1_layer`, data = filter(agg_pair, pair_type_unordered == "Gaba-Gaba"))$statistic
  ),
  p = c(
    kruskal.test(median_j ~ `cluster.1_layer`, data = filter(agg_pair, pair_type_unordered == "E-I"))$p.value,
    kruskal.test(prop_cont95 ~ `cluster.1_layer`, data = filter(agg_pair, pair_type_unordered == "E-I"))$p.value,
    kruskal.test(median_j ~ `cluster.1_layer`, data = filter(agg_pair, pair_type_unordered == "Glut-Glut"))$p.value,
    kruskal.test(prop_cont95 ~ `cluster.1_layer`, data = filter(agg_pair, pair_type_unordered == "Glut-Glut"))$p.value,
    kruskal.test(median_j ~ `cluster.1_layer`, data = filter(agg_pair, pair_type_unordered == "Gaba-Gaba"))$p.value,
    kruskal.test(prop_cont95 ~ `cluster.1_layer`, data = filter(agg_pair, pair_type_unordered == "Gaba-Gaba"))$p.value
  )
)
write_csv(layer_tests, file.path(path_out, "layer_kruskal_tests.csv"))

# -----------------------------
# key statistics
# -----------------------------
key_stats <- bind_rows(
  tibble(test = "paired Wilcoxon: slide-layer median Jaccard, Glut-Glut vs E-I",
         statistic = unname(wilcox_median_GG_vs_EI$statistic), p = wilcox_median_GG_vs_EI$p.value),
  tibble(test = "paired Wilcoxon: slide-layer median Jaccard, Gaba-Gaba vs E-I",
         statistic = unname(wilcox_median_II_vs_EI$statistic), p = wilcox_median_II_vs_EI$p.value),
  tibble(test = "paired Wilcoxon: slide-layer prop(max overlap>=0.95), E-I vs Glut-Glut",
         statistic = unname(wilcox_cont_EI_vs_GG$statistic), p = wilcox_cont_EI_vs_GG$p.value),
  tibble(test = "paired Wilcoxon: slide-layer prop(max overlap>=0.95), E-I vs Gaba-Gaba",
         statistic = unname(wilcox_cont_EI_vs_II$statistic), p = wilcox_cont_EI_vs_II$p.value),
  tibble(test = "Binomial: in E-I pairs, GABA is more contained than Glut (pair-level)",
         statistic = k_pair / n_pair, p = binom_pair$p.value),
  tibble(test = "Binomial: in E-I pairs, GABA is more contained than Glut (slide-layer median)",
         statistic = k_sl / n_sl, p = binom_slide_layer$p.value),
  tibble(test = "Spearman: log(Glut/GABA size ratio) vs delta containment",
         statistic = unname(cor_ei$estimate), p = cor_ei$p.value)
)
write_csv(key_stats, file.path(path_out, "key_statistics.csv"))

# ============================================================
# PLOT 1: 顶刊风格双 donut 图（彻底修正版）
# 不再用 coord_polar，改为手工绘制，保证标签位置稳定
# ============================================================

save_donut_pub <- function(filename, width = 12, height = 5.8, expr_draw) {
  pdf(file.path(path_out, paste0(filename, ".pdf")),
      width = width, height = height, family = "sans", useDingbats = FALSE)
  expr_draw()
  dev.off()
  
  png(file.path(path_out, paste0(filename, ".png")),
      width = width, height = height, units = "in", res = 600, bg = "white")
  expr_draw()
  dev.off()
}

draw_single_donut <- function(cx, cy, counts, groups, fills, center_title,
                              outer_r = 1.55, inner_r = 0.82,
                              line_r1 = 1.58, line_r2 = 1.95,
                              text_dx = 0.24,
                              cex_center = 1.55,
                              cex_label = 1.02) {
  total_n <- sum(counts)
  props <- counts / total_n
  cumprops <- c(0, cumsum(props))
  
  # 从顶部开始，顺时针
  theta_start <- pi/2 - 2*pi*cumprops[-length(cumprops)]
  theta_end   <- pi/2 - 2*pi*cumprops[-1]
  
  # 先画 donut
  for (i in seq_along(counts)) {
    th <- seq(theta_start[i], theta_end[i], length.out = 300)
    x_outer <- cx + outer_r * cos(th)
    y_outer <- cy + outer_r * sin(th)
    x_inner <- cx + inner_r * cos(rev(th))
    y_inner <- cy + inner_r * sin(rev(th))
    
    polygon(
      c(x_outer, x_inner),
      c(y_outer, y_inner),
      col = fills[i],
      border = "white",
      lwd = 2.2
    )
  }
  
  # 中心文字
  text(
    cx, cy,
    labels = paste0(center_title, "\nN = ", format(total_n, big.mark = ",")),
    cex = cex_center, font = 2, linespacing = 1.1
  )
  
  # 再画引导线和标签
  for (i in seq_along(counts)) {
    mid <- (theta_start[i] + theta_end[i]) / 2
    side <- ifelse(cos(mid) >= 0, 1, -1)
    
    # 线段起点、转折点、终点
    x0 <- cx + line_r1 * cos(mid)
    y0 <- cy + line_r1 * sin(mid)
    x1 <- cx + line_r2 * cos(mid)
    y1 <- cy + line_r2 * sin(mid)
    x2 <- x1 + side * 0.34
    y2 <- y1
    
    segments(x0, y0, x1, y1, col = "grey35", lwd = 1.3)
    segments(x1, y1, x2, y2, col = "grey35", lwd = 1.3)
    
    lab <- paste0(
      groups[i], "\n",
      "n = ", format(counts[i], big.mark = ","), " / ", format(total_n, big.mark = ","), "\n",
      sprintf("%.1f%%", props[i] * 100)
    )
    
    text(
      x = x2 + side * text_dx,
      y = y2,
      labels = lab,
      adj = c(ifelse(side > 0, 0, 1), 0.5),
      cex = cex_label,
      font = 2,
      linespacing = 1.05
    )
  }
}

save_donut_pub(
  filename = "01_EI_asymmetry_pair_vs_block",
  width = 12,
  height = 5.8,
  expr_draw = function() {
    par(mar = c(0.3, 0.3, 0.3, 0.3),
        oma = c(0.4, 0.4, 2.7, 0.4),
        xpd = NA)
    
    plot.new()
    plot.window(xlim = c(-6.6, 6.6), ylim = c(-2.5, 2.7), asp = 1)
    
    # 主标题
    mtext(
      "E-I spatial organization is strongly asymmetric",
      side = 3, outer = TRUE, line = 0.6, cex = 1.85, font = 2
    )
    
    # 副标题
    mtext(
      paste0(
        "Pair-level binomial p = ", fmt_p(binom_pair$p.value),
        "   |   Slide-layer median binomial p = ", fmt_p(binom_slide_layer$p.value)
      ),
      side = 3, outer = TRUE, line = -0.9, cex = 1.15, font = 1
    )
    
    # 左图：pair-level
    draw_single_donut(
      cx = -3.2, cy = -0.1,
      counts = c(k_pair, n_pair - k_pair),
      groups = c("GABA more contained", "Glut more contained / tie"),
      fills = c("#6A3D9A", "#D9D9D9"),
      center_title = "Unique E-I pairs",
      outer_r = 1.55,
      inner_r = 0.82,
      line_r1 = 1.60,
      line_r2 = 1.95,
      text_dx = 0.20,
      cex_center = 1.45,
      cex_label = 0.98
    )
    
    # 右图：slide-layer
    draw_single_donut(
      cx = 3.2, cy = -0.1,
      counts = c(k_sl, n_sl - k_sl),
      groups = c("GABA more contained", "Glut more contained / tie"),
      fills = c("#6A3D9A", "#D9D9D9"),
      center_title = "Slide-layer units",
      outer_r = 1.55,
      inner_r = 0.82,
      line_r1 = 1.60,
      line_r2 = 1.95,
      text_dx = 0.20,
      cex_center = 1.45,
      cex_label = 0.98
    )
  }
)

# ============================================================
# PLOT 2: size ratio vs asymmetry
# ============================================================
p2 <- ggplot(ei, aes(x = log_size_ratio, y = delta_gaba_minus_glut_cont)) +
  geom_hline(yintercept = 0, color = "grey75", linewidth = 0.5, linetype = "dashed") +
  geom_point(color = "#3B528B", alpha = 0.18, size = 1.0) +
  geom_smooth(method = "loess", se = TRUE, color = "#F03B20", fill = "#FDD0A2", linewidth = 1.0) +
  annotate(
    "text",
    x = Inf, y = Inf,
    hjust = 1.05, vjust = 1.5,
    label = paste0("Spearman rho = ", sprintf("%.3f", unname(cor_ei$estimate)),
                   "\np = ", fmt_p(cor_ei$p.value)),
    size = 4.2, fontface = "bold"
  ) +
  labs(
    title = "Asymmetry scales with excitatory-to-inhibitory size difference",
    subtitle = "Positive values indicate that the GABA cluster is more contained than the Glut cluster",
    x = "log(Glut cluster size / GABA cluster size)",
    y = "Containment asymmetry (GABA - Glut)"
  ) +
  theme_pub(12)
save_pub(p2, "02_EI_size_ratio_vs_asymmetry", width = 7.0, height = 5.4)

# ============================================================
# PLOT 3: slide-layer median Jaccard comparison
# ============================================================
plot_j <- pivot_pair %>%
  pivot_longer(cols = c(`Glut-Glut`, `E-I`, `Gaba-Gaba`),
               names_to = "pair_type", values_to = "value") %>%
  mutate(pair_type = factor(pair_type, levels = c("Glut-Glut", "E-I", "Gaba-Gaba")))

p3 <- ggplot(plot_j, aes(x = pair_type, y = value, color = pair_type)) +
  geom_line(aes(group = block_id), color = "grey82", alpha = 0.35, linewidth = 0.4) +
  geom_boxplot(aes(fill = pair_type), width = 0.50, alpha = 0.80, outlier.shape = NA, color = "black") +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3.2, fill = "white", color = "black") +
  scale_fill_manual(values = pair_cols) +
  scale_color_manual(values = pair_cols) +
  labs(
    title = "E-I pairs are not the most strongly overlapped under a symmetric metric",
    subtitle = paste0(
      "Paired Wilcoxon: Glut-Glut vs E-I p = ", fmt_p(wilcox_median_GG_vs_EI$p.value),
      "   |   Gaba-Gaba vs E-I p = ", fmt_p(wilcox_median_II_vs_EI$p.value)
    ),
    x = NULL,
    y = "Median Jaccard within each slide-layer block"
  ) +
  theme_pub(12) +
  theme(legend.position = "none")
save_pub(p3, "03_Jaccard_block_comparison", width = 7.6, height = 5.2)

# ============================================================
# PLOT 4: slide-layer near-complete containment comparison
# ============================================================
plot_c <- pivot_cont %>%
  pivot_longer(cols = c(`Glut-Glut`, `E-I`, `Gaba-Gaba`),
               names_to = "pair_type", values_to = "value") %>%
  mutate(pair_type = factor(pair_type, levels = c("Glut-Glut", "E-I", "Gaba-Gaba")))

p4 <- ggplot(plot_c, aes(x = pair_type, y = value, color = pair_type)) +
  geom_line(aes(group = block_id), color = "grey82", alpha = 0.35, linewidth = 0.4) +
  geom_boxplot(aes(fill = pair_type), width = 0.50, alpha = 0.80, outlier.shape = NA, color = "black") +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3.2, fill = "white", color = "black") +
  scale_fill_manual(values = pair_cols) +
  scale_color_manual(values = pair_cols) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "E-I pairs show the strongest near-complete containment",
    subtitle = paste0(
      "Paired Wilcoxon: E-I vs Glut-Glut p = ", fmt_p(wilcox_cont_EI_vs_GG$p.value),
      "   |   E-I vs Gaba-Gaba p = ", fmt_p(wilcox_cont_EI_vs_II$p.value)
    ),
    x = NULL,
    y = "Proportion with max containment >= 0.95"
  ) +
  theme_pub(12) +
  theme(legend.position = "none")
save_pub(p4, "04_Containment_block_comparison", width = 7.6, height = 5.2)

# ============================================================
# PLOT 5: canonical GABA-GABA heatmap
# ============================================================
canonical_order <- c("Vip","Lamp5","Sncg","Sst","Pvalb","Pvalb chandelier")

heat_df <- canonical_gaba_pair_summary %>%
  filter(subA %in% canonical_order, subB %in% canonical_order) %>%
  group_by(subA, subB) %>%
  summarise(median_j = median(median_j, na.rm = TRUE), .groups = "drop")

heat_df2 <- bind_rows(
  heat_df,
  heat_df %>% transmute(subA = subB, subB = subA, median_j = median_j)
) %>%
  group_by(subA, subB) %>%
  summarise(median_j = median(median_j, na.rm = TRUE), .groups = "drop") %>%
  complete(subA = canonical_order, subB = canonical_order, fill = list(median_j = NA_real_)) %>%
  mutate(
    subA = factor(subA, levels = canonical_order),
    subB = factor(subB, levels = rev(canonical_order)))

p5 <- ggplot(heat_df2, aes(x = subA, y = subB, fill = median_j)) +
  geom_tile(color = "white", linewidth = 0.55) +
  geom_text(
    data = heat_df2 %>% filter(!is.na(median_j)),
    aes(label = sprintf("%.2f", median_j)),
    size = 4.2, fontface = "bold", color = "black"
  ) +
  scale_fill_gradient(low = "#F7FBFF", high = "#1D91C0", na.value = "grey96") +
  labs(
    title = "Canonical inhibitory motifs do not show sharply selective spatial co-residence",
    subtitle = paste0(
      "Vip-Sst vs Vip-Pvalb p = ", fmt_p(canonical_tests$p[canonical_tests$comparison == "Vip-Sst vs Vip-Pvalb"]),
      "   |   Pvalb-Sst vs Pvalb-Vip p = ", fmt_p(canonical_tests$p[canonical_tests$comparison == "Pvalb-Sst vs Pvalb-Vip"]),
      "\nLamp5-Sst vs Lamp5-Vip p = ", fmt_p(canonical_tests$p[canonical_tests$comparison == "Lamp5-Sst vs Lamp5-Vip"])
    ),
    x = NULL, y = NULL, fill = "Median Jaccard"
  ) +
  theme_pub(11.5) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.text.y = element_text(angle = 0)
  )
save_pub(p5, "05_Canonical_GABA_heatmap", width = 8.0, height = 6.4)

# ============================================================
# PLOT 6: selected canonical GABA pair comparisons
# ============================================================
selected_plot_df <- bind_rows(
  canonical_piv %>% transmute(
    slide_layer,
    comparison = "Vip-centered comparison",
    pair_label = "Vip-Sst",
    median_j = `Sst||Vip`
  ),
  canonical_piv %>% transmute(
    slide_layer,
    comparison = "Vip-centered comparison",
    pair_label = "Vip-Pvalb",
    median_j = `Pvalb||Vip`
  ),
  canonical_piv %>% transmute(
    slide_layer,
    comparison = "Pvalb-centered comparison",
    pair_label = "Pvalb-Sst",
    median_j = `Pvalb||Sst`
  ),
  canonical_piv %>% transmute(
    slide_layer,
    comparison = "Pvalb-centered comparison",
    pair_label = "Pvalb-Vip",
    median_j = `Pvalb||Vip`
  ),
  canonical_piv %>% transmute(
    slide_layer,
    comparison = "Lamp5-centered comparison",
    pair_label = "Lamp5-Sst",
    median_j = `Lamp5||Sst`
  ),
  canonical_piv %>% transmute(
    slide_layer,
    comparison = "Lamp5-centered comparison",
    pair_label = "Lamp5-Vip",
    median_j = `Lamp5||Vip`
  )
) %>%
  filter(!is.na(median_j)) %>%
  mutate(
    comparison = factor(comparison,
                        levels = c("Vip-centered comparison",
                                   "Pvalb-centered comparison",
                                   "Lamp5-centered comparison")),
    pair_label = factor(pair_label,
                        levels = c("Vip-Sst","Vip-Pvalb","Pvalb-Sst","Pvalb-Vip","Lamp5-Sst","Lamp5-Vip"))
  )

ann_df <- tibble(
  comparison = factor(c("Vip-centered comparison",
                        "Pvalb-centered comparison",
                        "Lamp5-centered comparison"),
                      levels = c("Vip-centered comparison",
                                 "Pvalb-centered comparison",
                                 "Lamp5-centered comparison")),
  x = 1.5,
  y = 0.98,
  label = c(
    paste0("p = ", fmt_p(canonical_tests$p[canonical_tests$comparison == "Vip-Sst vs Vip-Pvalb"])),
    paste0("p = ", fmt_p(canonical_tests$p[canonical_tests$comparison == "Pvalb-Sst vs Pvalb-Vip"])),
    paste0("p = ", fmt_p(canonical_tests$p[canonical_tests$comparison == "Lamp5-Sst vs Lamp5-Vip"]))
  )
)

pair_fill <- c(
  "Vip-Sst" = "#3BAF92",
  "Vip-Pvalb" = "#E67E22",
  "Pvalb-Sst" = "#E84393",
  "Pvalb-Vip" = "#F4A261",
  "Lamp5-Sst" = "#7CB342",
  "Lamp5-Vip" = "#8DA0CB"
)

p6 <- ggplot(selected_plot_df, aes(x = pair_label, y = median_j, fill = pair_label)) +
  geom_boxplot(width = 0.62, alpha = 0.88, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.10, size = 1.0, alpha = 0.18, color = "black") +
  geom_text(
    data = ann_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 4.1
  ) +
  facet_wrap(~comparison, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = pair_fill) +
  coord_cartesian(ylim = c(0, 1.02), clip = "off") +
  labs(
    title = "Canonical inhibitory pairs show broad overlap rather than strong spatial exclusivity",
    subtitle = "Each panel corresponds exactly to one slide-layer paired comparison reported in the statistics",
    x = NULL,
    y = "Median Jaccard within each slide-layer block"
  ) +
  theme_pub(12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 1)
  )
save_pub(p6, "06_Canonical_GABA_selected_pairs", width = 11.2, height = 5.4)

# ============================================================
# PLOT 7: E-I layer dependence in median Jaccard
# ============================================================
ei_block_j <- agg_pair %>%
  filter(pair_type_unordered == "E-I", `cluster.1_layer` != "layer 6b") %>%
  mutate(`cluster.1_layer` = factor(`cluster.1_layer`,
                                    levels = c("layer 1","layer 2/3","layer 4","layer 5","layer 6a")))

kw_ei_j <- kruskal.test(median_j ~ `cluster.1_layer`, data = ei_block_j)

p7 <- ggplot(ei_block_j, aes(x = `cluster.1_layer`, y = median_j)) +
  geom_boxplot(width = 0.62, fill = "#B2DF8A", color = "black", alpha = 0.85, outlier.shape = NA) +
  geom_jitter(width = 0.14, size = 1.0, alpha = 0.25, color = "#4D4D4D") +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3.0, fill = "#1B7837", color = "black") +
  labs(
    title = "E-I overlap strength shows strong layer dependence",
    subtitle = paste0("Kruskal-Wallis p = ", fmt_p(kw_ei_j$p.value)),
    x = NULL, y = "Median Jaccard within each slide-layer block"
  ) +
  theme_pub(12)
save_pub(p7, "07_EI_layer_median_jaccard", width = 7.0, height = 5.0)

# ============================================================
# PLOT 8: E-I layer dependence in near-complete containment
# ============================================================
ei_block_c <- agg_pair %>%
  filter(pair_type_unordered == "E-I", `cluster.1_layer` != "layer 6b") %>%
  mutate(`cluster.1_layer` = factor(`cluster.1_layer`,
                                    levels = c("layer 1","layer 2/3","layer 4","layer 5","layer 6a")))

kw_ei_c <- kruskal.test(prop_cont95 ~ `cluster.1_layer`, data = ei_block_c)

layer_prop95_summary <- ei_block_c %>%
  group_by(`cluster.1_layer`) %>%
  summarise(
    mean_prop = mean(prop_cont95, na.rm = TRUE),
    median_prop = median(prop_cont95, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(layer_prop95_summary, file.path(path_out, "ei_layer_prop95_summary.csv"))

p8 <- ggplot(ei_block_c, aes(x = `cluster.1_layer`, y = prop_cont95)) +
  geom_boxplot(width = 0.62, fill = "#A6CEE3", color = "black", alpha = 0.85, outlier.shape = NA) +
  geom_jitter(width = 0.14, size = 1.0, alpha = 0.25, color = "#4D4D4D") +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3.0, fill = "#1F78B4", color = "black") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Near-complete E-I containment is strongest in superficial-middle layers",
    subtitle = paste0(
      "Kruskal-Wallis p = ", fmt_p(kw_ei_c$p.value),
      "   |   Highest block-level containment is expected in layer 2/3"
    ),
    x = NULL, y = "Proportion with max containment >= 0.95"
  ) +
  theme_pub(12)
save_pub(p8, "08_EI_layer_containment95", width = 7.0, height = 5.0)

# -----------------------------
# session info
# -----------------------------
sink(file.path(path_out, "sessionInfo.txt"))
print(sessionInfo())
sink()

cat("All analyses finished.\n")
cat("Output directory:\n", normalizePath(path_out), "\n")