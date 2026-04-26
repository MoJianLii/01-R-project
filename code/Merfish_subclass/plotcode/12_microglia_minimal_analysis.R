# 清空环境
rm(list = ls())
gc()
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
  library(scales)
  library(grid)
})

# avoid select/filter conflicts
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
arrange <- dplyr::arrange
distinct <- dplyr::distinct
pull <- dplyr::pull
count <- dplyr::count
group_by <- dplyr::group_by
summarise <- dplyr::summarise
transmute <- dplyr::transmute
ungroup <- dplyr::ungroup
rename <- dplyr::rename
left_join <- dplyr::left_join

# ============================================================
# 0. Paths
# ============================================================
project_root <- "E:/zaw/2603"
out_dir <- file.path(project_root, "R_output_12")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

find_first_existing <- function(paths) {
  hit <- paths[file.exists(paths)][1]
  if (is.na(hit)) stop("No valid input file found. Please check input path.")
  hit
}

input_neuron_nonneuron <- find_first_existing(c(
  file.path(project_root, "neuron-nonneuron-partner.csv"),
  file.path(project_root, "data", "neuron-nonneuron-partner.csv"),
  file.path(getwd(), "neuron-nonneuron-partner.csv")
))

# ============================================================
# 1. Helpers
# ============================================================
col_main <- c(
  "Micro->Gaba exact" = "#C44E52",
  "Micro->Glut exact" = "#4C78A8"
)
col_dir <- c("Enriched" = "#4C78A8", "Depleted" = "#E15759")

fmt_p <- function(p) {
  vapply(p, function(x) {
    if (is.na(x)) return(NA_character_)
    if (x < 1e-300) return("<1e-300")
    if (x < 1e-3) return(formatC(x, format = "e", digits = 2))
    sprintf("%.3f", x)
  }, character(1))
}
fmt_num <- function(x, digits = 3) sprintf(paste0("%.", digits, "f"), x)
wrap_title <- function(x) stringr::str_wrap(x, width = 62)
wrap_subtitle <- function(x) stringr::str_wrap(x, width = 86)

base_theme <- function() {
  theme_classic(base_size = 14) +
    theme(
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.title = element_text(size = 21, face = "bold", colour = "black", hjust = 0, lineheight = 1.02,
                                margin = margin(b = 8)),
      plot.subtitle = element_text(size = 12.8, colour = "#444444", hjust = 0, lineheight = 1.18,
                                   margin = margin(b = 12)),
      axis.title = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 12.5, colour = "black"),
      axis.line = element_line(linewidth = 0.8, colour = "black"),
      axis.ticks = element_line(linewidth = 0.8, colour = "black"),
      axis.ticks.length = unit(0.18, "cm"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      strip.background = element_rect(fill = "#F1F1F1", colour = NA),
      strip.text = element_text(size = 13.8, face = "bold"),
      panel.grid = element_blank(),
      plot.margin = margin(16, 36, 14, 26)
    )
}

save_pub <- function(plot, filename, width, height) {
  ggsave(file.path(fig_dir, filename), plot = plot, width = width, height = height,
         dpi = 600, bg = "white", limitsize = FALSE)
}

standardize_names <- function(df) {
  names(df) <- names(df) |>
    stringr::str_replace_all("\\.overlap\\.percent", "_overlap_percent") |>
    stringr::str_replace_all("\\.", "_")
  df
}

clean_subclass <- function(x) {
  x |>
    as.character() |>
    stringr::str_replace("^\\s*\\d+\\s*", "") |>
    stringr::str_replace("\\s+(Glut|Gaba|NN)\\s*$", "") |>
    stringr::str_trim()
}

pair_key <- function(a, b) ifelse(a <= b, paste0(a, "|||", b), paste0(b, "|||", a))
safe_logp <- function(p) -log10(pmax(p, 1e-300))
safe_div <- function(a, b) ifelse(b > 0, a / b, NA_real_)

fisher_tbl <- function(a, b, c, d) {
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))
  tibble(
    a = a, b = b, c = c, d = d,
    or = unname(ft$estimate),
    ci_low = unname(ft$conf.int[1]),
    ci_high = unname(ft$conf.int[2]),
    p = ft$p.value
  )
}

mw_compare <- function(df, group_col, focal, vars) {
  purrr::map_dfr(vars, function(v) {
    x <- df |> dplyr::filter(.data[[group_col]] == focal) |> dplyr::pull(.data[[v]])
    y <- df |> dplyr::filter(.data[[group_col]] != focal) |> dplyr::pull(.data[[v]])
    wt <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
    tibble(
      var = v,
      focal_median = median(x, na.rm = TRUE),
      other_median = median(y, na.rm = TRUE),
      delta = median(x, na.rm = TRUE) - median(y, na.rm = TRUE),
      p = wt$p.value,
      n_focal = sum(!is.na(x)),
      n_other = sum(!is.na(y))
    )
  }) |>
    mutate(fdr = p.adjust(p, method = "BH"))
}

topology_label <- function(non_overlap, neuron_overlap, tol = 1e-9) {
  dplyr::case_when(
    non_overlap >= 1 - tol & neuron_overlap < 1 - tol ~ "Non->Neuron exact",
    neuron_overlap >= 1 - tol & non_overlap < 1 - tol ~ "Neuron->Non exact",
    neuron_overlap >= 1 - tol & non_overlap >= 1 - tol ~ "Bidirectional exact",
    non_overlap >= 0.95 & non_overlap < 1 - tol & neuron_overlap < 0.95 ~ "Non->Neuron near-full",
    neuron_overlap >= 0.95 & neuron_overlap < 1 - tol & non_overlap < 0.95 ~ "Neuron->Non near-full",
    pmax(non_overlap, neuron_overlap) >= 0.95 ~ "High-overlap other",
    TRUE ~ "Other"
  )
}

# ============================================================
# 2. Data
# ============================================================
raw_nn_non <- readr::read_csv(input_neuron_nonneuron, show_col_types = FALSE) |>
  standardize_names()

dedup_nn_non <- raw_nn_non |>
  mutate(pair_key = pair_key(cluster_1_label, cluster_2_label)) |>
  distinct(pair_key, .keep_all = TRUE)

nn <- dedup_nn_non |>
  filter(
    (cluster_1_cell_Neruon_type == "NonNeuron" & cluster_2_cell_Neruon_type %in% c("Glut", "Gaba")) |
      (cluster_2_cell_Neruon_type == "NonNeuron" & cluster_1_cell_Neruon_type %in% c("Glut", "Gaba"))
  ) |>
  mutate(non_is_1 = cluster_1_cell_Neruon_type == "NonNeuron") |>
  transmute(
    pair_key,
    slide = if_else(non_is_1, cluster_2_slide, cluster_1_slide),
    layer = if_else(non_is_1, cluster_2_layer, cluster_1_layer),
    layer_simplified = stringr::str_remove(if_else(non_is_1, cluster_2_layer, cluster_1_layer), "^layer\\s+"),
    region = if_else(non_is_1, cluster_2_region, cluster_1_region),
    
    neuron_label = if_else(non_is_1, cluster_2_label, cluster_1_label),
    neuron_type = if_else(non_is_1, cluster_2_cell_Neruon_type, cluster_1_cell_Neruon_type),
    neuron_subclass = clean_subclass(if_else(non_is_1, cluster_2_subclass, cluster_1_subclass)),
    neuron_total = if_else(non_is_1, cluster_2_total_cell_num, cluster_1_total_cell_num),
    neuron_glut_n = if_else(non_is_1, cluster_2_Glut_Neruon_cell_ids_num, cluster_1_Glut_Neruon_cell_ids_num),
    neuron_gaba_n = if_else(non_is_1, cluster_2_GABA_Neruon_cell_ids_num, cluster_1_GABA_Neruon_cell_ids_num),
    neuron_p = if_else(non_is_1, cluster_2_cauchy_combination_p, cluster_1_cauchy_combination_p),
    neuron_overlap_pct = if_else(non_is_1, cluster_2_overlap_percent, cluster_1_overlap_percent),
    
    non_label = if_else(non_is_1, cluster_1_label, cluster_2_label),
    non_class = clean_subclass(if_else(non_is_1, cluster_1_subclass, cluster_2_subclass)),
    non_total = if_else(non_is_1, cluster_1_total_cell_num, cluster_2_total_cell_num),
    non_glut_n = if_else(non_is_1, cluster_1_Glut_Neruon_cell_ids_num, cluster_2_Glut_Neruon_cell_ids_num),
    non_gaba_n = if_else(non_is_1, cluster_1_GABA_Neruon_cell_ids_num, cluster_2_GABA_Neruon_cell_ids_num),
    non_p = if_else(non_is_1, cluster_1_cauchy_combination_p, cluster_2_cauchy_combination_p),
    non_overlap_pct = if_else(non_is_1, cluster_1_overlap_percent, cluster_2_overlap_percent),
    
    overlap_cell,
    union_cell,
    jaccard
  ) |>
  mutate(
    topology = topology_label(non_overlap_pct, neuron_overlap_pct),
    log10_neuron_non_ratio = log10((neuron_total + 1) / (non_total + 1)),
    containment_asym = non_overlap_pct - neuron_overlap_pct,
    non_glut_fraction = safe_div(non_glut_n, non_total),
    non_gaba_fraction = safe_div(non_gaba_n, non_total),
    non_non_fraction = 1 - safe_div(non_glut_n + non_gaba_n, non_total),
    neuron_logp = safe_logp(neuron_p),
    non_logp = safe_logp(non_p)
  )

# ============================================================
# 3. Figure 1 | geometry signature
# ============================================================
core_vars <- c("jaccard", "log10_neuron_non_ratio", "non_total")
metric_label_map <- c(
  jaccard = "Jaccard",
  log10_neuron_non_ratio = "log10(host/child size ratio)",
  non_total = "Microglia child size"
)

sig_glut_child <- nn |>
  filter(neuron_type == "Glut", topology == "Non->Neuron exact") |>
  mw_compare(group_col = "non_class", focal = "Microglia", vars = core_vars) |>
  mutate(context = "Micro->Glut exact")

sig_gaba_child <- nn |>
  filter(neuron_type == "Gaba", topology == "Non->Neuron exact") |>
  mw_compare(group_col = "non_class", focal = "Microglia", vars = core_vars) |>
  mutate(context = "Micro->Gaba exact")

geom_compare_summary <- bind_rows(sig_gaba_child, sig_glut_child) |>
  mutate(
    var_label = recode(var, !!!metric_label_map),
    focal_median_label = if_else(var == "non_total", as.character(round(focal_median)), fmt_num(focal_median, 3)),
    other_median_label = if_else(var == "non_total", as.character(round(other_median)), fmt_num(other_median, 3))
  ) |>
  dplyr::select(context, var, var_label, focal_median, other_median, delta, p, fdr, n_focal, n_other,
                focal_median_label, other_median_label)
readr::write_csv(geom_compare_summary, file.path(tab_dir, "table01_micro_exact_geometry_summary.csv"))

plot_df_a <- geom_compare_summary |>
  mutate(
    context = factor(context, levels = c("Micro->Gaba exact", "Micro->Glut exact")),
    var_label = factor(var_label, levels = c("Microglia child size", "log10(host/child size ratio)", "Jaccard")),
    effect_direction = case_when(
      context == "Micro->Gaba exact" & var == "jaccard" ~ "Higher in Micro",
      context == "Micro->Gaba exact" & var == "log10_neuron_non_ratio" ~ "Lower in Micro",
      context == "Micro->Gaba exact" & var == "non_total" ~ "Higher in Micro",
      context == "Micro->Glut exact" & fdr < 0.05 & delta > 0 ~ "Higher in Micro",
      context == "Micro->Glut exact" & fdr < 0.05 & delta < 0 ~ "Lower in Micro",
      TRUE ~ "No clear shift"
    )
  ) |>
  group_by(context) |>
  mutate(
    xmax_panel = max(c(focal_median, other_median), na.rm = TRUE),
    xmin_panel = min(c(focal_median, other_median), na.rm = TRUE),
    panel_range = pmax(xmax_panel - xmin_panel, xmax_panel * 0.35, 0.35),
    x_lower = pmin(0, xmin_panel - panel_range * 0.28),
    ann_x = xmax_panel + panel_range * 0.55,
    x_upper = ann_x + panel_range * 0.18,
    focal_lab_x = if_else(focal_median <= other_median,
                          pmax(focal_median + panel_range * 0.04, x_lower + panel_range * 0.02),
                          focal_median),
    other_lab_x = if_else(other_median <= focal_median,
                          pmax(other_median + panel_range * 0.04, x_lower + panel_range * 0.02),
                          other_median),
    fdr_lab = if_else(context == "Micro->Glut exact" & fdr >= 0.05 & var == "non_total", ">0.05", fmt_p(fdr))
  ) |>
  ungroup()

fig1 <- ggplot(plot_df_a, aes(y = var_label)) +
  geom_blank(aes(x = x_lower)) +
  geom_blank(aes(x = x_upper)) +
  geom_segment(aes(x = other_median, xend = focal_median, yend = var_label),
               linewidth = 1.25, colour = "#B8B8B8") +
  geom_point(aes(x = other_median), size = 4.2, colour = "#8F8F8F") +
  geom_point(aes(x = focal_median, colour = context), size = 4.4, show.legend = FALSE) +
  geom_text(aes(x = focal_lab_x, label = focal_median_label, colour = context),
            nudge_y = 0.22, hjust = 0, size = 3.9, show.legend = FALSE) +
  geom_text(aes(x = other_lab_x, label = other_median_label),
            nudge_y = -0.22, hjust = 0, size = 3.7, colour = "#4D4D4D") +
  geom_text(aes(x = ann_x, label = paste0(effect_direction, "\nFDR=", fdr_lab)),
            hjust = 0, size = 3.8, lineheight = 1.05, colour = "#333333") +
  facet_wrap(~context, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = col_main) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  coord_cartesian(clip = "off") +
  labs(
    title = wrap_title("Micro->Gaba exact is geometrically distinct, whereas Micro->Glut exact shows much weaker compactness-specific geometry"),
    subtitle = wrap_subtitle("Within exact nonneuron->neuron events, Microglia is compared against all other nonneuronal child classes. Positive evidence for a compact functional core in Micro->Gaba should appear as higher Jaccard, smaller host/child mismatch, and larger Microglia child size."),
    x = "Median value",
    y = NULL
  ) +
  base_theme() +
  theme(strip.text = element_text(size = 14.5, face = "bold"),
        plot.margin = margin(18, 64, 18, 34))

save_pub(fig1, "fig01_micro_exact_geometry_signature.png", 22.8, 6.3)

# ============================================================
# 4. Figure 2 | layer heatmap
# ============================================================
layer_tests <- nn |>
  filter(neuron_type %in% c("Glut", "Gaba")) |>
  group_by(neuron_type, layer_simplified) |>
  group_modify(~ {
    purrr::map_dfr(c("Non->Neuron exact", "Neuron->Non exact"), function(out) {
      a <- sum(.x$non_class == "Microglia" & .x$topology == out)
      b <- sum(.x$non_class == "Microglia" & .x$topology != out)
      c <- sum(.x$non_class != "Microglia" & .x$topology == out)
      d <- sum(.x$non_class != "Microglia" & .x$topology != out)
      fisher_tbl(a, b, c, d) |>
        mutate(
          outcome = out,
          rate_micro = safe_div(a, a + b),
          rate_other = safe_div(c, c + d)
        )
    })
  }) |>
  ungroup() |>
  group_by(neuron_type, outcome) |>
  mutate(fdr = p.adjust(p, method = "BH")) |>
  ungroup()

layer_glut_focus <- layer_tests |>
  filter(neuron_type == "Glut") |>
  mutate(
    outcome_label = recode(outcome,
                           "Non->Neuron exact" = "Micro->Glut exact",
                           "Neuron->Non exact" = "Glut->Micro exact"),
    outcome_label = factor(outcome_label, levels = c("Micro->Glut exact", "Glut->Micro exact")),
    layer_simplified = factor(layer_simplified, levels = c("1", "2/3", "4", "5", "6a", "6b")),
    tile_lab = paste0("OR=", fmt_num(or, 2), "\nFDR=", fmt_p(fdr))
  ) |>
  arrange(outcome_label, layer_simplified)

readr::write_csv(layer_tests, file.path(tab_dir, "table02_micro_layer_tests.csv"))
readr::write_csv(layer_glut_focus, file.path(tab_dir, "table03_micro_glut_layer_focus.csv"))

fig2 <- ggplot(layer_glut_focus, aes(x = layer_simplified, y = outcome_label, fill = log2(or))) +
  geom_tile(colour = "white", linewidth = 1.8) +
  geom_text(aes(label = tile_lab), size = 4.0, lineheight = 1.02) +
  scale_fill_gradient2(
    low = "#3B4CC0", mid = "#F7F7F7", high = "#C44E52", midpoint = 0,
    breaks = c(-1, 0, 1), labels = label_number(accuracy = 0.1)
  ) +
  labs(
    title = wrap_title("Microglia shows a pronounced layer bias in Glut pairs, with layer 1 being the strongest hotspot"),
    subtitle = wrap_subtitle("Odds ratios compare Microglia against all other nonneuronal classes within each layer. Red indicates enrichment, blue indicates depletion. The expected signature is simultaneous enrichment of Micro->Glut exact and depletion of Glut->Micro exact in shallow layers."),
    x = "Layer",
    y = NULL,
    fill = "Enrichment / depletion\n(log2 OR)"
  ) +
  base_theme() +
  theme(
    axis.text.y = element_text(face = "bold", size = 13.5),
    axis.text.x = element_text(size = 13.5),
    legend.position = "top",
    legend.key.width = unit(2.2, "cm"),
    plot.margin = margin(16, 18, 14, 18)
  )

save_pub(fig2, "fig02_micro_glut_layer_heatmap.png", 11.2, 5.8)

# ============================================================
# 5. Figure 3 | subclass forest
# ============================================================
subclass_tests <- purrr::map_dfr(c("Glut", "Gaba"), function(nt) {
  sub <- nn |> filter(neuron_type == nt)
  valid_subclasses <- sub |>
    count(neuron_subclass, name = "n_pairs") |>
    filter(n_pairs >= 50) |>
    pull(neuron_subclass)
  
  purrr::map_dfr(valid_subclasses, function(ns) {
    sub2 <- sub |> filter(neuron_subclass == ns)
    purrr::map_dfr(c("Non->Neuron exact", "Neuron->Non exact"), function(out) {
      a <- sum(sub2$non_class == "Microglia" & sub2$topology == out)
      b <- sum(sub2$non_class == "Microglia" & sub2$topology != out)
      c <- sum(sub2$non_class != "Microglia" & sub2$topology == out)
      d <- sum(sub2$non_class != "Microglia" & sub2$topology != out)
      fisher_tbl(a, b, c, d) |>
        mutate(
          neuron_type = nt,
          neuron_subclass = ns,
          outcome = out,
          rate_micro = safe_div(a, a + b),
          rate_other = safe_div(c, c + d),
          n_micro_pairs = a + b,
          n_total_pairs = a + b + c + d
        )
    })
  })
}) |>
  group_by(neuron_type, outcome) |>
  mutate(fdr = p.adjust(p, method = "BH")) |>
  ungroup()

readr::write_csv(subclass_tests, file.path(tab_dir, "table04_micro_subclass_tests.csv"))

focus_enriched <- c("L2/3 IT CTX", "CLA-EPd-CTX Car3", "L2/3 IT RSP")
focus_depleted <- c("L5 NP CTX", "L5 ET CTX", "L6 IT CTX")

plot_subclass_glut <- subclass_tests |>
  filter(neuron_type == "Glut", outcome == "Non->Neuron exact", neuron_subclass %in% c(focus_enriched, focus_depleted)) |>
  mutate(
    direction = if_else(neuron_subclass %in% focus_enriched, "Enriched", "Depleted"),
    neuron_subclass = factor(neuron_subclass, levels = rev(c(focus_depleted, focus_enriched))),
    ann = paste0("OR=", fmt_num(or, 2), "\nFDR=", fmt_p(fdr)),
    ann_x = if_else(or >= 1, ci_high * 1.12, ci_low / 1.12),
    ann_hjust = if_else(or >= 1, 0, 1)
  )

readr::write_csv(plot_subclass_glut, file.path(tab_dir, "table05_micro_glut_subclass_focus.csv"))

xmin3 <- min(plot_subclass_glut$ci_low, na.rm = TRUE)
xmax3 <- max(plot_subclass_glut$ci_high, na.rm = TRUE)

fig3 <- ggplot(plot_subclass_glut, aes(x = or, y = neuron_subclass, colour = direction)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.8, colour = "#7A7A7A") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.18, linewidth = 1.2, show.legend = FALSE) +
  geom_point(size = 4.0, show.legend = FALSE) +
  geom_text(aes(x = ann_x, label = ann, hjust = ann_hjust), size = 3.7, colour = "#444444", show.legend = FALSE) +
  scale_x_log10(labels = label_number(accuracy = 0.01), expand = expansion(mult = c(0.18, 0.18))) +
  scale_colour_manual(values = col_dir) +
  coord_cartesian(xlim = c(xmin3 / 1.25, xmax3 * 1.35), clip = "off") +
  labs(
    title = wrap_title("Micro->Glut exact is not uniform across excitatory subclasses"),
    subtitle = wrap_subtitle("Enrichment concentrates in shallow IT / association-like hosts, whereas several deeper subclasses are depleted."),
    x = "Odds ratio for Micro->Glut exact (log scale)",
    y = NULL
  ) +
  base_theme() +
  theme(plot.margin = margin(16, 58, 16, 26))

save_pub(fig3, "fig03_micro_glut_subclass_forest.png", 11.2, 6.4)

# ============================================================
# 6. Figure 4 | compact micro->gaba core summary
# ============================================================
core_gaba_summary <- sig_gaba_child |>
  mutate(
    metric = recode(var, !!!metric_label_map),
    effect_direction = case_when(
      var == "jaccard" ~ "Higher in Micro",
      var == "log10_neuron_non_ratio" ~ "Lower in Micro",
      var == "non_total" ~ "Higher in Micro",
      TRUE ~ NA_character_
    ),
    focal_median_label = if_else(var == "non_total", as.character(round(focal_median)), fmt_num(focal_median, 3)),
    other_median_label = if_else(var == "non_total", as.character(round(other_median)), fmt_num(other_median, 3))
  ) |>
  dplyr::select(metric, focal_median, other_median, focal_median_label, other_median_label, fdr, effect_direction) |>
  mutate(
    metric = factor(metric, levels = c("Microglia child size", "log10(host/child size ratio)", "Jaccard")),
    ann_x = c(52, 52, 52)
  )

readr::write_csv(core_gaba_summary, file.path(tab_dir, "table06_micro_gaba_compact_core_summary.csv"))

fig4 <- ggplot(core_gaba_summary, aes(y = metric)) +
  geom_segment(aes(x = other_median, xend = focal_median, yend = metric), linewidth = 1.35, colour = "#B8B8B8") +
  geom_point(aes(x = other_median), size = 4.4, colour = "#8F8F8F") +
  geom_point(aes(x = focal_median), size = 4.8, colour = "#C44E52") +
  geom_text(aes(x = other_median, label = other_median_label), nudge_y = -0.22, size = 4.0, colour = "#4D4D4D") +
  geom_text(aes(x = focal_median, label = focal_median_label), nudge_y = 0.22, size = 4.0, colour = "#C44E52") +
  geom_text(aes(x = ann_x, label = paste0(effect_direction, "\nFDR=", fmt_p(fdr))), hjust = 0, size = 3.9, lineheight = 1.05, colour = "#333333") +
  coord_cartesian(xlim = c(-3, 62), clip = "off") +
  labs(
    title = wrap_title("Micro->Gaba exact behaves like a compact functional core rather than a fragmented nonneuronal satellite"),
    subtitle = wrap_subtitle("Three geometric signals are highlighted: tighter overlap (higher Jaccard), smaller host-child mismatch, and larger Microglia child size, all relative to other nonneuronal child exact events in Gaba hosts."),
    x = "Median value",
    y = NULL
  ) +
  base_theme() +
  theme(plot.margin = margin(18, 62, 16, 28))

save_pub(fig4, "fig04_micro_gaba_core_summary.png", 11.0, 5.9)

# ============================================================
# 7. Console summary
# ============================================================
cat("\n==================== INPUT ====================\n")
cat("File:", input_neuron_nonneuron, "\n")
cat("Pairs used:", nrow(nn), "\n")

cat("\n==================== MICRO->GABA CORE ====================\n")
print(core_gaba_summary)

cat("\n==================== LAYER BIAS (GLUT) ====================\n")
print(layer_glut_focus |> dplyr::select(layer_simplified, outcome_label, or, ci_low, ci_high, p, fdr, rate_micro, rate_other))

cat("\n==================== GLUT SUBCLASS FOCUS ====================\n")
print(plot_subclass_glut |> dplyr::select(neuron_subclass, or, ci_low, ci_high, p, fdr, rate_micro, rate_other))

cat("\nOutputs saved to:\n")
cat("  ", fig_dir, "\n")
cat("  ", tab_dir, "\n")
