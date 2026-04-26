# 清空环境
rm(list = ls())
gc()
# ============================================================
# R_output_08
# Cross-dataset comparative analysis:
# neuron-neuron vs neuron-nonneuron spatial interaction data
# Focus:
# 1) GABA role is partner-dependent rather than fixed
# 2) Full multi-host containment obeys a distinct-host-subclass rule
# Input:
#   E:/zaw/2603/neuron-neuron-partner.csv
#   E:/zaw/2603/neuron-nonneuron-partner.csv
# Output:
#   E:/zaw/2603/R_output_08
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(forcats)
  library(scales)
  library(tibble)
})

options(scipen = 999)

# namespace guards
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

base_dir <- "E:/zaw/2603"
nn_file <- file.path(base_dir, "neuron-neuron-partner.csv")
nnon_file <- file.path(base_dir, "neuron-nonneuron-partner.csv")
out_dir <- file.path(base_dir, "R_output_08")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------
# Helper functions
# ---------------------------

dedup_unordered <- function(df) {
  a <- as.character(df$cluster.1_label)
  b <- as.character(df$cluster.2_label)
  pair_key <- ifelse(a < b, paste0(a, "||", b), paste0(b, "||", a))
  df %>%
    dplyr::mutate(pair_key = pair_key) %>%
    dplyr::distinct(pair_key, .keep_all = TRUE)
}

region_count <- function(x) {
  ifelse(
    is.na(x) | stringr::str_trim(x) == "",
    NA_integer_,
    lengths(stringr::str_split(x, ","))
  )
}

simplify_subclass <- function(x) {
  x %>%
    stringr::str_replace("^\\d+\\s+", "") %>%
    stringr::str_replace("\\s+(Glut|Gaba|NN)$", "")
}

cliffs_delta_fast <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  wt <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
  u <- unname(wt$statistic)
  (2 * u) / (length(x) * length(y)) - 1
}

paired_slide_test <- function(slide_summ, metric_name) {
  piv <- slide_summ %>%
    dplyr::select(slide, group, !!rlang::sym(metric_name)) %>%
    tidyr::pivot_wider(names_from = group, values_from = !!rlang::sym(metric_name))
  
  bind_rows(
    tibble(
      metric = metric_name,
      comparison = "GI vs GN",
      median_A = median(piv$GI, na.rm = TRUE),
      median_B = median(piv$GN, na.rm = TRUE),
      p_value = suppressWarnings(wilcox.test(piv$GI, piv$GN, paired = TRUE, exact = FALSE)$p.value)
    ),
    tibble(
      metric = metric_name,
      comparison = "GI vs IN",
      median_A = median(piv$GI, na.rm = TRUE),
      median_B = median(piv$IN, na.rm = TRUE),
      p_value = suppressWarnings(wilcox.test(piv$GI, piv$IN, paired = TRUE, exact = FALSE)$p.value)
    ),
    tibble(
      metric = metric_name,
      comparison = "GN vs IN",
      median_A = median(piv$GN, na.rm = TRUE),
      median_B = median(piv$IN, na.rm = TRUE),
      p_value = suppressWarnings(wilcox.test(piv$GN, piv$IN, paired = TRUE, exact = FALSE)$p.value)
    )
  )
}

wrap_plot_title <- function(x, width = 56) stringr::str_wrap(x, width = width)
wrap_plot_subtitle <- function(x, width = 95) stringr::str_wrap(x, width = width)

plot_theme_safe <- ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    plot.title.position = "plot",
    plot.margin = ggplot2::margin(20, 28, 14, 18),
    plot.title = ggplot2::element_text(face = "bold", size = 16.5, lineheight = 1.06, margin = ggplot2::margin(b = 8)),
    plot.subtitle = ggplot2::element_text(size = 11.6, colour = "#4D4D4D", lineheight = 1.10, margin = ggplot2::margin(b = 10)),
    axis.title = ggplot2::element_text(size = 13.2, colour = "black"),
    axis.text = ggplot2::element_text(size = 11.5, colour = "black"),
    legend.position = "top",
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 11.2)
  )

save_pub <- function(plot_obj, filename, width = 8, height = 5.2, dpi = 600) {
  ggsave(filename = file.path(fig_dir, paste0(filename, ".png")),
         plot = plot_obj, width = width, height = height, dpi = dpi,
         bg = "white", limitsize = FALSE)
  ggsave(filename = file.path(fig_dir, paste0(filename, ".pdf")),
         plot = plot_obj, width = width, height = height,
         bg = "white", device = cairo_pdf, limitsize = FALSE)
}

fmt_p_short <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-300) return("<1e-300")
  if (p < 1e-4) return(format(p, scientific = TRUE, digits = 2))
  sprintf("%.4f", p)
}

# ---------------------------
# Load and deduplicate
# ---------------------------

nn <- readr::read_csv(nn_file, show_col_types = FALSE)
nnon <- readr::read_csv(nnon_file, show_col_types = FALSE)

nn_u <- dedup_unordered(nn)
nnon_u <- dedup_unordered(nnon)

add_basic_features <- function(df) {
  df %>%
    dplyr::mutate(
      `cluster.1.region_n` = region_count(cluster.1_region),
      `cluster.2.region_n` = region_count(cluster.2_region),
      `cluster.1.log10p` = -log10(pmax(cluster.1_cauchy_combination_p, 1e-300)),
      `cluster.2.log10p` = -log10(pmax(cluster.2_cauchy_combination_p, 1e-300)),
      `cluster.1.log_size` = log10(cluster.1_total_cell_num + 1),
      `cluster.2.log_size` = log10(cluster.2_total_cell_num + 1),
      abs_asym = abs(cluster.1.overlap.percent - cluster.2.overlap.percent),
      c1_in_c2 = cluster.1.overlap.percent >= 0.999999,
      c2_in_c1 = cluster.2.overlap.percent >= 0.999999,
      any_containment = c1_in_c2 | c2_in_c1
    )
}

nn_u <- add_basic_features(nn_u)
nnon_u <- add_basic_features(nnon_u)

nn_u <- nn_u %>%
  dplyr::mutate(pair_class = purrr::map2_chr(
    cluster.1_cell_Neruon_type, cluster.2_cell_Neruon_type,
    ~ paste(sort(c(.x, .y)), collapse = "-")
  ))

nnon_u <- nnon_u %>%
  dplyr::mutate(pair_class = purrr::map2_chr(
    cluster.1_cell_Neruon_type, cluster.2_cell_Neruon_type,
    ~ paste(sort(c(.x, .y)), collapse = "-")
  ))

# ---------------------------
# Orient Gaba-Glut pairs (GI)
# ---------------------------

gi <- nn_u %>%
  dplyr::filter(pair_class == "Gaba-Glut") %>%
  dplyr::mutate(is_c1_glut = cluster.1_cell_Neruon_type == "Glut") %>%
  dplyr::transmute(
    pair_key,
    slide = if_else(is_c1_glut, cluster.1_slide, cluster.2_slide),
    layer = if_else(is_c1_glut, cluster.1_layer, cluster.2_layer),
    E_totalcellnum = if_else(is_c1_glut, cluster.1_total_cell_num, cluster.2_total_cell_num),
    I_totalcellnum = if_else(is_c1_glut, cluster.2_total_cell_num, cluster.1_total_cell_num),
    E_overlappercent = if_else(is_c1_glut, cluster.1.overlap.percent, cluster.2.overlap.percent),
    I_overlappercent = if_else(is_c1_glut, cluster.2.overlap.percent, cluster.1.overlap.percent),
    jaccard, overlap_cell, abs_asym,
    I_in_E = if_else(is_c1_glut, c2_in_c1, c1_in_c2),
    E_in_I = if_else(is_c1_glut, c1_in_c2, c2_in_c1)
  ) %>%
  dplyr::mutate(
    signed_asym = E_overlappercent - I_overlappercent,
    E_larger = E_totalcellnum >= I_totalcellnum
  )

# ---------------------------
# Orient neuron-nonneuron pairs: GN and IN
# ---------------------------

nnon_het <- nnon_u %>%
  dplyr::filter(pair_class %in% c("Glut-NonNeuron", "Gaba-NonNeuron")) %>%
  dplyr::mutate(is_c1_neuron = cluster.1_cell_Neruon_type %in% c("Glut", "Gaba")) %>%
  dplyr::transmute(
    pair_key,
    neuron_type = if_else(is_c1_neuron, cluster.1_cell_Neruon_type, cluster.2_cell_Neruon_type),
    neuron_label = if_else(is_c1_neuron, cluster.1_label, cluster.2_label),
    non_label = if_else(is_c1_neuron, cluster.2_label, cluster.1_label),
    slide = if_else(is_c1_neuron, cluster.1_slide, cluster.2_slide),
    layer = if_else(is_c1_neuron, cluster.1_layer, cluster.2_layer),
    neuron_totalcellnum = if_else(is_c1_neuron, cluster.1_total_cell_num, cluster.2_total_cell_num),
    non_totalcellnum = if_else(is_c1_neuron, cluster.2_total_cell_num, cluster.1_total_cell_num),
    neuron_overlappercent = if_else(is_c1_neuron, cluster.1.overlap.percent, cluster.2.overlap.percent),
    non_overlappercent = if_else(is_c1_neuron, cluster.2.overlap.percent, cluster.1.overlap.percent),
    neuron_in_non = if_else(is_c1_neuron, c1_in_c2, c2_in_c1),
    non_in_neuron = if_else(is_c1_neuron, c2_in_c1, c1_in_c2),
    jaccard, overlap_cell, abs_asym
  ) %>%
  dplyr::mutate(
    signed_asym = neuron_overlappercent - non_overlappercent,
    neuron_larger = neuron_totalcellnum >= non_totalcellnum
  )

gn <- nnon_het %>% dplyr::filter(neuron_type == "Glut")
inn <- nnon_het %>% dplyr::filter(neuron_type == "Gaba")

# ---------------------------
# 1) Partner-dependent directionality
# ---------------------------

direction_summary <- bind_rows(
  tibble(
    system = "GI",
    direction = c("Gaba ⊂ Glut", "Glut ⊂ Gaba"),
    count = c(sum(gi$I_in_E), sum(gi$E_in_I))
  ),
  tibble(
    system = "IN",
    direction = c("Gaba ⊂ NonNeuron", "NonNeuron ⊂ Gaba"),
    count = c(sum(inn$neuron_in_non), sum(inn$non_in_neuron))
  )
) %>%
  dplyr::group_by(system) %>%
  dplyr::mutate(prop = count / sum(count)) %>%
  dplyr::ungroup()

direction_tests <- bind_rows(
  tibble(
    system = "GI",
    total_oneway = sum(gi$I_in_E) + sum(gi$E_in_I),
    dominant_count = sum(gi$I_in_E),
    reverse_count = sum(gi$E_in_I),
    dominant_prop = sum(gi$I_in_E) / (sum(gi$I_in_E) + sum(gi$E_in_I)),
    p_value = binom.test(sum(gi$I_in_E), sum(gi$I_in_E) + sum(gi$E_in_I), p = 0.5)$p.value
  ),
  tibble(
    system = "IN",
    total_oneway = sum(inn$neuron_in_non) + sum(inn$non_in_neuron),
    dominant_count = sum(inn$neuron_in_non),
    reverse_count = sum(inn$non_in_neuron),
    dominant_prop = sum(inn$neuron_in_non) / (sum(inn$neuron_in_non) + sum(inn$non_in_neuron)),
    p_value = binom.test(sum(inn$neuron_in_non), sum(inn$neuron_in_non) + sum(inn$non_in_neuron), p = 0.5)$p.value
  )
)

# ---------------------------
# 2) Partial asymmetry not transferable from GI to IN
# ---------------------------

partial_gi <- gi %>% dplyr::filter(!(I_in_E | E_in_I))
partial_gn <- gn %>% dplyr::filter(!(neuron_in_non | non_in_neuron))
partial_in <- inn %>% dplyr::filter(!(neuron_in_non | non_in_neuron))

partial_pair_summary <- tibble(
  group = c("GI", "GN", "IN"),
  n_pairs = c(nrow(partial_gi), nrow(partial_gn), nrow(partial_in)),
  median_signed_asym = c(median(partial_gi$signed_asym), median(partial_gn$signed_asym), median(partial_in$signed_asym)),
  pair_level_p_vs_zero = c(
    ifelse(any(partial_gi$signed_asym != 0),
           wilcox.test(partial_gi$signed_asym[partial_gi$signed_asym != 0], mu = 0, exact = FALSE)$p.value, 1),
    ifelse(any(partial_gn$signed_asym != 0),
           wilcox.test(partial_gn$signed_asym[partial_gn$signed_asym != 0], mu = 0, exact = FALSE)$p.value, 1),
    ifelse(any(partial_in$signed_asym != 0),
           wilcox.test(partial_in$signed_asym[partial_in$signed_asym != 0], mu = 0, exact = FALSE)$p.value, 1)
  )
)

slide_partial <- bind_rows(
  partial_gi %>% dplyr::group_by(slide) %>% dplyr::summarise(med = median(signed_asym), .groups = "drop") %>% dplyr::mutate(group = "GI"),
  partial_gn %>% dplyr::group_by(slide) %>% dplyr::summarise(med = median(signed_asym), .groups = "drop") %>% dplyr::mutate(group = "GN"),
  partial_in %>% dplyr::group_by(slide) %>% dplyr::summarise(med = median(signed_asym), .groups = "drop") %>% dplyr::mutate(group = "IN")
)

slide_partial_summary <- slide_partial %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(
    n_slides = dplyr::n(),
    median_slide_signed_asym = median(med),
    slide_level_p_vs_zero = ifelse(any(med != 0), wilcox.test(med[med != 0], mu = 0, exact = FALSE)$p.value, 1),
    .groups = "drop"
  )

slide_partial_pairwise <- bind_rows(
  paired_slide_test(slide_partial %>% dplyr::rename(contain_rate = med), "contain_rate")
) %>%
  dplyr::mutate(metric = "median_signed_asym")

# ---------------------------
# 3) Distinct-host-subclass law in full containment vs overlap background
# ---------------------------

contain_edges <- bind_rows(
  nn_u %>%
    dplyr::filter(c1_in_c2) %>%
    dplyr::transmute(
      source = "neuron-neuron",
      guest_label = cluster.1_label,
      guest_type = cluster.1_cell_Neruon_type,
      host_label = cluster.2_label,
      host_type = cluster.2_cell_Neruon_type,
      host_subclass = cluster.2_subclass,
      edge_class = paste0(cluster.1_cell_Neruon_type, "_in_", cluster.2_cell_Neruon_type)
    ),
  nn_u %>%
    dplyr::filter(c2_in_c1) %>%
    dplyr::transmute(
      source = "neuron-neuron",
      guest_label = cluster.2_label,
      guest_type = cluster.2_cell_Neruon_type,
      host_label = cluster.1_label,
      host_type = cluster.1_cell_Neruon_type,
      host_subclass = cluster.1_subclass,
      edge_class = paste0(cluster.2_cell_Neruon_type, "_in_", cluster.1_cell_Neruon_type)
    ),
  nnon_u %>%
    dplyr::filter(c1_in_c2) %>%
    dplyr::transmute(
      source = "neuron-nonneuron",
      guest_label = cluster.1_label,
      guest_type = cluster.1_cell_Neruon_type,
      host_label = cluster.2_label,
      host_type = cluster.2_cell_Neruon_type,
      host_subclass = cluster.2_subclass,
      edge_class = paste0(cluster.1_cell_Neruon_type, "_in_", cluster.2_cell_Neruon_type)
    ),
  nnon_u %>%
    dplyr::filter(c2_in_c1) %>%
    dplyr::transmute(
      source = "neuron-nonneuron",
      guest_label = cluster.2_label,
      guest_type = cluster.2_cell_Neruon_type,
      host_label = cluster.1_label,
      host_type = cluster.1_cell_Neruon_type,
      host_subclass = cluster.1_subclass,
      edge_class = paste0(cluster.2_cell_Neruon_type, "_in_", cluster.1_cell_Neruon_type)
    )
) %>%
  dplyr::mutate(host_subclass_simple = simplify_subclass(host_subclass))

# helper long table for overlap neighborhoods
long_all <- bind_rows(
  nn_u %>%
    dplyr::transmute(
      source = "neuron-neuron",
      focal_label = cluster.1_label,
      focal_type = cluster.1_cell_Neruon_type,
      partner_label = cluster.2_label,
      partner_type = cluster.2_cell_Neruon_type,
      partner_subclass = cluster.2_subclass
    ),
  nn_u %>%
    dplyr::transmute(
      source = "neuron-neuron",
      focal_label = cluster.2_label,
      focal_type = cluster.2_cell_Neruon_type,
      partner_label = cluster.1_label,
      partner_type = cluster.1_cell_Neruon_type,
      partner_subclass = cluster.1_subclass
    ),
  nnon_u %>%
    dplyr::transmute(
      source = "neuron-nonneuron",
      focal_label = cluster.1_label,
      focal_type = cluster.1_cell_Neruon_type,
      partner_label = cluster.2_label,
      partner_type = cluster.2_cell_Neruon_type,
      partner_subclass = cluster.2_subclass
    ),
  nnon_u %>%
    dplyr::transmute(
      source = "neuron-nonneuron",
      focal_label = cluster.2_label,
      focal_type = cluster.2_cell_Neruon_type,
      partner_label = cluster.1_label,
      partner_type = cluster.1_cell_Neruon_type,
      partner_subclass = cluster.1_subclass
    )
) %>%
  dplyr::mutate(partner_subclass_simple = simplify_subclass(partner_subclass))

long_overlap_groups_for_edge <- function(edge_name) {
  parts <- stringr::str_split(edge_name, "_in_", simplify = TRUE)
  guest_type <- parts[1]
  host_type <- parts[2]
  src <- ifelse(
    edge_name %in% c("Gaba_in_Glut", "Glut_in_Gaba", "Glut_in_Glut", "Gaba_in_Gaba"),
    "neuron-neuron",
    "neuron-nonneuron"
  )
  
  long_all %>%
    dplyr::filter(source == src, focal_type == guest_type, partner_type == host_type) %>%
    dplyr::distinct(focal_label, partner_label, partner_subclass_simple) %>%
    dplyr::group_by(focal_label) %>%
    dplyr::summarise(
      n_partners = dplyr::n_distinct(partner_label),
      n_subcls = dplyr::n_distinct(partner_subclass_simple),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_partners >= 2) %>%
    dplyr::mutate(repeat_same_subclass = n_subcls < n_partners)
}

multihost_repeat_by_edge_class <- contain_edges %>%
  dplyr::group_by(edge_class) %>%
  dplyr::group_modify(function(.x, .y) {
    grp_cont <- .x %>%
      dplyr::distinct(guest_label, host_label, host_subclass_simple) %>%
      dplyr::group_by(guest_label) %>%
      dplyr::summarise(
        n_partners = dplyr::n_distinct(host_label),
        n_subcls = dplyr::n_distinct(host_subclass_simple),
        .groups = "drop"
      ) %>%
      dplyr::filter(n_partners >= 2) %>%
      dplyr::mutate(repeat_same_subclass = n_subcls < n_partners)
    
    grp_ov <- long_overlap_groups_for_edge(.y$edge_class)
    
    tibble::tibble(
      contain_multi_n = nrow(grp_cont),
      contain_repeat_rate = ifelse(nrow(grp_cont) > 0, mean(grp_cont$repeat_same_subclass), NA_real_),
      overlap_multi_n = nrow(grp_ov),
      overlap_repeat_rate = ifelse(nrow(grp_ov) > 0, mean(grp_ov$repeat_same_subclass), NA_real_),
      contain_repeat_count = ifelse(nrow(grp_cont) > 0, sum(grp_cont$repeat_same_subclass), 0),
      overlap_repeat_count = ifelse(nrow(grp_ov) > 0, sum(grp_ov$repeat_same_subclass), 0)
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(dplyr::desc(contain_multi_n))

contain_multi_n <- sum(multihost_repeat_by_edge_class$contain_multi_n, na.rm = TRUE)
contain_repeat_n <- sum(multihost_repeat_by_edge_class$contain_repeat_count, na.rm = TRUE)
overlap_multi_n <- sum(multihost_repeat_by_edge_class$overlap_multi_n, na.rm = TRUE)
overlap_repeat_n <- sum(multihost_repeat_by_edge_class$overlap_repeat_count, na.rm = TRUE)

multihost_repeat_global <- tibble(
  containment_multi_guest_n = contain_multi_n,
  containment_repeat_same_subclass_n = contain_repeat_n,
  containment_repeat_rate = contain_repeat_n / contain_multi_n,
  overlap_multi_guest_n = overlap_multi_n,
  overlap_repeat_same_subclass_n = overlap_repeat_n,
  overlap_repeat_rate = overlap_repeat_n / overlap_multi_n,
  fisher_p_value = fisher.test(matrix(c(
    contain_repeat_n, contain_multi_n - contain_repeat_n,
    overlap_repeat_n, overlap_multi_n - overlap_repeat_n
  ), 2, 2))$p.value
)

# ---------------------------
# Save core tables
# ---------------------------
write_csv(direction_summary, file.path(tab_dir, "table_01_direction_summary_GI_IN.csv"))
write_csv(direction_tests, file.path(tab_dir, "table_02_direction_tests_GI_IN.csv"))
write_csv(partial_pair_summary, file.path(tab_dir, "table_03_partial_pair_asymmetry_summary.csv"))
write_csv(slide_partial_summary, file.path(tab_dir, "table_04_partial_slide_asymmetry_summary.csv"))
write_csv(slide_partial_pairwise, file.path(tab_dir, "table_05_partial_slide_pairwise.csv"))
write_csv(multihost_repeat_global, file.path(tab_dir, "table_06_multihost_repeat_global.csv"))
write_csv(multihost_repeat_by_edge_class, file.path(tab_dir, "table_07_multihost_repeat_by_edge_class.csv"))

figure_manifest <- tibble(
  filename = c(
    "fig01_GI_vs_IN_direction_counts",
    "fig02_GI_vs_IN_direction_proportions",
    "fig03_partial_signed_asymmetry_boxplots",
    "fig04_partial_signed_asymmetry_slide",
    "fig05_multihost_repeat_global",
    "fig06_multihost_repeat_by_edge_class"
  ),
  summary = c(
    "Directional full-containment counts for GI and IN",
    "Directional composition within GI and IN exact one-way edges",
    "Pair-level partial signed asymmetry across GI, GN, IN",
    "Slide-level partial signed asymmetry across GI, GN, IN",
    "Global distinct-host-subclass law: full containment vs overlap background",
    "Distinct-host-subclass law by edge class"
  )
)
write_csv(figure_manifest, file.path(out_dir, "figure_manifest.csv"))

# ---------------------------
# Figures
# ---------------------------
col_system <- c("GI" = "#7B61A8", "IN" = "#C44E52")
col_group <- c("GI" = "#7B61A8", "GN" = "#4E79A7", "IN" = "#C44E52")
col_bg <- c("Full containment" = "#C44E52", "Overlap background" = "#9A9A9A")

# Fig1 counts
fig1 <- ggplot(direction_summary, aes(x = direction, y = count, fill = system)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64) +
  geom_text(
    aes(label = comma(count)),
    position = position_dodge(width = 0.72),
    vjust = -0.26, size = 4.2
  ) +
  scale_fill_manual(values = col_system) +
  scale_y_continuous(labels = comma_format(), expand = expansion(mult = c(0, 0.14))) +
  labs(
    title = wrap_plot_title("GABA shows a stable direction in GI but a partner-dependent flip in IN"),
    subtitle = wrap_plot_subtitle(
      paste0(
        "GI: Gaba ⊂ Glut = ", comma(sum(gi$I_in_E)),
        " vs Glut ⊂ Gaba = ", comma(sum(gi$E_in_I)),
        "; IN: Gaba ⊂ NonNeuron = ", comma(sum(inn$neuron_in_non)),
        " vs NonNeuron ⊂ Gaba = ", comma(sum(inn$non_in_neuron)), "."
      )
    ),
    x = NULL,
    y = "Number of one-way exact containment edges"
  ) +
  plot_theme_safe +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))
save_pub(fig1, "fig01_GI_vs_IN_direction_counts", 9.2, 5.8)

# Fig2 proportions
fig2 <- ggplot(direction_summary, aes(x = system, y = prop, fill = direction)) +
  geom_col(width = 0.64, position = "stack") +
  geom_text(
    aes(label = percent(prop, accuracy = 0.1)),
    position = position_stack(vjust = 0.5),
    colour = "white", fontface = "bold", size = 4.1
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = wrap_plot_title("The dominant containment direction flips once GABA interacts with NonNeuron partners"),
    subtitle = wrap_plot_subtitle(
      paste0(
        "GI binomial p = ", fmt_p_short(direction_tests$p_value[direction_tests$system == "GI"]),
        "; IN binomial p = ", fmt_p_short(direction_tests$p_value[direction_tests$system == "IN"]),
        "."
      )
    ),
    x = NULL,
    y = "Proportion within one-way exact edges"
  ) +
  plot_theme_safe
save_pub(fig2, "fig02_GI_vs_IN_direction_proportions", 7.0, 5.6)

# Fig3 pair-level partial asymmetry
partial_long <- bind_rows(
  tibble(group = "GI", value = partial_gi$signed_asym),
  tibble(group = "GN", value = partial_gn$signed_asym),
  tibble(group = "IN", value = partial_in$signed_asym)
) %>%
  dplyr::filter(is.finite(value))

med_labs <- partial_long %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(med = median(value), .groups = "drop")

fig3 <- ggplot(partial_long, aes(x = group, y = value, fill = group)) +
  geom_violin(alpha = 0.18, color = NA, trim = FALSE) +
  geom_boxplot(width = 0.24, outlier.shape = NA, linewidth = 0.7, alpha = 0.90) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.55, color = "#7F7F7F") +
  geom_text(
    data = med_labs,
    aes(x = group, y = max(partial_long$value, na.rm = TRUE) * 0.92, label = paste0("median=", sprintf("%.3f", med))),
    inherit.aes = FALSE, size = 4.0
  ) +
  scale_fill_manual(values = col_group) +
  labs(
    title = wrap_plot_title("IN partial overlap signed asymmetry is near zero rather than a simple extension of E-I geometry"),
    subtitle = wrap_plot_subtitle(
      paste0(
        "IN pair-level median = ", sprintf("%.3f", partial_pair_summary$median_signed_asym[partial_pair_summary$group == "IN"]),
        "; pair-level p = ", fmt_p_short(partial_pair_summary$pair_level_p_vs_zero[partial_pair_summary$group == "IN"]),
        "."
      )
    ),
    x = NULL,
    y = "Signed asymmetry"
  ) +
  plot_theme_safe +
  theme(legend.position = "none")
save_pub(fig3, "fig03_partial_signed_asymmetry_boxplots", 8.2, 5.8)

# Fig4 slide-level partial asymmetry
slide_meds <- slide_partial %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(med = median(med), .groups = "drop")

fig4 <- ggplot(slide_partial, aes(x = group, y = med, fill = group)) +
  geom_violin(alpha = 0.18, color = NA, trim = FALSE) +
  geom_boxplot(width = 0.24, outlier.shape = NA, linewidth = 0.7, alpha = 0.90) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.55, color = "#7F7F7F") +
  geom_text(
    data = slide_meds,
    aes(x = group, y = max(slide_partial$med, na.rm = TRUE) * 0.92, label = paste0("median=", sprintf("%.3f", med))),
    inherit.aes = FALSE, size = 4.0
  ) +
  scale_fill_manual(values = col_group) +
  labs(
    title = wrap_plot_title("At the slide level, IN still centers on zero asymmetry"),
    subtitle = wrap_plot_subtitle(
      paste0(
        "IN slide-level median = ", sprintf("%.3f", slide_partial_summary$median_slide_signed_asym[slide_partial_summary$group == "IN"]),
        "; slide-level p = ", fmt_p_short(slide_partial_summary$slide_level_p_vs_zero[slide_partial_summary$group == "IN"]),
        "."
      )
    ),
    x = NULL,
    y = "Slide-level median signed asymmetry"
  ) +
  plot_theme_safe +
  theme(legend.position = "none")
save_pub(fig4, "fig04_partial_signed_asymmetry_slide", 8.2, 5.8)

# Fig5 global multihost repeat law
fig5_df <- tibble(
  mode = c("Full containment", "Overlap background"),
  repeat_rate = c(multihost_repeat_global$containment_repeat_rate, multihost_repeat_global$overlap_repeat_rate),
  n_multi = c(multihost_repeat_global$containment_multi_guest_n, multihost_repeat_global$overlap_multi_guest_n),
  n_repeat = c(multihost_repeat_global$containment_repeat_same_subclass_n, multihost_repeat_global$overlap_repeat_same_subclass_n)
)

fig5 <- ggplot(fig5_df, aes(x = mode, y = repeat_rate, fill = mode)) +
  geom_col(width = 0.62) +
  geom_text(
    aes(label = paste0(comma(n_repeat), "/", comma(n_multi))),
    vjust = -0.24, size = 4.2
  ) +
  scale_fill_manual(values = col_bg) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.14))) +
  labs(
    title = wrap_plot_title("Repeated hosts of the same guest almost never come from the same subclass once full containment is imposed"),
    subtitle = wrap_plot_subtitle(
      paste0(
        "Global Fisher p = ", fmt_p_short(multihost_repeat_global$fisher_p_value),
        "; full containment repeat count = ", comma(multihost_repeat_global$containment_repeat_same_subclass_n),
        " across ", comma(multihost_repeat_global$containment_multi_guest_n),
        " multi-host guests."
      )
    ),
    x = NULL,
    y = "Rate of repeated same-host-subclass usage"
  ) +
  plot_theme_safe +
  theme(legend.position = "none")
save_pub(fig5, "fig05_multihost_repeat_global", 8.8, 5.6)

# Fig6 by edge class
mr_long <- multihost_repeat_by_edge_class %>%
  dplyr::select(edge_class, contain_repeat_rate, overlap_repeat_rate) %>%
  tidyr::pivot_longer(
    cols = c(contain_repeat_rate, overlap_repeat_rate),
    names_to = "kind", values_to = "rate"
  ) %>%
  dplyr::mutate(
    kind = dplyr::recode(
      kind,
      contain_repeat_rate = "Full containment",
      overlap_repeat_rate = "Overlap background"
    )
  )

fig6 <- ggplot(mr_long, aes(x = edge_class, y = rate, fill = kind)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64) +
  scale_fill_manual(values = col_bg) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.10))) +
  labs(
    title = wrap_plot_title("The distinct-host-subclass law holds across multi-host edge classes"),
    subtitle = wrap_plot_subtitle("Full-containment neighborhoods suppress same-subclass host repetition relative to ordinary overlap neighborhoods."),
    x = NULL,
    y = "Rate of repeated same-host-subclass usage"
  ) +
  plot_theme_safe +
  theme(axis.text.x = element_text(angle = 26, hjust = 1))
save_pub(fig6, "fig06_multihost_repeat_by_edge_class", 10.6, 6.0)

message("All outputs saved to: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
