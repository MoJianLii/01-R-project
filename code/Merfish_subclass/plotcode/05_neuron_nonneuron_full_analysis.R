# 清空环境
rm(list = ls())
gc()
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(scales)
})

options(stringsAsFactors = FALSE, scipen = 999)

# 防止与其它包的同名函数冲突
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
transmute <- dplyr::transmute
pull <- dplyr::pull

# =========================
# Paths
# =========================
root_dir <- "E:/zaw/2603"
input_file <- file.path(root_dir, "neuron-nonneuron-partner.csv")
out_dir <- file.path(root_dir, "R_output_05")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# =========================
# Helper functions
# =========================
clean_subclass <- function(x) {
  x %>%
    str_replace("^\\d+\\s+", "") %>%
    str_replace_all("_", " ")
}

safe_fisher <- function(a, b, c, d) {
  mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  ft <- fisher.test(mat)
  tibble(OR = unname(ft$estimate), p = ft$p.value)
}

pair_expand <- function(labels) {
  labels <- sort(unique(labels))
  if (length(labels) < 2) return(tibble(a = character(), b = character()))
  cmb <- t(combn(labels, 2))
  tibble(a = cmb[, 1], b = cmb[, 2])
}

same_subclass_multi <- function(df_pairs) {
  df_pairs %>%
    distinct(neuron_type, neuron_label, nn_label, nn_subclass_clean) %>%
    group_by(neuron_type, neuron_label) %>%
    summarise(
      n_partners = n(),
      n_subclasses = n_distinct(nn_subclass_clean),
      same_dup = n_subclasses < n_partners,
      .groups = "drop"
    ) %>%
    filter(n_partners >= 2)
}

nn_group_map <- function(x) {
  case_when(
    str_detect(x, "Astro") ~ "Astrocyte",
    str_detect(x, "Oligo|OPC") ~ "Oligodendro-lineage",
    str_detect(x, "SMC|VLMC|Peri|Endo") ~ "Vascular-associated",
    str_detect(x, "Micro|Macro|Immune") ~ "Immune/myeloid",
    TRUE ~ "Other"
  )
}

fmt_p <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 1e-4,
                format(p, scientific = TRUE, digits = 2),
                sprintf("%.4f", p)))
}

save_pub <- function(p, filename, width, height, dpi = 600) {
  ggsave(file.path(fig_dir, paste0(filename, ".png")), p,
         width = width, height = height, dpi = dpi, bg = "white")
  ggsave(file.path(fig_dir, paste0(filename, ".pdf")), p,
         width = width, height = height, bg = "white", device = cairo_pdf)
}

base_theme <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0),
    plot.subtitle = element_text(size = 11.5, colour = "#4D4D4D"),
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 11.5, colour = "black"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    plot.margin = margin(15, 20, 12, 12)
  )

# =========================
# Read data
# =========================
df <- read_csv(input_file, show_col_types = FALSE)

# Unique cluster table
clusters <- bind_rows(
  df %>% transmute(
    label = cluster.1_label,
    cell_Neruon_type = cluster.1_cell_Neruon_type,
    subclass = cluster.1_subclass,
    subclass_clean = clean_subclass(cluster.1_subclass),
    slide = cluster.1_slide,
    layer = cluster.1_layer,
    total_cell_num = cluster.1_total_cell_num
  ),
  df %>% transmute(
    label = cluster.2_label,
    cell_Neruon_type = cluster.2_cell_Neruon_type,
    subclass = cluster.2_subclass,
    subclass_clean = clean_subclass(cluster.2_subclass),
    slide = cluster.2_slide,
    layer = cluster.2_layer,
    total_cell_num = cluster.2_total_cell_num
  )
) %>% distinct(label, .keep_all = TRUE)

# Oriented neuron->nonneuron heterotypic pairs
het <- df %>%
  filter(cluster.1_cell_Neruon_type %in% c("Glut", "Gaba"),
         cluster.2_cell_Neruon_type == "NonNeuron") %>%
  mutate(
    interaction_class = case_when(
      cluster.1.overlap.percent >= 1 & cluster.2.overlap.percent >= 1 ~ "mutual_equal",
      cluster.1.overlap.percent >= 1 ~ "Neuron_in_NN",
      cluster.2.overlap.percent >= 1 ~ "NN_in_Neuron",
      TRUE ~ "Partial"
    ),
    neuron_type = cluster.1_cell_Neruon_type,
    neuron_subclass = cluster.1_subclass,
    neuron_subclass_clean = clean_subclass(cluster.1_subclass),
    nn_subclass = cluster.2_subclass,
    nn_subclass_clean = clean_subclass(cluster.2_subclass),
    neuron_size = cluster.1_total_cell_num,
    nn_size = cluster.2_total_cell_num,
    neuron_label = cluster.1_label,
    nn_label = cluster.2_label,
    layer = cluster.1_layer,
    slide = cluster.1_slide,
    neuron_region = cluster.1_region,
    nn_region = cluster.2_region,
    containment_event = interaction_class %in% c("Neuron_in_NN", "NN_in_Neuron", "mutual_equal"),
    size_ratio_nn_to_neuron = nn_size / neuron_size,
    nn_group = nn_group_map(clean_subclass(nn_subclass))
  )

# Unique NN-NN pairs for host-host lookups
nnn <- df %>%
  filter(cluster.1_cell_Neruon_type == "NonNeuron",
         cluster.2_cell_Neruon_type == "NonNeuron") %>%
  mutate(
    min_label = pmin(cluster.1_label, cluster.2_label),
    max_label = pmax(cluster.1_label, cluster.2_label)
  ) %>%
  distinct(min_label, max_label, .keep_all = TRUE)

nin <- het %>% filter(interaction_class == "Neuron_in_NN")
partial <- het %>% filter(interaction_class == "Partial")
nnin <- het %>% filter(interaction_class == "NN_in_Neuron")

# =========================
# Core summaries for the teacher's key claim
# =========================
mult_nin <- same_subclass_multi(nin)
mult_partial <- same_subclass_multi(partial)
mult_nnin <- same_subclass_multi(nnin)

overall_dup <- bind_rows(
  mult_nin %>% summarise(context = "Full host", n_total = n(), n_same = sum(same_dup), rate = mean(same_dup)),
  mult_partial %>% summarise(context = "Partial background", n_total = n(), n_same = sum(same_dup), rate = mean(same_dup)),
  mult_nnin %>% summarise(context = "Full guest", n_total = n(), n_same = sum(same_dup), rate = mean(same_dup))
)

dup_ft <- safe_fisher(
  sum(mult_nin$same_dup),
  sum(!mult_nin$same_dup),
  sum(mult_partial$same_dup),
  sum(!mult_partial$same_dup)
)

# host-host pairs for neurons with >=2 full hosts
host_pairs <- nin %>%
  distinct(neuron_type, neuron_label, nn_label, nn_subclass_clean) %>%
  group_by(neuron_type, neuron_label) %>%
  summarise(host_labels = list(sort(unique(nn_label))), .groups = "drop") %>%
  mutate(pair_tbl = map(host_labels, pair_expand)) %>%
  select(neuron_type, neuron_label, pair_tbl) %>%
  unnest(pair_tbl) %>%
  rename(host1 = a, host2 = b)

hosthost <- host_pairs %>%
  mutate(
    min_label = pmin(host1, host2),
    max_label = pmax(host1, host2)
  ) %>%
  left_join(nnn %>% select(min_label, max_label, jaccard), by = c("min_label", "max_label"))

hosthost_summary <- tibble(
  group = c("All NN-NN", "GABA multi-host", "Glut multi-host"),
  n_pairs = c(
    nrow(nnn),
    sum(hosthost$neuron_type == "Gaba" & !is.na(hosthost$jaccard)),
    sum(hosthost$neuron_type == "Glut" & !is.na(hosthost$jaccard))
  ),
  median_jaccard = c(
    median(nnn$jaccard, na.rm = TRUE),
    median(hosthost$jaccard[hosthost$neuron_type == "Gaba"], na.rm = TRUE),
    median(hosthost$jaccard[hosthost$neuron_type == "Glut"], na.rm = TRUE)
  ),
  mean_jaccard = c(
    mean(nnn$jaccard, na.rm = TRUE),
    mean(hosthost$jaccard[hosthost$neuron_type == "Gaba"], na.rm = TRUE),
    mean(hosthost$jaccard[hosthost$neuron_type == "Glut"], na.rm = TRUE)
  )
)

hosthost_comp <- bind_rows(
  {
    x <- hosthost$jaccard[hosthost$neuron_type == "Gaba"]
    y <- nnn$jaccard
    tibble(comparison = "GABA multi-host vs background", p = wilcox.test(x, y, exact = FALSE)$p.value)
  },
  {
    x <- hosthost$jaccard[hosthost$neuron_type == "Glut"]
    y <- nnn$jaccard
    tibble(comparison = "Glut multi-host vs background", p = wilcox.test(x, y, exact = FALSE)$p.value)
  }
)

# co-host enrichment
cohost_enrichment <- function(nin_df, nt, min_count = 20) {
  sub <- nin_df %>% filter(neuron_type == nt)
  subtype_counts <- sub %>% count(nn_subclass_clean, sort = TRUE) %>% filter(n >= min_count)
  subtypes <- subtype_counts$nn_subclass_clean
  pres <- sub %>%
    distinct(neuron_label, nn_subclass_clean) %>%
    filter(nn_subclass_clean %in% subtypes) %>%
    mutate(value = TRUE) %>%
    pivot_wider(names_from = nn_subclass_clean, values_from = value, values_fill = FALSE)
  if (nrow(pres) == 0) return(tibble())
  rows <- list(); idx <- 1
  for (i in seq_along(subtypes)) {
    for (j in seq((i + 1), length(subtypes))) {
      if (j > length(subtypes)) next
      a_name <- subtypes[i]; b_name <- subtypes[j]
      A <- pres[[a_name]]; B <- pres[[b_name]]
      both <- sum(A & B); a_only <- sum(A & !B); b_only <- sum(!A & B); neither <- sum(!A & !B)
      ft <- safe_fisher(both, a_only, b_only, neither)
      rows[[idx]] <- tibble(
        neuron_type = nt,
        subtype_a = a_name,
        subtype_b = b_name,
        both = both,
        OR = ft$OR,
        p = ft$p
      )
      idx <- idx + 1
    }
  }
  bind_rows(rows) %>% mutate(q = p.adjust(p, method = "BH"))
}

gaba_cohost <- cohost_enrichment(nin, "Gaba", min_count = 20)
glut_cohost <- cohost_enrichment(nin, "Glut", min_count = 20)
cohost_sig <- bind_rows(gaba_cohost, glut_cohost) %>%
  mutate(pair = paste(subtype_a, subtype_b, sep = " + "),
         log2OR = log2(OR)) %>%
  filter(q < 0.05) %>%
  group_by(neuron_type) %>%
  arrange(desc(log2OR), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

# representative motif enrichment for full host outcome
motifs <- map_dfr(c("Gaba", "Glut"), function(nt) {
  sub_nt <- het %>% filter(neuron_type == nt)
  valid <- sub_nt %>% count(neuron_subclass_clean, nn_subclass_clean) %>% filter(n >= 20)
  map_dfr(seq_len(nrow(valid)), function(i) {
    ns <- valid$neuron_subclass_clean[i]
    nn <- valid$nn_subclass_clean[i]
    sub <- sub_nt %>% filter(neuron_subclass_clean == ns, nn_subclass_clean == nn)
    a <- sum(sub$interaction_class == "Neuron_in_NN")
    b <- nrow(sub) - a
    c <- sum(sub_nt$interaction_class == "Neuron_in_NN") - a
    d <- nrow(sub_nt) - nrow(sub) - c
    ft <- safe_fisher(a, b, c, d)
    tibble(
      neuron_type = nt,
      motif = paste(ns, nn, sep = " × "),
      n_total = nrow(sub),
      n_outcome = a,
      rate = a / nrow(sub),
      OR = ft$OR,
      p = ft$p
    )
  })
}) %>%
  group_by(neuron_type) %>%
  mutate(q = p.adjust(p, method = "BH")) %>%
  ungroup()

motif_top <- motifs %>%
  filter(q < 0.05, OR > 1) %>%
  mutate(log2OR = log2(OR)) %>%
  group_by(neuron_type) %>%
  arrange(desc(log2OR), .by_group = TRUE) %>%
  slice_head(n = 8) %>%
  ungroup()

# =========================
# Save key tables
# =========================
write_csv(overall_dup, file.path(tab_dir, "table_01_same_subclass_duplication_summary.csv"))
write_csv(dup_ft, file.path(tab_dir, "table_02_same_subclass_duplication_fisher.csv"))
write_csv(hosthost_summary, file.path(tab_dir, "table_03_hosthost_jaccard_summary.csv"))
write_csv(hosthost_comp, file.path(tab_dir, "table_04_hosthost_jaccard_tests.csv"))
write_csv(bind_rows(gaba_cohost, glut_cohost), file.path(tab_dir, "table_05_cohost_enrichment.csv"))
write_csv(motifs, file.path(tab_dir, "table_06_fullhost_motif_enrichment.csv"))

figure_manifest <- tibble(
  figure = c(
    "fig01_same_subclass_duplication_core",
    "fig02_hosthost_jaccard_vs_background",
    "fig03_cohost_module_enrichment",
    "fig04_representative_fullhost_motifs"
  ),
  meaning = c(
    "核心结论：多full-host neuron中，同一nonneuron subclass重复为0，而partial背景中重复很常见",
    "多宿主neuron周围，不同nonneuron host彼此重叠显著高于普通NN-NN背景",
    "co-host组合并非随机，血管相关模块最稳定",
    "neuron-nonneuron full-host关系具有强烈subclass specificity"
  )
)
write_csv(figure_manifest, file.path(out_dir, "figure_manifest.csv"))

# =========================
# Figure 1: Core teacher result
# =========================
fig1_df <- overall_dup %>%
  filter(context %in% c("Full host", "Partial background")) %>%
  mutate(context = factor(context, levels = c("Full host", "Partial background")),
         label = paste0(scales::comma(n_same), "/", scales::comma(n_total)))

bracket_y <- max(fig1_df$rate) + 0.07
fig01 <- ggplot(fig1_df, aes(x = context, y = rate, fill = context)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = label), vjust = -0.22, size = 5.0) +
  geom_text(aes(y = pmax(rate * 0.45, 0.02), label = scales::percent(rate, accuracy = 0.1)),
            colour = "white", fontface = "bold", size = 4.4) +
  geom_segment(aes(x = 1, xend = 2, y = bracket_y, yend = bracket_y), inherit.aes = FALSE, linewidth = 0.6) +
  geom_segment(aes(x = 1, xend = 1, y = bracket_y, yend = bracket_y - 0.02), inherit.aes = FALSE, linewidth = 0.6) +
  geom_segment(aes(x = 2, xend = 2, y = bracket_y, yend = bracket_y - 0.02), inherit.aes = FALSE, linewidth = 0.6) +
  annotate("text", x = 1.5, y = bracket_y + 0.02,
           label = paste0("Fisher P = ", format(dup_ft$p, scientific = TRUE, digits = 2)),
           size = 4.4) +
  scale_fill_manual(values = c("Full host" = "#C44E52", "Partial background" = "#7E7E7E")) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, max(fig1_df$rate) + 0.14),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Repeated full hosts are never from the same nonneuron subclass",
    subtitle = "Among neurons with at least two complete nonneuron hosts, same-subclass duplication disappears, but common in the partial multi-partner background.",
    x = NULL,
    y = "Fraction with repeated same-subclass partners"
  ) +
  base_theme +
  theme(legend.position = "none")
save_pub(fig01, "fig01_same_subclass_duplication_core", 12.0, 5.6)

# =========================
# Figure 2: Host-host overlap structure
# =========================
fig2_df <- bind_rows(
  tibble(group = "All NN-NN background", neuron_type = "Background", jaccard = nnn$jaccard),
  tibble(group = "Multi-host around GABA", neuron_type = "Gaba", jaccard = hosthost$jaccard[hosthost$neuron_type == "Gaba"]),
  tibble(group = "Multi-host around Glut", neuron_type = "Glut", jaccard = hosthost$jaccard[hosthost$neuron_type == "Glut"])
) %>%
  filter(!is.na(jaccard)) %>%
  mutate(group = factor(group, levels = c("All NN-NN background", "Multi-host around GABA", "Multi-host around Glut")))

y_max_fig02 <- quantile(fig2_df$jaccard, 0.995, na.rm = TRUE)
med_labs <- fig2_df %>%
  group_by(group) %>%
  summarise(med = median(jaccard), .groups = "drop") %>%
  mutate(x = group,
         y = y_max_fig02 * 0.93,
         label = paste0("median=", sprintf("%.3f", med)))

fig02 <- ggplot(fig2_df, aes(x = group, y = jaccard, fill = group)) +
  geom_violin(trim = TRUE, alpha = 0.75, colour = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white", colour = "#222222", linewidth = 0.6) +
  geom_text(data = med_labs, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 4.2) +
  scale_fill_manual(values = c("All NN-NN background" = "#B9B9B9", "Multi-host around GABA" = "#CC6B6B", "Multi-host around Glut" = "#5B84B1")) +
  coord_cartesian(ylim = c(0, y_max_fig02)) +
  labs(
    title = "Different full hosts around the same neuron strongly overlap one another",
    subtitle = paste0(
      "Compared with the global NN-NN background, host-host Jaccard is strongly elevated around multi-host neurons. ",
      "GABA: P=", format(hosthost_comp$p[hosthost_comp$comparison == "GABA multi-host vs background"], scientific = TRUE, digits = 2),
      "; Glut: P=", format(hosthost_comp$p[hosthost_comp$comparison == "Glut multi-host vs background"], scientific = TRUE, digits = 2)
    ),
    x = NULL,
    y = "Nonneuron–nonneuron Jaccard"
  ) +
  base_theme +
  theme(legend.position = "none", axis.text.x = element_text(angle = 12, hjust = 1))
save_pub(fig02, "fig02_hosthost_jaccard_vs_background", 12.4, 5.8)

# =========================
# Figure 3: Significant co-host modules
# =========================
fig3_df <- cohost_sig %>%
  mutate(pair = fct_reorder(pair, log2OR))

fig03 <- ggplot(fig3_df, aes(x = log2OR, y = pair, size = both, colour = -log10(q))) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5, colour = "#808080") +
  geom_point(alpha = 0.95) +
  scale_colour_gradient(low = "#7EBCD2", high = "#C22E2E") +
  facet_wrap(~neuron_type, scales = "free_y", ncol = 1) +
  labs(
    title = "Co-host combinations are highly structured rather than random",
    subtitle = "The most stable enrichments are vascular-associated modules, especially SMC + VLMC and Endo + Peri.",
    x = "log2(odds ratio) for co-host enrichment among multi-host neurons",
    y = NULL,
    size = "# neurons",
    colour = "-log10(FDR)"
  ) +
  coord_cartesian(clip = "off") +
  base_theme +
  theme(
    plot.margin = margin(15, 28, 12, 12)
  )

save_pub(fig03, "fig03_cohost_module_enrichment", 12.6, 7.2)

# =========================
# Figure 4: Representative motif enrichments
# =========================
fig4_df <- motif_top %>%
  mutate(motif = fct_reorder(motif, log2OR))

fig04 <- ggplot(fig4_df, aes(x = log2OR, y = motif, size = n_total, colour = -log10(q))) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5, colour = "#808080") +
  geom_point(alpha = 0.95) +
  scale_colour_gradient(low = "#86BBD8", high = "#AE2012") +
  facet_wrap(~neuron_type, scales = "free_y", ncol = 1) +
  labs(
    title = "Full-host motifs are strongly subclass-specific",
    subtitle = "Representative enriched motifs show that neuron–nonneuron organization cannot be reduced to a simple GABA/Glut binary split.",
    x = "log2(odds ratio) for Neuron-in-Nonneuron enrichment",
    y = NULL,
    size = "# pairs",
    colour = "-log10(FDR)"
  ) +
  base_theme
save_pub(fig04, "fig04_representative_fullhost_motifs", 12.8, 7.4)

message("R_output_05 complete: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
