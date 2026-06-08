#!/usr/bin/env Rscript

# ============================================================
# Lamp5 inhibitory neuron layer-enrichment and spatial profile
# Core inputs follow Total_protocol_subclass_over3.R:
#   - E:/zaw/2511/mouseMerfish_zhuang_subclass/ws0.4_ss0.02/mouse_subclass_cluster_total_over3.txt
#   - E:/zaw/2511/mouseMerfish_zhuang_subclass/neocortex_new/*.txt
#   - E:/zaw/2511/mouseMerfish_zhuang_subclass/Merfish_mouse_neocortex_layer_region.txt
#
# Outputs:
#   E:/zaw/2511/mouseMerfish_zhuang_subclass/ws0.4_ss0.02/Lamp5_layer_enrichment/
# ============================================================

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# -----------------------------
# 0. Parameters
# -----------------------------
base_dir <- Sys.getenv("MERFISH_BASE_DIR", "E:/zaw/2511/mouseMerfish_zhuang_subclass")
window_size <- Sys.getenv("WINDOW_SIZE", "0.4")
step_size <- Sys.getenv("STEP_SIZE", "0.02")
subdir_tag <- paste0("ws", window_size, "_ss", step_size)
combo_root <- file.path(base_dir, subdir_tag)

out_dir <- file.path(combo_root, "Lamp5_layer_enrichment")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# Optional quick test, e.g. Sys.setenv(LAMP5_MAX_FILES=10) before source().
max_files_env <- Sys.getenv("LAMP5_MAX_FILES", "")
max_files <- suppressWarnings(as.integer(max_files_env))
if (is.na(max_files)) max_files <- Inf

top_spatial_slices <- 12L
layer_levels <- c("1", "2/3", "4", "5", "6a", "6b")
lamp5_pattern <- "\\bLamp5\\b"

cluster_file_candidates <- c(
  file.path(combo_root, "mouse_subclass_cluster_total_over3.txt"),
  file.path(combo_root, "mouse_subclass_cluster_total.txt"),
  file.path(base_dir, "mouse_subclass_cluster_total.txt")
)
cluster_file <- cluster_file_candidates[file.exists(cluster_file_candidates)][1]
if (is.na(cluster_file)) {
  stop("Cannot find mouse_subclass_cluster_total_over3.txt or fallback cluster_total file.", call. = FALSE)
}

neocortex_dir <- file.path(base_dir, "neocortex_new")
layer_map_file <- file.path(base_dir, "Merfish_mouse_neocortex_layer_region.txt")
cell_type_file <- file.path(base_dir, "Merfish_brain_cell_type.txt")

# -----------------------------
# 1. Helpers
# -----------------------------
clean_layer <- function(x) {
  y <- tolower(trimws(as.character(x)))
  y <- gsub("^layer\\s*", "", y)
  y <- gsub("\\s+", "", y)
  y <- gsub("^2[._-]?3$", "2/3", y)
  y
}

clean_subclass <- function(x) {
  y <- trimws(as.character(x))
  y <- gsub("\\s+", " ", y)
  y
}

safe_div <- function(a, b) {
  ifelse(is.na(b) | b <= 0, NA_real_, a / b)
}

safe_log10 <- function(x) {
  log10(pmax(as.numeric(x), 1))
}

parse_ids <- function(x) {
  if (length(x) == 0 || is.na(x) || !nzchar(x)) return(character())
  regmatches(as.character(x), gregexpr("[0-9]+", as.character(x), perl = TRUE))[[1]]
}

save_plot <- function(p, stem, width, height) {
  ggsave(file.path(fig_dir, paste0(stem, ".png")), p, width = width, height = height,
         dpi = 600, bg = "white", limitsize = FALSE)
  ggsave(file.path(fig_dir, paste0(stem, ".pdf")), p, width = width, height = height,
         device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE)
}

theme_lamp5 <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 4, hjust = 0),
      plot.subtitle = element_text(size = base_size, color = "#444444", hjust = 0),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "#111111"),
      axis.line = element_line(linewidth = 0.55, color = "#222222"),
      axis.ticks = element_line(linewidth = 0.45, color = "#222222"),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "#F3F3F3", color = NA),
      strip.text = element_text(face = "bold"),
      panel.grid.major.y = element_line(color = "#E8E8E8", linewidth = 0.35),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(10, 16, 10, 12)
    )
}

message_step <- function(...) {
  message(format(Sys.time(), "%H:%M:%S"), " | ", ...)
}

# -----------------------------
# 2. Cluster-level Lamp5 enrichment
# -----------------------------
message_step("Reading cluster table: ", cluster_file)
clu <- fread(cluster_file, sep = "\t", header = TRUE, data.table = TRUE, showProgress = TRUE)

required_cluster_cols <- c(
  "slide", "layer", "cell_Neruon_type", "subclass", "total_cell_num",
  "enrich_subclass_cell_ids_num", "enrich_subclass_cell_ids",
  "Glut_Neruon_cell_ids_num", "GABA_Neruon_cell_ids_num"
)
miss <- setdiff(required_cluster_cols, names(clu))
if (length(miss)) stop("Cluster table missing columns: ", paste(miss, collapse = ", "), call. = FALSE)

num_cols <- c("total_cell_num", "enrich_subclass_cell_ids_num",
              "Glut_Neruon_cell_ids_num", "GABA_Neruon_cell_ids_num")
for (nm in intersect(num_cols, names(clu))) {
  clu[, (nm) := as.numeric(get(nm))]
}

clu[, cluster_id := .I]
clu[, layer_f := clean_layer(layer)]
clu <- clu[layer_f %in% layer_levels]
clu[, layer_f := factor(layer_f, levels = layer_levels)]
clu[, subclass_clean := clean_subclass(subclass)]
clu[, is_lamp5 := grepl(lamp5_pattern, subclass_clean, ignore.case = TRUE)]
clu[, is_gaba_cluster := cell_Neruon_type %chin% c("Gaba", "GABA") |
      grepl("Gaba", subclass_clean, ignore.case = TRUE)]
clu[, ei_ratio := safe_div(Glut_Neruon_cell_ids_num, GABA_Neruon_cell_ids_num)]
clu[, lamp5_cluster_label := ifelse(is_lamp5, "Lamp5 enriched cluster", "Other cluster")]

lamp5_clu <- clu[is_lamp5 == TRUE]
if (nrow(lamp5_clu) == 0) {
  stop("No Lamp5 enriched clusters found in cluster table.", call. = FALSE)
}

cluster_layer_summary <- clu[, .(
  n_clusters = .N,
  n_gaba_clusters = sum(is_gaba_cluster, na.rm = TRUE),
  n_lamp5_clusters = sum(is_lamp5, na.rm = TRUE),
  total_cluster_cells = sum(total_cell_num, na.rm = TRUE),
  lamp5_region_cells = sum(total_cell_num[is_lamp5], na.rm = TRUE),
  lamp5_enriched_cells = sum(enrich_subclass_cell_ids_num[is_lamp5], na.rm = TRUE),
  median_lamp5_cluster_size = median(total_cell_num[is_lamp5], na.rm = TRUE),
  median_lamp5_enriched_cells = median(enrich_subclass_cell_ids_num[is_lamp5], na.rm = TRUE),
  median_lamp5_ei_ratio = median(ei_ratio[is_lamp5], na.rm = TRUE)
), by = layer_f][order(layer_f)]
cluster_layer_summary[, lamp5_pct_of_gaba_clusters := safe_div(n_lamp5_clusters, n_gaba_clusters)]
cluster_layer_summary[, lamp5_pct_of_all_clusters := safe_div(n_lamp5_clusters, n_clusters)]

lamp5_cluster_detail <- copy(lamp5_clu)[
  order(layer_f, slide, subclass_clean, -enrich_subclass_cell_ids_num),
  .(cluster_id, slide, layer = layer_f, region, merge_regions, cell_Neruon_type,
    subclass = subclass_clean, total_cell_num, enrich_subclass_cell_ids_num,
    Glut_Neruon_cell_ids_num, GABA_Neruon_cell_ids_num, ei_ratio)
]

lamp5_subclass_layer <- lamp5_clu[, .(
  n_clusters = .N,
  region_cells = sum(total_cell_num, na.rm = TRUE),
  enriched_cells = sum(enrich_subclass_cell_ids_num, na.rm = TRUE),
  median_ei_ratio = median(ei_ratio, na.rm = TRUE)
), by = .(layer_f, subclass_clean)][order(layer_f, -n_clusters, subclass_clean)]

fwrite(cluster_layer_summary, file.path(tab_dir, "Lamp5_cluster_layer_summary.tsv"),
       sep = "\t", quote = FALSE)
fwrite(lamp5_cluster_detail, file.path(tab_dir, "Lamp5_cluster_detail.tsv"),
       sep = "\t", quote = FALSE)
fwrite(lamp5_subclass_layer, file.path(tab_dir, "Lamp5_subclass_layer_summary.tsv"),
       sep = "\t", quote = FALSE)

# Unique enriched Lamp5 cells represented by merged enriched regions.
message_step("Unnesting Lamp5 enriched cell ids from cluster table")
lamp5_cluster_cells <- lamp5_clu[, {
  ids <- parse_ids(enrich_subclass_cell_ids)
  if (length(ids) == 0) {
    .(cell_label = character())
  } else {
    .(cell_label = ids)
  }
}, by = .(cluster_id, slide, layer_f, subclass_clean)]
lamp5_cluster_cells <- unique(lamp5_cluster_cells[nzchar(cell_label)])

lamp5_cluster_cell_summary <- lamp5_cluster_cells[, .(
  n_unique_lamp5_cells_in_enriched_regions = uniqueN(cell_label),
  n_lamp5_clusters_containing_unique_cells = uniqueN(cluster_id)
), by = layer_f][order(layer_f)]

fwrite(lamp5_cluster_cells, file.path(tab_dir, "Lamp5_enriched_region_unique_cells.tsv"),
       sep = "\t", quote = FALSE)
fwrite(lamp5_cluster_cell_summary, file.path(tab_dir, "Lamp5_enriched_region_unique_cell_summary.tsv"),
       sep = "\t", quote = FALSE)

# -----------------------------
# 3. Raw cell-level Lamp5 distribution
# -----------------------------
message_step("Reading raw neocortex cell files")
if (!dir.exists(neocortex_dir)) stop("Cannot find neocortex_new directory: ", neocortex_dir, call. = FALSE)
if (!file.exists(layer_map_file)) stop("Cannot find layer map file: ", layer_map_file, call. = FALSE)

layer_map <- fread(layer_map_file, sep = "\t", header = TRUE, data.table = TRUE)
if (!all(c("ccf_region_name", "layer") %in% names(layer_map))) {
  stop("Layer map must contain ccf_region_name and layer columns.", call. = FALSE)
}
layer_map[, layer_f := clean_layer(layer)]
layer_map <- unique(layer_map[, .(ccf_region_name, layer_f)])

cell_type_map <- NULL
if (file.exists(cell_type_file)) {
  cell_type_map <- fread(cell_type_file, sep = "\t", header = TRUE, data.table = TRUE)
  if (!all(c("class", "cell_Neuron_type") %in% names(cell_type_map))) {
    cell_type_map <- NULL
  } else {
    cell_type_map <- unique(cell_type_map[, .(class, cell_Neuron_type)])
  }
}

read_one_slice <- function(fp) {
  hdr <- names(fread(fp, sep = "\t", nrows = 0, header = TRUE))
  need <- c("brain_section_label", "cell_label", "ccf_region_name", "x", "y", "class", "subclass")
  sel <- intersect(need, hdr)
  if (!all(c("brain_section_label", "cell_label", "ccf_region_name", "x", "y", "subclass") %in% sel)) {
    warning("Skipping file with missing required columns: ", fp)
    return(NULL)
  }
  dt <- fread(fp, sep = "\t", select = sel, header = TRUE,
              colClasses = list(character = c("cell_label")), data.table = TRUE)
  setnames(dt, old = "brain_section_label", new = "slide", skip_absent = TRUE)
  if (!"class" %in% names(dt)) dt[, class := NA_character_]
  dt[, cell_label := as.character(cell_label)]
  dt[, x := as.numeric(x)]
  dt[, y := as.numeric(y)]
  dt <- dt[!is.na(x) & !is.na(y)]
  dt <- merge(dt, layer_map, by = "ccf_region_name", all.x = FALSE, all.y = FALSE)
  if (!is.null(cell_type_map) && "class" %in% names(dt)) {
    dt <- merge(dt, cell_type_map, by = "class", all.x = TRUE)
  } else {
    dt[, cell_Neuron_type := NA_character_]
  }
  dt[, subclass_clean := clean_subclass(subclass)]
  dt[, is_lamp5_cell := grepl(lamp5_pattern, subclass_clean, ignore.case = TRUE)]
  dt[, is_gaba_cell := cell_Neuron_type %chin% c("GABA", "Gaba") |
       grepl("Gaba", subclass_clean, ignore.case = TRUE) |
       grepl("Gaba", class, ignore.case = TRUE)]
  dt[, layer_f := factor(layer_f, levels = layer_levels)]
  dt[layer_f %in% layer_levels]
}

slice_files <- list.files(neocortex_dir, pattern = "\\.txt$", full.names = TRUE)
slice_files <- sort(slice_files)
if (length(slice_files) == 0) stop("No .txt files found in neocortex_new.", call. = FALSE)
if (is.finite(max_files)) slice_files <- head(slice_files, max_files)

cell_list <- vector("list", length(slice_files))
for (i in seq_along(slice_files)) {
  if (i %% 10L == 1L || i == length(slice_files)) {
    message_step("  raw slice ", i, "/", length(slice_files), ": ", basename(slice_files[i]))
  }
  cell_list[[i]] <- read_one_slice(slice_files[i])
}
cell_dt <- rbindlist(cell_list, use.names = TRUE, fill = TRUE)
if (nrow(cell_dt) == 0) stop("No raw cells were loaded after layer mapping.", call. = FALSE)

cell_layer_summary <- cell_dt[, .(
  n_cells = .N,
  n_gaba_cells = sum(is_gaba_cell, na.rm = TRUE),
  n_lamp5_cells = sum(is_lamp5_cell, na.rm = TRUE),
  n_slides = uniqueN(slide)
), by = layer_f][order(layer_f)]
cell_layer_summary[, lamp5_pct_all_cells := safe_div(n_lamp5_cells, n_cells)]
cell_layer_summary[, lamp5_pct_gaba_cells := safe_div(n_lamp5_cells, n_gaba_cells)]

cell_slide_layer_summary <- cell_dt[, .(
  n_cells = .N,
  n_gaba_cells = sum(is_gaba_cell, na.rm = TRUE),
  n_lamp5_cells = sum(is_lamp5_cell, na.rm = TRUE),
  mean_x_lamp5 = mean(x[is_lamp5_cell], na.rm = TRUE),
  mean_y_lamp5 = mean(y[is_lamp5_cell], na.rm = TRUE),
  sd_x_lamp5 = sd(x[is_lamp5_cell], na.rm = TRUE),
  sd_y_lamp5 = sd(y[is_lamp5_cell], na.rm = TRUE)
), by = .(slide, layer_f)][order(layer_f, slide)]
cell_slide_layer_summary[, lamp5_pct_all_cells := safe_div(n_lamp5_cells, n_cells)]
cell_slide_layer_summary[, lamp5_pct_gaba_cells := safe_div(n_lamp5_cells, n_gaba_cells)]

lamp5_spatial_feature <- cell_dt[is_lamp5_cell == TRUE, .(
  n_lamp5_cells = .N,
  x_centroid = mean(x, na.rm = TRUE),
  y_centroid = mean(y, na.rm = TRUE),
  x_sd = sd(x, na.rm = TRUE),
  y_sd = sd(y, na.rm = TRUE),
  x_min = min(x, na.rm = TRUE),
  x_max = max(x, na.rm = TRUE),
  y_min = min(y, na.rm = TRUE),
  y_max = max(y, na.rm = TRUE)
), by = .(slide, layer_f)][order(layer_f, -n_lamp5_cells)]

fwrite(cell_layer_summary, file.path(tab_dir, "Lamp5_raw_cell_layer_summary.tsv"),
       sep = "\t", quote = FALSE)
fwrite(cell_slide_layer_summary, file.path(tab_dir, "Lamp5_raw_cell_slide_layer_summary.tsv"),
       sep = "\t", quote = FALSE)
fwrite(lamp5_spatial_feature, file.path(tab_dir, "Lamp5_raw_cell_spatial_features_by_slide_layer.tsv"),
       sep = "\t", quote = FALSE)

# -----------------------------
# 4. Figures
# -----------------------------
lamp5_orange <- "#E65F3C"
lamp5_deep <- "#A92320"
gaba_blue <- "#266AA6"
neutral_gray <- "#B8B8B8"

plot_cluster_count <- ggplot(cluster_layer_summary, aes(x = layer_f, y = n_lamp5_clusters)) +
  geom_col(width = 0.72, fill = lamp5_deep, color = "black", linewidth = 0.25) +
  geom_text(aes(label = n_lamp5_clusters), vjust = -0.35, size = 3.8, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Lamp5 enriched clusters across cortical layers",
    subtitle = "Clusters are from mouse_subclass_cluster_total_over3; layer 2/3 is kept as one layer.",
    x = "Layer", y = "Number of Lamp5 enriched clusters"
  ) +
  theme_lamp5(13)
save_plot(plot_cluster_count, "Fig_Lamp5_01_cluster_count_by_layer", 7.2, 4.8)

plot_cluster_pct <- ggplot(cluster_layer_summary,
                           aes(x = layer_f, y = 100 * lamp5_pct_of_gaba_clusters)) +
  geom_col(width = 0.72, fill = lamp5_orange, color = "black", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * lamp5_pct_of_gaba_clusters)),
            vjust = -0.35, size = 3.7, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(
    title = "Lamp5 share among inhibitory enriched clusters",
    subtitle = "Denominator: all GABA-class enriched clusters in the same layer.",
    x = "Layer", y = "Lamp5 / GABA clusters (%)"
  ) +
  theme_lamp5(13)
save_plot(plot_cluster_pct, "Fig_Lamp5_02_cluster_fraction_among_gaba_by_layer", 7.2, 4.8)

plot_cell_count <- ggplot(cell_layer_summary, aes(x = layer_f, y = n_lamp5_cells)) +
  geom_col(width = 0.72, fill = lamp5_deep, color = "black", linewidth = 0.25) +
  geom_text(aes(label = format(n_lamp5_cells, big.mark = ",")),
            vjust = -0.35, size = 3.6, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.13))) +
  labs(
    title = "Raw Lamp5 cell distribution across layers",
    subtitle = "Counts are computed directly from neocortex_new cells after ccf-region to layer mapping.",
    x = "Layer", y = "Number of Lamp5 cells"
  ) +
  theme_lamp5(13)
save_plot(plot_cell_count, "Fig_Lamp5_03_raw_cell_count_by_layer", 7.2, 4.8)

plot_cell_pct <- ggplot(cell_layer_summary,
                        aes(x = layer_f, y = 100 * lamp5_pct_gaba_cells)) +
  geom_col(width = 0.72, fill = gaba_blue, color = "black", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * lamp5_pct_gaba_cells)),
            vjust = -0.35, size = 3.6, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.13))) +
  labs(
    title = "Lamp5 proportion within inhibitory cells",
    subtitle = "Denominator: raw GABA cells in the same mapped cortical layer.",
    x = "Layer", y = "Lamp5 / GABA cells (%)"
  ) +
  theme_lamp5(13)
save_plot(plot_cell_pct, "Fig_Lamp5_04_raw_cell_fraction_among_gaba_by_layer", 7.2, 4.8)

plot_cluster_size <- ggplot(clu[is_gaba_cluster == TRUE],
                            aes(x = layer_f, y = safe_log10(total_cell_num), fill = lamp5_cluster_label)) +
  geom_boxplot(width = 0.66, outlier.shape = NA, alpha = 0.88, color = "black", linewidth = 0.25) +
  geom_jitter(aes(color = lamp5_cluster_label), width = 0.18, size = 0.75, alpha = 0.35, show.legend = FALSE) +
  scale_fill_manual(values = c("Lamp5 enriched cluster" = lamp5_deep, "Other cluster" = neutral_gray)) +
  scale_color_manual(values = c("Lamp5 enriched cluster" = lamp5_deep, "Other cluster" = "#777777")) +
  labs(
    title = "Lamp5 cluster size relative to other inhibitory clusters",
    subtitle = "Y axis uses log10(total cells in merged region + 1).",
    x = "Layer", y = "log10(total cells + 1)", fill = "Cluster type"
  ) +
  theme_lamp5(13)
save_plot(plot_cluster_size, "Fig_Lamp5_05_cluster_size_vs_other_gaba", 8.4, 5.2)

ei_lims <- quantile(lamp5_clu$ei_ratio, probs = c(0.02, 0.98), na.rm = TRUE)
use_ei_lims <- all(is.finite(ei_lims)) && diff(ei_lims) > 0

plot_ei <- ggplot(lamp5_clu, aes(x = layer_f, y = ei_ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#333333", linewidth = 0.45) +
  geom_boxplot(width = 0.58, fill = "#F2B299", color = "black", outlier.shape = NA, linewidth = 0.3) +
  geom_jitter(width = 0.18, size = 1.1, alpha = 0.65, color = lamp5_deep) +
  labs(
    title = "E/I ratio inside Lamp5 enriched regions",
    subtitle = "E/I ratio = GLU cell number / GABA cell number within each merged enriched region.",
    x = "Layer", y = "E/I ratio"
  ) +
  theme_lamp5(13)
if (use_ei_lims) {
  plot_ei <- plot_ei + coord_cartesian(ylim = ei_lims)
}
save_plot(plot_ei, "Fig_Lamp5_06_lamp5_region_ei_ratio_by_layer", 7.8, 5)

plot_subclass_layer <- ggplot(lamp5_subclass_layer,
                              aes(x = layer_f, y = n_clusters, fill = subclass_clean)) +
  geom_col(width = 0.74, color = "black", linewidth = 0.18) +
  labs(
    title = "Lamp5 subtype composition by layer",
    subtitle = "If multiple Lamp5-labelled subclasses exist, this separates their enriched clusters.",
    x = "Layer", y = "Number of enriched clusters", fill = "Subclass"
  ) +
  theme_lamp5(12) +
  theme(legend.position = "right")
save_plot(plot_subclass_layer, "Fig_Lamp5_07_subclass_composition_by_layer", 8.8, 5.2)

# Top raw slices by Lamp5 cell count: gray all cells + red Lamp5 cells.
top_slice_tbl <- cell_dt[is_lamp5_cell == TRUE, .N, by = slide][order(-N)]
top_slices <- top_slice_tbl$slide[seq_len(min(top_spatial_slices, nrow(top_slice_tbl)))]
spatial_dt <- cell_dt[slide %chin% top_slices]
spatial_dt[, slide := factor(slide, levels = top_slices)]

plot_spatial <- ggplot() +
  geom_point(
    data = spatial_dt,
    aes(x = x, y = y),
    color = "#CFCFCF", size = 0.12, alpha = 0.35
  ) +
  geom_point(
    data = spatial_dt[is_lamp5_cell == TRUE],
    aes(x = x, y = y),
    color = lamp5_deep, size = 0.28, alpha = 0.95
  ) +
  facet_wrap(~ slide, ncol = 4) +
  coord_equal() +
  scale_y_reverse() +
  labs(
    title = "Spatial distribution of raw Lamp5 cells in representative slices",
    subtitle = "Top slices are selected by raw Lamp5 cell count; gray = all mapped neocortex cells, red = Lamp5.",
    x = NULL, y = NULL
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0),
    plot.subtitle = element_text(size = 11, color = "#444444", hjust = 0),
    strip.text = element_text(face = "bold", size = 10),
    panel.border = element_rect(fill = NA, color = "#222222", linewidth = 0.35),
    plot.margin = margin(10, 10, 10, 10)
  )
save_plot(plot_spatial, "Fig_Lamp5_08_raw_spatial_top_slices", 12, 9)

# Layer-resolved per-slice variability.
plot_slice_layer <- ggplot(cell_slide_layer_summary[n_lamp5_cells > 0],
                           aes(x = layer_f, y = 100 * lamp5_pct_gaba_cells)) +
  geom_violin(fill = "#F8D9CF", color = "#8F1E1B", linewidth = 0.35, trim = TRUE) +
  geom_boxplot(width = 0.18, fill = "white", outlier.shape = NA, linewidth = 0.35) +
  geom_jitter(width = 0.13, size = 0.75, alpha = 0.45, color = lamp5_deep) +
  labs(
    title = "Slice-to-slice variability of Lamp5 fraction",
    subtitle = "Each dot is one slice-layer pair; denominator is raw GABA cells in the same slice-layer.",
    x = "Layer", y = "Lamp5 / GABA cells (%)"
  ) +
  theme_lamp5(13)
save_plot(plot_slice_layer, "Fig_Lamp5_09_slice_layer_variability", 7.8, 5)

# Optional compact summary is intentionally skipped by default.
# Some local ggplot2/patchwork combinations fail during guide collection with:
#   Error in Ops.data.frame(guide_loc, panel_loc)
# The individual publication panels above are the authoritative outputs.
if (identical(Sys.getenv("LAMP5_MAKE_PATCHWORK", ""), "1") &&
    requireNamespace("patchwork", quietly = TRUE)) {
  tryCatch({
    combined <- patchwork::wrap_plots(
      plot_cluster_count + theme(legend.position = "none"),
      plot_cluster_pct + theme(legend.position = "none"),
      plot_cell_count + theme(legend.position = "none"),
      plot_cell_pct + theme(legend.position = "none"),
      plot_cluster_size + theme(legend.position = "none"),
      plot_ei + theme(legend.position = "none"),
      ncol = 2
    )
    save_plot(combined, "Fig_Lamp5_00_summary_multipanel", 14, 14)
  }, error = function(e) {
    message_step("Skip optional patchwork summary because patchwork failed: ", conditionMessage(e))
  })
}

# -----------------------------
# 5. Run metadata
# -----------------------------
run_info <- data.table(
  item = c(
    "base_dir", "combo_root", "cluster_file", "neocortex_dir",
    "n_cluster_rows_used", "n_lamp5_clusters", "n_raw_cells_used",
    "n_raw_lamp5_cells", "n_raw_files_read"
  ),
  value = c(
    base_dir, combo_root, cluster_file, neocortex_dir,
    as.character(nrow(clu)), as.character(nrow(lamp5_clu)),
    as.character(nrow(cell_dt)), as.character(sum(cell_dt$is_lamp5_cell, na.rm = TRUE)),
    as.character(length(slice_files))
  )
)
fwrite(run_info, file.path(tab_dir, "Lamp5_run_info.tsv"), sep = "\t", quote = FALSE)

message_step("Done.")
message_step("Figures: ", fig_dir)
message_step("Tables:  ", tab_dir)
