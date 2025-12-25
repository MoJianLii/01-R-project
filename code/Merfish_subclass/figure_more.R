rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(cowplot)
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
  library(ggforce)
  library(ggpubr)
  library(FSA)
  library(grid)
})

root_dir <- "E:/zaw/2511/mouseMerfish_zhuang_subclass"
fig_dir  <- file.path(root_dir, "figure_more")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

ws_ss_dirs <- c(
  "ws0.3_ss0.1",
  "ws0.3_ss0.02",
  "ws0.3_ss0.05",
  "ws0.4_ss0.1",
  "ws0.4_ss0.02",
  "ws0.4_ss0.05"
)

parse_ws_ss <- function(dir_name) {
  ws <- as.numeric(str_match(dir_name, "ws([0-9.]+)")[, 2])
  ss <- as.numeric(str_match(dir_name, "ss([0-9.]+)")[, 2])
  list(ws = ws, ss = ss)
}

find_subclass_col <- function(dt) {
  cand <- c("subclasstype", "subclass_type", "Subclass", "subclass", "sub_class")
  found <- intersect(cand, names(dt))
  if (length(found) == 0) {
    stop("no subclass")
  }
  found[1]
}


Table1_all <- rbindlist(
  lapply(ws_ss_dirs, function(d) {
    info <- parse_ws_ss(d)
    f    <- file.path(root_dir, d, "mouse_table", "table1", "Table_1_total_cell.txt")
    if (!file.exists(f)) stop("no file：", f)
    dt <- fread(f)
    dt[, window_size := info$ws]
    dt[, step_size   := info$ss]
    dt[, ws_ss       := d]
    dt[]
  }),
  use.names = TRUE, fill = TRUE
)

subcol <- find_subclass_col(Table1_all)
Table1_all[, subclass_raw   := get(subcol)]
Table1_all[, subclass_lower := tolower(subclass_raw)]

Table1_all[, subclass_group := "Other"]
Table1_all[grepl("pvalb|parvalb|parvalbumin", subclass_lower),
           subclass_group := "PVALB"]
Table1_all[grepl("som|sst|somatostatin",      subclass_lower),
           subclass_group := "SOM"]
Table1_all[grepl("\\bvip\\b",                 subclass_lower),
           subclass_group := "VIP"]

Table1_all[, subclass_group := factor(
  subclass_group,
  levels = c("PVALB", "SOM", "VIP", "Other")
)]

Table1_all[, ws_ss_label := sprintf("ws=%.2f\nss=%.3f", window_size, step_size)]
Table1_all[, ws_ss_label := factor(
  ws_ss_label,
  levels = unique(ws_ss_label[order(window_size, step_size)])
)]

psv_cols <- c(
  PVALB = "#4C72B0",
  SOM   = "#55A868",
  VIP   = "#C44E52"
)

Table1_all[ , .N, by = .(ws_ss, subclass_group)][order(ws_ss, subclass_group)]
psv_all <- Table1_all[subclass_group %in% c("PVALB", "SOM", "VIP")]

scale_psv <- psv_all[, .(
  n_hotspots   = .N,
  sum_cells    = sum(total_cell_num),
  mean_cells   = mean(total_cell_num),
  median_cells = median(total_cell_num)
), by = .(window_size, step_size, subclass_group)]

scale_psv[, step_factor := factor(step_size,
                                  levels = sort(unique(step_size)))]

step_levels <- levels(scale_psv$step_factor)
step_cols   <- setNames(
  brewer.pal(max(3, length(step_levels)), "Set1")[seq_along(step_levels)],
  step_levels
)

plot_psv_metric <- function(subclass_name,
                            metric = c("n_hotspots", "sum_cells")) {
  metric <- match.arg(metric)
  dt     <- scale_psv[subclass_group == subclass_name]
  
  y_lab <- if (metric == "n_hotspots") "Number of clusters" else "Total cells in hotspots"
  
  title_lab <- if (metric == "n_hotspots")
    sprintf("%s: Number of clusters", subclass_name)
  else
    sprintf("%s: Total cells in clusters", subclass_name)
  
  ggplot(dt,
         aes(x = window_size,
             y = .data[[metric]],
             colour = step_factor,
             group  = step_factor)) +
    geom_line(size = 1) +
    geom_point(size = 2.5) +
    scale_colour_manual(values = step_cols, name = "step_size") +
    scale_x_continuous(breaks = sort(unique(scale_psv$window_size))) +
    labs(x = "window_size", y = y_lab, title = title_lab) +
    theme_classic(base_size = 13) +
    theme(
      axis.text       = element_text(colour = "black"),
      plot.title      = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
}

p_PVALB_hot  <- plot_psv_metric("PVALB", "n_hotspots")
p_PVALB_cell <- plot_psv_metric("PVALB", "sum_cells")
p_SOM_hot    <- plot_psv_metric("SOM",   "n_hotspots")
p_SOM_cell   <- plot_psv_metric("SOM",   "sum_cells")
p_VIP_hot    <- plot_psv_metric("VIP",   "n_hotspots")
p_VIP_cell   <- plot_psv_metric("VIP",   "sum_cells")

legend_step <- cowplot::get_legend(
  p_PVALB_hot +
    theme(
      legend.position = "right",
      legend.title    = element_text(size = 11),
      legend.text     = element_text(size = 10)
    )
)

plots_grid <- cowplot::plot_grid(
  p_PVALB_hot  + theme(legend.position = "none"),
  p_PVALB_cell + theme(legend.position = "none"),
  p_SOM_hot    + theme(legend.position = "none"),
  p_SOM_cell   + theme(legend.position = "none"),
  p_VIP_hot    + theme(legend.position = "none"),
  p_VIP_cell   + theme(legend.position = "none"),
  nrow   = 3, ncol = 2,
  labels = c("A", "B", "C", "D", "E", "F"),
  label_size = 14,
  align  = "hv"
)

p2_combined <- cowplot::plot_grid(
  plots_grid,
  legend_step,
  ncol = 2,
  rel_widths = c(0.88, 0.12)
)

p2_combined

fig2_png <- file.path(fig_dir, "02_Fig_PSV_hotspot_scale_6panel_with_legend_right.png")
if (!file.exists(fig2_png)) {
  ggsave(fig2_png,
         p2_combined, width = 14, height = 9, dpi = 400)
}

main_ws <- 0.4
main_ss <- 0.02

main_tag     <- sprintf("ws%.2f_ss%.3f", main_ws, main_ss)
fig_dir_main <- file.path(fig_dir, main_tag)
dir.create(fig_dir_main, showWarnings = FALSE, recursive = TRUE)

dt_main <- Table1_all[
  abs(window_size - main_ws) < 1e-6 &
    abs(step_size   - main_ss) < 1e-6
]

if (!"layer" %in% names(dt_main)) {
  stop("no layer")
}

dt_main_psv <- dt_main[subclass_group %in% c("PVALB", "SOM", "VIP")]
layer_psv <- dt_main_psv[, .N, by = .(layer, subclass_group)]
layer_psv[, layer_total := sum(N), by = layer]
layer_psv[, frac        := N / layer_total]

if (!is.factor(layer_psv$layer)) {
  layer_psv[, layer := factor(layer, levels = sort(unique(layer)))]
}

p1_counts <- ggplot(layer_psv,
                    aes(x = layer, y = N, fill = subclass_group)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = psv_cols, name = "Subclass") +
  labs(
    x = "Layer",
    y = "Cell count in clusters",
    title = sprintf("PVALB/SOM/VIP cluster counts by layer (ws=%.2f, ss=%.3f)",
                    main_ws, main_ss)
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text  = element_text(colour = "black"),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p1_frac <- ggplot(layer_psv,
                  aes(x = layer, y = frac, fill = subclass_group)) +
  geom_col(position = "fill", width = 0.7,
           colour = "grey20", size = 0.2) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = psv_cols, name = "Subclass") +
  labs(
    x = "Layer",
    y = "Proportion in hotspots",
    title = sprintf("Composition of cluster cells (ws=%.2f, ss=%.3f)",
                    main_ws, main_ss)
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text  = element_text(colour = "black"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p1_combined <- plot_grid(p1_counts, p1_frac,
                         nrow = 1, rel_widths = c(1.2, 1))
p1_combined

ggsave(file.path(fig_dir_main, "01_Fig_PSV_layer_distribution_mainScale.png"),
       p1_combined, width = 15, height = 5, dpi = 400)


psv_main <- dt_main_psv
print(table(psv_main$subclass_group))

psv_main[, slide_suffix := tstrsplit(slide, "-", keep = 2)]
psv_main[, c("mouse_id", "section_raw") :=
           tstrsplit(slide_suffix, "\\.")]

psv_main[, section_index := as.numeric(section_raw)]
psv_main[is.na(section_index),
         section_index := frank(section_raw, ties.method = "dense"),
         by = mouse_id]

psv_main[, layer_f := factor(layer, levels = sort(unique(layer)))]


hot_layer_slide <- psv_main[
  ,
  .(
    n_hotspots = .N,
    sum_cells  = sum(total_cell_num)
  ),
  by = .(mouse_id, subclass_group, layer_f, section_index)
]

mouse_ids <- sort(unique(hot_layer_slide$mouse_id))
print(mouse_ids)


plot_mouse_heatmap <- function(dt_mouse) {
  mid <- unique(dt_mouse$mouse_id)
  
  ggplot(dt_mouse,
         aes(x = section_index, y = layer_f, fill = n_hotspots)) +
    geom_tile(colour = "grey90", size = 0.25) +
    scale_fill_gradientn(
      colours = colorRampPalette(c("#f7fbff", "#6baed6", "#08306b"))(100),
      name    = "Number of\nclusters"
    ) +
    facet_wrap(~ subclass_group, nrow = 3, strip.position = "right") +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
    labs(
      x = "Section index within mouse (A–P axis)",
      y = "Layer",
      title = paste0("Mouse ", mid)
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x    = element_text(colour = "black"),
      axis.text.y    = element_text(colour = "black"),
      strip.text     = element_text(face   = "bold"),
      panel.grid     = element_blank(),
      legend.position = "right",
      plot.title     = element_text(hjust  = 0.5, face = "bold")
    )
}

plots_mouse <- lapply(mouse_ids, function(mid) {
  plot_mouse_heatmap(hot_layer_slide[mouse_id == mid])
})

for (i in seq_along(mouse_ids)) {
  ggsave(
    filename = file.path(
      fig_dir_main,
      sprintf("03A_PSV_layer_slide_heatmap_mouse%s.png", mouse_ids[i])
    ),
    plot  = plots_mouse[[i]],
    width = 7.5, height = 5.5, dpi = 400
  )
}

p_layer_slide_4mice <- cowplot::plot_grid(
  plotlist   = plots_mouse,
  ncol       = 2,
  nrow       = 2,
  labels     = paste0("Mouse ", mouse_ids),
  label_size = 14
)

p_layer_slide_4mice

ggsave(
  filename = file.path(fig_dir_main,
                       "03A_PSV_layer_slide_heatmap_byMouse_2x2.png"),
  plot  = p_layer_slide_4mice,
  width = 12, height = 9, dpi = 400
)

split_region_smart <- function(reg_str) {
  if (is.null(reg_str) || length(reg_str) == 0L) return(character(0L))
  reg_str <- reg_str[1]
  
  toks <- strsplit(reg_str, ",")[[1]]
  toks <- trimws(toks)
  toks <- toks[toks != ""]
  if (!length(toks)) return(character(0L))
  
  layer_idx <- grep("^layer\\b", toks)
  if (!length(layer_idx)) return(paste(toks, collapse = ", "))
  
  res   <- character(0L)
  start <- 1L
  for (idx in layer_idx) {
    res   <- c(res, paste(toks[start:idx], collapse = ", "))
    start <- idx + 1L
  }
  res
}

region_long_main <- psv_main[
  ,
  {
    segs_list <- lapply(region, split_region_smart)
    .(region = unlist(segs_list, use.names = FALSE))
  },
  by = .(subclass_group)
]

region_long_main <- region_long_main[region != ""]

hot_region <- region_long_main[
  ,
  .(n_hotspots = .N),
  by = .(subclass_group, region)
]

hot_region[, rank := frank(-n_hotspots, ties.method = "min"),
           by = subclass_group]
hot_region_top <- hot_region[rank <= 10]

hot_region_top[, subclass_group :=
                 factor(subclass_group, levels = c("PVALB", "SOM", "VIP"))]

hot_region_top[, region_id := paste(subclass_group, region, sep = "__")]
setorder(hot_region_top, subclass_group, n_hotspots)
hot_region_top[, region_id := factor(region_id, levels = unique(region_id))]

p_region_bar <- ggplot(hot_region_top,
                       aes(x = region_id, y = n_hotspots,
                           fill = subclass_group)) +
  geom_col(width = 0.7, colour = "grey20", size = 0.25) +
  scale_fill_manual(values = psv_cols, name = "Subclass") +
  coord_flip() +
  facet_wrap(~ subclass_group, nrow = 1, scales = "free_y") +
  scale_x_discrete(labels = function(x) sub("^[^_]+__", "", x)) +
  labs(
    x     = "",
    y     = "Number of clusters",
    title = "Top anatomical regions containing clusters\n(ws = 0.40, step = 0.02)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text       = element_text(colour = "black"),
    strip.text      = element_text(face   = "bold"),
    legend.position = "none",
    plot.title      = element_text(hjust  = 0.5, face = "bold")
  )

p_region_bar

ggsave(file.path(fig_dir_main,
                 "03B_PSV_hotspot_region_bar_sorted.png"),
       p_region_bar, width = 20, height = 5.5, dpi = 400)

region_long_main[, region_no_layer :=
                   sub(",\\s*layer\\s+[^,]+$", "", region)]
region_long_main <- region_long_main[region_no_layer != ""]

hot_region_nl <- region_long_main[
  ,
  .(n_hotspots = .N),
  by = .(subclass_group, region_no_layer)
]

hot_region_nl[, rank := frank(-n_hotspots, ties.method = "min"),
              by = subclass_group]
hot_region_nl_top <- hot_region_nl[rank <= 10]

hot_region_nl_top[, subclass_group :=
                    factor(subclass_group, levels = c("PVALB", "SOM", "VIP"))]

hot_region_nl_top[, region_id :=
                    paste(subclass_group, region_no_layer, sep = "__")]
setorder(hot_region_nl_top, subclass_group, n_hotspots)
hot_region_nl_top[, region_id :=
                    factor(region_id, levels = unique(region_id))]

p_region_bar_nolayer <- ggplot(hot_region_nl_top,
                               aes(x = region_id, y = n_hotspots,
                                   fill = subclass_group)) +
  geom_col(width = 0.7, colour = "grey20", size = 0.25) +
  scale_fill_manual(values = psv_cols, name = "Subclass") +
  coord_flip() +
  facet_wrap(~ subclass_group, nrow = 1, scales = "free_y") +
  scale_x_discrete(labels = function(x) sub("^[^_]+__", "", x)) +
  labs(
    x     = "",
    y     = "Number of clusters",
    title = "Top anatomical regions (layers collapsed) containing clusters\n(ws = 0.40, step = 0.02)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text       = element_text(colour = "black"),
    strip.text      = element_text(face = "bold"),
    legend.position = "none",
    plot.title      = element_text(hjust  = 0.5, face  = "bold")
  )

p_region_bar_nolayer

ggsave(file.path(fig_dir_main,
                 "03C_PSV_hotspot_region_noLayer_bar_sorted.png"),
       p_region_bar_nolayer, width = 20, height = 5.5, dpi = 400)


slice_to_plot <- psv_main[, .N, by = slide][order(-N)]$slide[1]
message("Slice chosen for visualization: ", slice_to_plot)


psv_region_long <- psv_main[
  ,
  {
    segs_list <- lapply(region, split_region_smart)
    .(region = unlist(segs_list, use.names = FALSE))
  },
  by = .(subclass_group)
]

psv_region_long <- psv_region_long[region != ""]
psv_region_long[, region_no_layer :=
                  sub(",\\s*layer\\s+[^,]+$", "", region)]
psv_region_long <- psv_region_long[region_no_layer != ""]

hot_region_nl <- psv_region_long[
  ,
  .(n_hotspots = .N),
  by = .(subclass_group, region_no_layer)
]

hot_region_nl[, rank := frank(-n_hotspots, ties.method = "min"),
              by = subclass_group]
hot_region_nl_top <- hot_region_nl[rank <= 10]

top_regions_dt <- unique(hot_region_nl_top[, .(subclass_group, region_no_layer)])
setkey(top_regions_dt, subclass_group, region_no_layer)


base_rdata_dir <- file.path(root_dir, "ws0.4_ss0.02", "cell_window")
rdata_path     <- file.path(base_rdata_dir, paste0(slice_to_plot, ".RData"))

message("Loading cell-level RData: ", rdata_path)
if (!file.exists(rdata_path)) {
  stop("no RData:", rdata_path)
}

load(rdata_path)

if (exists("file_tmp_layer")) {
  cells <- as.data.table(file_tmp_layer)
} else if (exists("file_tmp")) {
  message("no file_tmp_layer")
  cells <- as.data.table(file_tmp)
} else {
  stop("no file_tmp_layer or file_tmp")
}

if (!"slide" %in% names(cells)) {
  if ("brain_section_label" %in% names(cells)) {
    cells[, slide := brain_section_label]
  } else {
    stop("no slide or brain_section_label")
  }
}

if (!"region" %in% names(cells)) {
  if ("ccf_region_name" %in% names(cells)) {
    cells[, region := ccf_region_name]
  } else {
    stop("no region or ccf_region_name")
  }
}

if (!"subclass_group" %in% names(cells)) {
  if (!"subclass" %in% names(cells)) {
    stop("no subclass_group or subclass")
  }
  cells[, subclass_lower := tolower(subclass)]
  cells[, subclass_group := fifelse(
    grepl("pvalb", subclass_lower), "PVALB",
    fifelse(
      grepl("vip",  subclass_lower), "VIP",
      fifelse(
        grepl("sst|somatostatin", subclass_lower), "SOM",
        NA_character_
      )
    )
  )]
}

if (!("x" %in% names(cells) && "y" %in% names(cells))) {
  x_col <- grep("(^x$|_x$|^X$|^pos_x$|^center_x$)",
                names(cells), ignore.case = TRUE, value = TRUE)[1]
  y_col <- grep("(^y$|_y$|^Y$|^pos_y$|^center_y$)",
                names(cells), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(x_col) || is.na(y_col)) {
    stop("no loc")
  }
  setnames(cells, c(x_col, y_col), c("x", "y"))
}

cells_slice <- cells[slide == slice_to_plot]

cells_slice[, region_no_layer :=
              sub(",\\s*layer\\s+[^,]+$", "", region)]

cells_psv <- cells_slice[subclass_group %in% c("PVALB", "SOM", "VIP")]

setkey(cells_psv, subclass_group, region_no_layer)
cells_psv[top_regions_dt, in_top_region := TRUE]
cells_psv[is.na(in_top_region), in_top_region := FALSE]

bg_n     <- 50000L
bg_cells <- cells_slice
if (nrow(bg_cells) > bg_n) {
  set.seed(1)
  bg_cells <- bg_cells[sample(.N, bg_n)]
}

plot_slice_single <- function(target) {
  ggplot() +
    geom_point(data = bg_cells,
               aes(x = x, y = y),
               colour = "grey95", size = 0.05, alpha = 0.5) +
    geom_point(data = cells_psv[subclass_group == target &
                                  in_top_region == FALSE],
               aes(x = x, y = y),
               colour = "grey70", size = 0.15, alpha = 0.3) +
    geom_point(data = cells_psv[subclass_group == target &
                                  in_top_region == TRUE],
               aes(x = x, y = y),
               colour = psv_cols[[target]],
               size = 0.35, alpha = 0.9) +
    coord_fixed() +
    scale_y_reverse() +
    labs(
      x = "x (µm)",
      y = "y (µm)",
      title = target
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text       = element_text(colour = "black"),
      plot.title      = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none",
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

p_slice_PVALB <- plot_slice_single("PVALB")
p_slice_SOM   <- plot_slice_single("SOM")
p_slice_VIP   <- plot_slice_single("VIP")

p_slice_all <- ggplot() +
  geom_point(data = bg_cells,
             aes(x = x, y = y),
             colour = "grey95", size = 0.05, alpha = 0.5) +
  geom_point(data = cells_psv[in_top_region == FALSE],
             aes(x = x, y = y),
             colour = "grey70", size = 0.15, alpha = 0.3) +
  geom_point(data = cells_psv[in_top_region == TRUE],
             aes(x = x, y = y, colour = subclass_group),
             size = 0.35, alpha = 0.9) +
  scale_colour_manual(values = psv_cols, name = "Subclass") +
  coord_fixed() +
  scale_y_reverse() +
  labs(
    x = "x (µm)",
    y = "y (µm)",
    title = paste0("PVALB SOM VIP — ", slice_to_plot)
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text       = element_text(colour = "black"),
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    plot.background = element_rect(fill = "white", colour = NA)
  )

p_slice_2x2 <- cowplot::plot_grid(
  p_slice_PVALB,
  p_slice_SOM,
  p_slice_VIP,
  p_slice_all + theme(legend.position = "none"),
  labels     = c("A", "B", "C", "D"),
  label_size = 14,
  ncol       = 2, nrow       = 2
)

p_slice_2x2

ggsave(
  file.path(fig_dir_main, "04_Fig_example_slice_PSV_topRegions_2x2.png"),
  p_slice_2x2,
  width = 10, height = 8, dpi = 400
)

ggsave(
  file.path(fig_dir_main, "04_Fig_example_slice_PSV_topRegions_allMerged.png"),
  p_slice_all,
  width = 5.5, height = 4.5, dpi = 400
)


all_regions_dt <- unique(hot_region_nl[, .(subclass_group, region_no_layer)])
setkey(all_regions_dt, subclass_group, region_no_layer)

cells_psv_hot <- copy(cells_psv)
setkey(cells_psv_hot, subclass_group, region_no_layer)
cells_psv_hot[all_regions_dt, in_any_hotspot_region := TRUE]
cells_psv_hot[is.na(in_any_hotspot_region), in_any_hotspot_region := FALSE]

cells_psv_hot <- cells_psv_hot[in_any_hotspot_region == TRUE]

plot_slice_single_F5 <- function(target) {
  ggplot() +
    geom_point(data = bg_cells,
               aes(x = x, y = y),
               colour = "grey95", size = 0.05, alpha = 0.5) +
    geom_point(data = cells_psv_hot[subclass_group == target],
               aes(x = x, y = y, alpha = in_top_region),
               colour = psv_cols[[target]],
               size = 0.35) +
    scale_alpha_manual(
      values = c("FALSE" = 0.3, "TRUE" = 1),
      name   = "Region type",
      labels = c("Non-top cluster regions", "Top cluster regions")
    ) +
    coord_fixed() +
    scale_y_reverse() +
    labs(
      x = "x (µm)",
      y = "y (µm)",
      title = target
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text       = element_text(colour = "black"),
      plot.title      = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

p5_slice_PVALB <- plot_slice_single_F5("PVALB")
p5_slice_SOM   <- plot_slice_single_F5("SOM")
p5_slice_VIP   <- plot_slice_single_F5("VIP")

p5_slice_all_simple <- ggplot() +
  geom_point(data = bg_cells,
             aes(x = x, y = y),
             colour = "grey95", size = 0.05, alpha = 0.5) +
  geom_point(data = cells_psv_hot,
             aes(x = x, y = y,
                 colour = subclass_group,
                 alpha  = in_top_region),
             size = 0.35) +
  scale_colour_manual(values = psv_cols, name = "Subclass") +
  scale_alpha_manual(
    values = c("FALSE" = 0.3, "TRUE" = 1),
    name   = "Region type",
    labels = c("Non-top cluster regions", "Top cluster regions"),
    guide  = guide_legend(override.aes = list(size = 3))
  ) +
  coord_fixed() +
  scale_y_reverse() +
  labs(
    x = "x (µm)",
    y = "y (µm)",
    title = paste0("All PVALB/SOM/VIP clusters — ", slice_to_plot)
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text       = element_text(colour = "black"),
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    plot.background = element_rect(fill = "white", colour = NA)
  )

p5_slice_2x2 <- cowplot::plot_grid(
  p5_slice_PVALB,
  p5_slice_SOM,
  p5_slice_VIP,
  p5_slice_all_simple + theme(legend.position = "none"),
  labels     = c("A", "B", "C", "D"),
  label_size = 14,
  ncol       = 2,
  nrow       = 2
)

p5_slice_2x2

ggsave(
  file.path(fig_dir_main, "05_Fig_example_slice_PSV_allClusters_2x2.png"),
  p5_slice_2x2,
  width = 10, height = 8, dpi = 400
)

cells_psv_hot[, region_type := ifelse(in_top_region,
                                      "Top cluster regions",
                                      "Non-top cluster regions")]

cells_psv_hot[, legend_group := factor(
  paste(subclass_group, region_type, sep = " | "),
  levels = c(
    "PVALB | Non-top cluster regions",
    "PVALB | Top cluster regions",
    "SOM | Non-top cluster regions",
    "SOM | Top cluster regions",
    "VIP | Non-top cluster regions",
    "VIP | Top cluster regions"
  )
)]

cols_6 <- c(
  "PVALB | Non-top cluster regions" = psv_cols[["PVALB"]],
  "PVALB | Top cluster regions"     = psv_cols[["PVALB"]],
  "SOM | Non-top cluster regions"   = psv_cols[["SOM"]],
  "SOM | Top cluster regions"       = psv_cols[["SOM"]],
  "VIP | Non-top cluster regions"   = psv_cols[["VIP"]],
  "VIP | Top cluster regions"       = psv_cols[["VIP"]]
)

alpha_6 <- c(
  "PVALB | Non-top cluster regions" = 0.3,
  "PVALB | Top cluster regions"     = 1.0,
  "SOM | Non-top cluster regions"   = 0.3,
  "SOM | Top cluster regions"       = 1.0,
  "VIP | Non-top cluster regions"   = 0.3,
  "VIP | Top cluster regions"       = 1.0
)

p5_slice_all_detailed <- ggplot() +
  geom_point(data = bg_cells,
             aes(x = x, y = y),
             colour = "grey95", size = 0.05, alpha = 0.5) +
  geom_point(data = cells_psv_hot,
             aes(x = x, y = y,
                 colour = legend_group,
                 alpha  = legend_group),
             size = 0.35) +
  scale_colour_manual(
    values = cols_6,
    name   = "Subclass & region type"
  ) +
  scale_alpha_manual(
    values = alpha_6,
    name   = "Subclass & region type",
    guide  = guide_legend(override.aes = list(size = 3))
  ) +
  coord_fixed() +
  scale_y_reverse() +
  labs(
    x = "x (µm)",
    y = "y (µm)",
    title = paste0("All PVALB/SOM/VIP clusters — ", slice_to_plot)
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text       = element_text(colour = "black"),
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    plot.background = element_rect(fill = "white", colour = NA)
  )

p5_slice_all_detailed

ggsave(
  file.path(fig_dir_main, "05_Fig_example_slice_PSV_allClusters_allMerged.png"),
  p5_slice_all_detailed,
  width = 7, height = 5.0, dpi = 400
)
