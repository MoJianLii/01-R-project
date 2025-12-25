rm(list = ls()); gc()

base_dir         <- "E:/zaw/2511/mouseMerfish_zhuang_subclass"

ws_env <- Sys.getenv("WINDOW_SIZE", "0.4")
ss_env <- Sys.getenv("STEP_SIZE",  "0.02")
window_size <- as.numeric(ws_env)
step_size   <- as.numeric(ss_env)

fmt <- function(x){
  s <- formatC(x, format = "f", digits = 6)
  s <- sub("0+$", "", s)
  s <- sub("\\.$", "", s)
  s
}
subdir_tag <- paste0("ws", fmt(window_size), "_ss", fmt(step_size))
combo_root <- file.path(base_dir, subdir_tag)
dir.create(combo_root, recursive = TRUE, showWarnings = FALSE)

setwd(base_dir)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(ggforce)
  library(ggpubr)
  library(FSA)
  library(scales)
  library(grid)
})

fig_dir <- file.path(combo_root, "figures1_subclass")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

mouse_subclass_cluster_total <- fread(
  file.path(combo_root, "mouse_subclass_cluster_total.txt"),
  sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE
)
mouse_subclass_cluster_total_superlow <- mouse_subclass_cluster_total[
  mouse_subclass_cluster_total$enrich_subclass_cell_ids_num < 3, 
]

mouse_subclass_cluster_total_3 <- mouse_subclass_cluster_total[
  mouse_subclass_cluster_total$enrich_subclass_cell_ids_num >= 3, 
]

fwrite(
  mouse_subclass_cluster_total_3,
  file = file.path(combo_root, "mouse_subclass_cluster_total_over3.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

mapping_path <- file.path(combo_root, "Merfish_brain_cell_type_subclass.txt")
Merfish_sub_map <- fread(
  mapping_path,
  sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE
)
colnames(Merfish_sub_map) <- c("subclass", "cell_Neruon_type", "subclass_sim")

Table_2_total_cell <- fread(
  file.path(combo_root, "mouse_table", "table2", "Table_2_total_cell.txt"),
  sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE
)

library(dplyr)
library(ggplot2)
library(scales)

type_cols <- c(
  "Gaba"      = "#6FB06F",
  "Glut"      = "#8C3D4A",
  "NonNeuron" = "#F2A94F"
)
type_labels <- c(
  "Gaba"      = "GABA",
  "Glut"      = "Glut",
  "NonNeuron" = "NonNeuron"
)

cluster_type_counts <- mouse_subclass_cluster_total %>%
  filter(!is.na(cell_Neruon_type)) %>%
  count(cell_Neruon_type, name = "cluster_n") %>%
  mutate(
    pct = cluster_n / sum(cluster_n),
    label_inner = sprintf(
      "%s clusters\n(%s)\n%s",
      type_labels[cell_Neruon_type],
      format(cluster_n, big.mark = ","),
      percent(pct, accuracy = 1)
    )
  ) %>%
  arrange(match(cell_Neruon_type, c("Gaba","Glut","NonNeuron")))

total_clusters <- sum(cluster_type_counts$cluster_n)

p1A_left <- ggplot(cluster_type_counts,
                   aes(x = 2, y = cluster_n, fill = cell_Neruon_type)) +
  geom_col(width = 1, color = "black") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  geom_text(
    aes(label = label_inner),
    position   = position_stack(vjust = 0.5),
    size       = 4,
    fontface   = "bold",
    color      = "black",
    lineheight = 1.1
  ) +
  annotate(
    "text",
    x = 0.5, y = 0,
    label    = paste0("Total clusters\n(",
                      format(total_clusters, big.mark = ","), ")"),
    size     = 4.2,
    fontface = "bold",
    lineheight = 1.1
  ) +
  scale_fill_manual(values = type_cols,
                    name   = "Type",
                    labels = type_labels) +
  theme_void(base_size = 14) +
  theme(
    legend.position = "none"
  )

ggsave(file.path(fig_dir, "Fig1A_left_clusterType_pie.png"),
       p1A_left, width = 5, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "Fig1A_left_clusterType_pie.pdf"),
       p1A_left, width = 5, height = 5)

subclass_counts <- mouse_subclass_cluster_total %>%
  filter(!is.na(subclass_sim),
         !is.na(cell_Neruon_type)) %>%
  count(subclass_sim, cell_Neruon_type, name = "cluster_n") %>%
  arrange(desc(cluster_n))

topN <- 15L
subclass_top <- subclass_counts %>%
  slice_head(n = topN) %>%
  mutate(
    subclass_sim = factor(subclass_sim, levels = subclass_sim)
  )

p1A_right <- ggplot(subclass_top,
                    aes(x = subclass_sim,
                        y = cluster_n,
                        fill = cell_Neruon_type)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 1000,
             linetype   = "dashed",
             color      = "#D55E00",
             linewidth  = 0.6) +
  geom_text(aes(label = cluster_n),
            vjust = -0.25,
            size  = 3.5,
            fontface = "bold") +
  scale_fill_manual(values = type_cols,
                    name   = "Type",
                    labels = type_labels) +
  labs(
    title = "Mouse subclass cluster",
    x     = NULL,
    y     = "Cluster number"
  ) +
  theme_classic(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.text.x  = element_text(angle = 35, hjust = 1, vjust = 1, size = 11),
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 15),
    legend.position = "top",
    legend.title    = element_text(face = "bold")
  )

ggsave(file.path(fig_dir, "Fig1A_right_topSubclass_bar.png"),
       p1A_right, width = 7.5, height = 4.5, dpi = 300)
ggsave(file.path(fig_dir, "Fig1A_right_topSubclass_bar.pdf"),
       p1A_right, width = 7.5, height = 4.5)


suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(ggforce)
  library(data.table)
})

top_sub_sim <- as.character(subclass_top$subclass_sim)
slice_dir   <- file.path(base_dir, "neocortex_new")
slice_files <- list.files(slice_dir, pattern = "\\.txt$", full.names = TRUE)

read_slice_cells <- function(fpath){
  dt <- data.table::fread(
    fpath,
    sep = "\t", header = TRUE,
    data.table = FALSE, check.names = FALSE
  )
  if ("cell_label" %in% names(dt)) {
    dt$cell_id <- as.character(dt$cell_label)
  } else if ("cell_id" %in% names(dt)) {
    dt$cell_id <- as.character(dt$cell_id)
  } else {
    return(NULL)
  }
  if (!"subclass" %in% names(dt)) return(NULL)
  
  dt$subclass <- as.character(dt$subclass)
  dt <- dt[, c("cell_id", "subclass"), drop = FALSE]
  
  dt <- dt %>%
    dplyr::left_join(
      Merfish_sub_map[, c("subclass", "subclass_sim")],
      by = "subclass"
    )
  dt
}

cell_list <- lapply(slice_files, read_slice_cells)
cells_all <- bind_rows(cell_list) %>%
  dplyr::filter(!is.na(subclass_sim)) %>%
  dplyr::distinct(cell_id, .keep_all = TRUE)

cells_top <- cells_all %>%
  dplyr::filter(subclass_sim %in% top_sub_sim)

total_ids_by_sub <- split(cells_top$cell_id, cells_top$subclass_sim)

if (!"label" %in% names(mouse_subclass_cluster_total)) {
  mouse_subclass_cluster_total$label <- with(
    mouse_subclass_cluster_total,
    paste(slide, layer, subclass, merge_regions, sep = "_")
  )
}
if (!"subclass_sim" %in% names(mouse_subclass_cluster_total)) {
  mouse_subclass_cluster_total <- mouse_subclass_cluster_total %>%
    dplyr::left_join(
      Merfish_sub_map[, c("subclass", "subclass_sim")],
      by = "subclass"
    )
}

split_ids <- function(x){
  if (is.null(x) || length(x) == 0) return(character(0))
  ids <- unlist(str_split(x, "[,]"))
  ids <- ids[ids != "" & !is.na(ids)]
  ids
}

donut_stats_list <- vector("list", length(top_sub_sim))
names(donut_stats_list) <- top_sub_sim

for (s in top_sub_sim) {
  
  total_ids <- unique(total_ids_by_sub[[s]])
  total_n   <- length(total_ids)
  
  if (is.null(total_ids) || total_n == 0) {
    donut_stats_list[[s]] <- data.frame(
      subclass_sim  = s,
      total_cells   = 0,
      cluster_cells = 0,
      stringsAsFactors = FALSE
    )
    next
  }
  
  clu_s <- mouse_subclass_cluster_total %>%
    dplyr::filter(subclass_sim == s)
  
  if (!nrow(clu_s) || !"enrich_subclass_cell_ids" %in% names(clu_s)) {
    donut_stats_list[[s]] <- data.frame(
      subclass_sim  = s,
      total_cells   = total_n,
      cluster_cells = 0,
      stringsAsFactors = FALSE
    )
    next
  }
  
  cluster_ids_list <- lapply(clu_s$enrich_subclass_cell_ids, split_ids)
  names(cluster_ids_list) <- clu_s$label
  
  cluster_cell_ids <- unique(unlist(cluster_ids_list))
  cluster_cell_ids <- intersect(cluster_cell_ids, total_ids)
  cluster_n <- length(cluster_cell_ids)
  
  donut_stats_list[[s]] <- data.frame(
    subclass_sim  = s,
    total_cells   = total_n,
    cluster_cells = cluster_n,
    stringsAsFactors = FALSE
  )
}

donut_stats <- bind_rows(donut_stats_list) %>%
  mutate(
    non_cluster_cells = pmax(total_cells - cluster_cells, 0),
    pct_cluster       = ifelse(total_cells == 0, NA, cluster_cells     / total_cells),
    pct_noncluster    = ifelse(total_cells == 0, NA, non_cluster_cells / total_cells)
  )

ring_df <- donut_stats %>%
  transmute(
    class = subclass_sim,
    total_cells,
    Clustered       = cluster_cells,
    `Non-clustered` = non_cluster_cells,
    pct_cluster,
    pct_noncluster
  ) %>%
  pivot_longer(
    cols      = c(Clustered, `Non-clustered`),
    names_to  = "group",
    values_to = "count"
  ) %>%
  mutate(
    group = factor(group, levels = c("Clustered", "Non-clustered")),
    pct   = ifelse(group == "Clustered", pct_cluster, pct_noncluster)
  ) %>%
  group_by(class) %>%
  arrange(group, .by_group = TRUE) %>%
  mutate(
    total_count = sum(count),
    frac        = ifelse(total_count > 0, count / total_count, 0),
    ymax        = cumsum(frac),
    ymin        = lag(ymax, default = 0),
    start       = 2 * pi * ymin,
    end         = 2 * pi * ymax,
    mid         = (start + end) / 2,
    r0          = 0.6,
    r           = 1.2,
    r_label_out = r + 0.25,
    lx          = r_label_out * sin(mid),
    ly          = r_label_out * cos(mid),
    hjust_lab   = ifelse(lx >= 0, 0, 1),
    label       = paste0(
      group, "\n",
      format(count, big.mark = ","), " cells\n",
      ifelse(is.na(pct),
             "NA",
             ifelse(pct < 0.01, "<1%", paste0(round(pct * 100), "%")))
    )
  ) %>%
  ungroup() %>%
  dplyr::filter(total_count > 0, count > 0)

center_df <- donut_stats %>%
  transmute(
    class        = subclass_sim,
    center_label = paste0(format(total_cells, big.mark = ","), " cells")
  )

fill_colors <- c(
  "Clustered"      = "#56B4E9",
  "Non-clustered"  = "#E69F00"
)

p1B <- ggplot() +
  geom_arc_bar(
    data = ring_df,
    aes(
      x0 = 0, y0 = 0,
      r0 = r0, r = r,
      start = start, end = end,
      fill  = group
    ),
    color     = "white",
    linewidth = 0.5
  ) +
  geom_text(
    data = ring_df,
    aes(x = lx, y = ly, label = label, hjust = hjust_lab),
    size       = 3,
    fontface   = "bold",
    lineheight = 1.05
  ) +
  geom_text(
    data = center_df,
    aes(x = 0, y = 0, label = center_label),
    size       = 3.0,
    fontface   = "bold",
    lineheight = 1.05
  ) +
facet_wrap(~ class, ncol = 3) +
  scale_fill_manual(
    values = fill_colors,
    breaks = names(fill_colors),
    name   = NULL
  ) +
  coord_fixed(xlim = c(-1.7, 1.7), ylim = c(-1.7, 1.7), clip = "off") +
  theme_void(base_size = 12) +
  theme(
    strip.text      = element_text(size = 11, face = "bold"),
    legend.position = "top",
    legend.text     = element_text(size = 10),
    panel.spacing.x = grid::unit(6, "lines"),
    panel.spacing.y = grid::unit(1.5, "lines"),
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.margin     = margin(5, 10, 5, 10)
  ) +
  labs(title = "Mouse class cell composition")

ggsave(file.path(fig_dir, "Fig1B_donut_topSubclasses_cluster_vs_noncluster.png"),
       p1B, width = 16, height = 9, dpi = 400)
ggsave(file.path(fig_dir, "Fig1B_donut_topSubclasses_cluster_vs_noncluster.pdf"),
       p1B, width = 16, height = 9)



library(dplyr)
library(ggplot2)
library(scales)

type_cols <- c(
  "Gaba"      = "#6FB06F",
  "Glut"      = "#8C3D4A",
  "NonNeuron" = "#F2A94F"
)
type_labels <- c(
  "Gaba"      = "GABA",
  "Glut"      = "Glut",
  "NonNeuron" = "NonNeuron"
)

cluster_type_counts <- mouse_subclass_cluster_total_3 %>%
  filter(!is.na(cell_Neruon_type)) %>%
  count(cell_Neruon_type, name = "cluster_n") %>%
  mutate(
    pct = cluster_n / sum(cluster_n),
    label_inner = sprintf(
      "%s clusters\n(%s)\n%s",
      type_labels[cell_Neruon_type],
      format(cluster_n, big.mark = ","),
      percent(pct, accuracy = 1)
    )
  ) %>%
  arrange(match(cell_Neruon_type, c("Gaba","Glut","NonNeuron")))

total_clusters <- sum(cluster_type_counts$cluster_n)

p1A_left <- ggplot(cluster_type_counts,
                   aes(x = 2, y = cluster_n, fill = cell_Neruon_type)) +
  geom_col(width = 1, color = "black") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  geom_text(
    aes(label = label_inner),
    position   = position_stack(vjust = 0.5),
    size       = 4,
    fontface   = "bold",
    color      = "black",
    lineheight = 1.1
  ) +
  annotate(
    "text",
    x = 0.5, y = 0,
    label    = paste0("Total clusters\n(",
                      format(total_clusters, big.mark = ","), ")"),
    size     = 4.2,
    fontface = "bold",
    lineheight = 1.1
  ) +
  scale_fill_manual(values = type_cols,
                    name   = "Type",
                    labels = type_labels) +
  theme_void(base_size = 14) +
  theme(
    legend.position = "none"
  )

ggsave(file.path(fig_dir, "Fig1A_left_clusterType_pie_3.png"),
       p1A_left, width = 5, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "Fig1A_left_clusterType_pie_3.pdf"),
       p1A_left, width = 5, height = 5)

subclass_counts <- mouse_subclass_cluster_total_3 %>%
  filter(!is.na(subclass_sim),
         !is.na(cell_Neruon_type)) %>%
  count(subclass_sim, cell_Neruon_type, name = "cluster_n") %>%
  arrange(desc(cluster_n))

topN <- 15L
subclass_top <- subclass_counts %>%
  slice_head(n = topN) %>%
  mutate(
    subclass_sim = factor(subclass_sim, levels = subclass_sim)
  )

p1A_right <- ggplot(subclass_top,
                    aes(x = subclass_sim,
                        y = cluster_n,
                        fill = cell_Neruon_type)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 1000,
             linetype   = "dashed",
             color      = "#D55E00",
             linewidth  = 0.6) +
  geom_text(aes(label = cluster_n),
            vjust = -0.25,
            size  = 3.5,
            fontface = "bold") +
  scale_fill_manual(values = type_cols,
                    name   = "Type",
                    labels = type_labels) +
  labs(
    title = "Mouse subclass cluster",
    x     = NULL,
    y     = "Cluster number"
  ) +
  theme_classic(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.text.x  = element_text(angle = 35, hjust = 1, vjust = 1, size = 11),
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 15),
    legend.position = "top",
    legend.title    = element_text(face = "bold")
  )

ggsave(file.path(fig_dir, "Fig1A_right_topSubclass_bar_3.png"),
       p1A_right, width = 7.5, height = 4.5, dpi = 300)
ggsave(file.path(fig_dir, "Fig1A_right_topSubclass_bar_3.pdf"),
       p1A_right, width = 7.5, height = 4.5)


suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(ggforce)
  library(data.table)
})

top_sub_sim <- as.character(subclass_top$subclass_sim)
slice_dir   <- file.path(base_dir, "neocortex_new")
slice_files <- list.files(slice_dir, pattern = "\\.txt$", full.names = TRUE)

read_slice_cells <- function(fpath){
  dt <- data.table::fread(
    fpath,
    sep = "\t", header = TRUE,
    data.table = FALSE, check.names = FALSE
  )
  if ("cell_label" %in% names(dt)) {
    dt$cell_id <- as.character(dt$cell_label)
  } else if ("cell_id" %in% names(dt)) {
    dt$cell_id <- as.character(dt$cell_id)
  } else {
    return(NULL)
  }
  if (!"subclass" %in% names(dt)) return(NULL)
  
  dt$subclass <- as.character(dt$subclass)
  dt <- dt[, c("cell_id", "subclass"), drop = FALSE]
  
  dt <- dt %>%
    dplyr::left_join(
      Merfish_sub_map[, c("subclass", "subclass_sim")],
      by = "subclass"
    )
  dt
}

cell_list <- lapply(slice_files, read_slice_cells)
cells_all <- bind_rows(cell_list) %>%
  dplyr::filter(!is.na(subclass_sim)) %>%
  dplyr::distinct(cell_id, .keep_all = TRUE)

cells_top <- cells_all %>%
  dplyr::filter(subclass_sim %in% top_sub_sim)

total_ids_by_sub <- split(cells_top$cell_id, cells_top$subclass_sim)

if (!"label" %in% names(mouse_subclass_cluster_total_3)) {
  mouse_subclass_cluster_total_3$label <- with(
    mouse_subclass_cluster_total_3,
    paste(slide, layer, subclass, merge_regions, sep = "_")
  )
}
if (!"subclass_sim" %in% names(mouse_subclass_cluster_total_3)) {
  mouse_subclass_cluster_total_3 <- mouse_subclass_cluster_total_3 %>%
    dplyr::left_join(
      Merfish_sub_map[, c("subclass", "subclass_sim")],
      by = "subclass"
    )
}

split_ids <- function(x){
  if (is.null(x) || length(x) == 0) return(character(0))
  ids <- unlist(str_split(x, "[,]"))
  ids <- ids[ids != "" & !is.na(ids)]
  ids
}

donut_stats_list <- vector("list", length(top_sub_sim))
names(donut_stats_list) <- top_sub_sim

for (s in top_sub_sim) {
  
  total_ids <- unique(total_ids_by_sub[[s]])
  total_n   <- length(total_ids)
  
  if (is.null(total_ids) || total_n == 0) {
    donut_stats_list[[s]] <- data.frame(
      subclass_sim  = s,
      total_cells   = 0,
      cluster_cells = 0,
      stringsAsFactors = FALSE
    )
    next
  }
  
  clu_s <- mouse_subclass_cluster_total_3 %>%
    dplyr::filter(subclass_sim == s)
  
  if (!nrow(clu_s) || !"enrich_subclass_cell_ids" %in% names(clu_s)) {
    donut_stats_list[[s]] <- data.frame(
      subclass_sim  = s,
      total_cells   = total_n,
      cluster_cells = 0,
      stringsAsFactors = FALSE
    )
    next
  }
  
  cluster_ids_list <- lapply(clu_s$enrich_subclass_cell_ids, split_ids)
  names(cluster_ids_list) <- clu_s$label
  
  cluster_cell_ids <- unique(unlist(cluster_ids_list))
  cluster_cell_ids <- intersect(cluster_cell_ids, total_ids)
  cluster_n <- length(cluster_cell_ids)
  
  donut_stats_list[[s]] <- data.frame(
    subclass_sim  = s,
    total_cells   = total_n,
    cluster_cells = cluster_n,
    stringsAsFactors = FALSE
  )
}

donut_stats <- bind_rows(donut_stats_list) %>%
  mutate(
    non_cluster_cells = pmax(total_cells - cluster_cells, 0),
    pct_cluster       = ifelse(total_cells == 0, NA, cluster_cells     / total_cells),
    pct_noncluster    = ifelse(total_cells == 0, NA, non_cluster_cells / total_cells)
  )

ring_df <- donut_stats %>%
  transmute(
    class = subclass_sim,
    total_cells,
    Clustered       = cluster_cells,
    `Non-clustered` = non_cluster_cells,
    pct_cluster,
    pct_noncluster
  ) %>%
  pivot_longer(
    cols      = c(Clustered, `Non-clustered`),
    names_to  = "group",
    values_to = "count"
  ) %>%
  mutate(
    group = factor(group, levels = c("Clustered", "Non-clustered")),
    pct   = ifelse(group == "Clustered", pct_cluster, pct_noncluster)
  ) %>%
  group_by(class) %>%
  arrange(group, .by_group = TRUE) %>%
  mutate(
    total_count = sum(count),
    frac        = ifelse(total_count > 0, count / total_count, 0),
    ymax        = cumsum(frac),
    ymin        = lag(ymax, default = 0),
    start       = 2 * pi * ymin,
    end         = 2 * pi * ymax,
    mid         = (start + end) / 2,
    r0          = 0.6,
    r           = 1.2,
    r_label_out = r + 0.25,
    lx          = r_label_out * sin(mid),
    ly          = r_label_out * cos(mid),
    hjust_lab   = ifelse(lx >= 0, 0, 1),
    label       = paste0(
      group, "\n",
      format(count, big.mark = ","), " cells\n",
      ifelse(is.na(pct),
             "NA",
             ifelse(pct < 0.01, "<1%", paste0(round(pct * 100), "%")))
    )
  ) %>%
  ungroup() %>%
  dplyr::filter(total_count > 0, count > 0)

center_df <- donut_stats %>%
  transmute(
    class        = subclass_sim,
    center_label = paste0(format(total_cells, big.mark = ","), " cells")
  )

fill_colors <- c(
  "Clustered"      = "#56B4E9",
  "Non-clustered"  = "#E69F00"
)

p1B <- ggplot() +
  geom_arc_bar(
    data = ring_df,
    aes(
      x0 = 0, y0 = 0,
      r0 = r0, r = r,
      start = start, end = end,
      fill  = group
    ),
    color     = "white",
    linewidth = 0.5
  ) +
  geom_text(
    data = ring_df,
    aes(x = lx, y = ly, label = label, hjust = hjust_lab),
    size       = 3,
    fontface   = "bold",
    lineheight = 1.05
  ) +
  geom_text(
    data = center_df,
    aes(x = 0, y = 0, label = center_label),
    size       = 3.0,
    fontface   = "bold",
    lineheight = 1.05
  ) +
facet_wrap(~ class, ncol = 3) +
  scale_fill_manual(
    values = fill_colors,
    breaks = names(fill_colors),
    name   = NULL
  ) +
  coord_fixed(xlim = c(-1.7, 1.7), ylim = c(-1.7, 1.7), clip = "off") +
  theme_void(base_size = 12) +
  theme(
    strip.text      = element_text(size = 11, face = "bold"),
    legend.position = "top",
    legend.text     = element_text(size = 10),
    panel.spacing.x = grid::unit(6, "lines"),
    panel.spacing.y = grid::unit(1.5, "lines"),
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.margin     = margin(5, 10, 5, 10)
  ) +
  labs(title = "Mouse class cell composition")

ggsave(file.path(fig_dir, "Fig1B_donut_topSubclasses_cluster_vs_noncluster_3.png"),
       p1B, width = 16, height = 9, dpi = 400)
ggsave(file.path(fig_dir, "Fig1B_donut_topSubclasses_cluster_vs_noncluster_3.pdf"),
       p1B, width = 16, height = 9)

