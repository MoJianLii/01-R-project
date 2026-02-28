rm(list = ls()); gc()

## ===================== 0) 基础路径与参数 =====================
base_dir <- "E:/zaw/2511/mouseMerfish_zhuang_subclass"
setwd(base_dir)

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

fig_dir <- file.path(combo_root, "figures1_subclass")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## ===================== 1) 加载依赖 =====================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(ggtext)
  library(ggpubr)
  library(ggExtra)
  library(ggforce)
  library(viridis)
  library(grid)
  library(forcats)
})

## ===================== 2) 读入核心数据 =====================
mouse_subclass_cluster_total <- fread(
  file.path(combo_root, "mouse_subclass_cluster_total.txt"),
  sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE
)

mouse_subclass_cluster_total_3 <- mouse_subclass_cluster_total[
  mouse_subclass_cluster_total$enrich_subclass_cell_ids_num >= 3, 
]

fwrite(
  mouse_subclass_cluster_total_3,
  file = file.path(combo_root, "mouse_subclass_cluster_total_over3.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

mapping_path <- file.path(combo_root, "Merfish_brain_cell_type_subclass.txt")
Merfish_sub_map <- fread(mapping_path, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
colnames(Merfish_sub_map) <- c("subclass", "cell_Neruon_type", "subclass_sim")

tab <- fread(
  file.path(combo_root, "mouse_table", "table2", "Table_2_total_cell.txt"),
  sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE
)

## --------------------- 2.1 补齐 cluster_total 的 label/subclass_sim ---------------------
clu_tbl <- mouse_subclass_cluster_total_3

if (!"label" %in% names(clu_tbl)) {
  stopifnot(all(c("slide","layer","subclass","merge_regions") %in% names(clu_tbl)))
  clu_tbl$label <- with(clu_tbl, paste(slide, layer, subclass, merge_regions, sep = "_"))
}
if (!"subclass_sim" %in% names(clu_tbl)) {
  clu_tbl <- clu_tbl %>%
    left_join(Merfish_sub_map[, c("subclass","subclass_sim")], by = "subclass")
}

## ================================================================
## 通用主题
## ================================================================
theme_paper <- function(base=12){
  theme_classic(base_size = base) %+replace%
    theme(
      panel.border = element_rect(color="black", fill=NA, linewidth=0.6),
      axis.text    = element_text(color="black"),
      plot.title   = element_text(hjust=0.5, face="bold", margin=margin(b=6)),
      plot.subtitle= element_text(hjust=0.5, margin=margin(b=6)),
      legend.position = "top",
      legend.title = element_text(face="bold"),
      plot.margin = margin(6, 10, 6, 6)
    )
}

theme_topjournal <- function(base = 11){
  theme_minimal(base_size = base) %+replace%
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, margin = margin(0,0,6,0)),
      axis.title.x = element_text(margin = margin(t = 6)),
      axis.title.y = element_text(margin = margin(r = 6)),
      axis.text    = element_text(color = "grey20"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "top",
      legend.title = element_blank(),
      plot.margin = margin(6, 10, 6, 6)
    )
}

## ================================================================
## ============================== Fig1 ==============================
## ================================================================

## ---------- Fig1A：左饼图（3种 type cluster 占比） ----------
type_cols <- c("Gaba"="#6FB06F", "Glut"="#8C3D4A", "NonNeuron"="#F2A94F")
type_labels <- c("Gaba"="GABA", "Glut"="Glut", "NonNeuron"="NonNeuron")

cluster_type_counts <- clu_tbl %>%
  filter(!is.na(cell_Neruon_type)) %>%
  count(cell_Neruon_type, name = "cluster_n") %>%
  mutate(
    pct = cluster_n / sum(cluster_n),
    label_inner = sprintf("%s clusters\n(%s)\n%s",
                          type_labels[cell_Neruon_type],
                          format(cluster_n, big.mark = ","),
                          percent(pct, accuracy = 1))
  ) %>%
  arrange(match(cell_Neruon_type, c("Gaba","Glut","NonNeuron")))

total_clusters <- sum(cluster_type_counts$cluster_n)

p1A_left <- ggplot(cluster_type_counts, aes(x = 2, y = cluster_n, fill = cell_Neruon_type)) +
  geom_col(width = 1, color = "black") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  geom_text(aes(label = label_inner),
            position = position_stack(vjust = 0.5),
            size = 4, fontface = "bold", color = "black", lineheight = 1.1) +
  annotate("text", x = 0.5, y = 0,
           label = paste0("Total clusters\n(", format(total_clusters, big.mark=","), ")"),
           size = 4.2, fontface = "bold", lineheight = 1.1) +
  scale_fill_manual(values = type_cols, labels = type_labels) +
  theme_void(base_size = 14) +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "Fig1A_left_clusterType_pie_over3.png"), p1A_left, width = 5, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "Fig1A_left_clusterType_pie_over3.pdf"), p1A_left, width = 5, height = 5)

## ---------- Fig1A：右柱图（TopN subclass_sim，✅先按subclass汇总选TopN，再回填分type堆叠） ----------
topN <- 15L

subclass_counts <- clu_tbl %>%
  filter(!is.na(subclass_sim), !is.na(cell_Neruon_type)) %>%
  count(subclass_sim, cell_Neruon_type, name = "cluster_n")

subclass_rank <- subclass_counts %>%
  group_by(subclass_sim) %>%
  summarise(total_n = sum(cluster_n), .groups = "drop") %>%
  slice_max(total_n, n = topN, with_ties = FALSE)

subclass_top <- subclass_counts %>%
  filter(subclass_sim %in% subclass_rank$subclass_sim) %>%
  left_join(subclass_rank, by = "subclass_sim") %>%
  mutate(subclass_sim = factor(subclass_sim, levels = subclass_rank$subclass_sim[order(subclass_rank$total_n, decreasing = TRUE)]))

p1A_right <- ggplot(subclass_top,
                    aes(x = subclass_sim, y = cluster_n, fill = cell_Neruon_type)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 1000, linetype = "dashed", color = "#D55E00", linewidth = 0.6) +
  geom_text(aes(label = cluster_n),
            position = position_stack(vjust = 1.02),
            size = 3.2, fontface = "bold") +
  scale_fill_manual(values = type_cols, labels = type_labels) +
  labs(title = "Mouse subclass cluster (Top subclasses)", x = NULL, y = "Cluster number") +
  theme_paper(13) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 11),
        legend.position = "top")

ggsave(file.path(fig_dir, "Fig1A_right_topSubclass_bar_over3.png"), p1A_right, width = 7.8, height = 4.8, dpi = 300)
ggsave(file.path(fig_dir, "Fig1A_right_topSubclass_bar_over3.pdf"), p1A_right, width = 7.8, height = 4.8)

## ---------- Fig1B：TopN subclass 的 Clustered vs Non-clustered donut（✅cell_uid 全局唯一） ----------
top_sub_sim <- as.character(subclass_rank$subclass_sim)

slice_dir   <- file.path(base_dir, "neocortex_new")
slice_files <- list.files(slice_dir, pattern = "\\.txt$", full.names = TRUE)

read_slice_cells <- function(fpath){
  slice_id <- tools::file_path_sans_ext(basename(fpath))
  dt <- data.table::fread(fpath, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
  
  if ("cell_label" %in% names(dt)) {
    dt$cell_id <- as.character(dt$cell_label)
  } else if ("cell_id" %in% names(dt)) {
    dt$cell_id <- as.character(dt$cell_id)
  } else return(NULL)
  
  if (!"subclass" %in% names(dt)) return(NULL)
  dt$subclass <- as.character(dt$subclass)
  
  dt <- dt[, c("cell_id","subclass"), drop = FALSE]
  dt$cell_uid <- paste0(slice_id, "|", dt$cell_id)
  
  dt <- dt %>%
    left_join(Merfish_sub_map[, c("subclass","subclass_sim")], by = "subclass")
  
  dt[, c("cell_uid","subclass_sim"), drop = FALSE]
}

cell_list <- lapply(slice_files, read_slice_cells)
cells_all <- bind_rows(cell_list) %>%
  filter(!is.na(subclass_sim)) %>%
  distinct(cell_uid, .keep_all = TRUE)

cells_top <- cells_all %>% filter(subclass_sim %in% top_sub_sim)
total_ids_by_sub <- split(cells_top$cell_uid, cells_top$subclass_sim)

split_ids <- function(x){
  if (is.null(x) || length(x) == 0) return(character(0))
  ids <- unlist(str_split(as.character(x), "[,]"))
  ids <- ids[ids != "" & !is.na(ids)]
  ids
}

donut_stats_list <- vector("list", length(top_sub_sim))
names(donut_stats_list) <- top_sub_sim

for (s in top_sub_sim) {
  total_ids <- unique(total_ids_by_sub[[s]])
  total_n   <- length(total_ids)
  
  if (is.null(total_ids) || total_n == 0) {
    donut_stats_list[[s]] <- data.frame(subclass_sim=s, total_cells=0, cluster_cells=0)
    next
  }
  
  clu_s <- clu_tbl %>% filter(subclass_sim == s)
  if (!nrow(clu_s) || !"enrich_subclass_cell_ids" %in% names(clu_s)) {
    donut_stats_list[[s]] <- data.frame(subclass_sim=s, total_cells=total_n, cluster_cells=0)
    next
  }
  
  cluster_ids_list <- lapply(seq_len(nrow(clu_s)), function(i){
    ids0 <- split_ids(clu_s$enrich_subclass_cell_ids[i])
    slide_i <- as.character(clu_s$slide[i])
    ifelse(grepl("\\|", ids0), ids0, paste0(slide_i, "|", ids0))
  })
  
  cluster_cell_ids <- unique(unlist(cluster_ids_list))
  cluster_cell_ids <- intersect(cluster_cell_ids, total_ids)
  cluster_n <- length(cluster_cell_ids)
  
  donut_stats_list[[s]] <- data.frame(subclass_sim=s, total_cells=total_n, cluster_cells=cluster_n)
}

donut_stats <- bind_rows(donut_stats_list) %>%
  mutate(
    non_cluster_cells = pmax(total_cells - cluster_cells, 0),
    pct_cluster       = ifelse(total_cells == 0, NA, cluster_cells / total_cells),
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
  pivot_longer(cols = c(Clustered, `Non-clustered`), names_to = "group", values_to = "count") %>%
  mutate(group = factor(group, levels = c("Clustered","Non-clustered")),
         pct   = ifelse(group=="Clustered", pct_cluster, pct_noncluster)) %>%
  group_by(class) %>%
  arrange(group, .by_group = TRUE) %>%
  mutate(
    total_count = sum(count),
    frac  = ifelse(total_count>0, count/total_count, 0),
    ymax  = cumsum(frac),
    ymin  = lag(ymax, default = 0),
    start = 2*pi*ymin,
    end   = 2*pi*ymax,
    mid   = (start+end)/2,
    r0    = 0.6,
    r     = 1.2,
    r_lab = r + 0.25,
    lx    = r_lab * sin(mid),
    ly    = r_lab * cos(mid),
    hjust_lab = ifelse(lx >= 0, 0, 1),
    label = paste0(
      group, "\n",
      format(count, big.mark=","), " cells\n",
      ifelse(is.na(pct), "NA",
             ifelse(pct < 0.01, "<1%", paste0(round(pct*100), "%")))
    )
  ) %>%
  ungroup() %>%
  filter(total_count > 0, count > 0)

center_df <- donut_stats %>%
  transmute(class = subclass_sim,
            center_label = paste0(format(total_cells, big.mark=","), " cells"))

fill_colors <- c("Clustered"="#56B4E9", "Non-clustered"="#E69F00")

p1B <- ggplot() +
  geom_arc_bar(
    data = ring_df,
    aes(x0=0, y0=0, r0=r0, r=r, start=start, end=end, fill=group),
    color="white", linewidth=0.5
  ) +
  geom_text(data = ring_df,
            aes(x=lx, y=ly, label=label, hjust=hjust_lab),
            size=3, fontface="bold", lineheight=1.05) +
  geom_text(data = center_df, aes(x=0, y=0, label=center_label),
            size=3.0, fontface="bold", lineheight=1.05) +
  facet_wrap(~ class, ncol = 3) +
  scale_fill_manual(values = fill_colors, name = NULL) +
  coord_fixed(xlim=c(-1.7,1.7), ylim=c(-1.7,1.7), clip="off") +
  theme_void(base_size = 12) +
  theme(
    strip.text      = element_text(size=11, face="bold"),
    legend.position = "top",
    legend.text     = element_text(size=10),
    panel.spacing.x = unit(6, "lines"),
    panel.spacing.y = unit(1.5, "lines"),
    plot.title      = element_text(hjust=0.5, face="bold", size=16),
    plot.margin     = margin(5,10,5,10)
  ) +
  labs(title = "Mouse class cell composition (Top subclasses)")

ggsave(file.path(fig_dir, "Fig1B_donut_topSubclasses_cluster_vs_noncluster_over3.png"), p1B, width=16, height=9, dpi=400)
ggsave(file.path(fig_dir, "Fig1B_donut_topSubclasses_cluster_vs_noncluster_over3.pdf"), p1B, width=16, height=9)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(cowplot)
  library(Cairo)
})

layer_levels <- c("1","2/3","4","5","6a","6b")
neuron_types <- c("Gaba","Glut")

clu_tbl <- if (exists("mouse_subclass_cluster_total_3")) mouse_subclass_cluster_total_3 else mouse_subclass_cluster_total
stopifnot(exists("Merfish_sub_map"))
stopifnot(exists("base_dir"))
stopifnot(exists("combo_root"))
if (!exists("fig_dir")) {
  fig_dir <- file.path(combo_root, "figures1_subclass")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
}

Merfish_sub_map_dt <- as.data.table(Merfish_sub_map)
Merfish_sub_map_dt <- Merfish_sub_map_dt[, .(subclass, cell_Neruon_type)]

dt_clu <- as.data.table(clu_tbl)
if (!"cell_Neruon_type" %in% names(dt_clu)) {
  dt_clu <- merge(dt_clu, Merfish_sub_map_dt, by = "subclass", all.x = TRUE)
}
dt_clu <- dt_clu[cell_Neruon_type %chin% neuron_types]
dt_clu <- dt_clu[!is.na(enrich_subclass_cell_ids) & enrich_subclass_cell_ids != ""]
stopifnot("slide" %in% names(dt_clu))

split_ids <- function(x){
  if (is.null(x) || length(x) == 0 || is.na(x) || x == "") return(character(0))
  ids <- unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE)
  ids <- ids[ids != "" & !is.na(ids)]
  ids
}

make_uid <- function(ids, slide){
  if (!length(ids)) return(character(0))
  if (any(grepl("\\|", ids))) {
    ids
  } else {
    paste0(slide, "|", ids)
  }
}

cluster_cell_uid <- unique(unlist(
  Map(function(x, s) make_uid(split_ids(x), s),
      dt_clu$enrich_subclass_cell_ids, dt_clu$slide),
  use.names = FALSE
))

slice_dir   <- file.path(base_dir, "neocortex_new")
slice_files <- list.files(slice_dir, pattern = "\\.txt$", full.names = TRUE)

get_layer_token_from_ccf <- function(ccf){
  ccf_clean <- str_squish(as.character(ccf))
  tok <- tolower(str_match(
    ccf_clean,
    regex("(?:,\\s*|/\\s*|\\s+)(?:layer\\s*)?([0-9]+(?:/[0-9]+)?[ab]?)\\s*$", ignore_case = TRUE)
  )[,2])
  tok <- ifelse(tok %in% layer_levels, tok, NA_character_)
  tok
}

read_slice_neuron <- function(fpath){
  dt <- data.table::fread(fpath, sep = "\t", header = TRUE, data.table = TRUE, check.names = FALSE, showProgress = FALSE)
  chip_id <- sub("\\.txt$", "", basename(fpath))
  
  if ("cell_label" %in% names(dt)) {
    dt[, cell_id := as.character(cell_label)]
  } else if ("cell_id" %in% names(dt)) {
    dt[, cell_id := as.character(cell_id)]
  } else {
    return(NULL)
  }
  
  if (!"subclass" %in% names(dt)) return(NULL)
  dt[, subclass := as.character(subclass)]
  
  if ("layer" %in% names(dt)) {
    dt[, layer_token := tolower(str_replace_all(as.character(layer), "^\\s*layer\\s*", ""))]
  } else if ("ccf_region_name" %in% names(dt)) {
    dt[, layer_token := get_layer_token_from_ccf(ccf_region_name)]
  } else {
    dt[, layer_token := NA_character_]
  }
  
  dt[, chip_id := chip_id]
  dt[, cell_uid := paste0(chip_id, "|", cell_id)]
  dt <- dt[, .(cell_uid, chip_id, layer_token, subclass)]
  dt
}

dt_cells <- rbindlist(lapply(slice_files, read_slice_neuron), use.names = TRUE, fill = TRUE)
dt_cells <- unique(dt_cells, by = "cell_uid")

dt_cells <- merge(dt_cells, Merfish_sub_map_dt, by = "subclass", all.x = TRUE)
dt_cells <- dt_cells[cell_Neruon_type %chin% neuron_types]
dt_cells <- dt_cells[layer_token %chin% layer_levels]
dt_cells[, layer_plot := factor(layer_token, levels = layer_levels)]

dt_cells[, in_cluster := cell_uid %chin% cluster_cell_uid]

total_neurons_all   <- nrow(dt_cells)
cluster_neurons_all <- dt_cells[, sum(in_cluster)]
noncluster_all      <- total_neurons_all - cluster_neurons_all

pie_dt <- data.table(
  Class = c("Cluster", "Non-cluster"),
  n     = c(cluster_neurons_all, noncluster_all)
)
pie_dt[, pct := n / sum(n)]
pie_dt[, label := sprintf("%.1f%%", pct * 100)]

g_cols <- c("Cluster" = "#F2A94F", "Non-cluster" = "#6FB06F")

p_g <- ggplot(pie_dt, aes(x = 2, y = n, fill = Class)) +
  geom_col(width = 1, color = "white", linewidth = 0.8) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5),
            size = 5, fontface = "bold", color = "black") +
  scale_fill_manual(values = g_cols, name = "Class") +
  labs(title = paste0("Neuronal cluster\n(", comma(total_neurons_all), " cells)")) +
  theme_void(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(size = 12)
  )

dt_chip_layer <- dt_cells[, .(
  total_neurons   = .N,
  cluster_neurons = sum(in_cluster),
  prop_cluster    = sum(in_cluster) / .N
), by = .(chip_id, layer_plot)]

dt_layer <- dt_chip_layer[, .(
  prop_mean = mean(prop_cluster, na.rm = TRUE),
  prop_sem  = sd(prop_cluster, na.rm = TRUE) / sqrt(sum(!is.na(prop_cluster))),
  total_neurons_allchips   = sum(total_neurons, na.rm = TRUE),
  cluster_neurons_allchips = sum(cluster_neurons, na.rm = TRUE)
), by = layer_plot]

dt_layer <- dt_layer[order(layer_plot)]

layer_cols <- c(
  "1"   = "#E76F51",
  "2/3" = "#4EA8DE",
  "4"   = "#6FB06F",
  "5"   = "#3A4E8C",
  "6a"  = "#F2A279",
  "6b"  = "#9AA5B1"
)

p_h <- ggplot(dt_layer, aes(x = layer_plot, y = prop_mean, fill = as.character(layer_plot))) +
  geom_col(width = 0.70, color = "black", linewidth = 0.8) +
  geom_errorbar(aes(ymin = pmax(prop_mean - prop_sem, 0), ymax = pmin(prop_mean + prop_sem, 1)),
                width = 0.18, linewidth = 0.9) +
  geom_text(aes(y = 0.03, label = comma(total_neurons_allchips)),
            angle = 90, vjust = 0.5, hjust = 0,
            size = 4.2, fontface = "bold", color = "black") +
  scale_fill_manual(values = layer_cols, guide = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0.02)) +
  labs(x = "Layer", y = "Prop. of cells in clusters") +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(linewidth = 1.0),
    axis.ticks = element_line(linewidth = 1.0),
    axis.text  = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    plot.margin = margin(6, 6, 6, 6)
  )

p_gh <- cowplot::plot_grid(p_g, p_h, ncol = 2, rel_widths = c(1, 1.35), labels = c("g", "h"),
                           label_fontface = "bold", label_size = 16)

fwrite(
  as.data.frame(dt_layer),
  file = file.path(fig_dir, "Fig1G-H_neuron_cluster_layer_stats.csv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

ggsave(file.path(fig_dir, "Fig1G_neuron_cluster_pie.png"), p_g, width = 5.2, height = 4.8, dpi = 400, bg = "white")
ggsave(file.path(fig_dir, "Fig1G_neuron_cluster_pie.pdf"), p_g, width = 5.2, height = 4.8)

ggsave(file.path(fig_dir, "Fig1H_layer_prop_in_clusters.png"), p_h, width = 6.2, height = 4.8, dpi = 400, bg = "white")
ggsave(file.path(fig_dir, "Fig1H_layer_prop_in_clusters.pdf"), p_h, width = 6.2, height = 4.8)

ggsave(file.path(fig_dir, "Fig1G-H_neuron_cluster_and_layer.png"), p_gh, width = 11.8, height = 4.8, dpi = 400, bg = "white")
Cairo::CairoPDF(file.path(fig_dir, "Fig1G-H_neuron_cluster_and_layer.pdf"), width = 11.8, height = 4.8)
print(p_gh)
dev.off()

cat("Total neurons (dt_cells): ", scales::comma(total_neurons_all), "\n")
cat("Cluster neurons:          ", scales::comma(cluster_neurons_all), "\n")
cat("Non-cluster neurons:      ", scales::comma(noncluster_all), "\n")
cat("Cluster proportion:       ", sprintf("%.4f", cluster_neurons_all / total_neurons_all), "\n")
cat("Check (cluster+non = total): ", (cluster_neurons_all + noncluster_all) == total_neurons_all, "\n\n")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(cowplot)
  library(scales)
  library(Cairo)
})

layer_levels <- c("1","2/3","4","5","6a","6b")
neuron_types <- c("Gaba","Glut")

clu_tbl <- if (exists("mouse_subclass_cluster_total_3")) mouse_subclass_cluster_total_3 else mouse_subclass_cluster_total
dt_clu0 <- as.data.table(clu_tbl)
stopifnot(exists("Merfish_sub_map"))
stopifnot(exists("base_dir"))
stopifnot(exists("fig_dir"))

Merfish_sub_map_dt <- as.data.table(Merfish_sub_map)[, .(subclass, cell_Neruon_type)]
if (!"cell_Neruon_type" %in% names(dt_clu0)) {
  dt_clu0 <- merge(dt_clu0, Merfish_sub_map_dt, by = "subclass", all.x = TRUE)
}
dt_clu0 <- dt_clu0[cell_Neruon_type %chin% neuron_types]
stopifnot(all(c("slide","layer","enrich_subclass_cell_ids") %in% names(dt_clu0)))

get_layer_token_from_ccf <- function(ccf){
  ccf_clean <- str_squish(as.character(ccf))
  tok <- tolower(str_match(
    ccf_clean,
    regex("(?:,\\s*|/\\s*|\\s+)(?:layer\\s*)?([0-9]+(?:/[0-9]+)?[ab]?)\\s*$", ignore_case = TRUE)
  )[,2])
  tok <- ifelse(tok %in% layer_levels, tok, NA_character_)
  tok
}

get_region_from_ccf <- function(ccf){
  x <- str_squish(as.character(ccf))
  if (all(is.na(x))) return(rep(NA_character_, length(x)))
  x2 <- str_trim(str_replace(
    x,
    regex("(?:,\\s*|/\\s*|\\s+)(?:layer\\s*)?[0-9]+(?:/[0-9]+)?[ab]?\\s*$", ignore_case = TRUE),
    ""
  ))
  reg <- str_trim(sub(",.*$", "", x2))
  reg <- str_squish(tolower(reg))
  reg
}

slice_dir   <- file.path(base_dir, "neocortex_new")
slice_files <- list.files(slice_dir, pattern = "\\.txt$", full.names = TRUE)

read_slice_map <- function(fpath){
  hdr <- names(data.table::fread(fpath, nrows = 0, data.table = FALSE, check.names = FALSE))
  sel <- intersect(hdr, c("cell_label","cell_id","subclass","ccf_region_name","layer"))
  dt  <- data.table::fread(fpath, sep = "\t", header = TRUE, data.table = TRUE, check.names = FALSE,
                           select = sel, showProgress = FALSE)
  chip_id <- sub("\\.txt$", "", basename(fpath))
  if ("cell_label" %in% names(dt)) dt[, cell_id := as.character(cell_label)]
  if (!"cell_id" %in% names(dt)) return(NULL)
  dt[, cell_id := as.character(cell_id)]
  if (!"subclass" %in% names(dt)) return(NULL)
  dt[, subclass := as.character(subclass)]
  if ("layer" %in% names(dt)) {
    dt[, layer_token := tolower(str_replace_all(as.character(layer), "^\\s*layer\\s*", ""))]
  } else if ("ccf_region_name" %in% names(dt)) {
    dt[, layer_token := get_layer_token_from_ccf(ccf_region_name)]
  } else {
    dt[, layer_token := NA_character_]
  }
  if ("ccf_region_name" %in% names(dt)) {
    dt[, region_key := get_region_from_ccf(ccf_region_name)]
  } else {
    dt[, region_key := NA_character_]
  }
  dt[, chip_id := chip_id]
  dt[, cell_uid := paste0(chip_id, "|", cell_id)]
  dt <- dt[, .(cell_uid, chip_id, layer_token, region_key, subclass)]
  dt
}

dt_map <- rbindlist(lapply(slice_files, read_slice_map), use.names = TRUE, fill = TRUE)
dt_map <- unique(dt_map, by = "cell_uid")
dt_map <- merge(dt_map, Merfish_sub_map_dt, by = "subclass", all.x = TRUE)
dt_map <- dt_map[cell_Neruon_type %chin% neuron_types]
dt_map <- dt_map[layer_token %chin% layer_levels]
dt_map <- dt_map[!is.na(region_key) & region_key != ""]

cell2region <- dt_map$region_key
names(cell2region) <- dt_map$cell_uid

split_ids <- function(x){
  if (is.null(x) || length(x) == 0 || is.na(x) || x == "") return(character(0))
  ids <- unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE)
  ids <- ids[ids != "" & !is.na(ids)]
  ids
}

dt_clu <- copy(dt_clu0)
dt_clu[, layer_token := tolower(str_replace_all(as.character(layer), "^\\s*layer\\s*", ""))]
dt_clu <- dt_clu[layer_token %chin% layer_levels]
dt_clu[, cluster_type := fifelse(cell_Neruon_type == "Glut", "Excitatory neuron cluster", "Inhibitory neuron cluster")]

dt_clu[, cell_uid_vec := Map(function(x, s){
  ids <- split_ids(x)
  if (!length(ids)) return(character(0))
  if (any(grepl("\\|", ids))) ids else paste0(s, "|", ids)
}, enrich_subclass_cell_ids, slide)]

dt_clu[, is_interregional := vapply(cell_uid_vec, function(uids){
  if (!length(uids)) return(FALSE)
  regs <- cell2region[uids]
  regs <- unique(regs[!is.na(regs) & regs != ""])
  length(regs) >= 2L
}, logical(1))]

dt_clu[, layer_plot := factor(layer_token, levels = layer_levels)]

dt_tot_type <- dt_clu[, .(
  total_clusters = .N,
  inter_clusters = sum(is_interregional)
), by = .(layer_plot, cluster_type)]

dt_tot_all <- dt_clu[, .(
  total_clusters = .N,
  inter_clusters = sum(is_interregional)
), by = layer_plot]
dt_tot_all[, prop_inter := ifelse(total_clusters == 0, NA_real_, inter_clusters / total_clusters)]

bar_cols <- c("Excitatory neuron cluster" = "#d62728", "Inhibitory neuron cluster" = "#1f77b4")

k_scale <- max(dt_tot_type$inter_clusters, na.rm = TRUE)
if (!is.finite(k_scale) || k_scale <= 0) k_scale <- 1

p_i <- ggplot() +
  geom_col(
    data = dt_tot_type,
    aes(x = layer_plot, y = inter_clusters, fill = cluster_type),
    width = 0.72, color = "black", linewidth = 0.35,
    position = position_dodge(width = 0.78)
  ) +
  geom_line(
    data = dt_tot_all,
    aes(x = layer_plot, y = prop_inter * k_scale, group = 1),
    linewidth = 1.1, color = "black"
  ) +
  geom_point(
    data = dt_tot_all,
    aes(x = layer_plot, y = prop_inter * k_scale),
    size = 2.6, color = "black"
  ) +
  scale_fill_manual(values = bar_cols, name = NULL) +
  scale_y_continuous(
    name = "Number of interregional clusters",
    sec.axis = sec_axis(~ . / k_scale, name = "Prop. of interregional clusters"),
    expand = c(0, 0.02)
  ) +
  labs(x = "Layer") +
  theme_classic(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.title.y.right = element_text(size = 12, face = "bold"),
    axis.title.y.left  = element_text(size = 12, face = "bold"),
    axis.title.x       = element_text(size = 12, face = "bold"),
    axis.text          = element_text(size = 11, face = "bold"),
    legend.position    = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background  = element_rect(fill = "white", color = "black", linewidth = 0.3),
    legend.text        = element_text(size = 10, face = "bold")
  )

fisher_or_ci <- function(a, b, c, d){
  m <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  ft <- fisher.test(m)
  or <- as.numeric(ft$estimate)
  ci <- as.numeric(ft$conf.int)
  p  <- as.numeric(ft$p.value)
  list(or = or, lo = ci[1], hi = ci[2], p = p)
}

dt_j_base <- dt_clu[, .(
  total_clusters = .N,
  inter_clusters = sum(is_interregional)
), by = .(cluster_type, layer_plot)]

dt_j_list <- list()
for (ct in unique(dt_j_base$cluster_type)) {
  tmp <- dt_j_base[cluster_type == ct]
  for (L in layer_levels) {
    a <- tmp[layer_plot == L, inter_clusters]
    b <- tmp[layer_plot == L, total_clusters - inter_clusters]
    if (length(a) == 0) a <- 0
    if (length(b) == 0) b <- 0
    c <- tmp[layer_plot != L, sum(inter_clusters)]
    d <- tmp[layer_plot != L, sum(total_clusters - inter_clusters)]
    st <- fisher_or_ci(a, b, c, d)
    dt_j_list[[length(dt_j_list) + 1L]] <- data.table(
      cluster_type = ct,
      layer_plot   = factor(L, levels = layer_levels),
      a = a, b = b, c = c, d = d,
      odds_ratio = st$or,
      ci_low = st$lo,
      ci_high = st$hi,
      p_value = st$p
    )
  }
}
dt_j <- rbindlist(dt_j_list, use.names = TRUE, fill = TRUE)
dt_j[, layer_lab := paste0("I", as.character(layer_plot))]

layer_cols <- c(
  "1"   = "#E76F51",
  "2/3" = "#4EA8DE",
  "4"   = "#6FB06F",
  "5"   = "#3A4E8C",
  "6a"  = "#F2A279",
  "6b"  = "#9AA5B1"
)

ymax_j <- max(dt_j$ci_high[is.finite(dt_j$ci_high)], na.rm = TRUE)
if (!is.finite(ymax_j) || ymax_j < 2.5) ymax_j <- 2.5
ymax_j <- min(ymax_j, 5)

p_j <- ggplot(dt_j, aes(x = layer_plot, y = odds_ratio, color = as.character(layer_plot))) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.9) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.18, linewidth = 0.9) +
  geom_point(size = 2.8) +
  facet_wrap(~ cluster_type, ncol = 1) +
  scale_color_manual(values = layer_cols, guide = "none") +
  scale_x_discrete(labels = function(x) paste0("I", x)) +
  coord_cartesian(ylim = c(0, ymax_j)) +
  labs(x = "Layer", y = "Odds ratio") +
  theme_classic(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 11, face = "bold")
  )

fwrite(
  as.data.frame(dt_tot_type),
  file = file.path(fig_dir, "Fig1I_interregional_clusters_by_layer_type.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
fwrite(
  as.data.frame(dt_j),
  file = file.path(fig_dir, "Fig1J_interregional_enrichment_OR.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

ggsave(file.path(fig_dir, "Fig1I_interregional_clusters_barline.png"), p_i, width = 9.3, height = 6.9, dpi = 400, bg = "white")
ggsave(file.path(fig_dir, "Fig1I_interregional_clusters_barline.pdf"), p_i, width = 9.3, height = 6.9)

ggsave(file.path(fig_dir, "Fig1J_interregional_enrichment_OR.png"), p_j, width = 6.2, height = 4.9, dpi = 400, bg = "white")
ggsave(file.path(fig_dir, "Fig1J_interregional_enrichment_OR.pdf"), p_j, width = 6.2, height = 4.9)

p_ij <- cowplot::plot_grid(p_i, p_j, ncol = 2, rel_widths = c(1.05, 1.0), labels = c("i","j"),
                           label_fontface = "bold", label_size = 16)
ggsave(file.path(fig_dir, "Fig1I-J_interregional_combined.png"), p_ij, width = 12.8, height = 4.9, dpi = 400, bg = "white")
Cairo::CairoPDF(file.path(fig_dir, "Fig1I-J_interregional_combined.pdf"), width = 12.8, height = 4.9)
print(p_ij)
dev.off()
