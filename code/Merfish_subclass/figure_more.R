## ================================================================
## 0. 环境准备 & 公共工具
## ================================================================
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

## ---------- 0.1 路径与 ws/ss 组合 ---------------------------------
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

## ---------- 0.2 从目录名解析 window_size / step_size ---------------
parse_ws_ss <- function(dir_name) {
  ws <- as.numeric(str_match(dir_name, "ws([0-9.]+)")[, 2])
  ss <- as.numeric(str_match(dir_name, "ss([0-9.]+)")[, 2])
  list(ws = ws, ss = ss)
}

## ---------- 0.3 自动识别 subclass 列名 -----------------------------
find_subclass_col <- function(dt) {
  cand <- c("subclasstype", "subclass_type", "Subclass", "subclass", "sub_class")
  found <- intersect(cand, names(dt))
  if (length(found) == 0) {
    stop("在 Table_1_total_cell 里找不到 subclass 相关列，请手工改下列名或代码。")
  }
  found[1]
}

## ================================================================
## 1. 合并所有 ws/ss 的 Table_1_total_cell，并统一 subclass
## ================================================================
Table1_all <- rbindlist(
  lapply(ws_ss_dirs, function(d) {
    info <- parse_ws_ss(d)
    f    <- file.path(root_dir, d, "mouse_table", "table1", "Table_1_total_cell.txt")
    if (!file.exists(f)) stop("文件不存在：", f)
    dt <- fread(f)
    dt[, window_size := info$ws]
    dt[, step_size   := info$ss]
    dt[, ws_ss       := d]
    dt[]
  }),
  use.names = TRUE, fill = TRUE
)

## ---------- 1.1 统一 subclass 名字为 PVALB / SOM / VIP / Other -----
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

## ---------- 1.2 为 ws/ss 组合做一个好看的标签 ----------------------
Table1_all[, ws_ss_label := sprintf("ws=%.2f\nss=%.3f", window_size, step_size)]
Table1_all[, ws_ss_label := factor(
  ws_ss_label,
  levels = unique(ws_ss_label[order(window_size, step_size)])
)]

## ---------- 1.3 一些通用颜色（给 P/S/V） ---------------------------
psv_cols <- c(
  PVALB = "#4C72B0",
  SOM   = "#55A868",
  VIP   = "#C44E52"
)

## （可选）检查各组合下每类窗口数量
Table1_all[ , .N, by = .(ws_ss, subclass_group)][order(ws_ss, subclass_group)]

## ================================================================
## 2. 图2（全局）：不同尺度 (ws/ss) 下 P/S/V 集群数量和细胞总数
##    —— 只画一次，统一放在 figure_more 根目录
## ================================================================
psv_all <- Table1_all[subclass_group %in% c("PVALB", "SOM", "VIP")]

## ---------- 2.1 计算不同尺度下的指标 -------------------------------
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

## ---------- 2.2 小工具：画单个 subclass 的曲线图 --------------------
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

## ---------- 2.3 六个子图 -------------------------------------------
p_PVALB_hot  <- plot_psv_metric("PVALB", "n_hotspots")
p_PVALB_cell <- plot_psv_metric("PVALB", "sum_cells")
p_SOM_hot    <- plot_psv_metric("SOM",   "n_hotspots")
p_SOM_cell   <- plot_psv_metric("SOM",   "sum_cells")
p_VIP_hot    <- plot_psv_metric("VIP",   "n_hotspots")
p_VIP_cell   <- plot_psv_metric("VIP",   "sum_cells")

## --------- 2.4 提取统一图例 & 合并成 3×2 panel ---------------------
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

## === Fig2 只需要画一次：如果文件已经存在，就不再覆盖 ===
fig2_png <- file.path(fig_dir, "02_Fig_PSV_hotspot_scale_6panel_with_legend_right.png")
if (!file.exists(fig2_png)) {
  ggsave(fig2_png,
         p2_combined, width = 14, height = 9, dpi = 400)
}

## ================================================================
## 3. 选择 main scale，并为其建立专属子目录
##    后续 Fig1 / Fig3 / Fig4 / Fig5 全部写入该子目录
## ================================================================
## 你只需修改这里的 main_ws / main_ss，就可以切换一个新的主尺度
main_ws <- 0.4
main_ss <- 0.02

## === 针对当前 main scale 建一个专属的输出子目录 ===
main_tag     <- sprintf("ws%.2f_ss%.3f", main_ws, main_ss)
fig_dir_main <- file.path(fig_dir, main_tag)
dir.create(fig_dir_main, showWarnings = FALSE, recursive = TRUE)

## 从 Table1_all 中筛出当前 main scale 的数据
dt_main <- Table1_all[
  abs(window_size - main_ws) < 1e-6 &
    abs(step_size   - main_ss) < 1e-6
]

if (!"layer" %in% names(dt_main)) {
  stop("当前数据中没有 'layer' 列，请把真实层列名替换成 layer。")
}

## 只看 PVALB / SOM / VIP
dt_main_psv <- dt_main[subclass_group %in% c("PVALB", "SOM", "VIP")]

## 每层 × subclass 的窗口数（N）以及该层中占比 frac
layer_psv <- dt_main_psv[, .N, by = .(layer, subclass_group)]
layer_psv[, layer_total := sum(N), by = layer]
layer_psv[, frac        := N / layer_total]

## layer 按自然顺序排成因子，画图更整齐
if (!is.factor(layer_psv$layer)) {
  layer_psv[, layer := factor(layer, levels = sort(unique(layer)))]
}

## ================================================================
## 4. 图1：主尺度下 P/S/V 的层分布
## ================================================================
## ---------- 4.1 图 1A：各层 cluster 数 ------------------------------
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

## ---------- 4.2 图 1B：各层中 P/S/V 的组成比例 ----------------------
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

## ---------- 4.3 图1 合并输出 ---------------------------------------
p1_combined <- plot_grid(p1_counts, p1_frac,
                         nrow = 1, rel_widths = c(1.2, 1))
p1_combined

ggsave(file.path(fig_dir_main, "01_Fig_PSV_layer_distribution_mainScale.png"),
       p1_combined, width = 15, height = 5, dpi = 400)

## ================================================================
## 5. 图3：沿 A–P 方向的层分布 & top 区域条形图（主尺度）
## ================================================================
## ---------- 5.1 为 main scale 的 P/S/V 定义 psv_main ----------
psv_main <- dt_main_psv  ## 就是 main_ws/main_ss 且 P/S/V 的子集

print(table(psv_main$subclass_group))

## ---------- 5.2 按小鼠拆分 slide，并在每只小鼠内部排序 ------------
# slide 形如 "C57BL6J-1.025"
psv_main[, slide_suffix := tstrsplit(slide, "-", keep = 2)]
psv_main[, c("mouse_id", "section_raw") :=
           tstrsplit(slide_suffix, "\\.")]

# 每只小鼠内部按 section_raw 的数值顺序给索引（代表 A–P 方向）
psv_main[, section_index := as.numeric(section_raw)]
psv_main[is.na(section_index),
         section_index := frank(section_raw, ties.method = "dense"),
         by = mouse_id]

# layer 作为因子，纵轴使用
psv_main[, layer_f := factor(layer, levels = sort(unique(layer)))]

## ---------- 5.3 统计每只小鼠×layer×section 的 cluster 数 ----------
hot_layer_slide <- psv_main[
  ,
  .(
    n_hotspots = .N,                     # cluster 数
    sum_cells  = sum(total_cell_num)     # 对应细胞总数（备用）
  ),
  by = .(mouse_id, subclass_group, layer_f, section_index)
]

mouse_ids <- sort(unique(hot_layer_slide$mouse_id))
print(mouse_ids)

## ---------- 5.4 小工具：画单只小鼠的 layer × section 热图 ----------
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

## 单独保存每只小鼠的 3A 图（可选）
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

## 合并 4 只小鼠的图为 2×2 大图
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

## ---------- 5.5 把 region 字符串拆成多个“区域+layer”片段 ------------
split_region_smart <- function(reg_str) {
  if (is.null(reg_str) || length(reg_str) == 0L) return(character(0L))
  reg_str <- reg_str[1]
  
  toks <- strsplit(reg_str, ",")[[1]]
  toks <- trimws(toks)
  toks <- toks[toks != ""]
  if (!length(toks)) return(character(0L))
  
  # 找出 "layer ..." 的位置，并把 [start:layer] 当成一个完整片段
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

## ---------- 5.6 3B：top 区域（包含 layer 信息） ----------------------
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

## ---------- 5.7 3C：合并 layer 后的 top 区域 ------------------------
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

## ================================================================
## 6. 图4：选择一张切片，显示 top 区域中的 P/S/V 细胞（主尺度）
## ================================================================
## ---------- 6.1 选择 hotspot 数最多的切片 --------------------------
slice_to_plot <- psv_main[, .N, by = slide][order(-N)]$slide[1]
message("Slice chosen for visualization: ", slice_to_plot)

## ---------- 6.2 准备“top 区域（去 layer）”列表，用于细胞级 join -------
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

## ---------- 6.3 载入该切片的细胞坐标 -------------------------------
base_rdata_dir <- file.path(root_dir, "ws0.4_ss0.02", "cell_window")
## 如果将来 main_ws/main_ss 改动，这里也相应改成对应路径即可

rdata_path     <- file.path(base_rdata_dir, paste0(slice_to_plot, ".RData"))

message("Loading cell-level RData: ", rdata_path)
if (!file.exists(rdata_path)) {
  stop("找不到 RData 文件：", rdata_path,
       "\n请确认该切片 slide 是否有对应的 .RData，或者手动指定 RData 路径。")
}

load(rdata_path)

if (exists("file_tmp_layer")) {
  cells <- as.data.table(file_tmp_layer)
} else if (exists("file_tmp")) {
  message("对象 file_tmp_layer 不存在，自动使用 file_tmp 作为细胞坐标表。")
  cells <- as.data.table(file_tmp)
} else {
  stop("在 RData 中找不到细胞坐标表：请确认存在 file_tmp_layer 或 file_tmp。")
}

## 补齐 slide / region / subclass_group 三列
if (!"slide" %in% names(cells)) {
  if ("brain_section_label" %in% names(cells)) {
    cells[, slide := brain_section_label]
  } else {
    stop("cells 表中找不到 'slide' 或 'brain_section_label' 列，请检查 RData 结构。")
  }
}

if (!"region" %in% names(cells)) {
  if ("ccf_region_name" %in% names(cells)) {
    cells[, region := ccf_region_name]
  } else {
    stop("cells 表中找不到 'region' 或 'ccf_region_name' 列，请检查 RData 结构。")
  }
}

if (!"subclass_group" %in% names(cells)) {
  if (!"subclass" %in% names(cells)) {
    stop("cells 表中没有 'subclass_group' 也没有 'subclass'，无法识别 PVALB / SOM / VIP。")
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

## 坐标列统一为 x/y
if (!("x" %in% names(cells) && "y" %in% names(cells))) {
  x_col <- grep("(^x$|_x$|^X$|^pos_x$|^center_x$)",
                names(cells), ignore.case = TRUE, value = TRUE)[1]
  y_col <- grep("(^y$|_y$|^Y$|^pos_y$|^center_y$)",
                names(cells), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(x_col) || is.na(y_col)) {
    stop("自动没有找到坐标列，请手动把 x_col / y_col 改成你实际的坐标列名。")
  }
  setnames(cells, c(x_col, y_col), c("x", "y"))
}

## ---------- 6.4 只保留该切片的细胞，并标注是否属于 top 区域 ----------
cells_slice <- cells[slide == slice_to_plot]

cells_slice[, region_no_layer :=
              sub(",\\s*layer\\s+[^,]+$", "", region)]

cells_psv <- cells_slice[subclass_group %in% c("PVALB", "SOM", "VIP")]

setkey(cells_psv, subclass_group, region_no_layer)
cells_psv[top_regions_dt, in_top_region := TRUE]
cells_psv[is.na(in_top_region), in_top_region := FALSE]

## 背景细胞用于勾勒切片轮廓
bg_n     <- 50000L
bg_cells <- cells_slice
if (nrow(bg_cells) > bg_n) {
  set.seed(1)
  bg_cells <- bg_cells[sample(.N, bg_n)]
}

## ---------- 6.5 小工具：单独画一个 subclass 的 top 区域 ------------
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

## 合并图：三种 subclass 一起，只高亮 top 区域
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

## 2×2：三个单独图 + 合并图（无图例）
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

## ================================================================
## 7. 图5：同一切片上显示“所有 cluster 区域”（主尺度）
##       深/浅色区分 top vs 非 top
## ================================================================
## ---------- 7.1 找出“所有出现过 cluster 的解剖区域” ----------------
all_regions_dt <- unique(hot_region_nl[, .(subclass_group, region_no_layer)])
setkey(all_regions_dt, subclass_group, region_no_layer)

## 在 cells_psv 中标记：是否落在任意 cluster 区域（不只 top）
cells_psv_hot <- copy(cells_psv)
setkey(cells_psv_hot, subclass_group, region_no_layer)
cells_psv_hot[all_regions_dt, in_any_hotspot_region := TRUE]
cells_psv_hot[is.na(in_any_hotspot_region), in_any_hotspot_region := FALSE]

## 只保留落在 cluster 区域中的 P/S/V 细胞（避免太多背景点）
cells_psv_hot <- cells_psv_hot[in_any_hotspot_region == TRUE]

## ---------- 7.2 2×2 图：PVALB / SOM / VIP / ALL --------------------
## 单个 subclass：颜色固定、alpha 区分 top vs 非 top
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
      legend.position = "right",   # 每个小图单独带图例
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

p5_slice_PVALB <- plot_slice_single_F5("PVALB")
p5_slice_SOM   <- plot_slice_single_F5("SOM")
p5_slice_VIP   <- plot_slice_single_F5("VIP")

## ALL 图：简化版，颜色区分 subclass，深浅区分 top / 非 top
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

## 2×2：前三个小图带图例，第四个 ALL 去掉图例
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

## ---------- 7.3 单独的“详细图例”合并图：6 个 legend 键 -------------
## legend_group = subclass × region_type（top / non-top）
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
