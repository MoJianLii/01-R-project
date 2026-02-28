rm(list = ls()); gc()
setwd('E:/zaw/2511/mouseMerfish_zhuang_subclass')
outdir <- "figures2_subclass"
dir.create(outdir, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(scales)
  library(ggrepel)
  library(ggtext)
  library(ggpubr)
})

tab <- read.delim('./mouse_table/table2/Table_2_total_cell.txt')

tab$pair_type <- paste(tab$cluster.1_cell_Neuron_type,
                       tab$cluster.2_cell_Neuron_type, sep = '-')
tab$pair_type[tab$pair_type == 'Gaba-Glut']      <- 'Glut-Gaba'
tab$pair_type[tab$pair_type == 'NonNeuron-Glut'] <- 'Glut-NonNeuron'
tab$pair_type[tab$pair_type == 'NonNeuron-Gaba'] <- 'Gaba-NonNeuron'

standardize_pair_type <- function(pair_str){
  parts <- strsplit(pair_str, " vs ")[[1]]
  paste(sort(parts), collapse = " vs ")
}
tab$subclass_pair_type <- paste(tab$cluster.1_subclass, tab$cluster.2_subclass, sep = " vs ")
tab$subclass_pair_type <- sapply(tab$subclass_pair_type, standardize_pair_type)

pair_keep <- c("Gaba-Gaba","Glut-Gaba","Glut-Glut")
dat_comp  <- tab %>%
  filter(pair_type %in% pair_keep) %>%
  count(pair_type, name = "n") %>%
  mutate(pct = n / sum(n),
         lbl = paste0(pair_type, "\n", n, " (", percent(pct, accuracy = 0.1), ")")) %>%
  arrange(match(pair_type, pair_keep))

dat_top <- tab %>%
  filter(pair_type %in% pair_keep) %>%
  count(pair_type, subclass_pair_type, name = "n") %>%
  arrange(desc(n))

# ===== 3) 颜色 + 主题 =====
pal_pair <- c(
  "Gaba-Gaba" = "#0072B2",
  "Glut-Gaba" = "#009E73",
  "Glut-Glut" = "#D55E00"
)

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

# =========================
# Fig2A 左
# =========================
total_n <- sum(dat_comp$n)
dat_comp <- dat_comp %>%
  dplyr::mutate(
    pct = n / total_n,
    lbl = sprintf("%s\n%.2f%%", pair_type, pct * 100)
  )

pA_left <- ggplot(dat_comp, aes(x = 2, y = n, fill = pair_type)) +
  geom_col(width = 0.8, color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.6) +
  
  geom_text(
    aes(x = 1.9, y = n, label = lbl),       # 显式提供 y=n
    position   = position_stack(vjust = 0.5),
    color      = "white",
    fontface   = "bold",
    size       = 3.2,
    lineheight = 0.95
  ) +
  annotate(
    "text", x = 0.5, y = 0,
    label    = paste0("Total\n", total_n),
    fontface = "bold",
    size     = 4.1,
    lineheight = 1.05
  ) +
  scale_fill_manual(
    values = pal_pair,
    breaks = c("Gaba-Gaba", "Glut-Gaba", "Glut-Glut")
  ) +
  labs(title = "Pair Type Composition", fill = NULL) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top"
  )


ggsave(file.path(outdir, "Fig2A_left_A1_donut.pdf"),
       pA_left, width = 6, height = 5, device = cairo_pdf)
ggsave(file.path(outdir, "Fig2A_left_A1_donut.png"),
       pA_left, width = 6, height = 5, dpi = 450)

# A2. 100% 堆叠条形
pA_left_100bar <- ggplot(dat_comp %>% mutate(pair_type=factor(pair_type, levels=pair_keep)),
                         aes(x = "All pairs", y = pct, fill = pair_type)) +
  geom_col(width = 0.5, color = "white") +
  geom_text(aes(label = percent(pct, accuracy = 0.1)),
            position = position_stack(vjust = 0.5),
            size = 4.5, color = "white", fontface = "bold") +
  scale_y_continuous(labels = percent_format(), expand = c(0,0), limits = c(0,1)) +
  scale_fill_manual(values = pal_pair) +
  theme_topjournal(12) +
  theme(axis.title = element_blank(), axis.text.x = element_blank(),
        panel.grid.major.y = element_line(color="grey90")) +
  labs(title = "Pair Type Composition (100%)")

ggsave(file.path(outdir, "Fig2A_left_A2_100bar.pdf"), pA_left_100bar, width = 5, height = 5, device = cairo_pdf)
ggsave(file.path(outdir, "Fig2A_left_A2_100bar.png"), pA_left_100bar, width = 5, height = 5, dpi = 450)

# =========================
# Fig2A 右 Top10
# =========================
topN <- 15L
pA_right_bar <- ggplot(dat_top %>%
                         slice_max(n, n=topN) %>%
                         mutate(subclass_pair_type=gsub("^\\d+\\s+", "", subclass_pair_type),
                                subclass_pair_type=gsub("vs\\s*\\d+\\s+", "vs ", subclass_pair_type),
                                subclass_pair_type=forcats::fct_reorder(subclass_pair_type, n),
                                pair_type=factor(pair_type, levels=pair_keep)),
                       aes(x=n, y=subclass_pair_type, fill=pair_type)) +
  geom_col(width=0.65) +
  geom_text(aes(label=scales::comma(n)), hjust=-0.15, size=4.2, fontface="bold") +
  scale_fill_manual(values=pal_pair) +
  scale_x_continuous(expand=expansion(mult=c(0,0.08))) +
  labs(x="Cluster number", y=NULL, title="Mouse subclass cluster (Top 10)") +
  theme_topjournal(11) +
  theme(axis.text.y=element_text(size=12, face="bold"))

ggsave(file.path(outdir, "Fig2A_right_B1_bar.pdf"), pA_right_bar, width=10, height=5.6, device=cairo_pdf)
ggsave(file.path(outdir, "Fig2A_right_B1_bar.png"), pA_right_bar, width=10, height=5.6, dpi=450)

# B3. 分面条形
K <- 10L
pA_right_facet <- ggplot(dat_top %>%
                           group_by(pair_type) %>%
                           slice_max(n, n=K, with_ties=FALSE) %>%
                           ungroup() %>%
                           mutate(subclass_pair_type=gsub("^\\d+\\s+", "", subclass_pair_type),
                                  subclass_pair_type=gsub("vs\\s*\\d+\\s+", "vs ", subclass_pair_type)) %>%
                           group_by(pair_type) %>%
                           mutate(subclass_pair_type=forcats::fct_reorder(subclass_pair_type, n)) %>%
                           ungroup(),
                         aes(x=n, y=subclass_pair_type, fill=pair_type)) +
  geom_col(width=0.62, show.legend=FALSE) +
  geom_text(aes(label=scales::comma(n)), hjust=-0.14, size=4, fontface="bold") +
  facet_wrap(~pair_type, scales="free_y", ncol=1) +
  scale_fill_manual(values=pal_pair) +
  scale_x_continuous(expand=expansion(mult=c(0,0.14))) +
  labs(x="Cluster number", y=NULL, title="Top subclass pairs by subclass") +
  theme_topjournal(11) +
  theme(strip.text=element_text(face="bold"),
        axis.text.y=element_text(size=12, face="bold"))

ggsave(file.path(outdir, "Fig2A_right_B3_facet.pdf"), pA_right_facet, width=7.2, height=10, device=cairo_pdf)
ggsave(file.path(outdir, "Fig2A_right_B3_facet.png"), pA_right_facet, width=7.2, height=10, dpi=450)

message("Done. PDFs saved to: ", normalizePath(outdir))


##################
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(scales)
})

## === 读数据（已有 tab 就注释掉） ===
tab <- read.delim('./mouse_table/table2/Table_2_total_cell.txt')

## === 1) 四类配对（保留方向） ===
tab <- tab %>%
  mutate(pair_type_raw = paste(cluster.1_cell_Neuron_type,
                               cluster.2_cell_Neuron_type, sep = "-"))

pair_levels <- c("Gaba-Gaba", "Gaba-Glut", "Glut-Gaba", "Glut-Glut")

plot_data <- tab %>%
  filter(pair_type_raw %in% pair_levels) %>%
  mutate(
    pair_type  = factor(pair_type_raw, levels = pair_levels),
    overlap_bin = cut(
      100 * cluster.2.overlap.percent,
      breaks = c(0, 20, 40, 60, 80, 100),
      include.lowest = TRUE,
      labels = c("0–20%", "20–40%", "40–60%", "60–80%", "80–100%")
    )
  )

## === 2) 统计比例 + 标签 ===
dfB <- plot_data %>%
  count(pair_type, overlap_bin, name = "n") %>%
  group_by(pair_type) %>%
  mutate(pct = n / sum(n),
         show_lab = pct >= 0.03,
         pct_lab  = paste0(round(100 * pct, 1), "%")) %>%
  ungroup()

## === 3) 颜色方案（深色 = 高重叠度 80–100%） ===
bin_pal <- c(
  "0–20%"   = "#EAE2BD",   # 浅色
  "20–40%"  = "#C7E2D5",
  "40–60%"  = "#59A8A2",
  "60–80%"  = "#2E74A6",
  "80–100%" = "#234663"    # 深色
)
lab_col <- c("0–20%"="#222222","20–40%"="#222222","40–60%"="#FFFFFF",
             "60–80%"="#FFFFFF","80–100%"="#FFFFFF")

## === 4) 主题（白底，去掉空隙） ===
theme_cleanstack <- function(base = 13){
  theme_minimal(base_size = base) %+replace%
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.6),
      axis.title         = element_blank(),
      axis.text.x        = element_text(margin = margin(t = 5), face = "bold"),
      plot.title         = element_text(face = "bold", hjust = 0.5, margin = margin(b = 6)),
      legend.position    = "right",
      legend.title       = element_blank(),
      legend.key.height  = unit(10, "pt"),
      plot.margin        = margin(6, 12, 6, 6)
    )
}

## === 5) 作图 ===
pB_right_100 <- ggplot(dfB, aes(x = pair_type, y = pct, fill = overlap_bin)) +
  geom_col(width = 0.68, color = "black", linewidth = 0.25) +
  geom_text(
    data = dplyr::filter(dfB, show_lab),
    aes(label = pct_lab, color = overlap_bin),
    position = position_stack(vjust = 0.5),
    size = 4, fontface = "bold", show.legend = FALSE
  ) +
  scale_color_manual(values = lab_col, guide = "none") +
  scale_fill_manual(values = bin_pal) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = c(0, 0), breaks = seq(0, 1, by = 0.25)) +   # 没有上方空隙
  labs(title = "Overlap-bin composition across pair types") +
  theme_cleanstack()

## === 6) 保存 ===
ggsave(file.path(outdir, "Fig_overlapbin_100stack_four_types.pdf"),
       pB_right_100, width = 10.5, height = 5.6, device = cairo_pdf)
ggsave(file.path(outdir, "Fig_overlapbin_100stack_four_types.png"),
       pB_right_100, width = 10.5, height = 5.6, dpi = 300)

pB_right_100


#
# === 分组柱状图：与示例图一致（白底、分组并排、顶部标注百分比、配色沿用当前方案） ===
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(scales)
})


# 1) 四类配对（保留方向）
tab <- tab %>%
  mutate(pair_type_raw = paste(cluster.1_cell_Neuron_type,
                               cluster.2_cell_Neuron_type, sep = "-"))

pair_levels <- c("Gaba-Gaba", "Gaba-Glut", "Glut-Gaba", "Glut-Glut")
bin_levels  <- c("0–20%", "20–40%", "40–60%", "60–80%", "80–100%")

plot_data <- tab %>%
  filter(pair_type_raw %in% pair_levels) %>%
  mutate(
    pair_type   = factor(pair_type_raw, levels = pair_levels),
    overlap_bin = cut(100 * cluster.2.overlap.percent,                # 如需 cluster.1 改这里
                      breaks = c(0, 20, 40, 60, 80, 100),
                      include.lowest = TRUE, labels = bin_levels)
  )

# 2) 统计每个 pair_type 内各分箱占比
dfB <- plot_data %>%
  count(pair_type, overlap_bin, name = "n") %>%
  group_by(pair_type) %>%
  mutate(pct = n / sum(n), pct_lab = paste0(round(100*pct, 1), "%")) %>%
  ungroup()

# 3) 配色（沿用：低→浅，高→深）
bin_pal <- c(
  "0–20%"   = "#EAE2BD",   # 浅
  "20–40%"  = "#C7E2D5",
  "40–60%"  = "#59A8A2",
  "60–80%"  = "#2E74A6",
  "80–100%" = "#234663"    # 深
)

# 4) 主题（白底、边框、上方图例）
theme_groupbar <- function(base = 12){
  theme_classic(base_size = base) %+replace%
    theme(
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.6),
      axis.title.x      = element_blank(),
      axis.title.y      = element_text(margin = margin(r = 6)),
      axis.text.x       = element_text(face = "bold", margin = margin(t = 5)),
      legend.position   = "top",
      legend.title      = element_text(face = "bold"),
      legend.key.height = unit(10, "pt"),
      plot.title        = element_text(hjust = 0.5, face = "bold", margin = margin(b = 6)),
      plot.margin       = margin(6, 10, 6, 6)
    )
}

# 5) 作图（分组并排、顶部标注百分比）
# 计算顶部留白（让标签不被裁切）
y_max <- max(dfB$pct, na.rm = TRUE)
ylim_top <- min(0.56, y_max * 1.15)  # 给一点富余，防止超过 60%

p_group <- ggplot(
  dfB, aes(x = pair_type, y = pct, fill = overlap_bin)
) +
  geom_col(
    position = position_dodge(width = 0.72),
    width = 0.62, color = "black", linewidth = 0.25
  ) +
  geom_text(
    aes(label = pct_lab),
    position = position_dodge(width = 0.72),
    vjust = -0.35, size = 3.2, fontface = "bold"
  ) +
  scale_fill_manual(values = bin_pal,
                    breaks = bin_levels,
                    name   = "Sub Cluster Cell Overlap Percent") +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, ylim_top),
    expand = c(0, 0),
    breaks = seq(0, 0.5, by = 0.1)
  ) +
  labs(y = NULL, title = NULL) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, title.position = "top")) +
  theme_groupbar(12)

# 6) 保存（确保 outdir 已定义）
# outdir <- "figures2_final"; dir.create(outdir, showWarnings = FALSE)
ggsave(file.path(outdir, "Fig_overlapbin_grouped_four_types.pdf"),
       p_group, width = 10, height = 5, device = cairo_pdf)
ggsave(file.path(outdir, "Fig_overlapbin_grouped_four_types.png"),
       p_group, width = 10, height = 5, dpi = 300)

p_group

######################
# =========================
# Fig 2C —— 三种顶刊级方案（散点 + 回归 / 边际分布 / 二维密度）
# =========================
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(ggExtra)   # 方案二：边际分布
  library(viridis)   # 方案三：二维密度
})

# ---- 预处理：只取 Glut/Gaba 互配；可切换以 cluster2 为“子集” ----
SUB_IS_CLUSTER2 <- TRUE

dfC_raw <- tab %>%
  filter(
    `cluster.1_cell_Neuron_type` %in% c("Glut","Gaba"),
    `cluster.2_cell_Neuron_type` %in% c("Glut","Gaba"),
    `cluster.1_cell_Neuron_type` != `cluster.2_cell_Neuron_type`
  )

if (SUB_IS_CLUSTER2) {
  dfC <- dfC_raw %>%
    mutate(
      overlap_sub_pct = 100 * `cluster.2.overlap.percent`,
      log2_size_sub  = log2(pmax(`cluster.2_total_cell_num`, 1)),
      log2_size_main = log2(pmax(`cluster.1_total_cell_num`, 1)),
      log2_ei_sub    = log2(pmax(`cluster.2_E_I_Ratio`, 1e-6)),
      log2_ei_main   = log2(pmax(`cluster.1_E_I_Ratio`, 1e-6))
    )
} else {
  dfC <- dfC_raw %>%
    mutate(
      overlap_sub_pct = 100 * `cluster.1.overlap.percent`,
      log2_size_sub  = log2(pmax(`cluster.1_total_cell_num`, 1)),
      log2_size_main = log2(pmax(`cluster.2_total_cell_num`, 1)),
      log2_ei_sub    = log2(pmax(`cluster.1_E_I_Ratio`, 1e-6)),
      log2_ei_main   = log2(pmax(`cluster.2_E_I_Ratio`, 1e-6))
    )
}

bin_levels <- c("0–20%","20–40%","40–60%","60–80%","80–100%")
dfC$overlap_bin <- cut(
  dfC$overlap_sub_pct, breaks = c(0,20,40,60,80,100),
  include.lowest = TRUE, labels = bin_levels
)

# 与 B 图一致：低→浅，高→深
c_cols <- c("0–20%"="#EAE2BD","20–40%"="#C7E2D5","40–60%"="#59A8A2",
            "60–80%"="#2E74A6","80–100%"="#234663")

# 统一坐标范围（便于跨 bin 对比）
size_limits <- c(2.5, 12.5); size_breaks <- seq(2.5, 12.5, 2.5)
eir_limits  <- c(-2.5, 7.5);  eir_breaks  <- seq(-2.5, 7.5, 2.5)

# ---- 方案一：经典顶刊风格（散点 + 线性回归 + 显著性） ----
make_scatter_regline <- function(dsub, x, y, xlabel, ylabel, title, col,
                                 xlim, ylim, xbr, ybr){
  ggplot(dsub, aes({{x}}, {{y}})) +
    geom_point(alpha = 0.35, size = 0.8, color = col) +
    geom_smooth(method = "lm", se = TRUE, color = "black",
                linetype = 2, linewidth = 0.7) +
    stat_cor(method = "pearson",
             label.x.npc = "left", label.y.npc = "top", size = 4) +
    labs(x = xlabel, y = ylabel, title = title) +
    scale_x_continuous(limits = xlim, breaks = xbr) +
    scale_y_continuous(limits = ylim, breaks = ybr) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
    )
}

# ---- 方案二：边际分布（直方图）+ 散点 + 回归 ----
make_scatter_marginal <- function(dsub, x, y, xlabel, ylabel, title, col,
                                  xlim, ylim, xbr, ybr){
  p <- ggplot(dsub, aes({{x}}, {{y}})) +
    geom_point(alpha = 0.35, size = 0.8, color = col) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.7) +
    stat_cor(method = "pearson",
             label.x.npc = "left", label.y.npc = "top", size = 4) +
    labs(x = xlabel, y = ylabel, title = title) +
    scale_x_continuous(limits = xlim, breaks = xbr) +
    scale_y_continuous(limits = ylim, breaks = ybr) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  ggMarginal(p, type = "histogram", fill = col, alpha = 0.5, bins = 30)
}

# ---- 方案三：二维密度（高点数/拥挤时最清晰） ----
make_scatter_density <- function(dsub, x, y, xlabel, ylabel, title, col,
                                 xlim, ylim, xbr, ybr){
  ggplot(dsub, aes({{x}}, {{y}})) +
    stat_density_2d(aes(fill = ..level..), geom = "polygon",
                    color = NA, alpha = 0.8) +
    scale_fill_viridis(option = "C", direction = -1, guide = "none") +
    geom_point(alpha = 0.15, size = 0.5, color = col) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.7) +
    stat_cor(method = "pearson",
             label.x.npc = "left", label.y.npc = "top", size = 4) +
    labs(x = xlabel, y = ylabel, title = title) +
    scale_x_continuous(limits = xlim, breaks = xbr) +
    scale_y_continuous(limits = ylim, breaks = ybr) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
    )
}

# ---- 生成三套图（每个 bin 两张：Size & E/I），并拼成 2×5 面板 ----
build_grid_and_save <- function(fun_builder, tag){
  plist_size <- list(); plist_eir <- list()
  for (b in bin_levels) {
    dsub <- dfC %>% filter(overlap_bin == b)
    colb <- c_cols[b]
    
    p_size <- fun_builder(dsub,
                          log2_size_sub, log2_size_main,
                          "log2 Main Cluster Size", "log2 Sub Cluster Size",
                          b, colb, size_limits, size_limits,
                          size_breaks, size_breaks)
    p_eir  <- fun_builder(dsub,
                          log2_ei_main,  log2_ei_sub,
                          "log2 Main E–I Ratio", "log2 Sub E–I Ratio",
                          b, colb, eir_limits, eir_limits,
                          eir_breaks, eir_breaks)
    plist_size[[b]] <- p_size
    plist_eir[[b]]  <- p_eir
  }
  
  grid <- ggarrange(
    plotlist = c(plist_size, plist_eir),
    ncol = length(bin_levels), nrow = 2,
    labels = c(rep("Cluster Size", length(bin_levels)),
               rep("E–I Ratio",   length(bin_levels))),
    font.label = list(size = 10, face = "bold"),
    hjust = -0.2, vjust = 1.2
  )
  
  ggsave(file.path(outdir, paste0("Fig2C_", tag, "_grid.png")),
         grid, width = 14, height = 6.2, dpi = 450)
  ggsave(file.path(outdir, paste0("Fig2C_", tag, "_grid.pdf")),
         grid, width = 14, height = 6.2)
}

# 分别输出三种方案
build_grid_and_save(make_scatter_regline,  "regline")
build_grid_and_save(make_scatter_marginal, "marginal")
build_grid_and_save(make_scatter_density,  "density")

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(scales); library(ggpubr)
})

# ===== 读入与准备 =====
tab <- read.delim('./mouse_table/table2/Table_2_total_cell.txt')

# 与原图一致：四类、以 cluster.2.overlap.percent × 100
df_all <- tab %>%
  mutate(
    pair_type = paste(cluster.1_cell_Neuron_type, cluster.2_cell_Neuron_type, sep = "-"),
    overlap2_pct = 100 * cluster.2.overlap.percent
  ) %>%
  filter(pair_type %in% c("Gaba-Gaba","Gaba-Glut","Glut-Gaba","Glut-Glut")) %>%
  mutate(pair_type = factor(pair_type, levels = c("Gaba-Gaba","Gaba-Glut","Glut-Gaba","Glut-Glut")))

pair_pal4 <- c("Gaba-Gaba"="#0072B2","Gaba-Glut"="#56B4E9","Glut-Gaba"="#009E73","Glut-Glut"="#D55E00")

# 组内数值（四舍五入到 1 位）
summary_tbl <- df_all %>%
  group_by(pair_type) %>%
  summarise(
    n     = n(),
    med   = median(overlap2_pct, na.rm = TRUE),
    meanv = mean(overlap2_pct,    na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(lab = paste0("n = ", scales::comma(n),
                      "\nmedian = ", round(med, 1), "%\nmean = ", round(meanv, 1), "%"))

# —— 统一放置高度：把“数值标签”放在上三分之一区域，括号在更靠上的区域，互不遮挡 —— 
y_max   <- max(df_all$overlap2_pct, na.rm = TRUE)
y_lbl   <- y_max * 0.82     # 组内 n/median/mean 放这里
y_brk_1 <- y_max * 0.94     # 最低一层括号起点
step_increase <- 0.065      # 括号层间距

# 两两比较（四类）：显示“数值 p 值”
comparisons_all <- list(
  c("Glut-Gaba","Gaba-Glut"),
  c("Glut-Gaba","Glut-Glut"),
  c("Glut-Gaba","Gaba-Gaba"),
  c("Gaba-Glut","Glut-Glut"),
  c("Gaba-Glut","Gaba-Gaba"),
  c("Gaba-Gaba","Glut-Glut")
)

# ===== 图 1：四类全量，带数值且避免遮挡 =====
p_full <- ggplot(df_all, aes(x = pair_type, y = overlap2_pct, fill = pair_type)) +
  geom_violin(width = 0.85, alpha = 0.15, color = NA, trim = TRUE) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.95, color = "black") +
  
  # 组内 n/median/mean（白底小标签，避免与底图混淆）
  geom_label(
    data = transform(summary_tbl, y = y_lbl),
    aes(x = pair_type, y = y, label = lab),
    inherit.aes = FALSE,
    size = 3.8, fontface = "bold",
    label.size = 0.25, label.r = unit(2.5, "pt"),
    fill = "white", color = "black"
  ) +
  
  # 两两 Wilcoxon 数值 p 值（放到更靠上的区间）
  stat_compare_means(
    comparisons = comparisons_all, method = "wilcox.test",
    label = "p.format", tip.length = 0.01, bracket.size = 0.45,
    size = 4.1, hide.ns = TRUE,
    step.increase = step_increase,
    label.y = y_brk_1   # 首层起点，其余按 step_increase 叠加
  ) +
  
  scale_fill_manual(values = pair_pal4, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.28))) +  # 顶部给足空白
  labs(title = "Overlap distribution with pairwise Wilcoxon tests",
       x = NULL, y = "Cluster2 overlap (%)") +
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    axis.text.x  = element_text(face = "bold")
  )

ggsave(file.path(outdir, "Fig_box_overlaps_pairwise_wilcox_WITH_NUMBERS_neat.pdf"),
       p_full, width = 8.6, height = 6.9, device = cairo_pdf)
ggsave(file.path(outdir, "Fig_box_overlaps_pairwise_wilcox_WITH_NUMBERS_neat.png"),
       p_full, width = 8.6, height = 6.9, dpi = 450)

# ===== 图 2：只保留 Glut-Gaba 与 Gaba-Glut =====
df_2dir <- df_all %>%
  filter(pair_type %in% c("Glut-Gaba","Gaba-Glut")) %>%
  droplevels()

summary_2dir <- df_2dir %>%
  group_by(pair_type) %>%
  summarise(
    n     = n(),
    med   = median(overlap2_pct, na.rm = TRUE),
    meanv = mean(overlap2_pct,    na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(lab = paste0("n = ", scales::comma(n),
                      "\nmedian = ", round(med, 1), "%\nmean = ", round(meanv, 1), "%"))

y_max2 <- max(df_2dir$overlap2_pct, na.rm = TRUE)
y_lbl2 <- y_max2 * 0.80
y_brk2 <- y_max2 * 0.92

p_2dir <- ggplot(df_2dir, aes(x = pair_type, y = overlap2_pct, fill = pair_type)) +
  geom_violin(width = 0.85, alpha = 0.15, color = NA, trim = TRUE) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.95, color = "black") +
  geom_label(
    data = transform(summary_2dir, y = y_lbl2),
    aes(x = pair_type, y = y, label = lab),
    inherit.aes = FALSE,
    size = 4.0, fontface = "bold",
    label.size = 0.25, label.r = unit(2.5, "pt"),
    fill = "white", color = "black"
  ) +
  stat_compare_means(
    comparisons = list(c("Glut-Gaba","Gaba-Glut")),
    method = "wilcox.test", label = "p.format",
    tip.length = 0.01, bracket.size = 0.5, size = 4.6,
    hide.ns = TRUE, label.y = y_brk2
  ) +
  scale_fill_manual(values = pair_pal4[c("Gaba-Glut","Glut-Gaba")], guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.22))) +
  labs(title = "Overlap: Glut→Gaba vs Gaba→Glut (Wilcoxon)",
       x = NULL, y = "Cluster2 overlap (%)") +
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    axis.text.x  = element_text(face = "bold")
  )

ggsave(file.path(outdir, "Fig_box_overlaps_GlutGaba_vs_GabaGlut_WITH_NUMBERS.pdf"),
       p_2dir, width = 6.6, height = 6.2, device = cairo_pdf)
ggsave(file.path(outdir, "Fig_box_overlaps_GlutGaba_vs_GabaGlut_WITH_NUMBERS.png"),
       p_2dir, width = 6.6, height = 6.2, dpi = 450)

message("✅ Saved:\n - ", normalizePath(file.path(outdir, "Fig_box_overlaps_pairwise_wilcox_WITH_NUMBERS_neat.*")),
        "\n - ", normalizePath(file.path(outdir, "Fig_box_overlaps_GlutGaba_vs_GabaGlut_WITH_NUMBERS.*")))
