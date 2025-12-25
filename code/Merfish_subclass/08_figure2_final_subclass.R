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
    aes(x = 1.9, y = n, label = lbl),
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


suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(scales)
})

tab <- read.delim('./mouse_table/table2/Table_2_total_cell.txt')
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

dfB <- plot_data %>%
  count(pair_type, overlap_bin, name = "n") %>%
  group_by(pair_type) %>%
  mutate(pct = n / sum(n),
         show_lab = pct >= 0.03,
         pct_lab  = paste0(round(100 * pct, 1), "%")) %>%
  ungroup()

bin_pal <- c(
  "0–20%"   = "#EAE2BD",   # 浅色
  "20–40%"  = "#C7E2D5",
  "40–60%"  = "#59A8A2",
  "60–80%"  = "#2E74A6",
  "80–100%" = "#234663"    # 深色
)
lab_col <- c("0–20%"="#222222","20–40%"="#222222","40–60%"="#FFFFFF",
             "60–80%"="#FFFFFF","80–100%"="#FFFFFF")

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
                     expand = c(0, 0), breaks = seq(0, 1, by = 0.25)) +
  labs(title = "Overlap-bin composition across pair types") +
  theme_cleanstack()

ggsave(file.path(outdir, "Fig_overlapbin_100stack_four_types.pdf"),
       pB_right_100, width = 10.5, height = 5.6, device = cairo_pdf)
ggsave(file.path(outdir, "Fig_overlapbin_100stack_four_types.png"),
       pB_right_100, width = 10.5, height = 5.6, dpi = 300)

pB_right_100

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(scales)
})


tab <- tab %>%
  mutate(pair_type_raw = paste(cluster.1_cell_Neuron_type,
                               cluster.2_cell_Neuron_type, sep = "-"))

pair_levels <- c("Gaba-Gaba", "Gaba-Glut", "Glut-Gaba", "Glut-Glut")
bin_levels  <- c("0–20%", "20–40%", "40–60%", "60–80%", "80–100%")

plot_data <- tab %>%
  filter(pair_type_raw %in% pair_levels) %>%
  mutate(
    pair_type   = factor(pair_type_raw, levels = pair_levels),
    overlap_bin = cut(100 * cluster.2.overlap.percent,
                      breaks = c(0, 20, 40, 60, 80, 100),
                      include.lowest = TRUE, labels = bin_levels)
  )

dfB <- plot_data %>%
  count(pair_type, overlap_bin, name = "n") %>%
  group_by(pair_type) %>%
  mutate(pct = n / sum(n), pct_lab = paste0(round(100*pct, 1), "%")) %>%
  ungroup()

bin_pal <- c(
  "0–20%"   = "#EAE2BD",
  "20–40%"  = "#C7E2D5",
  "40–60%"  = "#59A8A2",
  "60–80%"  = "#2E74A6",
  "80–100%" = "#234663"
)

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

y_max <- max(dfB$pct, na.rm = TRUE)
ylim_top <- min(0.56, y_max * 1.15)

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


ggsave(file.path(outdir, "Fig_overlapbin_grouped_four_types.pdf"),
       p_group, width = 10, height = 5, device = cairo_pdf)
ggsave(file.path(outdir, "Fig_overlapbin_grouped_four_types.png"),
       p_group, width = 10, height = 5, dpi = 300)

p_group

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(ggExtra)
  library(viridis)
})

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

c_cols <- c("0–20%"="#EAE2BD","20–40%"="#C7E2D5","40–60%"="#59A8A2",
            "60–80%"="#2E74A6","80–100%"="#234663")

size_limits <- c(2.5, 12.5); size_breaks <- seq(2.5, 12.5, 2.5)
eir_limits  <- c(-2.5, 7.5);  eir_breaks  <- seq(-2.5, 7.5, 2.5)

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

build_grid_and_save(make_scatter_regline,  "regline")
build_grid_and_save(make_scatter_marginal, "marginal")
build_grid_and_save(make_scatter_density,  "density")

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(scales); library(ggpubr)
})


tab <- read.delim('./mouse_table/table2/Table_2_total_cell.txt')
df_all <- tab %>%
  mutate(
    pair_type = paste(cluster.1_cell_Neuron_type, cluster.2_cell_Neuron_type, sep = "-"),
    overlap2_pct = 100 * cluster.2.overlap.percent
  ) %>%
  filter(pair_type %in% c("Gaba-Gaba","Gaba-Glut","Glut-Gaba","Glut-Glut")) %>%
  mutate(pair_type = factor(pair_type, levels = c("Gaba-Gaba","Gaba-Glut","Glut-Gaba","Glut-Glut")))

pair_pal4 <- c("Gaba-Gaba"="#0072B2","Gaba-Glut"="#56B4E9","Glut-Gaba"="#009E73","Glut-Glut"="#D55E00")

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

y_max   <- max(df_all$overlap2_pct, na.rm = TRUE)
y_lbl   <- y_max * 0.82
y_brk_1 <- y_max * 0.94
step_increase <- 0.065

comparisons_all <- list(
  c("Glut-Gaba","Gaba-Glut"),
  c("Glut-Gaba","Glut-Glut"),
  c("Glut-Gaba","Gaba-Gaba"),
  c("Gaba-Glut","Glut-Glut"),
  c("Gaba-Glut","Gaba-Gaba"),
  c("Gaba-Gaba","Glut-Glut")
)

p_full <- ggplot(df_all, aes(x = pair_type, y = overlap2_pct, fill = pair_type)) +
  geom_violin(width = 0.85, alpha = 0.15, color = NA, trim = TRUE) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.95, color = "black") +
  
  geom_label(
    data = transform(summary_tbl, y = y_lbl),
    aes(x = pair_type, y = y, label = lab),
    inherit.aes = FALSE,
    size = 3.8, fontface = "bold",
    label.size = 0.25, label.r = unit(2.5, "pt"),
    fill = "white", color = "black"
  ) +
  
  stat_compare_means(
    comparisons = comparisons_all, method = "wilcox.test",
    label = "p.format", tip.length = 0.01, bracket.size = 0.45,
    size = 4.1, hide.ns = TRUE,
    step.increase = step_increase,
    label.y = y_brk_1
  ) +
  
  scale_fill_manual(values = pair_pal4, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.28))) +
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
