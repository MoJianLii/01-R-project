rm(list=ls()); gc()

## ===================== 路径（按你项目结构） =====================
base_dir <- "E:/zaw/2511/mouseMerfish_zhuang_subclass"
ws_env <- Sys.getenv("WINDOW_SIZE", "0.4")
ss_env <- Sys.getenv("STEP_SIZE",  "0.02")

fmt <- function(x){
  s <- formatC(as.numeric(x), format="f", digits=6)
  s <- sub("0+$","",s); s <- sub("\\.$","",s)
  s
}
combo_root <- file.path(base_dir, paste0("ws", fmt(ws_env), "_ss", fmt(ss_env)))

## 你想用 over3 就用这个；否则改成 mouse_subclass_cluster_total.txt
clu_path <- file.path(combo_root, "mouse_subclass_cluster_total_over3.txt")
if (!file.exists(clu_path)) clu_path <- file.path(combo_root, "mouse_subclass_cluster_total.txt")

map_path <- file.path(combo_root, "Merfish_brain_cell_type_subclass.txt")

## ===================== 依赖 =====================
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(ggplot2)
})

## ===================== 读数据（最少用这两个） =====================
clu <- fread(clu_path, sep="\t", header=TRUE, data.table=TRUE, check.names=FALSE)
map <- fread(map_path, sep="\t", header=TRUE, data.table=TRUE, check.names=FALSE)

## mapping 文件列名兜底（避免 setnames 因列数不对报错）
if (!all(c("subclass","cell_Neruon_type") %in% names(map))) {
  # 假设前三列依次是 subclass / cell_Neruon_type / subclass_sim（你原始就是这样）
  stopifnot(ncol(map) >= 2)
  setnames(map, old = names(map)[1:2], new = c("subclass","cell_Neruon_type"))
  if (ncol(map) >= 3) setnames(map, old = names(map)[3], new = "subclass_sim")
}

## 如果 cluster 表里没有 cell_Neruon_type，就补上
if (!"cell_Neruon_type" %in% names(clu)) {
  clu <- merge(clu, map[, .(subclass, cell_Neruon_type)], by="subclass", all.x=TRUE)
}

## ===================== 统一类型命名（Gaba/Glut） =====================
clu[, cell_Neruon_type := str_trim(as.character(cell_Neruon_type))]
clu[grepl("^glu",  cell_Neruon_type, ignore.case=TRUE), cell_Neruon_type := "Glut"]
clu[grepl("^gaba", cell_Neruon_type, ignore.case=TRUE), cell_Neruon_type := "Gaba"]

## ===================== cluster size：用 total_cell_num =====================
if (!"total_cell_num" %in% names(clu)) {
  stop("clu 表里没有 total_cell_num 列；请确认 mouse_subclass_cluster_total*.txt 的列名是否一致。")
}
clu[, cluster_size := as.integer(total_cell_num)]

## 只看 Gaba / Glut
dt <- clu[cell_Neruon_type %chin% c("Gaba","Glut") &
            !is.na(cluster_size) & cluster_size > 0,
          .(cell_Neruon_type, cluster_size)]

## ===================== 汇总对比（最直接） =====================
sum_dt <- dt[, .(
  n_clusters = .N,
  mean_size  = mean(cluster_size),
  sd_size    = sd(cluster_size),
  median_size= median(cluster_size),
  iqr_size   = IQR(cluster_size),
  min_size   = min(cluster_size),
  max_size   = max(cluster_size)
), by = cell_Neruon_type]

print(sum_dt)

## ===================== 分布差异检验（Wilcoxon，稳健） =====================
x <- dt[cell_Neruon_type=="Glut", cluster_size]
y <- dt[cell_Neruon_type=="Gaba", cluster_size]

wt <- wilcox.test(x, y, alternative="two.sided", exact=FALSE)
cat("\nWilcoxon test (Glut vs Gaba):\n")
print(wt)

cat("\nMedian(Glut) - Median(Gaba) = ",
    median(x) - median(y), "\n", sep="")

## ===================== 快速可视化（log10 轴） =====================
p <- ggplot(dt, aes(x=cell_Neruon_type, y=cluster_size, fill=cell_Neruon_type)) +
  geom_boxplot(outlier.size=0.5) +
  scale_y_continuous(trans="log10") +
  theme_classic(base_size=12) +
  labs(x=NULL, y="Cluster size (#cells, log10)", title="Cluster size (total_cell_num): Glut vs Gaba")

print(p)



#############################################

rm(list=ls()); gc()

## ===================== 0) 路径与输入 =====================
base_dir <- "E:/zaw/2511/mouseMerfish_zhuang_subclass"
ws_env <- Sys.getenv("WINDOW_SIZE", "0.4")
ss_env <- Sys.getenv("STEP_SIZE",  "0.02")

fmt <- function(x){
  s <- formatC(as.numeric(x), format="f", digits=6)
  s <- sub("0+$","",s); s <- sub("\\.$","",s)
  s
}
combo_root <- file.path(base_dir, paste0("ws", fmt(ws_env), "_ss", fmt(ss_env)))

tab2_path <- file.path(combo_root, "mouse_table", "table2", "Table_2_total_cell.txt")
map_path  <- file.path(combo_root, "Merfish_brain_cell_type_subclass.txt")

## cluster_total：用于拿 cluster 的注释信息（layer/subclass/merge_regions/total_cell_num）
clu_path <- file.path(combo_root, "mouse_subclass_cluster_total_over3.txt")
if (!file.exists(clu_path)) clu_path <- file.path(combo_root, "mouse_subclass_cluster_total.txt")

slice_dir <- file.path(base_dir, "neocortex_new")

outdir <- file.path(combo_root, "fig_glut_gaba_overlap100")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(tab2_path))
stopifnot(file.exists(map_path))
stopifnot(file.exists(clu_path))
stopifnot(dir.exists(slice_dir))

## ===================== 1) 依赖 =====================
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(Cairo)
})

## ===================== 2) 读入数据 =====================
tab2 <- fread(tab2_path, sep="\t", header=TRUE, data.table=TRUE, check.names=FALSE)
mp   <- fread(map_path,  sep="\t", header=TRUE, data.table=TRUE, check.names=FALSE)
if (!all(c("subclass","cell_Neruon_type") %in% names(mp))) {
  stopifnot(ncol(mp) >= 2)
  setnames(mp, old=names(mp)[1:2], new=c("subclass","cell_Neruon_type"))
  if (ncol(mp) >= 3) setnames(mp, old=names(mp)[3], new="subclass_sim")
}

## ===================== 3) 通用工具 =====================
guess_col <- function(nm, patterns){
  for (p in patterns){
    hit <- grep(p, nm, ignore.case=TRUE, value=TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

norm_type <- function(x){
  x <- str_trim(as.character(x))
  x[grepl("^glu",  x, ignore.case=TRUE)]  <- "Glut"
  x[grepl("^gaba", x, ignore.case=TRUE)]  <- "Gaba"
  x
}

## 从 label 兜底解析 layer：默认 label = slide_layer_subclass_merge_regions
get_layer_from_label <- function(lbl){
  if (is.na(lbl) || lbl=="") return(NA_character_)
  ss <- strsplit(lbl, "_", fixed=TRUE)[[1]]
  if (length(ss) < 2) return(NA_character_)
  as.character(ss[2])
}

collapse_uniq <- function(x){
  ux <- unique(na.omit(as.character(x)))
  if (length(ux) == 0) return(NA_character_)
  if (length(ux) == 1) return(ux[1])
  paste(ux, collapse="|")
}

safe_int_max <- function(x){
  v <- suppressWarnings(as.integer(x))
  if (length(v)==0 || all(is.na(v))) return(NA_integer_)
  as.integer(max(v, na.rm=TRUE))
}

theme_pub <- function(base=12){
  theme_classic(base_size=base) +
    theme(
      panel.border = element_rect(color="black", fill=NA, linewidth=0.7),
      axis.text  = element_text(face="bold"),
      axis.title = element_text(face="bold"),
      plot.title = element_text(hjust=0.5, face="bold")
    )
}

## ===================== 4) 从 tab2 抽取 Glut-Gaba 且 overlap=100%（Gaba 被 Glut 完全包含） =====================
col_type1 <- guess_col(names(tab2), c("^cluster\\.1.*Neruon_type$", "cluster\\.1.*cell_Neruon_type", "cluster1.*Neruon_type"))
col_type2 <- guess_col(names(tab2), c("^cluster\\.2.*Neruon_type$", "cluster\\.2.*cell_Neruon_type", "cluster2.*Neruon_type"))
col_ov    <- guess_col(names(tab2), c("cluster\\.2\\.overlap\\.percent", "cluster2\\.overlap\\.percent", "overlap\\.percent"))

if (is.na(col_type1) || is.na(col_type2) || is.na(col_ov)) {
  stop("tab2 里找不到必要列：cluster.1_cell_Neruon_type / cluster.2_cell_Neruon_type / cluster.2.overlap.percent（或类似命名）。请检查 Table_2_total_cell.txt 列名。")
}

tab2[, type1 := norm_type(get(col_type1))]
tab2[, type2 := norm_type(get(col_type2))]
tab2[, ov_raw := as.numeric(get(col_ov))]

## 自动判断 overlap percent 是 0-1 还是 0-100
ov_max <- suppressWarnings(max(tab2$ov_raw, na.rm=TRUE))
tab2[, ov := ov_raw]
if (is.finite(ov_max) && ov_max > 1.5) tab2[, ov := ov_raw / 100]

pairs100 <- tab2[type1=="Glut" & type2=="Gaba" & is.finite(ov) & ov >= 0.999999]
if (nrow(pairs100) == 0) {
  stop("没有找到 Glut-Gaba 且 cluster.2.overlap.percent=100% 的记录（ov>=0.999999）。")
}

## cluster label & slide（尽量自动找）
col_lab1  <- guess_col(names(pairs100), c("cluster\\.1.*label", "cluster1.*label", "^label1$", "cluster\\.1\\.label"))
col_lab2  <- guess_col(names(pairs100), c("cluster\\.2.*label", "cluster2.*label", "^label2$", "cluster\\.2\\.label"))
col_slide <- guess_col(names(pairs100), c("^slide$", "chip", "slice"))

pairs100[, label1 := if (!is.na(col_lab1)) as.character(get(col_lab1)) else NA_character_]
pairs100[, label2 := if (!is.na(col_lab2)) as.character(get(col_lab2)) else NA_character_]

pairs100[, slide := if (!is.na(col_slide)) as.character(get(col_slide)) else NA_character_]
pairs100[is.na(slide) | slide=="", slide := ifelse(!is.na(label1) & label1!="", sub("_.*$", "", label1), NA_character_)]

## ===================== 5) 读 cluster_total 作为注释表（layer/subclass/merge_regions/size） =====================
clu_anno <- fread(clu_path, sep="\t", header=TRUE, data.table=TRUE, check.names=FALSE)

if (!"label" %in% names(clu_anno)) {
  req <- c("slide","layer","subclass","merge_regions")
  stopifnot(all(req %in% names(clu_anno)))
  clu_anno[, label := paste(slide, layer, subclass, merge_regions, sep="_")]
}

keep_cols <- intersect(c("label","slide","layer","subclass","merge_regions","total_cell_num"), names(clu_anno))
clu_anno <- clu_anno[, ..keep_cols]
setkey(clu_anno, label)

## ===================== 6) pairs100 -> pairs_stat：补齐 cluster1/cluster2 的 layer/subclass/merge_regions/size =====================
layer_levels <- c("1","2/3","4","5","6a","6b")

pairs_stat <- copy(pairs100)
pairs_stat[, pair_id := .I]

## cluster2 = Gaba (被包含者)
pairs_stat <- merge(
  pairs_stat,
  clu_anno[, .(
    label,
    slide2 = as.character(slide),
    layer2 = as.character(layer),
    subclass2 = as.character(if ("subclass" %in% names(clu_anno)) subclass else NA_character_),
    merge_regions2 = as.character(if ("merge_regions" %in% names(clu_anno)) merge_regions else NA_character_),
    size2 = as.integer(if ("total_cell_num" %in% names(clu_anno)) total_cell_num else NA_integer_)
  )],
  by.x="label2", by.y="label", all.x=TRUE
)

## cluster1 = Glut (包含者)
pairs_stat <- merge(
  pairs_stat,
  clu_anno[, .(
    label,
    slide1 = as.character(slide),
    layer1 = as.character(layer),
    subclass1 = as.character(if ("subclass" %in% names(clu_anno)) subclass else NA_character_),
    merge_regions1 = as.character(if ("merge_regions" %in% names(clu_anno)) merge_regions else NA_character_),
    size1 = as.integer(if ("total_cell_num" %in% names(clu_anno)) total_cell_num else NA_integer_)
  )],
  by.x="label1", by.y="label", all.x=TRUE
)

## slide：优先 tab2 的 slide；缺失则用 slide2 补
pairs_stat[is.na(slide) | slide=="", slide := slide2]

## 规范 layer token
pairs_stat[, layer2 := tolower(gsub("^\\s*layer\\s*", "", as.character(layer2), ignore.case=TRUE))]
pairs_stat[, layer1 := tolower(gsub("^\\s*layer\\s*", "", as.character(layer1), ignore.case=TRUE))]

## layer 兜底解析
pairs_stat[is.na(layer2) | layer2=="", layer2 := vapply(label2, get_layer_from_label, character(1))]
pairs_stat[is.na(layer1) | layer1=="", layer1 := vapply(label1, get_layer_from_label, character(1))]

pairs_stat[, layer2 := ifelse(layer2 %chin% layer_levels, layer2, NA_character_)]
pairs_stat[, layer1 := ifelse(layer1 %chin% layer_levels, layer1, NA_character_)]

pairs_stat[, layer2_f := factor(layer2, levels=layer_levels)]
pairs_stat[, layer1_f := factor(layer1, levels=layer_levels)]

## ===================== 7) 整体统计 + 可视化（不逐个 pair 画图） =====================

## 7.1 每层 pair 数量（按 cluster2=Gaba 所在层）
dt_pair_by_layer <- pairs_stat[!is.na(layer2_f), .N, by=layer2_f][order(layer2_f)]
dt_pair_by_layer[, pct := N/sum(N)]
fwrite(dt_pair_by_layer, file.path(outdir, "pair_count_by_layer.tsv"), sep="\t", quote=FALSE)

p1 <- ggplot(dt_pair_by_layer, aes(x=layer2_f, y=N)) +
  geom_col(width=0.72, fill="#4EA8DE", color="black", linewidth=0.3) +
  geom_text(aes(label=comma(N)), vjust=-0.35, size=4.2, fontface="bold") +
  labs(x="Layer (Gaba cluster layer)", y="Number of Glut→Gaba (ov=100%) pairs",
       title="Pairs per layer (cluster2=Gaba)") +
  scale_y_continuous(expand=expansion(mult=c(0,0.12))) +
  theme_pub(13)

ggsave(file.path(outdir, "Pairs_per_layer_bar.png"), p1, width=8.5, height=4.8, dpi=350, bg="white")
ggsave(file.path(outdir, "Pairs_per_layer_bar.pdf"), p1, width=8.5, height=4.8, device=cairo_pdf)

## 7.2 每层 unique 的 Gaba cluster 数量（label2 去重）
dt_unique_gaba <- unique(pairs_stat[!is.na(layer2_f) & !is.na(label2) & label2!="", .(label2, layer2_f, subclass2)])
dt_unique_gaba2 <- dt_unique_gaba[, .(n_unique_gaba_clusters=.N), by=layer2_f][order(layer2_f)]
fwrite(dt_unique_gaba2, file.path(outdir, "unique_gaba_cluster_by_layer.tsv"), sep="\t", quote=FALSE)

p2 <- ggplot(dt_unique_gaba2, aes(x=layer2_f, y=n_unique_gaba_clusters)) +
  geom_col(width=0.72, fill="#6FB06F", color="black", linewidth=0.3) +
  geom_text(aes(label=comma(n_unique_gaba_clusters)), vjust=-0.35, size=4.2, fontface="bold") +
  labs(x="Layer (Gaba cluster layer)", y="Unique Gaba clusters (label2)",
       title="Unique Gaba clusters per layer") +
  scale_y_continuous(expand=expansion(mult=c(0,0.12))) +
  theme_pub(13)

ggsave(file.path(outdir, "UniqueGabaClusters_per_layer_bar.png"), p2, width=8.5, height=4.8, dpi=350, bg="white")
ggsave(file.path(outdir, "UniqueGabaClusters_per_layer_bar.pdf"), p2, width=8.5, height=4.8, device=cairo_pdf)

## 7.3 Top slides 的 pair 数量（横向柱图）
dt_pair_by_slide <- pairs_stat[!is.na(slide) & slide!="", .N, by=slide][order(-N)]
topK_slide <- min(25, nrow(dt_pair_by_slide))
dt_top_slide <- dt_pair_by_slide[1:topK_slide]
dt_top_slide[, slide_f := factor(slide, levels=rev(slide))]
fwrite(dt_pair_by_slide, file.path(outdir, "pair_count_by_slide.tsv"), sep="\t", quote=FALSE)

p3 <- ggplot(dt_top_slide, aes(x=slide_f, y=N)) +
  geom_col(width=0.72, fill="#4EA8DE", color="black", linewidth=0.3) +
  geom_text(aes(label=comma(N)), hjust=-0.15, size=3.6, fontface="bold") +
  coord_flip(clip="off") +
  scale_y_continuous(expand=expansion(mult=c(0,0.15))) +
  labs(x=NULL, y="Number of pairs", title=paste0("Top ", topK_slide, " slides by #pairs")) +
  theme_pub(12) +
  theme(plot.margin = margin(6, 30, 6, 6))

ggsave(file.path(outdir, "TopSlides_pairs_barh.png"), p3, width=8.8, height=6.2, dpi=350, bg="white")
ggsave(file.path(outdir, "TopSlides_pairs_barh.pdf"), p3, width=8.8, height=6.2, device=cairo_pdf)

## 7.4 Top slides × layer 热图
dt_heat <- pairs_stat[!is.na(layer2_f) & !is.na(slide) & slide!="", .N, by=.(slide, layer2_f)]
top_slides <- dt_pair_by_slide[1:topK_slide, slide]
dt_heat <- dt_heat[slide %chin% top_slides]
dt_heat[, slide_f := factor(slide, levels=rev(top_slides))]
dt_heat[, layer_f := factor(as.character(layer2_f), levels=layer_levels)]

p4 <- ggplot(dt_heat, aes(x=layer_f, y=slide_f, fill=N)) +
  geom_tile(color="black", linewidth=0.2) +
  scale_fill_gradient(low="white", high="red", labels=comma) +
  labs(x="Layer (Gaba cluster layer)", y=NULL, fill="Pairs",
       title=paste0("Heatmap: Top ", topK_slide, " slides × layer")) +
  theme_pub(11)

ggsave(file.path(outdir, "Heatmap_TopSlides_byLayer.png"), p4, width=7.8, height=6.6, dpi=350, bg="white")
ggsave(file.path(outdir, "Heatmap_TopSlides_byLayer.pdf"), p4, width=7.8, height=6.6, device=cairo_pdf)

## 7.5 额外图：Gaba clusters（unique label2）的 subclass 组成（Top 30 subclass）
dt_gaba_subclass <- dt_unique_gaba[!is.na(subclass2) & subclass2!="", .N, by=subclass2][order(-N)]
topS <- min(30, nrow(dt_gaba_subclass))
dt_gaba_subclass_top <- dt_gaba_subclass[1:topS]
dt_gaba_subclass_top[, subclass_f := factor(subclass2, levels=rev(subclass2))]

p5 <- ggplot(dt_gaba_subclass_top, aes(x=subclass_f, y=N)) +
  geom_col(width=0.72, fill="#9AA5B1", color="black", linewidth=0.3) +
  geom_text(aes(label=comma(N)), hjust=-0.15, size=3.2, fontface="bold") +
  coord_flip(clip="off") +
  scale_y_continuous(expand=expansion(mult=c(0,0.15))) +
  labs(x=NULL, y="Unique Gaba clusters", title=paste0("Unique Gaba clusters: Top ", topS, " subclasses")) +
  theme_pub(11) +
  theme(plot.margin = margin(6, 35, 6, 6))

ggsave(file.path(outdir, "UniqueGaba_subclass_top30.png"), p5, width=10.0, height=7.0, dpi=350, bg="white")
ggsave(file.path(outdir, "UniqueGaba_subclass_top30.pdf"), p5, width=10.0, height=7.0, device=cairo_pdf)

## 7.6 输出 summary 文本
summary_txt <- file.path(outdir, "summary_glut_contains_gaba_ov100.txt")
sink(summary_txt)
cat("Total pairs (Glut contains Gaba, ov=100%): ", nrow(pairs_stat), "\n", sep="")
cat("Unique Gaba clusters (label2): ", nrow(dt_unique_gaba), "\n", sep="")
cat("Unique Glut clusters (label1): ", uniqueN(pairs_stat$label1), "\n\n", sep="")
cat("Pairs per layer (cluster2=Gaba):\n"); print(dt_pair_by_layer); cat("\n")
cat("Unique Gaba clusters per layer:\n"); print(dt_unique_gaba2); cat("\n")
cat("Top slides:\n"); print(dt_top_slide); cat("\n")
sink()

## ======================================================================
## ===================== 8) 特殊情况统计：多对一 / 一对多（100% overlap） =====================
## 你新增诉求：Top25 具体 cluster 是什么 subclass（并导出、并画图）
## ======================================================================
ps <- as.data.table(pairs_stat)
ps <- ps[!is.na(label1) & label1!="" & !is.na(label2) & label2!=""]

## --------------------- 8.1 多对一：同一个 Gaba cluster 被多个 Glut cluster 完全包含 ---------------------
dt_gaba_in_many_glut <- ps[, .(
  n_glut        = as.integer(uniqueN(label1)),
  slide2        = collapse_uniq(slide2),
  layer2        = collapse_uniq(layer2),
  subclass2     = collapse_uniq(subclass2),
  merge_regions2= collapse_uniq(merge_regions2),
  size2         = safe_int_max(size2)
), by=.(label2)]

dt_gaba_in_many_glut[, multi := (n_glut >= 2L)]
dt_gaba_in_many_glut[, layer2 := ifelse(tolower(layer2) %chin% layer_levels, tolower(layer2), layer2)]
dt_gaba_in_many_glut[, layer2_f := factor(layer2, levels=layer_levels)]

## --------------------- 8.2 一对多：同一个 Glut cluster 完全包含多个 Gaba cluster ---------------------
dt_glut_contains_many_gaba <- ps[, .(
  n_gaba        = as.integer(uniqueN(label2)),
  slide1        = collapse_uniq(slide1),
  layer1        = collapse_uniq(layer1),
  subclass1     = collapse_uniq(subclass1),
  merge_regions1= collapse_uniq(merge_regions1),
  size1         = safe_int_max(size1)
), by=.(label1)]

dt_glut_contains_many_gaba[, multi := (n_gaba >= 2L)]
dt_glut_contains_many_gaba[, layer1 := ifelse(tolower(layer1) %chin% layer_levels, tolower(layer1), layer1)]
dt_glut_contains_many_gaba[, layer1_f := factor(layer1, levels=layer_levels)]

## --------------------- 8.3 总览统计输出 ---------------------
pct_gaba_multi <- if (nrow(dt_gaba_in_many_glut)==0) NA_real_ else as.numeric(sum(dt_gaba_in_many_glut$multi)/nrow(dt_gaba_in_many_glut))
pct_glut_multi <- if (nrow(dt_glut_contains_many_gaba)==0) NA_real_ else as.numeric(sum(dt_glut_contains_many_gaba$multi)/nrow(dt_glut_contains_many_gaba))

sum_multi <- data.table(
  metric = c(
    "Total pairs (ov100, Glut->Gaba)",
    "Unique Gaba clusters (label2)",
    "Unique Glut clusters (label1)",
    "Gaba clusters with >=2 Glut containers",
    "Glut clusters containing >=2 Gaba clusters"
  ),
  value = c(
    nrow(ps),
    nrow(dt_gaba_in_many_glut),
    nrow(dt_glut_contains_many_gaba),
    sum(dt_gaba_in_many_glut$multi),
    sum(dt_glut_contains_many_gaba$multi)
  ),
  pct = c(NA_real_, NA_real_, NA_real_, pct_gaba_multi, pct_glut_multi)
)
fwrite(sum_multi, file.path(outdir, "multi_summary.tsv"), sep="\t", quote=FALSE)

## --------------------- 8.4 按 layer 汇总（类型统一，避免 data.table 报错） ---------------------
dt_gaba_by_layer <- dt_gaba_in_many_glut[!is.na(layer2) & layer2!="", .(
  n_gaba_clusters = as.integer(.N),
  n_multi         = as.integer(sum(multi)),
  pct_multi       = as.numeric(sum(multi)/.N),
  mean_n_glut     = as.numeric(mean(as.numeric(n_glut))),
  median_n_glut   = as.numeric(median(as.numeric(n_glut))),
  max_n_glut      = as.numeric(max(as.numeric(n_glut)))
), by=.(layer2)][]
dt_gaba_by_layer[, layer2_f := factor(layer2, levels=layer_levels)]
setorder(dt_gaba_by_layer, layer2_f)
fwrite(dt_gaba_by_layer[, .(layer2, n_gaba_clusters, n_multi, pct_multi, mean_n_glut, median_n_glut, max_n_glut)],
       file.path(outdir, "multi_gaba_by_layer.tsv"), sep="\t", quote=FALSE)

dt_glut_by_layer <- dt_glut_contains_many_gaba[!is.na(layer1) & layer1!="", .(
  n_glut_clusters = as.integer(.N),
  n_multi         = as.integer(sum(multi)),
  pct_multi       = as.numeric(sum(multi)/.N),
  mean_n_gaba     = as.numeric(mean(as.numeric(n_gaba))),
  median_n_gaba   = as.numeric(median(as.numeric(n_gaba))),
  max_n_gaba      = as.numeric(max(as.numeric(n_gaba)))
), by=.(layer1)][]
dt_glut_by_layer[, layer1_f := factor(layer1, levels=layer_levels)]
setorder(dt_glut_by_layer, layer1_f)
fwrite(dt_glut_by_layer[, .(layer1, n_glut_clusters, n_multi, pct_multi, mean_n_gaba, median_n_gaba, max_n_gaba)],
       file.path(outdir, "multi_glut_by_layer.tsv"), sep="\t", quote=FALSE)

## --------------------- 8.5 分布图：n_glut per Gaba / n_gaba per Glut（并补充 CDF） ---------------------
capN <- 10L
dt_gaba_in_many_glut[, ncat := ifelse(n_glut >= capN, paste0(">=", capN), as.character(n_glut))]
dt_glut_contains_many_gaba[, ncat := ifelse(n_gaba >= capN, paste0(">=", capN), as.character(n_gaba))]
levs <- c(as.character(1:(capN-1)), paste0(">=", capN))
dt_gaba_in_many_glut[, ncat := factor(ncat, levels=levs)]
dt_glut_contains_many_gaba[, ncat := factor(ncat, levels=levs)]

dist_gaba <- dt_gaba_in_many_glut[, .N, by=ncat][order(ncat)]
dist_glut <- dt_glut_contains_many_gaba[, .N, by=ncat][order(ncat)]
fwrite(dist_gaba, file.path(outdir, "dist_nGlut_perGaba.tsv"), sep="\t", quote=FALSE)
fwrite(dist_glut, file.path(outdir, "dist_nGaba_perGlut.tsv"), sep="\t", quote=FALSE)

p_dist_gaba <- ggplot(dist_gaba, aes(x=ncat, y=N)) +
  geom_col(width=0.72, fill="#4EA8DE", color="black", linewidth=0.3) +
  geom_text(aes(label=comma(N)), vjust=-0.35, size=3.8, fontface="bold") +
  scale_y_continuous(expand=expansion(mult=c(0,0.15))) +
  labs(x="#Glut containers per Gaba (ov=100%)", y="#Gaba clusters",
       title="Distribution: #Glut containers per Gaba cluster") +
  theme_pub(12)

p_dist_glut <- ggplot(dist_glut, aes(x=ncat, y=N)) +
  geom_col(width=0.72, fill="#6FB06F", color="black", linewidth=0.3) +
  geom_text(aes(label=comma(N)), vjust=-0.35, size=3.8, fontface="bold") +
  scale_y_continuous(expand=expansion(mult=c(0,0.15))) +
  labs(x="#Gaba contained per Glut (ov=100%)", y="#Glut clusters",
       title="Distribution: #Gaba contained per Glut cluster") +
  theme_pub(12)

ggsave(file.path(outdir, "Dist_nGlut_perGaba.png"), p_dist_gaba, width=9.2, height=4.8, dpi=350, bg="white")
ggsave(file.path(outdir, "Dist_nGlut_perGaba.pdf"), p_dist_gaba, width=9.2, height=4.8, device=cairo_pdf)
ggsave(file.path(outdir, "Dist_nGaba_perGlut.png"), p_dist_glut, width=9.2, height=4.8, dpi=350, bg="white")
ggsave(file.path(outdir, "Dist_nGaba_perGlut.pdf"), p_dist_glut, width=9.2, height=4.8, device=cairo_pdf)

## CDF（更适合论文一句话描述“80% 的 cluster 在 n<=2”这种）
cdf_gaba <- dt_gaba_in_many_glut[!is.na(n_glut), .N, by=n_glut][order(n_glut)]
cdf_gaba[, cdf := cumsum(N)/sum(N)]
cdf_glut <- dt_glut_contains_many_gaba[!is.na(n_gaba), .N, by=n_gaba][order(n_gaba)]
cdf_glut[, cdf := cumsum(N)/sum(N)]

p_cdf_gaba <- ggplot(cdf_gaba, aes(x=n_glut, y=cdf)) +
  geom_step(linewidth=1.0) +
  geom_point(size=2.2) +
  scale_y_continuous(labels=percent_format(accuracy=1), limits=c(0,1)) +
  labs(x="n_glut (containers per Gaba)", y="CDF",
       title="CDF of #Glut containers per Gaba cluster") +
  theme_pub(12)

p_cdf_glut <- ggplot(cdf_glut, aes(x=n_gaba, y=cdf)) +
  geom_step(linewidth=1.0) +
  geom_point(size=2.2) +
  scale_y_continuous(labels=percent_format(accuracy=1), limits=c(0,1)) +
  labs(x="n_gaba (contained per Glut)", y="CDF",
       title="CDF of #Gaba contained per Glut cluster") +
  theme_pub(12)

ggsave(file.path(outdir, "CDF_nGlut_perGaba.png"), p_cdf_gaba, width=7.6, height=4.8, dpi=350, bg="white")
ggsave(file.path(outdir, "CDF_nGlut_perGaba.pdf"), p_cdf_gaba, width=7.6, height=4.8, device=cairo_pdf)
ggsave(file.path(outdir, "CDF_nGaba_perGlut.png"), p_cdf_glut, width=7.6, height=4.8, dpi=350, bg="white")
ggsave(file.path(outdir, "CDF_nGaba_perGlut.pdf"), p_cdf_glut, width=7.6, height=4.8, device=cairo_pdf)

## --------------------- 8.6 新增图：按 layer 的 n 分布（箱线图/散点） ---------------------
p_box_gaba <- ggplot(dt_gaba_in_many_glut[!is.na(layer2_f)], aes(x=layer2_f, y=n_glut)) +
  geom_boxplot(width=0.65, outlier.size=0.8) +
  geom_jitter(width=0.18, height=0, alpha=0.35, size=1.6) +
  labs(x="Layer (Gaba cluster layer)", y="n_glut (how many Glut contain it)",
       title="By layer: n_glut distribution (Gaba clusters, ov=100%)") +
  theme_pub(12)

p_box_glut <- ggplot(dt_glut_contains_many_gaba[!is.na(layer1_f)], aes(x=layer1_f, y=n_gaba)) +
  geom_boxplot(width=0.65, outlier.size=0.8) +
  geom_jitter(width=0.18, height=0, alpha=0.35, size=1.6) +
  labs(x="Layer (Glut cluster layer)", y="n_gaba (how many Gaba it contains)",
       title="By layer: n_gaba distribution (Glut clusters, ov=100%)") +
  theme_pub(12)

ggsave(file.path(outdir, "Box_nGlut_byLayer_Gaba.png"), p_box_gaba, width=8.6, height=4.8, dpi=350, bg="white")
ggsave(file.path(outdir, "Box_nGlut_byLayer_Gaba.pdf"), p_box_gaba, width=8.6, height=4.8, device=cairo_pdf)
ggsave(file.path(outdir, "Box_nGaba_byLayer_Glut.png"), p_box_glut, width=8.6, height=4.8, dpi=350, bg="white")
ggsave(file.path(outdir, "Box_nGaba_byLayer_Glut.pdf"), p_box_glut, width=8.6, height=4.8, device=cairo_pdf)

## --------------------- 8.7 关键新增：Top25 cluster + subclass（导出 + 更完备可视化） ---------------------
topK <- 25L

## Top25：Gaba clusters（被多个 Glut 完全包含）
top_gaba <- dt_gaba_in_many_glut[multi==TRUE][order(-n_glut, -size2)][1:min(topK, .N)]
if (nrow(top_gaba) > 0) {
  ## 展示名：尽量用 slide_layer_merge_regions；并把 subclass 放进去
  top_gaba[, show_core := ifelse(!is.na(merge_regions2) & merge_regions2!="",
                                 paste0(slide2, "_", layer2, "_", merge_regions2),
                                 label2)]
  top_gaba[, show := ifelse(!is.na(subclass2) & subclass2!="",
                            paste0(show_core, "\n", subclass2),
                            show_core)]
  top_gaba[, show_f := factor(show, levels=rev(show))]
  
  ## 导出（你要看的“top25具体cluster是什么subclass类型”就在这张表）
  fwrite(top_gaba[, .(label2, slide2, layer2, merge_regions2, subclass2, size2, n_glut)],
         file.path(outdir, "Top25_Gaba_multiContained_withSubclass.tsv"),
         sep="\t", quote=FALSE)
  
  p_top_gaba <- ggplot(top_gaba, aes(x=show_f, y=n_glut)) +
    geom_col(width=0.72, fill="#1f77b4", color="black", linewidth=0.3) +
    geom_text(aes(label=paste0("n=", n_glut)), hjust=-0.15, size=3.6, fontface="bold") +
    coord_flip(clip="off") +
    scale_y_continuous(expand=expansion(mult=c(0,0.18))) +
    labs(x=NULL,
         y="#Glut clusters (ov100) containing this Gaba cluster",
         title=paste0("Top ", min(topK, nrow(top_gaba)), " Gaba clusters (multi-contained) + subclass")) +
    theme_pub(11) +
    theme(plot.margin = margin(6, 40, 6, 6))
  
  ggsave(file.path(outdir, "Top25_Gaba_multiContained.png"), p_top_gaba, width=11.2, height=7.2, dpi=350, bg="white")
  ggsave(file.path(outdir, "Top25_Gaba_multiContained.pdf"), p_top_gaba, width=11.2, height=7.2, device=cairo_pdf)
}

## Top25：Glut clusters（包含多个 Gaba）
top_glut <- dt_glut_contains_many_gaba[multi==TRUE][order(-n_gaba, -size1)][1:min(topK, .N)]
if (nrow(top_glut) > 0) {
  top_glut[, show_core := ifelse(!is.na(merge_regions1) & merge_regions1!="",
                                 paste0(slide1, "_", layer1, "_", merge_regions1),
                                 label1)]
  top_glut[, show := ifelse(!is.na(subclass1) & subclass1!="",
                            paste0(show_core, "\n", subclass1),
                            show_core)]
  top_glut[, show_f := factor(show, levels=rev(show))]
  
  fwrite(top_glut[, .(label1, slide1, layer1, merge_regions1, subclass1, size1, n_gaba)],
         file.path(outdir, "Top25_Glut_containsMany_withSubclass.tsv"),
         sep="\t", quote=FALSE)
  
  p_top_glut <- ggplot(top_glut, aes(x=show_f, y=n_gaba)) +
    geom_col(width=0.72, fill="#d62728", color="black", linewidth=0.3) +
    geom_text(aes(label=paste0("n=", n_gaba)), hjust=-0.15, size=3.6, fontface="bold") +
    coord_flip(clip="off") +
    scale_y_continuous(expand=expansion(mult=c(0,0.18))) +
    labs(x=NULL,
         y="#Gaba clusters (ov100) contained in this Glut cluster",
         title=paste0("Top ", min(topK, nrow(top_glut)), " Glut clusters (contain many) + subclass")) +
    theme_pub(11) +
    theme(plot.margin = margin(6, 40, 6, 6))
  
  ggsave(file.path(outdir, "Top25_Glut_containsMany.png"), p_top_glut, width=11.2, height=7.2, dpi=350, bg="white")
  ggsave(file.path(outdir, "Top25_Glut_containsMany.pdf"), p_top_glut, width=11.2, height=7.2, device=cairo_pdf)
}

## --------------------- 8.8 新增图：multi clusters 的 subclass 组成（Top30） ---------------------
gaba_multi_sub <- dt_gaba_in_many_glut[multi==TRUE & !is.na(subclass2) & subclass2!="", .N, by=subclass2][order(-N)]
glut_multi_sub <- dt_glut_contains_many_gaba[multi==TRUE & !is.na(subclass1) & subclass1!="", .N, by=subclass1][order(-N)]

topS2 <- min(30, nrow(gaba_multi_sub))
if (topS2 > 0) {
  gaba_multi_sub2 <- gaba_multi_sub[1:topS2]
  gaba_multi_sub2[, subclass_f := factor(subclass2, levels=rev(subclass2))]
  p_sub_gaba <- ggplot(gaba_multi_sub2, aes(x=subclass_f, y=N)) +
    geom_col(width=0.72, fill="#1f77b4", color="black", linewidth=0.3) +
    geom_text(aes(label=comma(N)), hjust=-0.15, size=3.2, fontface="bold") +
    coord_flip(clip="off") +
    scale_y_continuous(expand=expansion(mult=c(0,0.15))) +
    labs(x=NULL, y="#Gaba multi-contained clusters",
         title=paste0("Multi-contained Gaba clusters: Top ", topS2, " subclasses")) +
    theme_pub(11) +
    theme(plot.margin = margin(6, 40, 6, 6))
  
  ggsave(file.path(outdir, "Subclass_multiGaba_top30.png"), p_sub_gaba, width=10.8, height=7.0, dpi=350, bg="white")
  ggsave(file.path(outdir, "Subclass_multiGaba_top30.pdf"), p_sub_gaba, width=10.8, height=7.0, device=cairo_pdf)
  fwrite(gaba_multi_sub, file.path(outdir, "Subclass_multiGaba_all.tsv"), sep="\t", quote=FALSE)
}

topS3 <- min(30, nrow(glut_multi_sub))
if (topS3 > 0) {
  glut_multi_sub2 <- glut_multi_sub[1:topS3]
  glut_multi_sub2[, subclass_f := factor(subclass1, levels=rev(subclass1))]
  p_sub_glut <- ggplot(glut_multi_sub2, aes(x=subclass_f, y=N)) +
    geom_col(width=0.72, fill="#d62728", color="black", linewidth=0.3) +
    geom_text(aes(label=comma(N)), hjust=-0.15, size=3.2, fontface="bold") +
    coord_flip(clip="off") +
    scale_y_continuous(expand=expansion(mult=c(0,0.15))) +
    labs(x=NULL, y="#Glut clusters (contain many Gaba)",
         title=paste0("Glut clusters containing many Gaba: Top ", topS3, " subclasses")) +
    theme_pub(11) +
    theme(plot.margin = margin(6, 40, 6, 6))
  
  ggsave(file.path(outdir, "Subclass_multiGlut_top30.png"), p_sub_glut, width=10.8, height=7.0, dpi=350, bg="white")
  ggsave(file.path(outdir, "Subclass_multiGlut_top30.pdf"), p_sub_glut, width=10.8, height=7.0, device=cairo_pdf)
  fwrite(glut_multi_sub, file.path(outdir, "Subclass_multiGlut_all.tsv"), sep="\t", quote=FALSE)
}

## --------------------- 8.9 新增图：multi clusters 的 layer × subclass 热图（只取 top20 subclass，避免太大） ---------------------
topHeat <- 20L
if (nrow(gaba_multi_sub) > 0) {
  top_subs_gaba <- gaba_multi_sub[1:min(topHeat, .N), subclass2]
  heat_gaba <- dt_gaba_in_many_glut[multi==TRUE & subclass2 %chin% top_subs_gaba & !is.na(layer2_f),
                                    .N, by=.(layer2_f, subclass2)]
  heat_gaba[, subclass_f := factor(subclass2, levels=rev(top_subs_gaba))]
  p_hg <- ggplot(heat_gaba, aes(x=layer2_f, y=subclass_f, fill=N)) +
    geom_tile(color="black", linewidth=0.2) +
    scale_fill_gradient(low="white", high="red", labels=comma) +
    labs(x="Layer", y="Subclass (top)", fill="#clusters",
         title=paste0("Multi-contained Gaba: layer × subclass (Top ", min(topHeat, length(top_subs_gaba)), ")")) +
    theme_pub(11)
  ggsave(file.path(outdir, "Heatmap_multiGaba_layer_subclass.png"), p_hg, width=8.6, height=6.8, dpi=350, bg="white")
  ggsave(file.path(outdir, "Heatmap_multiGaba_layer_subclass.pdf"), p_hg, width=8.6, height=6.8, device=cairo_pdf)
}

if (nrow(glut_multi_sub) > 0) {
  top_subs_glut <- glut_multi_sub[1:min(topHeat, .N), subclass1]
  heat_glut <- dt_glut_contains_many_gaba[multi==TRUE & subclass1 %chin% top_subs_glut & !is.na(layer1_f),
                                          .N, by=.(layer1_f, subclass1)]
  heat_glut[, subclass_f := factor(subclass1, levels=rev(top_subs_glut))]
  p_hl <- ggplot(heat_glut, aes(x=layer1_f, y=subclass_f, fill=N)) +
    geom_tile(color="black", linewidth=0.2) +
    scale_fill_gradient(low="white", high="red", labels=comma) +
    labs(x="Layer", y="Subclass (top)", fill="#clusters",
         title=paste0("Glut contains-many: layer × subclass (Top ", min(topHeat, length(top_subs_glut)), ")")) +
    theme_pub(11)
  ggsave(file.path(outdir, "Heatmap_multiGlut_layer_subclass.png"), p_hl, width=8.6, height=6.8, dpi=350, bg="white")
  ggsave(file.path(outdir, "Heatmap_multiGlut_layer_subclass.pdf"), p_hl, width=8.6, height=6.8, device=cairo_pdf)
}

## --------------------- 8.10 导出详细表（全量） ---------------------
fwrite(dt_gaba_in_many_glut[order(-n_glut, -size2)],
       file.path(outdir, "detail_gaba_in_many_glut.tsv"),
       sep="\t", quote=FALSE)

fwrite(dt_glut_contains_many_gaba[order(-n_gaba, -size1)],
       file.path(outdir, "detail_glut_contains_many_gaba.tsv"),
       sep="\t", quote=FALSE)


## ===================== 9) 严格1对1完全包含pairs筛选（ov=100%） =====================
stopifnot(exists("pairs_stat"), exists("outdir"))

library(data.table)
library(ggplot2)
library(scales)
library(Cairo)

ps <- as.data.table(pairs_stat)
ps <- ps[!is.na(label1) & label1 != "" & !is.na(label2) & label2 != ""]

## 若前面没生成这两个表，就在这里现算（保证可独立运行）
if (!exists("dt_gaba_in_many_glut")) {
  dt_gaba_in_many_glut <- ps[, .(
    n_glut = as.integer(uniqueN(label1))
  ), by = .(label2)]
}
if (!exists("dt_glut_contains_many_gaba")) {
  dt_glut_contains_many_gaba <- ps[, .(
    n_gaba = as.integer(uniqueN(label2))
  ), by = .(label1)]
}

## 1) 找到“只被1个Glut完全包含”的Gaba clusters
gaba_1only <- dt_gaba_in_many_glut[n_glut == 1L, unique(label2)]

## 2) 找到“只完全包含1个Gaba”的Glut clusters
glut_1only <- dt_glut_contains_many_gaba[n_gaba == 1L, unique(label1)]

## 3) 同时满足两边唯一性的pairs（严格1对1）
pairs_1to1 <- ps[
  label1 %chin% glut_1only & label2 %chin% gaba_1only,
  .(
    label1, label2,
    slide  = if ("slide"  %in% names(ps)) as.character(slide)  else NA_character_,
    slide1 = if ("slide1" %in% names(ps)) as.character(slide1) else NA_character_,
    slide2 = if ("slide2" %in% names(ps)) as.character(slide2) else NA_character_,
    layer1 = if ("layer1" %in% names(ps)) as.character(layer1) else NA_character_,
    layer2 = if ("layer2" %in% names(ps)) as.character(layer2) else NA_character_,
    layer1_f = if ("layer1_f" %in% names(ps)) as.character(layer1_f) else NA_character_,
    layer2_f = if ("layer2_f" %in% names(ps)) as.character(layer2_f) else NA_character_,
    subclass1 = if ("subclass1" %in% names(ps)) as.character(subclass1) else NA_character_,
    subclass2 = if ("subclass2" %in% names(ps)) as.character(subclass2) else NA_character_,
    merge_regions1 = if ("merge_regions1" %in% names(ps)) as.character(merge_regions1) else NA_character_,
    merge_regions2 = if ("merge_regions2" %in% names(ps)) as.character(merge_regions2) else NA_character_,
    size1 = if ("size1" %in% names(ps)) as.integer(size1) else NA_integer_,
    size2 = if ("size2" %in% names(ps)) as.integer(size2) else NA_integer_,
    ov = if ("ov" %in% names(ps)) as.numeric(ov) else 1
  )
]

## 去重（理论上严格1对1应当每对只有1行，但保险起见）
pairs_1to1 <- unique(pairs_1to1, by = c("label1","label2"))

## 4) 输出汇总
cat("严格1对1pairs数量：", nrow(pairs_1to1), "\n", sep="")
cat("占全部ov100pairs比例：", round(100 * nrow(pairs_1to1) / nrow(ps), 2), "%\n", sep="")

## 5) 写出结果表
fwrite(pairs_1to1, file.path(outdir, "pairs_ov100_strict_1to1.tsv"),
       sep = "\t", quote = FALSE)

## ===================== 9.1 统计图：按layer2(=Gaba所在层)的1对1pairs数量 =====================
if (!exists("layer_levels")) layer_levels <- c("1","2/3","4","5","6a","6b")

dt_1to1_by_layer <- pairs_1to1[
  !is.na(layer2_f) & layer2_f != "",
  .N, by = .(layer2_f)
]
dt_1to1_by_layer[, layer2_f := factor(layer2_f, levels = layer_levels)]
setorder(dt_1to1_by_layer, layer2_f)

## 如果你前面定义过layer_cols就直接用；否则给一个兜底
if (!exists("layer_cols")) {
  layer_cols <- c("1"="#E76F51","2/3"="#4EA8DE","4"="#6FB06F","5"="#3A4E8C","6a"="#F2A279","6b"="#9AA5B1")
}

p_1to1_layer <- ggplot(dt_1to1_by_layer, aes(x=layer2_f, y=N, fill=as.character(layer2_f))) +
  geom_col(width=0.72, color="black", linewidth=0.3) +
  geom_text(aes(label=comma(N)), vjust=-0.35, size=4.2, fontface="bold") +
  scale_fill_manual(values=layer_cols, guide="none") +
  labs(x="Layer(Gaba cluster layer)", y="Number ofstrict1-to-1pairs",
       title="Strict1-to-1Glut→Gaba(ov=100%)pairsperlayer") +
  theme_classic(base_size=13) +
  theme(
    plot.title = element_text(hjust=0.5, face="bold"),
    axis.text  = element_text(face="bold"),
    axis.title = element_text(face="bold"),
    panel.border = element_rect(color="black", fill=NA, linewidth=0.7)
  )

ggsave(file.path(outdir, "Pairs_1to1_per_layer_bar.png"), p_1to1_layer,
       width=8.5, height=4.8, dpi=350, bg="white")
ggsave(file.path(outdir, "Pairs_1to1_per_layer_bar.pdf"), p_1to1_layer,
       width=8.5, height=4.8, device=cairo_pdf)

## ===================== 9.2 统计图：Topslides(1对1pairs) =====================
dt_1to1_by_slide <- pairs_1to1[
  !is.na(slide) & slide != "",
  .N, by = slide
][order(-N)]

topK <- min(25, nrow(dt_1to1_by_slide))
dt_1to1_top_slide <- dt_1to1_by_slide[1:topK]
dt_1to1_top_slide[, slide_f := factor(slide, levels=rev(slide))]

p_1to1_slide <- ggplot(dt_1to1_top_slide, aes(x=slide_f, y=N)) +
  geom_col(width=0.72, fill="#4EA8DE", color="black", linewidth=0.3) +
  geom_text(aes(label=comma(N)), hjust=-0.15, size=3.6, fontface="bold") +
  coord_flip(clip="off") +
  scale_y_continuous(expand=expansion(mult=c(0,0.12))) +
  labs(x=NULL, y="Number ofstrict1-to-1pairs",
       title=paste0("Top", topK, "slidesbystrict1-to-1pairs")) +
  theme_classic(base_size=12) +
  theme(
    plot.title = element_text(hjust=0.5, face="bold"),
    axis.text  = element_text(face="bold"),
    axis.title = element_text(face="bold"),
    panel.border = element_rect(color="black", fill=NA, linewidth=0.7),
    plot.margin = margin(6, 30, 6, 6)
  )

ggsave(file.path(outdir, "TopSlides_1to1_pairs_barh.png"), p_1to1_slide,
       width=8.8, height=6.2, dpi=350, bg="white")
ggsave(file.path(outdir, "TopSlides_1to1_pairs_barh.pdf"), p_1to1_slide,
       width=8.8, height=6.2, device=cairo_pdf)

## ===================== 9.3 可选：size1vs size2散点(若size列存在) =====================
if ("size1" %in% names(pairs_1to1) && "size2" %in% names(pairs_1to1) &&
    any(!is.na(pairs_1to1$size1)) && any(!is.na(pairs_1to1$size2))) {
  
  p_1to1_size <- ggplot(pairs_1to1[!is.na(size1) & !is.na(size2)], aes(x=size1, y=size2)) +
    geom_point(alpha=0.35, size=1.6) +
    scale_x_continuous(trans="log10", labels=comma) +
    scale_y_continuous(trans="log10", labels=comma) +
    labs(x="Glutclustersize(total_cell_num,log10)",
         y="Gabaclustersize(total_cell_num,log10)",
         title="Strict1-to-1pairs:sizerelationship") +
    theme_classic(base_size=12) +
    theme(
      plot.title = element_text(hjust=0.5, face="bold"),
      axis.text  = element_text(face="bold"),
      axis.title = element_text(face="bold"),
      panel.border = element_rect(color="black", fill=NA, linewidth=0.7)
    )
  
  ggsave(file.path(outdir, "Scatter_1to1_size1_vs_size2.png"), p_1to1_size,
         width=6.2, height=5.6, dpi=350, bg="white")
  ggsave(file.path(outdir, "Scatter_1to1_size1_vs_size2.pdf"), p_1to1_size,
         width=6.2, height=5.6, device=cairo_pdf)
}

## ===================== 9.4 输出一个summary文本 =====================
summary_txt2 <- file.path(outdir, "summary_strict_1to1.txt")
sink(summary_txt2)
cat("Strict1-to-1pairs(ov=100%): ", nrow(pairs_1to1), "\n", sep="")
cat("Allov100pairs: ", nrow(ps), "\n", sep="")
cat("Percent: ", round(100*nrow(pairs_1to1)/nrow(ps), 3), "%\n\n", sep="")
cat("Pairsperlayer2(Gaba layer):\n")
print(dt_1to1_by_layer)
cat("\nTopslides:\n")
print(dt_1to1_top_slide)
sink()


##
## =====================================================================
## 要求1：Scatter_1to1_size1_vs_size2.png添加相关性线+标注r与p
## 说明：用log10(size1/size2)空间做散点与线性拟合；相关性同样在log10空间算
## =====================================================================
stopifnot(exists("pairs_1to1"), exists("outdir"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(Cairo)
})

fmt_p <- function(p){
  if (is.na(p)) return("NA")
  if (p < 2.2e-16) return("<2.2e-16")
  formatC(p, format="e", digits=2)
}

df_sc <- as.data.table(pairs_1to1)
df_sc <- df_sc[!is.na(size1) & !is.na(size2) & size1 > 0 & size2 > 0]
df_sc[, lx := log10(as.numeric(size1))]
df_sc[, ly := log10(as.numeric(size2))]

## 相关性：默认Pearson(线性对应lm直线)；如果你想更稳健可改method="spearman"
ct <- cor.test(df_sc$lx, df_sc$ly, method="pearson")
r_val <- unname(ct$estimate)
p_val <- ct$p.value
lab_txt <- paste0("Pearson r=", sprintf("%.3f", r_val), "\n", "p=", fmt_p(p_val))

## 注释位置：左上角(按分位数自适应)
x_anno <- quantile(df_sc$lx, 0.02, na.rm=TRUE)
y_anno <- quantile(df_sc$ly, 0.98, na.rm=TRUE)

p_1to1_size_cor <- ggplot(df_sc, aes(x=lx, y=ly)) +
  geom_point(alpha=0.35, size=1.6) +
  geom_smooth(method="lm", se=FALSE, linewidth=1.0) +
  annotate("text", x=x_anno, y=y_anno, label=lab_txt,
           hjust=0, vjust=1, size=4.2, fontface="bold") +
  labs(
    x="log10(Glutclustersize(total_cell_num))",
    y="log10(Gabaclustersize(total_cell_num))",
    title="Strict1-to-1pairs:sizerelationship(withcorrelation)"
  ) +
  theme_classic(base_size=12) +
  theme(
    plot.title = element_text(hjust=0.5, face="bold"),
    axis.text  = element_text(face="bold"),
    axis.title = element_text(face="bold"),
    panel.border = element_rect(color="black", fill=NA, linewidth=0.7)
  )

ggsave(file.path(outdir, "Scatter_1to1_size1_vs_size2.png"), p_1to1_size_cor,
       width=6.2, height=5.6, dpi=350, bg="white")
ggsave(file.path(outdir, "Scatter_1to1_size1_vs_size2.pdf"), p_1to1_size_cor,
       width=6.2, height=5.6, device=cairo_pdf)



## =====================================================================
## 要求2：保留Pairs_1to1_per_layer_bar.png不动，新增一张综合对比图
## 目标：同一张图里给出每层总pairs(Pairs_per_layer_bar) + 1to1数量 + 占比
## 输出：Pairs_1to1_ratio_vs_total_combo.png/pdf
## =====================================================================
stopifnot(exists("dt_pair_by_layer"), exists("dt_1to1_by_layer"))

dt_total <- copy(dt_pair_by_layer)            # 来自7.1
setnames(dt_total, "N", "total_pairs")
dt_one <- copy(dt_1to1_by_layer)              # 来自9.1
setnames(dt_one, "N", "pairs_1to1")

## 统一key并补0
dt_cmp <- merge(dt_total[, .(layer2_f, total_pairs)],
                dt_one[,   .(layer2_f, pairs_1to1)],
                by="layer2_f", all.x=TRUE)
dt_cmp[is.na(pairs_1to1), pairs_1to1 := 0L]
dt_cmp[, pct_1to1 := fifelse(total_pairs > 0, pairs_1to1 / total_pairs, NA_real_)]

## 做长表用于双柱(总数vs1to1)
dt_long <- melt(
  dt_cmp,
  id.vars = "layer2_f",
  measure.vars = c("total_pairs", "pairs_1to1"),
  variable.name = "kind",
  value.name = "count"
)
dt_long[, kind := factor(kind, levels=c("total_pairs","pairs_1to1"),
                         labels=c("Totalpairs(ov100)","Strict1-to-1pairs"))]

## 用同一主轴画柱+线：比例线映射到主轴(用sec_axis还原成百分比)
scale_factor <- max(dt_cmp$total_pairs, na.rm=TRUE)
if (!is.finite(scale_factor) || scale_factor <= 0) scale_factor <- 1

p_layer_combo <- ggplot() +
  geom_col(
    data = dt_long,
    aes(x=layer2_f, y=count, fill=kind),
    width=0.65, position=position_dodge(width=0.75),
    color="black", linewidth=0.3
  ) +
  geom_text(
    data = dt_long,
    aes(x=layer2_f, y=count, label=comma(count), group=kind),
    position=position_dodge(width=0.75),
    vjust=-0.35, size=3.6, fontface="bold"
  ) +
  geom_line(
    data = dt_cmp,
    aes(x=layer2_f, y=pct_1to1 * scale_factor, group=1),
    linewidth=1.0
  ) +
  geom_point(
    data = dt_cmp,
    aes(x=layer2_f, y=pct_1to1 * scale_factor),
    size=2.4
  ) +
  geom_text(
    data = dt_cmp,
    aes(x=layer2_f, y=pct_1to1 * scale_factor,
        label=percent(pct_1to1, accuracy=0.1)),
    vjust=-0.8, size=3.4, fontface="bold"
  ) +
  scale_y_continuous(
    labels=comma,
    expand=expansion(mult=c(0,0.15)),
    sec.axis = sec_axis(~ . / scale_factor, labels=percent_format(accuracy=1), name="Strict1-to-1ratio")
  ) +
  labs(
    x="Layer(Gabaclusterlayer)",
    y="Numberofpairs",
    fill=NULL,
    title="Totalov100pairs vs Strict1-to-1pairs per layer(withratio)"
  ) +
  theme_classic(base_size=12) +
  theme(
    plot.title = element_text(hjust=0.5, face="bold"),
    axis.text  = element_text(face="bold"),
    axis.title = element_text(face="bold"),
    panel.border = element_rect(color="black", fill=NA, linewidth=0.7),
    legend.position = "top"
  )

ggsave(file.path(outdir, "Pairs_1to1_ratio_vs_total_combo.png"), p_layer_combo,
       width=9.2, height=5.2, dpi=350, bg="white")
ggsave(file.path(outdir, "Pairs_1to1_ratio_vs_total_combo.pdf"), p_layer_combo,
       width=9.2, height=5.2, device=cairo_pdf)

## 同时导出综合表
fwrite(dt_cmp[, .(layer=as.character(layer2_f), total_pairs, pairs_1to1, pct_1to1)],
       file.path(outdir, "pairs_total_vs_1to1_by_layer.tsv"),
       sep="\t", quote=FALSE)



## =====================================================================
## 要求3(改进版)：按C57BL6J后数字分4只小鼠；分别输出4张热图；颜色对比更显著
## 输出：
##   Heatmap_Mouse1_1to1_bySlide.png/pdf
##   Heatmap_Mouse2_1to1_bySlide.png/pdf
##   Heatmap_Mouse3_1to1_bySlide.png/pdf
##   Heatmap_Mouse4_1to1_bySlide.png/pdf
##   mouse_sum_1to1.tsv / allslides_1to1_count_byMouse.tsv
## =====================================================================
stopifnot(exists("clu_anno"), exists("dt_1to1_by_slide"))

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(Cairo)
})

## 1)构建“所有切片”列表(包含0次的切片)
all_slides <- sort(unique(na.omit(as.character(clu_anno$slide))))
if (length(all_slides) == 0 && exists("pairs_stat")) {
  all_slides <- sort(unique(na.omit(as.character(pairs_stat$slide))))
}
stopifnot(length(all_slides) > 0)

dt_all <- data.table(slide = all_slides)

dt_slide_cnt <- copy(dt_1to1_by_slide)
setnames(dt_slide_cnt, "N", "n_1to1")

dt_all <- merge(dt_all, dt_slide_cnt, by="slide", all.x=TRUE)
dt_all[is.na(n_1to1), n_1to1 := 0L]

## 2)解析mouse原始编号：取C57BL6J后面的第一段数字
get_mouse_raw <- function(x){
  m <- stringr::str_match(x, "C57BL6J\\D*([0-9]+)")
  suppressWarnings(as.integer(m[,2]))
}
dt_all[, mouse_raw := get_mouse_raw(slide)]

## 3)把mouse_raw映射成4组(Mouse1-4)，保证一定得到1:4
u_mouse <- sort(unique(na.omit(dt_all$mouse_raw)))

if (length(u_mouse) == 4) {
  ## 恰好4个编号：按数值大小映射成Mouse1-4(不丢数据)
  dt_all[, mouse_id := match(mouse_raw, u_mouse)]
  dt_all[, mouse_tag := paste0("C57BL6J", mouse_raw)]
} else {
  ## 否则兜底：按mouse_raw大小切成4组(确保分4只)
  dt_all[!is.na(mouse_raw), rk := frank(mouse_raw, ties.method="first")]
  dt_all[!is.na(rk), mouse_id := as.integer(cut(
    rk,
    breaks = quantile(rk, probs=seq(0,1,0.25), na.rm=TRUE),
    include.lowest=TRUE, labels=1:4
  ))]
  dt_all[, mouse_tag := ifelse(is.na(mouse_raw), NA_character_, paste0("C57BL6J", mouse_raw))]
}

## 只保留成功分到1-4组的切片
dt_all <- dt_all[mouse_id %chin% 1:4]

## 4)按mouse统计1to1总次数，做成标签
dt_mouse_sum <- dt_all[, .(sum_1to1 = sum(n_1to1), n_slides=.N), by=.(mouse_id)][order(mouse_id)]
dt_mouse_sum[, mouse_lab := paste0("Mouse", mouse_id, "(", sum_1to1, "pairs,", n_slides, "slides)")]
fwrite(dt_mouse_sum, file.path(outdir, "mouse_sum_1to1.tsv"), sep="\t", quote=FALSE)

dt_all <- merge(dt_all, dt_mouse_sum[, .(mouse_id, mouse_lab)], by="mouse_id", all.x=TRUE)

## 5)给每个mouse内部做“更合理”的slide排序：优先按slide里最后一个数字排序(没有就按字符串)
get_last_int <- function(x){
  m <- stringr::str_match(x, "([0-9]+)(?!.*[0-9])")
  suppressWarnings(as.integer(m[,2]))
}
dt_all[, slide_lastnum := get_last_int(slide)]
setorder(dt_all, mouse_id, slide_lastnum, slide)

## 6)统一色阶(跨4张图可比) + 更显著颜色对比 + sqrt变换增强低值差异
global_max <- suppressWarnings(max(dt_all$n_1to1, na.rm=TRUE))
if (!is.finite(global_max) || global_max < 1) global_max <- 1

fill_scale_strong <- scale_fill_gradientn(
  colours = c("white", "#ffffcc", "#ffeda0", "#feb24c", "#fd8d3c", "#f03b20", "#bd0026"),
  limits  = c(0, global_max),
  oob     = scales::squish,
  breaks  = pretty(c(0, global_max), n=5),
  labels  = scales::comma,
  trans   = "sqrt",
  name    = "1to1count"
)

## 7)分别输出4张热图(每只鼠一张)：做成“单行热图(strip heatmap)”
for (mid in 1:4) {
  dtm <- dt_all[mouse_id == mid]
  if (nrow(dtm) == 0) next
  
  ## 单行y
  dtm[, y := dt_mouse_sum[mouse_id == mid, mouse_lab][1]]
  dtm[, slide_f := factor(slide, levels=unique(slide))]
  
  p_hm <- ggplot(dtm, aes(x=slide_f, y=y, fill=n_1to1)) +
    geom_tile(color="black", linewidth=0.18) +
    fill_scale_strong +
    labs(
      x=NULL, y=NULL,
      title=paste0("Heatmap:Strict1-to-1pairs per slide - Mouse", mid)
    ) +
    theme_classic(base_size=12) +
    theme(
      plot.title  = element_text(hjust=0.5, face="bold"),
      axis.text.y = element_text(face="bold"),
      axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=7),
      panel.border = element_rect(color="black", fill=NA, linewidth=0.8),
      plot.margin = margin(6, 6, 6, 6),
      legend.position = "right"
    )
  
  ## 动态宽度：切片多就加宽
  w_hm <- max(10, min(80, 0.30 * nlevels(dtm$slide_f)))
  ggsave(file.path(outdir, paste0("Heatmap_Mouse", mid, "_1to1_bySlide.png")),
         p_hm, width=w_hm, height=2.6, dpi=350, bg="white")
  ggsave(file.path(outdir, paste0("Heatmap_Mouse", mid, "_1to1_bySlide.pdf")),
         p_hm, width=w_hm, height=2.6, device=cairo_pdf)
}

## 8)导出明细表(含0次)
fwrite(dt_all[, .(slide, mouse_raw, mouse_id, mouse_tag, mouse_lab, n_1to1)],
       file.path(outdir, "allslides_1to1_count_byMouse.tsv"),
       sep="\t", quote=FALSE)

## =====================================================================
## Strict1-to-1(ov100) pairs：pair_type(subclass1→subclass2)×mouse×layer 精细统计+可视化
## 依赖：pairs_1to1, outdir 已存在
## 输出：
##  1) 多个TSV统计表
##  2) Heatmap_TopPairTypes_byLayer_byMouse.png/pdf   (核心：pairtype×layer，按mouse分面)
##  3) Heatmap_TopPairTypes_byMouse.png/pdf          (核心：pairtype×mouse)
## =====================================================================
stopifnot(exists("pairs_1to1"), exists("outdir"))

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(Cairo)
})

## -------------------- 参数：你可以按需要调 --------------------
topK_pair <- 35L        # 热图只展示TopK的pair_type，其余合并为Other（35通常够清晰）
label_thres <- 5L       # tile里显示数字的阈值：N>=5才标注，避免太挤
layer_levels <- c("1","2/3","4","5","6a","6b")

## -------------------- 工具函数 --------------------
get_mouse_raw <- function(x){
  m <- stringr::str_match(x, "C57BL6J\\D*([0-9]+)")
  suppressWarnings(as.integer(m[,2]))
}

norm_layer <- function(x){
  x <- tolower(gsub("^\\s*layer\\s*", "", as.character(x), ignore.case=TRUE))
  x <- gsub("\\s+", "", x)
  x[x %chin% tolower(layer_levels)] <- layer_levels[match(x[x %chin% tolower(layer_levels)], tolower(layer_levels))]
  x
}

strong_fill <- function(maxv){
  scale_fill_gradientn(
    colours = c("white", "#ffffcc", "#ffeda0", "#feb24c", "#fd8d3c", "#f03b20", "#bd0026"),
    limits  = c(0, maxv),
    oob     = scales::squish,
    trans   = "sqrt",
    breaks  = pretty(c(0, maxv), n=6),
    labels  = comma,
    name    = "Count"
  )
}

theme_pub <- function(base=12){
  theme_classic(base_size=base) +
    theme(
      plot.title = element_text(hjust=0.5, face="bold"),
      axis.text  = element_text(face="bold"),
      axis.title = element_text(face="bold"),
      panel.border = element_rect(color="black", fill=NA, linewidth=0.8),
      legend.title = element_text(face="bold"),
      strip.text = element_text(face="bold")
    )
}

## -------------------- 1) 整理pairs_1to1：补slide/mouse/layer/pair_type --------------------
dt <- as.data.table(pairs_1to1)

## slide优先级：slide > slide2 > slide1；再不行就从label1截取
dt[, slide_use := fifelse(!is.na(slide)  & slide  != "", as.character(slide),
                          fifelse(!is.na(slide2) & slide2 != "", as.character(slide2),
                                  fifelse(!is.na(slide1) & slide1 != "", as.character(slide1),
                                          fifelse(!is.na(label1) & label1 != "", sub("_.*$", "", as.character(label1)), NA_character_))))]

## layer：默认用被包含者(Gaba)的layer2_f/layer2；缺失再用layer1_f/layer1兜底
dt[, layer_use := fifelse(!is.na(layer2_f) & as.character(layer2_f) != "", as.character(layer2_f),
                          fifelse(!is.na(layer2)   & as.character(layer2)   != "", as.character(layer2),
                                  fifelse(!is.na(layer1_f) & as.character(layer1_f) != "", as.character(layer1_f),
                                          fifelse(!is.na(layer1)   & as.character(layer1)   != "", as.character(layer1), NA_character_))))]

dt[, layer_use := norm_layer(layer_use)]
dt[, layer_f := factor(layer_use, levels=layer_levels)]

## mouse_raw + 映射成Mouse1-4
dt[, mouse_raw := get_mouse_raw(slide_use)]
u_mouse <- sort(unique(na.omit(dt$mouse_raw)))

if (length(u_mouse) == 4) {
  dt[, mouse_id := match(mouse_raw, u_mouse)]   # 按数值大小映射成1-4
} else {
  ## 兜底：按mouse_raw排序切4组（保证分成4只）
  dt[!is.na(mouse_raw), rk := frank(mouse_raw, ties.method="first")]
  dt[!is.na(rk), mouse_id := as.integer(cut(
    rk,
    breaks = quantile(rk, probs=seq(0,1,0.25), na.rm=TRUE),
    include.lowest=TRUE, labels=1:4
  ))]
}

dt <- dt[mouse_id %chin% 1:4]
dt[, mouse_f := factor(paste0("Mouse", mouse_id), levels=paste0("Mouse", 1:4))]

## pair_type：subclass1→subclass2（若缺失给NA标签避免丢）
dt[, subclass1 := fifelse(!is.na(subclass1) & subclass1!="", as.character(subclass1), "NA")]
dt[, subclass2 := fifelse(!is.na(subclass2) & subclass2!="", as.character(subclass2), "NA")]
dt[, pair_type := paste0(subclass1, "→", subclass2)]

## sanity
cat("Strict1-to-1pairs used:", nrow(dt), "\n")

## -------------------- 2) 导出“完美清晰”的统计表 --------------------
## 2.1 overall：pair_type分布
tab_pair_overall <- dt[, .N, by=pair_type][order(-N)]
fwrite(tab_pair_overall, file.path(outdir, "Strict1to1_pairtype_overall.tsv"), sep="\t", quote=FALSE)

## 2.2 mouse×layer×pair_type（最精细）
tab_pair_mouse_layer <- dt[, .N, by=.(mouse_f, layer_f, pair_type)][order(mouse_f, layer_f, -N)]
fwrite(tab_pair_mouse_layer, file.path(outdir, "Strict1to1_pairtype_byMouse_byLayer.tsv"), sep="\t", quote=FALSE)

## 2.3 mouse×pair_type
tab_pair_mouse <- dt[, .N, by=.(mouse_f, pair_type)][order(mouse_f, -N)]
fwrite(tab_pair_mouse, file.path(outdir, "Strict1to1_pairtype_byMouse.tsv"), sep="\t", quote=FALSE)

## 2.4 layer×pair_type
tab_pair_layer <- dt[!is.na(layer_f), .N, by=.(layer_f, pair_type)][order(layer_f, -N)]
fwrite(tab_pair_layer, file.path(outdir, "Strict1to1_pairtype_byLayer.tsv"), sep="\t", quote=FALSE)

## 2.5 边缘分布：Glut端subclass1、Gaba端subclass2（mouse×layer）
tab_sub1_mouse_layer <- dt[!is.na(layer_f), .N, by=.(mouse_f, layer_f, subclass1)][order(mouse_f, layer_f, -N)]
tab_sub2_mouse_layer <- dt[!is.na(layer_f), .N, by=.(mouse_f, layer_f, subclass2)][order(mouse_f, layer_f, -N)]
fwrite(tab_sub1_mouse_layer, file.path(outdir, "Strict1to1_subclass1_byMouse_byLayer.tsv"), sep="\t", quote=FALSE)
fwrite(tab_sub2_mouse_layer, file.path(outdir, "Strict1to1_subclass2_byMouse_byLayer.tsv"), sep="\t", quote=FALSE)

## 2.6 每只mouse的总pairs数（用于标题/说明）
tab_mouse_sum <- dt[, .(pairs_1to1=.N, n_slides=uniqueN(slide_use)), by=mouse_f][order(mouse_f)]
fwrite(tab_mouse_sum, file.path(outdir, "Strict1to1_mouse_summary.tsv"), sep="\t", quote=FALSE)

## -------------------- 3) 核心可视化A：TopPairTypes×Layer，按Mouse分面（强烈推荐） --------------------
top_pairs <- tab_pair_overall[1:min(topK_pair, .N), pair_type]
dt_hm <- dt[!is.na(layer_f), .N, by=.(mouse_f, layer_f, pair_type)]
dt_hm[, pair_show := fifelse(pair_type %chin% top_pairs, pair_type, "Other")]
dt_hm <- dt_hm[, .(N=sum(N)), by=.(mouse_f, layer_f, pair_show)]

pair_levels <- c(top_pairs, "Other")
dt_hm[, pair_f := factor(pair_show, levels=rev(pair_levels))]

maxv1 <- suppressWarnings(max(dt_hm$N, na.rm=TRUE)); if (!is.finite(maxv1) || maxv1 < 1) maxv1 <- 1
dt_hm[, label := ifelse(N >= label_thres, as.character(N), "")]

p_hm1 <- ggplot(dt_hm, aes(x=layer_f, y=pair_f, fill=N)) +
  geom_tile(color="black", linewidth=0.18) +
  geom_text(aes(label=label), size=2.8, fontface="bold") +
  facet_wrap(~mouse_f, ncol=1, scales="free_y") +
  strong_fill(maxv1) +
  labs(
    x="Layer(Gaba所在层)", y="PairType(subclass1→subclass2)",
    title=paste0("Strict1to1配对分布:Top", min(topK_pair, length(top_pairs)), "PairType×Layer(按Mouse分面)")
  ) +
  theme_pub(11) +
  theme(
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(size=7),
    strip.background = element_rect(color="black", fill="white", linewidth=0.7),
    plot.margin = margin(6, 6, 6, 6)
  )

h1 <- max(8, 0.22 * length(pair_levels) + 6.5)  # 动态高度
ggsave(file.path(outdir, "Heatmap_TopPairTypes_byLayer_byMouse.png"), p_hm1,
       width=8.8, height=h1, dpi=350, bg="white")
ggsave(file.path(outdir, "Heatmap_TopPairTypes_byLayer_byMouse.pdf"), p_hm1,
       width=8.8, height=h1, device=cairo_pdf)

## -------------------- 4) 核心可视化B：TopPairTypes×Mouse（更“总览”） --------------------
dt_hm2 <- dt[, .N, by=.(mouse_f, pair_type)]
dt_hm2[, pair_show := fifelse(pair_type %chin% top_pairs, pair_type, "Other")]
dt_hm2 <- dt_hm2[, .(N=sum(N)), by=.(mouse_f, pair_show)]
dt_hm2[, pair_f := factor(pair_show, levels=rev(pair_levels))]
dt_hm2[, label := ifelse(N >= label_thres, as.character(N), "")]

maxv2 <- suppressWarnings(max(dt_hm2$N, na.rm=TRUE)); if (!is.finite(maxv2) || maxv2 < 1) maxv2 <- 1

p_hm2 <- ggplot(dt_hm2, aes(x=mouse_f, y=pair_f, fill=N)) +
  geom_tile(color="black", linewidth=0.18) +
  geom_text(aes(label=label), size=3.0, fontface="bold") +
  strong_fill(maxv2) +
  labs(
    x=NULL, y="PairType(subclass1→subclass2)",
    title=paste0("Strict1to1配对分布:Top", min(topK_pair, length(top_pairs)), "PairType×Mouse")
  ) +
  theme_pub(12) +
  theme(
    axis.text.y = element_text(size=8),
    plot.margin = margin(6, 6, 6, 6)
  )

h2 <- max(6, 0.20 * length(pair_levels) + 2.5)
ggsave(file.path(outdir, "Heatmap_TopPairTypes_byMouse.png"), p_hm2,
       width=6.8, height=h2, dpi=350, bg="white")
ggsave(file.path(outdir, "Heatmap_TopPairTypes_byMouse.pdf"), p_hm2,
       width=6.8, height=h2, device=cairo_pdf)

cat("Done. Outputs written to:", outdir, "\n")



rm(list=ls()); gc()

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(Cairo)
})

## ===================== 0)路径 =====================
base_dir <- "E:/zaw/2511/mouseMerfish_zhuang_subclass"
ws_env <- Sys.getenv("WINDOW_SIZE", "0.4")
ss_env <- Sys.getenv("STEP_SIZE",  "0.02")

fmt <- function(x){
  s <- formatC(as.numeric(x), format="f", digits=6)
  s <- sub("0+$","",s); s <- sub("\\.$","",s)
  s
}
combo_root <- file.path(base_dir, paste0("ws", fmt(ws_env), "_ss", fmt(ss_env)))
slice_dir  <- file.path(base_dir, "neocortex_new")

pairs_path <- file.path(combo_root, "fig_glut_gaba_overlap100", "pairs_ov100_strict_1to1.tsv")
stopifnot(file.exists(pairs_path))
stopifnot(dir.exists(slice_dir))

clu_path <- file.path(combo_root, "mouse_subclass_cluster_total_over3.txt")
if (!file.exists(clu_path)) clu_path <- file.path(combo_root, "mouse_subclass_cluster_total.txt")
stopifnot(file.exists(clu_path))

outdir <- file.path(combo_root, "Mouse1_strict1to1_PLOTS_ROI_MATCHED")
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(outdir, "mouse1_all_strict1to1"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(outdir, "tables"), showWarnings=FALSE, recursive=TRUE)

## 绘图抽样（可调）
bg_max_n <- 250000L
hi_max_n <- 400000L
max_pairs_to_plot <- Inf

sample_dt <- function(dt, max_n=200000L, seed=1){
  if (nrow(dt) <= max_n) return(dt)
  set.seed(seed)
  dt[sample.int(.N, max_n)]
}

## ===================== 1)读pairs_1to1，强制只保留Mouse1且ov=1.0 =====================
pairs <- fread(pairs_path, sep="\t", header=TRUE, data.table=TRUE, check.names=FALSE)

for (cc in intersect(c("label1","label2","slide","slide1","slide2","subclass1","subclass2",
                       "ov"), names(pairs))){
  pairs[, (cc) := str_trim(as.character(get(cc)))]
}

pairs[, slide_use := fifelse(!is.na(slide) & slide!="", slide,
                             fifelse(!is.na(slide2) & slide2!="", slide2,
                                     fifelse(!is.na(slide1) & slide1!="", slide1,
                                             fifelse(!is.na(label1) & label1!="", sub("_.*$","",label1), NA_character_))))]

## 只保留Mouse1
pairs <- pairs[!is.na(slide_use) & grepl("^C57BL6J-1\\.", slide_use)]
stopifnot(nrow(pairs) > 0)

## ov归一化到0~1，并严格筛选1.0
stopifnot("ov" %in% names(pairs))
pairs[, ov01 := as.numeric(ov)]
ov_max <- suppressWarnings(max(pairs$ov01, na.rm=TRUE))
if (is.finite(ov_max) && ov_max > 1.5) pairs[, ov01 := ov01/100]
pairs <- pairs[is.finite(ov01) & ov01 >= 0.999999]
stopifnot(nrow(pairs) > 0)
stopifnot(all(pairs$ov01 >= 0.999999))
stopifnot(all(c("label1","label2") %in% names(pairs)))

pairs[, pair_id := .I]
pairs[, pair_type := paste0(subclass1, "→", subclass2)]
pairs <- unique(pairs, by=c("slide_use","label1","label2"))
if (is.finite(max_pairs_to_plot)) pairs <- pairs[1:min(.N, as.integer(max_pairs_to_plot))]
fwrite(pairs, file.path(outdir, "tables", "Mouse1_pairs_strict_ov1.tsv"), sep="\t", quote=FALSE)
cat("Mouse1严格ov=1pairs:", nrow(pairs), "\n")

## ===================== 2)读clu：构建“label -> region所有细胞ID”映射 =====================
clu <- fread(clu_path, sep="\t", header=TRUE, data.table=TRUE, check.names=FALSE)

stopifnot(all(c("slide","layer","subclass","merge_regions",
                "Glut_Neruon_cell_ids","GABA_Neruon_cell_ids","Non_Neruon_cell_ids") %in% names(clu)))

clu[, slide := str_trim(as.character(slide))]
clu[, layer := str_trim(as.character(layer))]
clu[, subclass := str_trim(as.character(subclass))]
clu[, merge_regions := str_trim(as.character(merge_regions))]
if (!"label" %in% names(clu)) {
  clu[, label := paste(slide, layer, subclass, merge_regions, sep="_")]
}
clu[, label := str_trim(as.character(label))]

## 注意：cell_label很长，必须用character匹配，不能转integer
extract_ids <- function(s){
  s <- as.character(s)
  if (is.na(s) || s=="") return(character())
  unique(str_extract_all(s, "\\d+")[[1]])
}

## 关键：直接用label做key，避免merge_regions解析差异导致匹配失败
allcell_map <- clu[
  !is.na(label) & label != "",
  .(all_ids = list(unique(c(
    unlist(lapply(Glut_Neruon_cell_ids, extract_ids)),
    unlist(lapply(GABA_Neruon_cell_ids, extract_ids)),
    unlist(lapply(Non_Neruon_cell_ids,  extract_ids))
  )))),
  by=.(label)
]
setkey(allcell_map, label)

fwrite(allcell_map[, .(label,n_all=lengths(all_ids))],
       file.path(outdir, "tables", "allcell_map_size.tsv"), sep="\t", quote=FALSE)

## ===================== 3)读切片坐标 =====================
file_try <- function(sid){
  f1 <- file.path(slice_dir, paste0(sid, ".txt"))
  f2 <- file.path(slice_dir, paste0(sid, ".tsv"))
  if (file.exists(f1)) return(f1)
  if (file.exists(f2)) return(f2)
  NA_character_
}

read_slice <- function(fp){
  dt <- fread(fp, sep="\t", header=TRUE, data.table=TRUE, showProgress=FALSE,
              select=c("cell_label","x","y","subclass"))
  dt[, cell_label := str_trim(as.character(cell_label))]
  dt[, x := as.numeric(x)]
  dt[, y := as.numeric(y)]
  dt[, subclass := str_trim(as.character(subclass))]
  dt
}

## ===================== 4)轮廓(可选)：给ROI画一个凸包边界线 =====================
get_hull <- function(dt_xy){
  dt_xy <- unique(dt_xy[is.finite(x) & is.finite(y), .(x,y)])
  if (nrow(dt_xy) < 3) return(NULL)
  idx <- chull(dt_xy$x, dt_xy$y)
  hull <- dt_xy[idx]
  rbind(hull, hull[1])
}

collapse_subclass <- function(dt_roi, max_levels=25L){
  tab <- dt_roi[, .N, by=subclass][order(-N)]
  if (nrow(tab) <= max_levels) {
    dt_roi[, subclass_show := subclass]
    return(dt_roi)
  }
  keep <- tab[1:max_levels, subclass]
  dt_roi[, subclass_show := ifelse(subclass %in% keep, subclass, "Other")]
  dt_roi[, subclass_show := factor(subclass_show, levels=c(keep, "Other"))]
  dt_roi
}

## ===================== 5)画图：图1双cluster色块(体现包含)；图2同ROI按subclass =====================
plot_two_regions <- function(dt_slice, ids_glu, ids_gab, title_main){
  bg <- sample_dt(dt_slice[is.finite(x) & is.finite(y), .(x,y)], bg_max_n, seed=1)
  
  hi <- dt_slice[cell_label %in% c(ids_glu, ids_gab), .(x,y,cell_label)]
  if (nrow(hi)==0) return(NULL)
  
  hi[, in_glu := cell_label %in% ids_glu]
  hi[, in_gab := cell_label %in% ids_gab]
  hi[, group := fifelse(in_glu & in_gab, "Overlap",
                        fifelse(in_gab, "GabaOnly",
                                fifelse(in_glu, "GlutOnly", NA_character_)))]
  hi <- hi[!is.na(group)]
  hi[, group := factor(group, levels=c("GlutOnly","Overlap","GabaOnly"))]
  hi <- hi[, .SD[sample.int(.N, min(.N, hi_max_n))], by=group]
  
  ggplot() +
    geom_point(data=bg, aes(x=x, y=y), color="grey80", alpha=0.25, size=0.55) +
    geom_point(data=hi, aes(x=x, y=y, color=group), alpha=0.95, size=0.80) +
    coord_equal() +
    scale_color_manual(values=c("GlutOnly"="#2C7FB8", "Overlap"="#7A5195", "GabaOnly"="#E31A1C"), name=NULL) +
    labs(title=title_main, x=NULL, y=NULL) +
    theme_void(base_size=14) +
    theme(
      plot.background = element_rect(fill="white", color=NA),
      legend.text = element_text(face="bold"),
      plot.title = element_text(face="bold", hjust=0.5, size=18),
      plot.margin = margin(10,10,10,10)
    )
}

plot_roi_subclass_sameROI <- function(dt_slice, roi_ids, title_main){
  bg <- sample_dt(dt_slice[is.finite(x) & is.finite(y), .(x,y)], bg_max_n, seed=1)
  
  roi <- dt_slice[cell_label %in% roi_ids, .(x,y,subclass)]
  if (nrow(roi)==0) return(NULL)
  
  roi <- collapse_subclass(roi, max_levels=25L)
  roi <- roi[, .SD[sample.int(.N, min(.N, hi_max_n))], by=subclass_show]
  
  ## 轮廓线只做展示，不参与筛选（保证ROI=图1一致）
  hull <- get_hull(roi[,.(x,y)])
  
  ggplot() +
    geom_point(data=bg, aes(x=x, y=y), color="grey80", alpha=0.25, size=0.55) +
    geom_point(data=roi, aes(x=x, y=y, color=subclass_show), alpha=0.95, size=0.80) +
    {if (!is.null(hull)) geom_path(data=hull, aes(x=x, y=y), linewidth=0.6, color="black", alpha=0.6)} +
    coord_equal() +
    labs(title=title_main, x=NULL, y=NULL, color=NULL) +
    theme_void(base_size=14) +
    theme(
      plot.background = element_rect(fill="white", color=NA),
      legend.text = element_text(face="bold", size=8),
      plot.title = element_text(face="bold", hjust=0.5, size=18),
      plot.margin = margin(10,10,10,10),
      legend.position="right"
    ) +
    guides(color=guide_legend(ncol=2, override.aes=list(size=2.5, alpha=1)))
}

## ===================== 6)主循环：逐pair逐切片出图，并输出解释用QC =====================
run_one_set <- function(dt_pairs, subdir, tag){
  if (nrow(dt_pairs)==0) return(NULL)
  
  log_dt <- list()
  
  for (sid in unique(dt_pairs$slide_use)){
    fp <- file_try(sid)
    if (is.na(fp)) next
    dt_slice <- read_slice(fp)
    
    rows <- dt_pairs[slide_use==sid]
    for (i in seq_len(nrow(rows))){
      rr <- rows[i]
      
      ## 用label直接拿“region所有细胞ID”(三类并集)
      g_glu <- allcell_map[J(as.character(rr$label1)), all_ids]
      g_gab <- allcell_map[J(as.character(rr$label2)), all_ids]
      
      if (length(g_glu)==0 || length(g_gab)==0) {
        log_dt[[length(log_dt)+1]] <- data.table(
          pair_id=rr$pair_id, slide=sid, note="MissingAllIDs",
          label1=as.character(rr$label1), label2=as.character(rr$label2)
        )
        next
      }
      
      ids_glu <- g_glu[[1]]
      ids_gab <- g_gab[[1]]
      roi_ids <- union(ids_glu, ids_gab)
      
      ## 图1：双cluster色块（注意：图2会用同一个roi_ids）
      title1 <- paste0(sid," | ",tag," | pair_id=",rr$pair_id,
                       " | ov=",sprintf("%.3f", rr$ov01),
                       " | nGlu=",length(ids_glu)," nGab=",length(ids_gab))
      p1 <- plot_two_regions(dt_slice, ids_glu, ids_gab, title1)
      if (!is.null(p1)) {
        ggsave(file.path(outdir, subdir, paste0(sid,"_pair",rr$pair_id,"_A_twoRegions.png")),
               p1, width=10.5, height=10.5, dpi=350, bg="white")
        ggsave(file.path(outdir, subdir, paste0(sid,"_pair",rr$pair_id,"_A_twoRegions.pdf")),
               p1, width=10.5, height=10.5, device=cairo_pdf)
      }
      
      ## 图2：同一ROI(roi_ids)按subclass上色 —— 这张图的区域必然与图1一致
      title2 <- paste0(sid," | ",tag," ROI-subclass(SameROI) | pair_id=",rr$pair_id,
                       " | ROIcells=",length(roi_ids))
      p2 <- plot_roi_subclass_sameROI(dt_slice, roi_ids, title2)
      if (!is.null(p2)) {
        ggsave(file.path(outdir, subdir, paste0(sid,"_pair",rr$pair_id,"_B_ROI_subclass.png")),
               p2, width=11.5, height=10.5, dpi=350, bg="white")
        ggsave(file.path(outdir, subdir, paste0(sid,"_pair",rr$pair_id,"_B_ROI_subclass.pdf")),
               p2, width=11.5, height=10.5, device=cairo_pdf)
      }
      
      log_dt[[length(log_dt)+1]] <- data.table(
        pair_id=rr$pair_id, slide=sid, note="OK",
        label1=as.character(rr$label1), label2=as.character(rr$label2),
        n_glu=length(ids_glu), n_gab=length(ids_gab),
        ROIcells=length(roi_ids), ov=rr$ov01
      )
    }
  }
  
  if (length(log_dt) == 0) {
    return(data.table(
      pair_id=integer(), slide=character(), note=character(),
      label1=character(), label2=character(),
      n_glu=integer(), n_gab=integer(), ROIcells=integer(), ov=numeric()
    ))
  }
  rbindlist(log_dt, use.names=TRUE, fill=TRUE)
}

log_all <- run_one_set(pairs, "mouse1_all_strict1to1", "Mouse1 strict1to1 Glut-contains-Gaba")

fwrite(log_all, file.path(outdir, "tables", "draw_log.tsv"), sep="\t", quote=FALSE)
cat("完成。输出目录：", outdir, "\n", sep="")

## ===================== 7) Mouse1: 每个strict pair的ROI内subclass统计 + 最少subclass筛选 =====================
top_n_base <- 10L
top_n_cap  <- 25L   # 可放宽上限：含并列时最多保留这么多

pair_subclass_list <- list()
pair_summary_list  <- list()
pair_issue_list    <- list()

for (sid in unique(pairs$slide_use)) {
  fp <- file_try(sid)
  rows <- pairs[slide_use == sid]
  
  if (is.na(fp)) {
    pair_issue_list[[length(pair_issue_list) + 1L]] <- rows[, .(
      pair_id, slide=slide_use, label1, label2, issue="SliceFileNotFound"
    )]
    next
  }
  
  dt_slice <- read_slice(fp)
  dt_slice[is.na(subclass) | subclass == "", subclass := "Unknown"]
  
  for (i in seq_len(nrow(rows))) {
    rr <- rows[i]
    
    g_glu <- allcell_map[J(as.character(rr$label1)), all_ids]
    g_gab <- allcell_map[J(as.character(rr$label2)), all_ids]
    if (length(g_glu) == 0 || length(g_gab) == 0) {
      pair_issue_list[[length(pair_issue_list) + 1L]] <- data.table(
        pair_id=rr$pair_id, slide=sid, label1=as.character(rr$label1), label2=as.character(rr$label2),
        issue="MissingAllIDs"
      )
      next
    }
    
    ids_glu <- g_glu[[1]]
    ids_gab <- g_gab[[1]]
    roi_ids <- union(ids_glu, ids_gab)
    dt_roi  <- dt_slice[cell_label %in% roi_ids, .(subclass)]
    
    if (nrow(dt_roi) == 0) {
      pair_issue_list[[length(pair_issue_list) + 1L]] <- data.table(
        pair_id=rr$pair_id, slide=sid, label1=as.character(rr$label1), label2=as.character(rr$label2),
        issue="NoROICellsInSlice"
      )
      next
    }
    
    sub_tab <- dt_roi[, .N, by=subclass][order(-N, subclass)]
    n_roi <- nrow(dt_roi)
    n_sub <- nrow(sub_tab)
    sub_tab[, `:=`(
      pair_id    = as.integer(rr$pair_id),
      slide      = sid,
      label1     = as.character(rr$label1),
      label2     = as.character(rr$label2),
      n_roi_cells= as.integer(n_roi),
      n_subclass = as.integer(n_sub),
      pct        = as.numeric(N / n_roi)
    )]
    setcolorder(sub_tab, c("pair_id","slide","label1","label2","n_roi_cells","n_subclass","subclass","N","pct"))
    
    pair_subclass_list[[length(pair_subclass_list) + 1L]] <- sub_tab
    pair_summary_list[[length(pair_summary_list) + 1L]] <- data.table(
      pair_id      = as.integer(rr$pair_id),
      slide        = sid,
      label1       = as.character(rr$label1),
      label2       = as.character(rr$label2),
      subclass1    = as.character(rr$subclass1),
      subclass2    = as.character(rr$subclass2),
      ov01         = as.numeric(rr$ov01),
      n_glu_cells  = as.integer(length(ids_glu)),
      n_gab_cells  = as.integer(length(ids_gab)),
      n_roi_cells  = as.integer(n_roi),
      n_subclass   = as.integer(n_sub)
    )
  }
}

pair_subclass_detail <- if (length(pair_subclass_list) > 0) {
  rbindlist(pair_subclass_list, use.names=TRUE, fill=TRUE)
} else {
  data.table(pair_id=integer(), slide=character(), label1=character(), label2=character(),
             n_roi_cells=integer(), n_subclass=integer(), subclass=character(), N=integer(), pct=numeric())
}

pair_subclass_summary <- if (length(pair_summary_list) > 0) {
  rbindlist(pair_summary_list, use.names=TRUE, fill=TRUE)
} else {
  data.table(pair_id=integer(), slide=character(), label1=character(), label2=character(),
             subclass1=character(), subclass2=character(), ov01=numeric(),
             n_glu_cells=integer(), n_gab_cells=integer(), n_roi_cells=integer(), n_subclass=integer())
}

pair_issues <- if (length(pair_issue_list) > 0) {
  rbindlist(pair_issue_list, use.names=TRUE, fill=TRUE)
} else {
  data.table(pair_id=integer(), slide=character(), label1=character(), label2=character(), issue=character())
}

fwrite(pair_subclass_detail,
       file.path(outdir, "tables", "Mouse1_pair_roi_subclass_detail.tsv"),
       sep="\t", quote=FALSE)
fwrite(pair_subclass_summary[order(n_subclass, n_roi_cells, -ov01, pair_id)],
       file.path(outdir, "tables", "Mouse1_pair_roi_subclass_summary.tsv"),
       sep="\t", quote=FALSE)
fwrite(pair_issues,
       file.path(outdir, "tables", "Mouse1_pair_roi_subclass_issues.tsv"),
       sep="\t", quote=FALSE)

## 按“subclass种类最少”筛候选：前10起步，含并列，最多放宽到top_n_cap
if (nrow(pair_subclass_summary) > 0) {
  cand_src <- pair_subclass_summary[order(n_subclass, n_roi_cells, -ov01, pair_id)]
  k <- min(top_n_base, nrow(cand_src))
  cut_sub <- cand_src[k, n_subclass]
  cand <- cand_src[n_subclass <= cut_sub]
  if (nrow(cand) > top_n_cap) cand <- cand[1:top_n_cap]
  cand[, rank_min_subclass := .I]
  
  fwrite(cand,
         file.path(outdir, "tables", "Mouse1_pair_candidates_min_subclass.tsv"),
         sep="\t", quote=FALSE)
  
  cand_detail <- merge(
    pair_subclass_detail,
    cand[, .(pair_id, rank_min_subclass)],
    by="pair_id",
    all.y=TRUE
  )[order(rank_min_subclass, pair_id, -N, subclass)]
  
  fwrite(cand_detail,
         file.path(outdir, "tables", "Mouse1_pair_candidates_min_subclass_detail.tsv"),
         sep="\t", quote=FALSE)
  
  ## 候选pairs单独绘图（独立目录 + 清晰命名）
  cand_plot_dir <- file.path(outdir, "mouse1_minSubclass_candidates_plots")
  dir.create(cand_plot_dir, showWarnings=FALSE, recursive=TRUE)
  
  safe_token <- function(x){
    x <- as.character(x)
    x[is.na(x) | x == ""] <- "NA"
    x <- gsub("[^A-Za-z0-9._-]+", "_", x)
    x <- gsub("_+", "_", x)
    x
  }
  
  cand_draw_log <- list()
  for (sid in unique(cand$slide)) {
    fp <- file_try(sid)
    rows <- cand[slide == sid]
    
    if (is.na(fp)) {
      cand_draw_log[[length(cand_draw_log) + 1L]] <- rows[, .(
        pair_id, rank_min_subclass, slide, n_subclass, note="SliceFileNotFound"
      )]
      next
    }
    
    dt_slice <- read_slice(fp)
    dt_slice[is.na(subclass) | subclass == "", subclass := "Unknown"]
    
    for (ii in seq_len(nrow(rows))) {
      rr <- rows[ii]
      g_glu <- allcell_map[J(as.character(rr$label1)), all_ids]
      g_gab <- allcell_map[J(as.character(rr$label2)), all_ids]
      
      if (length(g_glu) == 0 || length(g_gab) == 0) {
        cand_draw_log[[length(cand_draw_log) + 1L]] <- data.table(
          pair_id=rr$pair_id, rank_min_subclass=rr$rank_min_subclass, slide=sid, n_subclass=rr$n_subclass,
          note="MissingAllIDs"
        )
        next
      }
      
      ids_glu <- g_glu[[1]]
      ids_gab <- g_gab[[1]]
      roi_ids <- union(ids_glu, ids_gab)
      if (length(roi_ids) == 0) {
        cand_draw_log[[length(cand_draw_log) + 1L]] <- data.table(
          pair_id=rr$pair_id, rank_min_subclass=rr$rank_min_subclass, slide=sid, n_subclass=rr$n_subclass,
          note="EmptyROI"
        )
        next
      }
      
      prefix <- paste0(
        "rank", sprintf("%02d", as.integer(rr$rank_min_subclass)),
        "_pair", sprintf("%04d", as.integer(rr$pair_id)),
        "_nsub", sprintf("%02d", as.integer(rr$n_subclass)),
        "_", safe_token(sid)
      )
      
      t1 <- paste0(
        "[Candidate] rank=", rr$rank_min_subclass,
        " | pair=", rr$pair_id,
        " | nSub=", rr$n_subclass,
        " | ", sid
      )
      p1 <- plot_two_regions(dt_slice, ids_glu, ids_gab, t1)
      if (!is.null(p1)) {
        ggsave(file.path(cand_plot_dir, paste0(prefix, "_A_twoRegions.png")),
               p1, width=10.5, height=10.5, dpi=350, bg="white")
        ggsave(file.path(cand_plot_dir, paste0(prefix, "_A_twoRegions.pdf")),
               p1, width=10.5, height=10.5, device=cairo_pdf)
      }
      
      t2 <- paste0(
        "[Candidate] rank=", rr$rank_min_subclass,
        " | pair=", rr$pair_id,
        " | nSub=", rr$n_subclass,
        " | ROI subclass | ", sid
      )
      p2 <- plot_roi_subclass_sameROI(dt_slice, roi_ids, t2)
      if (!is.null(p2)) {
        ggsave(file.path(cand_plot_dir, paste0(prefix, "_B_ROI_subclass.png")),
               p2, width=11.5, height=10.5, dpi=350, bg="white")
        ggsave(file.path(cand_plot_dir, paste0(prefix, "_B_ROI_subclass.pdf")),
               p2, width=11.5, height=10.5, device=cairo_pdf)
      }
      
      cand_draw_log[[length(cand_draw_log) + 1L]] <- data.table(
        pair_id=rr$pair_id, rank_min_subclass=rr$rank_min_subclass, slide=sid, n_subclass=rr$n_subclass,
        n_glu=length(ids_glu), n_gab=length(ids_gab), ROIcells=length(roi_ids), note="OK"
      )
    }
  }
  
  cand_draw_log_dt <- if (length(cand_draw_log) > 0) {
    rbindlist(cand_draw_log, use.names=TRUE, fill=TRUE)
  } else {
    data.table(
      pair_id=integer(), rank_min_subclass=integer(), slide=character(), n_subclass=integer(),
      n_glu=integer(), n_gab=integer(), ROIcells=integer(), note=character()
    )
  }
  fwrite(cand_draw_log_dt,
         file.path(outdir, "tables", "Mouse1_pair_candidates_min_subclass_draw_log.tsv"),
         sep="\t", quote=FALSE)
  
  cat("ROI subclass统计完成：有效pairs=", nrow(pair_subclass_summary),
      "；候选(少subclass)=", nrow(cand),
      "；cut_n_subclass=", cut_sub,
      "；候选图目录=", cand_plot_dir, "\n", sep="")
} else {
  cat("ROI subclass统计完成，但没有可用pair（请检查issues表）。\n")
}
