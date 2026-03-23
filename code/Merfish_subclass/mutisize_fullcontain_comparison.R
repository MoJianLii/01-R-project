rm(list = ls()); gc()

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

## ===================== Config =====================
base_dir <- "E:/zaw/2511/mouseMerfish_zhuang_subclass"
combo_ids <- c(
  "ws0.4_ss0.02",
  "ws0.4_ss0.05",
  "ws0.3_ss0.05",
  "ws0.3_ss0.02",
  "ws0.2_ss0.05",
  "ws0.2_ss0.1"
)
rel_pairs_file <- file.path(
  "Mouse1_strict1to1_PLOTS_ROI_MATCHED",
  "tables",
  "Mouse1_strict_pairs_excludeLayer1.tsv"
)

outdir <- file.path(base_dir, "multisize_fullcontain_comparison")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "tables"), showWarnings = FALSE, recursive = TRUE)

## ===================== Helpers =====================
extract_slide <- function(dt) {
  n <- nrow(dt)
  if ("slide_use" %in% names(dt)) return(as.character(dt$slide_use))
  if ("slide" %in% names(dt)) return(as.character(dt$slide))
  if ("slide2" %in% names(dt)) return(as.character(dt$slide2))
  if ("slide1" %in% names(dt)) return(as.character(dt$slide1))
  if ("label1" %in% names(dt)) return(sub("_.*$", "", as.character(dt$label1)))
  rep(NA_character_, n)
}

normalize_overlap01 <- function(x) {
  v <- suppressWarnings(as.numeric(x))
  vmax <- suppressWarnings(max(v, na.rm = TRUE))
  if (is.finite(vmax) && vmax > 1.5) v <- v / 100
  v
}

read_one_combo <- function(combo_id) {
  fp <- file.path(base_dir, combo_id, rel_pairs_file)
  if (!file.exists(fp)) {
    stop(paste0("Missing required file: ", fp))
  }

  dt <- fread(fp, sep = "\t", header = TRUE, data.table = TRUE, check.names = FALSE)
  n_raw <- nrow(dt)

  dt[, slide_use := trimws(extract_slide(dt))]
  dt[slide_use == "", slide_use := NA_character_]

  if ("ov01" %in% names(dt)) {
    dt[, ov01_norm := normalize_overlap01(ov01)]
    dt <- dt[is.finite(ov01_norm) & ov01_norm >= 0.999999]
  }
  n_after_ov <- nrow(dt)

  dt <- dt[!is.na(slide_use)]
  n_after_slide <- nrow(dt)

  dedup_by <- intersect(c("slide_use", "label1", "label2"), names(dt))
  if (length(dedup_by) > 0) {
    dt <- unique(dt, by = dedup_by)
  }
  n_final <- nrow(dt)

  dt[, combo := combo_id]

  stat <- data.table(
    combo = combo_id,
    file_path = fp,
    n_raw = as.integer(n_raw),
    n_after_ov_filter = as.integer(n_after_ov),
    n_after_slide_filter = as.integer(n_after_slide),
    n_pairs_final = as.integer(n_final),
    n_unique_slices = as.integer(uniqueN(dt$slide_use))
  )

  list(data = dt, stat = stat)
}

## ===================== Read six files =====================
res <- lapply(combo_ids, read_one_combo)
pair_list <- lapply(res, `[[`, "data")
stat_list <- lapply(res, `[[`, "stat")

all_pairs <- rbindlist(pair_list, use.names = TRUE, fill = TRUE)
combo_stat <- rbindlist(stat_list, use.names = TRUE, fill = TRUE)
combo_stat <- combo_stat[match(combo_ids, combo)]

fwrite(combo_stat, file.path(outdir, "tables", "combo_input_and_counts.tsv"), sep = "\t", quote = FALSE)

## ===================== Slice presence comparison =====================
slice_presence_long <- unique(all_pairs[, .(combo, slide_use)])
slice_presence_long[, present := 1L]

slice_presence_wide <- dcast(
  slice_presence_long,
  slide_use ~ combo,
  value.var = "present",
  fill = 0L
)

for (cc in combo_ids) {
  if (!cc %in% names(slice_presence_wide)) {
    slice_presence_wide[, (cc) := 0L]
  }
}
setcolorder(slice_presence_wide, c("slide_use", combo_ids))
slice_presence_wide[, n_combos_present := rowSums(.SD), .SDcols = combo_ids]
setorder(slice_presence_wide, -n_combos_present, slide_use)

fwrite(slice_presence_wide, file.path(outdir, "tables", "slice_presence_matrix.tsv"), sep = "\t", quote = FALSE)

## ===================== Slices common to all six =====================
n_combo <- length(combo_ids)
common_slices <- slice_presence_wide[n_combos_present == n_combo, .(slide_use)]
setorder(common_slices, slide_use)

fwrite(common_slices, file.path(outdir, "tables", "slices_present_in_all6.tsv"), sep = "\t", quote = FALSE)

common_pairs <- all_pairs[slide_use %chin% common_slices$slide_use]
setorder(common_pairs, slide_use, combo)
fwrite(common_pairs, file.path(outdir, "tables", "fullcontain_pairs_on_common_slices.tsv"), sep = "\t", quote = FALSE)

common_slice_summary <- common_pairs[, .(
  n_pairs = .N,
  n_unique_label1 = if ("label1" %in% names(common_pairs)) uniqueN(label1) else NA_integer_,
  n_unique_label2 = if ("label2" %in% names(common_pairs)) uniqueN(label2) else NA_integer_
), by = .(slide_use, combo)][order(slide_use, combo)]

fwrite(common_slice_summary, file.path(outdir, "tables", "common_slices_pair_summary_by_combo.tsv"), sep = "\t", quote = FALSE)

## ===================== Slice set overlap between combos =====================
slice_sets <- lapply(combo_ids, function(cc) unique(all_pairs[combo == cc, slide_use]))
names(slice_sets) <- combo_ids

overlap_dt <- rbindlist(lapply(combo_ids, function(c1) {
  s1 <- slice_sets[[c1]]
  rbindlist(lapply(combo_ids, function(c2) {
    s2 <- slice_sets[[c2]]
    n_intersection <- length(intersect(s1, s2))
    n_union <- length(union(s1, s2))
    data.table(
      combo_1 = c1,
      combo_2 = c2,
      n_intersection = as.integer(n_intersection),
      n_union = as.integer(n_union),
      jaccard = ifelse(n_union > 0, n_intersection / n_union, NA_real_)
    )
  }))
}))

fwrite(overlap_dt, file.path(outdir, "tables", "slice_overlap_pairwise.tsv"), sep = "\t", quote = FALSE)

## ===================== Text summary =====================
summary_file <- file.path(outdir, "summary_multisize_fullcontain.txt")
sink(summary_file)
cat("Input combos:\n")
cat(paste0(" - ", combo_ids), sep = "\n")
cat("\n")
cat("Output directory:\n")
cat(outdir, "\n\n")

cat("Per-combo counts:\n")
print(combo_stat)
cat("\n")

cat("Slices present in all 6 combos:\n")
cat("n_common_slices = ", nrow(common_slices), "\n", sep = "")
if (nrow(common_slices) > 0) {
  print(common_slices)
}
cat("\n")
sink()

cat("Done.\n")
cat("Output directory: ", outdir, "\n", sep = "")
cat("Key files:\n")
cat(" - ", file.path(outdir, "tables", "combo_input_and_counts.tsv"), "\n", sep = "")
cat(" - ", file.path(outdir, "tables", "slice_presence_matrix.tsv"), "\n", sep = "")
cat(" - ", file.path(outdir, "tables", "slices_present_in_all6.tsv"), "\n", sep = "")
cat(" - ", file.path(outdir, "tables", "fullcontain_pairs_on_common_slices.tsv"), "\n", sep = "")

## ===================== Visualization & detailed saves =====================
plot_dir <- file.path(outdir, "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

common_slice_dir <- file.path(outdir, "tables", "common_slices_all6_details")
dir.create(common_slice_dir, showWarnings = FALSE, recursive = TRUE)

safe_token <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "NA"
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(nchar(x) == 0, "NA", x)
}

theme_pub <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
      axis.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

save_fig <- function(filename_base, plot_obj, width, height) {
  ggsave(file.path(plot_dir, paste0(filename_base, ".png")),
         plot_obj, width = width, height = height, dpi = 350, bg = "white")
  ggsave(file.path(plot_dir, paste0(filename_base, ".pdf")),
         plot_obj, width = width, height = height, bg = "white")
}

## ---- 图1: 每个size的pair数量 ----
p_pairs <- ggplot(
  combo_stat,
  aes(x = factor(combo, levels = combo_ids), y = n_pairs_final, fill = combo)
) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25, show.legend = FALSE) +
  geom_text(aes(label = comma(n_pairs_final)), vjust = -0.35, size = 4.0, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)), labels = comma) +
  labs(
    x = "Window/Step size",
    y = "Strict full-contain pairs",
    title = "Pair count across six size settings"
  ) +
  theme_pub(13)
save_fig("Fig01_pairs_per_size", p_pairs, width = 9.5, height = 5.2)

## ---- 图2: 每个size的unique切片数量 ----
p_slices <- ggplot(
  combo_stat,
  aes(x = factor(combo, levels = combo_ids), y = n_unique_slices, fill = combo)
) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25, show.legend = FALSE) +
  geom_text(aes(label = comma(n_unique_slices)), vjust = -0.35, size = 4.0, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)), labels = comma) +
  labs(
    x = "Window/Step size",
    y = "Unique slices",
    title = "Unique slices across six size settings"
  ) +
  theme_pub(13)
save_fig("Fig02_unique_slices_per_size", p_slices, width = 9.5, height = 5.2)

## ---- 图3: 切片出现次数分布（1~6）----
dt_ncombo <- slice_presence_wide[, .N, by = n_combos_present][order(n_combos_present)]
dt_ncombo <- dt_ncombo[data.table(n_combos_present = 1:length(combo_ids)), on = "n_combos_present"]
dt_ncombo[is.na(N), N := 0L]

p_ncombo <- ggplot(dt_ncombo, aes(x = factor(n_combos_present), y = N, fill = factor(n_combos_present))) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25, show.legend = FALSE) +
  geom_text(aes(label = comma(N)), vjust = -0.35, size = 4.0, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)), labels = comma) +
  labs(
    x = "In how many size settings a slice appears",
    y = "Number of slices",
    title = "Slice recurrence across six size settings"
  ) +
  theme_pub(13)
save_fig("Fig03_slice_recurrence_distribution", p_ncombo, width = 9.0, height = 5.2)

## ---- 图4: size两两切片重叠Jaccard热图 ----
hm_j <- copy(overlap_dt)
hm_j[, combo_1 := factor(combo_1, levels = combo_ids)]
hm_j[, combo_2 := factor(combo_2, levels = combo_ids)]
hm_j[, label := sprintf("%.2f", jaccard)]

p_jaccard <- ggplot(hm_j, aes(x = combo_1, y = combo_2, fill = jaccard)) +
  geom_tile(color = "black", linewidth = 0.25) +
  geom_text(aes(label = label), size = 3.8, fontface = "bold") +
  scale_fill_gradientn(
    colours = c("#F7FBFF", "#9ECAE1", "#3182BD", "#08519C"),
    limits = c(0, 1),
    labels = label_number(accuracy = 0.01)
  ) +
  labs(
    x = "Size setting",
    y = "Size setting",
    fill = "Jaccard",
    title = "Pairwise slice overlap (Jaccard)"
  ) +
  theme_pub(12)
save_fig("Fig04_jaccard_heatmap", p_jaccard, width = 8.0, height = 6.8)

## ---- 图5: 切片出现矩阵热图（按出现次数排序，最多展示前300个）----
show_n <- min(300L, nrow(slice_presence_wide))
if (show_n > 0) {
  hm_slice <- copy(slice_presence_wide[1:show_n])
  hm_slice[, slice_rank := sprintf("S%03d", .I)]
  
  hm_slice_long <- melt(
    hm_slice,
    id.vars = c("slide_use", "n_combos_present", "slice_rank"),
    measure.vars = combo_ids,
    variable.name = "combo",
    value.name = "present"
  )
  hm_slice_long[, combo := factor(as.character(combo), levels = combo_ids)]
  hm_slice_long[, slice_rank := factor(slice_rank, levels = rev(unique(slice_rank)))]
  
  p_slice_heat <- ggplot(hm_slice_long, aes(x = combo, y = slice_rank, fill = present)) +
    geom_tile(color = "grey85", linewidth = 0.1) +
    scale_fill_gradientn(
      colours = c("white", "#2B8CBE"),
      values = c(0, 1),
      limits = c(0, 1),
      breaks = c(0, 1),
      labels = c("Absent", "Present")
    ) +
    labs(
      x = "Size setting",
      y = paste0("Slice rank (top ", show_n, " by recurrence)"),
      fill = NULL,
      title = "Slice presence matrix"
    ) +
    theme_pub(10) +
    theme(axis.text.y = element_text(size = 6))
  save_fig("Fig05_slice_presence_heatmap_top300", p_slice_heat, width = 8.6, height = 10.5)
  
  fwrite(hm_slice[, .(slice_rank, slide_use, n_combos_present)], 
         file.path(outdir, "tables", "slice_presence_heatmap_top300_index.tsv"),
         sep = "\t", quote = FALSE)
}

## ---- 图6: 6个size共同切片在各size中的pair数热图 ----
if (nrow(common_slice_summary) > 0) {
  dt_common_total <- common_slice_summary[, .(total_pairs = sum(n_pairs)), by = slide_use][order(-total_pairs, slide_use)]
  show_common_n <- min(200L, nrow(dt_common_total))
  top_common <- dt_common_total[1:show_common_n, slide_use]

  hm_common <- common_slice_summary[slide_use %chin% top_common]
  hm_common <- hm_common[data.table(slide_use = top_common), on = "slide_use"]
  hm_common[, combo := factor(combo, levels = combo_ids)]
  hm_common[, slide_use := factor(slide_use, levels = rev(top_common))]
  hm_common[, label := comma(n_pairs)]

  p_common_heat <- ggplot(hm_common, aes(x = combo, y = slide_use, fill = n_pairs)) +
    geom_tile(color = "grey88", linewidth = 0.1) +
    geom_text(aes(label = label), size = 2.7, fontface = "bold") +
    scale_fill_gradientn(
      colours = c("#FFF5F0", "#FCAE91", "#FB6A4A", "#CB181D"),
      trans = "sqrt",
      labels = comma
    ) +
    labs(
      x = "Size setting",
      y = paste0("Common slices (top ", show_common_n, " by total pairs)"),
      fill = "Pairs",
      title = "Pair counts on slices present in all six sizes"
    ) +
    theme_pub(10) +
    theme(axis.text.y = element_text(size = 6))
  save_fig("Fig06_common_slices_paircount_heatmap_top200", p_common_heat, width = 10.2, height = 11.0)
}

## ---- 将“6组均出现切片”按切片逐个保存 ----
if (nrow(common_slices) > 0) {
  for (sid in common_slices$slide_use) {
    dt_one <- common_pairs[slide_use == sid][order(combo)]
    fwrite(
      dt_one,
      file.path(common_slice_dir, paste0("all6_", safe_token(sid), "_pairs.tsv")),
      sep = "\t", quote = FALSE
    )
  }
}

cat("Visualization done.\n")
cat("Plot directory: ", plot_dir, "\n", sep = "")
cat("Common-slice detail directory: ", common_slice_dir, "\n", sep = "")

## ===================== Collect all6 slice images from combined_with_cluster_stats =====================
all6_img_outdir <- file.path(outdir, "all6_common_slice_combined_with_cluster_stats")
dir.create(all6_img_outdir, showWarnings = FALSE, recursive = TRUE)

rel_combined_stats_dir <- file.path(
  "Mouse1_strict1to1_PLOTS_ROI_MATCHED",
  "mouse1_excludeLayer1_refstyle_from_local",
  "combined_with_cluster_stats"
)

scan_log_list <- list()
copy_log_list <- list()

if (nrow(common_slices) > 0) {
  for (cc in combo_ids) {
    src_dir <- file.path(base_dir, cc, rel_combined_stats_dir)
    dst_combo_dir <- file.path(all6_img_outdir, cc)
    dir.create(dst_combo_dir, showWarnings = FALSE, recursive = TRUE)
    
    if (!dir.exists(src_dir)) {
      scan_log_list[[length(scan_log_list) + 1L]] <- data.table(
        combo = cc,
        slide_use = NA_character_,
        source_dir = src_dir,
        n_files_matched = NA_integer_,
        status = "SourceDirMissing"
      )
      next
    }
    
    all_img_files <- list.files(src_dir, full.names = TRUE, recursive = FALSE)
    all_img_files <- all_img_files[grepl("\\.(png|pdf)$", all_img_files, ignore.case = TRUE)]
    if (length(all_img_files) == 0) {
      scan_log_list[[length(scan_log_list) + 1L]] <- data.table(
        combo = cc,
        slide_use = NA_character_,
        source_dir = src_dir,
        n_files_matched = 0L,
        status = "NoImageFiles"
      )
      next
    }
    
    for (sid in common_slices$slide_use) {
      sid <- as.character(sid)
      sid_tok <- safe_token(sid)
      bn <- basename(all_img_files)
      
      hit <- all_img_files[
        startsWith(bn, paste0(sid, "_")) |
          startsWith(bn, paste0(sid_tok, "_"))
      ]
      
      scan_log_list[[length(scan_log_list) + 1L]] <- data.table(
        combo = cc,
        slide_use = sid,
        source_dir = src_dir,
        n_files_matched = as.integer(length(hit)),
        status = ifelse(length(hit) > 0, "OK", "NoMatchedImage")
      )
      
      if (length(hit) > 0) {
        dst_slice_dir <- file.path(dst_combo_dir, safe_token(sid))
        dir.create(dst_slice_dir, showWarnings = FALSE, recursive = TRUE)
        
        dst <- file.path(dst_slice_dir, basename(hit))
        ok <- file.copy(hit, dst, overwrite = TRUE, copy.date = TRUE)
        
        copy_log_list[[length(copy_log_list) + 1L]] <- data.table(
          combo = cc,
          slide_use = sid,
          source_file = hit,
          dest_file = dst,
          copied = as.logical(ok)
        )
      }
    }
  }
}

scan_log <- if (length(scan_log_list) > 0) {
  rbindlist(scan_log_list, use.names = TRUE, fill = TRUE)
} else {
  data.table(
    combo = character(),
    slide_use = character(),
    source_dir = character(),
    n_files_matched = integer(),
    status = character()
  )
}

copy_log <- if (length(copy_log_list) > 0) {
  rbindlist(copy_log_list, use.names = TRUE, fill = TRUE)
} else {
  data.table(
    combo = character(),
    slide_use = character(),
    source_file = character(),
    dest_file = character(),
    copied = logical()
  )
}

fwrite(
  scan_log,
  file.path(outdir, "tables", "all6_common_slices_combined_with_cluster_stats_scan.tsv"),
  sep = "\t", quote = FALSE
)
fwrite(
  copy_log,
  file.path(outdir, "tables", "all6_common_slices_combined_with_cluster_stats_copy_log.tsv"),
  sep = "\t", quote = FALSE
)

img_summary <- scan_log[!is.na(slide_use), .(
  n_common_slices = uniqueN(slide_use),
  n_slices_with_image = sum(n_files_matched > 0, na.rm = TRUE),
  n_matched_files = sum(n_files_matched, na.rm = TRUE)
), by = combo]

if (nrow(copy_log) > 0) {
  copied_summary <- copy_log[, .(
    n_copied_files = sum(copied, na.rm = TRUE),
    n_copy_failed = sum(!copied, na.rm = TRUE)
  ), by = combo]
  img_summary <- merge(img_summary, copied_summary, by = "combo", all = TRUE)
}

img_summary <- img_summary[match(combo_ids, combo)]
fwrite(
  img_summary,
  file.path(outdir, "tables", "all6_common_slices_combined_with_cluster_stats_summary.tsv"),
  sep = "\t", quote = FALSE
)

cat("All6 combined_with_cluster_stats image collection done.\n")
cat("Image output directory: ", all6_img_outdir, "\n", sep = "")

## ===================== Keep slices present in >=2 parameter settings =====================
stopifnot(exists("slice_presence_wide"), exists("all_pairs"), exists("combo_ids"))

atleast2_outdir <- file.path(outdir, "atleast2_params_slice_results")
dir.create(atleast2_outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(atleast2_outdir, "tables"), showWarnings = FALSE, recursive = TRUE)

atleast2_slices <- copy(slice_presence_wide[n_combos_present >= 2L])
only1_slices <- copy(slice_presence_wide[n_combos_present == 1L])

setorder(atleast2_slices, -n_combos_present, slide_use)
setorder(only1_slices, slide_use)

fwrite(
  atleast2_slices,
  file.path(atleast2_outdir, "tables", "slices_present_in_atleast2_of6.tsv"),
  sep = "\t", quote = FALSE
)
fwrite(
  only1_slices,
  file.path(atleast2_outdir, "tables", "slices_present_in_only1_of6.tsv"),
  sep = "\t", quote = FALSE
)

atleast2_pairs <- all_pairs[slide_use %chin% atleast2_slices$slide_use]
setorder(atleast2_pairs, slide_use, combo)
fwrite(
  atleast2_pairs,
  file.path(atleast2_outdir, "tables", "fullcontain_pairs_on_slices_atleast2_of6.tsv"),
  sep = "\t", quote = FALSE
)

atleast2_summary_by_combo <- atleast2_pairs[, .(
  n_pairs = .N,
  n_unique_slices = uniqueN(slide_use),
  n_unique_label1 = if ("label1" %in% names(atleast2_pairs)) uniqueN(label1) else NA_integer_,
  n_unique_label2 = if ("label2" %in% names(atleast2_pairs)) uniqueN(label2) else NA_integer_
), by = combo]
atleast2_summary_by_combo <- atleast2_summary_by_combo[match(combo_ids, combo)]

fwrite(
  atleast2_summary_by_combo,
  file.path(atleast2_outdir, "tables", "summary_atleast2_by_combo.tsv"),
  sep = "\t", quote = FALSE
)

atleast2_meta <- data.table(
  filter_rule = "n_combos_present >= 2",
  total_combos = length(combo_ids),
  total_slices_before = nrow(slice_presence_wide),
  slices_kept = nrow(atleast2_slices),
  slices_removed_only1 = nrow(only1_slices),
  pair_rows_kept = nrow(atleast2_pairs)
)
fwrite(
  atleast2_meta,
  file.path(atleast2_outdir, "tables", "atleast2_filter_metadata.tsv"),
  sep = "\t", quote = FALSE
)

atleast2_summary_txt <- file.path(atleast2_outdir, "summary_atleast2_of6.txt")
sink(atleast2_summary_txt)
cat("Filter rule: keep slices with n_combos_present >= 2 (remove slices appearing in only one parameter setting).\n\n")
cat("Total slices in matrix: ", nrow(slice_presence_wide), "\n", sep = "")
cat("Slices kept (>=2): ", nrow(atleast2_slices), "\n", sep = "")
cat("Slices removed (=1): ", nrow(only1_slices), "\n\n", sep = "")
cat("Pair rows kept after filter: ", nrow(atleast2_pairs), "\n\n", sep = "")
cat("Summary by combo:\n")
print(atleast2_summary_by_combo)
sink()

cat("AtLeast2 filter done.\n")
cat("AtLeast2 output directory: ", atleast2_outdir, "\n", sep = "")
## ===================== End of AtLeast2 filter module =====================
