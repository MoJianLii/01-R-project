# ====== Rcpp版本 - 修复浮点数精度问题 ======
# 关键修复：使用sprintf处理大浮点数cell_label，避免精度损失

# 设置OpenMP线程数
Sys.setenv("OMP_NUM_THREADS" = "60")

Sys.setenv("PKG_CXXFLAGS" = "-std=c++11 -D_GLIBCXX_USE_CXX11_ABI=0 -fopenmp")
Sys.setenv("PKG_LIBS" = "-fopenmp")
Sys.setenv("PKG_CPPFLAGS" = "-I/usr/include")

suppressPackageStartupMessages({
  library(stringr)
  library(data.table)
  library(Rcpp)
})

message(sprintf("OpenMP threads: %s", Sys.getenv("OMP_NUM_THREADS")))

message("Compiling C++ code...")
sourceCpp("permutation_test_corrected.cpp")
message("C++ compilation complete!\n")

# ====== 辅助函数：安全地将cell_label转换为字符串 ======
# 使用sprintf确保大浮点数不丢失精度
safe_cell_label_to_str <- function(x) {
  sprintf("%.0f", x)
}

# ====== 配置 ======
base_dir        <- "/picb/neurosys/chenrenrui/Mouse_ST/A molecularly defined and spatially resolved cell atlas of the whole mouse brain/"
subclass_in_dir    <- file.path(base_dir, "subclass_10K/neocortex_sample/cell_window_layer")
data_sample_dir <- file.path(base_dir, 'subclass_10K/neocortex_sample_subclass/data_sample/')
pvalue_out_dir  <- file.path(base_dir, "subclass_10K/neocortex_sample_subclass/cell_pvalue")
num_iterations  <- 10000L

safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

setwd(base_dir)
safe_dir_create(pvalue_out_dir)

message("╔════════════════════════════════════════════════════════════╗")
message("║    Rcpp版本 - 修复浮点数精度问题（使用sprintf）            ║")
message("╚════════════════════════════════════════════════════════════╝\n")

start_time_all <- Sys.time()

chip <- list.files(file.path(base_dir, "data_filter/neocortex_new"))
chip <- str_split(chip, ".txt", simplify = TRUE)[, 1]

for (chip.1 in chip[101:155]) {
  message(sprintf("\n==> Processing chip: %s", chip.1))

  fn_dir <- file.path(subclass_in_dir, chip.1)
  file_names <- list.files(fn_dir)
  safe_dir_create(file.path(pvalue_out_dir, chip.1))

  for (sample.1 in 1:500) {
    sample_start <- Sys.time()
    message(sprintf("\n[Sample %d/500]", sample.1))

    for (file_name in file_names) {
      file_chose_names <- str_split(file_name, ".RD", simplify = TRUE)[, 1]
      file_chip        <- str_split(file_chose_names, "_layer", simplify = TRUE)[, 1]

      # 加载数据
      load(file.path(subclass_in_dir, file_chip, file_name))
      setDT(file_tmp_layer)
      setDT(cell_windows_layer_d_layer)

      # ====== 数据预处理 ======
      # 保存原始subclass名称用于后续输出
      file_tmp_layer[, subclass_original := subclass]
      
      # 创建subclass的字符串到整数的映射
      subclass_levels <- sort(unique(as.character(file_tmp_layer$subclass)))
      file_tmp_layer[, subclass_factor := factor(as.character(subclass), levels = subclass_levels)]
      file_tmp_layer[, subclass_int := as.integer(subclass_factor)]

      # ★★★ 关键修复：使用sprintf处理cell_label ★★★
      # 创建cell_label字符串（使用sprintf避免浮点数精度损失）
      file_tmp_layer[, cell_label_str := safe_cell_label_to_str(cell_label)]
      cell_windows_layer_d_layer[, cell_label_str := safe_cell_label_to_str(cell_label)]
      
      # 映射cell_label_str到整数ID
      unique_cells_str <- unique(file_tmp_layer$cell_label_str)
      cell_id_map <- seq_along(unique_cells_str)
      names(cell_id_map) <- unique_cells_str
      
      file_tmp_layer[, cell_id := cell_id_map[cell_label_str]]
      cell_windows_layer_d_layer[, cell_id := cell_id_map[cell_label_str]]
      
      # 检查是否有NA
      n_na_cells <- sum(is.na(file_tmp_layer$cell_id))
      n_na_windows <- sum(is.na(cell_windows_layer_d_layer$cell_id))
      if (n_na_cells > 0 || n_na_windows > 0) {
        message(sprintf("  ⚠ Warning: NA in cell_id: cells=%d, windows=%d", n_na_cells, n_na_windows))
      }
      
      # 移除无法映射的窗口
      cell_windows_layer_d_layer <- cell_windows_layer_d_layer[!is.na(cell_id)]

      # 映射loc到整数ID
      unique_locs <- unique(cell_windows_layer_d_layer$loc)
      loc_id_map <- seq_along(unique_locs)
      names(loc_id_map) <- as.character(unique_locs)
      loc_id_reverse <- names(loc_id_map)
      names(loc_id_reverse) <- as.character(loc_id_map)
      cell_windows_layer_d_layer[, loc_id := loc_id_map[as.character(loc)]]

      message(sprintf("  Loaded: %d cells, %d subclasses, %d layers",
                      nrow(file_tmp_layer),
                      length(subclass_levels),
                      length(unique(file_tmp_layer$layer))))

      # ★★★ 打乱class_int ★★★
      file_tmp_layer[, subclass_int := sample(subclass_int, .N)]
      
      # ★★★ 根据打乱后的class_int，反查class名称 ★★★
      file_tmp_layer[, subclass_shuffled := subclass_levels[subclass_int]]

      # ★★★ 保存打乱后的class（使用原始cell_label，不是字符串版本）★★★
      safe_dir_create(file.path(data_sample_dir, chip.1))
      data_sample_file <- file.path(data_sample_dir, chip.1,
                                    paste0('sample_', sample.1, '_', file_chose_names, '.txt'))
      fwrite(file_tmp_layer[, .(cell_label, subclass = subclass_shuffled, layer)],
             file = data_sample_file)
      message(sprintf("  Saved shuffled data_sample: %s", basename(data_sample_file)))

      # ====== 逐个subclass处理 ======
      for (subclass_idx in 1:length(subclass_levels)) {
        subclass_name <- subclass_levels[subclass_idx]
        message(sprintf("  subclass %d/%d: %s", subclass_idx, length(subclass_levels), subclass_name))

        subclass_start <- Sys.time()
        all_results <- list()

        for (layer.1 in sort(unique(file_tmp_layer$layer))) {
          # 筛选当前layer的细胞数据
          file_tmp_sample <- file_tmp_layer[layer == layer.1]

          # ★★★ 获取当前layer的所有窗口（不预筛选！）★★★
          cell_windows_sample <- cell_windows_layer_d_layer[layer == layer.1]

          # 只保留cell_id在file_tmp_sample中的窗口
          cell_windows_sample <- cell_windows_sample[cell_id %in% file_tmp_sample$cell_id]

          if (nrow(cell_windows_sample) == 0) {
            next
          }

          # 准备C++函数所需的向量
          cell_labels_vec <- as.integer(file_tmp_sample$cell_id)
          subclass_vec_int <- as.integer(file_tmp_sample$subclass_int)
          window_cell_labels_vec <- as.integer(cell_windows_sample$cell_id)
          window_locs_vec <- as.integer(cell_windows_sample$loc_id)
          target_subclass_int <- as.integer(subclass_idx)
          n_iter_int <- as.integer(num_iterations)

          # 调用C++
          result <- fast_permutation_test_cpp(
            cell_labels = cell_labels_vec,
            subclass_vec = subclass_vec_int,
            window_cell_labels = window_cell_labels_vec,
            window_locs = window_locs_vec,
            target_subclass = target_subclass_int,
            n_iterations = n_iter_int
          )

          # 处理结果
          setDT(result)
          if (nrow(result) > 0) {
            result[, loc_original := loc_id_reverse[as.character(loc)]]
            result[, subclass := subclass_name]
            result[, layer := layer.1]
            all_results[[as.character(layer.1)]] <-
              result[, .(subclass, loc = loc_original, true_sum, b, p, layer)]
          }
        }

        # 写出结果
        if (length(all_results) > 0) {
          res <- rbindlist(all_results, use.names = TRUE, fill = TRUE)

          out_file <- file.path(pvalue_out_dir, chip.1,
                                paste0('sample_', sample.1, '_', file_chose_names, "_",
                                       gsub("/", ".", subclass_name),
                                       "_cell_sliding_window_result_p.csv"))
          fwrite(res, file = out_file)

          subclass_time <- as.numeric(difftime(Sys.time(), subclass_start, units = "secs"))
          message(sprintf("    Completed in %.1fs (%d windows)", subclass_time, nrow(res)))
        }
      }

      gc()
    }

    sample_time <- as.numeric(difftime(Sys.time(), sample_start, units = "mins"))
    message(sprintf("\n[Sample %d/500] Done: %.1f min", sample.1, sample_time))

    if (sample.1 == 1) {
      est_total <- sample_time * 500 / 60
      message(sprintf("\n⏱  Estimate: %.1f hours (%.1f days)\n",
                      est_total, est_total/24))
    }
  }
}

total_time <- as.numeric(difftime(Sys.time(), start_time_all, units = "hours"))
message("\n╔════════════════════════════════════════════════════════════╗")
message(sprintf("║  DONE: %.1f hours (%.1f days)  ║",
                total_time, total_time/24))
message("╚════════════════════════════════════════════════════════════╝")
