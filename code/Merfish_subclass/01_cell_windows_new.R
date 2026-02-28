rm(list = ls()); gc()
setwd('E:/zaw/2511/mouseMerfish_zhuang_subclass')

suppressPackageStartupMessages({
  library(vroom)
  library(stringr)
  library(data.table)
  library(doParallel)
  library(foreach)
  library(dplyr)
})

start_time.t <- Sys.time()

## ================== 参数与工具 ==================
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

cat("[DEBUG] Reading Merfish reference files once...\n")
Merfish_mouse_neocortex_layer_region <- read.delim('./Merfish_mouse_neocortex_layer_region.txt',
                                                   header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
Merfish_brain_cell_type              <- read.delim('./Merfish_brain_cell_type.txt',
                                                   header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

in_dir  <- './neocortex_new/'
file_names <- list.files(in_dir)
cat("[INFO] Found", length(file_names), "files in", in_dir, "\n")

out_dir_base <- 'E:/zaw/2511/mouseMerfish_zhuang_subclass/cell_window'
out_dir <- file.path(out_dir_base, subdir_tag)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat("[INFO] Output dir: ", out_dir, "\n")
cat(sprintf("[INFO] window_size=%s step_size=%s\n", fmt(window_size), fmt(step_size)))

## ================== 开启并行（按切片并行） ==================
n_cores <- max(1L, min(parallel::detectCores(logical = TRUE) - 1L, length(file_names)))
cat("[INFO] Using", n_cores, "cores via doParallel\n")

cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

## 为了避免 data.table 在每个 worker 内再开多线程导致过度竞争，这里限制为 1 线程
parallel::clusterEvalQ(cl, {
  data.table::setDTthreads(1L)
  NULL
})

## ================== 并行处理每个文件 ==================
res_log <- foreach(file_name = file_names,
                   .packages = c("data.table", "stringr")) %dopar% {
                     ## 每个 worker 自己计时
                     start_time <- Sys.time()
                     msg_prefix <- paste0("[", file_name, "] ")
                     
                     fpath <- file.path(in_dir, file_name)
                     if (!file.exists(fpath)) {
                       return(paste0(msg_prefix, "ERROR: File not found: ", fpath))
                     }
                     
                     cat("\n[INFO] Processing file:", file_name, "\n")
                     
                     ## ========= 先只读表头，构造 colClasses，确保 cell_label / cell_id 从源头就是字符 =========
                     header_dt <- read.delim(
                       fpath,
                       nrows            = 0,
                       stringsAsFactors = FALSE,
                       check.names      = FALSE
                     )
                     col_nms <- names(header_dt)
                     col_cls <- rep(NA_character_, length(col_nms))
                     names(col_cls) <- col_nms
                     
                     if ("cell_label" %in% col_nms) {
                       col_cls["cell_label"] <- "character"
                     }
                     if ("cell_id" %in% col_nms) {
                       col_cls["cell_id"] <- "character"
                     }
                     
                     cat("[DEBUG]", msg_prefix, "Reading main data file with colClasses...\n")
                     file_tmp <- read.delim(
                       fpath,
                       stringsAsFactors = FALSE,
                       check.names      = FALSE,
                       colClasses       = col_cls
                     )
                     cat("[DEBUG]", msg_prefix, "file_tmp loaded:",
                         nrow(file_tmp), "rows x", ncol(file_tmp), "cols\n")
                     
                     ## 再保险：强制 cell_label / cell_id 为字符
                     if ("cell_label" %in% names(file_tmp)) {
                       file_tmp$cell_label <- as.character(file_tmp$cell_label)
                     }
                     if ("cell_id" %in% names(file_tmp)) {
                       file_tmp$cell_id <- as.character(file_tmp$cell_id)
                     }
                     
                     ## ================== 合并注释表 ==================
                     cat("[DEBUG]", msg_prefix, "Merging Merfish_mouse_neocortex_layer_region...\n")
                     file_tmp <- merge(file_tmp,
                                       Merfish_mouse_neocortex_layer_region,
                                       by = 'ccf_region_name')
                     cat("[DEBUG]", msg_prefix, "After first merge:", dim(file_tmp)[1], "x", dim(file_tmp)[2], "\n")
                     
                     cat("[DEBUG]", msg_prefix, "Merging Merfish_brain_cell_type...\n")
                     file_tmp <- merge(file_tmp,
                                       Merfish_brain_cell_type,
                                       by = 'class')
                     cat("[DEBUG]", msg_prefix, "After second merge:", dim(file_tmp)[1], "x", dim(file_tmp)[2], "\n")
                     
                     ## 合并之后再次确保 cell_label / cell_id 仍是字符
                     if ("cell_label" %in% names(file_tmp)) {
                       file_tmp$cell_label <- as.character(file_tmp$cell_label)
                     }
                     if ("cell_id" %in% names(file_tmp)) {
                       file_tmp$cell_id <- as.character(file_tmp$cell_id)
                     }
                     
                     file_chose_names <- sub("\\.txt$", "", file_name)
                     
                     ## 转成 data.table
                     setDT(file_tmp)
                     
                     ## ================== 定义滑窗函数（依赖 file_tmp_layer / min_y / step_size / window_size） ==================
                     cell_windows_function <- function(y_idx) {
                       y_low  <- min_y + y_idx * step_size
                       y_high <- y_low + window_size
                       
                       temp1 <- file_tmp_layer[y >= y_low & y < y_high,
                                               .(cell_label, x)]
                       
                       temp1[, cell_label := as.character(cell_label)]
                       
                       n1 <- nrow(temp1)
                       if (n1 <= 10L) return(NULL)
                       
                       o          <- order(temp1$x)
                       x_sorted   <- temp1$x[o]
                       lab_sorted <- as.character(temp1$cell_label[o])
                       
                       min_x <- x_sorted[1L]
                       max_x <- x_sorted[n1]
                       
                       gx_max <- floor((max_x - min_x) / step_size)
                       left   <- 1L
                       right  <- 0L
                       
                       res_list <- vector("list", gx_max + 1L)
                       k <- 0L
                       
                       for (gx in 0:gx_max) {
                         x_low  <- min_x + gx * step_size
                         x_high <- x_low + window_size
                         
                         while (left <= n1 && x_sorted[left] < x_low)          left  <- left + 1L
                         while (right < n1 && x_sorted[right + 1L] < x_high)   right <- right + 1L
                         
                         if (right - left + 1L > 1L) {
                           loc <- paste(x_low, y_low, sep = "_")
                           k   <- k + 1L
                           res_list[[k]] <- data.table(
                             cell_label = lab_sorted[left:right],
                             loc        = loc
                           )
                         }
                         if (left > n1) break
                       }
                       
                       if (k) rbindlist(res_list[1:k]) else NULL
                     }
                     
                     ## ================== 按 layer 做滑窗 ==================
                     layer_number      <- 1L
                     cell_windows_list <- list()
                     
                     for (layer.1 in sort(unique(file_tmp$layer))) {
                       cat("[DEBUG]", msg_prefix, "Processing layer:", layer.1, "\n")
                       file_tmp_layer <- file_tmp[layer == layer.1, ]
                       
                       if (nrow(file_tmp_layer) == 0) {
                         cat("[WARN]", msg_prefix, "Layer", layer.1, "has no data.\n")
                         next
                       }
                       
                       ## 每个 layer 内再次保证 cell_label 为字符
                       if ("cell_label" %in% names(file_tmp_layer)) {
                         file_tmp_layer[, cell_label := as.character(cell_label)]
                       }
                       
                       min_y <- min(file_tmp_layer$y)
                       max_y <- max(file_tmp_layer$y)
                       
                       y_values   <- 0:floor((max_y - min_y) / step_size)
                       all_result <- lapply(y_values, cell_windows_function)
                       all_result <- rbindlist(all_result, fill = TRUE)
                       
                       if (nrow(all_result) > 0) {
                         cat("[DEBUG]", msg_prefix, "Layer", layer.1,
                             "window results:", nrow(all_result), "rows\n")
                         all_result[, layer := layer.1]
                         cell_windows_list[[layer_number]] <- all_result
                         names(cell_windows_list)[layer_number] <- layer.1
                         layer_number <- layer_number + 1L
                       } else {
                         cat("[WARN]", msg_prefix, "Layer", layer.1, "produced no window results.\n")
                       }
                     }
                     
                     ## ================== 汇总所有 layer 的窗口 ==================
                     cat("[DEBUG]", msg_prefix, "Combining layer results...\n")
                     if (length(cell_windows_list) == 0) {
                       warn_msg <- paste0(msg_prefix, "ERROR: No layers produced results! Skip saving.")
                       cat(warn_msg, "\n")
                       return(warn_msg)
                     }
                     
                     cell_windows_layer_d <- rbindlist(cell_windows_list, fill = TRUE)
                     cat("[DEBUG]", msg_prefix, "cell_windows_layer_d dimensions:",
                         nrow(cell_windows_layer_d), "x", ncol(cell_windows_layer_d), "\n")
                     
                     if (ncol(cell_windows_layer_d) == 0) {
                       warn_msg <- paste0(msg_prefix, "ERROR: cell_windows_layer_d has 0 columns! Skip rename/save.")
                       cat(warn_msg, "\n")
                       return(warn_msg)
                     }
                     
                     setnames(cell_windows_layer_d, 1, "cell_label")
                     cell_windows_layer_d[, cell_label := as.character(cell_label)]
                     
                     end_time <- Sys.time()
                     execution_time <- difftime(end_time, start_time, units = "mins")
                     
                     save(file_tmp, cell_windows_layer_d,
                          file = file.path(out_dir, paste0(file_chose_names, '.RData')))
                     
                     log_msg <- paste0(
                       "[DONE] ", file_name, " in ", round(as.numeric(execution_time), 2), " mins | ",
                       "real=",   length(unique(file_tmp$cell_label)),
                       " sample=", length(unique(cell_windows_layer_d$cell_label))
                     )
                     cat(log_msg, "\n")
                     log_msg
                   }

## 关闭并行
parallel::stopCluster(cl)

end_time.t <- Sys.time()
execution_time <- difftime(end_time.t, start_time.t, units = "mins")
cat("[INFO] Total execution time:", execution_time, "mins\n")
cat("[INFO] Script finished.\n")
print(res_log)
