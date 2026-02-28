rm(list = ls()); gc()
suppressPackageStartupMessages({
  library(vroom)
  library(stringr)
  library(data.table)
  library(doParallel)
  library(foreach)
  library(dplyr)
  library(parallel)
})

t0 <- Sys.time()

pval_root   <- "E:/zaw/2511/mouseMerfish_zhuang_subclass/cell_pvalue_new"
cellwin_dir <- "E:/zaw/2511/mouseMerfish_zhuang_subclass/cell_window/ws0.4_ss0.02"
out_root    <- "E:/zaw/2511/mouseMerfish_zhuang_subclass/merge_region"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

file_names <- list.files(pval_root, pattern = "sig_high\\.csv$", full.names = FALSE)

## 去掉可能多余的第一列（你之前的 dummy-col 逻辑）
strip_dummy_col_if_needed <- function(dt){
  need <- c("subclass","loc","p","layer")
  if (!all(need %in% names(dt)) && ncol(dt) >= 2L) {
    dt2 <- dt[, -1, drop = FALSE]
    if (all(need %in% names(dt2))) return(dt2)
  }
  dt
}

WS <- 0.4

# --- 仅用于“文件名”的斜杠替换 ---
safe_name <- function(x) gsub("[/\\\\]", "&", x)

## ========= 并行设置 =========
workers <- max(1L, parallel::detectCores() - 1L)   # 视自己机器情况也可以手动写 8、12 等
cat("[INFO] Using workers:", workers, "\n")

cl <- makeCluster(workers)
doParallel::registerDoParallel(cl)

## 每个 worker 里关闭 BLAS/OMP 多线程 + data.table 线程，避免过度超线程
clusterEvalQ(cl, {
  Sys.setenv(
    OPENBLAS_NUM_THREADS = "1",
    OMP_NUM_THREADS      = "1",
    MKL_NUM_THREADS      = "1",
    GOTO_NUM_THREADS     = "1"
  )
  data.table::setDTthreads(1L)
  NULL
})

## ========= 并行按切片处理 =========
foreach(
  file_name = file_names,
  .packages = c("data.table","stringr"),
  .export   = c("pval_root","cellwin_dir","out_root","WS",
                "strip_dummy_col_if_needed","safe_name")
) %dopar% {
  message("[INFO] file: ", file_name)
  
  ## ========= 读入 p 值文件（强制字符串） =========
  p_csv <- file.path(pval_root, file_name)
  mouse_sig_p.total <- read.csv(
    p_csv,
    check.names      = FALSE,
    stringsAsFactors = FALSE
  )
  mouse_sig_p.total <- strip_dummy_col_if_needed(mouse_sig_p.total)
  stopifnot(all(c("subclass","loc","p","layer") %in% names(mouse_sig_p.total)))
  
  ## subclass / loc / layer 统一为字符
  mouse_sig_p.total$subclass <- as.character(mouse_sig_p.total$subclass)
  mouse_sig_p.total$loc      <- as.character(mouse_sig_p.total$loc)
  mouse_sig_p.total$layer    <- as.character(mouse_sig_p.total$layer)
  
  cell_p_file_name_layer <- unique(mouse_sig_p.total$layer)
  
  ## ========= 读入对应切片的 cell_window RData（强制 cell_label/cell_id 为字符） =========
  chip_id    <- str_split(file_name, "[_]", simplify = TRUE)[,1]
  rdata_path <- file.path(cellwin_dir, paste0(chip_id, ".RData"))
  if (!file.exists(rdata_path)) {
    message("[WARN] RData not found for chip: ", chip_id)
    return(NULL)
  }
  load(rdata_path)  # 通常包含 file_tmp, cell_windows_layer_d
  
  ## 源数据统一字符串
  if (exists("file_tmp")) {
    if ("cell_label" %in% names(file_tmp)) {
      file_tmp$cell_label <- as.character(file_tmp$cell_label)
    }
    if ("cell_id" %in% names(file_tmp)) {
      file_tmp$cell_id <- as.character(file_tmp$cell_id)
    }
    data.table::setDT(file_tmp)
  }
  if (exists("cell_windows_layer_d")) {
    if ("cell_label" %in% names(cell_windows_layer_d)) {
      cell_windows_layer_d$cell_label <- as.character(cell_windows_layer_d$cell_label)
    }
    if ("cell_id" %in% names(cell_windows_layer_d)) {
      cell_windows_layer_d$cell_id <- as.character(cell_windows_layer_d$cell_id)
    }
    data.table::setDT(cell_windows_layer_d)
  }
  
  ## loc_1 用 layer + loc 拼接
  cell_windows_layer_d[, layer := as.character(layer)]
  cell_windows_layer_d[, loc   := as.character(loc)]
  cell_windows_layer_d[, loc_1 := paste(layer, loc, sep = "_")]
  
  ## ========= 按 layer 处理 =========
  for (cell_p_file_name_layer.1 in cell_p_file_name_layer) {
    ## 当前 layer 的显著窗口结果
    cell_sliding_windows_result_p <- mouse_sig_p.total[mouse_sig_p.total$layer %in% cell_p_file_name_layer.1, ]
    if (!nrow(cell_sliding_windows_result_p)) next
    
    ## loc 拆分 x/y
    xy_mat <- str_split(cell_sliding_windows_result_p$loc, "[_]", simplify = TRUE)
    cell_sliding_windows_result_p$xstart <- as.numeric(xy_mat[,1])
    cell_sliding_windows_result_p$ystart <- as.numeric(xy_mat[,2])
    cell_sliding_windows_result_p$xend   <- cell_sliding_windows_result_p$xstart + WS
    cell_sliding_windows_result_p$yend   <- cell_sliding_windows_result_p$ystart + WS
    
    ## 只保留 p<0.05 的富集窗口
    Neruon_enrich <- cell_sliding_windows_result_p[cell_sliding_windows_result_p$p < 0.05, ]
    if (!nrow(Neruon_enrich)) next
    
    Neruon_enrich$layer    <- as.character(Neruon_enrich$layer)
    Neruon_enrich$loc      <- as.character(Neruon_enrich$loc)
    Neruon_enrich$loc_1    <- paste(Neruon_enrich$layer, Neruon_enrich$loc, sep = "_")
    Neruon_enrich$subclass <- as.character(Neruon_enrich$subclass)
    
    Neruon_subclass <- unique(Neruon_enrich$subclass)
    
    for (cls in Neruon_subclass) {
      temp_Neruon_enrich <- Neruon_enrich[Neruon_enrich$subclass %in% cls, ]
      if (!nrow(temp_Neruon_enrich)) next
      
      ## 当前 subclass 在窗口中的所有细胞
      temp_Neruon_enrich_cell_id <- cell_windows_layer_d[
        cell_windows_layer_d$loc_1 %in% temp_Neruon_enrich$loc_1, ]
      if (!nrow(temp_Neruon_enrich_cell_id)) next
      
      ## 再保险：cell_label 全部字符
      temp_Neruon_enrich_cell_id$cell_label <- as.character(temp_Neruon_enrich_cell_id$cell_label)
      
      ## =============== 步骤1：把每个窗口的 cell_label 列表收集起来 ===============
      temp_Neruon_enrich_cell_list <- list()
      a <- 1L
      for (i in unique(temp_Neruon_enrich_cell_id$loc_1)) {
        temp1 <- temp_Neruon_enrich_cell_id[temp_Neruon_enrich_cell_id$loc_1 %in% i, ]
        temp1$cell_label <- as.character(temp1$cell_label)
        temp_Neruon_enrich_cell_list[[a]] <- temp1$cell_label
        names(temp_Neruon_enrich_cell_list)[a] <- unique(paste(temp1$layer, temp1$loc, sep = "_"))
        a <- a + 1L
      }
      
      ## =============== 步骤2：按重叠合并窗口，形成 merged_regions ===============
      merged_regions <- list()
      for (i in seq_along(temp_Neruon_enrich_cell_list)) {
        region <- as.character(temp_Neruon_enrich_cell_list[[i]])  # 字符
        overlap <- sapply(merged_regions, function(x) any(intersect(x, region) != 0))
        overlap.loc <- which(overlap == TRUE)
        if (any(overlap)) {
          if (length(overlap.loc) > 1L) {
            tmpu <- union(unlist(merged_regions[overlap.loc]), region)
            merged_regions <- c(merged_regions, list(as.character(tmpu)))
            merged_regions <- merged_regions[-c(overlap.loc)]
          } else {
            merged_regions[[overlap.loc]] <- as.character(
              union(unlist(merged_regions[[overlap.loc]]), region)
            )
          }
        } else {
          merged_regions <- c(merged_regions, list(region))
        }
      }
      
      ## =============== 步骤3：逐个 merged_region 写 merged_rgions_loc（loc 文件） ===============
      temp_Neruon_merged_regions_table <- list()
      a <- 1L
      for (i in seq_along(merged_regions)) {
        ids  <- as.character(merged_regions[[i]])
        wins <- unique(
          temp_Neruon_enrich_cell_id$loc_1[
            temp_Neruon_enrich_cell_id$cell_label %in% ids
          ]
        )
        tmp  <- temp_Neruon_enrich[temp_Neruon_enrich$loc_1 %in% wins, ]
        tmp$merge_regions <- paste("regions", i, sep = "_")
        temp_Neruon_merged_regions_table[[a]] <- tmp
        names(temp_Neruon_merged_regions_table)[a] <- unique(tmp$layer)
        a <- a + 1L
      }
      temp_Neruon_merged_regions_table <- data.table::rbindlist(
        temp_Neruon_merged_regions_table, use.names = TRUE, fill = TRUE
      )
      
      ## 文件名中的 subclass/layer 做 safe_name 替换
      out_loc_file <- file.path(
        out_root,
        paste0(chip_id, "_", safe_name(cls), "_",
               safe_name(cell_p_file_name_layer.1), "_merged_rgions_loc.txt")
      )
      loc_cols <- c("subclass","layer","merge_regions","loc",
                    "xstart","ystart","xend","yend","p","sum_value")
      loc_cols <- loc_cols[loc_cols %in% names(temp_Neruon_merged_regions_table)]
      data.table::fwrite(temp_Neruon_merged_regions_table[, ..loc_cols],
                         out_loc_file, sep = "\t", quote = FALSE)
      
      ## =============== 步骤4：根据 merged_regions 生成 cell_id 聚合表（所有 id 都是字符串） ===============
      temp_Neruon_merged_regions_table_cell_id <- list()
      a <- 1L
      ## 源 file_tmp 中的 cell_label 也确保字符
      if ("cell_label" %in% names(file_tmp)) {
        file_tmp$cell_label <- as.character(file_tmp$cell_label)
      }
      if ("cell_id" %in% names(file_tmp)) {
        file_tmp$cell_id <- as.character(file_tmp$cell_id)
      }
      
      for (i in seq_along(merged_regions)) {
        ids <- as.character(merged_regions[[i]])
        tmpx <- file_tmp[file_tmp$cell_label %in% ids, ]
        tmpx$cell_label <- as.character(tmpx$cell_label)
        
        tmpx <- data.frame(
          slide = as.character(unique(tmpx$brain_section_label)),
          layer = as.character(unique(tmpx$layer)),
          region = paste(unique(as.character(tmpx$ccf_region_name)), collapse = ","),
          merge_regions = paste("regions", i, sep = "_"),
          total_cell_num = nrow(tmpx),
          subclass = as.character(cls),
          enrich_subclass_cell_ids_num =
            length(tmpx$cell_label[tmpx$subclass == cls]),
          enrich_subclass_cell_ids =
            paste(as.character(tmpx$cell_label[tmpx$subclass == cls]),
                  collapse = ","),
          Glut_Neruon_cell_ids_num =
            length(tmpx$cell_label[tmpx$cell_Neuron_type == "Glut"]),
          Glut_Neruon_cell_ids =
            paste(as.character(tmpx$cell_label[tmpx$cell_Neuron_type == "Glut"]),
                  collapse = ","),
          GABA_Neruon_cell_ids_num =
            length(tmpx$cell_label[tmpx$cell_Neuron_type == "GABA"]),
          GABA_Neruon_cell_ids =
            paste(as.character(tmpx$cell_label[tmpx$cell_Neuron_type == "GABA"]),
                  collapse = ","),
          Non_Neruon_cell_ids_num =
            length(tmpx$cell_label[tmpx$cell_Neuron_type == "NonNeuron"]),
          Non_Neruon_cell_ids =
            paste(as.character(tmpx$cell_label[tmpx$cell_Neuron_type == "NonNeuron"]),
                  collapse = ","),
          stringsAsFactors = FALSE
        )
        temp_Neruon_merged_regions_table_cell_id[[a]] <- tmpx
        names(temp_Neruon_merged_regions_table_cell_id)[a] <- unique(tmpx$layer)
        a <- a + 1L
      }
      temp_Neruon_merged_regions_table_cell_id <- data.table::rbindlist(
        temp_Neruon_merged_regions_table_cell_id, use.names = TRUE, fill = TRUE
      )
      
      out_id_file <- file.path(
        out_root,
        paste0(chip_id, "_", safe_name(cls), "_",
               safe_name(cell_p_file_name_layer.1), "_merged_regions_table_cell_id.txt")
      )
      data.table::fwrite(temp_Neruon_merged_regions_table_cell_id,
                         out_id_file, sep = "\t", quote = FALSE)
      
      ## 只清除本层循环临时对象，避免误删主数据
      rm(list = ls()[grep(
        "^tmp$|^tmpx$|^ids$|^wins$|^temp_Neruon_merged_regions_table$|^temp_Neruon_enrich_cell_list$",
        ls(), perl = TRUE
      )], inherits = FALSE)
      gc()
    }
  }
  
  message("[INFO] Finished chip: ", chip_id)
  chip_id
}

## ========= 收尾 =========
stopCluster(cl)

message(sprintf(
  "[INFO] total_minutes: %.2f",
  as.numeric(difftime(Sys.time(), t0, units = "mins"))
))
