rm(list = ls()); gc()

setwd('E:/zaw/2511/mouseMerfish_zhuang_subclass')

suppressPackageStartupMessages({
  library(vroom)
  library(stringr)
  library(data.table)
  library(doParallel)
  library(foreach)
  library(dplyr)
  library(parallel)
})

start_time_global <- Sys.time()
cat("Start time (cell_window): ", format(start_time_global), "\n")

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

base_dir <- "E:/zaw/2511/mouseMerfish_zhuang_subclass"
combo_root <- file.path(base_dir, subdir_tag)
dir.create(combo_root, recursive = TRUE, showWarnings = FALSE)

Merfish_mouse_neocortex_layer_region <- read.delim(
  file.path(base_dir, 'Merfish_mouse_neocortex_layer_region.txt'),
  header = TRUE, stringsAsFactors = FALSE, check.names = FALSE
)
Merfish_brain_cell_type <- read.delim(
  file.path(base_dir, 'Merfish_brain_cell_type.txt'),
  header = TRUE, stringsAsFactors = FALSE, check.names = FALSE
)

in_dir  <- file.path(base_dir, 'neocortex_new')
file_names <- list.files(in_dir)

out_dir <- file.path(combo_root, "cell_window")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_cores <- max(1L, min(parallel::detectCores(logical = TRUE) - 1L, length(file_names)))

cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)
parallel::clusterEvalQ(cl, {
  data.table::setDTthreads(1L)
  NULL
})

res_log <- foreach(file_name = file_names,
                   .packages = c("data.table", "stringr"),
                   .combine = c) %dopar% {
                     fpath <- file.path(in_dir, file_name)
                     if (!file.exists(fpath)) {
                       return(NA_character_)
                     }
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
                     file_tmp <- read.delim(
                       fpath,
                       stringsAsFactors = FALSE,
                       check.names      = FALSE,
                       colClasses       = col_cls
                     )
                     if ("cell_label" %in% names(file_tmp)) {
                       file_tmp$cell_label <- as.character(file_tmp$cell_label)
                     }
                     if ("cell_id" %in% names(file_tmp)) {
                       file_tmp$cell_id <- as.character(file_tmp$cell_id)
                     }
                     file_tmp <- merge(file_tmp,
                                       Merfish_mouse_neocortex_layer_region,
                                       by = 'ccf_region_name')
                     file_tmp <- merge(file_tmp,
                                       Merfish_brain_cell_type,
                                       by = 'class')
                     if ("cell_label" %in% names(file_tmp)) {
                       file_tmp$cell_label <- as.character(file_tmp$cell_label)
                     }
                     if ("cell_id" %in% names(file_tmp)) {
                       file_tmp$cell_id <- as.character(file_tmp$cell_id)
                     }
                     file_chose_names <- sub("\\.txt$", "", file_name)
                     setDT(file_tmp)
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
                     layer_number      <- 1L
                     cell_windows_list <- list()
                     for (layer.1 in sort(unique(file_tmp$layer))) {
                       file_tmp_layer <- file_tmp[layer == layer.1, ]
                       if (nrow(file_tmp_layer) == 0) {
                         next
                       }
                       if ("cell_label" %in% names(file_tmp_layer)) {
                         file_tmp_layer[, cell_label := as.character(cell_label)]
                       }
                       min_y <- min(file_tmp_layer$y)
                       max_y <- max(file_tmp_layer$y)
                       y_values   <- 0:floor((max_y - min_y) / step_size)
                       all_result <- lapply(y_values, cell_windows_function)
                       all_result <- rbindlist(all_result, fill = TRUE)
                       if (nrow(all_result) > 0) {
                         all_result[, layer := layer.1]
                         cell_windows_list[[layer_number]] <- all_result
                         names(cell_windows_list)[layer_number] <- layer.1
                         layer_number <- layer_number + 1L
                       }
                     }
                     if (length(cell_windows_list) == 0) {
                       return(file_name)
                     }
                     cell_windows_layer_d <- rbindlist(cell_windows_list, fill = TRUE)
                     setnames(cell_windows_layer_d, 1, "cell_label")
                     cell_windows_layer_d[, cell_label := as.character(cell_label)]
                     save(file_tmp, cell_windows_layer_d,
                          file = file.path(out_dir, paste0(file_chose_names, '.RData')))
                     file_name
                   }

parallel::stopCluster(cl)

end_time_global <- Sys.time()
cat("End time (cell_window): ", format(end_time_global), "\n")
cat("Elapsed minutes (cell_window): ",
    as.numeric(difftime(end_time_global, start_time_global, units = "mins")), "\n")

rm(list = ls()); gc()

setwd('E:/zaw/2511/mouseMerfish_zhuang_subclass')

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

start_time_global <- Sys.time()
cat("Start time (cell_pvalue): ", format(start_time_global), "\n")

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
base_dir <- "E:/zaw/2511/mouseMerfish_zhuang_subclass"
combo_root <- file.path(base_dir, subdir_tag)
dir.create(combo_root, recursive = TRUE, showWarnings = FALSE)

out_cell_sample_dir <- file.path(combo_root, "cell_sample_new")
out_cell_pvalue_dir  <- file.path(combo_root, "cell_pvalue_new")
dir.create(out_cell_sample_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_cell_pvalue_dir,  recursive = TRUE, showWarnings = FALSE)

file_names <- list.files('./neocortex_new/')

cell_window_dir <- file.path(combo_root, "cell_window")

for (file_name in file_names) {
  file_chose_names <- tools::file_path_sans_ext(file_name)
  load(file.path(cell_window_dir, paste0(file_chose_names, '.RData')))
  setDT(file_tmp)
  setDT(cell_windows_layer_d)
  layer_ids <- sort(unique(file_tmp$layer))
  cell_Sliding_window_result_p <- vector("list", length(layer_ids))
  for (i in seq_along(layer_ids)) {
    lyr <- layer_ids[i]
    ft_layer  <- file_tmp[layer == lyr, .(cell_label, subclass)]
    cwd_layer <- cell_windows_layer_d[layer == lyr, .(cell_label, loc)]
    N <- uniqueN(ft_layer$cell_label)
    Kc <- ft_layer[, .(K = .N), by = subclass]
    nl <- cwd_layer[, .(n = uniqueN(cell_label)), by = loc]
    obs <- cwd_layer[ft_layer, on = .(cell_label), nomatch = 0L][
      , .(sum_value = uniqueN(cell_label)), by = .(subclass, loc)]
    if (nrow(obs) == 0L) {
      cell_Sliding_window_result_p[[i]] <- data.table(subclass=character(), loc=character(), sum_value=integer(), p=numeric(), layer=character())
      next
    }
    res_layer <- obs[Kc, on = .(subclass)][nl, on = .(loc)][
      , p := phyper(sum_value - 1L, K, N - K, n, lower.tail = FALSE)][
        , .(subclass, loc, sum_value, p)][]
    res_layer[, layer := lyr]
    cell_Sliding_window_result_p[[i]] <- res_layer
  }
  cell_Sliding_window_result_p_d <- rbindlist(cell_Sliding_window_result_p, use.names = TRUE)
  fwrite(cell_Sliding_window_result_p_d,
         file = file.path(out_cell_pvalue_dir,
                          paste0(file_chose_names, '_', 'cell_sliding_window_result_p.csv')))
  sig_high <- cell_Sliding_window_result_p_d[p < 0.05]
  fwrite(sig_high,
         file = file.path(out_cell_pvalue_dir,
                          paste0(file_chose_names, '_', 'cell_sliding_window_result_p_sig_high.csv')))
  rm(file_tmp, cell_windows_layer_d, cell_Sliding_window_result_p, cell_Sliding_window_result_p_d, sig_high)
  gc()
}

end_time_global <- Sys.time()
cat("End time (cell_pvalue): ", format(end_time_global), "\n")
cat("Elapsed minutes (cell_pvalue): ",
    as.numeric(difftime(end_time_global, start_time_global, units = "mins")), "\n")

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

setwd('E:/zaw/2511/mouseMerfish_zhuang_subclass')

start_time_global <- Sys.time()
cat("Start time (merge_region): ", format(start_time_global), "\n")

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

base_dir  <- "E:/zaw/2511/mouseMerfish_zhuang_subclass"
combo_root <- file.path(base_dir, subdir_tag)
dir.create(combo_root, recursive = TRUE, showWarnings = FALSE)

pval_root   <- file.path(combo_root, "cell_pvalue_new")
cellwin_dir <- file.path(combo_root, "cell_window")
out_root    <- file.path(combo_root, "merge_region")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

file_names <- list.files(pval_root, pattern = "sig_high\\.csv$", full.names = FALSE)

WS <- window_size

strip_dummy_col_if_needed <- function(dt){
  need <- c("subclass","loc","p","layer")
  if (!all(need %in% names(dt)) && ncol(dt) >= 2L) {
    dt2 <- dt[, -1, drop = FALSE]
    if (all(need %in% names(dt2))) return(dt2)
  }
  dt
}

safe_name <- function(x) gsub("[/\\\\]", "&", x)

workers <- max(1L, parallel::detectCores() - 1L)
cat("[INFO] Using workers:", workers, "\n")

cl <- makeCluster(workers)
doParallel::registerDoParallel(cl)

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

foreach(
  file_name = file_names,
  .packages = c("data.table","stringr"),
  .export   = c("pval_root","cellwin_dir","out_root","WS",
                "strip_dummy_col_if_needed","safe_name")
) %dopar% {
  message("[INFO] file: ", file_name)
  p_csv <- file.path(pval_root, file_name)
  mouse_sig_p.total <- read.csv(
    p_csv,
    check.names      = FALSE,
    stringsAsFactors = FALSE
  )
  mouse_sig_p.total <- strip_dummy_col_if_needed(mouse_sig_p.total)
  stopifnot(all(c("subclass","loc","p","layer") %in% names(mouse_sig_p.total)))
  
  mouse_sig_p.total$subclass <- as.character(mouse_sig_p.total$subclass)
  mouse_sig_p.total$loc      <- as.character(mouse_sig_p.total$loc)
  mouse_sig_p.total$layer    <- as.character(mouse_sig_p.total$layer)
  
  cell_p_file_name_layer <- unique(mouse_sig_p.total$layer)
  
  chip_id    <- str_split(file_name, "[_]", simplify = TRUE)[,1]
  rdata_path <- file.path(cellwin_dir, paste0(chip_id, ".RData"))
  if (!file.exists(rdata_path)) {
    message("[WARN] RData not found for chip: ", chip_id)
    return(NULL)
  }
  load(rdata_path)
  
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
  
  cell_windows_layer_d[, layer := as.character(layer)]
  cell_windows_layer_d[, loc   := as.character(loc)]
  cell_windows_layer_d[, loc_1 := paste(layer, loc, sep = "_")]
  
  for (cell_p_file_name_layer.1 in cell_p_file_name_layer) {
    cell_sliding_windows_result_p <- mouse_sig_p.total[mouse_sig_p.total$layer %in% cell_p_file_name_layer.1, ]
    if (!nrow(cell_sliding_windows_result_p)) next
    
    xy_mat <- str_split(cell_sliding_windows_result_p$loc, "[_]", simplify = TRUE)
    cell_sliding_windows_result_p$xstart <- as.numeric(xy_mat[,1])
    cell_sliding_windows_result_p$ystart <- as.numeric(xy_mat[,2])
    cell_sliding_windows_result_p$xend   <- cell_sliding_windows_result_p$xstart + WS
    cell_sliding_windows_result_p$yend   <- cell_sliding_windows_result_p$ystart + WS
    
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
      
      temp_Neruon_enrich_cell_id <- cell_windows_layer_d[
        cell_windows_layer_d$loc_1 %in% temp_Neruon_enrich$loc_1, ]
      if (!nrow(temp_Neruon_enrich_cell_id)) next
      
      temp_Neruon_enrich_cell_id$cell_label <- as.character(temp_Neruon_enrich_cell_id$cell_label)
      temp_Neruon_enrich_cell_list <- list()
      a <- 1L
      for (i in unique(temp_Neruon_enrich_cell_id$loc_1)) {
        temp1 <- temp_Neruon_enrich_cell_id[temp_Neruon_enrich_cell_id$loc_1 %in% i, ]
        temp1$cell_label <- as.character(temp1$cell_label)
        temp_Neruon_enrich_cell_list[[a]] <- temp1$cell_label
        names(temp_Neruon_enrich_cell_list)[a] <- unique(paste(temp1$layer, temp1$loc, sep = "_"))
        a <- a + 1L
      }
      
      merged_regions <- list()
      for (i in seq_along(temp_Neruon_enrich_cell_list)) {
        region <- as.character(temp_Neruon_enrich_cell_list[[i]])
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
      
      temp_Neruon_merged_regions_table_cell_id <- list()
      a <- 1L
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
stopCluster(cl)

end_time_global <- Sys.time()
cat("End time (merge_region): ", format(end_time_global), "\n")
cat("Elapsed minutes (merge_region): ",
    as.numeric(difftime(end_time_global, start_time_global, units = "mins")), "\n")


rm(list = ls()); gc()

library(vroom)
library(stringr)
library(data.table)
library(doParallel)
library(foreach)
library(dplyr)
library(parallel)

base_dir <- "E:/zaw/2511/mouseMerfish_zhuang_subclass"

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

merge_region_dir <- file.path(combo_root, "merge_region")
mapping_path     <- file.path(combo_root, "Merfish_brain_cell_type_subclass.txt")
setwd(base_dir)

file_names <- list.files(
  merge_region_dir,
  pattern = "table_cell_id\\.txt$",
  full.names = FALSE
)

workers <- max(1L, min(length(file_names), parallel::detectCores() - 1L))
cl <- makeCluster(workers)
registerDoParallel(cl)

clusterEvalQ(cl, {
  library(data.table)
  data.table::setDTthreads(1L)
  NULL
})
data.table::setDTthreads(1L)

mouse_list <- foreach(
  fn = file_names,
  .packages = c("data.table"),
  .combine  = "c"
) %dopar% {
  fpath <- file.path(merge_region_dir, fn)
  dt <- read.delim(
    fpath,
    check.names      = FALSE,
    stringsAsFactors = FALSE
  )
  data.table::setDT(dt)
  chr_cols <- intersect(
    c("slide", "layer", "region", "merge_regions", "subclass"),
    names(dt)
  )
  for (cc in chr_cols) {
    dt[[cc]] <- as.character(dt[[cc]])
  }
  list(dt)
}

stopCluster(cl)

mouse_subclass_cluster_total <- rbindlist(
  mouse_list,
  use.names = TRUE,
  fill      = TRUE
)
setDT(mouse_subclass_cluster_total)

if ("subclass" %in% names(mouse_subclass_cluster_total)) {
  mouse_subclass_cluster_total[, subclass := as.character(subclass)]
}

sub_vec <- unique(mouse_subclass_cluster_total$subclass)
df <- data.frame(subclass = sub_vec, stringsAsFactors = FALSE)
last_token <- stringr::word(df$subclass, -1)

df$cell_Neuron_type <- dplyr::case_when(
  stringr::str_detect(last_token, "Glut$") ~ "Glut",
  stringr::str_detect(last_token, "Gaba$") ~ "Gaba",
  TRUE                                     ~ "NonNeuron"
)

tmp <- stringr::str_remove(df$subclass, "^\\S+\\s+")
df$subclass_sim <- tmp

idx_glut <- df$cell_Neuron_type == "Glut"
idx_gaba <- df$cell_Neuron_type == "Gaba"

df$subclass_sim[idx_glut] <- stringr::str_replace(
  df$subclass_sim[idx_glut],
  "([ -])Glut$",
  ""
)

df$subclass_sim[idx_gaba] <- stringr::str_replace(
  df$subclass_sim[idx_gaba],
  "([ -])Gaba$",
  ""
)

out_df <- df[, c("subclass", "cell_Neuron_type", "subclass_sim")]

write.table(
  out_df,
  file      = mapping_path,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

Merfish_brain_cell_type <- read.delim(
  mapping_path,
  check.names      = FALSE,
  stringsAsFactors = FALSE
)
setDT(Merfish_brain_cell_type)
setnames(
  Merfish_brain_cell_type,
  old = names(Merfish_brain_cell_type),
  new = c("subclass", "cell_Neruon_type", "subclass_sim")
)
Merfish_brain_cell_type[, subclass         := as.character(subclass)]
Merfish_brain_cell_type[, cell_Neruon_type := as.character(cell_Neruon_type)]
Merfish_brain_cell_type[, subclass_sim     := as.character(subclass_sim)]

mouse_subclass_cluster_total <- merge(
  mouse_subclass_cluster_total,
  Merfish_brain_cell_type,
  by    = "subclass",
  all.x = TRUE, all.y = FALSE
)

front_cols <- c(
  "slide", "layer", "region", "merge_regions",
  "cell_Neruon_type", "subclass", "subclass_sim"
)
front_cols <- intersect(front_cols, names(mouse_subclass_cluster_total))
other_cols <- setdiff(names(mouse_subclass_cluster_total), front_cols)
mouse_subclass_cluster_total <- mouse_subclass_cluster_total[
  , .SD, .SDcols = c(front_cols, other_cols)
]

write.table(
  mouse_subclass_cluster_total,
  file      = file.path(combo_root, "mouse_subclass_cluster_total.txt"),
  quote     = FALSE,
  sep       = "\t",
  row.names = FALSE
)

rm(list = ls()); gc()

library(vroom)
library(stringr)
library(data.table)
library(doParallel)
library(foreach)
library(dplyr)
library(parallel)

Sys.setenv(
  OPENBLAS_NUM_THREADS = "1",
  OMP_NUM_THREADS      = "1",
  MKL_NUM_THREADS      = "1",
  GOTO_NUM_THREADS     = "1"
)
data.table::setDTthreads(1L)

setwd('E:/zaw/2511/mouseMerfish_zhuang_subclass')

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
base_dir <- "E:/zaw/2511/mouseMerfish_zhuang_subclass"
combo_root <- file.path(base_dir, subdir_tag)
dir.create(combo_root, recursive = TRUE, showWarnings = FALSE)

file_names <- list.files(file.path(base_dir, 'neocortex_new'))

merged_regions_total <- list.files(file.path(combo_root, 'merge_region'))
merged_regions_total <- merged_regions_total[grep('cell_id.txt', merged_regions_total)]

workers <- max(1L, min(length(file_names), detectCores() - 1L))
cl1 <- makeCluster(workers)
registerDoParallel(cl1)

foreach(i = file_names,
        .combine   = c,
        .packages  = c("stringr")) %dopar% {
          file_name <- str_split(i, "\\.txt", simplify = TRUE)[,1]
          merged_regions <- merged_regions_total[grep(file_name, merged_regions_total)]
          temp1.cell.id.total <- data.frame(stringsAsFactors = FALSE)
          if (length(merged_regions) > 0) {
            for (x in merged_regions) {
              temp1.cell.id <- read.delim(
                file.path(combo_root, 'merge_region', x),
                sep             = '\t',
                header          = TRUE,
                stringsAsFactors = FALSE,
                check.names      = FALSE
              )
              chr_cols <- intersect(
                c("slide","layer","region","merge_regions","subclass",
                  "Glut_Neruon_cell_ids","GABA_Neruon_cell_ids","Non_Neruon_cell_ids"),
                names(temp1.cell.id)
              )
              for (cc in chr_cols) {
                temp1.cell.id[[cc]] <- as.character(temp1.cell.id[[cc]])
              }
              temp1.cell.id.total <- rbind(temp1.cell.id.total, temp1.cell.id)
            }
          }
          merge_region_cell_id <- data.frame(stringsAsFactors = FALSE)
          if (nrow(temp1.cell.id.total) > 0) {
            for (x in 1:nrow(temp1.cell.id.total)) {
              temp2 <- temp1.cell.id.total[x,]
              temp2.cell.id.GLU <- unique(as.character(
                str_split(temp2$Glut_Neruon_cell_ids, "[,]", simplify = TRUE)
              ))
              temp2.cell.id.GABA <- unique(as.character(
                str_split(temp2$GABA_Neruon_cell_ids,  "[,]", simplify = TRUE)
              ))
              temp2.cell.id.NonNeuron <- unique(as.character(
                str_split(temp2$Non_Neruon_cell_ids,   "[,]", simplify = TRUE)
              ))
              temp2.cell.id <- Reduce(union,
                                      list(temp2.cell.id.GLU,
                                           temp2.cell.id.GABA,
                                           temp2.cell.id.NonNeuron))
              temp3 <- data.frame(cell_id = as.character(temp2.cell.id),
                                  stringsAsFactors = FALSE)
              temp3 <- cbind(temp2[, c(1:4, 6, 5)], temp3)
              merge_region_cell_id <- rbind(merge_region_cell_id, temp3)
            }
          }
          out_dir_merge_cell <- file.path(combo_root, 'nearby_merge_regions_new', 'merge_cell_id')
          dir.create(out_dir_merge_cell, recursive = TRUE, showWarnings = FALSE)
          write.table(
            merge_region_cell_id,
            file = file.path(out_dir_merge_cell,
                             paste0(file_name, '.merge_region_cell_id.txt')),
            sep       = '\t',
            quote     = FALSE,
            row.names = FALSE
          )
          file_name
        }

stopCluster(cl1)

rm(list = ls()); gc()

library(vroom)
library(stringr)
library(data.table)
library(doParallel)
library(foreach)
library(dplyr)
library(parallel)

Sys.setenv(
  OPENBLAS_NUM_THREADS = "1",
  OMP_NUM_THREADS      = "1",
  MKL_NUM_THREADS      = "1",
  GOTO_NUM_THREADS     = "1"
)
data.table::setDTthreads(1L)

setwd('E:/zaw/2511/mouseMerfish_zhuang_subclass')

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
base_dir <- "E:/zaw/2511/mouseMerfish_zhuang_subclass"
combo_root <- file.path(base_dir, subdir_tag)
dir.create(combo_root, recursive = TRUE, showWarnings = FALSE)

file_names <- list.files(file.path(base_dir, 'neocortex_new'))

workers2 <- max(1L, min(length(file_names), detectCores() - 1L))
cl2 <- makeCluster(workers2)
registerDoParallel(cl2)

foreach(i = file_names,
        .combine  = c,
        .packages = c("stringr","data.table")) %dopar% {
          file_name <- str_split(i, "\\.txt", simplify = TRUE)[, 1]
          fpath_tmp <- file.path(base_dir, 'neocortex_new', i)
          header_dt <- read.delim(
            fpath_tmp,
            nrows           = 0,
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
          file_tmp <- read.delim(
            fpath_tmp,
            stringsAsFactors = FALSE,
            check.names      = FALSE,
            colClasses       = col_cls
          )
          if ("cell_label" %in% names(file_tmp)) {
            file_tmp$cell_label <- as.character(file_tmp$cell_label)
          }
          if ("cell_id" %in% names(file_tmp)) {
            file_tmp$cell_id <- as.character(file_tmp$cell_id)
          }
          fpath_merge <- file.path(combo_root, 'nearby_merge_regions_new', 'merge_cell_id',
                                   paste0(file_name, '.merge_region_cell_id.txt'))
          header_m <- read.delim(
            fpath_merge,
            nrows           = 0,
            stringsAsFactors = FALSE,
            check.names      = FALSE
          )
          col_nms_m <- names(header_m)
          col_cls_m <- rep(NA_character_, length(col_nms_m))
          names(col_cls_m) <- col_nms_m
          if ("cell_id" %in% col_nms_m) {
            col_cls_m["cell_id"] <- "character"
          }
          if ("layer" %in% col_nms_m) {
            col_cls_m["layer"] <- "character"
          }
          if ("subclass" %in% col_nms_m) {
            col_cls_m["subclass"] <- "character"
          }
          if ("merge_regions" %in% col_nms_m) {
            col_cls_m["merge_regions"] <- "character"
          }
          temp1 <- read.delim(
            fpath_merge,
            sep              = '\t',
            header           = TRUE,
            stringsAsFactors = FALSE,
            check.names      = FALSE,
            colClasses       = col_cls_m
          )
          if ("cell_id" %in% names(temp1)) {
            temp1$cell_id <- as.character(temp1$cell_id)
          }
          if ("layer" %in% names(temp1)) {
            temp1$layer <- as.character(temp1$layer)
          }
          if ("subclass" %in% names(temp1)) {
            temp1$subclass <- as.character(temp1$subclass)
          }
          if ("merge_regions" %in% names(temp1)) {
            temp1$merge_regions <- as.character(temp1$merge_regions)
          }
          if (nrow(temp1) == 0) {
            merge_regions_nearby_dist     <- data.frame()
            merge_regions_nearby_dist_sim <- data.frame()
          } else {
            temp1$loc <- paste(temp1$layer, temp1$subclass, temp1$merge_regions, sep = '_')
            merge_regions_nearby_dist <- list()
            a <- 1L
            for (layer.1 in unique(temp1$layer)) {
              temp1.layer <- temp1[temp1$layer %in% layer.1, , drop = FALSE]
              loc_vec <- unique(temp1.layer$loc)
              for (x_loc in loc_vec) {
                temp2 <- temp1.layer[temp1.layer$loc %in% x_loc, , drop = FALSE]
                temp2.cell <- file_tmp[file_tmp$cell_label %in% temp2$cell_id, , drop = FALSE]
                if (nrow(temp2.cell) == 0) next
                colnames(temp2.cell)[2] <- 'cell_id'
                temp2.cell <- merge(temp2,
                                    temp2.cell[, c('cell_id','x','y')],
                                    by = 'cell_id')
                if (nrow(temp2.cell) == 0) next
                temp2.cell$name <- paste(x_loc, rownames(temp2.cell), sep = '_')
                rownames(temp2.cell) <- temp2.cell$name
                for (y_loc in loc_vec) {
                  temp3 <- temp1.layer[temp1.layer$loc %in% y_loc, , drop = FALSE]
                  temp3.cell <- file_tmp[file_tmp$cell_label %in% temp3$cell_id, , drop = FALSE]
                  if (nrow(temp3.cell) == 0) next
                  colnames(temp3.cell)[2] <- 'cell_id'
                  temp3.cell <- merge(temp3,
                                      temp3.cell[, c('cell_id','x','y')],
                                      by = 'cell_id')
                  if (nrow(temp3.cell) == 0) next
                  temp3.cell$name <- paste(y_loc, rownames(temp3.cell), sep = '_')
                  rownames(temp3.cell) <- temp3.cell$name
                  temp.overlap.cell <- intersect(temp2.cell$cell_id, temp3.cell$cell_id)
                  if (length(temp.overlap.cell) != 0) {
                    temp4.min.row <- temp2.cell[temp2.cell$cell_id %in% temp.overlap.cell, , drop = FALSE]
                    temp4.min.col <- temp3.cell[temp3.cell$cell_id %in% temp.overlap.cell, , drop = FALSE]
                    temp5 <- cbind(temp4.min.row, temp4.min.col)
                    temp5$nearby_dist <- 0
                  } else {
                    num_temp2 <- nrow(temp2.cell)
                    num_temp3 <- nrow(temp3.cell)
                    if (num_temp2 == 0 || num_temp3 == 0) next
                    temp4.dist <- as.matrix(
                      dist(rbind(temp2.cell[, c("x","y")],
                                 temp3.cell[, c("x","y")]),
                           method = 'euclidean')
                    )
                    temp4.distances <- temp4.dist[1:num_temp2,
                                                  (num_temp2 + 1):(num_temp2 + num_temp3),
                                                  drop = FALSE]
                    temp4.nearest_indices <- apply(temp4.distances, 1, which.min)
                    temp4.min_distances  <- apply(temp4.distances, 1, min)
                    temp4.dist.total <- data.frame(
                      temp2_cell_id         = temp2.cell$cell_id,
                      nearest_temp3_cell_id = temp3.cell$cell_id[temp4.nearest_indices],
                      min_distance          = temp4.min_distances,
                      stringsAsFactors      = FALSE
                    )
                    which_min <- which.min(temp4.dist.total$min_distance)
                    temp4.min.row <- temp2.cell[temp2.cell$cell_id %in%
                                                  temp4.dist.total$temp2_cell_id[which_min],
                                                , drop = FALSE]
                    temp4.min.col <- temp3.cell[temp3.cell$cell_id %in%
                                                  temp4.dist.total$nearest_temp3_cell_id[which_min],
                                                , drop = FALSE]
                    temp5 <- cbind(temp4.min.row, temp4.min.col)
                    temp5$nearby_dist <- temp4.dist.total$min_distance[which_min]
                  }
                  merge_regions_nearby_dist[[a]] <- temp5
                  a <- a + 1L
                }
              }
            }
            if (length(merge_regions_nearby_dist) > 0) {
              merge_regions_nearby_dist <- data.table::rbindlist(merge_regions_nearby_dist, fill = TRUE)
              merge_regions_nearby_dist <- merge_regions_nearby_dist[, c(2,3,8,1,19,12,23)]
              colnames(merge_regions_nearby_dist)[3:6] <- c('cluster.1','cluster.1.cell_id','cluster.2','cluster.2.cell_id')
              merge_regions_nearby_dist <- merge_regions_nearby_dist[
                merge_regions_nearby_dist$cluster.1 != merge_regions_nearby_dist$cluster.2, ]
              merge_regions_nearby_dist_sim <- merge_regions_nearby_dist[, c(1:3,5,7)]
              merge_regions_nearby_dist_sim <- merge_regions_nearby_dist_sim[!duplicated(merge_regions_nearby_dist_sim), ]
            } else {
              merge_regions_nearby_dist     <- data.frame()
              merge_regions_nearby_dist_sim <- data.frame()
            }
          }
          out_dir_near_dist <- file.path(combo_root, 'nearby_merge_regions_new', 'merge_regions_nearby_dist')
          dir.create(out_dir_near_dist, recursive = TRUE, showWarnings = FALSE)
          write.table(
            merge_regions_nearby_dist,
            file = file.path(out_dir_near_dist,
                             paste0(file_name, '.merge_regions_nearby_dist.txt')),
            sep       = '\t',
            quote     = FALSE,
            row.names = FALSE
          )
          write.table(
            merge_regions_nearby_dist_sim,
            file = file.path(out_dir_near_dist,
                             paste0(file_name, '.merge_regions_nearby_dist_sim.txt')),
            sep       = '\t',
            quote     = FALSE,
            row.names = FALSE
          )
          file_name
        }

stopCluster(cl2)

rm(list = ls()); gc()
setwd('E:/zaw/2511/mouseMerfish_zhuang_subclass')

suppressPackageStartupMessages({
  library(stringr)
  library(data.table)
  library(doParallel)
  library(foreach)
  library(parallel)
})

Sys.setenv(
  OPENBLAS_NUM_THREADS = "1",
  OMP_NUM_THREADS      = "1",
  MKL_NUM_THREADS      = "1",
  GOTO_NUM_THREADS     = "1"
)
data.table::setDTthreads(1L)

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
base_dir <- "E:/zaw/2511/mouseMerfish_zhuang_subclass"
combo_root <- file.path(base_dir, subdir_tag)
dir.create(combo_root, recursive = TRUE, showWarnings = FALSE)

near_dist_dir <- file.path(combo_root, 'nearby_merge_regions_new', 'merge_regions_nearby_dist')

file_names <- list.files(near_dist_dir)
file_names <- file_names[grep('dist.txt$', file_names)]

workers <- max(1L, min(length(file_names), detectCores() - 1L))
cl1 <- makeCluster(workers)
registerDoParallel(cl1)

dt_rbind <- function(...) data.table::rbindlist(list(...), use.names = TRUE, fill = TRUE)

error_list_all <- foreach(
  fn = file_names,
  .combine      = dt_rbind,
  .multicombine = TRUE,
  .packages     = c("stringr","data.table")
) %dopar% {
  dist_path <- file.path(near_dist_dir, fn)
  if (!file.exists(dist_path) || file.info(dist_path)$size == 0L) {
    return(data.table())
  }
  
  file_temp <- read.delim(
    dist_path,
    sep              = '\t',
    header           = TRUE,
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )
  if (!nrow(file_temp)) return(data.table())
  
  id_cols1 <- grep("cell_label|cell_id", names(file_temp),
                   ignore.case = TRUE, value = TRUE)
  for (cc in id_cols1) {
    file_temp[[cc]] <- as.character(file_temp[[cc]])
  }
  
  slice_id  <- str_split(fn, "\\.merge", simplify = TRUE)[, 1]
  cell_path <- file.path(base_dir, "neocortex_new", paste0(slice_id, ".txt"))
  
  if (!file.exists(cell_path) || file.info(cell_path)$size == 0L) {
    return(data.table())
  }
  
  header_dt <- read.delim(
    cell_path,
    nrows            = 0,
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )
  col_nms <- names(header_dt)
  col_cls <- rep(NA_character_, length(col_nms))
  names(col_cls) <- col_nms
  if ("cell_label" %in% col_nms) col_cls["cell_label"] <- "character"
  if ("cell_id"    %in% col_nms) col_cls["cell_id"]    <- "character"
  
  file_temp_cell <- read.delim(
    cell_path,
    sep              = '\t',
    header           = TRUE,
    stringsAsFactors = FALSE,
    check.names      = FALSE,
    colClasses       = col_cls
  )
  if ("cell_label" %in% names(file_temp_cell)) {
    file_temp_cell$cell_label <- as.character(file_temp_cell$cell_label)
  }
  if ("cell_id" %in% names(file_temp_cell)) {
    file_temp_cell$cell_id <- as.character(file_temp_cell$cell_id)
  }
  
  n_row <- nrow(file_temp)
  if (n_row == 0L) return(data.table())
  
  iter <- min(1000L, n_row)
  err_list <- vector("list", iter)
  k <- 0L
  
  for (i in seq_len(iter)) {
    temp1 <- file_temp[sample.int(n_row, 1L), , drop = FALSE]
    
    temp1.x <- file_temp_cell[file_temp_cell$cell_label %in% temp1$cluster.1.cell_id, ]
    temp1.y <- file_temp_cell[file_temp_cell$cell_label %in% temp1$cluster.2.cell_id, ]
    
    if (nrow(temp1.x) == 0L || nrow(temp1.y) == 0L) next
    
    temp_xy <- rbind(
      temp1.x[, c("x","y")],
      temp1.y[, c("x","y")]
    )
    temp1.dist <- dist(as.matrix(temp_xy), method = "euclidean")
    
    a <- all.equal(as.numeric(temp1.dist[1]), as.numeric(temp1$nearby_dist))
    if (!isTRUE(a)) {
      k <- k + 1L
      err_row <- as.data.table(temp1)
      err_row[, src_file := fn]
      err_list[[k]] <- err_row
    }
  }
  
  if (k > 0L) {
    data.table::rbindlist(err_list[1:k], use.names = TRUE, fill = TRUE)
  } else {
    data.table()
  }
}


stopCluster(cl1)

error_list <- error_list_all

rm(list = ls()); gc()
setwd('E:/zaw/2511/mouseMerfish_zhuang_subclass')

suppressPackageStartupMessages({
  library(stringr)
  library(data.table)
  library(doParallel)
  library(foreach)
  library(parallel)
})

Sys.setenv(
  OPENBLAS_NUM_THREADS = "1",
  OMP_NUM_THREADS      = "1",
  MKL_NUM_THREADS      = "1",
  GOTO_NUM_THREADS     = "1"
)
data.table::setDTthreads(1L)

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
base_dir <- "E:/zaw/2511/mouseMerfish_zhuang_subclass"
combo_root <- file.path(base_dir, subdir_tag)
dir.create(combo_root, recursive = TRUE, showWarnings = FALSE)

near_dist_dir <- file.path(combo_root, 'nearby_merge_regions_new', 'merge_regions_nearby_dist')

file_names <- list.files(near_dist_dir)
file_names <- file_names[grep('sim.txt$', file_names)]

workers2 <- max(1L, min(length(file_names), detectCores() - 1L))
cl2 <- makeCluster(workers2)
registerDoParallel(cl2)

dt_rbind <- function(...) data.table::rbindlist(list(...), use.names = TRUE, fill = TRUE)

mouse_cluster_pairship_list <- foreach(fn = file_names,
                                       .combine      = dt_rbind,
                                       .multicombine = TRUE,
                                       .packages     = c("stringr","data.table")) %dopar% {
                                         file_temp <- read.delim(
                                           file.path(near_dist_dir, fn),
                                           sep              = '\t',
                                           header           = TRUE,
                                           stringsAsFactors = FALSE,
                                           check.names      = FALSE
                                         )
                                         if (!nrow(file_temp)) return(data.table())
                                         id_cols2 <- grep("cell_label|cell_id", names(file_temp), ignore.case = TRUE, value = TRUE)
                                         for (cc in id_cols2) {
                                           file_temp[[cc]] <- as.character(file_temp[[cc]])
                                         }
                                         as.data.table(file_temp)
                                       }

stopCluster(cl2)

mouse_cluster_pairship <- mouse_cluster_pairship_list
rm(mouse_cluster_pairship_list); gc()

n <- nrow(mouse_cluster_pairship)
mouse_cluster_pairship_new <- vector("list", n)
for (i in seq_len(n)) {
  temp1 <- mouse_cluster_pairship[i, ]
  temp1.label <- union(
    as.character(temp1$cluster.1),
    as.character(temp1$cluster.2)
  )
  temp1.label <- temp1.label[order(temp1.label)]
  temp1$label <- paste(temp1.label[1], temp1.label[2], sep = '-')
  mouse_cluster_pairship_new[[i]] <- temp1
}
mouse_cluster_pairship_new <- data.table::rbindlist(mouse_cluster_pairship_new, use.names = TRUE, fill = TRUE)

write.table(
  mouse_cluster_pairship_new,
  file = file.path(combo_root, 'mouse_cluster_pairship_new_total.txt'),
  quote = FALSE, row.names = FALSE, sep = '\t'
)

mouse_cluster_pairship_new.1 <- mouse_cluster_pairship_new[
  !duplicated(mouse_cluster_pairship_new[, c(1,2,5,6)])
]

write.table(
  mouse_cluster_pairship_new.1,
  file = file.path(combo_root, 'mouse_cluster_pairship_new.txt'),
  quote = FALSE, row.names = FALSE, sep = '\t'
)

rm(list = ls()); gc()

suppressPackageStartupMessages({
  library(stringr)
  library(data.table)
  library(doParallel)
  library(foreach)
  library(parallel)
})

Sys.setenv(
  OPENBLAS_NUM_THREADS = "1",
  OMP_NUM_THREADS      = "1",
  MKL_NUM_THREADS      = "1",
  GOTO_NUM_THREADS     = "1"
)
data.table::setDTthreads(1L)

setwd('E:/zaw/2511/mouseMerfish_zhuang_subclass')

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
base_dir <- "E:/zaw/2511/mouseMerfish_zhuang_subclass"
combo_root <- file.path(base_dir, subdir_tag)
dir.create(combo_root, recursive = TRUE, showWarnings = FALSE)

inp_cluster_file <- file.path(combo_root, "mouse_subclass_cluster_total.txt")
loc_dir         <- file.path(combo_root, "merge_region")

force_id_cols_char <- function(df) {
  id_cols <- names(df)[grepl("cell_label$|cell_id$|cell_ids$", names(df))]
  for (cc in id_cols) {
    df[[cc]] <- as.character(df[[cc]])
  }
  df
}

stopifnot(file.exists(inp_cluster_file))
mouse_subclass_cluster_total <- fread(
  inp_cluster_file,
  sep         = "\t",
  header      = TRUE,
  data.table  = FALSE,
  check.names = FALSE
)
mouse_subclass_cluster_total <- force_id_cols_char(mouse_subclass_cluster_total)

num_cols <- c("total_cell_num",
              "Glut_Neruon_cell_ids_num",
              "GABA_Neruon_cell_ids_num",
              "enrich_subclass_cell_ids_num")
for (nm in num_cols) {
  if (nm %in% names(mouse_subclass_cluster_total)) {
    mouse_subclass_cluster_total[[nm]] <- as.numeric(mouse_subclass_cluster_total[[nm]])
  }
}

need_cols_cluster <- c("slide","layer","region","cell_Neruon_type","subclass",
                       "total_cell_num","Glut_Neruon_cell_ids_num","GABA_Neruon_cell_ids_num",
                       "Glut_Neruon_cell_ids","GABA_Neruon_cell_ids","Non_Neruon_cell_ids",
                       "merge_regions","enrich_subclass_cell_ids_num")
miss_cluster <- setdiff(need_cols_cluster, names(mouse_subclass_cluster_total))
if (length(miss_cluster)) {
  stop(sprintf("missing columns in mouse_subclass_cluster_total: %s",
               paste(miss_cluster, collapse = ",")))
}

mouse_subclass_cluster_total <- mouse_subclass_cluster_total[
  mouse_subclass_cluster_total$enrich_subclass_cell_ids_num >= 3,
  ,
  drop = FALSE
]

fwrite(
  mouse_subclass_cluster_total,
  file = file.path(combo_root, "mouse_subclass_cluster_total_over3.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cauchy_combination_test <- function(pvals, weights = NULL) {
  pvals <- as.numeric(pvals)
  pvals[pvals < 1e-15]     <- 1e-15
  pvals[pvals > 1 - 1e-15] <- 1 - 1e-15
  if (is.null(weights)) weights <- rep(1 / length(pvals), length(pvals))
  t_stat <- sum(weights * tan((0.5 - pvals) * pi))
  0.5 - atan(t_stat) / pi
}

stopifnot(dir.exists(loc_dir))
file_names <- list.files(loc_dir,
                         pattern = "merged_rgions_loc\\.txt$",
                         full.names = FALSE)

temp.total <- data.table()

if (length(file_names) > 0) {
  workers <- max(1L, min(length(file_names), detectCores() - 1L))
  cl <- makeCluster(workers)
  registerDoParallel(cl)
  temp.list <- foreach(file_name = file_names,
                       .packages = c("data.table", "stringr"),
                       .export   = c("cauchy_combination_test")) %dopar% {
                         fp <- file.path(loc_dir, file_name)
                         dt <- tryCatch(
                           fread(fp, sep = "\t", header = TRUE,
                                 data.table = FALSE, check.names = FALSE),
                           error = function(e) NULL
                         )
                         if (is.null(dt) || !nrow(dt)) {
                           return(data.frame(
                             slide = character(0),
                             layer = character(0),
                             subclass = character(0),
                             merge_regions = character(0),
                             cauchy_combination_p = numeric(0),
                             stringsAsFactors = FALSE
                           ))
                         }
                         req_cols <- c("merge_regions","layer","subclass","p")
                         if (!all(req_cols %in% names(dt))) {
                           return(data.frame(
                             slide = character(0),
                             layer = character(0),
                             subclass = character(0),
                             merge_regions = character(0),
                             cauchy_combination_p = numeric(0),
                             stringsAsFactors = FALSE
                           ))
                         }
                         slide_id <- str_split(file_name, "[_]", simplify = TRUE)[,1]
                         lyr      <- unique(dt$layer)
                         cls      <- unique(dt$subclass)
                         if (length(lyr) != 1L) lyr <- lyr[1]
                         if (length(cls) != 1L) cls <- cls[1]
                         regs <- sort(unique(dt$merge_regions))
                         if (!length(regs)) {
                           return(data.frame(
                             slide = character(0),
                             layer = character(0),
                             subclass = character(0),
                             merge_regions = character(0),
                             cauchy_combination_p = numeric(0),
                             stringsAsFactors = FALSE
                           ))
                         }
                         out_list <- vector("list", length(regs))
                         k <- 1L
                         for (rg in regs) {
                           sub <- dt[dt$merge_regions %in% rg, , drop = FALSE]
                           if (!nrow(sub)) next
                           pc <- cauchy_combination_test(sub$p)
                           out_list[[k]] <- data.frame(
                             slide = slide_id,
                             layer = lyr,
                             subclass = cls,
                             merge_regions = rg,
                             cauchy_combination_p = pc,
                             stringsAsFactors = FALSE
                           )
                           k <- k + 1L
                         }
                         if (k == 1L) {
                           data.frame(
                             slide = character(0),
                             layer = character(0),
                             subclass = character(0),
                             merge_regions = character(0),
                             cauchy_combination_p = numeric(0),
                             stringsAsFactors = FALSE
                           )
                         } else {
                           do.call(rbind, out_list[1:(k-1)])
                         }
                       }
  stopCluster(cl)
  temp.total <- if (length(temp.list)) {
    rbindlist(temp.list, use.names = TRUE, fill = TRUE)
  } else data.table()
}

if (nrow(temp.total)) {
  temp.total$label <- paste(
    temp.total$slide,
    temp.total$layer,
    temp.total$subclass,
    temp.total$merge_regions, sep = "_"
  )
} else {
  temp.total$label <- character(0)
}

mouse_subclass_cluster_total$label <- paste(
  mouse_subclass_cluster_total$slide,
  mouse_subclass_cluster_total$layer,
  mouse_subclass_cluster_total$subclass,
  mouse_subclass_cluster_total$merge_regions, sep = "_"
)

keep_cols <- c("label","cauchy_combination_p")
mouse_subclass_cluster_total_cauchy_combination_test <-
  merge(mouse_subclass_cluster_total,
        temp.total[, ..keep_cols],
        by = "label", all.x = FALSE, all.y = FALSE)

dir.create(file.path(combo_root, "mouse_table"), showWarnings = FALSE)
dir.create(file.path(combo_root, "mouse_table", "table1"), showWarnings = FALSE)

tab1_cols <- c("label","slide","layer","region","cell_Neruon_type","subclass",
               "total_cell_num","Glut_Neruon_cell_ids_num",
               "GABA_Neruon_cell_ids_num","cauchy_combination_p")
Table_1_total_cell <- mouse_subclass_cluster_total_cauchy_combination_test[, tab1_cols, drop = FALSE]

Table_1_total_cell$E_I_Ratio <- with(
  Table_1_total_cell,
  Glut_Neruon_cell_ids_num / ifelse(GABA_Neruon_cell_ids_num == 0, NA, GABA_Neruon_cell_ids_num)
)

fwrite(Table_1_total_cell,
       file.path(combo_root, "mouse_table", "table1", "Table_1_total_cell.txt"),
       sep = "\t", quote = FALSE, row.names = FALSE)

mouse_cluster_number <- as.data.frame(
  table(mouse_subclass_cluster_total$subclass),
  stringsAsFactors = FALSE
)
colnames(mouse_cluster_number) <- c("subclass","cluster_number")
mouse_cluster_number <- mouse_cluster_number[
  order(mouse_cluster_number$cluster_number, decreasing = TRUE), ]
top_subclasses <- head(as.character(mouse_cluster_number$subclass), 8)

library(data.table)
library(stringr)

table1_dir <- file.path(combo_root, "mouse_table", "table1")
dir.create(table1_dir, recursive = TRUE, showWarnings = FALSE)

safe_subclass_name <- function(x) {
  x <- as.character(x)
  x <- gsub("[/\\\\]", "-", x)
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("[^[:alnum:]_\\-]+", "", x)
  x
}

for (cl in top_subclasses) {
  
  sub1 <- Table_1_total_cell[Table_1_total_cell$subclass %in% cl, , drop = FALSE]
  
  if (nrow(sub1) == 0) {
    message("Warning: no rows for subclass: ", cl)
    next
  }

  cl_safe <- safe_subclass_name(cl)
  
  outp <- file.path(
    table1_dir,
    sprintf("Table_1_%s_cell.txt", cl_safe)
  )
  
  fwrite(sub1, outp, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Written: ", outp)
}


pairship_file <- file.path(combo_root, "mouse_cluster_pairship_new_total.txt")
stopifnot(file.exists(pairship_file))
mouse_cluster_pairship_new_total <- fread(
  pairship_file,
  sep         = "\t",
  header      = TRUE,
  data.table  = FALSE,
  check.names = FALSE
)
mouse_cluster_pairship_new_total <- force_id_cols_char(mouse_cluster_pairship_new_total)

mouse_cluster_pairship_new_total <- mouse_cluster_pairship_new_total[
  mouse_cluster_pairship_new_total$nearby_dist == 0, , drop = FALSE
]

mouse_cluster_pairship_new_total$label.1 <- paste(
  mouse_cluster_pairship_new_total$slide,
  mouse_cluster_pairship_new_total$cluster.1, sep = "_"
)
mouse_cluster_pairship_new_total$label.2 <- paste(
  mouse_cluster_pairship_new_total$slide,
  mouse_cluster_pairship_new_total$cluster.2, sep = "_"
)
mouse_cluster_pairship_new_total <- mouse_cluster_pairship_new_total[, c("label.1","label.2"), drop = FALSE]

mouse_subclass_cluster_total_cauchy_combination_test$label <- paste(
  mouse_subclass_cluster_total_cauchy_combination_test$slide,
  mouse_subclass_cluster_total_cauchy_combination_test$layer,
  mouse_subclass_cluster_total_cauchy_combination_test$subclass,
  mouse_subclass_cluster_total_cauchy_combination_test$merge_regions, sep = "_"
)
mouse_subclass_cluster_total$label <- paste(
  mouse_subclass_cluster_total$slide,
  mouse_subclass_cluster_total$layer,
  mouse_subclass_cluster_total$subclass,
  mouse_subclass_cluster_total$merge_regions, sep = "_"
)

dir.create(file.path(combo_root, "mouse_table", "table2"), showWarnings = FALSE)

n_pairs <- nrow(mouse_cluster_pairship_new_total)

if (n_pairs > 0) {
  workers2 <- max(1L, min(n_pairs, detectCores() - 1L))
  cl2 <- makeCluster(workers2)
  registerDoParallel(cl2)
  pair_overlap_dt <- foreach(i = seq_len(n_pairs),
                             .packages = c("stringr","data.table"),
                             .export   = c("mouse_subclass_cluster_total")) %dopar% {
                               pr <- mouse_cluster_pairship_new_total[i, ]
                               xlab <- pr$label.1
                               ylab <- pr$label.2
                               xrow <- mouse_subclass_cluster_total[
                                 mouse_subclass_cluster_total$label %in% xlab, , drop = FALSE]
                               yrow <- mouse_subclass_cluster_total[
                                 mouse_subclass_cluster_total$label %in% ylab, , drop = FALSE]
                               x_cells <- Reduce(
                                 union,
                                 list(
                                   unique(as.character(str_split(xrow$Glut_Neruon_cell_ids, "[,]", simplify = TRUE))),
                                   unique(as.character(str_split(xrow$GABA_Neruon_cell_ids,  "[,]", simplify = TRUE))),
                                   unique(as.character(str_split(xrow$Non_Neruon_cell_ids,   "[,]", simplify = TRUE)))
                                 )
                               )
                               x_cells <- x_cells[x_cells != "" & !is.na(x_cells)]
                               y_cells <- Reduce(
                                 union,
                                 list(
                                   unique(as.character(str_split(yrow$Glut_Neruon_cell_ids, "[,]", simplify = TRUE))),
                                   unique(as.character(str_split(yrow$GABA_Neruon_cell_ids,  "[,]", simplify = TRUE))),
                                   unique(as.character(str_split(yrow$Non_Neruon_cell_ids,   "[,]", simplify = TRUE)))
                                 )
                               )
                               y_cells <- y_cells[y_cells != "" & !is.na(y_cells)]
                               overlap <- length(intersect(x_cells, y_cells))
                               un      <- length(union(x_cells, y_cells))
                               data.table(
                                 label.1      = xlab,
                                 label.2      = ylab,
                                 overlap_cell = overlap,
                                 union_cell   = un
                               )
                             }
  stopCluster(cl2)
  pair_overlap <- rbindlist(pair_overlap_dt, use.names = TRUE, fill = TRUE)
} else {
  pair_overlap <- data.table()
}

Table_1_min <- Table_1_total_cell

colnames(pair_overlap)[1] <- "label"
Table_2_total_cell <- merge(Table_1_min, pair_overlap, by = "label")
colnames(Table_2_total_cell)[1:11] <- paste0("cluster.1_", colnames(Table_2_total_cell)[1:11])

colnames(Table_2_total_cell)[12] <- "label"
Table_2_total_cell <- merge(Table_2_total_cell, Table_1_min, by = "label")
colnames(Table_2_total_cell)[c(1, 15:24)] <- paste0("cluster.2_", colnames(Table_2_total_cell)[c(1, 15:24)])

Table_2_total_cell <- Table_2_total_cell[, c(2:12, 1, 15:24, 13:14)]
Table_2_total_cell$jaccard <- Table_2_total_cell$overlap_cell /
  ifelse(Table_2_total_cell$union_cell == 0, NA, Table_2_total_cell$union_cell)
Table_2_total_cell$cluster.1.overlap.percent <-
  Table_2_total_cell$overlap_cell /
  ifelse(Table_2_total_cell$cluster.1_total_cell_num == 0, NA,
         Table_2_total_cell$cluster.1_total_cell_num)
Table_2_total_cell$cluster.2.overlap.percent <-
  Table_2_total_cell$overlap_cell /
  ifelse(Table_2_total_cell$cluster.2_total_cell_num == 0, NA,
         Table_2_total_cell$cluster.2_total_cell_num)

fwrite(Table_2_total_cell,
       file.path(combo_root, "mouse_table", "table2", "Table_2_total_cell.txt"),
       sep = "\t", quote = FALSE, row.names = FALSE)


rm(list = ls()); gc()
setwd('E:/zaw/2511/mouseMerfish_zhuang_subclass')
outdir <- NULL

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(scales)
  library(ggrepel)
  library(ggtext)
  library(ggpubr)
})

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
base_dir <- "E:/zaw/2511/mouseMerfish_zhuang_subclass"
combo_root <- file.path(base_dir, subdir_tag)
dir.create(combo_root, recursive = TRUE, showWarnings = FALSE)

outdir <- file.path(combo_root, "figures2_subclass")
dir.create(outdir, showWarnings = FALSE)

tab <- read.delim(file.path(combo_root, "mouse_table", "table2", "Table_2_total_cell.txt"))

tab$pair_type <- paste(tab$cluster.1_cell_Neruon_type,
                       tab$cluster.2_cell_Neruon_type, sep = '-')
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

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(scales)
})

tab <- read.delim(file.path(combo_root, "mouse_table", "table2", "Table_2_total_cell.txt"))

tab <- tab %>%
  mutate(pair_type_raw = paste(cluster.1_cell_Neruon_type,
                               cluster.2_cell_Neruon_type, sep = "-"))

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
  "0–20%"   = "#EAE2BD",
  "20–40%"  = "#C7E2D5",
  "40–60%"  = "#59A8A2",
  "60–80%"  = "#2E74A6",
  "80–100%" = "#234663"
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

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(scales)
})

tab <- tab %>%
  mutate(pair_type_raw = paste(cluster.1_cell_Neruon_type,
                               cluster.2_cell_Neruon_type, sep = "-"))

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

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(ggExtra)
  library(viridis)
})

SUB_IS_CLUSTER2 <- TRUE

dfC_raw <- tab %>%
  filter(
    `cluster.1_cell_Neruon_type` %in% c("Glut","Gaba"),
    `cluster.2_cell_Neruon_type` %in% c("Glut","Gaba"),
    `cluster.1_cell_Neruon_type` != `cluster.2_cell_Neruon_type`
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

tab <- read.delim(file.path(combo_root, "mouse_table", "table2", "Table_2_total_cell.txt"))

df_all <- tab %>%
  mutate(
    pair_type = paste(cluster.1_cell_Neruon_type, cluster.2_cell_Neruon_type, sep = "-"),
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


rm(list = ls()); gc()

base_dir <- "E:/zaw/2511/mouseMerfish_zhuang_subclass"

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
