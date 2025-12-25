rm(list = ls()); gc()

setwd('E:/zaw/2511/mouseMerfish_zhuang_subclass')

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

out_base <- "E:/zaw/2511/mouseMerfish_zhuang_subclass"
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
out_cell_sample_dir <- file.path(out_base, "cell_sample_new")
out_cell_pvalue_dir  <- file.path(out_base, "cell_pvalue_new")
dir.create(out_cell_sample_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_cell_pvalue_dir,  recursive = TRUE, showWarnings = FALSE)

start_time.t <- Sys.time()
file_names <- list.files('./neocortex_new/')

for (file_name in file_names) {
  file_chose_names <- tools::file_path_sans_ext(file_name)
  load(paste0('./cell_window/ws0.4_ss0.02/', file_chose_names, '.RData'))

  setDT(file_tmp)
  setDT(cell_windows_layer_d)
  
  start_time <- Sys.time()
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
  
  end_time <- Sys.time()
  message(sprintf("[%s] done in %.1f mins",
                  file_chose_names, as.numeric(difftime(end_time, start_time, units = "mins"))))
  
  rm(file_tmp, cell_windows_layer_d, cell_Sliding_window_result_p, cell_Sliding_window_result_p_d, sig_high)
  gc()
}

end_time.t <- Sys.time()
execution_time <- difftime(end_time.t, start_time.t, units = "mins")
print(paste('total', execution_time, sep = ' '))
