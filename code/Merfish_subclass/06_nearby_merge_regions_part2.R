###########################
## Part 1: 距离结果抽查校验
###########################
rm(list = ls()); gc()
setwd('E:/zaw/2511/mouseMerfish_zhuang_subclass')

suppressPackageStartupMessages({
  library(stringr)
  library(data.table)
  library(doParallel)
  library(foreach)
  library(parallel)
})

## 限制 BLAS/OMP 线程，避免和并行抢核
Sys.setenv(
  OPENBLAS_NUM_THREADS = "1",
  OMP_NUM_THREADS      = "1",
  MKL_NUM_THREADS      = "1",
  GOTO_NUM_THREADS     = "1"
)
data.table::setDTthreads(1L)

file_names <- list.files('./nearby_merge_regions_new/merge_regions_nearby_dist/')
file_names <- file_names[grep('dist.txt$', file_names)]

cat("[INFO] dist files:", length(file_names), "\n")

## 并行设置
workers <- max(1L, min(length(file_names), detectCores() - 1L))
cat("[INFO] Part1 workers:", workers, "\n")
cl1 <- makeCluster(workers)
registerDoParallel(cl1)

## data.table 安全合并函数
dt_rbind <- function(...) data.table::rbindlist(list(...), use.names = TRUE, fill = TRUE)

## 每个文件跑 1000 次随机抽查，返回本文件的 error rows
error_list_all <- foreach(fn = file_names,
                          .combine      = dt_rbind,
                          .multicombine = TRUE,
                          .packages     = c("stringr","data.table")) %dopar% {
                            
                            ## ---------- 读距离结果 ----------
                            dist_path <- file.path('./nearby_merge_regions_new/merge_regions_nearby_dist/', fn)
                            file_temp <- read.delim(
                              dist_path,
                              sep              = '\t',
                              header           = TRUE,
                              stringsAsFactors = FALSE,
                              check.names      = FALSE
                            )
                            if (!nrow(file_temp)) return(data.table())
                            
                            ## 仅对 *cell_id* / *cell_label* 相关列强制转字符
                            id_cols1 <- grep("cell_label|cell_id", names(file_temp), ignore.case = TRUE, value = TRUE)
                            for (cc in id_cols1) {
                              file_temp[[cc]] <- as.character(file_temp[[cc]])
                            }
                            
                            ## ---------- 读对应切片坐标文件，保证 cell_label/cell_id 为字符 ----------
                            slice_id <- str_split(fn, "\\.merge", simplify = TRUE)[, 1]
                            cell_path <- paste0(
                              'E:/zaw/2511/mouseMerfish_zhuang_subclass/neocortex_new/',
                              slice_id, '.txt'
                            )
                            
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
                            ## cell_label / cell_id 强制字符
                            if ("cell_label" %in% names(file_temp_cell)) {
                              file_temp_cell$cell_label <- as.character(file_temp_cell$cell_label)
                            }
                            if ("cell_id" %in% names(file_temp_cell)) {
                              file_temp_cell$cell_id <- as.character(file_temp_cell$cell_id)
                            }
                            
                            ## ---------- 随机抽查 1000 行 ----------
                            n_row <- nrow(file_temp)
                            if (n_row == 0L) return(data.table())
                            
                            iter <- min(1000L, n_row)  ## 行很少时不超出
                            err_list <- vector("list", iter)
                            k <- 0L
                            
                            for (i in seq_len(iter)) {
                              temp1 <- file_temp[sample.int(n_row, 1L), , drop = FALSE]
                              
                              temp1.x <- file_temp_cell[file_temp_cell$cell_label %in% temp1$cluster.1.cell_id, ]
                              temp1.y <- file_temp_cell[file_temp_cell$cell_label %in% temp1$cluster.2.cell_id, ]
                              
                              ## 如果有一侧没取到坐标，直接跳过这一轮
                              if (nrow(temp1.x) == 0L || nrow(temp1.y) == 0L) next
                              
                              temp_xy <- rbind(temp1.x[, c('x','y')],
                                               temp1.y[, c('x','y')])
                              temp1.dist <- dist(as.matrix(temp_xy), method = 'euclidean')
                              
                              a <- all.equal(as.numeric(temp1.dist[1]), as.numeric(temp1$nearby_dist))
                              if (!isTRUE(a)) {
                                k <- k + 1L
                                err_row <- as.data.table(temp1)
                                err_row[, src_file := fn]
                                err_list[[k]] <- err_row
                              }
                            }
                            
                            if (k > 0L) data.table::rbindlist(err_list[1:k], use.names = TRUE, fill = TRUE)
                            else data.table()
                          }

stopCluster(cl1)

error_list <- error_list_all
cat("[INFO] Part1 done. Error rows:", nrow(error_list), "\n")

##########################################
## Part 2: 聚合 sim 结果并生成 pairship 表
##########################################
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

file_names <- list.files('./nearby_merge_regions_new/merge_regions_nearby_dist/')
file_names <- file_names[grep('sim.txt$', file_names)]

cat("[INFO] sim files:", length(file_names), "\n")

workers2 <- max(1L, min(length(file_names), detectCores() - 1L))
cat("[INFO] Part2 workers:", workers2, "\n")
cl2 <- makeCluster(workers2)
registerDoParallel(cl2)

dt_rbind <- function(...) data.table::rbindlist(list(...), use.names = TRUE, fill = TRUE)

## 并行读入所有 sim.txt
mouse_cluster_pairship_list <- foreach(fn = file_names,
                                       .combine      = dt_rbind,
                                       .multicombine = TRUE,
                                       .packages     = c("stringr","data.table")) %dopar% {
                                         
                                         file_temp <- read.delim(
                                           file.path('./nearby_merge_regions_new/merge_regions_nearby_dist/', fn),
                                           sep              = '\t',
                                           header           = TRUE,
                                           stringsAsFactors = FALSE,
                                           check.names      = FALSE
                                         )
                                         
                                         if (!nrow(file_temp)) return(data.table())
                                         
                                         ## 仅对 cell_id / cell_label 相关列强制转字符
                                         id_cols2 <- grep("cell_label|cell_id", names(file_temp), ignore.case = TRUE, value = TRUE)
                                         for (cc in id_cols2) {
                                           file_temp[[cc]] <- as.character(file_temp[[cc]])
                                         }
                                         
                                         as.data.table(file_temp)
                                       }

stopCluster(cl2)

mouse_cluster_pairship <- mouse_cluster_pairship_list
rm(mouse_cluster_pairship_list); gc()

cat("[INFO] nrow(mouse_cluster_pairship) =", nrow(mouse_cluster_pairship), "\n")

## ========= 保持你原来逻辑：逐行 union / 排序 / 拼 label =========
n <- nrow(mouse_cluster_pairship)
mouse_cluster_pairship_new <- vector("list", n)

for (i in seq_len(n)) {
  temp1 <- mouse_cluster_pairship[i, ]
  ## 注意：这里转字符只在临时变量里，不改原数据列类型
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
  file = './mouse_cluster_pairship_new_total.txt',
  quote = FALSE, row.names = FALSE, sep = '\t'
)

## 去重逻辑保持不变（按前 4 列去重）
mouse_cluster_pairship_new.1 <- mouse_cluster_pairship_new[
  !duplicated(mouse_cluster_pairship_new[, c(1,2,5,6)])
]

write.table(
  mouse_cluster_pairship_new.1,
  file = './mouse_cluster_pairship_new.txt',
  quote = FALSE, row.names = FALSE, sep = '\t'
)

cat("[INFO] Part2 done. Unique rows:", nrow(mouse_cluster_pairship_new.1), "\n")

