rm(list = ls()); gc()

library(vroom)
library(stringr)
library(data.table)
library(doParallel)
library(foreach)
library(dplyr)
library(parallel)

## 限制 BLAS / OMP 线程，避免和并行抢核
Sys.setenv(
  OPENBLAS_NUM_THREADS = "1",
  OMP_NUM_THREADS      = "1",
  MKL_NUM_THREADS      = "1",
  GOTO_NUM_THREADS     = "1"
)
data.table::setDTthreads(1L)

setwd('E:/zaw/2511/mouseMerfish_zhuang_subclass')

start_time.t <- Sys.time()

## ================== 第一段：merge_region_cell_id 生成 ==================
file_names <- list.files('E:/zaw/2511/mouseMerfish_zhuang_subclass/neocortex_new/')

merged_regions_total <- list.files('E:/zaw/2511/mouseMerfish_zhuang_subclass/merge_region')
merged_regions_total <- merged_regions_total[grep('cell_id.txt', merged_regions_total)]

## 并行 backend
workers <- max(1L, min(length(file_names), detectCores() - 1L))
cat("[INFO] Part1 workers:", workers, "\n")
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
              ## 读入 merge_region 结果，强制字符串
              temp1.cell.id <- read.delim(
                file.path('E:/zaw/2511/mouseMerfish_zhuang_subclass/merge_region', x),
                sep             = '\t',
                header          = TRUE,
                stringsAsFactors = FALSE,
                check.names      = FALSE
              )
              
              ## 关键列保证字符（包括 ID 串）
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
              
              ## 三类细胞 ID 都按字符串拆分
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
              
              ## cell_id 保证字符
              temp3 <- data.frame(cell_id = as.character(temp2.cell.id),
                                  stringsAsFactors = FALSE)
              
              ## 其它列顺序不变
              temp3 <- cbind(temp2[, c(1:4, 6, 5)], temp3)
              merge_region_cell_id <- rbind(merge_region_cell_id, temp3)
            }
          }
          
          ## 写出 merge_region_cell_id（cell_id 为字符串）
          write.table(
            merge_region_cell_id,
            file = paste(
              'E:/zaw/2511/mouseMerfish_zhuang_subclass/nearby_merge_regions_new/merge_cell_id/',
              file_name, '.merge_region_cell_id.txt', sep = ''
            ),
            sep       = '\t',
            quote     = FALSE,
            row.names = FALSE
          )
          print(file_name)
          
          file_name
        }

stopCluster(cl1)

## ================== 清环境，再跑第二段（保持你原来的习惯） ==================
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

file_names <- list.files('E:/zaw/2511/mouseMerfish_zhuang_subclass/neocortex_new/')

workers2 <- max(1L, min(length(file_names), detectCores() - 1L))
cat("[INFO] Part2 workers:", workers2, "\n")
cl2 <- makeCluster(workers2)
registerDoParallel(cl2)

foreach(i = file_names,
        .combine  = c,
        .packages = c("stringr","data.table")) %dopar% {
          
          file_name <- str_split(i, "\\.txt", simplify = TRUE)[, 1]
          
          ## ---------- 读 neocortex_new，强制 cell_label 为字符 ----------
          fpath_tmp <- paste(
            'E:/zaw/2511/mouseMerfish_zhuang_subclass/neocortex_new/', i, sep = ''
          )
          
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
          
          ## ---------- 读 merge_region_cell_id，强制 cell_id 为字符 ----------
          fpath_merge <- paste(
            'E:/zaw/2511/mouseMerfish_zhuang_subclass/nearby_merge_regions_new/merge_cell_id/',
            file_name, '.merge_region_cell_id.txt', sep = ''
          )
          
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
          
          ## 如果这个切片根本没有任何 merge_region，直接写空结果并返回
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
                
                ## 这里 cell_label 和 cell_id 都是字符串
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
                  print(paste(file_name, ':', x_loc, ': ', a, '/', length(loc_vec), sep = ''))
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
          
          ## 写结果
          write.table(
            merge_regions_nearby_dist,
            file = paste(
              'E:/zaw/2511/mouseMerfish_zhuang_subclass/nearby_merge_regions_new/merge_regions_nearby_dist/',
              file_name, '.merge_regions_nearby_dist.txt', sep = ''
            ),
            sep       = '\t',
            quote     = FALSE,
            row.names = FALSE
          )
          
          write.table(
            merge_regions_nearby_dist_sim,
            file = paste(
              'E:/zaw/2511/mouseMerfish_zhuang_subclass/nearby_merge_regions_new/merge_regions_nearby_dist/',
              file_name, '.merge_regions_nearby_dist_sim.txt', sep = ''
            ),
            sep       = '\t',
            quote     = FALSE,
            row.names = FALSE
          )
          
          file_name
        }

stopCluster(cl2)
