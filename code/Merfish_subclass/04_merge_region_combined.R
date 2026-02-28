rm(list = ls()); gc()

library(vroom)
library(stringr)
library(data.table)
library(doParallel)
library(foreach)
library(dplyr)
library(parallel)

base_dir         <- "E:/zaw/2511/mouseMerfish_zhuang_subclass/"
merge_region_dir <- file.path(base_dir, "merge_region")
mapping_path     <- file.path(base_dir, "Merfish_brain_cell_type_subclass.txt")
setwd(base_dir)

## ========== 找到所有 merge_region 结果文件 ==========
file_names <- list.files(
  merge_region_dir,
  pattern = "table_cell_id\\.txt$",
  full.names = FALSE
)

cat("[INFO] Found", length(file_names), "files in", merge_region_dir, "\n")

## ========== 并行读取所有 table_cell_id.txt 并保证字符串 ==========
workers <- max(1L, min(length(file_names), parallel::detectCores() - 1L))
cat("[INFO] Using workers:", workers, "\n")

cl <- makeCluster(workers)
registerDoParallel(cl)

## 每个 worker 内 data.table 设为单线程
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
  
  ## 保证关键标识列为字符
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

## 合并为一个大表
mouse_subclass_cluster_total <- rbindlist(
  mouse_list,
  use.names = TRUE,
  fill      = TRUE
)
setDT(mouse_subclass_cluster_total)

## 再保险：subclass 一定是字符
if ("subclass" %in% names(mouse_subclass_cluster_total)) {
  mouse_subclass_cluster_total[, subclass := as.character(subclass)]
}

cat("[INFO] Total rows in mouse_subclass_cluster_total:",
    nrow(mouse_subclass_cluster_total), "\n")

## ========== 下面保持你原来的逻辑不变 ==========

print(length(unique(mouse_subclass_cluster_total$subclass)))

sub_vec <- unique(mouse_subclass_cluster_total$subclass)

## 构建 subclass → cell_Neuron_type, subclass_sim 的映射
df <- data.frame(subclass = sub_vec, stringsAsFactors = FALSE)

# 取最后一个 token（例如 "Dopa-Gaba"、"Glut"、"Gaba"）
last_token <- stringr::word(df$subclass, -1)

# 1）根据最后一个 token 是否以 Glut / Gaba 结尾判断神经元类型
df$cell_Neuron_type <- dplyr::case_when(
  stringr::str_detect(last_token, "Glut$") ~ "Glut",
  stringr::str_detect(last_token, "Gaba$") ~ "Gaba",
  TRUE                                     ~ "NonNeuron"
)

# 2）构建 subclass_sim：
#    先去掉开头编号（例如 "035 "）
tmp <- stringr::str_remove(df$subclass, "^\\S+\\s+")

# 先默认简写 = 去掉编号后的字符串
df$subclass_sim <- tmp

# 对 Glut / Gaba 的行，把末尾的 " Glut" / " Gaba" / "-Glut" / "-Gaba" 去掉
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

## 写出 mapping 文件
write.table(
  out_df,
  file      = mapping_path,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

## 读回 mapping 文件并合并
Merfish_brain_cell_type <- read.delim(
  mapping_path,
  check.names      = FALSE,
  stringsAsFactors = FALSE
)
setDT(Merfish_brain_cell_type)  # ★ 关键：变成 data.table

setnames(
  Merfish_brain_cell_type,
  old = names(Merfish_brain_cell_type),
  new = c("subclass", "cell_Neuron_type", "subclass_sim")
)

## 保证 mapping 这三列也是字符串
Merfish_brain_cell_type[, subclass         := as.character(subclass)]
Merfish_brain_cell_type[, cell_Neuron_type := as.character(cell_Neuron_type)]
Merfish_brain_cell_type[, subclass_sim     := as.character(subclass_sim)]

mouse_subclass_cluster_total <- merge(
  mouse_subclass_cluster_total,
  Merfish_brain_cell_type,
  by    = "subclass",
  all.x = TRUE, all.y = FALSE
)

## 调整列顺序
front_cols <- c(
  "slide", "layer", "region", "merge_regions",
  "cell_Neuron_type", "subclass", "subclass_sim"
)
front_cols <- intersect(front_cols, names(mouse_subclass_cluster_total))
other_cols <- setdiff(names(mouse_subclass_cluster_total), front_cols)
mouse_subclass_cluster_total <- mouse_subclass_cluster_total[
  , .SD, .SDcols = c(front_cols, other_cols)
]

print(table(mouse_subclass_cluster_total$cell_Neuron_type))

## 最终写出总表
write.table(
  mouse_subclass_cluster_total,
  file      = file.path(base_dir, "mouse_subclass_cluster_total.txt"),
  quote     = FALSE,
  sep       = "\t",
  row.names = FALSE
)

cat("[INFO] Done.\n")

