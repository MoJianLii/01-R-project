# Mouse MERFISH 分层滑窗采样

核心任务：对 `./neocortex_new/` 目录下的每个切片文件进行并行处理，合并注释表后，按 **layer 分层**在二维空间（y 方向先滑窗，再在每个 y 窗内对 x 方向滑窗）生成窗口采样结果，并将每个切片的结果保存为 `.RData`。

------------------------------------------------------------------------

## 1. 脚本目标

1.  **一次性读取两个参考注释表**（region/layer 与 cell type 映射）。
2.  **遍历切片文件**（每个切片一个 `.txt`），对每个文件：
    -   读取数据（强制 `cell_label`、`cell_id` 为字符型，避免类型问题）。
    -   与参考注释表进行两次合并（merge）。
    -   按 `layer` 分层做滑窗采样：
        -   先沿 **y** 方向以 `STEP_SIZE` 递进生成 y-window；
        -   在每个 y-window 内，对点按 **x** 排序，再沿 **x** 方向以 `STEP_SIZE` 滑动生成 x-window；
        -   为每个有效窗口生成一个位置标识 `loc = "{x_low}_{y_low}"`，并记录该窗口中出现的 `cell_label`。
3.  将每个切片的合并后数据 `file_tmp` 以及窗口采样结果 `cell_windows_layer_d` 保存为 `.RData` 文件。

------------------------------------------------------------------------

## 2. 输入文件与目录结构

### 2.1 参考注释表

-   `./Merfish_mouse_neocortex_layer_region.txt`\
    用于按 `ccf_region_name` 合并 region/layer 等注释信息。
-   `./Merfish_brain_cell_type.txt`\
    用于按 `subclass` 合并 cell type 等注释信息。

### 2.2 切片数据目录

-   `./neocortex_new/`\
    该目录下的所有文件都会被 `list.files()` 读取并进入并行处理。

> 注意：脚本使用 `setwd('E:/zaw/2511/mouseMerfish_zhuang_subclass')`，因此以上相对路径都相对于该工作目录。

------------------------------------------------------------------------

## 3. 输出结果

### 3.1 输出目录

输出目录由两部分构成： - 基础输出目录：`E:/zaw/2511/mouseMerfish_zhuang_subclass/cell_window` - 参数子目录：`ws{WINDOW_SIZE}_ss{STEP_SIZE}`

最终输出目录类似： `E:/zaw/2511/mouseMerfish_zhuang_subclass/cell_window/ws0.4_ss0.02/`

### 3.2 输出文件

对每个切片 `file_name`，输出： - `{file_name 去掉 .txt}.RData`

每个 `.RData` 内包含两个对象： - `file_tmp`：切片数据与注释表合并后的完整数据（data.table） - `cell_windows_layer_d`：滑窗采样结果（data.table），包含： - `cell_label`：窗口内的细胞标签 - `loc`：窗口位置标识（x_low 与 y_low 拼接） - `layer`：该窗口所属 layer

------------------------------------------------------------------------

## 4. 关键参数

脚本通过环境变量读取参数：

-   `WINDOW_SIZE`：滑窗窗口的边长（默认 `0.4`）
-   `STEP_SIZE`：滑动步长（默认 `0.02`）

这两个值会共同影响： - 窗口覆盖范围（越大越粗，越小越细） - 窗口数量与计算量（步长越小窗口越密，计算越慢） - 输出目录名（用于区分不同参数实验的结果）

------------------------------------------------------------------------

## 5. 并行策略

脚本采用 **按切片文件并行**： - 每个 worker 负责一个 `file_name` - 优点：切片之间互相独立，天然适合并行

同时，脚本显式设置： - `data.table::setDTthreads(1L)`\
目的是避免 `data.table` 在每个 worker 内再开启多线程，导致 CPU 线程过度竞争。

------------------------------------------------------------------------

## 6. 算法逻辑

### 6.1 分层（Layer-wise）

对 `file_tmp$layer` 的每一个唯一值： - 取出该层数据 `file_tmp_layer <- file_tmp[layer == layer.1, ]` - 在该 layer 内单独计算 y 方向范围 `min_y`、`max_y` - 然后在这个 layer 内做滑窗

意义： 不同 layer 的细胞空间分布可能不同，分层处理可以避免跨层混杂。

### 6.2 先 y-window 再 x-window

对每个 y-window：

1\. 计算 y 范围： - `y_low  <- min_y + y_idx * step_size` - `y_high <- y_low + window_size`

2\. 在该 y 范围内取点： - `temp1 <- file_tmp_layer[y >= y_low & y < y_high, .(cell_label, x)]`

3\. 若该 y-window 内点太少（`n1 <= 10`），直接跳过（返回 `NULL`）

4\. 对剩下点按 x 排序后，用双指针方式沿 x 做滑窗： - 对每个 `x_low`，维护窗口 `[x_low, x_high)` 内的点的左右边界 - 若窗口内点数足够（`right - left + 1L > 1L`），就输出这些点对应的 `cell_label` 与窗口位置 `loc`

最终把每个 layer 的所有窗口结果合并，再把所有 layer 合并成 `cell_windows_layer_d`。

------------------------------------------------------------------------

## 7. 完整脚本

------------------------------------------------------------------------

### 7.1 初始化与依赖加载

这一段做三件事：\
1）清空环境并 `gc()`，避免上一次会话残留变量影响结果；\
2）设置工作目录；\
3）加载依赖包，并记录脚本开始时间。

``` r
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

start_time.t <- Sys.time()  # Script start time
```

------------------------------------------------------------------------

### 7.2 参数读取与输出目录标签

这一段用于从环境变量读取窗口参数：\
- `WINDOW_SIZE`：窗口大小\
- `STEP_SIZE`：步长

并通过 `fmt()` 把参数格式化成目录标签（例如 `ws0.4_ss0.02`），避免出现 `0.400000` 这种不美观的名字。

``` r
## ================== Parameters and utilities ==================
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
```

------------------------------------------------------------------------

### 7.3 读取参考注释表与输入切片列表

这一段做两件事：\
1）读取两张参考注释表；\
2）读取 `./neocortex_new/` 下所有切片文件名，后续按文件并行处理。

``` r
cat("[DEBUG] Reading Merfish reference files once...\n")

Merfish_mouse_neocortex_layer_region <- read.delim(
  './Merfish_mouse_neocortex_layer_region.txt',
  header = TRUE, stringsAsFactors = FALSE, check.names = FALSE
)

Merfish_brain_cell_type <- read.delim(
  './Merfish_brain_cell_type.txt',
  header = TRUE, stringsAsFactors = FALSE, check.names = FALSE
)

in_dir  <- './neocortex_new/'
file_names <- list.files(in_dir)
cat("[INFO] Found", length(file_names), "files in", in_dir, "\n")
```

------------------------------------------------------------------------

### 7.4 创建输出目录

输出目录结构是：\
- 基础目录：`.../cell_window/`\
- 参数子目录：`ws*_ss*`

这样不同窗口参数的结果不会混在一起。

``` r
out_dir_base <- 'E:/zaw/2511/mouseMerfish_zhuang_subclass/cell_window'
out_dir <- file.path(out_dir_base, subdir_tag)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("[INFO] Output dir: ", out_dir, "\n")
cat(sprintf("[INFO] window_size=%s step_size=%s\n", fmt(window_size), fmt(step_size)))
```

------------------------------------------------------------------------

### 7.5 并行设置

这里的策略是"一个切片文件 = 一个并行任务"。\
同时把每个 worker 内的 `data.table` 线程数设为 1，避免"并行 + 多线程"叠加造成 CPU 竞争。

``` r
n_cores <- max(1L, min(parallel::detectCores(logical = TRUE) - 1L, length(file_names)))
cat("[INFO] Using", n_cores, "cores via doParallel\n")

cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

## To avoid over-subscription (data.table threads inside workers), force 1 thread per worker
parallel::clusterEvalQ(cl, {
  data.table::setDTthreads(1L)
  NULL
})
```

------------------------------------------------------------------------

### 7.6 并行处理每个切片文件（核心）

每个切片文件的处理步骤是：

1）**读表头**构造 `colClasses`，保证 `cell_label` / `cell_id` 从读入开始就是字符型，这是为了防止后面 merge 或 data.table 转换后类型变掉；\
2）读主表 → 再次强制字符型；\
3）两次 merge：\
- `by = ccf_region_name` 合并 layer/region 注释\
- `by = class` 合并 cell type 注释\
4）按 `layer` 分层：对每一层做滑窗采样（先 y 再 x），并生成窗口位置 `loc = "{x_low}_{y_low}"`；\
5）把 `file_tmp`（完整合并后的表）和 `cell_windows_layer_d`保存为 `.RData`。

``` r
res_log <- foreach(
  file_name = file_names,
  .packages = c("data.table", "stringr")
) %dopar% {

  start_time <- Sys.time()
  msg_prefix <- paste0("[", file_name, "] ")

  fpath <- file.path(in_dir, file_name)
  if (!file.exists(fpath)) {
    return(paste0(msg_prefix, "ERROR: File not found: ", fpath))
  }

  cat("\n[INFO] Processing file:", file_name, "\n")

  header_dt <- read.delim(
    fpath,
    nrows            = 0,
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )
  col_nms <- names(header_dt)
  col_cls <- rep(NA_character_, length(col_nms))
  names(col_cls) <- col_nms

  if ("cell_label" %in% col_nms) col_cls["cell_label"] <- "character"
  if ("cell_id"    %in% col_nms) col_cls["cell_id"]    <- "character"

  cat("[DEBUG]", msg_prefix, "Reading main data file with colClasses...\n")
  file_tmp <- read.delim(
    fpath,
    stringsAsFactors = FALSE,
    check.names      = FALSE,
    colClasses       = col_cls
  )
  cat("[DEBUG]", msg_prefix, "file_tmp loaded:",
      nrow(file_tmp), "rows x", ncol(file_tmp), "cols\n")

  if ("cell_label" %in% names(file_tmp)) file_tmp$cell_label <- as.character(file_tmp$cell_label)
  if ("cell_id"    %in% names(file_tmp)) file_tmp$cell_id    <- as.character(file_tmp$cell_id)

  cat("[DEBUG]", msg_prefix, "Merging Merfish_mouse_neocortex_layer_region...\n")
  file_tmp <- merge(file_tmp, Merfish_mouse_neocortex_layer_region, by = 'ccf_region_name')
  cat("[DEBUG]", msg_prefix, "After first merge:", dim(file_tmp)[1], "x", dim(file_tmp)[2], "\n")

  cat("[DEBUG]", msg_prefix, "Merging Merfish_brain_cell_type...\n")
  file_tmp <- merge(file_tmp, Merfish_brain_cell_type, by = 'class')
  cat("[DEBUG]", msg_prefix, "After second merge:", dim(file_tmp)[1], "x", dim(file_tmp)[2], "\n")

  if ("cell_label" %in% names(file_tmp)) file_tmp$cell_label <- as.character(file_tmp$cell_label)
  if ("cell_id"    %in% names(file_tmp)) file_tmp$cell_id    <- as.character(file_tmp$cell_id)

  file_chose_names <- sub("\\.txt$", "", file_name)

  setDT(file_tmp)

  cell_windows_function <- function(y_idx) {
    y_low  <- min_y + y_idx * step_size
    y_high <- y_low + window_size

    temp1 <- file_tmp_layer[y >= y_low & y < y_high, .(cell_label, x)]
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
    cat("[DEBUG]", msg_prefix, "Processing layer:", layer.1, "\n")
    file_tmp_layer <- file_tmp[layer == layer.1, ]

    if (nrow(file_tmp_layer) == 0) {
      cat("[WARN]", msg_prefix, "Layer", layer.1, "has no data.\n")
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

  save(
    file_tmp,
    cell_windows_layer_d,
    file = file.path(out_dir, paste0(file_chose_names, '.RData'))
  )

  log_msg <- paste0(
    "[DONE] ", file_name, " in ", round(as.numeric(execution_time), 2), " mins | ",
    "real=",   length(unique(file_tmp$cell_label)),
    " sample=", length(unique(cell_windows_layer_d$cell_label))
  )
  cat(log_msg, "\n")
  log_msg
}
```

------------------------------------------------------------------------

### 7.7 收尾：关闭并行并输出总耗时

这部分负责释放并行资源，并输出总耗时与日志。

``` r
parallel::stopCluster(cl)

end_time.t <- Sys.time()
execution_time <- difftime(end_time.t, start_time.t, units = "mins")
cat("[INFO] Total execution time:", execution_time, "mins\n")
cat("[INFO] Script finished.\n")
print(res_log)
```
