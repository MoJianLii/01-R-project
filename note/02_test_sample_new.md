# Mouse MERFISH 滑窗富集显著性计算

核心任务：对每个切片对应的 `.RData`（包含 `file_tmp` 与 `cell_windows_layer_d`）进行逐层统计，在每个滑窗位置（`loc`）上计算各 `subclass` 的富集显著性（超几何检验），输出每个切片的 p 值结果表与显著窗口（`p < 0.05`）子集。

------------------------------------------------------------------------

## 1. 脚本目标

对 `./neocortex_new/` 下的每个切片文件：

1. 读取第一步输出的 `.RData`（`./cell_window/ws0.4_ss0.02/{slice}.RData`）。
2. 按 `layer` 分层，对每一层分别计算：
   - `N`：该层总细胞数（unique `cell_label`）；
   - `K`：每个 `subclass` 在该层的细胞数；
   - `n`：每个滑窗位置 `loc` 的窗口样本量（窗口内 unique `cell_label`）；
   - `sum_value`：每个窗口内某个 `subclass` 的观察到数量（unique `cell_label`）；
   - 用 `phyper()` 计算超几何检验右尾 p 值：
     - `p = P(X >= sum_value)`（通过 `lower.tail = FALSE` 且输入 `sum_value - 1` 实现）。
3. 输出每个切片对应的两个 CSV：
   - 全量结果：`{slice}_cell_sliding_window_result_p.csv`
   - 显著窗口：`{slice}_cell_sliding_window_result_p_sig_high.csv`（`p < 0.05`）

------------------------------------------------------------------------

## 2. 输入文件与目录结构

### 2.1 输入数据

- `./neocortex_new/`  
  用于获取切片文件列表（文件名用于定位对应 `.RData`）。

- `./cell_window/ws0.4_ss0.02/{slice}.RData`  
  每个切片的滑窗采样结果（来自上一步脚本），要求至少包含两个对象：
  - `file_tmp`：带注释的原始细胞表（至少包含 `cell_label, layer, subclass`）
  - `cell_windows_layer_d`：滑窗采样表（至少包含 `cell_label, layer, loc`）

> 注意：脚本使用 `setwd('E:/zaw/2511/mouseMerfish_zhuang_subclass')`，因此以上相对路径均相对于该工作目录。

------------------------------------------------------------------------

## 3. 输出结果

### 3.1 输出目录

- `E:/zaw/2511/mouseMerfish_zhuang_subclass/cell_pvalue_new/`

（脚本里也创建了 `cell_sample_new/`，但本脚本实际写文件仅使用 `cell_pvalue_new/`。）

### 3.2 输出文件

对每个切片 `file_name`（去掉扩展名后记为 `{slice}`），输出：

- `{slice}_cell_sliding_window_result_p.csv`  
  全量结果，字段包含：
  - `subclass`
  - `loc`
  - `sum_value`
  - `p`
  - `layer`

- `{slice}_cell_sliding_window_result_p_sig_high.csv`  
  显著窗口子集（`p < 0.05`），字段同上。

------------------------------------------------------------------------

## 4. 关键参数与统计量定义

本脚本没有通过环境变量读取参数；关键依赖来自第一步生成的滑窗结果（其窗口大小/步长已体现在 `cell_window/ws0.4_ss0.02/` 路径中）。

在每个 `layer` 内，超几何检验的统计量含义如下：

- `N`：该层总细胞数（unique `cell_label`）
- `K`：某个 `subclass` 在该层的细胞数
- `n`：某个滑窗 `loc` 内的细胞数（unique `cell_label`）
- `sum_value`：该滑窗 `loc` 内属于某个 `subclass` 的细胞数（unique `cell_label`）

检验问题：总体 `N` 中有 `K` 个“成功”，窗口相当于抽取 `n` 个，观察成功数是否至少为 `sum_value`。  
p 值右尾计算公式：

- `phyper(sum_value - 1, K, N - K, n, lower.tail = FALSE)`

------------------------------------------------------------------------

## 5. 算法逻辑（按 layer 分层计算）

对每个切片：

1. `ft_layer`：提取该层的 `(cell_label, subclass)`（定义总体与每类规模）。
2. `cwd_layer`：提取该层窗口数据 `(cell_label, loc)`（定义每个窗口抽样集合）。
3. 计算：
   - `N <- uniqueN(ft_layer$cell_label)`
   - `Kc <- ft_layer[, .(K = .N), by = subclass]`
   - `nl <- cwd_layer[, .(n = uniqueN(cell_label)), by = loc]`
4. 通过窗口与总体按 `cell_label` 连接，得到 `obs`：每个 `(subclass, loc)` 的观察到数量 `sum_value`。
5. 拼接 `obs + Kc + nl`，并计算 `p`，最终得到该层结果。
6. 合并所有层结果并写出 CSV；再筛选 `p < 0.05` 写出显著子集。

------------------------------------------------------------------------

## 6. 完整脚本

------------------------------------------------------------------------

### 6.1 初始化与依赖加载

这一段做三件事：\
1）清空环境并 gc()，避免上一次会话残留变量影响结果；\
2）设置工作目录；\
3）加载依赖包（data.table 用于高性能 join/汇总），并初始化环境。

``` r
rm(list = ls()); gc()

setwd('E:/zaw/2511/mouseMerfish_zhuang_subclass')

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})
```

------------------------------------------------------------------------

### 6.2 参数读取与输出目录标签

这一段用于创建输出目录：\

- cell_pvalue_new/：保存每个切片的全量 p 值表、以及 p < 0.05 的显著窗口子集；\

- cell_sample_new/：本脚本中未写出文件，但保留目录创建，方便项目结构统一（后续若扩展“窗口抽样/细胞抽样结果”可直接用）。

``` r
out_base <- "E:/zaw/2511/mouseMerfish_zhuang_subclass"
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

out_cell_sample_dir <- file.path(out_base, "cell_sample_new")
out_cell_pvalue_dir <- file.path(out_base, "cell_pvalue_new")
dir.create(out_cell_sample_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_cell_pvalue_dir,  recursive = TRUE, showWarnings = FALSE)
```

------------------------------------------------------------------------

### 6.3 获取切片列表与全局计时

这一段做两件事：\
1）记录脚本总耗时起点 start_time.t；\
2）读取 ./neocortex_new/ 目录下所有切片文件名（这里只是用文件名来定位对应的 .RData，不直接读取 .txt）。

``` r
start_time.t <- Sys.time()
file_names <- list.files('./neocortex_new/')
```

------------------------------------------------------------------------

### 6.4 逐切片处理：读取 .RData 并按 layer 计算富集 p 值（核心）

每个切片文件的处理步骤是：\
1）用切片名（去扩展名）作为输出文件前缀；\
2）加载上一步滑窗脚本输出的 .RData（包含 file_tmp 与 cell_windows_layer_d）；\
3）按 layer 分层，分别构造超几何检验的统计量：N/K/n/sum_value；\
4）用 phyper() 计算右尾 p 值：P(X >= sum_value)；\
5）输出全量结果与显著窗口子集（p < 0.05），并做日志与内存清理。

``` r
for (file_name in file_names) {

  #切片命名 + 载入 Step 1 输出
  ##file_chose_names：作为输出文件名前缀
  file_chose_names <- tools::file_path_sans_ext(file_name)

  #载入 Step 1 的滑窗结果（每切片一个 .RData）
  ##包含：file_tmp（带注释的细胞表）与 cell_windows_layer_d（滑窗采样表）
  load(paste0('./cell_window/ws0.4_ss0.02/', file_chose_names, '.RData'))

  #转为 data.table，便于后续 join 与分组统计
  setDT(file_tmp)
  setDT(cell_windows_layer_d)

  #单切片计时
  start_time <- Sys.time()

  #获取该切片的 layer 列表 + 结果容器
  layer_ids <- sort(unique(file_tmp$layer))
  cell_Sliding_window_result_p <- vector("list", length(layer_ids))

  #按 layer 分层计算：每个 (subclass, loc) 做超几何富集检验
  for (i in seq_along(layer_ids)) {

    lyr <- layer_ids[i]

    #背景总体（该 layer 内每个细胞的 subclass
    ##ft_layer：定义总体 N，并可统计每个 subclass 的总数 K
    ft_layer  <- file_tmp[layer == lyr, .(cell_label, subclass)]

    #窗口成员关系（该 layer 内每个细胞属于哪个 loc
    ##cwd_layer：定义每个窗口 loc 的抽样大小 n
    cwd_layer <- cell_windows_layer_d[layer == lyr, .(cell_label, loc)]

    #\超几何检验的四个核心统计量：N / K / n / sum_value
    ## N：该 layer 总细胞数（unique cell_label）
    N <- uniqueN(ft_layer$cell_label)

    #K：每个 subclass 在该 layer 的总细胞数
    Kc <- ft_layer[, .(K = .N), by = subclass]

    #n：每个窗口 loc 的细胞数（unique cell_label）
    nl <- cwd_layer[, .(n = uniqueN(cell_label)), by = loc]

    #sum_value：每个窗口 loc 内属于某 subclass 的细胞数（unique cell_label）
    ##实现：把窗口表按 cell_label 映射到 subclass，再对 (subclass, loc) 计数
    obs <- cwd_layer[ft_layer, on = .(cell_label), nomatch = 0L][
      , .(sum_value = uniqueN(cell_label)), by = .(subclass, loc)
    ]

    #若该层在窗口与背景之间完全无交集，则写入空表
    if (nrow(obs) == 0L) {
      cell_Sliding_window_result_p[[i]] <- data.table(
        subclass  = character(),
        loc       = character(),
        sum_value = integer(),
        p         = numeric(),
        layer     = character()
      )
      next
    }

    #合并统计量并计算右尾超几何 p 值
    ##右尾：p = P(X >= sum_value)
    ###等价写法：phyper(sum_value - 1, K, N - K, n, lower.tail = FALSE)
    res_layer <- obs[Kc, on = .(subclass)][nl, on = .(loc)][
      , p := phyper(sum_value - 1L, K, N - K, n, lower.tail = FALSE)
    ][
      , .(subclass, loc, sum_value, p)
    ][
      ][]

    #标记 layer，便于追踪来源
    res_layer[, layer := lyr]

    cell_Sliding_window_result_p[[i]] <- res_layer
  }

  #合并该切片所有 layer 的结果
  cell_Sliding_window_result_p_d <- rbindlist(cell_Sliding_window_result_p, use.names = TRUE)

  #写出该切片全量结果
  fwrite(
    cell_Sliding_window_result_p_d,
    file = file.path(
      out_cell_pvalue_dir,
      paste0(file_chose_names, '_', 'cell_sliding_window_result_p.csv')
    )
  )

  #写出该切片显著窗口子集（p < 0.05）
  sig_high <- cell_Sliding_window_result_p_d[p < 0.05]
  fwrite(
    sig_high,
    file = file.path(
      out_cell_pvalue_dir,
      paste0(file_chose_names, '_', 'cell_sliding_window_result_p_sig_high.csv')
    )
  )

  #输出单切片耗时
  end_time <- Sys.time()
  message(sprintf(
    "[%s] done in %.1f mins",
    file_chose_names,
    as.numeric(difftime(end_time, start_time, units = "mins"))
  ))

  #清理内存：避免长循环占用过大
  rm(file_tmp, cell_windows_layer_d,
     cell_Sliding_window_result_p, cell_Sliding_window_result_p_d, sig_high)
  gc()
}
```

------------------------------------------------------------------------

### 6.5 收尾：输出总耗时

这一段负责输出整个脚本的总耗时，便于记录与复现。

``` r
end_time.t <- Sys.time()
execution_time <- difftime(end_time.t, start_time.t, units = "mins")
print(paste('total', execution_time, sep = ' '))
```

