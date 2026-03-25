#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

PAIR_ORDER <- c("Gaba-Gaba", "Gaba-Glut", "Glut-Gaba", "Glut-Glut")
BINS <- c("0-20%", "20-40%", "40-60%", "60-80%", "80-100%")
BIN_COLORS <- c(
  "0-20%" = "#d9d2ad",
  "20-40%" = "#b7d4ca",
  "40-60%" = "#57a9a5",
  "60-80%" = "#3d80ad",
  "80-100%" = "#294f72"
)

get_script_dir <- function() {
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_all, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

normalize_bin_label <- function(x) {
  x <- as.character(x)
  x <- gsub("[[:space:]]+", "", x)
  x <- gsub("[\u2012\u2013\u2014\u2212]", "-", x)
  x
}

parse_bool <- function(x) {
  x <- tolower(trimws(as.character(x)))
  if (x %in% c("1", "true", "t", "yes", "y")) return(TRUE)
  if (x %in% c("0", "false", "f", "no", "n")) return(FALSE)
  stop(sprintf("Invalid boolean value: %s", x), call. = FALSE)
}

parse_args <- function(defaults) {
  args <- commandArgs(trailingOnly = TRUE)
  out <- defaults
  i <- 1L
  while (i <= length(args)) {
    key <- args[i]
    if (!startsWith(key, "--")) {
      stop(sprintf("Invalid argument format: %s", key), call. = FALSE)
    }
    key_name <- sub("^--", "", key)
    key_name <- gsub("-", "_", key_name, fixed = TRUE)
    if (!key_name %in% names(out)) {
      stop(sprintf("Unknown argument: %s", key), call. = FALSE)
    }
    if (i == length(args)) {
      stop(sprintf("Missing value for argument: %s", key), call. = FALSE)
    }
    out[[key_name]] <- args[i + 1L]
    i <- i + 2L
  }
  out
}

resolve_python <- function(python_arg) {
  if (nzchar(trimws(python_arg))) {
    return(python_arg)
  }
  py3 <- Sys.which("python3")
  if (nzchar(py3)) return(py3)
  py <- Sys.which("python")
  if (nzchar(py)) return(py)
  stop("No python executable found. Set --python explicitly.", call. = FALSE)
}

run_cmd <- function(command, args) {
  message(sprintf("[RUN] %s %s", command, paste(shQuote(args), collapse = " ")))
  status <- system2(command, args = args, stdout = "", stderr = "")
  if (!identical(as.integer(status), 0L)) {
    stop(sprintf("Command failed with exit code %s", as.character(status)), call. = FALSE)
  }
}

link_or_copy <- function(src, dst) {
  if (file.exists(dst)) return(FALSE)

  ok_link <- FALSE
  suppressWarnings({
    ok_link <- isTRUE(file.link(src, dst))
  })
  if (ok_link && file.exists(dst)) return(TRUE)

  ok_copy <- isTRUE(file.copy(src, dst, overwrite = FALSE))
  if (!ok_copy || !file.exists(dst)) {
    stop(sprintf("Failed to materialize file: %s -> %s", src, dst), call. = FALSE)
  }
  TRUE
}

collect_sample_sources <- function(merge_region_root, start_folder, end_folder) {
  chip_dirs <- list.dirs(merge_region_root, recursive = FALSE, full.names = TRUE)
  if (length(chip_dirs) == 0L) {
    return(data.table())
  }

  all_rows <- vector("list", length(chip_dirs))
  for (i in seq_along(chip_dirs)) {
    chip_dir <- chip_dirs[i]
    chip_name <- basename(chip_dir)
    sample_dirs <- list.dirs(chip_dir, recursive = FALSE, full.names = TRUE)
    sample_names <- basename(sample_dirs)
    keep_name <- grepl("^sample_[0-9]+$", sample_names)
    sample_dirs <- sample_dirs[keep_name]
    sample_names <- sample_names[keep_name]
    if (length(sample_dirs) == 0L) {
      all_rows[[i]] <- NULL
      next
    }

    sample_idx <- as.integer(sub("^sample_", "", sample_names))
    keep_range <- !is.na(sample_idx) & sample_idx >= start_folder & (is.na(end_folder) | sample_idx <= end_folder)
    sample_dirs <- sample_dirs[keep_range]
    sample_names <- sample_names[keep_range]
    sample_idx <- sample_idx[keep_range]

    if (length(sample_dirs) == 0L) {
      all_rows[[i]] <- NULL
      next
    }

    all_rows[[i]] <- data.table(
      chip = chip_name,
      sample_name = sample_names,
      sample_idx = sample_idx,
      sample_dir = sample_dirs
    )
  }

  rbindlist(all_rows, use.names = TRUE, fill = TRUE)
}

prepare_merge_samples <- function(merge_region_root, merge_samples_dir, start_folder, end_folder, rebuild) {
  if (!dir.exists(merge_region_root)) {
    stop(sprintf("merge_region root not found: %s", merge_region_root), call. = FALSE)
  }

  if (rebuild && dir.exists(merge_samples_dir)) {
    message(sprintf("[INFO] Removing existing directory: %s", merge_samples_dir))
    unlink(merge_samples_dir, recursive = TRUE, force = TRUE)
  }
  dir.create(merge_samples_dir, recursive = TRUE, showWarnings = FALSE)

  src_tbl <- collect_sample_sources(merge_region_root, start_folder, end_folder)
  if (nrow(src_tbl) == 0L) {
    stop("No sample_* directories found in selected range under merge_region root.", call. = FALSE)
  }

  linked <- 0L
  skipped <- 0L
  sample_names <- unique(src_tbl$sample_name)
  sample_names <- sample_names[order(as.integer(sub("^sample_", "", sample_names)))]

  for (sname in sample_names) {
    target_dir <- file.path(merge_samples_dir, sname)
    dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)

    rows <- src_tbl[sample_name == sname]
    for (k in seq_len(nrow(rows))) {
      src_dir <- rows$sample_dir[k]
      chip <- rows$chip[k]
      files <- list.files(src_dir, pattern = "_merged_regions_table_cell_id\\.txt$", full.names = TRUE)
      if (length(files) == 0L) next

      for (src_file in files) {
        dst <- file.path(target_dir, basename(src_file))
        if (file.exists(dst)) {
          dst <- file.path(target_dir, paste0(chip, "__", basename(src_file)))
        }
        if (file.exists(dst)) {
          skipped <- skipped + 1L
          next
        }
        if (link_or_copy(src_file, dst)) {
          linked <- linked + 1L
        } else {
          skipped <- skipped + 1L
        }
      }
    }
  }

  list(
    samples = length(sample_names),
    files_materialized = linked,
    files_skipped = skipped
  )
}

compute_null_mean_table <- function(sample_csv, out_csv) {
  dt <- fread(sample_csv)
  need <- c("sample", "pair_type", "overlap_bin", "percent")
  if (!all(need %in% names(dt))) {
    stop(sprintf("Missing columns in %s: %s", sample_csv, paste(setdiff(need, names(dt)), collapse = ", ")), call. = FALSE)
  }

  dt[, sample := as.character(sample)]
  dt[, pair_type := as.character(pair_type)]
  dt[, overlap_bin := normalize_bin_label(overlap_bin)]
  dt[, percent := as.numeric(percent)]

  dt <- dt[pair_type %in% PAIR_ORDER & overlap_bin %in% BINS]

  summary_dt <- dt[
    ,
    .(
      sample_n = uniqueN(sample),
      mean_percent = mean(percent, na.rm = TRUE),
      sd_percent = if (.N > 1L) sd(percent, na.rm = TRUE) else NA_real_
    ),
    by = .(pair_type, overlap_bin)
  ]

  full_grid <- CJ(pair_type = PAIR_ORDER, overlap_bin = BINS, unique = TRUE)
  summary_dt <- merge(full_grid, summary_dt, by = c("pair_type", "overlap_bin"), all.x = TRUE, sort = FALSE)
  summary_dt[is.na(sample_n), sample_n := 0L]

  fwrite(summary_dt, out_csv)
  summary_dt
}

plot_null_mean_base <- function(summary_dt, out_svg) {
  plot_dt <- copy(summary_dt)
  plot_dt[, pair_type := factor(pair_type, levels = PAIR_ORDER)]
  plot_dt[, overlap_bin := factor(overlap_bin, levels = BINS)]

  cast_mean <- dcast(plot_dt, overlap_bin ~ pair_type, value.var = "mean_percent")
  cast_sd <- dcast(plot_dt, overlap_bin ~ pair_type, value.var = "sd_percent")

  mat <- as.matrix(cast_mean[, -1, with = FALSE])
  rownames(mat) <- as.character(cast_mean$overlap_bin)
  err <- as.matrix(cast_sd[, -1, with = FALSE])
  err[is.na(err)] <- 0

  y_max <- max(0, mat + err, na.rm = TRUE)
  if (!is.finite(y_max) || y_max <= 0) y_max <- 1

  svg(filename = out_svg, width = 12, height = 7)
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  bp <- barplot(
    mat,
    beside = TRUE,
    col = BIN_COLORS[rownames(mat)],
    border = "#222222",
    ylim = c(0, y_max * 1.15),
    las = 1,
    ylab = "Mean Overlap Percent (%)",
    xlab = "Pair Type",
    main = "Permutation-500 Null Control Mean"
  )

  xs <- as.numeric(bp)
  means <- as.numeric(mat)
  sds <- as.numeric(err)
  arrows(xs, pmax(0, means - sds), xs, means + sds, angle = 90, code = 3, length = 0.03, lwd = 1.2)

  axis(1, at = colMeans(bp), labels = colnames(mat), tick = FALSE, line = 0.5)
  legend("topright", legend = rownames(mat), fill = BIN_COLORS[rownames(mat)], bty = "n", title = "Overlap Bin")
}

plot_null_mean_gg <- function(summary_dt, out_svg) {
  suppressPackageStartupMessages(library(ggplot2))

  plot_dt <- copy(summary_dt)
  plot_dt[, pair_type := factor(pair_type, levels = PAIR_ORDER)]
  plot_dt[, overlap_bin := factor(overlap_bin, levels = BINS)]

  p <- ggplot(plot_dt, aes(x = pair_type, y = mean_percent, fill = overlap_bin)) +
    geom_col(position = position_dodge(width = 0.82), width = 0.72, color = "#222222", linewidth = 0.22) +
    geom_errorbar(
      aes(
        ymin = pmax(0, mean_percent - fifelse(is.na(sd_percent), 0, sd_percent)),
        ymax = mean_percent + fifelse(is.na(sd_percent), 0, sd_percent)
      ),
      width = 0.22,
      linewidth = 0.28,
      position = position_dodge(width = 0.82)
    ) +
    scale_fill_manual(values = BIN_COLORS, breaks = BINS, limits = BINS) +
    scale_x_discrete(limits = PAIR_ORDER) +
    labs(
      title = "Permutation-500 Null Control Mean",
      x = "Pair Type",
      y = "Mean Overlap Percent (%)",
      fill = "Overlap Bin"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "top"
    )

  ggsave(filename = out_svg, plot = p, width = 12, height = 7, units = "in", dpi = 300)
}

main <- function() {
  script_dir <- get_script_dir()
  project_dir <- script_dir

  default_base_dir <- "/picb/neurosys/chenrenrui/Mouse_ST/A molecularly defined and spatially resolved cell atlas of the whole mouse brain/"
  default_merge_region_root <- file.path(default_base_dir, "subclass_10K/neocortex_sample_subclass/merge_region")
  default_merge_samples_dir <- file.path(default_base_dir, "subclass_10K/neocortex_sample_subclass/merge_samples")
  default_overlap_dir <- file.path(default_base_dir, "subclass_10K/neocortex_sample_subclass/overlap_results")
  default_true_percent <- file.path(default_base_dir, "true_data_overlap_results/true_data_pair_overlap_percent.csv")
  default_true_csv_direct <- file.path(default_base_dir, "true_data_pair_overlap_percent_nonzero_pair.csv")
  default_mapping <- file.path(default_base_dir, "Merfish_brain_cell_type_subclass.txt")

  default_sample_table <- file.path(default_overlap_dir, "sample_pair_overlap_percent.csv")
  default_mean_csv <- file.path(default_overlap_dir, "sample500_null_control_mean.csv")
  default_mean_svg <- file.path(default_overlap_dir, "sample500_null_control_mean.svg")
  default_direct_svg <- file.path(default_overlap_dir, "sample_vs_true_overlap_plot.direct.svg")

  args <- parse_args(list(
    base_dir = default_base_dir,
    merge_region_root = default_merge_region_root,
    merge_samples_dir = default_merge_samples_dir,
    overlap_dir = default_overlap_dir,
    true_percent = default_true_percent,
    true_csv_direct = default_true_csv_direct,
    mapping = default_mapping,
    sample_table = default_sample_table,
    mean_csv = default_mean_csv,
    mean_svg = default_mean_svg,
    direct_svg = default_direct_svg,
    workers = "20",
    start_folder = "1",
    end_folder = "",
    python = "",
    rebuild_merge_samples = "true",
    run_overlap = "true",
    run_plot = "true",
    run_direct_plot = "false",
    build_mean_control = "true"
  ))

  base_dir_norm <- normalizePath(args$base_dir, mustWork = FALSE)
  default_base_dir_norm <- normalizePath(default_base_dir, mustWork = FALSE)
  if (!identical(base_dir_norm, default_base_dir_norm)) {
    if (identical(args$merge_region_root, default_merge_region_root)) {
      args$merge_region_root <- file.path(args$base_dir, "subclass_10K/neocortex_sample_subclass/merge_region")
    }
    if (identical(args$merge_samples_dir, default_merge_samples_dir)) {
      args$merge_samples_dir <- file.path(args$base_dir, "subclass_10K/neocortex_sample_subclass/merge_samples")
    }
    if (identical(args$overlap_dir, default_overlap_dir)) {
      args$overlap_dir <- file.path(args$base_dir, "subclass_10K/neocortex_sample_subclass/overlap_results")
    }
    if (identical(args$true_percent, default_true_percent)) {
      args$true_percent <- file.path(args$base_dir, "true_data_overlap_results/true_data_pair_overlap_percent.csv")
    }
    if (identical(args$true_csv_direct, default_true_csv_direct)) {
      args$true_csv_direct <- file.path(args$base_dir, "true_data_pair_overlap_percent_nonzero_pair.csv")
    }
    if (identical(args$mapping, default_mapping)) {
      args$mapping <- file.path(args$base_dir, "Merfish_brain_cell_type_subclass.txt")
    }
    if (identical(args$sample_table, default_sample_table)) {
      args$sample_table <- file.path(args$overlap_dir, "sample_pair_overlap_percent.csv")
    }
    if (identical(args$mean_csv, default_mean_csv)) {
      args$mean_csv <- file.path(args$overlap_dir, "sample500_null_control_mean.csv")
    }
    if (identical(args$mean_svg, default_mean_svg)) {
      args$mean_svg <- file.path(args$overlap_dir, "sample500_null_control_mean.svg")
    }
    if (identical(args$direct_svg, default_direct_svg)) {
      args$direct_svg <- file.path(args$overlap_dir, "sample_vs_true_overlap_plot.direct.svg")
    }
  }

  merge_region_root <- normalizePath(args$merge_region_root, mustWork = FALSE)
  merge_samples_dir <- normalizePath(args$merge_samples_dir, mustWork = FALSE)
  overlap_dir <- normalizePath(args$overlap_dir, mustWork = FALSE)
  true_percent <- normalizePath(args$true_percent, mustWork = FALSE)
  true_csv_direct <- normalizePath(args$true_csv_direct, mustWork = FALSE)
  mapping <- normalizePath(args$mapping, mustWork = FALSE)
  sample_table <- normalizePath(args$sample_table, mustWork = FALSE)
  mean_csv <- normalizePath(args$mean_csv, mustWork = FALSE)
  mean_svg <- normalizePath(args$mean_svg, mustWork = FALSE)
  direct_svg <- normalizePath(args$direct_svg, mustWork = FALSE)

  workers <- as.integer(args$workers)
  if (is.na(workers) || workers <= 0L) stop("--workers must be a positive integer.", call. = FALSE)

  start_folder <- as.integer(args$start_folder)
  if (is.na(start_folder) || start_folder <= 0L) stop("--start_folder must be >= 1.", call. = FALSE)

  end_folder <- NA_integer_
  if (nzchar(trimws(args$end_folder))) {
    end_folder <- as.integer(args$end_folder)
    if (is.na(end_folder) || end_folder < start_folder) {
      stop("--end_folder must be empty or >= --start_folder.", call. = FALSE)
    }
  }

  rebuild_merge_samples <- parse_bool(args$rebuild_merge_samples)
  run_overlap <- parse_bool(args$run_overlap)
  run_plot <- parse_bool(args$run_plot)
  run_direct_plot <- parse_bool(args$run_direct_plot)
  build_mean_control <- parse_bool(args$build_mean_control)

  python_exec <- resolve_python(args$python)

  compute_py <- file.path(project_dir, "compute_data1_overlap_500.py")
  plot_py <- file.path(project_dir, "plot_sample_overlap_permutation_500.py")
  direct_plot_py <- file.path(project_dir, "plot_sample_overlap_permutation_direct_500.py")

  if (run_overlap && !file.exists(compute_py)) {
    stop(sprintf("Missing script: %s", compute_py), call. = FALSE)
  }
  if (run_plot && !file.exists(plot_py)) {
    stop(sprintf("Missing script: %s", plot_py), call. = FALSE)
  }
  if (run_direct_plot && !file.exists(direct_plot_py)) {
    stop(sprintf("Missing script: %s", direct_plot_py), call. = FALSE)
  }

  message("[STEP 1] Prepare merge_samples intermediate folders")
  prep <- prepare_merge_samples(
    merge_region_root = merge_region_root,
    merge_samples_dir = merge_samples_dir,
    start_folder = start_folder,
    end_folder = end_folder,
    rebuild = rebuild_merge_samples
  )
  message(sprintf("[INFO] merge_samples ready: samples=%d, files_materialized=%d, files_skipped=%d",
                  prep$samples, prep$files_materialized, prep$files_skipped))

  dir.create(overlap_dir, recursive = TRUE, showWarnings = FALSE)

  if (run_overlap) {
    message("[STEP 2] Compute overlap summaries from merge_samples")
    overlap_args <- c(
      compute_py,
      "--data-dir", merge_samples_dir,
      "--output-dir", overlap_dir,
      "--workers", as.character(workers),
      "--start-folder", as.character(start_folder)
    )
    if (!is.na(end_folder)) {
      overlap_args <- c(overlap_args, "--end-folder", as.character(end_folder))
    }
    run_cmd(python_exec, overlap_args)
  } else {
    message("[STEP 2] Skipped overlap computation (--run_overlap=false)")
  }

  if (run_plot) {
    message("[STEP 3] Build sample distribution table and permutation plot")
    if (!file.exists(true_percent)) {
      stop(sprintf("Missing true-percent file: %s", true_percent), call. = FALSE)
    }
    if (!file.exists(mapping)) {
      stop(sprintf("Missing mapping file: %s", mapping), call. = FALSE)
    }

    plot_args <- c(
      plot_py,
      "--overlap-dir", overlap_dir,
      "--true-percent", true_percent,
      "--mapping", mapping,
      "--out-dir", overlap_dir,
      "--jobs", as.character(workers)
    )
    run_cmd(python_exec, plot_args)
  } else {
    message("[STEP 3] Skipped permutation plotting (--run_plot=false)")
  }

  if (run_direct_plot) {
    message("[STEP 4] Build direct plot from precomputed CSV")
    sample_csv_for_direct <- file.path(overlap_dir, "sample_pair_overlap_percent.csv")
    pvalue_csv_for_direct <- file.path(overlap_dir, "sample_vs_true_permutation_pvalues.csv")

    if (!file.exists(sample_csv_for_direct)) {
      stop(sprintf("Missing sample CSV for direct plot: %s", sample_csv_for_direct), call. = FALSE)
    }
    if (!file.exists(true_csv_direct)) {
      stop(sprintf("Missing true CSV for direct plot: %s", true_csv_direct), call. = FALSE)
    }

    direct_args <- c(
      direct_plot_py,
      "--sample-csv", sample_csv_for_direct,
      "--true-csv", true_csv_direct,
      "--out-svg", direct_svg
    )
    if (file.exists(pvalue_csv_for_direct)) {
      direct_args <- c(direct_args, "--pvalue-csv", pvalue_csv_for_direct)
    }
    run_cmd(python_exec, direct_args)
  } else {
    message("[STEP 4] Skipped direct plot (--run_direct_plot=false)")
  }

  if (build_mean_control) {
    message("[STEP 5] Build 500-sample mean control table/plot")
    if (!file.exists(sample_table)) {
      sample_table <- file.path(overlap_dir, "sample_pair_overlap_percent.csv")
    }
    if (!file.exists(sample_table)) {
      stop(sprintf("Missing sample table for mean control: %s", sample_table), call. = FALSE)
    }

    summary_dt <- compute_null_mean_table(sample_table, mean_csv)
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      plot_null_mean_gg(summary_dt, mean_svg)
    } else {
      message("[INFO] ggplot2 not found, using base plot fallback.")
      plot_null_mean_base(summary_dt, mean_svg)
    }

    n_samples <- unique(summary_dt$sample_n)
    n_samples <- n_samples[!is.na(n_samples)]
    n_samples <- max(n_samples, 0)
    message(sprintf("[INFO] mean control done: sample_n=%d", n_samples))
    message(sprintf("[INFO] mean CSV: %s", mean_csv))
    message(sprintf("[INFO] mean SVG: %s", mean_svg))
  } else {
    message("[STEP 5] Skipped mean control build (--build_mean_control=false)")
  }

  message("[DONE] sample500 pipeline finished.")
}

main()
