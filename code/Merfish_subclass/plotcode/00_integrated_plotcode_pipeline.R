#!/usr/bin/env Rscript

# ============================================================
# Integrated runner + figure organizer for the 12 plotcode scripts.
#
# Goal
# 1) Keep original 12 scripts unchanged.
# 2) Optionally run all scripts in one place.
# 3) Collect all generated figures.
# 4) Group similar analyses together.
# 5) Rename figures as continuous main index + panel letters:
#    figure_01a, figure_01b, ..., figure_02a, ...
#
# Inputs (shared by all analyses, placed in base_dir):
# - neuron-neuron-partner.csv (or txt in some scripts)
# - neuron-nonneuron-partner.csv (or txt in some scripts)
#
# Output:
# - plotcode/R_output_integrated/figures
# - plotcode/R_output_integrated/metadata
# ============================================================

options(stringsAsFactors = FALSE, scipen = 999)

get_script_dir <- function() {
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_all, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

patch_base_assignment <- function(lines, var_name, new_value) {
  # Replace e.g. base_dir <- "E:/zaw/2603"
  pattern <- paste0("(^[[:space:]]*", var_name, "[[:space:]]*<-[[:space:]]*)(\"[^\"]*\"|'[^']*')")
  replacement <- paste0("\\1\"", new_value, "\"")
  sub(pattern, replacement, lines, perl = TRUE)
}

slugify <- function(x, max_len = 70L) {
  y <- tolower(as.character(x))
  y <- gsub("[/\\\\]+", "_", y)
  y <- gsub("[^a-z0-9]+", "_", y)
  y <- gsub("^_+|_+$", "", y)
  if (!nzchar(y)) y <- "plot"
  if (nchar(y, type = "chars") > max_len) y <- substr(y, 1L, max_len)
  y
}

index_to_letters <- function(idx) {
  out <- character(length(idx))
  for (k in seq_along(idx)) {
    n <- as.integer(idx[k])
    if (is.na(n) || n <= 0L) {
      out[k] <- "a"
      next
    }
    s <- ""
    while (n > 0L) {
      n <- n - 1L
      s <- paste0(letters[(n %% 26L) + 1L], s)
      n <- n %/% 26L
    }
    out[k] <- s
  }
  out
}

infer_topic <- function(stem_text) {
  s <- tolower(stem_text)

  # Topic priority from specific to general.
  if (grepl("micro|glia", s)) return("Microglia Specific")
  if (grepl("priority|causal|precursor|tier", s)) return("Causal Priority")
  if (grepl("cross|dataset|gi_vs_in|gaba_in_glu|shared_source|targetclass", s)) return("Cross-Dataset Comparison")
  if (grepl("model|forest|regression|logit|adjusted|coef|contrast", s)) return("Adjusted Model Effects")
  if (grepl("layer|l23|l4|heatmap", s)) return("Layer Effect")
  if (grepl("direction|asym|flip|signed", s)) return("Directionality")
  if (grepl("subclass|subpair|pairclass|pair_class|motif", s)) return("Subclass/Pair Preference")
  if (grepl("jaccard|overlap|compact|cohesion|hosthost", s)) return("Overlap/Compactness")
  if (grepl("host|guest|multihost|multiplicity|cohost|fullhost", s)) return("Host-Guest Topology")
  if (grepl("contain|nest|exact|full", s)) return("Containment Pattern")
  if (grepl("size|geometry|scatter", s)) return("Geometry/Scale")
  "Other"
}

run_single_script <- function(script_path, base_dir_for_run) {
  code <- readLines(script_path, warn = FALSE, encoding = "UTF-8")

  # Most scripts use one of these names for project root.
  code <- patch_base_assignment(code, "base_dir", base_dir_for_run)
  code <- patch_base_assignment(code, "root_dir", base_dir_for_run)
  code <- patch_base_assignment(code, "project_root", base_dir_for_run)

  # Run in isolated env so rm(list = ls()) inside module scripts is harmless.
  env <- new.env(parent = globalenv())
  eval(parse(text = code), envir = env)
  invisible(TRUE)
}

collect_module_figures <- function(module_row, base_dir_for_outputs, keep_ext) {
  fig_dir <- file.path(base_dir_for_outputs, module_row[["fig_rel"]])
  if (!dir.exists(fig_dir)) {
    return(data.frame())
  }

  files <- list.files(fig_dir, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
  if (length(files) == 0L) {
    return(data.frame())
  }

  ext <- tolower(tools::file_ext(files))
  keep <- ext %in% keep_ext
  files <- files[keep]
  ext <- ext[keep]

  if (length(files) == 0L) {
    return(data.frame())
  }

  rel <- substring(files, nchar(fig_dir) + 2L)
  stem_rel <- tools::file_path_sans_ext(rel)

  data.frame(
    module_id = module_row[["module_id"]],
    script_file = module_row[["script_file"]],
    fig_dir = fig_dir,
    source_file = files,
    rel_path = rel,
    stem_rel = stem_rel,
    ext = ext,
    stringsAsFactors = FALSE
  )
}

main <- function() {
  script_dir <- get_script_dir()

  # -------------------------
  # User config
  # -------------------------
  base_dir <- "E:/zaw/2603"
  run_modules <- FALSE
  overwrite_existing <- TRUE
  keep_extensions <- c("png", "pdf", "svg")

  # If TRUE, only run selected modules.
  run_only <- integer(0)

  # -------------------------
  # Module registry (12 scripts)
  # -------------------------
  modules <- data.frame(
    module_id = 1:12,
    script_file = c(
      "01_ei_spatial_cluster_analysis.R",
      "02_neuron_interaction_analysis.R",
      "03_neuron_containment_extended_analysis.R",
      "04_neuron_containment_extension_analysis.R",
      "05_neuron_nonneuron_full_analysis.R",
      "06_neuron_nonneuron_interaction_analysis.R",
      "07_neuron_nonneuron_interaction_analysis_full.R",
      "08_cross_dataset_comparative_analysis.R",
      "09_gaba_in_glu_cross_interaction_full_analysis.R",
      "10_neuron_ie_special_cross_dataset_analysis_full.R",
      "11_causal_priority_analysis.R",
      "12_microglia_minimal_analysis.R"
    ),
    fig_rel = c(
      "R_output_01",
      "R_output_02/figures",
      "R_output_03/figures",
      "R_output_04/figures",
      "R_output_05/figures",
      "R_output_06/figures",
      "R_output_07/figures",
      "R_output_08/figures",
      "R_output_09/figures",
      "R_output_10/figures",
      "R_output_11/figures",
      "R_output_12/figures"
    ),
    stringsAsFactors = FALSE
  )

  if (length(run_only) > 0L) {
    modules <- modules[modules$module_id %in% run_only, , drop = FALSE]
  }

  output_root <- file.path(script_dir, "R_output_integrated")
  output_fig_dir <- file.path(output_root, "figures")
  output_meta_dir <- file.path(output_root, "metadata")
  dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_meta_dir, recursive = TRUE, showWarnings = FALSE)

  # -------------------------
  # 1) Optional: run modules
  # -------------------------
  run_log <- data.frame(
    module_id = modules$module_id,
    script_file = modules$script_file,
    status = "SKIPPED",
    message = "",
    stringsAsFactors = FALSE
  )

  if (isTRUE(run_modules)) {
    for (i in seq_len(nrow(modules))) {
      script_path <- file.path(script_dir, modules$script_file[i])
      if (!file.exists(script_path)) {
        run_log$status[i] <- "MISSING"
        run_log$message[i] <- paste("Script not found:", script_path)
        next
      }

      message(sprintf("[RUN] Module %02d: %s", modules$module_id[i], modules$script_file[i]))
      ok <- TRUE
      err_msg <- ""
      tryCatch(
        {
          run_single_script(script_path, base_dir_for_run = base_dir)
        },
        error = function(e) {
          ok <<- FALSE
          err_msg <<- conditionMessage(e)
        }
      )

      run_log$status[i] <- if (ok) "OK" else "FAILED"
      run_log$message[i] <- err_msg
    }
  }

  # -------------------------
  # 2) Collect all figures
  # -------------------------
  fig_list <- vector("list", nrow(modules))
  for (i in seq_len(nrow(modules))) {
    fig_list[[i]] <- collect_module_figures(
      module_row = modules[i, , drop = FALSE],
      base_dir_for_outputs = base_dir,
      keep_ext = keep_extensions
    )
  }
  fig_df <- do.call(rbind, fig_list)

  if (is.null(fig_df) || nrow(fig_df) == 0L) {
    write.csv(run_log, file.path(output_meta_dir, "module_run_log.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    stop("No figures found. Check base_dir and module output folders.")
  }

  # -------------------------
  # 3) Group same analyses together
  # -------------------------
  panel_df <- unique(fig_df[, c("module_id", "script_file", "stem_rel")])
  panel_df$topic <- vapply(panel_df$stem_rel, infer_topic, character(1))

  topic_order <- c(
    "Containment Pattern",
    "Directionality",
    "Layer Effect",
    "Subclass/Pair Preference",
    "Overlap/Compactness",
    "Host-Guest Topology",
    "Geometry/Scale",
    "Cross-Dataset Comparison",
    "Adjusted Model Effects",
    "Causal Priority",
    "Microglia Specific",
    "Other"
  )

  panel_df$topic_rank <- match(panel_df$topic, topic_order)
  panel_df$topic_rank[is.na(panel_df$topic_rank)] <- length(topic_order) + 1L
  panel_df <- panel_df[order(panel_df$topic_rank, panel_df$module_id, panel_df$stem_rel), , drop = FALSE]

  ordered_topics <- unique(panel_df$topic)
  topic_to_main <- setNames(seq_along(ordered_topics), ordered_topics)
  panel_df$main_idx <- unname(topic_to_main[panel_df$topic])

  panel_df$sub_idx <- 0L
  for (mid in unique(panel_df$main_idx)) {
    idx <- which(panel_df$main_idx == mid)
    panel_df$sub_idx[idx] <- seq_along(idx)
  }
  panel_df$sub_label <- index_to_letters(panel_df$sub_idx)
  panel_df$figure_label <- paste0(panel_df$main_idx, panel_df$sub_label)

  topic_slug <- vapply(panel_df$topic, slugify, character(1), max_len = 40L)
  stem_slug <- vapply(panel_df$stem_rel, slugify, character(1), max_len = 70L)
  panel_df$new_stem <- sprintf(
    "figure_%02d%s__%s__%s",
    panel_df$main_idx,
    panel_df$sub_label,
    topic_slug,
    stem_slug
  )
  panel_df$new_stem <- make.unique(panel_df$new_stem, sep = "_dup")

  # -------------------------
  # 4) Copy + rename figures
  # -------------------------
  panel_key <- paste(panel_df$module_id, panel_df$stem_rel, sep = "::")
  fig_df$panel_key <- paste(fig_df$module_id, fig_df$stem_rel, sep = "::")
  hit <- match(fig_df$panel_key, panel_key)

  fig_df$topic <- panel_df$topic[hit]
  fig_df$main_idx <- panel_df$main_idx[hit]
  fig_df$sub_label <- panel_df$sub_label[hit]
  fig_df$figure_label <- panel_df$figure_label[hit]
  fig_df$new_stem <- panel_df$new_stem[hit]
  fig_df$new_name <- paste0(fig_df$new_stem, ".", fig_df$ext)
  fig_df$output_file <- file.path(output_fig_dir, fig_df$new_name)

  copied <- logical(nrow(fig_df))
  for (i in seq_len(nrow(fig_df))) {
    copied[i] <- isTRUE(file.copy(fig_df$source_file[i], fig_df$output_file[i], overwrite = overwrite_existing))
  }
  fig_df$copied <- copied

  # -------------------------
  # 5) Save manifests
  # -------------------------
  fig_df <- fig_df[order(fig_df$main_idx, fig_df$sub_label, fig_df$module_id, fig_df$ext), , drop = FALSE]
  panel_df <- panel_df[order(panel_df$main_idx, panel_df$sub_label, panel_df$module_id), , drop = FALSE]

  write.csv(fig_df, file.path(output_meta_dir, "figure_files_manifest.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(panel_df, file.path(output_meta_dir, "figure_panels_manifest.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(run_log, file.path(output_meta_dir, "module_run_log.csv"), row.names = FALSE, fileEncoding = "UTF-8")

  topic_summary <- aggregate(
    x = list(n_panels = panel_df$figure_label),
    by = list(main_idx = panel_df$main_idx, topic = panel_df$topic),
    FUN = length
  )
  topic_summary <- topic_summary[order(topic_summary$main_idx), , drop = FALSE]
  write.csv(topic_summary, file.path(output_meta_dir, "topic_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")

  message("Integrated figure organization complete.")
  message("Output figures: ", normalizePath(output_fig_dir, mustWork = FALSE))
  message("Manifest folder: ", normalizePath(output_meta_dir, mustWork = FALSE))
}

main()

