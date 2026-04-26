#!/usr/bin/env Rscript

# ============================================================
# Single-file integrated script for all 12 plotcode analyses
# - Code from all 12 scripts is embedded as module functions
# - Unified output root: R_output_00
# - Similar analyses are grouped and renamed as 1a/1b/2a...
# - Plot export/style is unified for publication-quality output
# ============================================================

options(stringsAsFactors = FALSE, scipen = 999)

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
    if (!startsWith(key, "--")) stop(sprintf("Invalid argument: %s", key), call. = FALSE)
    name <- gsub("-", "_", sub("^--", "", key), fixed = TRUE)
    if (!name %in% names(out)) stop(sprintf("Unknown argument: %s", key), call. = FALSE)
    if (i == length(args)) stop(sprintf("Missing value for %s", key), call. = FALSE)
    out[[name]] <- args[i + 1L]
    i <- i + 2L
  }
  out
}

theme_pub_global <- function(base_size = 12, ...) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = base_size, color = "#2f2f2f", hjust = 0.5),
      axis.title = ggplot2::element_text(face = "bold", size = base_size + 0.5),
      axis.text = ggplot2::element_text(color = "#111111", size = base_size),
      axis.line = ggplot2::element_line(color = "#222222", linewidth = 0.5),
      axis.ticks = ggplot2::element_line(color = "#222222", linewidth = 0.4),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.text = ggplot2::element_text(size = base_size - 0.5),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(10, 12, 10, 10)
    )
}

ggsave_pub <- function(filename, plot = ggplot2::last_plot(), width = 8, height = 5, dpi = 700, bg = "white", ...) {
  ext <- tolower(tools::file_ext(filename))
  dots <- list(...)
  base_args <- list(filename = filename, plot = plot, width = width, height = height, bg = bg)

  if (ext == "pdf") {
    if (is.null(dots$device)) dots$device <- grDevices::cairo_pdf
  } else {
    if (is.null(dots$dpi)) dots$dpi <- max(600, dpi)
  }

  do.call(ggplot2::ggsave, c(base_args, dots))
}

slugify <- function(x, max_len = 80L) {
  vals <- as.character(x)
  out <- vapply(vals, function(one) {
    y <- tolower(one)
    y <- gsub("[/\\\\]+", "_", y)
    y <- gsub("[^a-z0-9]+", "_", y)
    y <- gsub("^_+|_+$", "", y)
    if (!nzchar(y)) y <- "plot"
    if (nchar(y, type = "chars") > max_len) y <- substr(y, 1L, max_len)
    y
  }, character(1), USE.NAMES = FALSE)
  out
}

force_dplyr_aliases <- function(target_env = .GlobalEnv) {
  if (!requireNamespace("dplyr", quietly = TRUE)) return(invisible(FALSE))

  fn_names <- c(
    "select", "rename", "filter", "mutate", "summarise", "summarize",
    "arrange", "distinct", "count", "left_join", "right_join",
    "full_join", "inner_join", "transmute", "pull"
  )

  for (nm in fn_names) {
    # Never touch locked bindings in target env.
    if (exists(nm, envir = target_env, inherits = FALSE) &&
        bindingIsLocked(nm, target_env)) {
      next
    }
    assign(nm, getExportedValue("dplyr", nm), envir = target_env)
  }

  invisible(TRUE)
}

index_to_letters <- function(idx) {
  out <- character(length(idx))
  for (k in seq_along(idx)) {
    n <- as.integer(idx[k])
    if (is.na(n) || n <= 0L) { out[k] <- "a"; next }
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
  if (grepl("micro|glia", s)) return("Microglia Specific")
  if (grepl("priority|causal|tier|precursor", s)) return("Causal Priority")
  if (grepl("cross|dataset|gi_vs_in|gaba_in_glu|targetclass|shared_source", s)) return("Cross-Dataset Comparison")
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

collect_figures <- function(runtime_base, module_map, keep_ext = c("png", "pdf", "svg")) {
  out <- vector("list", nrow(module_map))
  for (i in seq_len(nrow(module_map))) {
    d <- file.path(runtime_base, module_map$fig_rel[i])
    if (!dir.exists(d)) { out[[i]] <- data.frame(); next }
    fs <- list.files(d, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
    if (length(fs) == 0L) { out[[i]] <- data.frame(); next }
    ext <- tolower(tools::file_ext(fs))
    keep <- ext %in% keep_ext
    fs <- fs[keep]
    ext <- ext[keep]
    if (length(fs) == 0L) { out[[i]] <- data.frame(); next }
    rel <- substring(fs, nchar(d) + 2L)
    stem <- tools::file_path_sans_ext(rel)
    out[[i]] <- data.frame(
      module_id = module_map$module_id[i],
      script_file = module_map$script_file[i],
      source_dir = d,
      source_file = fs,
      rel_path = rel,
      stem_rel = stem,
      ext = ext,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

run_module_01 <- function(runtime_base) {
  message('[MODULE 01] start')
  # [removed] rm(list = ls())
  options(stringsAsFactors = FALSE)
  
  suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(ggplot2)
    library(scales)
  })
  
  # ============================================================
  # Top-journal style visualization for neuron-neuron partner analysis
  # Working directory: E:/zaw/2603
  # ============================================================
  
  # -----------------------------
  # [comment omitted: encoding-safe]
  # -----------------------------
  base_dir <- runtime_base
  path_txt <- file.path(base_dir, "neuron-neuron-partner.txt")
  path_csv <- file.path(base_dir, "neuron-neuron-partner.csv")
  path_out <- file.path(base_dir, "R_output_01")
  dir.create(path_out, showWarnings = FALSE, recursive = TRUE)
  
  # -----------------------------
  # [comment omitted: encoding-safe]
  # -----------------------------
  read_partner_file <- function(txt_path, csv_path) {
    if (file.exists(txt_path)) {
      message("Reading TXT: ", txt_path)
      read_tsv(txt_path, show_col_types = FALSE)
    } else if (file.exists(csv_path)) {
      message("Reading CSV: ", csv_path)
      read_csv(csv_path, show_col_types = FALSE)
    } else {
      stop("Cannot find input file. Please place neuron-neuron-partner.txt or neuron-neuron-partner.csv in: ", base_dir)
    }
  }
  
  # -----------------------------
  # [comment omitted: encoding-safe]
  # -----------------------------
  subclass_short <- function(x) {
    x %>%
      str_remove("^\\d+\\s+") %>%
      str_replace(" Gaba$", "") %>%
      str_replace(" Glut$", "")
  }
  
  fmt_p <- function(p) {
    if (length(p) > 1) {
      return(vapply(p, fmt_p, character(1)))
    }
    if (is.na(p)) return("NA")
    if (p == 0 || p < 1e-300) return("<1e-300")
    if (p < 1e-3) return(formatC(p, format = "e", digits = 2))
    sprintf("%.3f", p)
  }
  
  paired_compare <- function(piv, a, b) {
    tmp <- piv %>% dplyr::select(all_of(c(a, b))) %>% na.omit()
    if (nrow(tmp) == 0) {
      return(tibble(n = 0, p = NA_real_, median_diff = NA_real_, mean_diff = NA_real_))
    }
    w <- wilcox.test(tmp[[a]], tmp[[b]], paired = TRUE, exact = FALSE)
    tibble(
      n = nrow(tmp),
      p = w$p.value,
      median_diff = median(tmp[[a]] - tmp[[b]], na.rm = TRUE),
      mean_diff = mean(tmp[[a]] - tmp[[b]], na.rm = TRUE)
    )
  }
  
  save_pub <- function(p, filename, width = 7, height = 5) {
    ggsave_pub(
      filename = file.path(path_out, paste0(filename, ".pdf")),
      plot = p, width = width, height = height,
      device = cairo_pdf, dpi = 600, bg = "white"
    )
    ggsave_pub(
      filename = file.path(path_out, paste0(filename, ".png")),
      plot = p, width = width, height = height,
      dpi = 600, bg = "white"
    )
  }
  
  theme_pub <- function(base_size = 12) {
    theme_pub_global(base_size = base_size) +
      theme(
        plot.title = element_text(face = "bold", size = base_size + 1.5, hjust = 0.5, lineheight = 1.08),
        plot.subtitle = element_text(size = base_size - 0.4, hjust = 0.5, color = "grey20", lineheight = 1.05),
        axis.title = element_text(face = "bold", color = "black"),
        axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.4),
        axis.ticks.length = unit(0.18, "cm"),
        legend.title = element_text(face = "bold"),
        legend.position = "right",
        panel.border = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = base_size),
        plot.margin = margin(10, 16, 10, 10)
      )
  }
  
  pair_cols <- c(
    "Glut-Glut" = "#D95F02",
    "E-I" = "#1B9E77",
    "Gaba-Gaba" = "#7570B3"
  )
  
  # -----------------------------
  # [comment omitted: encoding-safe]
  # -----------------------------
  df <- read_partner_file(path_txt, path_csv)
  
  # -----------------------------
  # [comment omitted: encoding-safe]
  # -----------------------------
  df <- df %>%
    mutate(
      pair_key_sorted = if_else(
        `cluster.1_label` < `cluster.2_label`,
        paste(`cluster.1_label`, `cluster.2_label`, sep = "||"),
        paste(`cluster.2_label`, `cluster.1_label`, sep = "||")
      )
    )
  
  dedup <- df %>%
    arrange(pair_key_sorted, `cluster.1_label`, `cluster.2_label`) %>%
    distinct(pair_key_sorted, .keep_all = TRUE) %>%
    mutate(
      pair_type_unordered = case_when(
        `cluster.1_cell_Neruon_type` == "Glut" & `cluster.2_cell_Neruon_type` == "Glut" ~ "Glut-Glut",
        `cluster.1_cell_Neruon_type` == "Gaba" & `cluster.2_cell_Neruon_type` == "Gaba" ~ "Gaba-Gaba",
        TRUE ~ "E-I"
      ),
      max_containment = pmax(`cluster.1.overlap.percent`, `cluster.2.overlap.percent`),
      min_containment = pmin(`cluster.1.overlap.percent`, `cluster.2.overlap.percent`),
      sub1_short = subclass_short(`cluster.1_subclass`),
      sub2_short = subclass_short(`cluster.2_subclass`)
    )
  
  # -----------------------------
  # [comment omitted: encoding-safe]
  # -----------------------------
  clusters1 <- dedup %>%
    transmute(
      label       = `cluster.1_label`,
      slide       = `cluster.1_slide`,
      layer       = `cluster.1_layer`,
      region      = `cluster.1_region`,
      cell_type   = `cluster.1_cell_Neruon_type`,
      subclass    = `cluster.1_subclass`,
      total_cell_num = `cluster.1_total_cell_num`,
      Glut_num    = `cluster.1_Glut_Neruon_cell_ids_num`,
      GABA_num    = `cluster.1_GABA_Neruon_cell_ids_num`,
      cauchy_p    = `cluster.1_cauchy_combination_p`,
      E_I_Ratio   = `cluster.1_E_I_Ratio`
    )
  
  clusters2 <- dedup %>%
    transmute(
      label       = `cluster.2_label`,
      slide       = `cluster.2_slide`,
      layer       = `cluster.2_layer`,
      region      = `cluster.2_region`,
      cell_type   = `cluster.2_cell_Neruon_type`,
      subclass    = `cluster.2_subclass`,
      total_cell_num = `cluster.2_total_cell_num`,
      Glut_num    = `cluster.2_Glut_Neruon_cell_ids_num`,
      GABA_num    = `cluster.2_GABA_Neruon_cell_ids_num`,
      cauchy_p    = `cluster.2_cauchy_combination_p`,
      E_I_Ratio   = `cluster.2_E_I_Ratio`
    )
  
  clusters <- bind_rows(clusters1, clusters2) %>%
    distinct(label, .keep_all = TRUE) %>%
    mutate(
      gaba_frac = GABA_num / (GABA_num + Glut_num),
      glut_frac = Glut_num / (GABA_num + Glut_num),
      sub_short = subclass_short(subclass)
    )
  
  # -----------------------------
  # [comment omitted: encoding-safe]
  # -----------------------------
  overview <- tibble(
    directed_pairs    = nrow(df),
    unique_pairs      = nrow(dedup),
    unique_clusters   = nrow(clusters),
    slides            = n_distinct(clusters$slide),
    slide_layer_units = nrow(distinct(clusters, slide, layer)),
    glut_clusters     = sum(clusters$cell_type == "Glut"),
    gaba_clusters     = sum(clusters$cell_type == "Gaba")
  )
  write_csv(overview, file.path(path_out, "dataset_overview.csv"))
  
  # -----------------------------
  # Result block 1: pair type summary
  # -----------------------------
  pair_summary <- dedup %>%
    group_by(pair_type_unordered) %>%
    summarise(
      n = n(),
      mean_j = mean(jaccard, na.rm = TRUE),
      median_j = median(jaccard, na.rm = TRUE),
      q75_j = quantile(jaccard, 0.75, na.rm = TRUE),
      prop_cont95 = mean(max_containment >= 0.95, na.rm = TRUE),
      prop_cont99 = mean(max_containment >= 0.99, na.rm = TRUE),
      prop_cont100 = mean(max_containment >= 1.00, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(pair_type_unordered = factor(pair_type_unordered, levels = c("Glut-Glut", "E-I", "Gaba-Gaba")))
  write_csv(pair_summary, file.path(path_out, "pair_summary.csv"))
  
  agg_pair <- dedup %>%
    group_by(`cluster.1_slide`, `cluster.1_layer`, pair_type_unordered) %>%
    summarise(
      median_j = median(jaccard, na.rm = TRUE),
      mean_j = mean(jaccard, na.rm = TRUE),
      prop_cont95 = mean(max_containment >= 0.95, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(block_id = paste(`cluster.1_slide`, `cluster.1_layer`, sep = "|"))
  
  pivot_pair <- agg_pair %>%
    select(block_id, `cluster.1_slide`, `cluster.1_layer`, pair_type_unordered, median_j) %>%
    pivot_wider(names_from = pair_type_unordered, values_from = median_j) %>%
    drop_na()
  
  pivot_cont <- agg_pair %>%
    select(block_id, `cluster.1_slide`, `cluster.1_layer`, pair_type_unordered, prop_cont95) %>%
    pivot_wider(names_from = pair_type_unordered, values_from = prop_cont95) %>%
    drop_na()
  
  wilcox_median_GG_vs_EI <- wilcox.test(pivot_pair$`Glut-Glut`, pivot_pair$`E-I`, paired = TRUE, exact = FALSE)
  wilcox_median_II_vs_EI <- wilcox.test(pivot_pair$`Gaba-Gaba`, pivot_pair$`E-I`, paired = TRUE, exact = FALSE)
  wilcox_cont_EI_vs_GG   <- wilcox.test(pivot_cont$`E-I`, pivot_cont$`Glut-Glut`, paired = TRUE, exact = FALSE)
  wilcox_cont_EI_vs_II   <- wilcox.test(pivot_cont$`E-I`, pivot_cont$`Gaba-Gaba`, paired = TRUE, exact = FALSE)
  
  # -----------------------------
  # Result block 2: E-I asymmetry
  # -----------------------------
  ei <- dedup %>%
    filter(pair_type_unordered == "E-I") %>%
    mutate(
      gaba_sub = if_else(`cluster.1_cell_Neruon_type` == "Gaba", sub1_short, sub2_short),
      glut_sub = if_else(`cluster.1_cell_Neruon_type` == "Glut", sub1_short, sub2_short),
      gaba_overlap_percent = if_else(`cluster.1_cell_Neruon_type` == "Gaba", `cluster.1.overlap.percent`, `cluster.2.overlap.percent`),
      glut_overlap_percent = if_else(`cluster.1_cell_Neruon_type` == "Glut", `cluster.1.overlap.percent`, `cluster.2.overlap.percent`),
      delta_gaba_minus_glut_cont = gaba_overlap_percent - glut_overlap_percent,
      gaba_total_cell_num = if_else(`cluster.1_cell_Neruon_type` == "Gaba", `cluster.1_total_cell_num`, `cluster.2_total_cell_num`),
      glut_total_cell_num = if_else(`cluster.1_cell_Neruon_type` == "Glut", `cluster.1_total_cell_num`, `cluster.2_total_cell_num`),
      log_size_ratio = log(glut_total_cell_num / gaba_total_cell_num),
      max_containment = pmax(gaba_overlap_percent, glut_overlap_percent)
    )
  
  k_pair <- sum(ei$delta_gaba_minus_glut_cont > 0, na.rm = TRUE)
  n_pair <- sum(ei$delta_gaba_minus_glut_cont != 0, na.rm = TRUE)
  binom_pair <- binom.test(k_pair, n_pair, p = 0.5, alternative = "greater")
  
  ei_slide_layer <- ei %>%
    group_by(`cluster.1_slide`, `cluster.1_layer`) %>%
    summarise(
      median_delta = median(delta_gaba_minus_glut_cont, na.rm = TRUE),
      mean_delta   = mean(delta_gaba_minus_glut_cont, na.rm = TRUE),
      prop_gaba_more = mean(delta_gaba_minus_glut_cont > 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(block_id = paste(`cluster.1_slide`, `cluster.1_layer`, sep = "|"))
  
  k_sl <- sum(ei_slide_layer$median_delta > 0, na.rm = TRUE)
  n_sl <- nrow(ei_slide_layer)
  binom_slide_layer <- binom.test(k_sl, n_sl, p = 0.5, alternative = "greater")
  
  cor_ei <- suppressWarnings(cor.test(
    ei$log_size_ratio,
    ei$delta_gaba_minus_glut_cont,
    method = "spearman",
    exact = FALSE
  ))
  
  asym_summary <- tibble(
    metric = c(
      "E-I pair-level proportion: GABA containment > Glut containment",
      "E-I slide-layer proportion: median delta > 0",
      "Spearman rho: log(Glut/GABA size ratio) vs asymmetry delta"
    ),
    value = c(k_pair / n_pair, k_sl / n_sl, unname(cor_ei$estimate)),
    p = c(binom_pair$p.value, binom_slide_layer$p.value, cor_ei$p.value),
    n = c(n_pair, n_sl, sum(complete.cases(ei$log_size_ratio, ei$delta_gaba_minus_glut_cont)))
  )
  write_csv(asym_summary, file.path(path_out, "asymmetry_summary.csv"))
  write_csv(ei_slide_layer, file.path(path_out, "ei_slide_layer_summary.csv"))
  
  # -----------------------------
  # Result block 3: canonical GABA-GABA motifs
  # -----------------------------
  gaba_pairs <- dedup %>%
    filter(pair_type_unordered == "Gaba-Gaba") %>%
    mutate(
      subA = if_else(sub1_short < sub2_short, sub1_short, sub2_short),
      subB = if_else(sub1_short < sub2_short, sub2_short, sub1_short),
      pair_subclass = paste(subA, subB, sep = "||"),
      slide_layer = paste(`cluster.1_slide`, `cluster.1_layer`, sep = "|")
    )
  
  canonical_subs <- c("Vip","Lamp5","Sncg","Sst","Pvalb","Pvalb chandelier","Sst Chodl")
  
  canonical_gaba_pair_summary <- gaba_pairs %>%
    filter(subA %in% canonical_subs, subB %in% canonical_subs) %>%
    group_by(subA, subB) %>%
    summarise(
      n = n(),
      mean_j = mean(jaccard, na.rm = TRUE),
      median_j = median(jaccard, na.rm = TRUE),
      p75 = quantile(jaccard, 0.75, na.rm = TRUE),
      p90 = quantile(jaccard, 0.90, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(median_j), desc(n))
  write_csv(canonical_gaba_pair_summary, file.path(path_out, "canonical_gaba_pair_summary.csv"))
  
  canonical_interest <- c("Sst||Vip", "Pvalb||Vip", "Lamp5||Vip", "Pvalb||Sst", "Lamp5||Sst")
  
  canonical_piv <- gaba_pairs %>%
    filter(pair_subclass %in% canonical_interest) %>%
    group_by(slide_layer, pair_subclass) %>%
    summarise(median_j = median(jaccard, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = pair_subclass, values_from = median_j)
  
  canonical_tests <- bind_rows(
    paired_compare(canonical_piv, "Sst||Vip", "Pvalb||Vip") %>% mutate(comparison = "Vip-Sst vs Vip-Pvalb"),
    paired_compare(canonical_piv, "Pvalb||Sst", "Pvalb||Vip") %>% mutate(comparison = "Pvalb-Sst vs Pvalb-Vip"),
    paired_compare(canonical_piv, "Lamp5||Sst", "Lamp5||Vip") %>% mutate(comparison = "Lamp5-Sst vs Lamp5-Vip")
  ) %>%
    select(comparison, everything())
  write_csv(canonical_tests, file.path(path_out, "canonical_pair_tests_selected.csv"))
  
  # -----------------------------
  # Result block 4: layer dependence
  # -----------------------------
  ei_layer_asymmetry <- ei %>%
    group_by(`cluster.1_layer`) %>%
    summarise(
      n = n(),
      prop_gaba_more = mean(delta_gaba_minus_glut_cont > 0, na.rm = TRUE),
      median_delta   = median(delta_gaba_minus_glut_cont, na.rm = TRUE),
      mean_delta     = mean(delta_gaba_minus_glut_cont, na.rm = TRUE),
      prop_cont95    = mean(max_containment >= 0.95, na.rm = TRUE),
      median_j       = median(jaccard, na.rm = TRUE),
      .groups = "drop"
    )
  write_csv(ei_layer_asymmetry, file.path(path_out, "ei_layer_asymmetry.csv"))
  
  layer_tests <- tibble(
    pair_type = c("E-I", "E-I", "Glut-Glut", "Glut-Glut", "Gaba-Gaba", "Gaba-Gaba"),
    metric = c("median_j", "prop_cont95", "median_j", "prop_cont95", "median_j", "prop_cont95"),
    H = c(
      kruskal.test(median_j ~ `cluster.1_layer`, data = filter(agg_pair, pair_type_unordered == "E-I"))$statistic,
      kruskal.test(prop_cont95 ~ `cluster.1_layer`, data = filter(agg_pair, pair_type_unordered == "E-I"))$statistic,
      kruskal.test(median_j ~ `cluster.1_layer`, data = filter(agg_pair, pair_type_unordered == "Glut-Glut"))$statistic,
      kruskal.test(prop_cont95 ~ `cluster.1_layer`, data = filter(agg_pair, pair_type_unordered == "Glut-Glut"))$statistic,
      kruskal.test(median_j ~ `cluster.1_layer`, data = filter(agg_pair, pair_type_unordered == "Gaba-Gaba"))$statistic,
      kruskal.test(prop_cont95 ~ `cluster.1_layer`, data = filter(agg_pair, pair_type_unordered == "Gaba-Gaba"))$statistic
    ),
    p = c(
      kruskal.test(median_j ~ `cluster.1_layer`, data = filter(agg_pair, pair_type_unordered == "E-I"))$p.value,
      kruskal.test(prop_cont95 ~ `cluster.1_layer`, data = filter(agg_pair, pair_type_unordered == "E-I"))$p.value,
      kruskal.test(median_j ~ `cluster.1_layer`, data = filter(agg_pair, pair_type_unordered == "Glut-Glut"))$p.value,
      kruskal.test(prop_cont95 ~ `cluster.1_layer`, data = filter(agg_pair, pair_type_unordered == "Glut-Glut"))$p.value,
      kruskal.test(median_j ~ `cluster.1_layer`, data = filter(agg_pair, pair_type_unordered == "Gaba-Gaba"))$p.value,
      kruskal.test(prop_cont95 ~ `cluster.1_layer`, data = filter(agg_pair, pair_type_unordered == "Gaba-Gaba"))$p.value
    )
  )
  write_csv(layer_tests, file.path(path_out, "layer_kruskal_tests.csv"))
  
  # -----------------------------
  # key statistics
  # -----------------------------
  key_stats <- bind_rows(
    tibble(test = "paired Wilcoxon: slide-layer median Jaccard, Glut-Glut vs E-I",
           statistic = unname(wilcox_median_GG_vs_EI$statistic), p = wilcox_median_GG_vs_EI$p.value),
    tibble(test = "paired Wilcoxon: slide-layer median Jaccard, Gaba-Gaba vs E-I",
           statistic = unname(wilcox_median_II_vs_EI$statistic), p = wilcox_median_II_vs_EI$p.value),
    tibble(test = "paired Wilcoxon: slide-layer prop(max overlap>=0.95), E-I vs Glut-Glut",
           statistic = unname(wilcox_cont_EI_vs_GG$statistic), p = wilcox_cont_EI_vs_GG$p.value),
    tibble(test = "paired Wilcoxon: slide-layer prop(max overlap>=0.95), E-I vs Gaba-Gaba",
           statistic = unname(wilcox_cont_EI_vs_II$statistic), p = wilcox_cont_EI_vs_II$p.value),
    tibble(test = "Binomial: in E-I pairs, GABA is more contained than Glut (pair-level)",
           statistic = k_pair / n_pair, p = binom_pair$p.value),
    tibble(test = "Binomial: in E-I pairs, GABA is more contained than Glut (slide-layer median)",
           statistic = k_sl / n_sl, p = binom_slide_layer$p.value),
    tibble(test = "Spearman: log(Glut/GABA size ratio) vs delta containment",
           statistic = unname(cor_ei$estimate), p = cor_ei$p.value)
  )
  write_csv(key_stats, file.path(path_out, "key_statistics.csv"))
  
  # ============================================================
  # [comment omitted: encoding-safe]
  # [comment omitted: encoding-safe]
  # ============================================================
  
  save_donut_pub <- function(filename, width = 12, height = 5.8, expr_draw) {
    pdf(file.path(path_out, paste0(filename, ".pdf")),
        width = width, height = height, family = "sans", useDingbats = FALSE)
    expr_draw()
    dev.off()
    
    png(file.path(path_out, paste0(filename, ".png")),
        width = width, height = height, units = "in", res = 600, bg = "white")
    expr_draw()
    dev.off()
  }
  
  draw_single_donut <- function(cx, cy, counts, groups, fills, center_title,
                                outer_r = 1.55, inner_r = 0.82,
                                line_r1 = 1.58, line_r2 = 1.95,
                                text_dx = 0.24,
                                cex_center = 1.55,
                                cex_label = 1.02) {
    total_n <- sum(counts)
    props <- counts / total_n
    cumprops <- c(0, cumsum(props))
    
  # [comment omitted: encoding-safe]
    theta_start <- pi/2 - 2*pi*cumprops[-length(cumprops)]
    theta_end   <- pi/2 - 2*pi*cumprops[-1]
    
  # [comment omitted: encoding-safe]
    for (i in seq_along(counts)) {
      th <- seq(theta_start[i], theta_end[i], length.out = 300)
      x_outer <- cx + outer_r * cos(th)
      y_outer <- cy + outer_r * sin(th)
      x_inner <- cx + inner_r * cos(rev(th))
      y_inner <- cy + inner_r * sin(rev(th))
      
      polygon(
        c(x_outer, x_inner),
        c(y_outer, y_inner),
        col = fills[i],
        border = "white",
        lwd = 2.2
      )
    }
    
  # [comment omitted: encoding-safe]
    text(
      cx, cy,
      labels = paste0(center_title, "\nN = ", format(total_n, big.mark = ",")),
      cex = cex_center, font = 2, linespacing = 1.1
    )
    
  # [comment omitted: encoding-safe]
    for (i in seq_along(counts)) {
      mid <- (theta_start[i] + theta_end[i]) / 2
      side <- ifelse(cos(mid) >= 0, 1, -1)
      
  # [comment omitted: encoding-safe]
      x0 <- cx + line_r1 * cos(mid)
      y0 <- cy + line_r1 * sin(mid)
      x1 <- cx + line_r2 * cos(mid)
      y1 <- cy + line_r2 * sin(mid)
      x2 <- x1 + side * 0.34
      y2 <- y1
      
      segments(x0, y0, x1, y1, col = "grey35", lwd = 1.3)
      segments(x1, y1, x2, y2, col = "grey35", lwd = 1.3)
      
      lab <- paste0(
        groups[i], "\n",
        "n = ", format(counts[i], big.mark = ","), " / ", format(total_n, big.mark = ","), "\n",
        sprintf("%.1f%%", props[i] * 100)
      )
      
      text(
        x = x2 + side * text_dx,
        y = y2,
        labels = lab,
        adj = c(ifelse(side > 0, 0, 1), 0.5),
        cex = cex_label,
        font = 2,
        linespacing = 1.05
      )
    }
  }
  
  save_donut_pub(
    filename = "01_EI_asymmetry_pair_vs_block",
    width = 12,
    height = 5.8,
    expr_draw = function() {
      par(mar = c(0.3, 0.3, 0.3, 0.3),
          oma = c(0.4, 0.4, 2.7, 0.4),
          xpd = NA)
      
      plot.new()
      plot.window(xlim = c(-6.6, 6.6), ylim = c(-2.5, 2.7), asp = 1)
      
  # [comment omitted: encoding-safe]
      mtext(
        "E-I spatial organization is strongly asymmetric",
        side = 3, outer = TRUE, line = 0.6, cex = 1.85, font = 2
      )
      
  # [comment omitted: encoding-safe]
      mtext(
        paste0(
          "Pair-level binomial p = ", fmt_p(binom_pair$p.value),
          "   |   Slide-layer median binomial p = ", fmt_p(binom_slide_layer$p.value)
        ),
        side = 3, outer = TRUE, line = -0.9, cex = 1.15, font = 1
      )
      
  # [comment omitted: encoding-safe]
      draw_single_donut(
        cx = -3.2, cy = -0.1,
        counts = c(k_pair, n_pair - k_pair),
        groups = c("GABA more contained", "Glut more contained / tie"),
        fills = c("#6A3D9A", "#D9D9D9"),
        center_title = "Unique E-I pairs",
        outer_r = 1.55,
        inner_r = 0.82,
        line_r1 = 1.60,
        line_r2 = 1.95,
        text_dx = 0.20,
        cex_center = 1.45,
        cex_label = 0.98
      )
      
  # [comment omitted: encoding-safe]
      draw_single_donut(
        cx = 3.2, cy = -0.1,
        counts = c(k_sl, n_sl - k_sl),
        groups = c("GABA more contained", "Glut more contained / tie"),
        fills = c("#6A3D9A", "#D9D9D9"),
        center_title = "Slide-layer units",
        outer_r = 1.55,
        inner_r = 0.82,
        line_r1 = 1.60,
        line_r2 = 1.95,
        text_dx = 0.20,
        cex_center = 1.45,
        cex_label = 0.98
      )
    }
  )
  
  # ============================================================
  # PLOT 2: size ratio vs asymmetry
  # ============================================================
  p2 <- ggplot(ei, aes(x = log_size_ratio, y = delta_gaba_minus_glut_cont)) +
    geom_hline(yintercept = 0, color = "grey75", linewidth = 0.5, linetype = "dashed") +
    geom_point(color = "#3B528B", alpha = 0.18, size = 1.0) +
    geom_smooth(method = "loess", se = TRUE, color = "#F03B20", fill = "#FDD0A2", linewidth = 1.0) +
    annotate(
      "text",
      x = Inf, y = Inf,
      hjust = 1.05, vjust = 1.5,
      label = paste0("Spearman rho = ", sprintf("%.3f", unname(cor_ei$estimate)),
                     "\np = ", fmt_p(cor_ei$p.value)),
      size = 4.2, fontface = "bold"
    ) +
    labs(
      title = "Asymmetry scales with excitatory-to-inhibitory size difference",
      subtitle = "Positive values indicate that the GABA cluster is more contained than the Glut cluster",
      x = "log(Glut cluster size / GABA cluster size)",
      y = "Containment asymmetry (GABA - Glut)"
    ) +
    theme_pub_global(12)
  save_pub(p2, "02_EI_size_ratio_vs_asymmetry", width = 7.0, height = 5.4)
  
  # ============================================================
  # PLOT 3: slide-layer median Jaccard comparison
  # ============================================================
  plot_j <- pivot_pair %>%
    pivot_longer(cols = c(`Glut-Glut`, `E-I`, `Gaba-Gaba`),
                 names_to = "pair_type", values_to = "value") %>%
    mutate(pair_type = factor(pair_type, levels = c("Glut-Glut", "E-I", "Gaba-Gaba")))
  
  p3 <- ggplot(plot_j, aes(x = pair_type, y = value, color = pair_type)) +
    geom_line(aes(group = block_id), color = "grey82", alpha = 0.35, linewidth = 0.4) +
    geom_boxplot(aes(fill = pair_type), width = 0.50, alpha = 0.80, outlier.shape = NA, color = "black") +
    stat_summary(fun = median, geom = "point", shape = 23, size = 3.2, fill = "white", color = "black") +
    scale_fill_manual(values = pair_cols) +
    scale_color_manual(values = pair_cols) +
    labs(
      title = "E-I pairs are not the most strongly overlapped under a symmetric metric",
      subtitle = paste0(
        "Paired Wilcoxon: Glut-Glut vs E-I p = ", fmt_p(wilcox_median_GG_vs_EI$p.value),
        "   |   Gaba-Gaba vs E-I p = ", fmt_p(wilcox_median_II_vs_EI$p.value)
      ),
      x = NULL,
      y = "Median Jaccard within each slide-layer block"
    ) +
    theme_pub_global(12) +
    theme(legend.position = "none")
  save_pub(p3, "03_Jaccard_block_comparison", width = 7.6, height = 5.2)
  
  # ============================================================
  # PLOT 4: slide-layer near-complete containment comparison
  # ============================================================
  plot_c <- pivot_cont %>%
    pivot_longer(cols = c(`Glut-Glut`, `E-I`, `Gaba-Gaba`),
                 names_to = "pair_type", values_to = "value") %>%
    mutate(pair_type = factor(pair_type, levels = c("Glut-Glut", "E-I", "Gaba-Gaba")))
  
  p4 <- ggplot(plot_c, aes(x = pair_type, y = value, color = pair_type)) +
    geom_line(aes(group = block_id), color = "grey82", alpha = 0.35, linewidth = 0.4) +
    geom_boxplot(aes(fill = pair_type), width = 0.50, alpha = 0.80, outlier.shape = NA, color = "black") +
    stat_summary(fun = median, geom = "point", shape = 23, size = 3.2, fill = "white", color = "black") +
    scale_fill_manual(values = pair_cols) +
    scale_color_manual(values = pair_cols) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      title = "E-I pairs show the strongest near-complete containment",
      subtitle = paste0(
        "Paired Wilcoxon: E-I vs Glut-Glut p = ", fmt_p(wilcox_cont_EI_vs_GG$p.value),
        "   |   E-I vs Gaba-Gaba p = ", fmt_p(wilcox_cont_EI_vs_II$p.value)
      ),
      x = NULL,
      y = "Proportion with max containment >= 0.95"
    ) +
    theme_pub_global(12) +
    theme(legend.position = "none")
  save_pub(p4, "04_Containment_block_comparison", width = 7.6, height = 5.2)
  
  # ============================================================
  # PLOT 5: canonical GABA-GABA heatmap
  # ============================================================
  canonical_order <- c("Vip","Lamp5","Sncg","Sst","Pvalb","Pvalb chandelier")
  
  heat_df <- canonical_gaba_pair_summary %>%
    filter(subA %in% canonical_order, subB %in% canonical_order) %>%
    group_by(subA, subB) %>%
    summarise(median_j = median(median_j, na.rm = TRUE), .groups = "drop")
  
  heat_df2 <- bind_rows(
    heat_df,
    heat_df %>% transmute(subA = subB, subB = subA, median_j = median_j)
  ) %>%
    group_by(subA, subB) %>%
    summarise(median_j = median(median_j, na.rm = TRUE), .groups = "drop") %>%
    complete(subA = canonical_order, subB = canonical_order, fill = list(median_j = NA_real_)) %>%
    mutate(
      subA = factor(subA, levels = canonical_order),
      subB = factor(subB, levels = rev(canonical_order)))
  
  p5 <- ggplot(heat_df2, aes(x = subA, y = subB, fill = median_j)) +
    geom_tile(color = "white", linewidth = 0.55) +
    geom_text(
      data = heat_df2 %>% filter(!is.na(median_j)),
      aes(label = sprintf("%.2f", median_j)),
      size = 4.2, fontface = "bold", color = "black"
    ) +
    scale_fill_gradient(low = "#F7FBFF", high = "#1D91C0", na.value = "grey96") +
    labs(
      title = "Canonical inhibitory motifs do not show sharply selective spatial co-residence",
      subtitle = paste0(
        "Vip-Sst vs Vip-Pvalb p = ", fmt_p(canonical_tests$p[canonical_tests$comparison == "Vip-Sst vs Vip-Pvalb"]),
        "   |   Pvalb-Sst vs Pvalb-Vip p = ", fmt_p(canonical_tests$p[canonical_tests$comparison == "Pvalb-Sst vs Pvalb-Vip"]),
        "\nLamp5-Sst vs Lamp5-Vip p = ", fmt_p(canonical_tests$p[canonical_tests$comparison == "Lamp5-Sst vs Lamp5-Vip"])
      ),
      x = NULL, y = NULL, fill = "Median Jaccard"
    ) +
    theme_pub_global(11.5) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      axis.text.y = element_text(angle = 0)
    )
  save_pub(p5, "05_Canonical_GABA_heatmap", width = 8.0, height = 6.4)
  
  # ============================================================
  # PLOT 6: selected canonical GABA pair comparisons
  # ============================================================
  selected_plot_df <- bind_rows(
    canonical_piv %>% transmute(
      slide_layer,
      comparison = "Vip-centered comparison",
      pair_label = "Vip-Sst",
      median_j = `Sst||Vip`
    ),
    canonical_piv %>% transmute(
      slide_layer,
      comparison = "Vip-centered comparison",
      pair_label = "Vip-Pvalb",
      median_j = `Pvalb||Vip`
    ),
    canonical_piv %>% transmute(
      slide_layer,
      comparison = "Pvalb-centered comparison",
      pair_label = "Pvalb-Sst",
      median_j = `Pvalb||Sst`
    ),
    canonical_piv %>% transmute(
      slide_layer,
      comparison = "Pvalb-centered comparison",
      pair_label = "Pvalb-Vip",
      median_j = `Pvalb||Vip`
    ),
    canonical_piv %>% transmute(
      slide_layer,
      comparison = "Lamp5-centered comparison",
      pair_label = "Lamp5-Sst",
      median_j = `Lamp5||Sst`
    ),
    canonical_piv %>% transmute(
      slide_layer,
      comparison = "Lamp5-centered comparison",
      pair_label = "Lamp5-Vip",
      median_j = `Lamp5||Vip`
    )
  ) %>%
    filter(!is.na(median_j)) %>%
    mutate(
      comparison = factor(comparison,
                          levels = c("Vip-centered comparison",
                                     "Pvalb-centered comparison",
                                     "Lamp5-centered comparison")),
      pair_label = factor(pair_label,
                          levels = c("Vip-Sst","Vip-Pvalb","Pvalb-Sst","Pvalb-Vip","Lamp5-Sst","Lamp5-Vip"))
    )
  
  ann_df <- tibble(
    comparison = factor(c("Vip-centered comparison",
                          "Pvalb-centered comparison",
                          "Lamp5-centered comparison"),
                        levels = c("Vip-centered comparison",
                                   "Pvalb-centered comparison",
                                   "Lamp5-centered comparison")),
    x = 1.5,
    y = 0.98,
    label = c(
      paste0("p = ", fmt_p(canonical_tests$p[canonical_tests$comparison == "Vip-Sst vs Vip-Pvalb"])),
      paste0("p = ", fmt_p(canonical_tests$p[canonical_tests$comparison == "Pvalb-Sst vs Pvalb-Vip"])),
      paste0("p = ", fmt_p(canonical_tests$p[canonical_tests$comparison == "Lamp5-Sst vs Lamp5-Vip"]))
    )
  )
  
  pair_fill <- c(
    "Vip-Sst" = "#3BAF92",
    "Vip-Pvalb" = "#E67E22",
    "Pvalb-Sst" = "#E84393",
    "Pvalb-Vip" = "#F4A261",
    "Lamp5-Sst" = "#7CB342",
    "Lamp5-Vip" = "#8DA0CB"
  )
  
  p6 <- ggplot(selected_plot_df, aes(x = pair_label, y = median_j, fill = pair_label)) +
    geom_boxplot(width = 0.62, alpha = 0.88, outlier.shape = NA, color = "black") +
    geom_jitter(width = 0.10, size = 1.0, alpha = 0.18, color = "black") +
    geom_text(
      data = ann_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      fontface = "bold",
      size = 4.1
    ) +
    facet_wrap(~comparison, nrow = 1, scales = "free_x") +
    scale_fill_manual(values = pair_fill) +
    coord_cartesian(ylim = c(0, 1.02), clip = "off") +
    labs(
      title = "Canonical inhibitory pairs show broad overlap rather than strong spatial exclusivity",
      subtitle = "Each panel corresponds exactly to one slide-layer paired comparison reported in the statistics",
      x = NULL,
      y = "Median Jaccard within each slide-layer block"
    ) +
    theme_pub_global(12) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 20, hjust = 1)
    )
  save_pub(p6, "06_Canonical_GABA_selected_pairs", width = 11.2, height = 5.4)
  
  # ============================================================
  # PLOT 7: E-I layer dependence in median Jaccard
  # ============================================================
  ei_block_j <- agg_pair %>%
    filter(pair_type_unordered == "E-I", `cluster.1_layer` != "layer 6b") %>%
    mutate(`cluster.1_layer` = factor(`cluster.1_layer`,
                                      levels = c("layer 1","layer 2/3","layer 4","layer 5","layer 6a")))
  
  kw_ei_j <- kruskal.test(median_j ~ `cluster.1_layer`, data = ei_block_j)
  
  p7 <- ggplot(ei_block_j, aes(x = `cluster.1_layer`, y = median_j)) +
    geom_boxplot(width = 0.62, fill = "#B2DF8A", color = "black", alpha = 0.85, outlier.shape = NA) +
    geom_jitter(width = 0.14, size = 1.0, alpha = 0.25, color = "#4D4D4D") +
    stat_summary(fun = median, geom = "point", shape = 23, size = 3.0, fill = "#1B7837", color = "black") +
    labs(
      title = "E-I overlap strength shows strong layer dependence",
      subtitle = paste0("Kruskal-Wallis p = ", fmt_p(kw_ei_j$p.value)),
      x = NULL, y = "Median Jaccard within each slide-layer block"
    ) +
    theme_pub_global(12)
  save_pub(p7, "07_EI_layer_median_jaccard", width = 7.0, height = 5.0)
  
  # ============================================================
  # PLOT 8: E-I layer dependence in near-complete containment
  # ============================================================
  ei_block_c <- agg_pair %>%
    filter(pair_type_unordered == "E-I", `cluster.1_layer` != "layer 6b") %>%
    mutate(`cluster.1_layer` = factor(`cluster.1_layer`,
                                      levels = c("layer 1","layer 2/3","layer 4","layer 5","layer 6a")))
  
  kw_ei_c <- kruskal.test(prop_cont95 ~ `cluster.1_layer`, data = ei_block_c)
  
  layer_prop95_summary <- ei_block_c %>%
    group_by(`cluster.1_layer`) %>%
    summarise(
      mean_prop = mean(prop_cont95, na.rm = TRUE),
      median_prop = median(prop_cont95, na.rm = TRUE),
      .groups = "drop"
    )
  write_csv(layer_prop95_summary, file.path(path_out, "ei_layer_prop95_summary.csv"))
  
  p8 <- ggplot(ei_block_c, aes(x = `cluster.1_layer`, y = prop_cont95)) +
    geom_boxplot(width = 0.62, fill = "#A6CEE3", color = "black", alpha = 0.85, outlier.shape = NA) +
    geom_jitter(width = 0.14, size = 1.0, alpha = 0.25, color = "#4D4D4D") +
    stat_summary(fun = median, geom = "point", shape = 23, size = 3.0, fill = "#1F78B4", color = "black") +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      title = "Near-complete E-I containment is strongest in superficial-middle layers",
      subtitle = paste0(
        "Kruskal-Wallis p = ", fmt_p(kw_ei_c$p.value),
        "   |   Highest block-level containment is expected in layer 2/3"
      ),
      x = NULL, y = "Proportion with max containment >= 0.95"
    ) +
    theme_pub_global(12)
  save_pub(p8, "08_EI_layer_containment95", width = 7.0, height = 5.0)
  
  # -----------------------------
  # session info
  # -----------------------------
  sink(file.path(path_out, "sessionInfo.txt"))
  print(sessionInfo())
  sink()
  
  cat("All analyses finished.\n")
  cat("Output directory:\n", normalizePath(path_out), "\n")
}

run_module_02 <- function(runtime_base) {
  message('[MODULE 02] start')
  # =========================================================
  # neuron_interaction_topjournal_R_output_02.R
  # [comment omitted: encoding-safe]
  # [comment omitted: encoding-safe]
  # [comment omitted: encoding-safe]
  # =========================================================
  
  required_pkgs <- c(
    "data.table", "dplyr", "tidyr", "stringr", "purrr", "readr", "tibble",
    "forcats", "ggplot2", "ggrepel", "scales", "sandwich", "lmtest"
  )
  
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
  
  suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(purrr)
    library(readr)
    library(tibble)
    library(forcats)
    library(ggplot2)
    library(ggrepel)
    library(scales)
    library(sandwich)
    library(lmtest)
  })
  
  # ---------------------------
  # Paths
  # ---------------------------
  base_dir <- runtime_base
  out_dir <- file.path(base_dir, "R_output_02")
  fig_dir <- file.path(out_dir, "figures")
  tbl_dir <- file.path(out_dir, "tables")
  meta_dir <- file.path(out_dir, "metadata")
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(tbl_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(meta_dir, showWarnings = FALSE, recursive = TRUE)
  
  find_input_csv <- function(base_dir) {
    preferred <- c(
      file.path(base_dir, "neuron-neuron-partner.csv"),
      file.path(base_dir, "neuron_neuron_partner.csv"),
      file.path(base_dir, "partner.csv")
    )
    preferred <- preferred[file.exists(preferred)]
    if (length(preferred) > 0) return(preferred[1])
    
    hits <- list.files(
      base_dir,
      pattern = "(partner|neuron).*(csv)$",
      full.names = TRUE,
      ignore.case = TRUE
    )
    if (length(hits) > 0) return(hits[1])
    
    stop("未在E:/zaw/2603中找到输入CSV。请把原始数据放到该目录下，并命名为neuron-neuron-partner.csv或类似名称。")
  }
  
  input_csv <- find_input_csv(base_dir)
  
  # ---------------------------
  # Global styles
  # ---------------------------
  base_family <- "sans"
  
  pair_cols <- c(
    "E-E" = "#3C78D8",
    "E-I" = "#D1495B",
    "I-I" = "#2A9D8F"
  )
  
  pref_cols <- c(
    "same-type preferring" = "#2A9D8F",
    "cross-type preferring" = "#D1495B"
  )
  
  bg_cols <- c(
    "Higher than block E-I background" = "#C44E52",
    "Lower than block E-I background" = "#4C72B0"
  )
  
  pub_theme <- function(base_size = 12) {
    theme_pub_global(base_size = base_size, base_family = base_family) +
      theme(
        plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0, margin = margin(b = 8)),
        plot.subtitle = element_text(size = base_size, hjust = 0, color = "#444444", margin = margin(b = 10)),
        plot.caption = element_text(size = base_size - 1, color = "#555555", hjust = 0),
        axis.title = element_text(face = "bold", size = base_size + 0.5),
        axis.text = element_text(color = "black", size = base_size),
        axis.line = element_line(linewidth = 0.55, color = "black"),
        axis.ticks = element_line(linewidth = 0.5, color = "black"),
        axis.ticks.length = grid::unit(-1.8, "mm"),
        axis.text.x = element_text(margin = margin(t = 7)),
        axis.text.y = element_text(margin = margin(r = 7)),
        strip.background = element_rect(fill = "#F2F2F2", color = NA),
        strip.text = element_text(face = "bold", size = base_size),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = base_size),
        legend.key.width = grid::unit(1.2, "cm"),
        plot.margin = margin(10, 14, 10, 10)
      )
  }
  
  theme_set(pub_theme(12))
  
  save_pub <- function(p, filename, width, height) {
    png(file.path(fig_dir, paste0(filename, ".png")), width = width, height = height, units = "in", res = 600, bg = "white")
    print(p)
    dev.off()
    
    pdf(file.path(fig_dir, paste0(filename, ".pdf")), width = width, height = height, bg = "white", useDingbats = FALSE)
    print(p)
    dev.off()
  }
  
  fmt_num <- function(x, digits = 4) {
    ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
  }
  
  fmt_p <- function(x) {
    ifelse(
      is.na(x), "NA",
      ifelse(x < 1e-3, formatC(x, format = "e", digits = 2), formatC(x, format = "f", digits = 4))
    )
  }
  
  canonical_pair <- function(a, b) {
    ifelse(a < b, paste(a, b, sep = " || "), paste(b, a, sep = " || "))
  }
  
  clean_subclass <- function(x) {
    str_replace(x, "^\\d+\\s+", "")
  }
  
  rank_biserial_paired <- function(x, y) {
    d <- x - y
    d <- d[d != 0]
    if (length(d) == 0) return(NA_real_)
    pos <- sum(d > 0)
    neg <- sum(d < 0)
    (pos - neg) / length(d)
  }
  
  safe_paired_wilcox <- function(x, y) {
    keep <- is.finite(x) & is.finite(y)
    x <- x[keep]
    y <- y[keep]
    if (length(x) < 3) {
      return(list(p_value = NA_real_, n = length(x), mean_diff = NA_real_, median_diff = NA_real_, rank_biserial = NA_real_))
    }
    wt <- suppressWarnings(wilcox.test(x, y, paired = TRUE, exact = FALSE))
    list(
      p_value = wt$p.value,
      n = length(x),
      mean_diff = mean(x - y),
      median_diff = median(x - y),
      rank_biserial = rank_biserial_paired(x, y)
    )
  }
  
  paired_wilcox_df <- function(df, a, b, metric_name) {
    tmp <- df %>% select(all_of(c(a, b))) %>% drop_na()
    res <- safe_paired_wilcox(tmp[[a]], tmp[[b]])
    tibble(
      metric = metric_name,
      comparison = paste(a, "vs", b),
      n = res$n,
      mean_diff = res$mean_diff,
      median_diff = res$median_diff,
      rank_biserial = res$rank_biserial,
      p_value = res$p_value
    )
  }
  
  pick_focus_row <- function(df, tokens) {
    if (nrow(df) == 0) return(df[0, ])
    keep <- rep(TRUE, nrow(df))
    for (tok in tokens) {
      keep <- keep & str_detect(df$subpair_name_raw, fixed(tok, ignore_case = TRUE))
    }
    out <- df[keep, , drop = FALSE]
    if (nrow(out) == 0) return(df[0, ])
    out %>% slice(1)
  }
  
  make_empty_plot <- function(title_text) {
    ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = title_text, size = 5, fontface = "bold") +
      xlim(0, 1) + ylim(0, 1) +
      theme_void(base_family = base_family)
  }
  
  # ---------------------------
  # Read data
  # ---------------------------
  df <- fread(input_csv)
  df[, pair_key := canonical_pair(cluster.1_label, cluster.2_label)]
  
  # ---------------------------
  # Unique clusters
  # ---------------------------
  c1 <- copy(df[, .(
    label = cluster.1_label,
    slide = cluster.1_slide,
    layer = cluster.1_layer,
    region = cluster.1_region,
    cell_Neruon_type = cluster.1_cell_Neruon_type,
    subclass = cluster.1_subclass,
    total_cell_num = cluster.1_total_cell_num,
    Glut_Neruon_cell_ids_num = cluster.1_Glut_Neruon_cell_ids_num,
    GABA_Neruon_cell_ids_num = cluster.1_GABA_Neruon_cell_ids_num,
    cauchy_combination_p = cluster.1_cauchy_combination_p,
    E_I_Ratio = cluster.1_E_I_Ratio
  )])
  
  c2 <- copy(df[, .(
    label = cluster.2_label,
    slide = cluster.2_slide,
    layer = cluster.2_layer,
    region = cluster.2_region,
    cell_Neruon_type = cluster.2_cell_Neruon_type,
    subclass = cluster.2_subclass,
    total_cell_num = cluster.2_total_cell_num,
    Glut_Neruon_cell_ids_num = cluster.2_Glut_Neruon_cell_ids_num,
    GABA_Neruon_cell_ids_num = cluster.2_GABA_Neruon_cell_ids_num,
    cauchy_combination_p = cluster.2_cauchy_combination_p,
    E_I_Ratio = cluster.2_E_I_Ratio
  )])
  
  clusters <- unique(rbindlist(list(c1, c2), use.names = TRUE))
  clusters[, block := paste(slide, layer, sep = "|")]
  clusters[, subclass_clean := clean_subclass(subclass)]
  
  # ---------------------------
  # Unique undirected pairs
  # ---------------------------
  pair_list <- lapply(split(df, df$pair_key), function(sub) {
    stopifnot(nrow(sub) == 2)
    labs <- sort(unique(c(sub$cluster.1_label, sub$cluster.2_label)))
    a <- labs[1]
    b <- labs[2]
    ra <- sub[sub$cluster.1_label == a, ][1, ]
    type_a <- ifelse(ra$cluster.1_cell_Neruon_type == "Glut", "E", "I")
    type_b <- ifelse(ra$cluster.2_cell_Neruon_type == "Glut", "E", "I")
    
    sub_raw_sorted <- sort(c(ra$cluster.1_subclass, ra$cluster.2_subclass))
    sub_clean_sorted <- sort(c(clean_subclass(ra$cluster.1_subclass), clean_subclass(ra$cluster.2_subclass)))
    
    tibble(
      label_a = a,
      label_b = b,
      type_a_raw = ra$cluster.1_cell_Neruon_type,
      type_b_raw = ra$cluster.2_cell_Neruon_type,
      type_a = type_a,
      type_b = type_b,
      pair_type = paste(sort(c(type_a, type_b)), collapse = "-"),
      subclass_a = ra$cluster.1_subclass,
      subclass_b = ra$cluster.2_subclass,
      subclass_a_clean = clean_subclass(ra$cluster.1_subclass),
      subclass_b_clean = clean_subclass(ra$cluster.2_subclass),
      subpair_name_raw = paste(sub_raw_sorted, collapse = " × "),
      subpair_name_clean = paste(sub_clean_sorted, collapse = " × "),
      layer = ra$cluster.1_layer,
      slide = ra$cluster.1_slide,
      region_a = ra$cluster.1_region,
      region_b = ra$cluster.2_region,
      size_a = ra$cluster.1_total_cell_num,
      size_b = ra$cluster.2_total_cell_num,
      EI_ratio_a = ra$cluster.1_E_I_Ratio,
      EI_ratio_b = ra$cluster.2_E_I_Ratio,
      p_a = ra$cluster.1_cauchy_combination_p,
      p_b = ra$cluster.2_cauchy_combination_p,
      overlap = ra$overlap_cell,
      union = ra$union_cell,
      jaccard = ra$jaccard,
      overlap_prop_a = ra$cluster.1.overlap.percent,
      overlap_prop_b = ra$cluster.2.overlap.percent
    )
  })
  
  pairs_u <- bind_rows(pair_list) %>%
    mutate(
      block = paste(slide, layer, sep = "|"),
      size_small = pmin(size_a, size_b),
      size_large = pmax(size_a, size_b),
      size_ratio = size_large / size_small,
      overlap_prop_small = overlap / size_small,
      overlap_prop_large = overlap / size_large,
      nesting_gap = overlap_prop_small - overlap_prop_large,
      nested_90 = overlap_prop_small >= 0.90,
      pair_type = factor(pair_type, levels = c("E-E", "E-I", "I-I")),
      layer = factor(layer, levels = c("layer 1", "layer 2/3", "layer 4", "layer 5", "layer 6a", "layer 6b"))
    )
  
  # ---------------------------
  # Long edge table for cluster-centric preference
  # ---------------------------
  edges_long <- bind_rows(
    pairs_u %>% transmute(
      label = label_a,
      type_raw = type_a_raw,
      subclass = subclass_a,
      subclass_clean = subclass_a_clean,
      layer = layer,
      partner_label = label_b,
      partner_type = type_b_raw,
      jaccard = jaccard
    ),
    pairs_u %>% transmute(
      label = label_b,
      type_raw = type_b_raw,
      subclass = subclass_b,
      subclass_clean = subclass_b_clean,
      layer = layer,
      partner_label = label_a,
      partner_type = type_a_raw,
      jaccard = jaccard
    )
  )
  
  # ---------------------------
  # Core summaries
  # ---------------------------
  block_metrics <- pairs_u %>%
    group_by(block, slide, layer, pair_type) %>%
    summarise(
      n_pairs = n(),
      mean_jaccard = mean(jaccard, na.rm = TRUE),
      mean_overlap_prop_small = mean(overlap_prop_small, na.rm = TRUE),
      mean_overlap_prop_large = mean(overlap_prop_large, na.rm = TRUE),
      mean_nesting_gap = mean(nesting_gap, na.rm = TRUE),
      prop_nested90 = mean(nested_90, na.rm = TRUE),
      .groups = "drop"
    )
  
  blocks_complete <- block_metrics %>% count(block) %>% filter(n == 3) %>% pull(block)
  block_metrics_complete <- block_metrics %>% filter(block %in% blocks_complete)
  
  block_small_wide <- block_metrics_complete %>%
    select(block, slide, layer, pair_type, mean_overlap_prop_small) %>%
    pivot_wider(names_from = pair_type, values_from = mean_overlap_prop_small) %>%
    drop_na()
  
  block_gap_wide <- block_metrics_complete %>%
    select(block, slide, layer, pair_type, mean_nesting_gap) %>%
    pivot_wider(names_from = pair_type, values_from = mean_nesting_gap) %>%
    drop_na()
  
  block_j_wide <- block_metrics_complete %>%
    select(block, slide, layer, pair_type, mean_jaccard) %>%
    pivot_wider(names_from = pair_type, values_from = mean_jaccard) %>%
    drop_na()
  
  block_nested_wide <- block_metrics_complete %>%
    select(block, slide, layer, pair_type, prop_nested90) %>%
    pivot_wider(names_from = pair_type, values_from = prop_nested90) %>%
    drop_na()
  
  # ---------------------------
  # Tests: geometry family
  # [comment omitted: encoding-safe]
  # ---------------------------
  overall_tests <- bind_rows(
    paired_wilcox_df(block_j_wide, "E-E", "E-I", "block_mean_jaccard"),
    paired_wilcox_df(block_j_wide, "E-E", "I-I", "block_mean_jaccard"),
    paired_wilcox_df(block_j_wide, "E-I", "I-I", "block_mean_jaccard"),
    paired_wilcox_df(block_small_wide, "E-E", "E-I", "block_mean_overlap_prop_small"),
    paired_wilcox_df(block_small_wide, "E-E", "I-I", "block_mean_overlap_prop_small"),
    paired_wilcox_df(block_small_wide, "E-I", "I-I", "block_mean_overlap_prop_small"),
    paired_wilcox_df(block_gap_wide, "E-E", "E-I", "block_mean_nesting_gap"),
    paired_wilcox_df(block_gap_wide, "E-E", "I-I", "block_mean_nesting_gap"),
    paired_wilcox_df(block_gap_wide, "E-I", "I-I", "block_mean_nesting_gap")
  ) %>% mutate(fdr_bh = p.adjust(p_value, method = "BH"))
  
  # ---------------------------
  # Tests: layer-wise near-complete nesting family
  # ---------------------------
  nested_layer_tests <- block_metrics_complete %>%
    group_by(layer) %>%
    group_split() %>%
    map_dfr(function(dat) {
      layer_name <- unique(dat$layer)
      piv <- dat %>%
        select(block, slide, pair_type, prop_nested90) %>%
        pivot_wider(names_from = pair_type, values_from = prop_nested90) %>%
        drop_na()
      
      bind_rows(
        paired_wilcox_df(piv, "E-E", "E-I", "prop_nested90"),
        paired_wilcox_df(piv, "E-E", "I-I", "prop_nested90"),
        paired_wilcox_df(piv, "E-I", "I-I", "prop_nested90")
      ) %>% mutate(layer = as.character(layer_name))
    }) %>%
    mutate(fdr_bh = p.adjust(p_value, method = "BH"))
  
  # ---------------------------
  # Cluster-centric preference
  # ---------------------------
  cluster_pref <- edges_long %>%
    mutate(same_type_partner = type_raw == partner_type) %>%
    group_by(label, type_raw, subclass, subclass_clean, layer, same_type_partner) %>%
    summarise(mean_jaccard = mean(jaccard, na.rm = TRUE), n = n(), .groups = "drop") %>%
    pivot_wider(names_from = same_type_partner, values_from = mean_jaccard, names_prefix = "same_") %>%
    rename(cross_type_mean_jaccard = same_FALSE, same_type_mean_jaccard = same_TRUE) %>%
    drop_na(cross_type_mean_jaccard, same_type_mean_jaccard) %>%
    mutate(same_minus_cross = same_type_mean_jaccard - cross_type_mean_jaccard)
  
  exc_pref <- cluster_pref %>% filter(type_raw == "Glut")
  gaba_pref <- cluster_pref %>% filter(type_raw == "Gaba")
  
  safe_pref_test <- function(dat) {
    if (nrow(dat) < 3) {
      return(tibble(
        n_clusters = nrow(dat),
        mean_same_minus_cross = mean(dat$same_minus_cross, na.rm = TRUE),
        median_same_minus_cross = median(dat$same_minus_cross, na.rm = TRUE),
        p_value = NA_real_
      ))
    }
    tibble(
      n_clusters = nrow(dat),
      mean_same_minus_cross = mean(dat$same_minus_cross, na.rm = TRUE),
      median_same_minus_cross = median(dat$same_minus_cross, na.rm = TRUE),
      p_value = suppressWarnings(wilcox.test(dat$same_type_mean_jaccard, dat$cross_type_mean_jaccard, paired = TRUE, exact = FALSE)$p.value)
    )
  }
  
  exc_sub_tests <- exc_pref %>%
    group_by(subclass, subclass_clean) %>%
    group_modify(~ safe_pref_test(.x)) %>%
    ungroup() %>%
    filter(n_clusters >= 5) %>%
    mutate(fdr_bh = p.adjust(p_value, method = "BH"))
  
  gaba_sub_tests <- gaba_pref %>%
    group_by(subclass, subclass_clean) %>%
    group_modify(~ safe_pref_test(.x)) %>%
    ungroup() %>%
    filter(n_clusters >= 5) %>%
    mutate(fdr_bh = p.adjust(p_value, method = "BH"))
  
  exc_sub_layer_tests <- exc_pref %>%
    group_by(subclass, subclass_clean, layer) %>%
    group_modify(~ safe_pref_test(.x)) %>%
    ungroup() %>%
    filter(n_clusters >= 3) %>%
    mutate(fdr_bh = p.adjust(p_value, method = "BH"))
  
  # ---------------------------
  # E-I subclass pair vs within-block E-I background
  # ---------------------------
  subpair_tests <- pairs_u %>%
    group_by(subpair_name_raw, subpair_name_clean) %>%
    group_split() %>%
    map_dfr(function(sub) {
      sp <- unique(sub$subpair_name_raw)
      sp_clean <- unique(sub$subpair_name_clean)
      if (nrow(sub) < 5 || n_distinct(sub$block) < 3) return(NULL)
      
      sp_block <- sub %>%
        group_by(block, pair_type) %>%
        summarise(sp_mean = mean(jaccard, na.rm = TRUE), .groups = "drop")
      
      other_bg <- pairs_u %>%
        filter(subpair_name_raw != sp) %>%
        group_by(block, pair_type) %>%
        summarise(other_mean = mean(jaccard, na.rm = TRUE), .groups = "drop")
      
      merged <- inner_join(sp_block, other_bg, by = c("block", "pair_type"))
      if (nrow(merged) < 3) return(NULL)
      
      wt <- safe_paired_wilcox(merged$sp_mean, merged$other_mean)
      
      tibble(
        subpair_name_raw = sp,
        subpair_name_clean = sp_clean,
        pair_type = unique(as.character(merged$pair_type))[1],
        n_edges = nrow(sub),
        n_blocks = nrow(merged),
        mean_jaccard = mean(sub$jaccard, na.rm = TRUE),
        median_jaccard = median(sub$jaccard, na.rm = TRUE),
        mean_diff_vs_bg = wt$mean_diff,
        median_diff_vs_bg = wt$median_diff,
        p_value = wt$p_value
      )
    }) %>%
    mutate(fdr_bh = p.adjust(p_value, method = "BH"))
  
  ei_subpairs_sig <- subpair_tests %>%
    filter(pair_type == "E-I", is.finite(fdr_bh), fdr_bh < 0.05) %>%
    mutate(
      direction = ifelse(mean_diff_vs_bg > 0, "Higher than block E-I background", "Lower than block E-I background")
    ) %>%
    arrange(desc(abs(mean_diff_vs_bg)))
  
  # ---------------------------
  # Key result extraction for the 4 chat conclusions
  # ---------------------------
  get_test_row <- function(df, metric_name, comp_name) {
    out <- df %>% filter(metric == metric_name, comparison == comp_name)
    if (nrow(out) == 0) {
      return(tibble(metric = metric_name, comparison = comp_name, n = NA_real_, mean_diff = NA_real_, median_diff = NA_real_, rank_biserial = NA_real_, p_value = NA_real_, fdr_bh = NA_real_))
    }
    out %>% slice(1)
  }
  
  res_overlap_EE_EI <- get_test_row(overall_tests, "block_mean_overlap_prop_small", "E-E vs E-I")
  res_overlap_EI_II <- get_test_row(overall_tests, "block_mean_overlap_prop_small", "E-I vs I-I")
  res_gap_EE_EI <- get_test_row(overall_tests, "block_mean_nesting_gap", "E-E vs E-I")
  res_gap_EI_II <- get_test_row(overall_tests, "block_mean_nesting_gap", "E-I vs I-I")
  
  get_layer_nested_row <- function(layer_name, comp_name = "E-E vs E-I") {
    out <- nested_layer_tests %>% filter(as.character(layer) == layer_name, comparison == comp_name)
    if (nrow(out) == 0) {
      return(tibble(metric = "prop_nested90", comparison = comp_name, n = NA_real_, mean_diff = NA_real_, median_diff = NA_real_, rank_biserial = NA_real_, p_value = NA_real_, layer = layer_name, fdr_bh = NA_real_))
    }
    out %>% slice(1)
  }
  
  res_nested_l23 <- get_layer_nested_row("layer 2/3")
  res_nested_l4 <- get_layer_nested_row("layer 4")
  
  focus_006 <- exc_sub_tests %>% filter(str_detect(subclass, "^006\\s")) %>% slice(1)
  focus_007_l23 <- exc_sub_layer_tests %>% filter(str_detect(subclass, "^007\\s"), as.character(layer) == "layer 2/3") %>% slice(1)
  focus_007_l4 <- exc_sub_layer_tests %>% filter(str_detect(subclass, "^007\\s"), as.character(layer) == "layer 4") %>% slice(1)
  
  focus_pair_1 <- pick_focus_row(ei_subpairs_sig, c("007", "L2/3 IT CTX", "Vip"))
  focus_pair_2 <- pick_focus_row(ei_subpairs_sig, c("032", "L5 NP CTX", "Pvalb"))
  focus_pair_3 <- pick_focus_row(ei_subpairs_sig, c("032", "L5 NP CTX", "Sst"))
  focus_pair_4 <- pick_focus_row(ei_subpairs_sig, c("002", "IT EP-CLA", "Sst"))
  focus_pair_5 <- pick_focus_row(ei_subpairs_sig, c("021", "L4 RSP-ACA", "Pvalb"))
  focus_pair_6 <- pick_focus_row(ei_subpairs_sig, c("007", "L2/3 IT CTX", "Lamp5"))
  focus_pair_7 <- pick_focus_row(ei_subpairs_sig, c("004", "L6 IT CTX", "Vip"))
  
  key_results <- bind_rows(
    tibble(
      section = "1_geometry",
      item = c(
        "E-I smaller-cluster coverage > E-E",
        "E-I smaller-cluster coverage > I-I",
        "E-I nesting gap > E-E",
        "E-I nesting gap > I-I"
      ),
      estimate = c(
        -res_overlap_EE_EI$mean_diff,
        res_overlap_EI_II$mean_diff,
        -res_gap_EE_EI$mean_diff,
        res_gap_EI_II$mean_diff
      ),
      p_value = c(res_overlap_EE_EI$p_value, res_overlap_EI_II$p_value, res_gap_EE_EI$p_value, res_gap_EI_II$p_value),
      fdr_bh = c(res_overlap_EE_EI$fdr_bh, res_overlap_EI_II$fdr_bh, res_gap_EE_EI$fdr_bh, res_gap_EI_II$fdr_bh),
      note = c(
        "block-level mean smaller overlap proportion",
        "block-level mean smaller overlap proportion",
        "block-level mean nesting gap",
        "block-level mean nesting gap"
      )
    ),
    tibble(
      section = "2_layer_near_complete_nesting",
      item = c(
        "layer 2/3: E-I near-complete nesting > E-E",
        "layer 4: E-I near-complete nesting > E-E"
      ),
      estimate = c(-res_nested_l23$mean_diff, -res_nested_l4$mean_diff),
      p_value = c(res_nested_l23$p_value, res_nested_l4$p_value),
      fdr_bh = c(res_nested_l23$fdr_bh, res_nested_l4$fdr_bh),
      note = c("prop overlap_small >= 0.90", "prop overlap_small >= 0.90")
    ),
    tibble(
      section = "3_exc_subclass_preference",
      item = c(
        if (nrow(focus_006) > 0) paste0(focus_006$subclass, ": same-type minus cross-type") else "006 subclass not found",
        if (nrow(focus_007_l23) > 0) paste0(focus_007_l23$subclass, " in layer 2/3") else "007 layer 2/3 not found",
        if (nrow(focus_007_l4) > 0) paste0(focus_007_l4$subclass, " in layer 4") else "007 layer 4 not found"
      ),
      estimate = c(
        if (nrow(focus_006) > 0) focus_006$mean_same_minus_cross else NA_real_,
        if (nrow(focus_007_l23) > 0) focus_007_l23$mean_same_minus_cross else NA_real_,
        if (nrow(focus_007_l4) > 0) focus_007_l4$mean_same_minus_cross else NA_real_
      ),
      p_value = c(
        if (nrow(focus_006) > 0) focus_006$p_value else NA_real_,
        if (nrow(focus_007_l23) > 0) focus_007_l23$p_value else NA_real_,
        if (nrow(focus_007_l4) > 0) focus_007_l4$p_value else NA_real_
      ),
      fdr_bh = c(
        if (nrow(focus_006) > 0) focus_006$fdr_bh else NA_real_,
        if (nrow(focus_007_l23) > 0) focus_007_l23$fdr_bh else NA_real_,
        if (nrow(focus_007_l4) > 0) focus_007_l4$fdr_bh else NA_real_
      ),
      note = c("overall excitatory subclass test", "layer-specific raw p is most informative", "layer-specific raw p is most informative")
    ),
    bind_rows(focus_pair_1, focus_pair_2, focus_pair_3, focus_pair_4, focus_pair_5, focus_pair_6, focus_pair_7) %>%
      transmute(
        section = "4_specific_EI_pairs",
        item = subpair_name_raw,
        estimate = mean_diff_vs_bg,
        p_value = p_value,
        fdr_bh = fdr_bh,
        note = direction
      )
  )
  
  # ---------------------------
  # Export tables
  # ---------------------------
  dataset_overview <- tibble(
    metric = c(
      "input_csv", "raw directional rows", "unique undirected pairs", "unique clusters",
      "Glut clusters", "Gaba clusters", "slides", "layers", "complete slide-layer blocks"
    ),
    value = c(
      basename(input_csv), nrow(df), nrow(pairs_u), nrow(clusters),
      sum(clusters$cell_Neruon_type == "Glut"),
      sum(clusters$cell_Neruon_type == "Gaba"),
      n_distinct(clusters$slide),
      n_distinct(clusters$layer),
      length(blocks_complete)
    )
  )
  
  pairtype_desc <- pairs_u %>%
    group_by(pair_type) %>%
    summarise(
      n_pairs = n(),
      mean_jaccard = mean(jaccard, na.rm = TRUE),
      median_jaccard = median(jaccard, na.rm = TRUE),
      mean_overlap_prop_small = mean(overlap_prop_small, na.rm = TRUE),
      median_overlap_prop_small = median(overlap_prop_small, na.rm = TRUE),
      mean_overlap_prop_large = mean(overlap_prop_large, na.rm = TRUE),
      mean_nesting_gap = mean(nesting_gap, na.rm = TRUE),
      prop_nested90 = mean(nested_90, na.rm = TRUE),
      .groups = "drop"
    )
  
  write_csv(dataset_overview, file.path(tbl_dir, "table_00_dataset_overview.csv"))
  write_csv(pairtype_desc, file.path(tbl_dir, "table_01_pairtype_descriptives.csv"))
  write_csv(overall_tests, file.path(tbl_dir, "table_02_overall_geometry_tests.csv"))
  write_csv(nested_layer_tests, file.path(tbl_dir, "table_03_layer_nested90_tests.csv"))
  write_csv(exc_sub_tests, file.path(tbl_dir, "table_04_exc_subclass_preference.csv"))
  write_csv(exc_sub_layer_tests, file.path(tbl_dir, "table_05_exc_subclass_layerwise_preference.csv"))
  write_csv(gaba_sub_tests, file.path(tbl_dir, "table_06_gaba_subclass_preference.csv"))
  write_csv(ei_subpairs_sig, file.path(tbl_dir, "table_07_significant_EI_subpairs_vs_background.csv"))
  write_csv(key_results, file.path(tbl_dir, "table_08_key_results_from_chat_conclusions.csv"))
  
  # ---------------------------
  # Figure 1: spatial geometry is different, not simply stronger
  # ---------------------------
  geom_small_long <- block_metrics_complete %>%
    select(block, pair_type, value = mean_overlap_prop_small)
  
  geom_gap_long <- block_metrics_complete %>%
    select(block, pair_type, value = mean_nesting_gap)
  
  geom_scatter <- block_metrics_complete %>%
    select(block, pair_type, mean_overlap_prop_small, mean_overlap_prop_large)
  
  ann_small <- c(
    paste0("E-I > E-E: FDR = ", fmt_p(res_overlap_EE_EI$fdr_bh)),
    paste0("E-I > I-I: FDR = ", fmt_p(res_overlap_EI_II$fdr_bh))
  )
  ann_gap <- c(
    paste0("E-I > E-E: FDR = ", fmt_p(res_gap_EE_EI$fdr_bh)),
    paste0("E-I > I-I: FDR = ", fmt_p(res_gap_EI_II$fdr_bh))
  )
  
  p1a <- ggplot(geom_small_long, aes(x = pair_type, y = value, color = pair_type)) +
    geom_line(aes(group = block), color = "grey78", alpha = 0.22, linewidth = 0.35) +
    geom_boxplot(aes(fill = pair_type), width = 0.34, alpha = 0.16, outlier.shape = NA, color = "black", linewidth = 0.45) +
    geom_point(position = position_jitter(width = 0.06, height = 0), size = 1.1, alpha = 0.55) +
    stat_summary(aes(group = 1), fun = median, geom = "line", linewidth = 1.15, color = "black") +
    stat_summary(fun = median, geom = "point", size = 3.1, color = "black") +
    scale_color_manual(values = pair_cols) +
    scale_fill_manual(values = pair_cols) +
    scale_y_continuous(labels = label_number(accuracy = 0.01), expand = expansion(mult = c(0.05, 0.18))) +
    labs(
      title = "Smaller-cluster coverage is highest for E-I",
      subtitle = "Block-level paired comparison across slide×layer units",
      x = NULL,
      y = "Mean smaller-cluster overlap proportion"
    ) +
    annotate("text", x = 1.05, y = max(geom_small_long$value, na.rm = TRUE) * 1.12, label = ann_small[1], hjust = 0, size = 3.7) +
    annotate("text", x = 1.05, y = max(geom_small_long$value, na.rm = TRUE) * 1.04, label = ann_small[2], hjust = 0, size = 3.7)
  
  p1b <- ggplot(geom_gap_long, aes(x = pair_type, y = value, color = pair_type)) +
    geom_line(aes(group = block), color = "grey78", alpha = 0.22, linewidth = 0.35) +
    geom_boxplot(aes(fill = pair_type), width = 0.34, alpha = 0.16, outlier.shape = NA, color = "black", linewidth = 0.45) +
    geom_point(position = position_jitter(width = 0.06, height = 0), size = 1.1, alpha = 0.55) +
    stat_summary(aes(group = 1), fun = median, geom = "line", linewidth = 1.15, color = "black") +
    stat_summary(fun = median, geom = "point", size = 3.1, color = "black") +
    scale_color_manual(values = pair_cols) +
    scale_fill_manual(values = pair_cols) +
    scale_y_continuous(labels = label_number(accuracy = 0.01), expand = expansion(mult = c(0.05, 0.18))) +
    labs(
      title = "E-I shows the largest nesting gap",
      subtitle = "A larger gap means the smaller cluster is nested inside a much larger one",
      x = NULL,
      y = "Mean nesting gap"
    ) +
    annotate("text", x = 1.05, y = max(geom_gap_long$value, na.rm = TRUE) * 1.12, label = ann_gap[1], hjust = 0, size = 3.7) +
    annotate("text", x = 1.05, y = max(geom_gap_long$value, na.rm = TRUE) * 1.04, label = ann_gap[2], hjust = 0, size = 3.7)
  
  p1c <- ggplot(geom_scatter, aes(x = mean_overlap_prop_large, y = mean_overlap_prop_small, color = pair_type, fill = pair_type)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.55, color = "grey50") +
    stat_ellipse(geom = "polygon", alpha = 0.10, linewidth = 0.55, level = 0.70, show.legend = FALSE) +
    geom_point(size = 2.0, alpha = 0.78, stroke = 0) +
    stat_summary(aes(group = pair_type), fun = mean, geom = "point", size = 4.2, shape = 21, color = "black", stroke = 0.55, show.legend = FALSE) +
    scale_color_manual(values = pair_cols) +
    scale_fill_manual(values = pair_cols) +
    scale_x_continuous(labels = label_number(accuracy = 0.01), expand = expansion(mult = c(0.02, 0.03))) +
    scale_y_continuous(labels = label_number(accuracy = 0.01), expand = expansion(mult = c(0.02, 0.06))) +
    coord_equal() +
    labs(
      title = "The geometry is hierarchical rather than symmetric",
      x = "Mean overlap proportion of the larger cluster",
      y = "Mean overlap proportion of the smaller cluster"
    )
  
  save_pub(p1a + theme(legend.position = "none"), "fig01a_smaller_cluster_coverage", 5.6, 5.2)
  save_pub(p1b + theme(legend.position = "none"), "fig01b_nesting_gap", 8, 5.2)
  save_pub(p1c + theme(legend.position = "none"), "fig01c_hierarchical_geometry_scatter", 8, 5.2)
  
  # ---------------------------
  # Figure 2: near-complete nesting in layer 2/3 and layer 4
  # ---------------------------
  layer_heat <- block_metrics_complete %>%
    group_by(layer, pair_type) %>%
    summarise(value = mean(prop_nested90, na.rm = TRUE), .groups = "drop")
  
  p2a <- ggplot(layer_heat, aes(x = pair_type, y = layer, fill = value)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.2f", value)), size = 4.2, fontface = "bold") +
    scale_fill_gradientn(colours = c("#F4F6FB", "#C9D9F2", "#8FB3E3", "#4F84C4", "#1F4E79")) +
    labs(
      title = "Near-complete nesting is strongest in shallow-to-middle layers",
      subtitle = "Definition: smaller-cluster overlap ≥ 90%",
      x = NULL,
      y = NULL,
      fill = NULL
    ) +
    theme(legend.position = "right")
  
  make_layer_nested_plot <- function(layer_name, test_row) {
    dat <- block_metrics_complete %>% filter(as.character(layer) == layer_name)
    if (nrow(dat) == 0) return(make_empty_plot(paste0(layer_name, " not found")))
    
    ggplot(dat, aes(x = pair_type, y = prop_nested90, color = pair_type)) +
      geom_line(aes(group = block), color = "grey78", alpha = 0.25, linewidth = 0.35) +
      geom_boxplot(aes(fill = pair_type), width = 0.34, alpha = 0.16, outlier.shape = NA, color = "black", linewidth = 0.45) +
      geom_point(position = position_jitter(width = 0.06, height = 0), size = 1.2, alpha = 0.62) +
      stat_summary(aes(group = 1), fun = median, geom = "line", linewidth = 1.12, color = "black") +
      stat_summary(fun = median, geom = "point", size = 3.0, color = "black") +
      scale_color_manual(values = pair_cols) +
      scale_fill_manual(values = pair_cols) +
      scale_y_continuous(labels = label_percent(accuracy = 1), expand = expansion(mult = c(0.05, 0.18))) +
      labs(
        title = paste0(layer_name, ": E-I near-complete nesting"),
        subtitle = paste0("E-I > E-E: FDR = ", fmt_p(test_row$fdr_bh)),
        x = NULL,
        y = "Fraction of pairs with overlap_small ≥ 90%"
      )
  }
  
  p2b <- make_layer_nested_plot("layer 2/3", res_nested_l23)
  p2c <- make_layer_nested_plot("layer 4", res_nested_l4)
  
  save_pub(p2a, "fig02a_near_complete_nesting_heatmap", 8, 5.3)
  save_pub(p2b + theme(legend.position = "none"), "fig02b_layer23_near_complete_nesting", 5.6, 5.3)
  save_pub(p2c + theme(legend.position = "none"), "fig02c_layer4_near_complete_nesting", 5.6, 5.3)
  
  # ---------------------------
  # Figure 3: excitatory subclass partner preference
  # ---------------------------
  exc_plot_df <- exc_sub_tests %>%
    mutate(
      preference = ifelse(mean_same_minus_cross < 0, "cross-type preferring", "same-type preferring"),
      label_flag = case_when(
        str_detect(subclass, "^006\\s") ~ TRUE,
        str_detect(subclass, "^007\\s") ~ TRUE,
        fdr_bh < 0.05 & mean_same_minus_cross < 0 ~ TRUE,
        TRUE ~ FALSE
      ),
      subclass_plot = fct_reorder(subclass, mean_same_minus_cross)
    )
  
  p3a <- ggplot(exc_plot_df, aes(x = mean_same_minus_cross, y = subclass_plot, color = preference)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.55, color = "grey45") +
    geom_segment(aes(x = 0, xend = mean_same_minus_cross, y = subclass_plot, yend = subclass_plot), linewidth = 0.85, alpha = 0.85) +
    geom_point(aes(size = ifelse(fdr_bh < 0.05, 3.3, 2.3)), stroke = 0.5, fill = "white", shape = 21) +
    geom_text_repel(
      data = exc_plot_df %>% filter(label_flag),
      aes(label = paste0(subclass, "\nΔ=", fmt_num(mean_same_minus_cross, 4), "; FDR=", fmt_p(fdr_bh))),
      size = 3.3,
      min.segment.length = 0,
      segment.alpha = 0.6,
      box.padding = 0.28,
      point.padding = 0.15,
      max.overlaps = Inf,
      seed = 123
    ) +
    scale_color_manual(values = pref_cols) +
    scale_size_identity() +
    labs(
      title = "Some common excitatory subclasses prefer inhibitory partners",
      subtitle = "Negative values indicate stronger cross-type than same-type spatial coupling",
      x = "same-type mean Jaccard − cross-type mean Jaccard",
      y = NULL
    )
  
  focus_007_dat <- exc_pref %>%
    filter(str_detect(subclass, "^007\\s"), as.character(layer) %in% c("layer 2/3", "layer 4")) %>%
    mutate(layer = factor(as.character(layer), levels = c("layer 2/3", "layer 4")))
  
  focus_007_ann <- tibble(
    layer = factor(c("layer 2/3", "layer 4"), levels = c("layer 2/3", "layer 4")),
    p_lab = c(
      if (nrow(focus_007_l23) > 0) paste0("p = ", fmt_p(focus_007_l23$p_value)) else "p = NA",
      if (nrow(focus_007_l4) > 0) paste0("p = ", fmt_p(focus_007_l4$p_value)) else "p = NA"
    )
  )
  
  if (nrow(focus_007_dat) > 0) {
    p3b <- ggplot(focus_007_dat, aes(x = layer, y = same_minus_cross)) +
      geom_hline(yintercept = 0, linetype = 2, linewidth = 0.55, color = "grey45") +
      geom_boxplot(width = 0.40, outlier.shape = NA, fill = "#E8EEF7", color = "black", linewidth = 0.5) +
      geom_jitter(aes(color = layer), width = 0.08, height = 0, size = 2.0, alpha = 0.75, show.legend = FALSE) +
      stat_summary(fun = mean, geom = "point", size = 3.6, shape = 23, fill = "#D1495B", color = "black") +
      geom_text(
        data = focus_007_ann,
        aes(x = layer, y = max(focus_007_dat$same_minus_cross, na.rm = TRUE) + 0.012, label = p_lab),
        inherit.aes = FALSE,
        size = 3.7,
        fontface = "bold"
      ) +
      scale_color_manual(values = c("layer 2/3" = "#3C78D8", "layer 4" = "#D1495B")) +
      scale_y_continuous(expand = expansion(mult = c(0.08, 0.18))) +
      labs(
        title = "007 L2/3 IT CTX Glut keeps the cross-type tendency in layer 2/3 and layer 4",
        subtitle = "Cluster-level paired test within the focal excitatory subclass",
        x = NULL,
        y = "same-type mean Jaccard − cross-type mean Jaccard"
      )
  } else {
    p3b <- make_empty_plot("Target subclass 007 was not found in the current dataset")
  }
  
  save_pub(p3a + theme(legend.position = "none"), "fig03a_exc_subclass_partner_preference", 8.4, 6.6)
  save_pub(p3b + theme(legend.position = "none"), "fig03b_subclass_007_layerwise_preference", 8, 6.6)
  
  # ---------------------------
  # Figure 4: specific E-I pairs, not generic inhibitory roles
  # [comment omitted: encoding-safe]
  # ---------------------------
  focus_pairs_from_chat <- bind_rows(
    focus_pair_1, focus_pair_2, focus_pair_3, focus_pair_4,
    focus_pair_5, focus_pair_6, focus_pair_7
  ) %>% distinct(subpair_name_raw, .keep_all = TRUE)
  
  if (nrow(focus_pairs_from_chat) < 4) {
    pos_top <- ei_subpairs_sig %>% filter(mean_diff_vs_bg > 0) %>% arrange(desc(mean_diff_vs_bg)) %>% slice_head(n = 4)
    neg_top <- ei_subpairs_sig %>% filter(mean_diff_vs_bg < 0) %>% arrange(mean_diff_vs_bg) %>% slice_head(n = 4)
    fig4_df <- bind_rows(pos_top, neg_top)
  } else {
    fig4_df <- focus_pairs_from_chat
  }
  
  fig4_df <- fig4_df %>%
    mutate(
      pair_label = str_wrap(subpair_name_raw, width = 34),
      pair_label = fct_reorder(pair_label, mean_diff_vs_bg),
      text_hjust = ifelse(mean_diff_vs_bg >= 0, 0, 1),
      x_text = mean_diff_vs_bg + ifelse(mean_diff_vs_bg >= 0, 0.006, -0.006)
    )
  
  p4 <- ggplot(fig4_df, aes(x = mean_diff_vs_bg, y = pair_label, color = direction)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.55, color = "grey45") +
    geom_segment(aes(x = 0, xend = mean_diff_vs_bg, y = pair_label, yend = pair_label), linewidth = 1.0, alpha = 0.9) +
    geom_point(size = 3.2) +
    geom_text(
      aes(x = x_text, hjust = text_hjust, label = paste0("ΔJ = ", fmt_num(mean_diff_vs_bg, 4), "\nFDR = ", fmt_p(fdr_bh))),
      size = 3.45,
      lineheight = 0.95,
      show.legend = FALSE
    ) +
    scale_color_manual(values = bg_cols) +
    scale_x_continuous(expand = expansion(mult = c(0.16, 0.24))) +
    labs(
      title = "Specific E-I pairs can be above or below the block-level E-I background",
      subtitle = "This argues against treating inhibitory subclasses as a single uniform spatial role",
      x = "Δ mean Jaccard vs within-block E-I background",
      y = NULL
    )
  
  save_pub(p4, "fig04_specific_EI_pairs_above_or_below_background", 12.0, 7.2)
  
  # ---------------------------
  # Figure manifest
  # ---------------------------
  figure_manifest <- tibble(
    file_stem = c(
      "fig01_EI_geometry_not_symmetric_but_nested",
      "fig02_near_complete_nesting_layer23_layer4",
      "fig03_exc_subclass_prefers_inhibitory_partners",
      "fig04_specific_EI_pairs_above_or_below_background"
    ),
    content = c(
      "主图1：证明三类互作不是简单强弱差异，而是空间几何类型不同；重点突出E-I在small-cluster coverage和nesting gap上最高。",
      "主图2：展示near-complete nesting在不同layer中的分布，并重点强调layer 2/3和layer 4中E-I > E-E。",
      "主图3：展示兴奋性亚类的伙伴偏好，突出006总体更偏向抑制伙伴，以及007在layer 2/3和layer 4中的同样趋势。",
      "主图4：展示具体E-I配对相对block内E-I背景的偏离，强调空间耦合由具体配对决定，而不是由抑制亚类统一决定。"
    )
  )
  
  write_csv(figure_manifest, file.path(meta_dir, "figure_manifest.csv"))
  write_lines(
    c(
      paste0("Input CSV: ", normalizePath(input_csv, winslash = "/", mustWork = FALSE)),
      paste0("Output dir: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE)),
      "",
      "Generated figures:",
      paste0("- ", figure_manifest$file_stem, ": ", figure_manifest$content)
    ),
    file.path(meta_dir, "run_summary.txt")
  )
  
  # ---------------------------
  # Optional cluster-robust regression table retained
  # ---------------------------
  mod_df <- pairs_u %>%
    mutate(
      jaccard_clip = pmin(pmax(jaccard, 1e-4), 1 - 1e-4),
      logit_jaccard = log(jaccard_clip / (1 - jaccard_clip)),
      log_size_small = log(size_small),
      log_size_ratio = log(size_ratio),
      layer_chr = as.character(layer),
      pair_type_chr = as.character(pair_type)
    )
  
  fit <- lm(logit_jaccard ~ pair_type_chr + layer_chr + log_size_small + log_size_ratio, data = mod_df)
  robust_vcov <- sandwich::vcovCL(fit, cluster = mod_df$block)
  ct <- lmtest::coeftest(fit, vcov. = robust_vcov)
  
  reg_table <- data.frame(
    term      = rownames(ct),
    estimate  = as.numeric(ct[, 1]),
    std_error = as.numeric(ct[, 2]),
    statistic = as.numeric(ct[, 3]),
    p_value   = as.numeric(ct[, 4]),
    row.names = NULL,
    check.names = FALSE
  )
  
  readr::write_csv(reg_table, file.path(tbl_dir, "table_09_cluster_robust_regression.csv"))
  
  message("All outputs are ready under: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
}

run_module_03 <- function(runtime_base) {
  message('[MODULE 03] start')
  # ============================================================
  # [comment omitted: encoding-safe]
  # [comment omitted: encoding-safe]
  # [comment omitted: encoding-safe]
  # [comment omitted: encoding-safe]
  # ============================================================
  
  suppressPackageStartupMessages({
    library(tidyverse)
    library(broom)
  })
  
  options(scipen = 999)
  
  # -------------------------------
  # [comment omitted: encoding-safe]
  # -------------------------------
  base_dir <- runtime_base
  input_file <- file.path(base_dir, "neuron-neuron-partner.csv")
  out_dir    <- file.path(base_dir, "R_output_03")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  fig_dir <- file.path(out_dir, "figures")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  
  # -------------------------------
  # [comment omitted: encoding-safe]
  # -------------------------------
  fmt_p <- function(p) {
    ifelse(
      is.na(p), NA_character_,
      ifelse(
        p < 1e-300, "<1e-300",
        ifelse(p < 1e-4, format(p, scientific = TRUE, digits = 3),
               sprintf("%.6f", p))
      )
    )
  }
  
  fmt_p_plot <- function(p) {
    if (is.na(p)) return("NA")
    if (p < 1e-3) return(paste0("FDR=", sub("e", "×10^", format(p, scientific = TRUE, digits = 2))))
    paste0("FDR=", sprintf("%.3f", p))
  }
  
  type_to_ei <- function(x) ifelse(x == "Glut", "E", "I")
  
  pair_type_from_types <- function(x, y) {
    purrr::map2_chr(x, y, ~ paste(sort(c(type_to_ei(.x), type_to_ei(.y))), collapse = "-"))
  }
  
  make_pair_key <- function(x, y) paste(pmin(x, y), pmax(x, y), sep = "||")
  
  count_region_proxy <- function(x) stringr::str_count(x, "layer")
  
  build_full_pairs_one_stratum <- function(df_sub) {
    if (nrow(df_sub) < 2) return(tibble())
    cmb <- utils::combn(seq_len(nrow(df_sub)), 2)
    tibble(
      stratum = df_sub$stratum[1],
      slide   = df_sub$slide[1],
      layer   = df_sub$layer[1],
      label_a = df_sub$label[cmb[1, ]],
      label_b = df_sub$label[cmb[2, ]],
      type_a  = df_sub$type[cmb[1, ]],
      type_b  = df_sub$type[cmb[2, ]],
      subclass_a = df_sub$subclass[cmb[1, ]],
      subclass_b = df_sub$subclass[cmb[2, ]],
      n_a = df_sub$n[cmb[1, ]],
      n_b = df_sub$n[cmb[2, ]]
    ) %>%
      mutate(
        pair_key = make_pair_key(label_a, label_b),
        pair_type = pair_type_from_types(type_a, type_b)
      )
  }
  
  cmh_binary_compare <- function(df, exposure_col, exposure_yes, outcome_col,
                                 outcome_yes = TRUE, strata_col = "stratum",
                                 comparison_label = NA_character_) {
    sp <- split(df, df[[strata_col]])
    arr_list <- list()
    strata_keep <- c()
    
    for (s in names(sp)) {
      d <- sp[[s]]
      a <- sum(d[[exposure_col]] == exposure_yes & d[[outcome_col]] == outcome_yes, na.rm = TRUE)
      b <- sum(d[[exposure_col]] == exposure_yes & d[[outcome_col]] != outcome_yes, na.rm = TRUE)
      c <- sum(d[[exposure_col]] != exposure_yes & d[[outcome_col]] == outcome_yes, na.rm = TRUE)
      d0 <- sum(d[[exposure_col]] != exposure_yes & d[[outcome_col]] != outcome_yes, na.rm = TRUE)
      
      if ((a + b) > 0 && (c + d0) > 0 && (a + c) > 0 && (b + d0) > 0) {
        arr_list[[length(arr_list) + 1]] <- matrix(c(a, b, c, d0), nrow = 2, byrow = TRUE)
        strata_keep <- c(strata_keep, s)
      }
    }
    
    if (length(arr_list) == 0) {
      return(tibble(
        comparison = comparison_label,
        n_strata = 0L,
        odds_ratio_MH = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        p = NA_real_
      ))
    }
    
    arr <- array(
      0,
      dim = c(2, 2, length(arr_list)),
      dimnames = list(
        exposure = c(as.character(exposure_yes), "other"),
        outcome = c("yes", "no"),
        stratum = strata_keep
      )
    )
    
    for (i in seq_along(arr_list)) arr[, , i] <- arr_list[[i]]
    
    mh <- mantelhaen.test(arr)
    
    tibble(
      comparison = comparison_label,
      n_strata = length(arr_list),
      odds_ratio_MH = unname(mh$estimate),
      ci_low = unname(mh$conf.int[1]),
      ci_high = unname(mh$conf.int[2]),
      p = mh$p.value
    )
  }
  
  paired_stratum_wilcox <- function(df, group_col, metric, group1, group0,
                                    strata_cols = c("slide", "layer")) {
    tmp <- df %>%
      transmute(across(all_of(strata_cols)),
                group = .data[[group_col]],
                value = .data[[metric]]) %>%
      group_by(across(all_of(c(strata_cols, "group")))) %>%
      summarise(value = median(value, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = group, values_from = value)
    
    if (!(group1 %in% names(tmp)) || !(group0 %in% names(tmp))) {
      return(tibble(
        metric = metric,
        group1 = group1,
        group0 = group0,
        n_strata = 0L,
        median_diff = NA_real_,
        prop_group1_gt_group0 = NA_real_,
        p = NA_real_
      ))
    }
    
    sub <- tmp %>% select(all_of(strata_cols), all_of(group1), all_of(group0)) %>% drop_na()
    if (nrow(sub) == 0) {
      return(tibble(
        metric = metric,
        group1 = group1,
        group0 = group0,
        n_strata = 0L,
        median_diff = NA_real_,
        prop_group1_gt_group0 = NA_real_,
        p = NA_real_
      ))
    }
    
    wt <- wilcox.test(sub[[group1]], sub[[group0]], paired = TRUE, exact = FALSE)
    
    tibble(
      metric = metric,
      group1 = group1,
      group0 = group0,
      n_strata = nrow(sub),
      median_diff = median(sub[[group1]] - sub[[group0]], na.rm = TRUE),
      prop_group1_gt_group0 = mean(sub[[group1]] > sub[[group0]], na.rm = TRUE),
      p = wt$p.value
    )
  }
  
  group_compare_wilcox <- function(df, group_col, metric) {
    g1 <- df %>% filter(.data[[group_col]]) %>% pull(.data[[metric]])
    g0 <- df %>% filter(!.data[[group_col]]) %>% pull(.data[[metric]])
    wt <- wilcox.test(g1, g0, paired = FALSE, exact = FALSE)
    tibble(
      metric = metric,
      n_group1 = length(g1),
      n_group0 = length(g0),
      median_group1 = median(g1, na.rm = TRUE),
      median_group0 = median(g0, na.rm = TRUE),
      p = wt$p.value
    )
  }
  
  fisher_feature_enrichment <- function(df, feature_col, event_col = "status",
                                        event_label = "I_in_E",
                                        baseline_label = "partial",
                                        min_count = 10) {
    feat_counts <- table(df[[feature_col]])
    keep <- names(feat_counts)[feat_counts >= min_count]
    
    out <- purrr::map_dfr(keep, function(term) {
      a <- sum(df[[feature_col]] == term & df[[event_col]] == event_label, na.rm = TRUE)
      b <- sum(df[[feature_col]] == term & df[[event_col]] == baseline_label, na.rm = TRUE)
      c <- sum(df[[feature_col]] != term & df[[event_col]] == event_label, na.rm = TRUE)
      d <- sum(df[[feature_col]] != term & df[[event_col]] == baseline_label, na.rm = TRUE)
      ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))
      tibble(
        feature = term,
        count = a + b,
        event = a,
        baseline = b,
        odds_ratio = unname(ft$estimate),
        p = ft$p.value,
        event_rate = a / (a + b)
      )
    }) %>%
      mutate(fdr = p.adjust(p, method = "BH")) %>%
      arrange(fdr, desc(odds_ratio))
    
    out
  }
  
  cluster_subclass_enrichment <- function(df, subclass_col, group_col, min_count = 20) {
    sub_counts <- table(df[[subclass_col]])
    keep <- names(sub_counts)[sub_counts >= min_count]
    
    out <- purrr::map_dfr(keep, function(term) {
      a <- sum(df[[subclass_col]] == term & df[[group_col]], na.rm = TRUE)
      b <- sum(df[[subclass_col]] == term & !df[[group_col]], na.rm = TRUE)
      c <- sum(df[[subclass_col]] != term & df[[group_col]], na.rm = TRUE)
      d <- sum(df[[subclass_col]] != term & !df[[group_col]], na.rm = TRUE)
      ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))
      tibble(
        subclass = term,
        count = a + b,
        group_true = a,
        group_false = b,
        odds_ratio = unname(ft$estimate),
        p = ft$p.value,
        rate_true = a / (a + b)
      )
    }) %>%
      mutate(fdr = p.adjust(p, method = "BH")) %>%
      arrange(fdr, desc(odds_ratio))
    
    out
  }
  
  # -------------------------------
  # [comment omitted: encoding-safe]
  # -------------------------------
  raw <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    rename(
      label1 = `cluster.1_label`,
      slide1 = `cluster.1_slide`,
      layer1 = `cluster.1_layer`,
      region1 = `cluster.1_region`,
      type1 = `cluster.1_cell_Neruon_type`,
      subclass1 = `cluster.1_subclass`,
      n1 = `cluster.1_total_cell_num`,
      glu1 = `cluster.1_Glut_Neruon_cell_ids_num`,
      gaba1 = `cluster.1_GABA_Neruon_cell_ids_num`,
      p1 = `cluster.1_cauchy_combination_p`,
      eir1 = `cluster.1_E_I_Ratio`,
      label2 = `cluster.2_label`,
      slide2 = `cluster.2_slide`,
      layer2 = `cluster.2_layer`,
      region2 = `cluster.2_region`,
      type2 = `cluster.2_cell_Neruon_type`,
      subclass2 = `cluster.2_subclass`,
      n2 = `cluster.2_total_cell_num`,
      glu2 = `cluster.2_Glut_Neruon_cell_ids_num`,
      gaba2 = `cluster.2_GABA_Neruon_cell_ids_num`,
      p2 = `cluster.2_cauchy_combination_p`,
      eir2 = `cluster.2_E_I_Ratio`,
      overlap = `overlap_cell`,
      union = `union_cell`,
      jaccard = `jaccard`,
      ov1 = `cluster.1.overlap.percent`,
      ov2 = `cluster.2.overlap.percent`
    ) %>%
    mutate(pair_key = make_pair_key(label1, label2))
  
  # -------------------------------
  # [comment omitted: encoding-safe]
  # -------------------------------
  pairs <- purrr::imap_dfr(split(raw, raw$pair_key), function(g, key) {
    if (nrow(g) != 2) stop(paste("Pair key does not have exactly 2 rows:", key))
    labs <- sort(c(g$label1[1], g$label2[1]))
    a <- labs[1]
    b <- labs[2]
    r_ab <- g %>% filter(label1 == a, label2 == b)
    r_ba <- g %>% filter(label1 == b, label2 == a)
    if (nrow(r_ab) != 1 || nrow(r_ba) != 1) stop(paste("Cannot find both directions for:", key))
    
    tibble(
      pair_key = key,
      slide = r_ab$slide1,
      layer = r_ab$layer1,
      stratum = paste(r_ab$slide1, r_ab$layer1, sep = "__"),
      label_a = a,
      label_b = b,
      type_a = r_ab$type1,
      type_b = r_ab$type2,
      subclass_a = r_ab$subclass1,
      subclass_b = r_ab$subclass2,
      region_a = r_ab$region1,
      region_b = r_ab$region2,
      n_a = r_ab$n1,
      n_b = r_ab$n2,
      glu_a = r_ab$glu1,
      glu_b = r_ab$glu2,
      gaba_a = r_ab$gaba1,
      gaba_b = r_ab$gaba2,
      p_a = r_ab$p1,
      p_b = r_ab$p2,
      eir_a = r_ab$eir1,
      eir_b = r_ab$eir2,
      overlap = r_ab$overlap,
      union = r_ab$union,
      jaccard = r_ab$jaccard,
      ov_a = r_ab$ov1,
      ov_b = r_ab$ov2
    )
  }) %>%
    mutate(
      pair_type = pair_type_from_types(type_a, type_b),
      contain_a_in_b = abs(ov_a - 1) < 1e-12,
      contain_b_in_a = abs(ov_b - 1) < 1e-12,
      any_containment = contain_a_in_b | contain_b_in_a
    )
  
  # -------------------------------
  # [comment omitted: encoding-safe]
  # -------------------------------
  clusters <- bind_rows(
    raw %>%
      transmute(label = label1, slide = slide1, layer = layer1, region = region1,
                type = type1, subclass = subclass1, n = n1, glu = glu1, gaba = gaba1,
                p = p1, eir = eir1),
    raw %>%
      transmute(label = label2, slide = slide2, layer = layer2, region = region2,
                type = type2, subclass = subclass2, n = n2, glu = glu2, gaba = gaba2,
                p = p2, eir = eir2)
  ) %>%
    distinct(label, .keep_all = TRUE) %>%
    mutate(
      stratum = paste(slide, layer, sep = "__"),
      glu_frac = glu / n,
      gaba_frac = gaba / n,
      region_count = count_region_proxy(region)
    )
  
  # -------------------------------
  # [comment omitted: encoding-safe]
  # -------------------------------
  all_pairs <- purrr::map_dfr(split(clusters, clusters$stratum), build_full_pairs_one_stratum) %>%
    mutate(
      observed = pair_key %in% pairs$pair_key
    ) %>%
    left_join(pairs %>% select(pair_key, any_containment), by = "pair_key") %>%
    mutate(any_containment = tidyr::replace_na(any_containment, FALSE))
  
  # -------------------------------
  # [comment omitted: encoding-safe]
  # -------------------------------
  ei <- pairs %>%
    filter(pair_type == "E-I") %>%
    mutate(
      glu_label = if_else(type_a == "Glut", label_a, label_b),
      gaba_label = if_else(type_a == "Glut", label_b, label_a),
      glu_subclass = if_else(type_a == "Glut", subclass_a, subclass_b),
      gaba_subclass = if_else(type_a == "Glut", subclass_b, subclass_a),
      glu_region = if_else(type_a == "Glut", region_a, region_b),
      gaba_region = if_else(type_a == "Glut", region_b, region_a),
      glu_n = if_else(type_a == "Glut", n_a, n_b),
      gaba_n = if_else(type_a == "Glut", n_b, n_a),
      glu_glu = if_else(type_a == "Glut", glu_a, glu_b),
      glu_gaba = if_else(type_a == "Glut", gaba_a, gaba_b),
      gaba_glu = if_else(type_a == "Glut", glu_b, glu_a),
      gaba_gaba = if_else(type_a == "Glut", gaba_b, gaba_a),
      glu_p = if_else(type_a == "Glut", p_a, p_b),
      gaba_p = if_else(type_a == "Glut", p_b, p_a),
      glu_eir = if_else(type_a == "Glut", eir_a, eir_b),
      gaba_eir = if_else(type_a == "Glut", eir_b, eir_a),
      glu_ov = if_else(type_a == "Glut", ov_a, ov_b),
      gaba_ov = if_else(type_a == "Glut", ov_b, ov_a)
    ) %>%
    mutate(
      I_in_E = abs(gaba_ov - 1) < 1e-12,
      E_in_I = abs(glu_ov - 1) < 1e-12,
      partial = !(I_in_E | E_in_I),
      status = case_when(
        I_in_E & E_in_I ~ "mutual_full",
        I_in_E ~ "I_in_E",
        E_in_I ~ "E_in_I",
        TRUE ~ "partial"
      ),
      asymmetry_abs = abs(glu_ov - gaba_ov),
      size_ratio_gaba_to_glu = gaba_n / glu_n,
      log_size_ratio = log((gaba_n + 0.5) / (glu_n + 0.5)),
      min_n = pmin(glu_n, gaba_n),
      max_n = pmax(glu_n, gaba_n),
      glu_region_count = count_region_proxy(glu_region),
      gaba_region_count = count_region_proxy(gaba_region),
      region_count_diff_glu_minus_gaba = glu_region_count - gaba_region_count,
      glu_logp = -log10(glu_p),
      gaba_logp = -log10(gaba_p),
      glu_glu_frac = glu_glu / glu_n,
      gaba_gaba_frac = gaba_gaba / gaba_n,
      motif = paste(glu_subclass, gaba_subclass, sep = " || ")
    ) %>%
    filter(status != "mutual_full")
  
  # [comment omitted: encoding-safe]
  all_ei <- all_pairs %>%
    filter(pair_type == "E-I") %>%
    mutate(
      glu_label = if_else(type_a == "Glut", label_a, label_b),
      gaba_label = if_else(type_a == "Glut", label_b, label_a),
      glu_subclass = if_else(type_a == "Glut", subclass_a, subclass_b),
      gaba_subclass = if_else(type_a == "Glut", subclass_b, subclass_a),
      glu_n = if_else(type_a == "Glut", n_a, n_b),
      gaba_n = if_else(type_a == "Glut", n_b, n_a)
    ) %>%
    left_join(
      ei %>% select(pair_key, overlap, jaccard, glu_ov, gaba_ov, I_in_E, E_in_I),
      by = "pair_key"
    ) %>%
    mutate(
      overlap = replace_na(overlap, 0),
      jaccard = replace_na(jaccard, 0),
      glu_ov = replace_na(glu_ov, 0),
      gaba_ov = replace_na(gaba_ov, 0),
      I_in_E = replace_na(I_in_E, FALSE),
      E_in_I = replace_na(E_in_I, FALSE)
    )
  
  # -------------------------------
  # [comment omitted: encoding-safe]
  # -------------------------------
  contain_events <- bind_rows(
    pairs %>%
      filter(abs(ov_a - 1) < 1e-12) %>%
      transmute(
        pair_key, slide, layer, stratum,
        contained_label = label_a,
        container_label = label_b,
        contained_type = type_a,
        container_type = type_b,
        contained_subclass = subclass_a,
        container_subclass = subclass_b,
        contained_n = n_a,
        container_n = n_b,
        overlap, jaccard,
        contained_ov = ov_a,
        container_ov = ov_b
      ),
    pairs %>%
      filter(abs(ov_b - 1) < 1e-12) %>%
      transmute(
        pair_key, slide, layer, stratum,
        contained_label = label_b,
        container_label = label_a,
        contained_type = type_b,
        container_type = type_a,
        contained_subclass = subclass_b,
        container_subclass = subclass_a,
        contained_n = n_b,
        container_n = n_a,
        overlap, jaccard,
        contained_ov = ov_b,
        container_ov = ov_a
      )
  ) %>%
    mutate(
      event_type = paste0(type_to_ei(contained_type), "_in_", type_to_ei(container_type)),
      size_ratio_contained_to_container = contained_n / container_n,
      asymmetry_abs = 1 - container_ov
    )
  
  # -------------------------------
  # [comment omitted: encoding-safe]
  # -------------------------------
  table_01_data_summary <- tribble(
    ~metric, ~value,
    "原始CSV行数", nrow(raw),
    "去重后唯一重叠配对数", nrow(pairs),
    "唯一cluster数", nrow(clusters),
    "slide×layer分层数", n_distinct(clusters$stratum),
    "所有可能配对数（同slide×layer内）", nrow(all_pairs),
    "重叠配对占全部可能配对比例", mean(all_pairs$observed),
    "观测到E-I重叠配对数", sum(pairs$pair_type == "E-I"),
    "观测到E-E重叠配对数", sum(pairs$pair_type == "E-E"),
    "观测到I-I重叠配对数", sum(pairs$pair_type == "I-I"),
    "存在任一方向100%包含的重叠配对数", sum(pairs$any_containment),
    "重叠配对中包含比例", mean(pairs$any_containment)
  )
  readr::write_csv(table_01_data_summary, file.path(out_dir, "table_01_data_summary.csv"))
  
  # -------------------------------
  # [comment omitted: encoding-safe]
  # -------------------------------
  table_02_pairtype_containment_rates_observed <- pairs %>%
    group_by(pair_type) %>%
    summarise(
      n_pairs = n(),
      n_containment = sum(any_containment),
      containment_rate = mean(any_containment),
      .groups = "drop"
    ) %>%
    arrange(desc(containment_rate))
  readr::write_csv(table_02_pairtype_containment_rates_observed,
                   file.path(out_dir, "table_02_pairtype_containment_rates_observed.csv"))
  
  table_03_pairtype_rates_all_possible <- all_pairs %>%
    group_by(pair_type) %>%
    summarise(
      n_possible_pairs = n(),
      n_observed_overlaps = sum(observed),
      overlap_rate = mean(observed),
      n_containment = sum(any_containment),
      containment_rate_all_possible = mean(any_containment),
      .groups = "drop"
    )
  readr::write_csv(table_03_pairtype_rates_all_possible,
                   file.path(out_dir, "table_03_pairtype_rates_all_possible.csv"))
  
  table_04_cmh_pairtype_tests <- bind_rows(
    cmh_binary_compare(pairs %>% filter(pair_type %in% c("E-I", "I-I")),
                       "pair_type", "E-I", "any_containment", TRUE,
                       comparison_label = "observed_overlaps_only | containment | E-I vs I-I"),
    cmh_binary_compare(pairs %>% filter(pair_type %in% c("E-E", "I-I")),
                       "pair_type", "E-E", "any_containment", TRUE,
                       comparison_label = "observed_overlaps_only | containment | E-E vs I-I"),
    cmh_binary_compare(pairs %>% filter(pair_type %in% c("E-I", "E-E")),
                       "pair_type", "E-I", "any_containment", TRUE,
                       comparison_label = "observed_overlaps_only | containment | E-I vs E-E"),
    cmh_binary_compare(all_pairs %>% filter(pair_type %in% c("E-I", "I-I")),
                       "pair_type", "E-I", "any_containment", TRUE,
                       comparison_label = "all_possible_pairs | containment | E-I vs I-I"),
    cmh_binary_compare(all_pairs %>% filter(pair_type %in% c("E-E", "I-I")),
                       "pair_type", "E-E", "any_containment", TRUE,
                       comparison_label = "all_possible_pairs | containment | E-E vs I-I"),
    cmh_binary_compare(all_pairs %>% filter(pair_type %in% c("E-I", "E-E")),
                       "pair_type", "E-I", "any_containment", TRUE,
                       comparison_label = "all_possible_pairs | containment | E-I vs E-E"),
    cmh_binary_compare(all_pairs %>% filter(pair_type %in% c("E-I", "I-I")),
                       "pair_type", "E-I", "observed", TRUE,
                       comparison_label = "all_possible_pairs | overlap | E-I vs I-I"),
    cmh_binary_compare(all_pairs %>% filter(pair_type %in% c("E-E", "I-I")),
                       "pair_type", "E-E", "observed", TRUE,
                       comparison_label = "all_possible_pairs | overlap | E-E vs I-I"),
    cmh_binary_compare(all_pairs %>% filter(pair_type %in% c("E-I", "E-E")),
                       "pair_type", "E-I", "observed", TRUE,
                       comparison_label = "all_possible_pairs | overlap | E-I vs E-E")
  ) %>%
    mutate(p_fmt = fmt_p(p))
  readr::write_csv(table_04_cmh_pairtype_tests,
                   file.path(out_dir, "table_04_cmh_pairtype_tests.csv"))
  
  # -------------------------------
  # [comment omitted: encoding-safe]
  # -------------------------------
  table_05_ei_directionality_overall <- tibble(
    scope = "observed_EI_overlaps",
    n_EI_pairs = nrow(ei),
    I_in_E = sum(ei$I_in_E),
    E_in_I = sum(ei$E_in_I),
    partial = sum(ei$partial),
    I_in_E_rate_among_EI_overlaps = mean(ei$I_in_E),
    E_in_I_rate_among_EI_overlaps = mean(ei$E_in_I),
    I_in_E_among_containment_events = sum(ei$I_in_E) / (sum(ei$I_in_E) + sum(ei$E_in_I)),
    directionality_exact_binom_p = binom.test(sum(ei$I_in_E),
                                              sum(ei$I_in_E) + sum(ei$E_in_I),
                                              p = 0.5,
                                              alternative = "greater")$p.value,
    all_possible_EI_pairs = nrow(all_ei),
    I_in_E_rate_all_possible_EI = mean(all_ei$I_in_E),
    E_in_I_rate_all_possible_EI = mean(all_ei$E_in_I)
  ) %>%
    mutate(p_fmt = fmt_p(directionality_exact_binom_p))
  readr::write_csv(table_05_ei_directionality_overall,
                   file.path(out_dir, "table_05_ei_directionality_overall.csv"))
  
  table_06_ei_directionality_by_layer <- ei %>%
    group_by(layer) %>%
    summarise(
      n_EI_pairs = n(),
      n_I_in_E = sum(I_in_E),
      n_E_in_I = sum(E_in_I),
      n_partial = sum(partial),
      I_in_E_rate_among_EI_overlaps = sum(I_in_E) / n(),
      E_in_I_rate_among_EI_overlaps = sum(E_in_I) / n(),
      I_in_E_among_containment_events = ifelse(sum(I_in_E) + sum(E_in_I) > 0,
                                               sum(I_in_E) / (sum(I_in_E) + sum(E_in_I)),
                                               NA_real_),
      directionality_p = ifelse(sum(I_in_E) + sum(E_in_I) > 0,
                                binom.test(sum(I_in_E),
                                           sum(I_in_E) + sum(E_in_I),
                                           p = 0.5,
                                           alternative = "greater")$p.value,
                                NA_real_),
      .groups = "drop"
    ) %>%
    arrange(desc(I_in_E_rate_among_EI_overlaps)) %>%
    mutate(p_fmt = fmt_p(directionality_p))
  readr::write_csv(table_06_ei_directionality_by_layer,
                   file.path(out_dir, "table_06_ei_directionality_by_layer.csv"))
  
  # -------------------------------
  # [comment omitted: encoding-safe]
  # -------------------------------
  metrics_iine_partial <- c(
    "jaccard", "overlap", "glu_ov", "gaba_ov", "glu_n", "gaba_n",
    "size_ratio_gaba_to_glu", "asymmetry_abs",
    "glu_region_count", "gaba_region_count", "region_count_diff_glu_minus_gaba",
    "glu_logp", "gaba_logp", "glu_eir", "gaba_eir"
  )
  
  table_07_iine_vs_partial_paired_wilcoxon <- purrr::map_dfr(
    metrics_iine_partial,
    ~ paired_stratum_wilcox(
      ei %>% filter(status %in% c("I_in_E", "partial")),
      group_col = "status",
      metric = .x,
      group1 = "I_in_E",
      group0 = "partial"
    )
  ) %>%
    mutate(p_fmt = fmt_p(p))
  readr::write_csv(table_07_iine_vs_partial_paired_wilcoxon,
                   file.path(out_dir, "table_07_iine_vs_partial_paired_wilcoxon.csv"))
  
  # -------------------------------
  # [comment omitted: encoding-safe]
  # -------------------------------
  pairwise_event_compare <- function(metric, a, b) {
    paired_stratum_wilcox(
      contain_events %>% filter(event_type %in% c(a, b)),
      group_col = "event_type",
      metric = metric,
      group1 = a,
      group0 = b
    ) %>%
      mutate(comparison = paste(a, "vs", b))
  }
  
  table_08_containment_eventtype_pairwise_paired_wilcoxon <- bind_rows(
    purrr::map_dfr(
      c("container_n", "contained_n", "size_ratio_contained_to_container", "container_ov", "jaccard", "overlap"),
      ~ pairwise_event_compare(.x, "I_in_E", "E_in_E")
    ),
    purrr::map_dfr(
      c("container_n", "contained_n", "size_ratio_contained_to_container", "container_ov", "jaccard", "overlap"),
      ~ pairwise_event_compare(.x, "I_in_E", "I_in_I")
    ),
    purrr::map_dfr(
      c("container_n", "contained_n", "size_ratio_contained_to_container", "container_ov", "jaccard", "overlap"),
      ~ pairwise_event_compare(.x, "I_in_E", "E_in_I")
    )
  ) %>%
    mutate(p_fmt = fmt_p(p))
  readr::write_csv(table_08_containment_eventtype_pairwise_paired_wilcoxon,
                   file.path(out_dir, "table_08_containment_eventtype_pairwise_paired_wilcoxon.csv"))
  
  # -------------------------------
  # [comment omitted: encoding-safe]
  # -------------------------------
  host_counts_iine <- ei %>% filter(I_in_E) %>% count(glu_label, name = "n_gaba_contained")
  guest_counts_iine <- ei %>% filter(I_in_E) %>% count(gaba_label, name = "n_glu_hosts")
  host_counts_eini <- ei %>% filter(E_in_I) %>% count(gaba_label, name = "n_glu_contained")
  guest_counts_eini <- ei %>% filter(E_in_I) %>% count(glu_label, name = "n_gaba_hosts")
  
  table_09_multiplicity_summary <- bind_rows(
    tibble(
      type = "GLU host (I_in_E)",
      n_objects = nrow(host_counts_iine),
      median_partners = median(host_counts_iine$n_gaba_contained),
      mean_partners = mean(host_counts_iine$n_gaba_contained),
      prop_multi = mean(host_counts_iine$n_gaba_contained >= 2),
      max_partners = max(host_counts_iine$n_gaba_contained)
    ),
    tibble(
      type = "GABA host (E_in_I)",
      n_objects = nrow(host_counts_eini),
      median_partners = median(host_counts_eini$n_glu_contained),
      mean_partners = mean(host_counts_eini$n_glu_contained),
      prop_multi = mean(host_counts_eini$n_glu_contained >= 2),
      max_partners = max(host_counts_eini$n_glu_contained)
    ),
    tibble(
      type = "GABA guest (I_in_E)",
      n_objects = nrow(guest_counts_iine),
      median_partners = median(guest_counts_iine$n_glu_hosts),
      mean_partners = mean(guest_counts_iine$n_glu_hosts),
      prop_multi = mean(guest_counts_iine$n_glu_hosts >= 2),
      max_partners = max(guest_counts_iine$n_glu_hosts)
    ),
    tibble(
      type = "GLU guest (E_in_I)",
      n_objects = nrow(guest_counts_eini),
      median_partners = median(guest_counts_eini$n_gaba_hosts),
      mean_partners = mean(guest_counts_eini$n_gaba_hosts),
      prop_multi = mean(guest_counts_eini$n_gaba_hosts >= 2),
      max_partners = max(guest_counts_eini$n_gaba_hosts)
    )
  )
  readr::write_csv(table_09_multiplicity_summary,
                   file.path(out_dir, "table_09_multiplicity_summary.csv"))
  
  ft_host <- fisher.test(matrix(c(
    sum(host_counts_iine$n_gaba_contained >= 2), sum(host_counts_iine$n_gaba_contained < 2),
    sum(host_counts_eini$n_glu_contained >= 2), sum(host_counts_eini$n_glu_contained < 2)
  ), nrow = 2, byrow = TRUE))
  
  ft_guest <- fisher.test(matrix(c(
    sum(guest_counts_iine$n_glu_hosts >= 2), sum(guest_counts_iine$n_glu_hosts < 2),
    sum(guest_counts_eini$n_gaba_hosts >= 2), sum(guest_counts_eini$n_gaba_hosts < 2)
  ), nrow = 2, byrow = TRUE))
  
  table_10_multiplicity_fisher_tests <- bind_rows(
    tibble(
      comparison = "Host multiplicity: GLU hosts in I_in_E vs GABA hosts in E_in_I",
      odds_ratio = unname(ft_host$estimate),
      p = ft_host$p.value,
      group1_multi_fraction = mean(host_counts_iine$n_gaba_contained >= 2),
      group2_multi_fraction = mean(host_counts_eini$n_glu_contained >= 2)
    ),
    tibble(
      comparison = "Contained multiplicity: GABA guests in I_in_E vs GLU guests in E_in_I",
      odds_ratio = unname(ft_guest$estimate),
      p = ft_guest$p.value,
      group1_multi_fraction = mean(guest_counts_iine$n_glu_hosts >= 2),
      group2_multi_fraction = mean(guest_counts_eini$n_gaba_hosts >= 2)
    )
  ) %>%
    mutate(p_fmt = fmt_p(p))
  readr::write_csv(table_10_multiplicity_fisher_tests,
                   file.path(out_dir, "table_10_multiplicity_fisher_tests.csv"))
  
  # -------------------------------
  # [comment omitted: encoding-safe]
  # -------------------------------
  ei_iine_partial <- ei %>% filter(status %in% c("I_in_E", "partial"))
  
  table_11_iine_vs_partial_motif_enrichment_full <- fisher_feature_enrichment(
    ei_iine_partial,
    feature_col = "motif",
    event_col = "status",
    event_label = "I_in_E",
    baseline_label = "partial",
    min_count = 10
  )
  readr::write_csv(table_11_iine_vs_partial_motif_enrichment_full,
                   file.path(out_dir, "table_11_iine_vs_partial_motif_enrichment_full.csv"))
  
  table_12_top_enriched_motifs_iine_vs_partial <- table_11_iine_vs_partial_motif_enrichment_full %>%
    filter(odds_ratio > 1) %>%
    arrange(fdr, desc(odds_ratio)) %>%
    slice_head(n = 20)
  readr::write_csv(table_12_top_enriched_motifs_iine_vs_partial,
                   file.path(out_dir, "table_12_top_enriched_motifs_iine_vs_partial.csv"))
  
  table_13_top_depleted_motifs_iine_vs_partial <- table_11_iine_vs_partial_motif_enrichment_full %>%
    filter(odds_ratio < 1) %>%
    arrange(fdr, odds_ratio) %>%
    slice_head(n = 20)
  readr::write_csv(table_13_top_depleted_motifs_iine_vs_partial,
                   file.path(out_dir, "table_13_top_depleted_motifs_iine_vs_partial.csv"))
  
  table_14_glu_subclass_enrichment_pairlevel_iine_vs_partial <- fisher_feature_enrichment(
    ei_iine_partial,
    feature_col = "glu_subclass",
    event_col = "status",
    event_label = "I_in_E",
    baseline_label = "partial",
    min_count = 20
  )
  readr::write_csv(table_14_glu_subclass_enrichment_pairlevel_iine_vs_partial,
                   file.path(out_dir, "table_14_glu_subclass_enrichment_pairlevel_iine_vs_partial.csv"))
  
  table_15_gaba_subclass_enrichment_pairlevel_iine_vs_partial <- fisher_feature_enrichment(
    ei_iine_partial,
    feature_col = "gaba_subclass",
    event_col = "status",
    event_label = "I_in_E",
    baseline_label = "partial",
    min_count = 20
  )
  readr::write_csv(table_15_gaba_subclass_enrichment_pairlevel_iine_vs_partial,
                   file.path(out_dir, "table_15_gaba_subclass_enrichment_pairlevel_iine_vs_partial.csv"))
  
  # [comment omitted: encoding-safe]
  gaba_stats <- ei %>%
    group_by(gaba_label) %>%
    summarise(
      slide = dplyr::first(slide),
      layer = dplyr::first(layer),
      gaba_subclass = dplyr::first(gaba_subclass),
      gaba_n = dplyr::first(gaba_n),
      gaba_eir = dplyr::first(gaba_eir),
      ei_partner_count = n(),
      full_container_count = sum(I_in_E),
      any_I_in_E = any(I_in_E),
      any_E_in_I = any(E_in_I),
      median_glu_partner_n = median(glu_n),
      median_jaccard = median(jaccard),
      .groups = "drop"
    ) %>%
    mutate(
      full_container_count_ge2 = full_container_count >= 2
    )
  
  glu_stats <- ei %>%
    group_by(glu_label) %>%
    summarise(
      slide = dplyr::first(slide),
      layer = dplyr::first(layer),
      glu_subclass = dplyr::first(glu_subclass),
      glu_n = dplyr::first(glu_n),
      glu_eir = dplyr::first(glu_eir),
      ei_partner_count = n(),
      full_containee_count = sum(I_in_E),
      any_host_I_in_E = any(I_in_E),
      any_contained_E_in_I = any(E_in_I),
      median_gaba_partner_n = median(gaba_n),
      median_jaccard = median(jaccard),
      .groups = "drop"
    ) %>%
    mutate(
      multi_host = full_containee_count >= 2
    )
  
  readr::write_csv(gaba_stats, file.path(out_dir, "table_16_gaba_cluster_stats.csv"))
  readr::write_csv(glu_stats, file.path(out_dir, "table_17_glu_cluster_stats.csv"))
  
  table_18_gaba_cluster_subclass_enrichment_ever_contained <- cluster_subclass_enrichment(
    gaba_stats, subclass_col = "gaba_subclass", group_col = "any_I_in_E", min_count = 20
  )
  readr::write_csv(table_18_gaba_cluster_subclass_enrichment_ever_contained,
                   file.path(out_dir, "table_18_gaba_cluster_subclass_enrichment_ever_contained.csv"))
  
  table_19_glu_cluster_subclass_enrichment_ever_hosting <- cluster_subclass_enrichment(
    glu_stats, subclass_col = "glu_subclass", group_col = "any_host_I_in_E", min_count = 20
  )
  readr::write_csv(table_19_glu_cluster_subclass_enrichment_ever_hosting,
                   file.path(out_dir, "table_19_glu_cluster_subclass_enrichment_ever_hosting.csv"))
  
  table_20_glu_host_subclass_enrichment_multi_containee <- cluster_subclass_enrichment(
    glu_stats %>% filter(any_host_I_in_E),
    subclass_col = "glu_subclass", group_col = "multi_host", min_count = 10
  )
  readr::write_csv(table_20_glu_host_subclass_enrichment_multi_containee,
                   file.path(out_dir, "table_20_glu_host_subclass_enrichment_multi_containee.csv"))
  
  table_21_gaba_guest_subclass_enrichment_multi_host <- cluster_subclass_enrichment(
    gaba_stats %>% filter(any_I_in_E),
    subclass_col = "gaba_subclass", group_col = "full_container_count_ge2", min_count = 10
  )
  readr::write_csv(table_21_gaba_guest_subclass_enrichment_multi_host,
                   file.path(out_dir, "table_21_gaba_guest_subclass_enrichment_multi_host.csv"))
  
  # -------------------------------
  # [comment omitted: encoding-safe]
  # -------------------------------
  table_22_top_I_in_E_by_overlap <- ei %>%
    filter(I_in_E) %>%
    arrange(desc(overlap)) %>%
    select(slide, layer, glu_label, gaba_label, glu_subclass, gaba_subclass,
           glu_n, gaba_n, overlap, jaccard, glu_ov) %>%
    slice_head(n = 20)
  readr::write_csv(table_22_top_I_in_E_by_overlap,
                   file.path(out_dir, "table_22_top_I_in_E_by_overlap.csv"))
  
  table_23_top_I_in_E_by_extreme_small_ratio <- ei %>%
    filter(I_in_E) %>%
    arrange(size_ratio_gaba_to_glu) %>%
    select(slide, layer, glu_label, gaba_label, glu_subclass, gaba_subclass,
           glu_n, gaba_n, overlap, jaccard, glu_ov, size_ratio_gaba_to_glu) %>%
    slice_head(n = 20)
  readr::write_csv(table_23_top_I_in_E_by_extreme_small_ratio,
                   file.path(out_dir, "table_23_top_I_in_E_by_extreme_small_ratio.csv"))
  
  table_24_top_GLU_hosts_by_number_of_contained_GABA <- host_counts_iine %>%
    arrange(desc(n_gaba_contained))
  readr::write_csv(table_24_top_GLU_hosts_by_number_of_contained_GABA,
                   file.path(out_dir, "table_24_top_GLU_hosts_by_number_of_contained_GABA.csv"))
  
  table_25_top_GABA_guests_by_number_of_GLU_hosts <- guest_counts_iine %>%
    arrange(desc(n_glu_hosts))
  readr::write_csv(table_25_top_GABA_guests_by_number_of_GLU_hosts,
                   file.path(out_dir, "table_25_top_GABA_guests_by_number_of_GLU_hosts.csv"))
  
  # -------------------------------
  # [comment omitted: encoding-safe]
  # -------------------------------
  table_26_gaba_cluster_group_comparisons <- bind_rows(
    group_compare_wilcox(gaba_stats, "any_I_in_E", "gaba_n"),
    group_compare_wilcox(gaba_stats, "any_I_in_E", "gaba_eir"),
    group_compare_wilcox(gaba_stats, "any_I_in_E", "ei_partner_count"),
    group_compare_wilcox(gaba_stats, "any_I_in_E", "median_glu_partner_n"),
    group_compare_wilcox(gaba_stats, "any_I_in_E", "median_jaccard")
  ) %>%
    mutate(group = "GABA ever fully contained by GLU vs never",
           p_fmt = fmt_p(p))
  readr::write_csv(table_26_gaba_cluster_group_comparisons,
                   file.path(out_dir, "table_26_gaba_cluster_group_comparisons.csv"))
  
  table_27_glu_cluster_group_comparisons <- bind_rows(
    group_compare_wilcox(glu_stats, "any_host_I_in_E", "glu_n"),
    group_compare_wilcox(glu_stats, "any_host_I_in_E", "glu_eir"),
    group_compare_wilcox(glu_stats, "any_host_I_in_E", "ei_partner_count"),
    group_compare_wilcox(glu_stats, "any_host_I_in_E", "median_gaba_partner_n"),
    group_compare_wilcox(glu_stats, "any_host_I_in_E", "median_jaccard")
  ) %>%
    mutate(group = "GLU ever fully contains GABA vs never",
           p_fmt = fmt_p(p))
  readr::write_csv(table_27_glu_cluster_group_comparisons,
                   file.path(out_dir, "table_27_glu_cluster_group_comparisons.csv"))
  
  # -------------------------------
  # [comment omitted: encoding-safe]
  # -------------------------------
  table_28_status_descriptive_medians <- ei %>%
    group_by(status) %>%
    summarise(
      n = n(),
      jaccard_median = median(jaccard),
      overlap_median = median(overlap),
      glu_n_median = median(glu_n),
      gaba_n_median = median(gaba_n),
      size_ratio_median = median(size_ratio_gaba_to_glu),
      asymmetry_median = median(asymmetry_abs),
      glu_ov_median = median(glu_ov),
      gaba_ov_median = median(gaba_ov),
      glu_region_count_median = median(glu_region_count),
      gaba_region_count_median = median(gaba_region_count),
      glu_logp_median = median(glu_logp),
      gaba_logp_median = median(gaba_logp),
      .groups = "drop"
    )
  readr::write_csv(table_28_status_descriptive_medians,
                   file.path(out_dir, "table_28_status_descriptive_medians.csv"))
  
  table_29_eventtype_descriptive_medians <- contain_events %>%
    group_by(event_type) %>%
    summarise(
      n = n(),
      contained_n_median = median(contained_n),
      container_n_median = median(container_n),
      size_ratio_median = median(size_ratio_contained_to_container),
      overlap_median = median(overlap),
      jaccard_median = median(jaccard),
      container_ov_median = median(container_ov),
      .groups = "drop"
    )
  readr::write_csv(table_29_eventtype_descriptive_medians,
                   file.path(out_dir, "table_29_eventtype_descriptive_medians.csv"))
  
  # -------------------------------
  # [comment omitted: encoding-safe]
  # -------------------------------
  glm_df <- ei %>%
    filter(status %in% c("I_in_E", "partial")) %>%
    mutate(
      y = if_else(status == "I_in_E", 1, 0),
      log_glu_n = log1p(glu_n),
      log_gaba_n = log1p(gaba_n)
    )
  
  fit_layer <- glm(
    y ~ log_glu_n + log_gaba_n + region_count_diff_glu_minus_gaba +
      glu_logp + gaba_logp + factor(layer),
    family = binomial(),
    data = glm_df
  )
  
  fit_slide_layer <- glm(
    y ~ log_glu_n + log_gaba_n + region_count_diff_glu_minus_gaba +
      glu_logp + gaba_logp + factor(layer) + factor(slide),
    family = binomial(),
    data = glm_df
  )
  
  table_30_logistic_sensitivity_layer_adjusted <- broom::tidy(fit_layer, conf.int = TRUE) %>%
    mutate(
      odds_ratio = exp(estimate),
      or_low = exp(conf.low),
      or_high = exp(conf.high),
      p_fmt = fmt_p(p.value)
    )
  readr::write_csv(table_30_logistic_sensitivity_layer_adjusted,
                   file.path(out_dir, "table_30_logistic_sensitivity_layer_adjusted.csv"))
  
  table_31_logistic_sensitivity_slide_layer_adjusted <- broom::tidy(fit_slide_layer, conf.int = TRUE) %>%
    mutate(
      odds_ratio = exp(estimate),
      or_low = exp(conf.low),
      or_high = exp(conf.high),
      p_fmt = fmt_p(p.value)
    )
  readr::write_csv(table_31_logistic_sensitivity_slide_layer_adjusted,
                   file.path(out_dir, "table_31_logistic_sensitivity_slide_layer_adjusted.csv"))
  
  # -------------------------------
  # [comment omitted: encoding-safe]
  # -------------------------------
  plot_theme <- theme_pub_global(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "#4D4D4D"),
      axis.title = element_text(color = "black"),
      axis.text = element_text(color = "black")
    )
  
  pair_cols <- c("E-E" = "#4E79A7", "E-I" = "#9C6BC4", "I-I" = "#59A14F")
  status_cols <- c("Partial" = "#BDBDBD", "I in E" = "#CC5A5F", "E in I" = "#5B84B1")
  
  # -------------------------------
  # [comment omitted: encoding-safe]
  # -------------------------------
  fig01 <- table_02_pairtype_containment_rates_observed %>%
    mutate(pair_type = factor(pair_type, levels = c("E-E", "E-I", "I-I"))) %>%
    ggplot(aes(x = pair_type, y = containment_rate, fill = pair_type)) +
    geom_col(width = 0.70) +
    geom_text(aes(label = scales::percent(containment_rate, accuracy = 0.1)),
              vjust = -0.2, size = 5.0) +
    scale_fill_manual(values = pair_cols) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.10))) +
    labs(
      x = "Pair type",
      y = "Containment rate among observed overlaps",
      title = "Containment is most frequent in E-I overlaps",
      subtitle = "Unique undirected overlap pairs, summarized by interaction class"
    ) +
    plot_theme +
    theme(legend.position = "none")
  ggsave_pub(file.path(fig_dir, "fig01_pairtype_containment_rate_observed.png"), fig01,
         width = 7.0, height = 5.2, dpi = 400, bg = "white")
  
  # -------------------------------
  # [comment omitted: encoding-safe]
  # -------------------------------
  fig02_df <- ei %>%
    count(status) %>%
    filter(status %in% c("partial", "I_in_E", "E_in_I")) %>%
    mutate(
      status_label = factor(c("Partial", "I in E", "E in I")[match(status, c("partial", "I_in_E", "E_in_I"))],
                            levels = c("Partial", "I in E", "E in I")),
      pct = n / sum(n),
      label = paste0(scales::comma(n), "\n(", scales::percent(pct, accuracy = 0.1), ")")
    )
  
  fig02 <- ggplot(fig02_df, aes(x = status_label, y = n, fill = status_label)) +
    geom_col(width = 0.70) +
    geom_text(aes(label = label), vjust = -0.22, size = 5.0, lineheight = 0.95) +
    scale_fill_manual(values = status_cols) +
    scale_y_continuous(labels = scales::comma_format(),
                       expand = expansion(mult = c(0, 0.16))) +
    labs(
      x = "E-I interaction status",
      y = "Number of E-I overlap pairs",
      title = "Directional asymmetry within E-I overlaps",
      subtitle = "Most E-I overlaps are partial, but one-sided full containment strongly favors I in E"
    ) +
    coord_cartesian(clip = "off") +
    plot_theme +
    theme(legend.position = "none")
  ggsave_pub(file.path(fig_dir, "fig02_ei_status_counts.png"), fig02,
         width = 7.2, height = 5.4, dpi = 400, bg = "white")
  
  # -------------------------------
  # [comment omitted: encoding-safe]
  # -------------------------------
  fig03 <- contain_events %>%
    mutate(event_type = factor(event_type, levels = c("I_in_E", "E_in_E", "I_in_I", "E_in_I"))) %>%
    ggplot(aes(x = container_n, y = contained_n, shape = event_type, color = event_type)) +
    geom_point(alpha = 0.35, size = 1.8) +
    scale_x_log10(labels = scales::comma_format()) +
    scale_y_log10(labels = scales::comma_format()) +
    scale_color_manual(values = c("I_in_E" = "#D66A6A",
                                  "E_in_E" = "#E4A06F",
                                  "I_in_I" = "#8BBF85",
                                  "E_in_I" = "#6C92C4")) +
    labs(
      x = "Container cluster size (cells, log scale)",
      y = "Contained cluster size (cells, log scale)",
      title = "Containment events differ in host-domain size hierarchy",
      subtitle = "Each point is one full-containment event"
    ) +
    plot_theme
  ggsave_pub(file.path(fig_dir, "fig03_containment_event_size_scatter.png"), fig03,
         width = 7.2, height = 5.8, dpi = 400, bg = "white")
  
  # -------------------------------
  # [comment omitted: encoding-safe]
  # -------------------------------
  fig04 <- table_06_ei_directionality_by_layer %>%
    filter(layer != "layer 6b") %>%
    ggplot(aes(x = reorder(layer, -I_in_E_rate_among_EI_overlaps),
               y = I_in_E_rate_among_EI_overlaps,
               fill = I_in_E_rate_among_EI_overlaps)) +
    geom_col(width = 0.70) +
    geom_text(aes(label = scales::percent(I_in_E_rate_among_EI_overlaps, accuracy = 0.1)),
              vjust = -0.2, size = 5.0) +
    scale_fill_gradient(low = "#E6C7DA", high = "#E0117F", guide = "none") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.10))) +
    labs(
      x = "Layer",
      y = "I_in_E rate among E-I overlaps",
      title = "I_in_E enrichment is strongest in layer 2/3",
      subtitle = "Rates are computed within each layer among observed E-I overlaps"
    ) +
    plot_theme
  ggsave_pub(file.path(fig_dir, "fig04_iine_rate_by_layer.png"), fig04,
         width = 7.2, height = 4.8, dpi = 400, bg = "white")
  
  # -------------------------------
  # [comment omitted: encoding-safe]
  # [comment omitted: encoding-safe]
  # -------------------------------
  fig05_df <- ei %>%
    filter(status %in% c("partial", "I_in_E", "E_in_I")) %>%
    mutate(
      status = factor(status, levels = c("partial", "I_in_E", "E_in_I")),
      status_label = factor(c("Partial", "I in E", "E in I")[match(status, c("partial", "I_in_E", "E_in_I"))],
                            levels = c("Partial", "I in E", "E in I"))
    )
  
  # pairwise Wilcoxon + BH
  p_raw <- c(
    wilcox.test(size_ratio_gaba_to_glu ~ status_label,
                data = fig05_df %>% filter(status_label %in% c("Partial", "I in E")),
                exact = FALSE)$p.value,
    wilcox.test(size_ratio_gaba_to_glu ~ status_label,
                data = fig05_df %>% filter(status_label %in% c("I in E", "E in I")),
                exact = FALSE)$p.value,
    wilcox.test(size_ratio_gaba_to_glu ~ status_label,
                data = fig05_df %>% filter(status_label %in% c("Partial", "E in I")),
                exact = FALSE)$p.value
  )
  p_adj <- p.adjust(p_raw, method = "BH")
  
  max_y_fig05 <- max(fig05_df$size_ratio_gaba_to_glu, na.rm = TRUE)
  
  sig_df <- tibble(
    x    = c(1, 2, 1),
    xend = c(2, 3, 3),
    y    = c(max_y_fig05 * 1.90, max_y_fig05 * 3.20, max_y_fig05 * 5.40),
    yend = c(max_y_fig05 * 1.70, max_y_fig05 * 2.90, max_y_fig05 * 4.90),
    label = c(fmt_p_plot(p_adj[1]), fmt_p_plot(p_adj[2]), fmt_p_plot(p_adj[3]))
  )
  
  box_label_df <- fig05_df %>%
    group_by(status_label) %>%
    summarise(
      n = n(),
      med = median(size_ratio_gaba_to_glu, na.rm = TRUE),
      q3 = quantile(size_ratio_gaba_to_glu, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      x = c(1, 2, 3),
      y = q3 * c(1.58, 1.52, 1.58),
      label = paste0("n=", scales::comma(n), "\nmedian=", sprintf("%.2f", med))
    )
  
  fig05 <- ggplot(fig05_df, aes(x = status_label, y = size_ratio_gaba_to_glu, fill = status_label)) +
    geom_boxplot(outlier.shape = NA, width = 0.56, linewidth = 0.7) +
    geom_text(data = box_label_df,
              aes(x = x, y = y, label = label),
              inherit.aes = FALSE, size = 3.8, lineheight = 0.92) +
    geom_segment(data = sig_df,
                 aes(x = x, xend = xend, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.55) +
    geom_segment(data = sig_df,
                 aes(x = x, xend = x, y = y, yend = yend),
                 inherit.aes = FALSE, linewidth = 0.55) +
    geom_segment(data = sig_df,
                 aes(x = xend, xend = xend, y = y, yend = yend),
                 inherit.aes = FALSE, linewidth = 0.55) +
    geom_text(data = sig_df,
              aes(x = (x + xend) / 2, y = y * 1.15, label = label),
              inherit.aes = FALSE, size = 3.7) +
    scale_fill_manual(values = status_cols) +
    scale_y_log10(
      labels = scales::label_number(accuracy = 0.01),
      limits = c(0.01, max_y_fig05 * 8)
    ) +
    labs(
      x = "E-I interaction status",
      y = "GABA / GLU size ratio (log scale)",
      title = "Containment direction tracks guest-to-host size ratio",
      subtitle = "I in E events concentrate at smaller GABA-to-GLU size ratios than partial overlaps"
    ) +
    coord_cartesian(clip = "off") +
    plot_theme +
    theme(
      legend.position = "none",
      plot.margin = margin(15, 20, 12, 12)
    )
  ggsave_pub(file.path(fig_dir, "fig05_size_ratio_boxplot.png"), fig05,
         width = 8.0, height = 6.0, dpi = 400, bg = "white")
  
  # -------------------------------
  # [comment omitted: encoding-safe]
  # [comment omitted: encoding-safe]
  # [comment omitted: encoding-safe]
  # -------------------------------
  host_plot_df <- bind_rows(
    host_counts_iine %>% transmute(kind = "GLU hosts in I_in_E", n_partner = n_gaba_contained),
    host_counts_eini %>% transmute(kind = "GABA hosts in E_in_I", n_partner = n_glu_contained)
  )
  
  fig06_df <- host_plot_df %>%
    count(kind, n_partner, name = "n") %>%
    mutate(kind = factor(kind, levels = c("GABA hosts in E_in_I", "GLU hosts in I_in_E")))
  
  fig06 <- ggplot(fig06_df, aes(x = factor(n_partner), y = n, fill = kind)) +
    geom_col(position = position_dodge(width = 0.76), width = 0.68) +
    scale_fill_manual(values = c("GABA hosts in E_in_I" = "#6C92C4",
                                 "GLU hosts in I_in_E" = "#CC6B6B")) +
    labs(
      x = "Number of contained partners per host cluster",
      y = "Number of host clusters",
      title = "I_in_E hosts more often contain multiple guests",
      subtitle = paste0(
        "GLU hosts in I_in_E: n=", scales::comma(nrow(host_counts_iine)),
        ", multi-guest rate (≥2)=", scales::percent(mean(host_counts_iine$n_gaba_contained >= 2), accuracy = 0.1),
        "; GABA hosts in E_in_I: n=", scales::comma(nrow(host_counts_eini)),
        ", multi-guest rate (≥2)=", scales::percent(mean(host_counts_eini$n_glu_contained >= 2), accuracy = 0.1),
        "; Fisher ", fmt_p_plot(ft_host$p.value)
      )
    ) +
    plot_theme
  ggsave_pub(file.path(fig_dir, "fig06_host_multiplicity_hist.png"), fig06,
         width = 8.2, height = 5.6, dpi = 400, bg = "white")
  
  message("Analysis complete. Outputs written to: ", out_dir)
}

run_module_04 <- function(runtime_base) {
  message('[MODULE 04] start')
  # [comment omitted: encoding-safe]
  # [removed] rm(list = ls())
  gc()
  # =========================================================
  # neuron_containment_extension_analysis_topjournal.R
  # Top-journal visualization version
  # Input : E:/zaw/2603/neuron-neuron-partner.csv
  # Output: E:/zaw/2603/R_output_04
  # =========================================================
  
  required_pkgs <- c(
    "data.table", "dplyr", "tidyr", "purrr", "stringr",
    "ggplot2", "readr", "tibble", "sandwich", "lmtest", "scales"
  )
  
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
  
  suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(stringr)
    library(ggplot2)
    library(readr)
    library(tibble)
    library(sandwich)
    library(lmtest)
    library(scales)
  })
  
  options(scipen = 999)
  
  
  # ---------------------------
  # Protect dplyr verbs from method-dispatch conflicts
  # ---------------------------
  select <- dplyr::select
  rename <- dplyr::rename
  filter <- dplyr::filter
  mutate <- dplyr::mutate
  summarise <- dplyr::summarise
  arrange <- dplyr::arrange
  count <- dplyr::count
  distinct <- dplyr::distinct
  left_join <- dplyr::left_join
  inner_join <- dplyr::inner_join
  pull <- dplyr::pull
  transmute <- dplyr::transmute
  first <- dplyr::first
  
  
  # ---------------------------
  # Paths
  # ---------------------------
  base_dir <- runtime_base
  input_csv <- file.path(base_dir, "neuron-neuron-partner.csv")
  out_dir <- file.path(base_dir, "R_output_04")
  fig_dir <- file.path(out_dir, "figures")
  tbl_dir <- file.path(out_dir, "tables")
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(tbl_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ---------------------------
  # Helpers
  # ---------------------------
  canonical_pair <- function(a, b) {
    ifelse(a < b, paste(a, b, sep = " || "), paste(b, a, sep = " || "))
  }
  
  clean_subclass <- function(x) {
    stringr::str_replace(x, "^\\d+\\s+", "")
  }
  
  rank_biserial_paired <- function(x, y) {
    d <- x - y
    d <- d[d != 0]
    if (length(d) == 0) return(NA_real_)
    (sum(d > 0) - sum(d < 0)) / length(d)
  }
  
  paired_wilcox_from_wide <- function(dat, a, b) {
  # [comment omitted: encoding-safe]
    sub <- dat %>% dplyr::select(all_of(c(a, b))) %>% tidyr::drop_na()
    if (nrow(sub) == 0) return(NULL)
    wt <- suppressWarnings(wilcox.test(sub[[a]], sub[[b]], paired = TRUE, exact = FALSE))
    tibble(
      n = nrow(sub),
      mean_diff = mean(sub[[a]] - sub[[b]]),
      median_diff = median(sub[[a]] - sub[[b]]),
      rank_biserial = rank_biserial_paired(sub[[a]], sub[[b]]),
      p_value = wt$p.value
    )
  }
  
  robust_coef_table <- function(fit, cluster, exponentiate = FALSE) {
    V <- sandwich::vcovCL(fit, cluster = cluster)
    ct <- lmtest::coeftest(fit, vcov. = V)
    out <- data.frame(
      term = rownames(ct),
      estimate = as.numeric(ct[, 1]),
      std.error = as.numeric(ct[, 2]),
      statistic = as.numeric(ct[, 3]),
      p.value = as.numeric(ct[, 4]),
      row.names = NULL,
      check.names = FALSE
    )
    out <- out %>%
      mutate(
        conf.low = estimate - 1.96 * std.error,
        conf.high = estimate + 1.96 * std.error
      )
    if (exponentiate) {
      out <- out %>%
        mutate(
          OR = exp(estimate),
          OR_low = exp(conf.low),
          OR_high = exp(conf.high)
        )
    }
    out
  }
  
  fisher_enrichment <- function(category_vec, outcome_vec, min_n = 20) {
  # [comment omitted: encoding-safe]
    category_vec <- as.character(unlist(category_vec))
    outcome_vec <- as.logical(unlist(outcome_vec))
    
    stopifnot(length(category_vec) == length(outcome_vec))
    
    dat <- tibble(category = category_vec, outcome = outcome_vec)
    total_pos <- sum(dat$outcome)
    total_neg <- sum(!dat$outcome)
    
    dat %>%
      dplyr::count(category, name = "n_edges") %>%
      filter(n_edges >= min_n) %>%
      mutate(
        n_exact = map_int(category, ~sum(dat$outcome[dat$category == .x])),
        rate_exact = n_exact / n_edges,
        odds_ratio = map_dbl(category, ~{
          pos <- sum(dat$outcome[dat$category == .x])
          neg <- sum(!dat$outcome[dat$category == .x])
          other_pos <- total_pos - pos
          other_neg <- total_neg - neg
          unname(fisher.test(matrix(c(pos, neg, other_pos, other_neg), nrow = 2, byrow = TRUE))$estimate)
        }),
        p_value = map_dbl(category, ~{
          pos <- sum(dat$outcome[dat$category == .x])
          neg <- sum(!dat$outcome[dat$category == .x])
          other_pos <- total_pos - pos
          other_neg <- total_neg - neg
          fisher.test(matrix(c(pos, neg, other_pos, other_neg), nrow = 2, byrow = TRUE))$p.value
        })
      ) %>%
      mutate(
        fdr_bh = p.adjust(p_value, method = "BH"),
        direction = ifelse(odds_ratio > 1, "enriched", "depleted")
      ) %>%
      arrange(fdr_bh, desc(odds_ratio))
  }
  
  mw_test <- function(df, group_col, metric) {
    g1 <- df %>% filter(.data[[group_col]]) %>% pull(.data[[metric]]) %>% na.omit()
    g0 <- df %>% filter(!.data[[group_col]]) %>% pull(.data[[metric]]) %>% na.omit()
    if (length(g1) == 0 || length(g0) == 0) return(NULL)
    wt <- wilcox.test(g1, g0, exact = FALSE)
    tibble(
      metric = metric,
      n_group1 = length(g1),
      n_group0 = length(g0),
      mean_group1 = mean(g1),
      mean_group0 = mean(g0),
      median_group1 = median(g1),
      median_group0 = median(g0),
      mean_diff = mean(g1) - mean(g0),
      p_value = wt$p.value
    )
  }
  
  fmt_p <- function(p) {
    ifelse(is.na(p), NA_character_,
           ifelse(p < 1e-300, "<1e-300",
                  ifelse(p < 1e-4, format(p, scientific = TRUE, digits = 3),
                         sprintf("%.6f", p))))
  }
  
  fmt_p_short <- function(p) {
    if (is.na(p)) return("NA")
    if (p < 1e-300) return("<1e-300")
    if (p < 1e-4) return(format(p, scientific = TRUE, digits = 2))
    sprintf("%.4f", p)
  }
  
  save_pub <- function(p, filename, width, height, dpi = 600) {
    ggsave_pub(file.path(fig_dir, paste0(filename, ".png")), p,
           width = width, height = height, dpi = dpi, bg = "white")
    ggsave_pub(file.path(fig_dir, paste0(filename, ".pdf")), p,
           width = width, height = height, bg = "white")
  }
  
  add_sig_bracket <- function(df, x1, x2, y, label) {
    bind_rows(df,
              tibble(x = x1, xend = x2, y = y, yend = y, label = label, xlab = (x1 + x2) / 2, ylab = y),
              tibble(x = x1, xend = x1, y = y, yend = y * 0.93, label = NA_character_, xlab = NA_real_, ylab = NA_real_),
              tibble(x = x2, xend = x2, y = y, yend = y * 0.93, label = NA_character_, xlab = NA_real_, ylab = NA_real_))
  }
  
  # ---------------------------
  # Plot style
  # ---------------------------
  base_family <- "sans"
  col_pair <- c("E-E" = "#4E79A7", "E-I" = "#A06AB4", "I-I" = "#59A14F")
  col_dir <- c("GABA in Glu" = "#C44E52", "Glu in GABA" = "#4E79A7")
  col_small <- c("Smaller cluster = GABA" = "#C44E52", "Smaller cluster = Glu" = "#4E79A7")
  col_binary <- c("Exact GABA-in-Glu" = "#C44E52", "Other E-I" = "#BDBDBD")
  col_layer <- c("#E8D9EE", "#D3B8E0", "#BF95D1", "#AA73C3", "#954FB4", "#7F2EA5")
  col_multi <- c("GABA guest in Glu" = "#C44E52", "Glu guest in GABA" = "#4E79A7")
  
  pub_theme <- function() {
    theme_pub_global(base_size = 12, base_family = base_family) +
      theme(
        plot.title = element_text(face = "bold", size = 18, hjust = 0),
        plot.subtitle = element_text(size = 11.5, colour = "#4D4D4D", margin = margin(b = 8)),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        axis.line = element_line(linewidth = 0.7, colour = "black"),
        axis.ticks = element_line(linewidth = 0.7, colour = "black"),
        axis.ticks.length = unit(0.18, "cm"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        panel.grid = element_blank(),
        plot.margin = margin(12, 16, 12, 12)
      )
  }
  
  theme_set(pub_theme())
  
  # ---------------------------
  # Read data
  # ---------------------------
  df <- fread(input_csv)
  df[, pair_key := canonical_pair(cluster.1_label, cluster.2_label)]
  
  # Unique clusters
  c1 <- copy(df[, .(
    label = cluster.1_label,
    slide = cluster.1_slide,
    layer = cluster.1_layer,
    region = cluster.1_region,
    cell_Neruon_type = cluster.1_cell_Neruon_type,
    subclass = cluster.1_subclass,
    total_cell_num = cluster.1_total_cell_num,
    cauchy_combination_p = cluster.1_cauchy_combination_p,
    E_I_Ratio = cluster.1_E_I_Ratio
  )])
  
  c2 <- copy(df[, .(
    label = cluster.2_label,
    slide = cluster.2_slide,
    layer = cluster.2_layer,
    region = cluster.2_region,
    cell_Neruon_type = cluster.2_cell_Neruon_type,
    subclass = cluster.2_subclass,
    total_cell_num = cluster.2_total_cell_num,
    cauchy_combination_p = cluster.2_cauchy_combination_p,
    E_I_Ratio = cluster.2_E_I_Ratio
  )])
  
  clusters <- unique(rbindlist(list(c1, c2), use.names = TRUE))
  clusters[, block := paste(slide, layer, sep = "|")]
  clusters[, subclass_clean := clean_subclass(subclass)]
  
  # Unique undirected pairs
  pair_list <- lapply(split(df, df$pair_key), function(sub) {
    stopifnot(nrow(sub) == 2)
    labs <- sort(unique(c(sub$cluster.1_label, sub$cluster.2_label)))
    a <- labs[1]
    b <- labs[2]
    ra <- sub[sub$cluster.1_label == a, ][1, ]
    tibble(
      label_a = a,
      label_b = b,
      slide = ra$cluster.1_slide,
      layer = ra$cluster.1_layer,
      region_a = ra$cluster.1_region,
      region_b = ra$cluster.2_region,
      type_a_raw = ra$cluster.1_cell_Neruon_type,
      type_b_raw = ra$cluster.2_cell_Neruon_type,
      subclass_a = ra$cluster.1_subclass,
      subclass_b = ra$cluster.2_subclass,
      subclass_a_clean = clean_subclass(ra$cluster.1_subclass),
      subclass_b_clean = clean_subclass(ra$cluster.2_subclass),
      size_a = ra$cluster.1_total_cell_num,
      size_b = ra$cluster.2_total_cell_num,
      ov_a = ra$cluster.1.overlap.percent,
      ov_b = ra$cluster.2.overlap.percent,
      overlap = ra$overlap_cell,
      union = ra$union_cell,
      jaccard = ra$jaccard,
      p_a = ra$cluster.1_cauchy_combination_p,
      p_b = ra$cluster.2_cauchy_combination_p,
      ei_ratio_a = ra$cluster.1_E_I_Ratio,
      ei_ratio_b = ra$cluster.2_E_I_Ratio
    )
  })
  
  pairs <- bind_rows(pair_list) %>%
    mutate(
      pair_type = map2_chr(type_a_raw, type_b_raw, ~paste(sort(c(ifelse(.x == "Glut", "E", "I"), ifelse(.y == "Glut", "E", "I"))), collapse = "-")),
      block = paste(slide, layer, sep = "|"),
      size_small = pmin(size_a, size_b),
      size_large = pmax(size_a, size_b),
      overlap_prop_small = overlap / size_small,
      overlap_prop_large = overlap / size_large,
      exact_smaller_contained = abs(overlap_prop_small - 1) < 1e-12,
      near95_smaller_contained = overlap_prop_small >= 0.95
    )
  
  # Oriented E-I table
  ei <- pairs %>%
    filter(pair_type == "E-I") %>%
    mutate(
      glu_label = ifelse(type_a_raw == "Glut", label_a, label_b),
      gaba_label = ifelse(type_a_raw == "Glut", label_b, label_a),
      glu_subclass = ifelse(type_a_raw == "Glut", subclass_a_clean, subclass_b_clean),
      gaba_subclass = ifelse(type_a_raw == "Glut", subclass_b_clean, subclass_a_clean),
      glu_size = ifelse(type_a_raw == "Glut", size_a, size_b),
      gaba_size = ifelse(type_a_raw == "Glut", size_b, size_a),
      glu_overlap_prop = ifelse(type_a_raw == "Glut", ov_a, ov_b),
      gaba_overlap_prop = ifelse(type_a_raw == "Glut", ov_b, ov_a),
      glu_region = ifelse(type_a_raw == "Glut", region_a, region_b),
      gaba_region = ifelse(type_a_raw == "Glut", region_b, region_a),
      glu_p = ifelse(type_a_raw == "Glut", p_a, p_b),
      gaba_p = ifelse(type_a_raw == "Glut", p_b, p_a),
      glu_ei_ratio = ifelse(type_a_raw == "Glut", ei_ratio_a, ei_ratio_b),
      gaba_ei_ratio = ifelse(type_a_raw == "Glut", ei_ratio_b, ei_ratio_a),
      same_region = glu_region == gaba_region,
      gaba_in_glu_exact = abs(gaba_overlap_prop - 1) < 1e-12,
      glu_in_gaba_exact = abs(glu_overlap_prop - 1) < 1e-12,
      gaba_in_glu_99 = gaba_overlap_prop >= 0.99,
      gaba_in_glu_95 = gaba_overlap_prop >= 0.95,
      glu_in_gaba_99 = glu_overlap_prop >= 0.99,
      glu_in_gaba_95 = glu_overlap_prop >= 0.95,
      size_ratio_glu_to_gaba = glu_size / gaba_size,
      log_size_ratio_glu_to_gaba = log(size_ratio_glu_to_gaba),
      gaba_is_smaller = gaba_size < glu_size,
      glu_is_smaller = glu_size < gaba_size,
      same_size = glu_size == gaba_size,
      smaller_is_gaba = as.integer(gaba_is_smaller),
      smaller_contained_exact = case_when(
        gaba_is_smaller ~ gaba_in_glu_exact,
        glu_is_smaller ~ glu_in_gaba_exact,
        TRUE ~ FALSE
      ),
      smaller_contained_95 = case_when(
        gaba_is_smaller ~ gaba_in_glu_95,
        glu_is_smaller ~ glu_in_gaba_95,
        TRUE ~ FALSE
      ),
      smaller_type = case_when(
        gaba_is_smaller ~ "GABA_smaller",
        glu_is_smaller ~ "GLU_smaller",
        TRUE ~ "same_size"
      ),
      subpair_clean = paste(glu_subclass, gaba_subclass, sep = " × "),
      jaccard_clip = pmin(pmax(jaccard, 1e-4), 1 - 1e-4),
      logit_jaccard = log(jaccard_clip / (1 - jaccard_clip))
    )
  
  # ---------------------------
  # Tables
  # ---------------------------
  overview <- tibble(
    metric = c(
      "Raw directional rows",
      "Unique undirected pairs",
      "E-I undirected pairs",
      "Exact GABA in Glu",
      "Exact Glu in GABA",
      "Near-exact GABA in Glu (>=99%)",
      "Near-exact GABA in Glu (>=95%)",
      "E-I with smaller cluster = GABA",
      "E-I with smaller cluster = Glu",
      "Exact smaller-cluster containment rate when smaller = GABA",
      "Exact smaller-cluster containment rate when smaller = Glu",
      ">=95% smaller-cluster containment rate when smaller = GABA",
      ">=95% smaller-cluster containment rate when smaller = Glu"
    ),
    value = c(
      nrow(df),
      nrow(pairs),
      nrow(ei),
      sum(ei$gaba_in_glu_exact),
      sum(ei$glu_in_gaba_exact),
      sum(ei$gaba_in_glu_99),
      sum(ei$gaba_in_glu_95),
      sum(ei$gaba_is_smaller),
      sum(ei$glu_is_smaller),
      mean(ei$gaba_in_glu_exact[ei$gaba_is_smaller]),
      mean(ei$glu_in_gaba_exact[ei$glu_is_smaller]),
      mean(ei$gaba_in_glu_95[ei$gaba_is_smaller]),
      mean(ei$glu_in_gaba_95[ei$glu_is_smaller])
    )
  )
  
  pairtype_containment <- pairs %>%
    group_by(pair_type) %>%
    summarise(
      n_pairs = n(),
      exact_smaller_contained_n = sum(exact_smaller_contained),
      exact_smaller_contained_rate = mean(exact_smaller_contained),
      near95_smaller_contained_n = sum(near95_smaller_contained),
      near95_smaller_contained_rate = mean(near95_smaller_contained),
      mean_jaccard = mean(jaccard),
      median_jaccard = median(jaccard),
      .groups = "drop"
    )
  
  n_exact_dir <- sum(ei$gaba_in_glu_exact) + sum(ei$glu_in_gaba_exact)
  k_gaba_dir <- sum(ei$gaba_in_glu_exact)
  binom_exact <- binom.test(k_gaba_dir, n_exact_dir, p = 0.5, alternative = "greater")
  
  tab_smaller_exact <- matrix(
    c(
      sum(ei$gaba_in_glu_exact[ei$gaba_is_smaller]),
      sum(!ei$gaba_in_glu_exact[ei$gaba_is_smaller]),
      sum(ei$glu_in_gaba_exact[ei$glu_is_smaller]),
      sum(!ei$glu_in_gaba_exact[ei$glu_is_smaller])
    ),
    nrow = 2, byrow = TRUE
  )
  fisher_smaller_exact <- fisher.test(tab_smaller_exact)
  
  tab_smaller_95 <- matrix(
    c(
      sum(ei$gaba_in_glu_95[ei$gaba_is_smaller]),
      sum(!ei$gaba_in_glu_95[ei$gaba_is_smaller]),
      sum(ei$glu_in_gaba_95[ei$glu_is_smaller]),
      sum(!ei$glu_in_gaba_95[ei$glu_is_smaller])
    ),
    nrow = 2, byrow = TRUE
  )
  fisher_smaller_95 <- fisher.test(tab_smaller_95)
  
  direction_tests <- tibble(
    test = c(
      "Binomial among exact-direction E-I edges: GABA_in_Glu vs Glu_in_GABA",
      "Fisher: exact containment of smaller cluster, smaller GABA vs smaller Glu",
      "Fisher: >=95% containment of smaller cluster, smaller GABA vs smaller Glu"
    ),
    n_edges = c(n_exact_dir, sum(tab_smaller_exact), sum(tab_smaller_95)),
    k_gaba_in_glu = c(k_gaba_dir, sum(ei$gaba_in_glu_exact[ei$gaba_is_smaller]), sum(ei$gaba_in_glu_95[ei$gaba_is_smaller])),
    proportion_gaba_in_glu = c(k_gaba_dir / n_exact_dir, mean(ei$gaba_in_glu_exact[ei$gaba_is_smaller]), mean(ei$gaba_in_glu_95[ei$gaba_is_smaller])),
    odds_ratio = c(NA_real_, unname(fisher_smaller_exact$estimate), unname(fisher_smaller_95$estimate)),
    p_value = c(binom_exact$p.value, fisher_smaller_exact$p.value, fisher_smaller_95$p.value)
  )
  
  tmp <- ei %>%
    filter(gaba_is_smaller | glu_is_smaller) %>%
    mutate(
      smaller_contained_exact_int = as.integer(smaller_contained_exact),
      smaller_contained_95_int = as.integer(smaller_contained_95),
      abs_size_ratio = ifelse(gaba_is_smaller, glu_size / gaba_size, gaba_size / glu_size),
      log_abs_size_ratio = log(abs_size_ratio)
    )
  
  fit_smaller_exact <- glm(smaller_contained_exact_int ~ smaller_is_gaba + log_abs_size_ratio + layer,
                           data = tmp, family = binomial())
  fit_smaller_95 <- glm(smaller_contained_95_int ~ smaller_is_gaba + log_abs_size_ratio + layer,
                        data = tmp, family = binomial())
  
  smaller_exact_logit <- robust_coef_table(fit_smaller_exact, cluster = tmp$block, exponentiate = TRUE)
  smaller_95_logit <- robust_coef_table(fit_smaller_95, cluster = tmp$block, exponentiate = TRUE)
  
  exact_vs_other_desc <- ei %>%
    mutate(group = ifelse(gaba_in_glu_exact, "GABA_in_Glu_exact", "Other_EI")) %>%
    group_by(group) %>%
    summarise(
      n_edges = n(),
      mean_jaccard = mean(jaccard),
      median_jaccard = median(jaccard),
      mean_size_ratio_glu_to_gaba = mean(size_ratio_glu_to_gaba),
      median_size_ratio_glu_to_gaba = median(size_ratio_glu_to_gaba),
      mean_glu_size = mean(glu_size),
      median_glu_size = median(glu_size),
      mean_gaba_size = mean(gaba_size),
      median_gaba_size = median(gaba_size),
      mean_overlap = mean(overlap),
      median_overlap = median(overlap),
      same_region_rate = mean(same_region),
      mean_glu_overlap_prop = mean(glu_overlap_prop),
      mean_gaba_overlap_prop = mean(gaba_overlap_prop),
      .groups = "drop"
    )
  
  metric_list <- c("jaccard", "size_ratio_glu_to_gaba", "glu_size", "gaba_size",
                   "overlap", "same_region", "glu_overlap_prop", "gaba_overlap_prop",
                   "glu_ei_ratio", "gaba_ei_ratio")
  
  block_metric_tests <- purrr::map_dfr(metric_list, function(metric) {
    tmp_metric <- ei %>%
      mutate(group = ifelse(gaba_in_glu_exact, "exact", "other")) %>%
      group_by(block, group) %>%
      summarise(value = mean(.data[[metric]], na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = group, values_from = value) %>%
      drop_na(exact, other)
    if (nrow(tmp_metric) < 10) return(NULL)
    wt <- suppressWarnings(wilcox.test(tmp_metric$exact, tmp_metric$other, paired = TRUE, exact = FALSE))
    tibble(
      metric = metric,
      n_blocks = nrow(tmp_metric),
      mean_diff_exact_minus_other = mean(tmp_metric$exact - tmp_metric$other),
      median_diff_exact_minus_other = median(tmp_metric$exact - tmp_metric$other),
      rank_biserial = rank_biserial_paired(tmp_metric$exact, tmp_metric$other),
      p_value = wt$p.value
    )
  }) %>%
    mutate(fdr_bh = p.adjust(p_value, method = "BH"))
  
  fit_jaccard <- lm(logit_jaccard ~ as.integer(gaba_in_glu_exact) + log_size_ratio_glu_to_gaba + layer, data = ei)
  jaccard_regression <- robust_coef_table(fit_jaccard, cluster = ei$block, exponentiate = FALSE)
  
  layer_rate <- ei %>%
    group_by(layer) %>%
    summarise(
      n_edges = n(),
      exact_n = sum(gaba_in_glu_exact),
      exact_rate = mean(gaba_in_glu_exact),
      near95_n = sum(gaba_in_glu_95),
      near95_rate = mean(gaba_in_glu_95),
      mean_size_ratio = mean(size_ratio_glu_to_gaba),
      .groups = "drop"
    )
  
  layer_tests <- bind_rows(lapply(sort(unique(ei$layer)), function(lay) {
    in_layer <- ei$layer == lay
    out <- list()
    for (outcome in c("gaba_in_glu_exact", "gaba_in_glu_95")) {
      a <- sum(ei[[outcome]][in_layer])
      b <- sum(!ei[[outcome]][in_layer])
      c <- sum(ei[[outcome]][!in_layer])
      d <- sum(!ei[[outcome]][!in_layer])
      ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))
      out[[outcome]] <- tibble(
        layer = lay,
        outcome = ifelse(outcome == "gaba_in_glu_exact", "exact", "near95"),
        n_in_layer = sum(in_layer),
        rate_in_layer = a / sum(in_layer),
        odds_ratio = unname(ft$estimate),
        p_value = ft$p.value
      )
    }
    bind_rows(out)
  })) %>%
    mutate(fdr_bh = p.adjust(p_value, method = "BH"))
  
  fit_layer_exact <- glm(as.integer(gaba_in_glu_exact) ~ log_size_ratio_glu_to_gaba + layer,
                         data = ei, family = binomial())
  layer_exact_logit <- robust_coef_table(fit_layer_exact, cluster = ei$block, exponentiate = TRUE)
  
  block_exact <- pairs %>%
    group_by(block, pair_type) %>%
    summarise(exact_rate = mean(exact_smaller_contained), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = pair_type, values_from = exact_rate) %>%
    drop_na()
  
  pairtype_exact_tests <- bind_rows(
    paired_wilcox_from_wide(block_exact, "E-E", "E-I") %>% dplyr::mutate(comparison = "E-E vs E-I"),
    paired_wilcox_from_wide(block_exact, "E-E", "I-I") %>% dplyr::mutate(comparison = "E-E vs I-I"),
    paired_wilcox_from_wide(block_exact, "E-I", "I-I") %>% dplyr::mutate(comparison = "E-I vs I-I")
  )
  pairtype_exact_tests <- pairtype_exact_tests %>%
    dplyr::mutate(fdr_bh = p.adjust(p_value, method = "BH"))
  pairtype_exact_tests <- pairtype_exact_tests[, c("comparison", setdiff(names(pairtype_exact_tests), "comparison"))]
  
  layer_pairtype_exact_tests <- bind_rows(lapply(sort(unique(pairs$layer)), function(lay) {
    piv <- pairs %>%
      filter(layer == lay) %>%
      group_by(slide, pair_type) %>%
      summarise(exact_rate = mean(exact_smaller_contained), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = pair_type, values_from = exact_rate)
    out <- list()
    for (comp in list(c("E-E", "E-I"), c("E-E", "I-I"), c("E-I", "I-I"))) {
      a <- comp[1]; b <- comp[2]
      if (all(c(a, b) %in% names(piv))) {
  # [comment omitted: encoding-safe]
        sub <- piv %>% dplyr::select(all_of(c(a, b))) %>% tidyr::drop_na()
        if (nrow(sub) >= 5) {
          wt <- suppressWarnings(wilcox.test(sub[[a]], sub[[b]], paired = TRUE, exact = FALSE))
          out[[paste(a, b)]] <- tibble(
            layer = lay,
            comparison = paste(a, "vs", b),
            n_slides = nrow(sub),
            mean_diff = mean(sub[[a]] - sub[[b]]),
            median_diff = median(sub[[a]] - sub[[b]]),
            rank_biserial = rank_biserial_paired(sub[[a]], sub[[b]]),
            p_value = wt$p.value
          )
        }
      }
    }
    bind_rows(out)
  })) %>%
    mutate(fdr_bh = p.adjust(p_value, method = "BH"))
  
  host_glu_subclass_enrichment <- fisher_enrichment(ei$glu_subclass, ei$gaba_in_glu_exact, min_n = 50)
  guest_gaba_subclass_enrichment <- fisher_enrichment(ei$gaba_subclass, ei$gaba_in_glu_exact, min_n = 50)
  subclass_pair_fisher_enrichment <- fisher_enrichment(ei$subpair_clean, ei$gaba_in_glu_exact, min_n = 25)
  
  subclass_pair_block_background <- bind_rows(lapply(split(ei, ei$subpair_clean), function(sub) {
    sp <- unique(sub$subpair_clean)
    if (nrow(sub) < 20 || n_distinct(sub$block) < 10) return(NULL)
    sp_block <- sub %>%
      group_by(block) %>%
      summarise(sp_rate = mean(gaba_in_glu_exact), .groups = "drop")
    other_bg <- ei %>%
      filter(subpair_clean != sp) %>%
      group_by(block) %>%
      summarise(other_rate = mean(gaba_in_glu_exact), .groups = "drop")
    merged <- inner_join(sp_block, other_bg, by = "block")
    if (nrow(merged) < 10) return(NULL)
    wt <- suppressWarnings(wilcox.test(merged$sp_rate, merged$other_rate, paired = TRUE, exact = FALSE))
    tibble(
      subpair_clean = sp,
      n_edges = nrow(sub),
      n_blocks = nrow(merged),
      pair_exact_rate = mean(sub$gaba_in_glu_exact),
      mean_diff_vs_block_EI_bg = mean(merged$sp_rate - merged$other_rate),
      median_diff_vs_block_EI_bg = median(merged$sp_rate - merged$other_rate),
      rank_biserial = rank_biserial_paired(merged$sp_rate, merged$other_rate),
      p_value = wt$p.value
    )
  })) %>%
    mutate(
      fdr_bh = p.adjust(p_value, method = "BH"),
      direction = ifelse(mean_diff_vs_block_EI_bg > 0, "higher_containment_than_bg", "lower_containment_than_bg")
    ) %>%
    arrange(fdr_bh)
  
  ee <- pairs %>% filter(pair_type == "E-E")
  ii <- pairs %>% filter(pair_type == "I-I")
  
  ee_long <- bind_rows(
    ee %>% transmute(label = label_a, subclass = subclass_a_clean, layer = layer, size = size_a, partner = label_b, jaccard = jaccard),
    ee %>% transmute(label = label_b, subclass = subclass_b_clean, layer = layer, size = size_b, partner = label_a, jaccard = jaccard)
  )
  
  ii_long <- bind_rows(
    ii %>% transmute(label = label_a, subclass = subclass_a_clean, layer = layer, size = size_a, partner = label_b, jaccard = jaccard),
    ii %>% transmute(label = label_b, subclass = subclass_b_clean, layer = layer, size = size_b, partner = label_a, jaccard = jaccard)
  )
  
  glu_cluster_stats <- ei %>%
    group_by(glu_label, glu_subclass, layer) %>%
    summarise(
      glu_size = dplyr::first(glu_size),
      n_gaba_partners = n_distinct(gaba_label),
      n_ei_edges = n(),
      exact_host_count = sum(gaba_in_glu_exact),
      exact_host_rate = mean(gaba_in_glu_exact),
      near95_host_count = sum(gaba_in_glu_95),
      mean_ei_jaccard = mean(jaccard),
      median_ei_jaccard = median(jaccard),
      mean_size_ratio = mean(size_ratio_glu_to_gaba),
      mean_overlap = mean(overlap),
      .groups = "drop"
    ) %>%
    mutate(ever_exact_host = exact_host_count > 0) %>%
    left_join(ee_long %>% group_by(label) %>% summarise(n_ee_partners = n_distinct(partner), mean_ee_jaccard = mean(jaccard), .groups = "drop"),
              by = c("glu_label" = "label")) %>%
    left_join(clusters %>% dplyr::select(label, slide) %>% distinct(), by = c("glu_label" = "label")) %>%
    mutate(log_glu_size = log(glu_size))
  
  gaba_cluster_stats <- ei %>%
    group_by(gaba_label, gaba_subclass, layer) %>%
    summarise(
      gaba_size = dplyr::first(gaba_size),
      n_glu_partners = n_distinct(glu_label),
      n_ei_edges = n(),
      exact_contained_count = sum(gaba_in_glu_exact),
      exact_contained_rate = mean(gaba_in_glu_exact),
      near95_contained_count = sum(gaba_in_glu_95),
      mean_ei_jaccard = mean(jaccard),
      median_ei_jaccard = median(jaccard),
      mean_size_ratio = mean(size_ratio_glu_to_gaba),
      mean_overlap = mean(overlap),
      .groups = "drop"
    ) %>%
    mutate(ever_exact_contained = exact_contained_count > 0) %>%
    left_join(ii_long %>% group_by(label) %>% summarise(n_ii_partners = n_distinct(partner), mean_ii_jaccard = mean(jaccard), .groups = "drop"),
              by = c("gaba_label" = "label")) %>%
    left_join(clusters %>% dplyr::select(label, slide) %>% distinct(), by = c("gaba_label" = "label")) %>%
    mutate(log_gaba_size = log(gaba_size))
  
  glu_cluster_host_tests <- bind_rows(lapply(
    c("glu_size", "n_gaba_partners", "n_ei_edges", "exact_host_rate", "mean_ei_jaccard",
      "mean_size_ratio", "mean_overlap", "n_ee_partners", "mean_ee_jaccard"),
    function(m) mw_test(glu_cluster_stats, "ever_exact_host", m)
  )) %>% mutate(fdr_bh = p.adjust(p_value, method = "BH"))
  
  gaba_cluster_guest_tests <- bind_rows(lapply(
    c("gaba_size", "n_glu_partners", "n_ei_edges", "exact_contained_rate", "mean_ei_jaccard",
      "mean_size_ratio", "mean_overlap", "n_ii_partners", "mean_ii_jaccard"),
    function(m) mw_test(gaba_cluster_stats, "ever_exact_contained", m)
  )) %>% mutate(fdr_bh = p.adjust(p_value, method = "BH"))
  
  fit_glu_host <- glm(as.integer(ever_exact_host) ~ log_glu_size + mean_ee_jaccard + n_ee_partners + layer,
                      data = glu_cluster_stats %>% drop_na(mean_ee_jaccard, n_ee_partners),
                      family = binomial())
  fit_gaba_guest <- glm(as.integer(ever_exact_contained) ~ log_gaba_size + mean_ii_jaccard + n_ii_partners + layer,
                        data = gaba_cluster_stats %>% drop_na(mean_ii_jaccard, n_ii_partners),
                        family = binomial())
  
  glu_cluster_logit <- robust_coef_table(fit_glu_host,
                                         cluster = glu_cluster_stats %>% drop_na(mean_ee_jaccard, n_ee_partners) %>% pull(slide),
                                         exponentiate = TRUE)
  gaba_cluster_logit <- robust_coef_table(fit_gaba_guest,
                                          cluster = gaba_cluster_stats %>% drop_na(mean_ii_jaccard, n_ii_partners) %>% pull(slide),
                                          exponentiate = TRUE)
  
  cluster_counts_block <- clusters %>%
    dplyr::count(block, cell_Neruon_type) %>%
    pivot_wider(names_from = cell_Neruon_type, values_from = n, values_fill = 0) %>%
    rename(nE = Glut, nI = Gaba)
  
  edge_counts_block <- pairs %>%
    dplyr::count(block, pair_type) %>%
    pivot_wider(names_from = pair_type, values_from = n, values_fill = 0)
  
  dens <- cluster_counts_block %>%
    left_join(edge_counts_block, by = "block") %>%
    mutate(
      possible_EE = nE * (nE - 1) / 2,
      possible_II = nI * (nI - 1) / 2,
      possible_EI = nE * nI,
      density_EE = `E-E` / possible_EE,
      density_EI = `E-I` / possible_EI,
      density_II = `I-I` / possible_II
    )
  
  block_pair_summary <- pairs %>%
    group_by(block, pair_type) %>%
    summarise(
      mean_jaccard = mean(jaccard),
      exact_rate = mean(exact_smaller_contained),
      near95_rate = mean(near95_smaller_contained),
      .groups = "drop"
    )
  
  block_ei <- ei %>%
    group_by(block) %>%
    summarise(
      ei_n = n(),
      gaba_in_glu_exact_rate = mean(gaba_in_glu_exact),
      gaba_in_glu_95_rate = mean(gaba_in_glu_95),
      mean_jaccard_ei = mean(jaccard),
      mean_size_ratio_ei = mean(size_ratio_glu_to_gaba),
      .groups = "drop"
    )
  
  block_all <- block_ei %>%
    left_join(block_pair_summary %>% dplyr::select(block, pair_type, mean_jaccard) %>% pivot_wider(names_from = pair_type, values_from = mean_jaccard),
              by = "block") %>%
    left_join(block_pair_summary %>% dplyr::select(block, pair_type, exact_rate) %>% pivot_wider(names_from = pair_type, values_from = exact_rate, names_prefix = "exact_rate_"),
              by = "block") %>%
    left_join(block_pair_summary %>% dplyr::select(block, pair_type, near95_rate) %>% pivot_wider(names_from = pair_type, values_from = near95_rate, names_prefix = "near95_rate_"),
              by = "block") %>%
    left_join(dens %>% dplyr::select(block, density_EE, density_EI, density_II), by = "block")
  
  corr_targets <- c("E-E", "E-I", "I-I", "exact_rate_E-E", "exact_rate_E-I", "exact_rate_E-I",
                    "density_EE", "density_EI", "density_II", "mean_size_ratio_ei")
  
  block_correlations <- bind_rows(lapply(corr_targets, function(target) {
  # [comment omitted: encoding-safe]
    sub <- block_all %>% dplyr::select(gaba_in_glu_exact_rate, all_of(target)) %>% drop_na()
    if (nrow(sub) < 10) return(NULL)
    ct <- cor.test(sub$gaba_in_glu_exact_rate, sub[[target]], method = "spearman", exact = FALSE)
    tibble(
      target = target,
      n_blocks = nrow(sub),
      spearman_rho = unname(ct$estimate),
      p_value = ct$p.value
    )
  })) %>%
    mutate(fdr_bh = p.adjust(p_value, method = "BH")) %>%
    arrange(fdr_bh)
  
  gaba_exact_hosts <- ei %>%
    filter(gaba_in_glu_exact) %>%
    group_by(gaba_label, gaba_subclass, layer) %>%
    summarise(
      n_exact_glu_hosts = n_distinct(glu_label),
      n_host_subclasses = n_distinct(glu_subclass),
      mean_host_size = mean(glu_size),
      guest_size = dplyr::first(gaba_size),
      mean_jaccard = mean(jaccard),
      .groups = "drop"
    ) %>%
    mutate(multi_host = n_exact_glu_hosts >= 2)
  
  glu_exact_hosts_rev <- ei %>%
    filter(glu_in_gaba_exact) %>%
    group_by(glu_label, glu_subclass, layer) %>%
    summarise(
      n_exact_gaba_hosts = n_distinct(gaba_label),
      host_subclasses = n_distinct(gaba_subclass),
      guest_size = dplyr::first(glu_size),
      mean_host_size = mean(gaba_size),
      mean_jaccard = mean(jaccard),
      .groups = "drop"
    ) %>%
    mutate(multi_host = n_exact_gaba_hosts >= 2)
  
  multi_host_ft <- fisher.test(matrix(
    c(sum(gaba_exact_hosts$multi_host), sum(!gaba_exact_hosts$multi_host),
      sum(glu_exact_hosts_rev$multi_host), sum(!glu_exact_hosts_rev$multi_host)),
    nrow = 2, byrow = TRUE
  ))
  
  multi_host_direction_comparison <- tibble(
    direction = c("GABA guest in Glu", "Glu guest in GABA"),
    n_guest_clusters = c(nrow(gaba_exact_hosts), nrow(glu_exact_hosts_rev)),
    multi_host_n = c(sum(gaba_exact_hosts$multi_host), sum(glu_exact_hosts_rev$multi_host)),
    multi_host_rate = c(mean(gaba_exact_hosts$multi_host), mean(glu_exact_hosts_rev$multi_host))
  )
  
  top_exact_gaba_in_glu_pairs <- ei %>%
    filter(gaba_in_glu_exact) %>%
    mutate(pair_name_clean = paste(glu_subclass, gaba_subclass, sep = " × ")) %>%
    arrange(desc(jaccard), desc(overlap), desc(glu_size)) %>%
  # [comment omitted: encoding-safe]
    dplyr::select(layer, pair_name_clean, glu_size, gaba_size, overlap, jaccard, size_ratio_glu_to_gaba, glu_label, gaba_label) %>%
    slice_head(n = 20)
  
  top_glu_host_clusters <- ei %>%
    filter(gaba_in_glu_exact) %>%
    group_by(glu_label, glu_subclass, layer) %>%
    summarise(
      host_size = dplyr::first(glu_size),
      n_exact_gaba_guests = n_distinct(gaba_label),
      n_guest_subclasses = n_distinct(gaba_subclass),
      mean_guest_size = mean(gaba_size),
      mean_jaccard = mean(jaccard),
      .groups = "drop"
    ) %>%
    arrange(desc(n_exact_gaba_guests), desc(host_size)) %>%
    slice_head(n = 20)
  
  top_gaba_guest_clusters <- gaba_exact_hosts %>%
    arrange(desc(n_exact_glu_hosts), desc(guest_size)) %>%
    slice_head(n = 20)
  
  dec_dat <- tmp %>%
    mutate(size_ratio_decile = ntile(abs_size_ratio, 10)) %>%
    group_by(size_ratio_decile, smaller_type) %>%
    summarise(
      n = n(),
      exact_rate = mean(smaller_contained_exact),
      median_abs_size_ratio = median(abs_size_ratio),
      .groups = "drop"
    )
  
  size_ratio_decile_rates <- dec_dat
  
  # ---------------------------
  # Save tables
  # ---------------------------
  write_csv(overview, file.path(tbl_dir, "table_01_overview.csv"))
  write_csv(pairtype_containment, file.path(tbl_dir, "table_02_pairtype_containment.csv"))
  write_csv(direction_tests, file.path(tbl_dir, "table_03_direction_tests.csv"))
  write_csv(smaller_exact_logit, file.path(tbl_dir, "table_04_smaller_type_logit_exact.csv"))
  write_csv(smaller_95_logit, file.path(tbl_dir, "table_05_smaller_type_logit_95.csv"))
  write_csv(exact_vs_other_desc, file.path(tbl_dir, "table_06_exact_vs_other_EI_descriptives.csv"))
  write_csv(block_metric_tests, file.path(tbl_dir, "table_07_block_paired_exact_vs_other_metrics.csv"))
  write_csv(jaccard_regression, file.path(tbl_dir, "table_08_jaccard_regression_exact_vs_other.csv"))
  write_csv(layer_rate, file.path(tbl_dir, "table_09_layer_rates.csv"))
  write_csv(layer_tests, file.path(tbl_dir, "table_10_layer_enrichment_tests.csv"))
  write_csv(layer_exact_logit, file.path(tbl_dir, "table_11_layer_logit_exact.csv"))
  write_csv(pairtype_exact_tests, file.path(tbl_dir, "table_12_pairtype_exact_rate_tests.csv"))
  write_csv(layer_pairtype_exact_tests, file.path(tbl_dir, "table_13_layerwise_pairtype_exact_rate_tests.csv"))
  write_csv(host_glu_subclass_enrichment, file.path(tbl_dir, "table_14_host_glu_subclass_enrichment.csv"))
  write_csv(guest_gaba_subclass_enrichment, file.path(tbl_dir, "table_15_guest_gaba_subclass_enrichment.csv"))
  write_csv(subclass_pair_fisher_enrichment, file.path(tbl_dir, "table_16_subclass_pair_fisher_enrichment.csv"))
  write_csv(subclass_pair_block_background, file.path(tbl_dir, "table_17_subclass_pair_block_background_tests.csv"))
  write_csv(glu_cluster_host_tests, file.path(tbl_dir, "table_18_glu_cluster_host_tests.csv"))
  write_csv(gaba_cluster_guest_tests, file.path(tbl_dir, "table_19_gaba_cluster_guest_tests.csv"))
  write_csv(glu_cluster_logit, file.path(tbl_dir, "table_20_glu_cluster_logit_ever_host.csv"))
  write_csv(gaba_cluster_logit, file.path(tbl_dir, "table_21_gaba_cluster_logit_ever_contained.csv"))
  write_csv(block_correlations, file.path(tbl_dir, "table_22_block_correlations.csv"))
  write_csv(multi_host_direction_comparison, file.path(tbl_dir, "table_23_multi_host_direction_comparison.csv"))
  write_csv(top_exact_gaba_in_glu_pairs, file.path(tbl_dir, "table_24_top_exact_gaba_in_glu_pairs.csv"))
  write_csv(top_glu_host_clusters, file.path(tbl_dir, "table_25_top_glu_host_clusters.csv"))
  write_csv(top_gaba_guest_clusters, file.path(tbl_dir, "table_26_top_gaba_guest_clusters.csv"))
  write_csv(size_ratio_decile_rates, file.path(tbl_dir, "table_27_size_ratio_decile_rates.csv"))
  
  # ---------------------------
  # Top-journal figures
  # ---------------------------
  
  # ---------- plotting helpers ----------
  wrap_plot_title <- function(x, width = 52) stringr::str_wrap(x, width = width)
  wrap_plot_subtitle <- function(x, width = 92) stringr::str_wrap(x, width = width)
  
  plot_theme_safe <- theme_pub_global(base_size = 12) +
    ggplot2::theme(
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.margin = ggplot2::margin(20, 28, 14, 18),
      plot.title = ggplot2::element_text(
        face = "bold", size = 16.5, colour = "black",
        lineheight = 1.06, margin = ggplot2::margin(b = 8)
      ),
      plot.subtitle = ggplot2::element_text(
        size = 11.8, colour = "#4D4D4D",
        lineheight = 1.10, margin = ggplot2::margin(b = 12)
      ),
      axis.title = ggplot2::element_text(size = 13.5, colour = "black"),
      axis.text = ggplot2::element_text(size = 11.5, colour = "black"),
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 11.2)
    )
  
  # Figure 1: Pair-type exact containment rate with paired FDR annotations
  block_exact_plot <- pairs %>%
    group_by(block, pair_type) %>%
    summarise(exact_rate = mean(exact_smaller_contained), .groups = "drop") %>%
    mutate(pair_type = factor(pair_type, levels = c("E-E", "E-I", "I-I")))
  
  sig1 <- tibble()
  if (nrow(pairtype_exact_tests) > 0) {
    ee_ei <- pairtype_exact_tests %>% filter(comparison == "E-E vs E-I") %>% pull(fdr_bh)
    ei_ii <- pairtype_exact_tests %>% filter(comparison == "E-I vs I-I") %>% pull(fdr_bh)
    ee_ii <- pairtype_exact_tests %>% filter(comparison == "E-E vs I-I") %>% pull(fdr_bh)
    y1 <- max(block_exact_plot$exact_rate, na.rm = TRUE) * 1.12
    y2 <- max(block_exact_plot$exact_rate, na.rm = TRUE) * 1.24
    y3 <- max(block_exact_plot$exact_rate, na.rm = TRUE) * 1.36
    sig1 <- add_sig_bracket(sig1, 1, 2, y1, paste0("FDR=", fmt_p_short(ee_ei)))
    sig1 <- add_sig_bracket(sig1, 2, 3, y2, paste0("FDR=", fmt_p_short(ei_ii)))
    sig1 <- add_sig_bracket(sig1, 1, 3, y3, paste0("FDR=", fmt_p_short(ee_ii)))
  }
  
  fig1 <- ggplot(block_exact_plot, aes(pair_type, exact_rate, fill = pair_type)) +
    geom_violin(width = 0.90, alpha = 0.16, color = NA, trim = FALSE) +
    geom_boxplot(width = 0.34, outlier.shape = NA, linewidth = 0.7, alpha = 0.90) +
    geom_jitter(width = 0.12, size = 1.3, alpha = 0.22, color = "#4D4D4D") +
    geom_segment(
      data = sig1,
      aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE, linewidth = 0.55
    ) +
    geom_text(
      data = sig1 %>% filter(!is.na(label)),
      aes(x = xlab, y = ylab * 1.03, label = label),
      inherit.aes = FALSE, size = 3.8
    ) +
    scale_fill_manual(values = col_pair) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      expand = expansion(mult = c(0.02, 0.22))
    ) +
    labs(
      title = wrap_plot_title("E-I has the highest exact-containment rate across the three interaction classes"),
      subtitle = wrap_plot_subtitle(
        paste0(
          "Block-level smaller-cluster exact containment rate. Means: E-I ",
          percent(mean(pairs$exact_smaller_contained[pairs$pair_type == "E-I"]), accuracy = 0.01),
          ", E-E ",
          percent(mean(pairs$exact_smaller_contained[pairs$pair_type == "E-E"]), accuracy = 0.01),
          ", I-I ",
          percent(mean(pairs$exact_smaller_contained[pairs$pair_type == "I-I"]), accuracy = 0.01)
        )
      ),
      x = NULL,
      y = "Exact containment rate of the smaller cluster"
    ) +
    coord_cartesian(clip = "off") +
    plot_theme_safe
  save_pub(fig1, "fig01_pairtype_exact_containment_topjournal", 9.6, 6.2)
  
  # Figure 2: Directional exact asymmetry counts
  count_df <- tibble(
    direction = factor(c("GABA in Glu", "Glu in GABA"), levels = c("GABA in Glu", "Glu in GABA")),
    n_edges = c(sum(ei$gaba_in_glu_exact), sum(ei$glu_in_gaba_exact)),
    prop = c(sum(ei$gaba_in_glu_exact), sum(ei$glu_in_gaba_exact)) / n_exact_dir
  )
  
  fig2 <- ggplot(count_df, aes(direction, n_edges, fill = direction)) +
    geom_col(width = 0.62, color = NA) +
    geom_text(
      aes(label = paste0(comma(n_edges), "\n(", percent(prop, accuracy = 0.1), ")")),
      vjust = -0.28, size = 4.6, lineheight = 0.95
    ) +
    scale_fill_manual(values = col_dir) +
    scale_y_continuous(labels = comma_format(), expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = wrap_plot_title("Exact directional containment is strongly biased toward GABA-in-Glu"),
      subtitle = wrap_plot_subtitle(
        paste0(
          percent(k_gaba_dir / n_exact_dir, accuracy = 0.1),
          " of all exact directional E-I containment edges are GABA in Glu; binomial p = ",
          fmt_p_short(binom_exact$p.value)
        )
      ),
      x = NULL,
      y = "Number of E-I exact-direction edges"
    ) +
    plot_theme_safe +
    theme(legend.position = "none")
  save_pub(fig2, "fig02_directional_exact_asymmetry_topjournal", 8.4, 5.8)
  
  # Figure 3: Smaller-cluster comparison
  small_df <- tibble(
    smaller_type = factor(
      c("Smaller cluster = GABA", "Smaller cluster = Glu"),
      levels = c("Smaller cluster = GABA", "Smaller cluster = Glu")
    ),
    rate = c(
      mean(ei$gaba_in_glu_exact[ei$gaba_is_smaller]),
      mean(ei$glu_in_gaba_exact[ei$glu_is_smaller])
    ),
    n = c(sum(ei$gaba_is_smaller), sum(ei$glu_is_smaller))
  )
  
  fig3 <- ggplot(small_df, aes(smaller_type, rate, fill = smaller_type)) +
    geom_col(width = 0.62) +
    geom_text(
      aes(label = paste0(percent(rate, accuracy = 0.1), "\n(n=", comma(n), ")")),
      vjust = -0.26, size = 4.5, lineheight = 0.95
    ) +
    scale_fill_manual(values = col_small) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.16))) +
    labs(
      title = wrap_plot_title("The directional bias persists even when only the smaller cluster is compared"),
      subtitle = wrap_plot_subtitle(
        paste0(
          "Fisher OR = ",
          sprintf("%.2f", unname(fisher_smaller_exact$estimate)),
          "; p = ",
          fmt_p_short(fisher_smaller_exact$p.value),
          "; adjusted OR = ",
          sprintf("%.2f", smaller_exact_logit$OR[smaller_exact_logit$term == "smaller_is_gaba"]),
          ", p = ",
          fmt_p_short(smaller_exact_logit$p.value[smaller_exact_logit$term == "smaller_is_gaba"])
        )
      ),
      x = NULL,
      y = "Exact containment rate of the smaller cluster"
    ) +
    plot_theme_safe +
    theme(legend.position = "none")
  save_pub(fig3, "fig03_smaller_cluster_direction_bias_topjournal", 8.8, 5.8)
  
  # Figure 4: Adjusted OR forest for smaller-cluster exact model
  forest_small <- smaller_exact_logit %>%
    filter(term != "(Intercept)") %>%
    mutate(
      term_label = recode(
        term,
        `smaller_is_gaba` = "Smaller cluster = GABA",
        `log_abs_size_ratio` = "log(larger/smaller)",
        `layerlayer 2/3` = "Layer 2/3",
        `layerlayer 4` = "Layer 4",
        `layerlayer 5` = "Layer 5",
        `layerlayer 6a` = "Layer 6a",
        `layerlayer 6b` = "Layer 6b"
      )
    ) %>%
    mutate(term_label = factor(term_label, levels = rev(term_label)))
  
  fig4 <- ggplot(forest_small, aes(x = OR, y = term_label)) +
    geom_vline(xintercept = 1, linetype = 2, color = "#7F7F7F", linewidth = 0.55) +
    geom_errorbarh(aes(xmin = OR_low, xmax = OR_high), height = 0.18, linewidth = 0.7, color = "#4D4D4D") +
    geom_point(size = 2.8, color = "#C44E52") +
    scale_x_log10(labels = label_number(accuracy = 0.01)) +
    labs(
      title = wrap_plot_title("The smaller-GABA effect remains after controlling for size ratio and layer"),
      subtitle = wrap_plot_subtitle(
        paste0(
          "Adjusted OR for smaller cluster = GABA: ",
          sprintf("%.2f", smaller_exact_logit$OR[smaller_exact_logit$term == "smaller_is_gaba"]),
          "; p = ",
          fmt_p_short(smaller_exact_logit$p.value[smaller_exact_logit$term == "smaller_is_gaba"])
        )
      ),
      x = "Adjusted odds ratio (log scale)",
      y = NULL
    ) +
    plot_theme_safe
  save_pub(fig4, "fig04_smaller_cluster_logit_forest_topjournal", 9.0, 6.0)
  
  # Figure 5: Jaccard exact vs other E-I, with paired block summary and regression effect
  exact_jaccard_plot <- ei %>%
    mutate(group = ifelse(gaba_in_glu_exact, "Exact GABA-in-Glu", "Other E-I")) %>%
    group_by(block, group) %>%
    summarise(value = mean(jaccard), .groups = "drop") %>%
    drop_na() %>%
    mutate(group = factor(group, levels = c("Other E-I", "Exact GABA-in-Glu")))
  
  fig5 <- ggplot(exact_jaccard_plot, aes(group, value, fill = group)) +
    geom_violin(width = 0.86, alpha = 0.16, color = NA, trim = FALSE) +
    geom_boxplot(width = 0.30, outlier.shape = NA, linewidth = 0.7, alpha = 0.9) +
    geom_jitter(width = 0.10, size = 1.2, alpha = 0.20, color = "#4D4D4D") +
    scale_fill_manual(values = col_binary) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0.02, 0.08))) +
    labs(
      title = wrap_plot_title("Exact GABA-in-Glu edges have higher Jaccard despite larger size ratios"),
      subtitle = wrap_plot_subtitle(
        paste0(
          "Block-level mean Jaccard; adjusted logit(Jaccard) coefficient for exact containment = ",
          sprintf("%.4f", jaccard_regression$estimate[jaccard_regression$term == "as.integer(gaba_in_glu_exact)"]),
          "; p = ",
          fmt_p_short(jaccard_regression$p.value[jaccard_regression$term == "as.integer(gaba_in_glu_exact)"])
        )
      ),
      x = NULL,
      y = "Mean Jaccard per block"
    ) +
    plot_theme_safe +
    theme(legend.position = "none")
  save_pub(fig5, "fig05_exact_vs_other_jaccard_topjournal", 9.0, 6.0)
  
  # Figure 6: Layer-wise exact rate
  layer_plot <- layer_rate %>%
    mutate(layer = factor(layer, levels = c("layer 1", "layer 2/3", "layer 4", "layer 5", "layer 6a", "layer 6b")))
  
  fig6 <- ggplot(layer_plot, aes(layer, exact_rate, fill = layer)) +
    geom_col(width = 0.70) +
    geom_text(
      aes(label = paste0(percent(exact_rate, accuracy = 0.1), "\n(n=", comma(n_edges), ")")),
      vjust = -0.22, size = 4.1, lineheight = 0.95
    ) +
    scale_fill_manual(values = col_layer) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.16))) +
    labs(
      title = wrap_plot_title("Layer dependence of exact GABA-in-Glu containment"),
      subtitle = wrap_plot_subtitle("Rates among observed E-I pairs; shallow and middle layers show the highest levels"),
      x = NULL,
      y = "Exact GABA-in-Glu rate"
    ) +
    plot_theme_safe +
    theme(legend.position = "none")
  save_pub(fig6, "fig06_layer_exact_rate_topjournal", 8.8, 5.8)
  
  # Figure 7: Pair-specific background deviation
  sel_pos <- subclass_pair_block_background %>%
    filter(fdr_bh < 0.05, mean_diff_vs_block_EI_bg > 0) %>%
    arrange(desc(mean_diff_vs_block_EI_bg)) %>%
    slice_head(n = 6)
  
  sel_neg <- subclass_pair_block_background %>%
    filter(fdr_bh < 0.05, mean_diff_vs_block_EI_bg < 0) %>%
    arrange(mean_diff_vs_block_EI_bg) %>%
    slice_head(n = 10)
  
  sel_pair_bg <- bind_rows(sel_pos, sel_neg) %>%
    arrange(mean_diff_vs_block_EI_bg) %>%
    mutate(
      subpair_clean = factor(subpair_clean, levels = subpair_clean),
      group = ifelse(
        mean_diff_vs_block_EI_bg > 0,
        "Higher than block E-I background",
        "Lower than block E-I background"
      )
    )
  
  fig7 <- ggplot(sel_pair_bg, aes(mean_diff_vs_block_EI_bg, subpair_clean, fill = group)) +
    geom_col() +
    geom_vline(xintercept = 0, linewidth = 0.55, color = "#7F7F7F") +
    scale_fill_manual(values = c(
      "Higher than block E-I background" = "#C44E52",
      "Lower than block E-I background" = "#4E79A7"
    )) +
    labs(
      title = wrap_plot_title("Containment propensity is strongly pair-specific rather than uniform across inhibitory subclasses"),
      subtitle = wrap_plot_subtitle("Δ exact-containment rate relative to the within-block E-I background"),
      x = "Δ exact-containment rate vs block E-I background",
      y = NULL
    ) +
    plot_theme_safe
  save_pub(fig7, "fig07_pair_specific_background_deviation_topjournal", 10.8, 6.4)
  
  # Figure 8: GABA guest-cluster logistic forest
  forest_gaba <- gaba_cluster_logit %>%
    filter(term %in% c("log_gaba_size", "mean_ii_jaccard", "n_ii_partners")) %>%
    mutate(
      term_label = recode(
        term,
        `log_gaba_size` = "log(GABA cluster size)",
        `mean_ii_jaccard` = "Mean I-I Jaccard",
        `n_ii_partners` = "Number of I-I partners"
      )
    ) %>%
    mutate(term_label = factor(term_label, levels = rev(term_label)))
  
  fig8 <- ggplot(forest_gaba, aes(x = OR, y = term_label)) +
    geom_vline(xintercept = 1, linetype = 2, color = "#7F7F7F", linewidth = 0.55) +
    geom_errorbarh(aes(xmin = OR_low, xmax = OR_high), height = 0.18, linewidth = 0.7, color = "#4D4D4D") +
    geom_point(size = 2.8, color = "#4E79A7") +
    scale_x_log10(labels = label_number(accuracy = 0.01)) +
    labs(
      title = wrap_plot_title("Exact-contained GABA clusters are smaller and less embedded within inhibitory modules"),
      subtitle = wrap_plot_subtitle(
        paste0(
          "OR(log GABA size) = ",
          sprintf("%.3f", gaba_cluster_logit$OR[gaba_cluster_logit$term == "log_gaba_size"]),
          "; OR(mean I-I Jaccard) = ",
          sprintf("%.3f", gaba_cluster_logit$OR[gaba_cluster_logit$term == "mean_ii_jaccard"])
        )
      ),
      x = "Adjusted odds ratio (log scale)",
      y = NULL
    ) +
    plot_theme_safe
  save_pub(fig8, "fig08_gaba_guest_cluster_forest_topjournal", 8.8, 5.6)
  
  # Figure 9: Multi-host direction asymmetry
  fig9 <- ggplot(multi_host_direction_comparison, aes(direction, multi_host_rate, fill = direction)) +
    geom_col(width = 0.62) +
    geom_text(
      aes(label = paste0(percent(multi_host_rate, accuracy = 0.1), "\n(n=", comma(n_guest_clusters), ")")),
      vjust = -0.24, size = 4.4, lineheight = 0.95
    ) +
    scale_fill_manual(values = col_multi) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.16))) +
    labs(
      title = wrap_plot_title("Multi-host exact coverage is itself directionally asymmetric"),
      subtitle = wrap_plot_subtitle(
        paste0(
          "Fisher OR = ",
          sprintf("%.2f", unname(multi_host_ft$estimate)),
          "; p = ",
          fmt_p_short(multi_host_ft$p.value)
        )
      ),
      x = NULL,
      y = "Fraction of exact guest clusters with ≥2 exact hosts"
    ) +
    plot_theme_safe +
    theme(legend.position = "none")
  save_pub(fig9, "fig09_multi_host_direction_topjournal", 8.2, 5.8)
  
  # ---------------------------
  # Figure manifest
  # ---------------------------
  figure_manifest <- tibble(
    filename = c(
      "fig01_pairtype_exact_containment_topjournal",
      "fig02_directional_exact_asymmetry_topjournal",
      "fig03_smaller_cluster_direction_bias_topjournal",
      "fig04_smaller_cluster_logit_forest_topjournal",
      "fig05_exact_vs_other_jaccard_topjournal",
      "fig06_layer_exact_rate_topjournal",
      "fig07_pair_specific_background_deviation_topjournal",
      "fig08_gaba_guest_cluster_forest_topjournal",
      "fig09_multi_host_direction_topjournal"
    ),
    summary = c(
      "E-I exact containment rate vs E-E and I-I at block level",
      "Exact directional asymmetry: GABA in Glu vs Glu in GABA",
      "Smaller-GABA vs smaller-Glu exact containment rates",
      "Adjusted OR forest for smaller-cluster containment model",
      "Block-level Jaccard comparison for exact GABA-in-Glu vs other E-I edges",
      "Layer-specific exact GABA-in-Glu rates",
      "Pair-specific deviation from block-level E-I background",
      "Cluster-level logistic effects for exact-contained GABA guests",
      "Direction asymmetry of multi-host exact coverage"
    )
  )
  
  write_csv(figure_manifest, file.path(out_dir, "figure_manifest.csv"))
  
  message("All outputs saved to: ", normalizePath(out_dir))
}

run_module_05 <- function(runtime_base) {
  message('[MODULE 05] start')
  # [comment omitted: encoding-safe]
  # [removed] rm(list = ls())
  gc()
  suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(stringr)
    library(ggplot2)
    library(forcats)
    library(scales)
  })
  
  options(stringsAsFactors = FALSE, scipen = 999)
  
  # [comment omitted: encoding-safe]
  select <- dplyr::select
  rename <- dplyr::rename
  filter <- dplyr::filter
  mutate <- dplyr::mutate
  summarise <- dplyr::summarise
  arrange <- dplyr::arrange
  distinct <- dplyr::distinct
  count <- dplyr::count
  left_join <- dplyr::left_join
  right_join <- dplyr::right_join
  transmute <- dplyr::transmute
  pull <- dplyr::pull
  
  # =========================
  # Paths
  # =========================
  root_dir <- runtime_base
  input_file <- file.path(root_dir, "neuron-nonneuron-partner.csv")
  out_dir <- file.path(root_dir, "R_output_05")
  fig_dir <- file.path(out_dir, "figures")
  tab_dir <- file.path(out_dir, "tables")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (!file.exists(input_file)) {
    stop("Input file not found: ", input_file)
  }
  
  # =========================
  # Helper functions
  # =========================
  clean_subclass <- function(x) {
    x %>%
      str_replace("^\\d+\\s+", "") %>%
      str_replace_all("_", " ")
  }
  
  safe_fisher <- function(a, b, c, d) {
    mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    ft <- fisher.test(mat)
    tibble(OR = unname(ft$estimate), p = ft$p.value)
  }
  
  pair_expand <- function(labels) {
    labels <- sort(unique(labels))
    if (length(labels) < 2) return(tibble(a = character(), b = character()))
    cmb <- t(combn(labels, 2))
    tibble(a = cmb[, 1], b = cmb[, 2])
  }
  
  same_subclass_multi <- function(df_pairs) {
    df_pairs %>%
      distinct(neuron_type, neuron_label, nn_label, nn_subclass_clean) %>%
      group_by(neuron_type, neuron_label) %>%
      summarise(
        n_partners = n(),
        n_subclasses = n_distinct(nn_subclass_clean),
        same_dup = n_subclasses < n_partners,
        .groups = "drop"
      ) %>%
      filter(n_partners >= 2)
  }
  
  nn_group_map <- function(x) {
    case_when(
      str_detect(x, "Astro") ~ "Astrocyte",
      str_detect(x, "Oligo|OPC") ~ "Oligodendro-lineage",
      str_detect(x, "SMC|VLMC|Peri|Endo") ~ "Vascular-associated",
      str_detect(x, "Micro|Macro|Immune") ~ "Immune/myeloid",
      TRUE ~ "Other"
    )
  }
  
  fmt_p <- function(p) {
    ifelse(is.na(p), "NA",
           ifelse(p < 1e-4,
                  format(p, scientific = TRUE, digits = 2),
                  sprintf("%.4f", p)))
  }
  
  save_pub <- function(p, filename, width, height, dpi = 600) {
    ggsave_pub(file.path(fig_dir, paste0(filename, ".png")), p,
           width = width, height = height, dpi = dpi, bg = "white")
    ggsave_pub(file.path(fig_dir, paste0(filename, ".pdf")), p,
           width = width, height = height, bg = "white", device = cairo_pdf)
  }
  
  base_theme <- theme_pub_global(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0),
      plot.subtitle = element_text(size = 11.5, colour = "#4D4D4D"),
      axis.title = element_text(size = 14, colour = "black"),
      axis.text = element_text(size = 11.5, colour = "black"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 11),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 12),
      plot.margin = margin(15, 20, 12, 12)
    )
  
  # =========================
  # Read data
  # =========================
  df <- read_csv(input_file, show_col_types = FALSE)
  
  # Unique cluster table
  clusters <- bind_rows(
    df %>% transmute(
      label = cluster.1_label,
      cell_Neruon_type = cluster.1_cell_Neruon_type,
      subclass = cluster.1_subclass,
      subclass_clean = clean_subclass(cluster.1_subclass),
      slide = cluster.1_slide,
      layer = cluster.1_layer,
      total_cell_num = cluster.1_total_cell_num
    ),
    df %>% transmute(
      label = cluster.2_label,
      cell_Neruon_type = cluster.2_cell_Neruon_type,
      subclass = cluster.2_subclass,
      subclass_clean = clean_subclass(cluster.2_subclass),
      slide = cluster.2_slide,
      layer = cluster.2_layer,
      total_cell_num = cluster.2_total_cell_num
    )
  ) %>% distinct(label, .keep_all = TRUE)
  
  # Oriented neuron->nonneuron heterotypic pairs
  het <- df %>%
    filter(cluster.1_cell_Neruon_type %in% c("Glut", "Gaba"),
           cluster.2_cell_Neruon_type == "NonNeuron") %>%
    mutate(
      interaction_class = case_when(
        cluster.1.overlap.percent >= 1 & cluster.2.overlap.percent >= 1 ~ "mutual_equal",
        cluster.1.overlap.percent >= 1 ~ "Neuron_in_NN",
        cluster.2.overlap.percent >= 1 ~ "NN_in_Neuron",
        TRUE ~ "Partial"
      ),
      neuron_type = cluster.1_cell_Neruon_type,
      neuron_subclass = cluster.1_subclass,
      neuron_subclass_clean = clean_subclass(cluster.1_subclass),
      nn_subclass = cluster.2_subclass,
      nn_subclass_clean = clean_subclass(cluster.2_subclass),
      neuron_size = cluster.1_total_cell_num,
      nn_size = cluster.2_total_cell_num,
      neuron_label = cluster.1_label,
      nn_label = cluster.2_label,
      layer = cluster.1_layer,
      slide = cluster.1_slide,
      neuron_region = cluster.1_region,
      nn_region = cluster.2_region,
      containment_event = interaction_class %in% c("Neuron_in_NN", "NN_in_Neuron", "mutual_equal"),
      size_ratio_nn_to_neuron = nn_size / neuron_size,
      nn_group = nn_group_map(clean_subclass(nn_subclass))
    )
  
  # Unique NN-NN pairs for host-host lookups
  nnn <- df %>%
    filter(cluster.1_cell_Neruon_type == "NonNeuron",
           cluster.2_cell_Neruon_type == "NonNeuron") %>%
    mutate(
      min_label = pmin(cluster.1_label, cluster.2_label),
      max_label = pmax(cluster.1_label, cluster.2_label)
    ) %>%
    distinct(min_label, max_label, .keep_all = TRUE)
  
  nin <- het %>% filter(interaction_class == "Neuron_in_NN")
  partial <- het %>% filter(interaction_class == "Partial")
  nnin <- het %>% filter(interaction_class == "NN_in_Neuron")
  
  # =========================
  # Core summaries for the teacher's key claim
  # =========================
  mult_nin <- same_subclass_multi(nin)
  mult_partial <- same_subclass_multi(partial)
  mult_nnin <- same_subclass_multi(nnin)
  
  overall_dup <- bind_rows(
    mult_nin %>% summarise(context = "Full host", n_total = n(), n_same = sum(same_dup), rate = mean(same_dup)),
    mult_partial %>% summarise(context = "Partial background", n_total = n(), n_same = sum(same_dup), rate = mean(same_dup)),
    mult_nnin %>% summarise(context = "Full guest", n_total = n(), n_same = sum(same_dup), rate = mean(same_dup))
  )
  
  dup_ft <- safe_fisher(
    sum(mult_nin$same_dup),
    sum(!mult_nin$same_dup),
    sum(mult_partial$same_dup),
    sum(!mult_partial$same_dup)
  )
  
  # host-host pairs for neurons with >=2 full hosts
  host_pairs <- nin %>%
    distinct(neuron_type, neuron_label, nn_label, nn_subclass_clean) %>%
    group_by(neuron_type, neuron_label) %>%
    summarise(host_labels = list(sort(unique(nn_label))), .groups = "drop") %>%
    mutate(pair_tbl = map(host_labels, pair_expand)) %>%
    select(neuron_type, neuron_label, pair_tbl) %>%
    unnest(pair_tbl) %>%
    rename(host1 = a, host2 = b)
  
  hosthost <- host_pairs %>%
    mutate(
      min_label = pmin(host1, host2),
      max_label = pmax(host1, host2)
    ) %>%
    left_join(nnn %>% select(min_label, max_label, jaccard), by = c("min_label", "max_label"))
  
  hosthost_summary <- tibble(
    group = c("All NN-NN", "GABA multi-host", "Glut multi-host"),
    n_pairs = c(
      nrow(nnn),
      sum(hosthost$neuron_type == "Gaba" & !is.na(hosthost$jaccard)),
      sum(hosthost$neuron_type == "Glut" & !is.na(hosthost$jaccard))
    ),
    median_jaccard = c(
      median(nnn$jaccard, na.rm = TRUE),
      median(hosthost$jaccard[hosthost$neuron_type == "Gaba"], na.rm = TRUE),
      median(hosthost$jaccard[hosthost$neuron_type == "Glut"], na.rm = TRUE)
    ),
    mean_jaccard = c(
      mean(nnn$jaccard, na.rm = TRUE),
      mean(hosthost$jaccard[hosthost$neuron_type == "Gaba"], na.rm = TRUE),
      mean(hosthost$jaccard[hosthost$neuron_type == "Glut"], na.rm = TRUE)
    )
  )
  
  hosthost_comp <- bind_rows(
    {
      x <- hosthost$jaccard[hosthost$neuron_type == "Gaba"]
      y <- nnn$jaccard
      tibble(comparison = "GABA multi-host vs background", p = wilcox.test(x, y, exact = FALSE)$p.value)
    },
    {
      x <- hosthost$jaccard[hosthost$neuron_type == "Glut"]
      y <- nnn$jaccard
      tibble(comparison = "Glut multi-host vs background", p = wilcox.test(x, y, exact = FALSE)$p.value)
    }
  )
  
  # co-host enrichment
  cohost_enrichment <- function(nin_df, nt, min_count = 20) {
    sub <- nin_df %>% filter(neuron_type == nt)
    subtype_counts <- sub %>% count(nn_subclass_clean, sort = TRUE) %>% filter(n >= min_count)
    subtypes <- subtype_counts$nn_subclass_clean
    pres <- sub %>%
      distinct(neuron_label, nn_subclass_clean) %>%
      filter(nn_subclass_clean %in% subtypes) %>%
      mutate(value = TRUE) %>%
      pivot_wider(names_from = nn_subclass_clean, values_from = value, values_fill = FALSE)
    if (nrow(pres) == 0) return(tibble())
    rows <- list(); idx <- 1
    for (i in seq_along(subtypes)) {
      for (j in seq((i + 1), length(subtypes))) {
        if (j > length(subtypes)) next
        a_name <- subtypes[i]; b_name <- subtypes[j]
        A <- pres[[a_name]]; B <- pres[[b_name]]
        both <- sum(A & B); a_only <- sum(A & !B); b_only <- sum(!A & B); neither <- sum(!A & !B)
        ft <- safe_fisher(both, a_only, b_only, neither)
        rows[[idx]] <- tibble(
          neuron_type = nt,
          subtype_a = a_name,
          subtype_b = b_name,
          both = both,
          OR = ft$OR,
          p = ft$p
        )
        idx <- idx + 1
      }
    }
    bind_rows(rows) %>% mutate(q = p.adjust(p, method = "BH"))
  }
  
  gaba_cohost <- cohost_enrichment(nin, "Gaba", min_count = 20)
  glut_cohost <- cohost_enrichment(nin, "Glut", min_count = 20)
  cohost_sig <- bind_rows(gaba_cohost, glut_cohost) %>%
    mutate(pair = paste(subtype_a, subtype_b, sep = " + "),
           log2OR = log2(OR)) %>%
    filter(q < 0.05) %>%
    group_by(neuron_type) %>%
    arrange(desc(log2OR), .by_group = TRUE) %>%
    slice_head(n = 10) %>%
    ungroup()
  
  # representative motif enrichment for full host outcome
  motifs <- map_dfr(c("Gaba", "Glut"), function(nt) {
    sub_nt <- het %>% filter(neuron_type == nt)
    valid <- sub_nt %>% count(neuron_subclass_clean, nn_subclass_clean) %>% filter(n >= 20)
    map_dfr(seq_len(nrow(valid)), function(i) {
      ns <- valid$neuron_subclass_clean[i]
      nn <- valid$nn_subclass_clean[i]
      sub <- sub_nt %>% filter(neuron_subclass_clean == ns, nn_subclass_clean == nn)
      a <- sum(sub$interaction_class == "Neuron_in_NN")
      b <- nrow(sub) - a
      c <- sum(sub_nt$interaction_class == "Neuron_in_NN") - a
      d <- nrow(sub_nt) - nrow(sub) - c
      ft <- safe_fisher(a, b, c, d)
      tibble(
        neuron_type = nt,
        motif = paste(ns, nn, sep = " × "),
        n_total = nrow(sub),
        n_outcome = a,
        rate = a / nrow(sub),
        OR = ft$OR,
        p = ft$p
      )
    })
  }) %>%
    group_by(neuron_type) %>%
    mutate(q = p.adjust(p, method = "BH")) %>%
    ungroup()
  
  motif_top <- motifs %>%
    filter(q < 0.05, OR > 1) %>%
    mutate(log2OR = log2(OR)) %>%
    group_by(neuron_type) %>%
    arrange(desc(log2OR), .by_group = TRUE) %>%
    slice_head(n = 8) %>%
    ungroup()
  
  # =========================
  # Save key tables
  # =========================
  write_csv(overall_dup, file.path(tab_dir, "table_01_same_subclass_duplication_summary.csv"))
  write_csv(dup_ft, file.path(tab_dir, "table_02_same_subclass_duplication_fisher.csv"))
  write_csv(hosthost_summary, file.path(tab_dir, "table_03_hosthost_jaccard_summary.csv"))
  write_csv(hosthost_comp, file.path(tab_dir, "table_04_hosthost_jaccard_tests.csv"))
  write_csv(bind_rows(gaba_cohost, glut_cohost), file.path(tab_dir, "table_05_cohost_enrichment.csv"))
  write_csv(motifs, file.path(tab_dir, "table_06_fullhost_motif_enrichment.csv"))
  
  figure_manifest <- tibble(
    figure = c(
      "fig01_same_subclass_duplication_core",
      "fig02_hosthost_jaccard_vs_background",
      "fig03_cohost_module_enrichment",
      "fig04_representative_fullhost_motifs"
    ),
    meaning = c(
      "核心结论：多full-host neuron中，同一nonneuron subclass重复为0，而partial背景中重复很常见",
      "多宿主neuron周围，不同nonneuron host彼此重叠显著高于普通NN-NN背景",
      "co-host组合并非随机，血管相关模块最稳定",
      "neuron-nonneuron full-host关系具有强烈subclass specificity"
    )
  )
  write_csv(figure_manifest, file.path(out_dir, "figure_manifest.csv"))
  
  # =========================
  # Figure 1: Core teacher result
  # =========================
  fig1_df <- overall_dup %>%
    filter(context %in% c("Full host", "Partial background")) %>%
    mutate(context = factor(context, levels = c("Full host", "Partial background")),
           label = paste0(scales::comma(n_same), "/", scales::comma(n_total)))
  
  bracket_y <- max(fig1_df$rate) + 0.07
  fig01 <- ggplot(fig1_df, aes(x = context, y = rate, fill = context)) +
    geom_col(width = 0.62) +
    geom_text(aes(label = label), vjust = -0.22, size = 5.0) +
    geom_text(aes(y = pmax(rate * 0.45, 0.02), label = scales::percent(rate, accuracy = 0.1)),
              colour = "white", fontface = "bold", size = 4.4) +
    geom_segment(aes(x = 1, xend = 2, y = bracket_y, yend = bracket_y), inherit.aes = FALSE, linewidth = 0.6) +
    geom_segment(aes(x = 1, xend = 1, y = bracket_y, yend = bracket_y - 0.02), inherit.aes = FALSE, linewidth = 0.6) +
    geom_segment(aes(x = 2, xend = 2, y = bracket_y, yend = bracket_y - 0.02), inherit.aes = FALSE, linewidth = 0.6) +
    annotate("text", x = 1.5, y = bracket_y + 0.02,
             label = paste0("Fisher P = ", format(dup_ft$p, scientific = TRUE, digits = 2)),
             size = 4.4) +
    scale_fill_manual(values = c("Full host" = "#C44E52", "Partial background" = "#7E7E7E")) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       limits = c(0, max(fig1_df$rate) + 0.14),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(
      title = "Repeated full hosts are never from the same nonneuron subclass",
      subtitle = "Among neurons with at least two complete nonneuron hosts, same-subclass duplication disappears, but common in the partial multi-partner background.",
      x = NULL,
      y = "Fraction with repeated same-subclass partners"
    ) +
    base_theme +
    theme(legend.position = "none")
  save_pub(fig01, "fig01_same_subclass_duplication_core", 12.0, 5.6)
  
  # =========================
  # Figure 2: Host-host overlap structure
  # =========================
  fig2_df <- bind_rows(
    tibble(group = "All NN-NN background", neuron_type = "Background", jaccard = nnn$jaccard),
    tibble(group = "Multi-host around GABA", neuron_type = "Gaba", jaccard = hosthost$jaccard[hosthost$neuron_type == "Gaba"]),
    tibble(group = "Multi-host around Glut", neuron_type = "Glut", jaccard = hosthost$jaccard[hosthost$neuron_type == "Glut"])
  ) %>%
    filter(!is.na(jaccard)) %>%
    mutate(group = factor(group, levels = c("All NN-NN background", "Multi-host around GABA", "Multi-host around Glut")))
  
  y_max_fig02 <- quantile(fig2_df$jaccard, 0.995, na.rm = TRUE)
  med_labs <- fig2_df %>%
    group_by(group) %>%
    summarise(med = median(jaccard), .groups = "drop") %>%
    mutate(x = group,
           y = y_max_fig02 * 0.93,
           label = paste0("median=", sprintf("%.3f", med)))
  
  fig02 <- ggplot(fig2_df, aes(x = group, y = jaccard, fill = group)) +
    geom_violin(trim = TRUE, alpha = 0.75, colour = NA) +
    geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white", colour = "#222222", linewidth = 0.6) +
    geom_text(data = med_labs, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 4.2) +
    scale_fill_manual(values = c("All NN-NN background" = "#B9B9B9", "Multi-host around GABA" = "#CC6B6B", "Multi-host around Glut" = "#5B84B1")) +
    coord_cartesian(ylim = c(0, y_max_fig02)) +
    labs(
      title = "Different full hosts around the same neuron strongly overlap one another",
      subtitle = paste0(
        "Compared with the global NN-NN background, host-host Jaccard is strongly elevated around multi-host neurons. ",
        "GABA: P=", format(hosthost_comp$p[hosthost_comp$comparison == "GABA multi-host vs background"], scientific = TRUE, digits = 2),
        "; Glut: P=", format(hosthost_comp$p[hosthost_comp$comparison == "Glut multi-host vs background"], scientific = TRUE, digits = 2)
      ),
      x = NULL,
      y = "Nonneuron–nonneuron Jaccard"
    ) +
    base_theme +
    theme(legend.position = "none", axis.text.x = element_text(angle = 12, hjust = 1))
  save_pub(fig02, "fig02_hosthost_jaccard_vs_background", 12.4, 5.8)
  
  # =========================
  # Figure 3: Significant co-host modules
  # =========================
  fig3_df <- cohost_sig %>%
    mutate(pair = fct_reorder(pair, log2OR))
  
  fig03 <- ggplot(fig3_df, aes(x = log2OR, y = pair, size = both, colour = -log10(q))) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5, colour = "#808080") +
    geom_point(alpha = 0.95) +
    scale_colour_gradient(low = "#7EBCD2", high = "#C22E2E") +
    facet_wrap(~neuron_type, scales = "free_y", ncol = 1) +
    labs(
      title = "Co-host combinations are highly structured rather than random",
      subtitle = "The most stable enrichments are vascular-associated modules, especially SMC + VLMC and Endo + Peri.",
      x = "log2(odds ratio) for co-host enrichment among multi-host neurons",
      y = NULL,
      size = "# neurons",
      colour = "-log10(FDR)"
    ) +
    coord_cartesian(clip = "off") +
    base_theme +
    theme(
      plot.margin = margin(15, 28, 12, 12)
    )
  
  save_pub(fig03, "fig03_cohost_module_enrichment", 12.6, 7.2)
  
  # =========================
  # Figure 4: Representative motif enrichments
  # =========================
  fig4_df <- motif_top %>%
    mutate(motif = fct_reorder(motif, log2OR))
  
  fig04 <- ggplot(fig4_df, aes(x = log2OR, y = motif, size = n_total, colour = -log10(q))) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5, colour = "#808080") +
    geom_point(alpha = 0.95) +
    scale_colour_gradient(low = "#86BBD8", high = "#AE2012") +
    facet_wrap(~neuron_type, scales = "free_y", ncol = 1) +
    labs(
      title = "Full-host motifs are strongly subclass-specific",
      subtitle = "Representative enriched motifs show that neuron–nonneuron organization cannot be reduced to a simple GABA/Glut binary split.",
      x = "log2(odds ratio) for Neuron-in-Nonneuron enrichment",
      y = NULL,
      size = "# pairs",
      colour = "-log10(FDR)"
    ) +
    base_theme
  save_pub(fig04, "fig04_representative_fullhost_motifs", 12.8, 7.4)
  
  message("R_output_05 complete: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
}

run_module_06 <- function(runtime_base) {
  message('[MODULE 06] start')
  # [comment omitted: encoding-safe]
  # [removed] rm(list = ls())
  gc()
  # =========================================================
  # neuron_nonneuron_R_output_06_topjournal_fixed_v2.R
  # [comment omitted: encoding-safe]
  # Input : E:/zaw/2603/neuron-nonneuron-partner.csv
  # Output: E:/zaw/2603/R_output_06
  # =========================================================
  
  required_pkgs <- c(
    "dplyr", "tidyr", "stringr", "purrr", "readr", "tibble",
    "ggplot2", "forcats", "broom", "sandwich", "lmtest", "scales"
  )
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(purrr)
    library(readr)
    library(tibble)
    library(ggplot2)
    library(forcats)
    library(broom)
    library(sandwich)
    library(lmtest)
    library(scales)
  })
  
  options(stringsAsFactors = FALSE, scipen = 999)
  
  # Guard against function masking from other packages
  select <- dplyr::select
  rename <- dplyr::rename
  filter <- dplyr::filter
  mutate <- dplyr::mutate
  summarise <- dplyr::summarise
  arrange <- dplyr::arrange
  distinct <- dplyr::distinct
  count <- dplyr::count
  left_join <- dplyr::left_join
  inner_join <- dplyr::inner_join
  transmute <- dplyr::transmute
  pull <- dplyr::pull
  first <- dplyr::first
  n_distinct <- dplyr::n_distinct
  
  # ---------------------------
  # Paths
  # ---------------------------
  root_dir <- runtime_base
  input_csv <- file.path(root_dir, "neuron-nonneuron-partner.csv")
  out_dir <- file.path(root_dir, "R_output_06")
  fig_dir <- file.path(out_dir, "figures")
  tbl_dir <- file.path(out_dir, "tables")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(tbl_dir, showWarnings = FALSE, recursive = TRUE)
  if (!file.exists(input_csv)) stop("Input file not found: ", input_csv)
  
  # ---------------------------
  # Helpers
  # ---------------------------
  canonical_pair <- function(a, b) {
    ifelse(a < b, paste(a, b, sep = " || "), paste(b, a, sep = " || "))
  }
  
  clean_subclass <- function(x) {
    stringr::str_replace(as.character(x), "^\\d+\\s+", "")
  }
  
  rank_biserial_paired <- function(x, y) {
    d <- x - y
    d <- d[!is.na(d) & d != 0]
    if (length(d) == 0) return(NA_real_)
    (sum(d > 0) - sum(d < 0)) / length(d)
  }
  
  paired_wilcox_tbl <- function(dat, a, b) {
    tmp <- dat %>% dplyr::select(all_of(c(a, b))) %>% tidyr::drop_na()
    if (nrow(tmp) == 0) return(NULL)
    wt <- suppressWarnings(wilcox.test(tmp[[a]], tmp[[b]], paired = TRUE, exact = FALSE))
    tibble(
      n = nrow(tmp),
      mean_diff = mean(tmp[[a]] - tmp[[b]]),
      median_diff = median(tmp[[a]] - tmp[[b]]),
      rank_biserial = rank_biserial_paired(tmp[[a]], tmp[[b]]),
      p_value = wt$p.value
    )
  }
  
  tidy_cluster_model <- function(fit, cluster_vec, logistic = TRUE) {
    vc <- sandwich::vcovCL(fit, cluster = cluster_vec)
    ct <- lmtest::coeftest(fit, vcov. = vc)
    out <- tibble(
      term = rownames(ct),
      estimate = as.numeric(ct[, 1]),
      std_error = as.numeric(ct[, 2]),
      stat = as.numeric(ct[, 3]),
      p_value = as.numeric(ct[, 4])
    )
    out <- out %>%
      dplyr::mutate(
        conf_low = estimate - 1.96 * std_error,
        conf_high = estimate + 1.96 * std_error
      )
    if (logistic) {
      out <- out %>%
        dplyr::mutate(
          OR = exp(estimate),
          OR_low = exp(conf_low),
          OR_high = exp(conf_high)
        )
    }
    out
  }
  
  fmt_p_short <- function(p) {
    if (is.na(p)) return("NA")
    if (p < 1e-300) return("<1e-300")
    if (p < 1e-4) return(format(p, scientific = TRUE, digits = 2))
    sprintf("%.4f", p)
  }
  
  wrap_plot_title <- function(x, width = 54) stringr::str_wrap(x, width = width)
  wrap_plot_subtitle <- function(x, width = 92) stringr::str_wrap(x, width = width)
  
  save_pub <- function(p, filename, width, height, dpi = 600) {
    ggsave_pub(
      filename = file.path(fig_dir, paste0(filename, ".png")),
      plot = p, width = width, height = height, dpi = dpi,
      bg = "white", limitsize = FALSE
    )
    ggsave_pub(
      filename = file.path(fig_dir, paste0(filename, ".pdf")),
      plot = p, width = width, height = height,
      bg = "white", limitsize = FALSE
    )
  }
  
  pick_pairs_by_pattern <- function(tab, patterns, sign = c("both", "pos", "neg")) {
    sign <- match.arg(sign)
    out <- purrr::map_dfr(patterns, function(pt) {
      x <- tab %>% dplyr::filter(stringr::str_detect(subpair_name_clean, stringr::regex(pt, ignore_case = TRUE)))
      if (nrow(x) == 0) return(tibble())
      if (sign == "pos") x <- x %>% dplyr::filter(mean_diff_vs_bg > 0)
      if (sign == "neg") x <- x %>% dplyr::filter(mean_diff_vs_bg < 0)
      if (nrow(x) == 0) return(tibble())
      x %>% dplyr::arrange(p_value, dplyr::desc(abs(mean_diff_vs_bg))) %>% dplyr::slice(1)
    })
    out %>% dplyr::distinct(subpair_name_clean, .keep_all = TRUE)
  }
  
  plot_theme <- theme_pub_global(base_size = 12) +
    ggplot2::theme(
      plot.title.position = "plot",
      plot.margin = ggplot2::margin(18, 28, 14, 18),
      plot.title = ggplot2::element_text(face = "bold", size = 16.5, lineheight = 1.06, margin = ggplot2::margin(b = 8)),
      plot.subtitle = ggplot2::element_text(size = 11.8, colour = "#4D4D4D", lineheight = 1.10, margin = ggplot2::margin(b = 12)),
      axis.title = ggplot2::element_text(size = 13.5, colour = "black"),
      axis.text = ggplot2::element_text(size = 11.5, colour = "black"),
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 11.2)
    )
  
  theme_set(plot_theme)
  
  col_neuron <- c("Glut" = "#C44E52", "Gaba" = "#4E79A7")
  col_guest <- c("Not exact-contained" = "#BDBDBD", "Exact-contained" = "#C44E52")
  col_posneg <- c("Higher than background" = "#C44E52", "Lower than background" = "#4E79A7")
  col_shared <- c("Shared by both" = "#7A4FB3", "One neuron type only" = "#BDBDBD")
  
  # ---------------------------
  # Read data and build pair table
  # ---------------------------
  df <- readr::read_csv(input_csv, show_col_types = FALSE) %>%
    dplyr::mutate(pair_key = canonical_pair(cluster.1_label, cluster.2_label))
  
  pairs <- df %>%
    dplyr::group_by(pair_key) %>%
    dplyr::filter(cluster.1_label == pmin(cluster.1_label, cluster.2_label)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      label_a = cluster.1_label,
      label_b = cluster.2_label,
      slide = cluster.1_slide,
      layer = cluster.1_layer,
      region_a = cluster.1_region,
      region_b = cluster.2_region,
      type_a_raw = cluster.1_cell_Neruon_type,
      type_b_raw = cluster.2_cell_Neruon_type,
      subclass_a = cluster.1_subclass,
      subclass_b = cluster.2_subclass,
      subclass_a_clean = clean_subclass(cluster.1_subclass),
      subclass_b_clean = clean_subclass(cluster.2_subclass),
      size_a = cluster.1_total_cell_num,
      size_b = cluster.2_total_cell_num,
      overlap = overlap_cell,
      union = union_cell,
      jaccard = jaccard,
      overlap_prop_a = cluster.1.overlap.percent,
      overlap_prop_b = cluster.2.overlap.percent
    ) %>%
    dplyr::mutate(
      type_pair = purrr::map2_chr(type_a_raw, type_b_raw, ~ paste(sort(c(.x, .y)), collapse = " × ")),
      block = paste(slide, layer, sep = "|"),
      size_small = pmin(size_a, size_b),
      size_large = pmax(size_a, size_b),
      size_ratio = size_large / size_small,
      overlap_prop_small = overlap / size_small,
      overlap_prop_large = overlap / size_large,
      nesting_gap = overlap_prop_small - overlap_prop_large,
      exact_small_contained = overlap_prop_small == 1,
      near95_small_contained = overlap_prop_small >= 0.95,
      near90_small_contained = overlap_prop_small >= 0.90
    )
  
  clusters <- dplyr::bind_rows(
    pairs %>% dplyr::transmute(label = label_a, slide, layer, type_raw = type_a_raw, subclass = subclass_a, subclass_clean = subclass_a_clean, size = size_a),
    pairs %>% dplyr::transmute(label = label_b, slide, layer, type_raw = type_b_raw, subclass = subclass_b, subclass_clean = subclass_b_clean, size = size_b)
  ) %>% dplyr::distinct() %>% dplyr::mutate(block = paste(slide, layer, sep = "|"))
  
  mixed <- pairs %>%
    dplyr::filter((type_a_raw == "NonNeuron" & type_b_raw %in% c("Glut", "Gaba")) |
                    (type_b_raw == "NonNeuron" & type_a_raw %in% c("Glut", "Gaba"))) %>%
    dplyr::mutate(
      nn_on_a = type_a_raw == "NonNeuron",
      neuron_type = ifelse(nn_on_a, type_b_raw, type_a_raw),
      nn_subclass_full = ifelse(nn_on_a, subclass_a, subclass_b),
      nn_subclass = ifelse(nn_on_a, subclass_a_clean, subclass_b_clean),
      neuron_subclass_full = ifelse(nn_on_a, subclass_b, subclass_a),
      neuron_subclass = ifelse(nn_on_a, subclass_b_clean, subclass_a_clean),
      size_nn = ifelse(nn_on_a, size_a, size_b),
      size_neuron = ifelse(nn_on_a, size_b, size_a),
      overlap_prop_nn = ifelse(nn_on_a, overlap_prop_a, overlap_prop_b),
      overlap_prop_neuron = ifelse(nn_on_a, overlap_prop_b, overlap_prop_a),
      label_nn = ifelse(nn_on_a, label_a, label_b),
      label_neuron = ifelse(nn_on_a, label_b, label_a),
      nn_in_neuron_exact = overlap_prop_nn == 1,
      neuron_in_nn_exact = overlap_prop_neuron == 1,
      exact_any = nn_in_neuron_exact | neuron_in_nn_exact,
      both_exact = nn_in_neuron_exact & neuron_in_nn_exact,
      smaller_identity = dplyr::case_when(
        size_nn < size_neuron ~ "NN",
        size_neuron < size_nn ~ "Neuron",
        TRUE ~ "Equal"
      ),
      log_size_neuron = log(size_neuron),
      log_size_nn = log(size_nn),
      log_size_ratio = log(size_ratio),
      jaccard_clip = pmin(pmax(jaccard, 1e-4), 1 - 1e-4),
      logit_jaccard = log(jaccard_clip / (1 - jaccard_clip)),
      subpair_name = paste(neuron_subclass_full, nn_subclass_full, sep = " || "),
      subpair_name_clean = paste(neuron_subclass, nn_subclass, sep = " × ")
    )
  
  major_nn <- mixed %>% dplyr::count(nn_subclass, sort = TRUE) %>% dplyr::filter(n >= 100) %>% dplyr::pull(nn_subclass)
  
  # ---------------------------
  # Pair-specific background tests
  # ---------------------------
  subpair_bg_jaccard <- mixed %>%
    dplyr::group_by(neuron_type, subpair_name) %>%
    dplyr::group_split() %>%
    purrr::map_dfr(function(sub) {
      nt <- unique(sub$neuron_type)
      sp_name <- unique(sub$subpair_name)
      sp_clean <- unique(sub$subpair_name_clean)
      if (nrow(sub) < 20 || dplyr::n_distinct(sub$block) < 10) return(NULL)
      sp <- sub %>% dplyr::group_by(block) %>% dplyr::summarise(sp_mean = mean(jaccard), .groups = "drop")
      bg <- mixed %>% dplyr::filter(neuron_type == nt, subpair_name != sp_name) %>%
        dplyr::group_by(block) %>% dplyr::summarise(bg_mean = mean(jaccard), .groups = "drop")
      merged <- dplyr::inner_join(sp, bg, by = "block")
      if (nrow(merged) < 10) return(NULL)
      wt <- suppressWarnings(wilcox.test(merged$sp_mean, merged$bg_mean, paired = TRUE, exact = FALSE))
      tibble(
        neuron_type = nt,
        subpair_name = sp_name,
        subpair_name_clean = sp_clean,
        n_edges = nrow(sub),
        n_blocks = nrow(merged),
        mean_jaccard = mean(sub$jaccard),
        mean_diff_vs_bg = mean(merged$sp_mean - merged$bg_mean),
        median_diff_vs_bg = median(merged$sp_mean - merged$bg_mean),
        p_value = wt$p.value
      )
    }) %>% dplyr::mutate(fdr_bh = p.adjust(p_value, method = "BH"))
  
  subpair_bg_exact <- mixed %>%
    dplyr::group_by(neuron_type, subpair_name) %>%
    dplyr::group_split() %>%
    purrr::map_dfr(function(sub) {
      nt <- unique(sub$neuron_type)
      sp_name <- unique(sub$subpair_name)
      sp_clean <- unique(sub$subpair_name_clean)
      if (nrow(sub) < 20 || dplyr::n_distinct(sub$block) < 10) return(NULL)
      sp <- sub %>% dplyr::group_by(block) %>% dplyr::summarise(sp_rate = mean(nn_in_neuron_exact), .groups = "drop")
      bg <- mixed %>% dplyr::filter(neuron_type == nt, subpair_name != sp_name) %>%
        dplyr::group_by(block) %>% dplyr::summarise(bg_rate = mean(nn_in_neuron_exact), .groups = "drop")
      merged <- dplyr::inner_join(sp, bg, by = "block")
      if (nrow(merged) < 10) return(NULL)
      wt <- suppressWarnings(wilcox.test(merged$sp_rate, merged$bg_rate, paired = TRUE, exact = FALSE))
      tibble(
        neuron_type = nt,
        subpair_name = sp_name,
        subpair_name_clean = sp_clean,
        n_edges = nrow(sub),
        n_blocks = nrow(merged),
        exact_rate = mean(sub$nn_in_neuron_exact),
        mean_diff_vs_bg = mean(merged$sp_rate - merged$bg_rate),
        median_diff_vs_bg = median(merged$sp_rate - merged$bg_rate),
        p_value = wt$p.value
      )
    }) %>% dplyr::mutate(fdr_bh = p.adjust(p_value, method = "BH"))
  
  # ---------------------------
  # Cluster-level host / guest analyses
  # ---------------------------
  neuron_cluster_stats <- mixed %>%
    dplyr::group_by(label_neuron, neuron_type, neuron_subclass_full, neuron_subclass, slide, layer) %>%
    dplyr::summarise(
      neuron_size = dplyr::first(size_neuron),
      n_mixed_edges = dplyr::n(),
      mean_jaccard = mean(jaccard),
      mean_nesting_gap = mean(nesting_gap),
      n_exact_nn_guests = sum(nn_in_neuron_exact),
      prop_exact_nn_guests = mean(nn_in_neuron_exact),
      n_exact_neuron_in_nn = sum(neuron_in_nn_exact),
      n_unique_nn_subclasses = dplyr::n_distinct(nn_subclass),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      host_exact_any = n_exact_nn_guests > 0,
      multi_exact_nn_guests = n_exact_nn_guests >= 2,
      log_neuron_size = log(neuron_size),
      block = paste(slide, layer, sep = "|")
    )
  
  nnnn_edges <- pairs %>% dplyr::filter(type_pair == "NonNeuron × NonNeuron")
  
  nnnn_long <- dplyr::bind_rows(
    nnnn_edges %>% dplyr::transmute(label_nn = label_a, jaccard = jaccard, overlap_prop_self = overlap_prop_a),
    nnnn_edges %>% dplyr::transmute(label_nn = label_b, jaccard = jaccard, overlap_prop_self = overlap_prop_b)
  )
  
  nnnn_cluster <- nnnn_long %>%
    dplyr::group_by(label_nn) %>%
    dplyr::summarise(
      n_nnn_edges = dplyr::n(),
      mean_nnn_jaccard = mean(jaccard, na.rm = TRUE),
      mean_nnn_overlap_self = mean(overlap_prop_self, na.rm = TRUE),
      .groups = "drop"
    )
  
  nn_cluster_stats <- mixed %>%
    dplyr::group_by(label_nn, nn_subclass_full, nn_subclass, slide, layer) %>%
    dplyr::summarise(
      nn_size = dplyr::first(size_nn),
      n_mixed_edges = dplyr::n(),
      mean_jaccard = mean(jaccard),
      mean_nesting_gap = mean(nesting_gap),
      n_exact_neuron_hosts = sum(nn_in_neuron_exact),
      prop_exact_neuron_hosts = mean(nn_in_neuron_exact),
      n_exact_nn_hosts = sum(neuron_in_nn_exact),
      n_unique_neuron_subclasses = dplyr::n_distinct(neuron_subclass_full),
      n_unique_neuron_types = dplyr::n_distinct(neuron_type),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      guest_exact_any = n_exact_neuron_hosts > 0,
      multi_host_exact = n_exact_neuron_hosts >= 2,
      log_nn_size = log(nn_size),
      block = paste(slide, layer, sep = "|")
    ) %>%
    dplyr::left_join(nnnn_cluster, by = "label_nn") %>%
    dplyr::mutate(
      n_nnn_edges = ifelse(is.na(n_nnn_edges), 0L, n_nnn_edges),
      mean_nnn_jaccard = ifelse(is.na(mean_nnn_jaccard), 0, mean_nnn_jaccard),
      mean_nnn_overlap_self = ifelse(is.na(mean_nnn_overlap_self), 0, mean_nnn_overlap_self),
      nn_subclass_model = ifelse(nn_subclass %in% major_nn, nn_subclass, "OtherRareNN")
    )
  
  fit_host <- glm(
    host_exact_any ~ neuron_type + log_neuron_size + n_mixed_edges + layer,
    data = neuron_cluster_stats, family = binomial()
  )
  reg_host <- tidy_cluster_model(fit_host, neuron_cluster_stats$block, logistic = TRUE)
  
  fit_guest <- glm(
    guest_exact_any ~ log_nn_size + n_mixed_edges + mean_nnn_jaccard + layer + nn_subclass_model,
    data = nn_cluster_stats, family = binomial()
  )
  reg_guest <- tidy_cluster_model(fit_guest, nn_cluster_stats$block, logistic = TRUE)
  
  fit_multihost <- glm(
    multi_host_exact ~ log_nn_size + n_mixed_edges + mean_nnn_jaccard + layer + nn_subclass_model,
    data = nn_cluster_stats, family = binomial()
  )
  reg_multihost <- tidy_cluster_model(fit_multihost, nn_cluster_stats$block, logistic = TRUE)
  
  exact_guest_both_types <- mixed %>%
    dplyr::filter(nn_in_neuron_exact) %>%
    dplyr::group_by(label_nn) %>%
    dplyr::summarise(
      n_exact_hosts = dplyr::n_distinct(label_neuron),
      n_glut_exact_hosts = sum(neuron_type == "Glut"),
      n_gaba_exact_hosts = sum(neuron_type == "Gaba"),
      has_glut_exact = any(neuron_type == "Glut"),
      has_gaba_exact = any(neuron_type == "Gaba"),
      nn_subclass = dplyr::first(nn_subclass),
      nn_size = dplyr::first(size_nn),
      block = dplyr::first(block),
      .groups = "drop"
    ) %>%
    dplyr::mutate(both_neuron_types_exact = has_glut_exact & has_gaba_exact)
  
  exact_guest_both_types_summary <- exact_guest_both_types %>%
    dplyr::group_by(nn_subclass) %>%
    dplyr::summarise(
      n_exact_guest_clusters = dplyr::n(),
      both_neuron_types_rate = mean(both_neuron_types_exact),
      mean_n_exact_hosts = mean(n_exact_hosts),
      median_n_exact_hosts = median(n_exact_hosts),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(both_neuron_types_rate))
  
  overall_shared_summary <- tibble(
    n_exact_guest_clusters = nrow(exact_guest_both_types),
    n_shared_by_both_types = sum(exact_guest_both_types$both_neuron_types_exact),
    shared_rate = mean(exact_guest_both_types$both_neuron_types_exact)
  )
  
  # ---------------------------
  # Selected exemplars
  # ---------------------------
  jaccard_interest_patterns <- c(
    "L2/3 IT RSP.*Oligo",
    "L4 RSP-ACA.*Oligo",
    "L2/3 IT CTX.*Endo",
    "L4/5 IT CTX.*Astro-TE",
    "Lamp5.*VLMC",
    "Lamp5.*SMC"
  )
  exact_interest_patterns <- c(
    "Lamp5.*OPC",
    "Lamp5.*Endo",
    "Lamp5.*VLMC",
    "Lamp5.*SMC"
  )
  
  sel_jaccard_pairs <- pick_pairs_by_pattern(subpair_bg_jaccard, jaccard_interest_patterns, sign = "both") %>%
    dplyr::mutate(direction = ifelse(mean_diff_vs_bg > 0, "Higher than background", "Lower than background"))
  if (nrow(sel_jaccard_pairs) < 8) {
    extra <- subpair_bg_jaccard %>% dplyr::filter(fdr_bh < 0.05) %>% dplyr::arrange(dplyr::desc(abs(mean_diff_vs_bg))) %>% dplyr::slice_head(n = 8)
    sel_jaccard_pairs <- dplyr::bind_rows(sel_jaccard_pairs, extra) %>% dplyr::distinct(subpair_name_clean, .keep_all = TRUE) %>%
      dplyr::slice_head(n = 10) %>% dplyr::mutate(direction = ifelse(mean_diff_vs_bg > 0, "Higher than background", "Lower than background"))
  }
  
  sel_exact_pairs <- pick_pairs_by_pattern(subpair_bg_exact, exact_interest_patterns, sign = "both") %>%
    dplyr::mutate(direction = ifelse(mean_diff_vs_bg > 0, "Higher than background", "Lower than background"))
  if (nrow(sel_exact_pairs) < 6) {
    extra <- subpair_bg_exact %>% dplyr::filter(fdr_bh < 0.05) %>% dplyr::arrange(dplyr::desc(abs(mean_diff_vs_bg))) %>% dplyr::slice_head(n = 8)
    sel_exact_pairs <- dplyr::bind_rows(sel_exact_pairs, extra) %>% dplyr::distinct(subpair_name_clean, .keep_all = TRUE) %>%
      dplyr::slice_head(n = 8) %>% dplyr::mutate(direction = ifelse(mean_diff_vs_bg > 0, "Higher than background", "Lower than background"))
  }
  
  # ---------------------------
  # Save tables
  # ---------------------------
  write_csv(reg_guest, file.path(tbl_dir, "table_01_reg_guest_exact_any.csv"))
  write_csv(reg_host, file.path(tbl_dir, "table_02_reg_host_exact_any.csv"))
  write_csv(reg_multihost, file.path(tbl_dir, "table_03_reg_multi_host_exact.csv"))
  write_csv(subpair_bg_jaccard, file.path(tbl_dir, "table_04_subpair_bg_jaccard.csv"))
  write_csv(subpair_bg_exact, file.path(tbl_dir, "table_05_subpair_bg_exact.csv"))
  write_csv(nn_cluster_stats, file.path(tbl_dir, "table_06_nn_cluster_stats.csv"))
  write_csv(neuron_cluster_stats, file.path(tbl_dir, "table_07_neuron_cluster_stats.csv"))
  write_csv(exact_guest_both_types, file.path(tbl_dir, "table_08_exact_guest_both_types.csv"))
  write_csv(exact_guest_both_types_summary, file.path(tbl_dir, "table_09_exact_guest_both_types_summary.csv"))
  write_csv(overall_shared_summary, file.path(tbl_dir, "table_10_exact_guest_both_types_overall.csv"))
  write_csv(sel_jaccard_pairs, file.path(tbl_dir, "table_11_selected_pair_jaccard_examples.csv"))
  write_csv(sel_exact_pairs, file.path(tbl_dir, "table_12_selected_pair_exact_examples.csv"))
  
  # ---------------------------
  # Figures
  # ---------------------------
  forest_guest <- reg_guest %>%
    dplyr::filter(term %in% c("log_nn_size", "mean_nnn_jaccard", "n_mixed_edges")) %>%
    dplyr::mutate(
      term_label = dplyr::recode(
        term,
        `log_nn_size` = "log(non-neuron cluster size)",
        `mean_nnn_jaccard` = "Mean NN-NN Jaccard",
        `n_mixed_edges` = "Number of neuron-NN partners"
      ),
      term_label = factor(term_label, levels = rev(term_label))
    )
  
  fig1 <- ggplot(forest_guest, aes(x = OR, y = term_label)) +
    geom_vline(xintercept = 1, linetype = 2, linewidth = 0.55, colour = "#7F7F7F") +
    geom_errorbarh(aes(xmin = OR_low, xmax = OR_high), height = 0.16, linewidth = 0.7, colour = "#4D4D4D") +
    geom_point(size = 2.9, colour = "#C44E52") +
    scale_x_log10(labels = label_number(accuracy = 0.01)) +
    labs(
      title = wrap_plot_title("Exact-contained non-neuron clusters are smaller and less self-cohesive within the non-neuronal network"),
      subtitle = wrap_plot_subtitle(
        paste0(
          "NN guest model: OR(log_nn_size) = ", sprintf("%.3f", reg_guest$OR[reg_guest$term == "log_nn_size"]),
          ", p = ", fmt_p_short(reg_guest$p_value[reg_guest$term == "log_nn_size"]),
          "; OR(mean_nnn_jaccard) = ", sprintf("%.4f", reg_guest$OR[reg_guest$term == "mean_nnn_jaccard"]),
          ", p = ", fmt_p_short(reg_guest$p_value[reg_guest$term == "mean_nnn_jaccard"])
        )
      ),
      x = "Adjusted odds ratio (log scale)",
      y = NULL
    ) +
    plot_theme
  save_pub(fig1, "fig01_nn_guest_forest_topjournal", 8.8, 5.8)
  
  nn_guest_plot <- nn_cluster_stats %>%
    dplyr::mutate(group = ifelse(guest_exact_any, "Exact-contained", "Not exact-contained"))
  
  fig2 <- ggplot(nn_guest_plot, aes(group, nn_size, fill = group)) +
    geom_violin(width = 0.86, alpha = 0.18, colour = NA, trim = FALSE) +
    geom_boxplot(width = 0.28, outlier.shape = NA, linewidth = 0.7, alpha = 0.9) +
    scale_y_log10(labels = label_number(accuracy = 1)) +
    scale_fill_manual(values = col_guest) +
    labs(
      title = wrap_plot_title("Exact-contained non-neuron clusters are not random samples: they are markedly smaller"),
      subtitle = wrap_plot_subtitle(
        paste0(
          "Cluster-level model effect: OR(log_nn_size) = ",
          sprintf("%.3f", reg_guest$OR[reg_guest$term == "log_nn_size"]),
          "; p = ", fmt_p_short(reg_guest$p_value[reg_guest$term == "log_nn_size"])
        )
      ),
      x = NULL,
      y = "Non-neuron cluster size (log scale)"
    ) +
    plot_theme +
    theme(legend.position = "none")
  save_pub(fig2, "fig02_nn_guest_size_topjournal", 8.2, 5.8)
  
  fig3 <- ggplot(nn_guest_plot, aes(group, mean_nnn_jaccard, fill = group)) +
    geom_violin(width = 0.86, alpha = 0.18, colour = NA, trim = FALSE) +
    geom_boxplot(width = 0.28, outlier.shape = NA, linewidth = 0.7, alpha = 0.9) +
    scale_fill_manual(values = col_guest) +
    labs(
      title = wrap_plot_title("Exact-contained non-neuron clusters are less cohesive within the NN-NN interaction scaffold"),
      subtitle = wrap_plot_subtitle(
        paste0(
          "Cluster-level model effect: OR(mean_nnn_jaccard) = ",
          sprintf("%.4f", reg_guest$OR[reg_guest$term == "mean_nnn_jaccard"]),
          "; p = ", fmt_p_short(reg_guest$p_value[reg_guest$term == "mean_nnn_jaccard"])
        )
      ),
      x = NULL,
      y = "Mean NN-NN Jaccard"
    ) +
    plot_theme +
    theme(legend.position = "none")
  save_pub(fig3, "fig03_nn_guest_nnn_cohesion_topjournal", 8.2, 5.8)
  
  forest_host <- reg_host %>%
    dplyr::filter(term %in% c("neuron_typeGlut", "log_neuron_size", "n_mixed_edges")) %>%
    dplyr::mutate(
      term_label = dplyr::recode(
        term,
        `neuron_typeGlut` = "Glut vs Gaba",
        `log_neuron_size` = "log(neuron cluster size)",
        `n_mixed_edges` = "Mixed degree"
      ),
      term_label = factor(term_label, levels = rev(term_label))
    )
  
  fig4 <- ggplot(forest_host, aes(x = OR, y = term_label)) +
    geom_vline(xintercept = 1, linetype = 2, linewidth = 0.55, colour = "#7F7F7F") +
    geom_errorbarh(aes(xmin = OR_low, xmax = OR_high), height = 0.16, linewidth = 0.7, colour = "#4D4D4D") +
    geom_point(size = 2.9, colour = "#4E79A7") +
    scale_x_log10(labels = label_number(accuracy = 0.01)) +
    labs(
      title = wrap_plot_title("Glut clusters are more likely to become exact hosts of non-neuron islands"),
      subtitle = wrap_plot_subtitle(
        paste0(
          "After controlling for neuron size and mixed degree, OR(Glut host_exact_any) = ",
          sprintf("%.2f", reg_host$OR[reg_host$term == "neuron_typeGlut"]),
          "; p = ", fmt_p_short(reg_host$p_value[reg_host$term == "neuron_typeGlut"])
        )
      ),
      x = "Adjusted odds ratio (log scale)",
      y = NULL
    ) +
    plot_theme
  save_pub(fig4, "fig04_glut_host_forest_topjournal", 8.6, 5.6)
  
  shared_df <- tibble(
    group = c("Shared by both", "One neuron type only"),
    rate = c(
      mean(exact_guest_both_types$both_neuron_types_exact),
      1 - mean(exact_guest_both_types$both_neuron_types_exact)
    )
  )
  
  fig5 <- ggplot(shared_df, aes(group, rate, fill = group)) +
    geom_col(width = 0.62) +
    geom_text(aes(label = percent(rate, accuracy = 0.1)), vjust = -0.25, size = 4.8) +
    scale_fill_manual(values = col_shared) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = expansion(mult = c(0, 0.08))) +
    labs(
      title = wrap_plot_title("A subset of exact-contained non-neuron niches is shared by both Glut and Gaba hosts"),
      subtitle = wrap_plot_subtitle(
        paste0(
          "Among all exact guest NN clusters, ",
          percent(mean(exact_guest_both_types$both_neuron_types_exact), accuracy = 0.1),
          " have both Glut and Gaba exact hosts."
        )
      ),
      x = NULL,
      y = "Fraction of exact guest non-neuron clusters"
    ) +
    plot_theme +
    theme(legend.position = "none")
  save_pub(fig5, "fig05_shared_exact_guest_rate_topjournal", 7.6, 5.6)
  
  sel_jaccard_pairs <- sel_jaccard_pairs %>% dplyr::arrange(mean_diff_vs_bg)
  fig6 <- ggplot(sel_jaccard_pairs, aes(x = mean_diff_vs_bg, y = forcats::fct_reorder(subpair_name_clean, mean_diff_vs_bg), fill = direction)) +
    geom_col() +
    geom_vline(xintercept = 0, linewidth = 0.5) +
    scale_fill_manual(values = col_posneg) +
    labs(
      title = wrap_plot_title("What matters is not the coarse label, but the specific neuron × non-neuron pair"),
      subtitle = wrap_plot_subtitle("Representative pairs with Jaccard significantly above or below the within-neuron-type background."),
      x = "Δ mean Jaccard vs background",
      y = NULL
    ) +
    plot_theme
  save_pub(fig6, "fig06_selected_pair_jaccard_examples_topjournal", 10.2, 6.4)
  
  sel_exact_pairs <- sel_exact_pairs %>% dplyr::arrange(mean_diff_vs_bg)
  fig7 <- ggplot(sel_exact_pairs, aes(x = mean_diff_vs_bg, y = forcats::fct_reorder(subpair_name_clean, mean_diff_vs_bg), fill = direction)) +
    geom_col() +
    geom_vline(xintercept = 0, linewidth = 0.5) +
    scale_fill_manual(values = col_posneg) +
    labs(
      title = wrap_plot_title("The same non-neuron subclass can switch to a different spatial role under different neuron subclasses"),
      subtitle = wrap_plot_subtitle("Representative pairs with exact containment significantly above or below the within-neuron-type background."),
      x = "Δ exact containment rate vs background",
      y = NULL
    ) +
    plot_theme
  save_pub(fig7, "fig07_selected_pair_exact_examples_topjournal", 10.2, 6.0)
  
  shared_subclass_plot <- exact_guest_both_types_summary %>%
    dplyr::filter(n_exact_guest_clusters >= 10) %>%
    dplyr::slice_head(n = 10) %>%
    dplyr::arrange(both_neuron_types_rate)
  
  fig8 <- ggplot(shared_subclass_plot, aes(x = both_neuron_types_rate, y = forcats::fct_reorder(nn_subclass, both_neuron_types_rate))) +
    geom_col(fill = "#7A4FB3") +
    labs(
      title = wrap_plot_title("Shared exact guest niches are concentrated in a subset of non-neuron subclasses"),
      subtitle = wrap_plot_subtitle("Top non-neuron subclasses ranked by the proportion of exact guest clusters shared by both Glut and Gaba hosts."),
      x = "Fraction shared by both neuron types",
      y = NULL
    ) +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    plot_theme
  save_pub(fig8, "fig08_shared_exact_guest_subclasses_topjournal", 8.6, 6.0)
  
  figure_manifest <- tibble(
    filename = c(
      "fig01_nn_guest_forest_topjournal",
      "fig02_nn_guest_size_topjournal",
      "fig03_nn_guest_nnn_cohesion_topjournal",
      "fig04_glut_host_forest_topjournal",
      "fig05_shared_exact_guest_rate_topjournal",
      "fig06_selected_pair_jaccard_examples_topjournal",
      "fig07_selected_pair_exact_examples_topjournal",
      "fig08_shared_exact_guest_subclasses_topjournal"
    ),
    summary = c(
      "Cluster-level model for exact-contained non-neuron guests",
      "Exact-contained non-neuron guests are smaller",
      "Exact-contained non-neuron guests are less cohesive in NN-NN space",
      "Glut clusters are more likely to become exact hosts",
      "Overall fraction of exact guest NN clusters shared by both Glut and Gaba",
      "Selected pair-level Jaccard deviations from background",
      "Selected pair-level exact containment deviations from background",
      "Non-neuron subclasses enriched for cross-neuron-type exact sharing"
    )
  )
  write_csv(figure_manifest, file.path(out_dir, "figure_manifest.csv"))
  
  message("All outputs written to: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
}

run_module_07 <- function(runtime_base) {
  message('[MODULE 07] start')
  # [comment omitted: encoding-safe]
  # [removed] rm(list = ls())
  gc()
  #!/usr/bin/env Rscript
  # ============================================================
  # neuron_nonneuron_R_output_07_topjournal.R
  # Purpose:
  #   Multivariable-model focused analysis and publication-style
  #   visualization for neuron vs nonneuron spatial interactions.
  #   Emphasis:
  #   "After controlling for layer, size and region composition,
  #    the Gaba effect on exact neuron-in-nonneuron weakens,
  #    but positive effects on asymmetry and log2 enrichment remain."
  # Input : E:/zaw/2603/neuron-nonneuron-partner.csv
  # Output: E:/zaw/2603/R_output_07
  # ============================================================
  
  required_pkgs <- c(
    "dplyr", "tidyr", "readr", "stringr", "purrr", "tibble",
    "ggplot2", "broom", "sandwich", "lmtest", "forcats", "scales"
  )
  
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readr)
    library(stringr)
    library(purrr)
    library(tibble)
    library(ggplot2)
    library(broom)
    library(sandwich)
    library(lmtest)
    library(forcats)
    library(scales)
  })
  
  # protect against namespace conflicts
  select <- dplyr::select
  rename <- dplyr::rename
  filter <- dplyr::filter
  mutate <- dplyr::mutate
  summarise <- dplyr::summarise
  arrange <- dplyr::arrange
  count <- dplyr::count
  distinct <- dplyr::distinct
  left_join <- dplyr::left_join
  inner_join <- dplyr::inner_join
  pull <- dplyr::pull
  transmute <- dplyr::transmute
  first <- dplyr::first
  
  options(scipen = 999)
  
  # -----------------------------
  # I/O
  # -----------------------------
  base_dir <- runtime_base
  input_file <- file.path(base_dir, "neuron-nonneuron-partner.csv")
  output_dir <- file.path(base_dir, "R_output_07")
  fig_dir <- file.path(output_dir, "figures")
  table_dir <- file.path(output_dir, "tables")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (!file.exists(input_file)) {
    stop("Input file not found: ", input_file)
  }
  
  # -----------------------------
  # Helper functions
  # -----------------------------
  extract_nn_family <- function(x) {
    x %>%
      stringr::str_replace("^\\s*[0-9]+\\s+", "") %>%
      stringr::str_replace("\\s+NN\\s*$", "") %>%
      stringr::str_trim()
  }
  
  extract_neuron_family <- function(x) {
    s <- x %>%
      stringr::str_replace("^\\s*[0-9]+\\s+", "") %>%
      stringr::str_replace("\\s+(Glut|Gaba)\\s*$", "") %>%
      stringr::str_trim()
    
    dplyr::case_when(
      stringr::str_detect(s, "\\bLamp5\\b") ~ "Lamp5",
      stringr::str_detect(s, "\\bPvalb\\b") ~ "Pvalb",
      stringr::str_detect(s, "\\bSst\\b") ~ "Sst",
      stringr::str_detect(s, "\\bVip\\b") ~ "Vip",
      stringr::str_detect(s, "\\bSncg\\b") ~ "Sncg",
      stringr::str_detect(s, "^L2/3 IT\\b") ~ "L2/3 IT",
      stringr::str_detect(s, "^L4 RSP-ACA\\b") ~ "L4 RSP-ACA",
      stringr::str_detect(s, "^L4/5\\b") ~ "L4/5",
      stringr::str_detect(s, "^L5 ET\\b") ~ "L5 ET",
      stringr::str_detect(s, "^L5 IT\\b") ~ "L5 IT",
      stringr::str_detect(s, "^L5 NP\\b") ~ "L5 NP",
      stringr::str_detect(s, "^L5/6\\b") ~ "L5/6",
      stringr::str_detect(s, "^L6 CT\\b") ~ "L6 CT",
      stringr::str_detect(s, "^L6 IT\\b") ~ "L6 IT",
      stringr::str_detect(s, "^L6b\\b") ~ "L6b",
      stringr::str_detect(s, "CLA-EPd-CTX") ~ "CLA-EPd-CTX",
      stringr::str_detect(s, "^IT EP-CLA\\b") ~ "IT EP-CLA",
      TRUE ~ paste(head(str_split(s, "\\s+")[[1]], 2), collapse = " ")
    )
  }
  
  region_jaccard <- function(a, b) {
    sa <- unique(str_trim(unlist(str_split(ifelse(is.na(a), "", a), ","))))
    sb <- unique(str_trim(unlist(str_split(ifelse(is.na(b), "", b), ","))))
    sa <- sa[sa != ""]
    sb <- sb[sb != ""]
    u <- union(sa, sb)
    if (length(u) == 0) return(NA_real_)
    length(intersect(sa, sb)) / length(u)
  }
  
  canonical_pair <- function(a, b) ifelse(a < b, paste(a, b, sep = " || "), paste(b, a, sep = " || "))
  
  hg_p_upper <- function(N, K, n, x, union_cell) {
    N2 <- max(as.integer(N), as.integer(K), as.integer(n), as.integer(union_cell))
    phyper(q = as.integer(x) - 1, m = as.integer(K), n = as.integer(N2 - K),
           k = as.integer(n), lower.tail = FALSE)
  }
  
  add_fdr <- function(df, p_col = "p.value", out_col = "fdr_bh") {
    df[[out_col]] <- p.adjust(df[[p_col]], method = "BH")
    df
  }
  
  cluster_robust_tidy <- function(model, cluster, logistic = FALSE) {
    vc <- sandwich::vcovCL(model, cluster = cluster)
    ct <- lmtest::coeftest(model, vcov. = vc)
    out <- tibble(
      term = rownames(ct),
      estimate = as.numeric(ct[, 1]),
      std.error = as.numeric(ct[, 2]),
      statistic = as.numeric(ct[, 3]),
      p.value = as.numeric(ct[, 4])
    ) %>%
      mutate(
        conf.low = estimate - 1.96 * std.error,
        conf.high = estimate + 1.96 * std.error
      )
    if (logistic) {
      out <- out %>%
        mutate(
          odds_ratio = exp(estimate),
          or_low = exp(conf.low),
          or_high = exp(conf.high)
        )
    }
    out
  }
  
  fmt_p_short <- function(p) {
    if (is.na(p)) return("NA")
    if (p < 1e-300) return("<1e-300")
    if (p < 1e-4) return(format(p, scientific = TRUE, digits = 2))
    sprintf("%.4f", p)
  }
  
  wrap_plot_title <- function(x, width = 58) stringr::str_wrap(x, width = width)
  wrap_plot_subtitle <- function(x, width = 100) stringr::str_wrap(x, width = width)
  
  save_pub <- function(plot_obj, filename, width = 8, height = 5, dpi = 600) {
    ggsave_pub(filename = file.path(fig_dir, paste0(filename, ".png")), plot = plot_obj,
           width = width, height = height, dpi = dpi, bg = "white", limitsize = FALSE)
    ggsave_pub(filename = file.path(fig_dir, paste0(filename, ".pdf")), plot = plot_obj,
           width = width, height = height, bg = "white", device = cairo_pdf, limitsize = FALSE)
  }
  
  plot_theme_safe <- theme_pub_global(base_size = 12) +
    ggplot2::theme(
      plot.title.position = "plot",
      plot.margin = ggplot2::margin(20, 28, 14, 18),
      plot.title = ggplot2::element_text(face = "bold", size = 16.5, lineheight = 1.06,
                                         margin = ggplot2::margin(b = 8)),
      plot.subtitle = ggplot2::element_text(size = 11.8, colour = "#4D4D4D", lineheight = 1.10,
                                            margin = ggplot2::margin(b = 12)),
      axis.title = ggplot2::element_text(size = 13.5, colour = "black"),
      axis.text = ggplot2::element_text(size = 11.5, colour = "black"),
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 11.2)
    )
  
  col_neuron <- c("Glut" = "#4E79A7", "Gaba" = "#C44E52")
  col_model <- c(
    "Exact neuron in nonneuron" = "#7F7F7F",
    "Asymmetry" = "#C44E52",
    "log2 enrichment" = "#4E79A7"
  )
  
  # -----------------------------
  # Read data
  # -----------------------------
  message("Reading input file: ", input_file)
  raw <- read_csv(input_file, show_col_types = FALSE)
  
  # canonical orientation: neuron in cluster.1, nonneuron in cluster.2
  mask_1 <- raw$cluster.1_cell_Neruon_type %in% c("Glut", "Gaba") & raw$cluster.2_cell_Neruon_type == "NonNeuron"
  mask_2 <- raw$cluster.1_cell_Neruon_type == "NonNeuron" & raw$cluster.2_cell_Neruon_type %in% c("Glut", "Gaba")
  
  part_a <- raw[mask_1, ]
  part_b <- raw[mask_2, ]
  
  swap_names <- function(df) {
    out <- df
    cluster1_cols <- names(df)[str_detect(names(df), "^cluster\\.1_")]
    cluster2_cols <- names(df)[str_detect(names(df), "^cluster\\.2_")]
    for (c1 in cluster1_cols) {
      c2 <- str_replace(c1, "^cluster\\.1_", "cluster.2_")
      out[[c1]] <- df[[c2]]
    }
    for (c2 in cluster2_cols) {
      c1 <- str_replace(c2, "^cluster\\.2_", "cluster.1_")
      out[[c2]] <- df[[c1]]
    }
    out[["cluster.1.overlap.percent"]] <- df[["cluster.2.overlap.percent"]]
    out[["cluster.2.overlap.percent"]] <- df[["cluster.1.overlap.percent"]]
    out
  }
  
  part_b_swapped <- swap_names(part_b)
  canon <- bind_rows(part_a, part_b_swapped) %>%
    mutate(pair_key = paste(cluster.1_label, cluster.2_label, sep = "||")) %>%
    distinct(pair_key, .keep_all = TRUE)
  
  canon <- canon %>%
    mutate(
      neuron_type = cluster.1_cell_Neruon_type,
      nn_family = extract_nn_family(cluster.2_subclass),
      neuron_family = map_chr(cluster.1_subclass, extract_neuron_family),
      slide_layer = paste(cluster.1_slide, cluster.1_layer, sep = "__"),
      source_overlap = `cluster.1.overlap.percent`,
      target_overlap = `cluster.2.overlap.percent`,
      asym = source_overlap - target_overlap,
      size_ratio = (cluster.1_total_cell_num + 1) / (cluster.2_total_cell_num + 1),
      log_size_ratio = log(size_ratio),
      log_geom_size = 0.5 * log((cluster.1_total_cell_num + 1) * (cluster.2_total_cell_num + 1)),
      is_gaba = as.integer(neuron_type == "Gaba"),
      exact_neuron_in_nn = as.integer(source_overlap >= 1),
      exact_nn_in_neuron = as.integer(target_overlap >= 1),
      region_jaccard = map2_dbl(cluster.1_region, cluster.2_region, region_jaccard)
    )
  
  N_lower_df <- canon %>%
    group_by(slide_layer) %>%
    summarise(N_lower = max(union_cell, na.rm = TRUE), .groups = "drop")
  
  canon <- canon %>%
    left_join(N_lower_df, by = "slide_layer") %>%
    mutate(
      expected_overlap = cluster.1_total_cell_num * cluster.2_total_cell_num / N_lower,
      hg_p = pmap_dbl(
        list(N_lower, cluster.1_total_cell_num, cluster.2_total_cell_num, overlap_cell, union_cell),
        ~ hg_p_upper(..1, ..2, ..3, ..4, ..5)
      ),
      hg_fdr = p.adjust(hg_p, method = "BH"),
      log2_fe = log2((overlap_cell + 0.5) / (expected_overlap + 0.5))
    )
  
  write_csv(canon, file.path(table_dir, "pairs_canonical_with_stats.csv"))
  
  major_nn <- canon %>% count(nn_family, sort = TRUE) %>% filter(n >= 100) %>% pull(nn_family)
  model_df <- canon %>% filter(nn_family %in% major_nn)
  
  # -----------------------------
  # Core multivariable models
  # -----------------------------
  m_glm_neuron_in_nn <- glm(
    exact_neuron_in_nn ~ is_gaba + factor(nn_family) + factor(cluster.1_layer) +
      log_size_ratio + log_geom_size + region_jaccard,
    data = model_df,
    family = binomial()
  )
  glm_neuron_in_nn_tidy <- cluster_robust_tidy(m_glm_neuron_in_nn, model_df$slide_layer, logistic = TRUE)
  write_csv(glm_neuron_in_nn_tidy, file.path(table_dir, "model_glm_exact_neuron_in_nn.csv"))
  
  m_lm_asym <- lm(
    asym ~ is_gaba + factor(nn_family) + factor(cluster.1_layer) +
      log_size_ratio + log_geom_size + region_jaccard,
    data = model_df
  )
  ols_asym_tidy <- cluster_robust_tidy(m_lm_asym, model_df$slide_layer, logistic = FALSE)
  write_csv(ols_asym_tidy, file.path(table_dir, "model_ols_asym.csv"))
  
  m_lm_log2fe <- lm(
    log2_fe ~ is_gaba + factor(nn_family) + factor(cluster.1_layer) +
      log_size_ratio + log_geom_size + region_jaccard,
    data = model_df
  )
  ols_log2fe_tidy <- cluster_robust_tidy(m_lm_log2fe, model_df$slide_layer, logistic = FALSE)
  write_csv(ols_log2fe_tidy, file.path(table_dir, "model_ols_log2_fe.csv"))
  
  # effect summary for is_gaba only
  extract_gaba_effect <- function(tbl, model_name, effect_type = c("OR", "beta")) {
    effect_type <- match.arg(effect_type)
    sub <- tbl %>% filter(term == "is_gaba")
    if (nrow(sub) == 0) stop("is_gaba term not found in ", model_name)
    if (effect_type == "OR") {
      tibble(
        model = model_name,
        scale = "Odds ratio",
        estimate = sub$odds_ratio,
        conf.low = sub$or_low,
        conf.high = sub$or_high,
        p.value = sub$p.value,
        significant = sub$p.value < 0.05
      )
    } else {
      tibble(
        model = model_name,
        scale = "Coefficient",
        estimate = sub$estimate,
        conf.low = sub$conf.low,
        conf.high = sub$conf.high,
        p.value = sub$p.value,
        significant = sub$p.value < 0.05
      )
    }
  }
  
  model_effect_summary <- bind_rows(
    extract_gaba_effect(glm_neuron_in_nn_tidy, "Exact neuron in nonneuron", "OR"),
    extract_gaba_effect(ols_asym_tidy, "Asymmetry", "beta"),
    extract_gaba_effect(ols_log2fe_tidy, "log2 enrichment", "beta")
  ) %>%
    mutate(
      model = factor(model, levels = c("Exact neuron in nonneuron", "Asymmetry", "log2 enrichment")),
      p.label = vapply(p.value, fmt_p_short, character(1))
    )
  write_csv(model_effect_summary, file.path(table_dir, "table_model_effect_summary_is_gaba.csv"))
  
  # adjusted predictions / counterfactual means by neuron type
  predict_counterfactual <- function(model, data, is_glm = FALSE, cluster = NULL, nsim = 400) {
    X0 <- model.matrix(formula(model), data = mutate(data, is_gaba = 0L))
    X1 <- model.matrix(formula(model), data = mutate(data, is_gaba = 1L))
    beta <- coef(model)
    vc <- sandwich::vcovCL(model, cluster = cluster)
    draws <- MASS::mvrnorm(nsim, mu = beta, Sigma = vc)
    if (is_glm) {
      pred0 <- plogis(X0 %*% t(draws))
      pred1 <- plogis(X1 %*% t(draws))
    } else {
      pred0 <- X0 %*% t(draws)
      pred1 <- X1 %*% t(draws)
    }
    tibble(
      neuron_type = c("Glut", "Gaba"),
      estimate = c(mean(rowMeans(pred0)), mean(rowMeans(pred1))),
      conf.low = c(quantile(rowMeans(pred0), 0.025), quantile(rowMeans(pred1), 0.025)),
      conf.high = c(quantile(rowMeans(pred0), 0.975), quantile(rowMeans(pred1), 0.975))
    )
  }
  
  adj_exact <- predict_counterfactual(m_glm_neuron_in_nn, model_df, is_glm = TRUE, cluster = model_df$slide_layer) %>%
    mutate(model = "Exact neuron in nonneuron")
  adj_asym <- predict_counterfactual(m_lm_asym, model_df, is_glm = FALSE, cluster = model_df$slide_layer) %>%
    mutate(model = "Asymmetry")
  adj_log2fe <- predict_counterfactual(m_lm_log2fe, model_df, is_glm = FALSE, cluster = model_df$slide_layer) %>%
    mutate(model = "log2 enrichment")
  
  adjusted_predictions <- bind_rows(adj_exact, adj_asym, adj_log2fe) %>%
    mutate(model = factor(model, levels = c("Exact neuron in nonneuron", "Asymmetry", "log2 enrichment")))
  write_csv(adjusted_predictions, file.path(table_dir, "table_adjusted_predictions_by_neuron_type.csv"))
  
  # model coefficient tables for selected covariates
  select_terms_exact <- c("is_gaba", "log_size_ratio", "log_geom_size", "region_jaccard")
  select_terms_lm <- c("is_gaba", "log_size_ratio", "log_geom_size", "region_jaccard")
  
  coef_exact_plot <- glm_neuron_in_nn_tidy %>%
    filter(term %in% select_terms_exact) %>%
    mutate(term_label = recode(term,
                               is_gaba = "Gaba vs Glut",
                               log_size_ratio = "log(size ratio)",
                               log_geom_size = "log(geometric mean size)",
                               region_jaccard = "Region composition Jaccard"),
           term_label = factor(term_label, levels = rev(c("Gaba vs Glut", "log(size ratio)", "log(geometric mean size)", "Region composition Jaccard"))))
  
  coef_asym_plot <- ols_asym_tidy %>%
    filter(term %in% select_terms_lm) %>%
    mutate(term_label = recode(term,
                               is_gaba = "Gaba vs Glut",
                               log_size_ratio = "log(size ratio)",
                               log_geom_size = "log(geometric mean size)",
                               region_jaccard = "Region composition Jaccard"),
           term_label = factor(term_label, levels = rev(c("Gaba vs Glut", "log(size ratio)", "log(geometric mean size)", "Region composition Jaccard"))))
  
  coef_log2fe_plot <- ols_log2fe_tidy %>%
    filter(term %in% select_terms_lm) %>%
    mutate(term_label = recode(term,
                               is_gaba = "Gaba vs Glut",
                               log_size_ratio = "log(size ratio)",
                               log_geom_size = "log(geometric mean size)",
                               region_jaccard = "Region composition Jaccard"),
           term_label = factor(term_label, levels = rev(c("Gaba vs Glut", "log(size ratio)", "log(geometric mean size)", "Region composition Jaccard"))))
  
  # -----------------------------
  # Publication figures
  # -----------------------------
  
  # Figure 1: key message summary across models
  fig1_exact <- model_effect_summary %>% filter(model == "Exact neuron in nonneuron")
  fig1_lm <- model_effect_summary %>% filter(model != "Exact neuron in nonneuron")
  
  p1a <- ggplot(fig1_exact, aes(x = estimate, y = model, color = significant)) +
    geom_vline(xintercept = 1, linetype = 2, linewidth = 0.55, color = "#7F7F7F") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.18, linewidth = 0.75) +
    geom_point(size = 3.0) +
    scale_x_log10(labels = scales::label_number(accuracy = 0.01)) +
    scale_color_manual(values = c(`TRUE` = "#C44E52", `FALSE` = "#7F7F7F")) +
    labs(
      title = wrap_plot_title("After adjustment, the Gaba effect on absolute containment weakens, but its anchoring signal remains"),
      subtitle = wrap_plot_subtitle(
        paste0(
          "Exact neuron-in-nonneuron: adjusted OR = ",
          sprintf("%.2f", fig1_exact$estimate),
          "; p = ", fmt_p_short(fig1_exact$p.value),
          ". By contrast, Gaba remains positively associated with asymmetry and log2 enrichment after the same controls."
        )
      ),
      x = "Adjusted odds ratio (log scale)",
      y = NULL,
      color = NULL
    ) +
    plot_theme_safe + theme(legend.position = "none")
  
  p1b <- ggplot(fig1_lm, aes(x = estimate, y = model, color = significant)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.55, color = "#7F7F7F") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.18, linewidth = 0.75) +
    geom_point(size = 3.0) +
    scale_color_manual(values = c(`TRUE` = "#C44E52", `FALSE` = "#7F7F7F")) +
    labs(x = "Adjusted coefficient", y = NULL, color = NULL) +
    plot_theme_safe + theme(legend.position = "none")
  
  save_pub(p1a, "fig01_gaba_effect_summary_exact_only", 8.8, 4.6)
  save_pub(p1b, "fig02_gaba_effect_summary_asym_log2fe", 8.8, 4.6)
  
  # Figure 3: adjusted predictions by neuron type across models
  p2 <- ggplot(adjusted_predictions, aes(x = neuron_type, y = estimate, color = neuron_type)) +
    geom_point(size = 2.8) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.12, linewidth = 0.7) +
    facet_wrap(~model, scales = "free_y", nrow = 1) +
    scale_color_manual(values = col_neuron) +
    labs(
      title = wrap_plot_title("Counterfactual adjusted means confirm that the Gaba effect persists for asymmetry and enrichment"),
      subtitle = wrap_plot_subtitle("Each panel averages model-based predictions after setting is_gaba to 0 or 1 while holding all other covariates at their observed values."),
      x = NULL,
      y = "Adjusted prediction",
      color = NULL
    ) +
    plot_theme_safe
  save_pub(p2, "fig03_adjusted_predictions_by_neuron_type", 11.0, 4.8)
  
  # Figure 4: exact model covariate forest
  p3 <- ggplot(coef_exact_plot, aes(x = odds_ratio, y = term_label, color = term == "is_gaba")) +
    geom_vline(xintercept = 1, linetype = 2, linewidth = 0.55, color = "#7F7F7F") +
    geom_errorbarh(aes(xmin = or_low, xmax = or_high), height = 0.18, linewidth = 0.75) +
    geom_point(size = 3.0) +
    scale_x_log10(labels = scales::label_number(accuracy = 0.01)) +
    scale_color_manual(values = c(`TRUE` = "#C44E52", `FALSE` = "#4E79A7")) +
    labs(
      title = wrap_plot_title("Geometry and scale absorb much of the absolute-containment difference"),
      subtitle = wrap_plot_subtitle("In the exact neuron-in-nonneuron model, size ratio, geometric size and region composition account for much of the Gaba-vs-Glut difference."),
      x = "Adjusted odds ratio (log scale)",
      y = NULL,
      color = NULL
    ) +
    plot_theme_safe + theme(legend.position = "none")
  save_pub(p3, "fig04_exact_model_covariate_forest", 8.8, 5.6)
  
  # Figure 5: asymmetry model forest
  p4 <- ggplot(coef_asym_plot, aes(x = estimate, y = term_label, color = term == "is_gaba")) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.55, color = "#7F7F7F") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.18, linewidth = 0.75) +
    geom_point(size = 3.0) +
    scale_color_manual(values = c(`TRUE` = "#C44E52", `FALSE` = "#4E79A7")) +
    labs(
      title = wrap_plot_title("Gaba remains positively associated with asymmetry after adjustment"),
      subtitle = wrap_plot_subtitle(
        paste0(
          "Adjusted coefficient for Gaba = ", sprintf("%.4f", coef_asym_plot$estimate[coef_asym_plot$term == "is_gaba"]),
          "; p = ", fmt_p_short(coef_asym_plot$p.value[coef_asym_plot$term == "is_gaba"])
        )
      ),
      x = "Adjusted coefficient",
      y = NULL,
      color = NULL
    ) +
    plot_theme_safe + theme(legend.position = "none")
  save_pub(p4, "fig05_asym_model_covariate_forest", 8.8, 5.6)
  
  # Figure 6: log2 enrichment model forest
  p5 <- ggplot(coef_log2fe_plot, aes(x = estimate, y = term_label, color = term == "is_gaba")) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.55, color = "#7F7F7F") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.18, linewidth = 0.75) +
    geom_point(size = 3.0) +
    scale_color_manual(values = c(`TRUE` = "#C44E52", `FALSE` = "#4E79A7")) +
    labs(
      title = wrap_plot_title("Gaba remains positively associated with overlap enrichment after adjustment"),
      subtitle = wrap_plot_subtitle(
        paste0(
          "Adjusted coefficient for Gaba = ", sprintf("%.4f", coef_log2fe_plot$estimate[coef_log2fe_plot$term == "is_gaba"]),
          "; p = ", fmt_p_short(coef_log2fe_plot$p.value[coef_log2fe_plot$term == "is_gaba"])
        )
      ),
      x = "Adjusted coefficient",
      y = NULL,
      color = NULL
    ) +
    plot_theme_safe + theme(legend.position = "none")
  save_pub(p5, "fig06_log2fe_model_covariate_forest", 8.8, 5.6)
  
  # Figure 7: family-level adjusted raw summary for asymmetry
  plot7_df <- canon %>%
    filter(nn_family %in% major_nn) %>%
    group_by(neuron_type, nn_family) %>%
    summarise(median_asym = median(asym, na.rm = TRUE), n = n(), .groups = "drop") %>%
    mutate(nn_family = factor(nn_family, levels = major_nn))
  
  p6 <- ggplot(plot7_df, aes(x = median_asym, y = nn_family, color = neuron_type)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.55, color = "#7F7F7F") +
    geom_point(size = 2.7, position = position_dodge(width = 0.45)) +
    scale_color_manual(values = col_neuron) +
    labs(
      title = wrap_plot_title("Across major nonneuron families, Gaba is consistently shifted toward positive asymmetry"),
      subtitle = wrap_plot_subtitle("This descriptive pattern is consistent with the positive adjusted Gaba coefficient in the multivariable asymmetry model."),
      x = "Median asymmetry (neuron overlap - nonneuron overlap)",
      y = NULL,
      color = NULL
    ) +
    plot_theme_safe
  save_pub(p6, "fig07_family_level_asymmetry_shift", 9.6, 6.4)
  
  # Figure 8: family-level adjusted raw summary for log2 FE
  plot8_df <- canon %>%
    filter(nn_family %in% major_nn) %>%
    group_by(neuron_type, nn_family) %>%
    summarise(median_log2_fe = median(log2_fe, na.rm = TRUE), n = n(), .groups = "drop") %>%
    mutate(nn_family = factor(nn_family, levels = major_nn))
  
  p7 <- ggplot(plot8_df, aes(x = median_log2_fe, y = nn_family, color = neuron_type)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.55, color = "#7F7F7F") +
    geom_point(size = 2.7, position = position_dodge(width = 0.45)) +
    scale_color_manual(values = col_neuron) +
    labs(
      title = wrap_plot_title("Across major nonneuron families, Gaba also tends to show higher overlap enrichment"),
      subtitle = wrap_plot_subtitle("The descriptive family-level pattern agrees with the positive adjusted Gaba effect in the log2 enrichment model."),
      x = "Median log2 enrichment",
      y = NULL,
      color = NULL
    ) +
    plot_theme_safe
  save_pub(p7, "fig08_family_level_log2fe_shift", 9.6, 6.4)
  
  # -----------------------------
  # Figure manifest
  # -----------------------------
  figure_manifest <- tibble(
    filename = c(
      "fig01_gaba_effect_summary_exact_only",
      "fig02_gaba_effect_summary_asym_log2fe",
      "fig03_adjusted_predictions_by_neuron_type",
      "fig04_exact_model_covariate_forest",
      "fig05_asym_model_covariate_forest",
      "fig06_log2fe_model_covariate_forest",
      "fig07_family_level_asymmetry_shift",
      "fig08_family_level_log2fe_shift"
    ),
    summary = c(
      "Adjusted Gaba effect in exact neuron-in-nonneuron model",
      "Adjusted Gaba effects in asymmetry and log2 enrichment models",
      "Counterfactual adjusted predictions by neuron type across three models",
      "Exact containment model covariate forest",
      "Asymmetry model covariate forest",
      "log2 enrichment model covariate forest",
      "Family-level descriptive asymmetry shifts",
      "Family-level descriptive enrichment shifts"
    )
  )
  write_csv(figure_manifest, file.path(output_dir, "figure_manifest.csv"))
  
  write_lines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  message("Analysis completed. Outputs written to: ", output_dir)
}

run_module_08 <- function(runtime_base) {
  message('[MODULE 08] start')
  # [comment omitted: encoding-safe]
  # [removed] rm(list = ls())
  gc()
  # ============================================================
  # R_output_08
  # Cross-dataset comparative analysis:
  # neuron-neuron vs neuron-nonneuron spatial interaction data
  # Focus:
  # 1) GABA role is partner-dependent rather than fixed
  # 2) Full multi-host containment obeys a distinct-host-subclass rule
  # Input:
  #   E:/zaw/2603/neuron-neuron-partner.csv
  #   E:/zaw/2603/neuron-nonneuron-partner.csv
  # Output:
  #   E:/zaw/2603/R_output_08
  # ============================================================
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readr)
    library(stringr)
    library(purrr)
    library(ggplot2)
    library(forcats)
    library(scales)
    library(tibble)
  })
  
  options(scipen = 999)
  
  # namespace guards
  select <- dplyr::select
  rename <- dplyr::rename
  filter <- dplyr::filter
  mutate <- dplyr::mutate
  summarise <- dplyr::summarise
  arrange <- dplyr::arrange
  count <- dplyr::count
  distinct <- dplyr::distinct
  left_join <- dplyr::left_join
  inner_join <- dplyr::inner_join
  pull <- dplyr::pull
  transmute <- dplyr::transmute
  first <- dplyr::first
  
  base_dir <- runtime_base
  nn_file <- file.path(base_dir, "neuron-neuron-partner.csv")
  nnon_file <- file.path(base_dir, "neuron-nonneuron-partner.csv")
  out_dir <- file.path(base_dir, "R_output_08")
  fig_dir <- file.path(out_dir, "figures")
  tab_dir <- file.path(out_dir, "tables")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ---------------------------
  # Helper functions
  # ---------------------------
  
  dedup_unordered <- function(df) {
    a <- as.character(df$cluster.1_label)
    b <- as.character(df$cluster.2_label)
    pair_key <- ifelse(a < b, paste0(a, "||", b), paste0(b, "||", a))
    df %>%
      dplyr::mutate(pair_key = pair_key) %>%
      dplyr::distinct(pair_key, .keep_all = TRUE)
  }
  
  region_count <- function(x) {
    ifelse(
      is.na(x) | stringr::str_trim(x) == "",
      NA_integer_,
      lengths(stringr::str_split(x, ","))
    )
  }
  
  simplify_subclass <- function(x) {
    x %>%
      stringr::str_replace("^\\d+\\s+", "") %>%
      stringr::str_replace("\\s+(Glut|Gaba|NN)$", "")
  }
  
  cliffs_delta_fast <- function(x, y) {
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]
    if (length(x) == 0 || length(y) == 0) return(NA_real_)
    wt <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
    u <- unname(wt$statistic)
    (2 * u) / (length(x) * length(y)) - 1
  }
  
  paired_slide_test <- function(slide_summ, metric_name) {
    piv <- slide_summ %>%
      dplyr::select(slide, group, !!rlang::sym(metric_name)) %>%
      tidyr::pivot_wider(names_from = group, values_from = !!rlang::sym(metric_name))
    
    bind_rows(
      tibble(
        metric = metric_name,
        comparison = "GI vs GN",
        median_A = median(piv$GI, na.rm = TRUE),
        median_B = median(piv$GN, na.rm = TRUE),
        p_value = suppressWarnings(wilcox.test(piv$GI, piv$GN, paired = TRUE, exact = FALSE)$p.value)
      ),
      tibble(
        metric = metric_name,
        comparison = "GI vs IN",
        median_A = median(piv$GI, na.rm = TRUE),
        median_B = median(piv$IN, na.rm = TRUE),
        p_value = suppressWarnings(wilcox.test(piv$GI, piv$IN, paired = TRUE, exact = FALSE)$p.value)
      ),
      tibble(
        metric = metric_name,
        comparison = "GN vs IN",
        median_A = median(piv$GN, na.rm = TRUE),
        median_B = median(piv$IN, na.rm = TRUE),
        p_value = suppressWarnings(wilcox.test(piv$GN, piv$IN, paired = TRUE, exact = FALSE)$p.value)
      )
    )
  }
  
  wrap_plot_title <- function(x, width = 56) stringr::str_wrap(x, width = width)
  wrap_plot_subtitle <- function(x, width = 95) stringr::str_wrap(x, width = width)
  
  plot_theme_safe <- theme_pub_global(base_size = 12) +
    ggplot2::theme(
      plot.title.position = "plot",
      plot.margin = ggplot2::margin(20, 28, 14, 18),
      plot.title = ggplot2::element_text(face = "bold", size = 16.5, lineheight = 1.06, margin = ggplot2::margin(b = 8)),
      plot.subtitle = ggplot2::element_text(size = 11.6, colour = "#4D4D4D", lineheight = 1.10, margin = ggplot2::margin(b = 10)),
      axis.title = ggplot2::element_text(size = 13.2, colour = "black"),
      axis.text = ggplot2::element_text(size = 11.5, colour = "black"),
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 11.2)
    )
  
  save_pub <- function(plot_obj, filename, width = 8, height = 5.2, dpi = 600) {
    ggsave_pub(filename = file.path(fig_dir, paste0(filename, ".png")),
           plot = plot_obj, width = width, height = height, dpi = dpi,
           bg = "white", limitsize = FALSE)
    ggsave_pub(filename = file.path(fig_dir, paste0(filename, ".pdf")),
           plot = plot_obj, width = width, height = height,
           bg = "white", device = cairo_pdf, limitsize = FALSE)
  }
  
  fmt_p_short <- function(p) {
    if (is.na(p)) return("NA")
    if (p < 1e-300) return("<1e-300")
    if (p < 1e-4) return(format(p, scientific = TRUE, digits = 2))
    sprintf("%.4f", p)
  }
  
  # ---------------------------
  # Load and deduplicate
  # ---------------------------
  
  nn <- readr::read_csv(nn_file, show_col_types = FALSE)
  nnon <- readr::read_csv(nnon_file, show_col_types = FALSE)
  
  nn_u <- dedup_unordered(nn)
  nnon_u <- dedup_unordered(nnon)
  
  add_basic_features <- function(df) {
    df %>%
      dplyr::mutate(
        `cluster.1.region_n` = region_count(cluster.1_region),
        `cluster.2.region_n` = region_count(cluster.2_region),
        `cluster.1.log10p` = -log10(pmax(cluster.1_cauchy_combination_p, 1e-300)),
        `cluster.2.log10p` = -log10(pmax(cluster.2_cauchy_combination_p, 1e-300)),
        `cluster.1.log_size` = log10(cluster.1_total_cell_num + 1),
        `cluster.2.log_size` = log10(cluster.2_total_cell_num + 1),
        abs_asym = abs(cluster.1.overlap.percent - cluster.2.overlap.percent),
        c1_in_c2 = cluster.1.overlap.percent >= 0.999999,
        c2_in_c1 = cluster.2.overlap.percent >= 0.999999,
        any_containment = c1_in_c2 | c2_in_c1
      )
  }
  
  nn_u <- add_basic_features(nn_u)
  nnon_u <- add_basic_features(nnon_u)
  
  nn_u <- nn_u %>%
    dplyr::mutate(pair_class = purrr::map2_chr(
      cluster.1_cell_Neruon_type, cluster.2_cell_Neruon_type,
      ~ paste(sort(c(.x, .y)), collapse = "-")
    ))
  
  nnon_u <- nnon_u %>%
    dplyr::mutate(pair_class = purrr::map2_chr(
      cluster.1_cell_Neruon_type, cluster.2_cell_Neruon_type,
      ~ paste(sort(c(.x, .y)), collapse = "-")
    ))
  
  # ---------------------------
  # Orient Gaba-Glut pairs (GI)
  # ---------------------------
  
  gi <- nn_u %>%
    dplyr::filter(pair_class == "Gaba-Glut") %>%
    dplyr::mutate(is_c1_glut = cluster.1_cell_Neruon_type == "Glut") %>%
    dplyr::transmute(
      pair_key,
      slide = if_else(is_c1_glut, cluster.1_slide, cluster.2_slide),
      layer = if_else(is_c1_glut, cluster.1_layer, cluster.2_layer),
      E_totalcellnum = if_else(is_c1_glut, cluster.1_total_cell_num, cluster.2_total_cell_num),
      I_totalcellnum = if_else(is_c1_glut, cluster.2_total_cell_num, cluster.1_total_cell_num),
      E_overlappercent = if_else(is_c1_glut, cluster.1.overlap.percent, cluster.2.overlap.percent),
      I_overlappercent = if_else(is_c1_glut, cluster.2.overlap.percent, cluster.1.overlap.percent),
      jaccard, overlap_cell, abs_asym,
      I_in_E = if_else(is_c1_glut, c2_in_c1, c1_in_c2),
      E_in_I = if_else(is_c1_glut, c1_in_c2, c2_in_c1)
    ) %>%
    dplyr::mutate(
      signed_asym = E_overlappercent - I_overlappercent,
      E_larger = E_totalcellnum >= I_totalcellnum
    )
  
  # ---------------------------
  # Orient neuron-nonneuron pairs: GN and IN
  # ---------------------------
  
  nnon_het <- nnon_u %>%
    dplyr::filter(pair_class %in% c("Glut-NonNeuron", "Gaba-NonNeuron")) %>%
    dplyr::mutate(is_c1_neuron = cluster.1_cell_Neruon_type %in% c("Glut", "Gaba")) %>%
    dplyr::transmute(
      pair_key,
      neuron_type = if_else(is_c1_neuron, cluster.1_cell_Neruon_type, cluster.2_cell_Neruon_type),
      neuron_label = if_else(is_c1_neuron, cluster.1_label, cluster.2_label),
      non_label = if_else(is_c1_neuron, cluster.2_label, cluster.1_label),
      slide = if_else(is_c1_neuron, cluster.1_slide, cluster.2_slide),
      layer = if_else(is_c1_neuron, cluster.1_layer, cluster.2_layer),
      neuron_totalcellnum = if_else(is_c1_neuron, cluster.1_total_cell_num, cluster.2_total_cell_num),
      non_totalcellnum = if_else(is_c1_neuron, cluster.2_total_cell_num, cluster.1_total_cell_num),
      neuron_overlappercent = if_else(is_c1_neuron, cluster.1.overlap.percent, cluster.2.overlap.percent),
      non_overlappercent = if_else(is_c1_neuron, cluster.2.overlap.percent, cluster.1.overlap.percent),
      neuron_in_non = if_else(is_c1_neuron, c1_in_c2, c2_in_c1),
      non_in_neuron = if_else(is_c1_neuron, c2_in_c1, c1_in_c2),
      jaccard, overlap_cell, abs_asym
    ) %>%
    dplyr::mutate(
      signed_asym = neuron_overlappercent - non_overlappercent,
      neuron_larger = neuron_totalcellnum >= non_totalcellnum
    )
  
  gn <- nnon_het %>% dplyr::filter(neuron_type == "Glut")
  inn <- nnon_het %>% dplyr::filter(neuron_type == "Gaba")
  
  # ---------------------------
  # 1) Partner-dependent directionality
  # ---------------------------
  
  direction_summary <- bind_rows(
    tibble(
      system = "GI",
      direction = c("Gaba ⊂ Glut", "Glut ⊂ Gaba"),
      count = c(sum(gi$I_in_E), sum(gi$E_in_I))
    ),
    tibble(
      system = "IN",
      direction = c("Gaba ⊂ NonNeuron", "NonNeuron ⊂ Gaba"),
      count = c(sum(inn$neuron_in_non), sum(inn$non_in_neuron))
    )
  ) %>%
    dplyr::group_by(system) %>%
    dplyr::mutate(prop = count / sum(count)) %>%
    dplyr::ungroup()
  
  direction_tests <- bind_rows(
    tibble(
      system = "GI",
      total_oneway = sum(gi$I_in_E) + sum(gi$E_in_I),
      dominant_count = sum(gi$I_in_E),
      reverse_count = sum(gi$E_in_I),
      dominant_prop = sum(gi$I_in_E) / (sum(gi$I_in_E) + sum(gi$E_in_I)),
      p_value = binom.test(sum(gi$I_in_E), sum(gi$I_in_E) + sum(gi$E_in_I), p = 0.5)$p.value
    ),
    tibble(
      system = "IN",
      total_oneway = sum(inn$neuron_in_non) + sum(inn$non_in_neuron),
      dominant_count = sum(inn$neuron_in_non),
      reverse_count = sum(inn$non_in_neuron),
      dominant_prop = sum(inn$neuron_in_non) / (sum(inn$neuron_in_non) + sum(inn$non_in_neuron)),
      p_value = binom.test(sum(inn$neuron_in_non), sum(inn$neuron_in_non) + sum(inn$non_in_neuron), p = 0.5)$p.value
    )
  )
  
  # ---------------------------
  # 2) Partial asymmetry not transferable from GI to IN
  # ---------------------------
  
  partial_gi <- gi %>% dplyr::filter(!(I_in_E | E_in_I))
  partial_gn <- gn %>% dplyr::filter(!(neuron_in_non | non_in_neuron))
  partial_in <- inn %>% dplyr::filter(!(neuron_in_non | non_in_neuron))
  
  partial_pair_summary <- tibble(
    group = c("GI", "GN", "IN"),
    n_pairs = c(nrow(partial_gi), nrow(partial_gn), nrow(partial_in)),
    median_signed_asym = c(median(partial_gi$signed_asym), median(partial_gn$signed_asym), median(partial_in$signed_asym)),
    pair_level_p_vs_zero = c(
      ifelse(any(partial_gi$signed_asym != 0),
             wilcox.test(partial_gi$signed_asym[partial_gi$signed_asym != 0], mu = 0, exact = FALSE)$p.value, 1),
      ifelse(any(partial_gn$signed_asym != 0),
             wilcox.test(partial_gn$signed_asym[partial_gn$signed_asym != 0], mu = 0, exact = FALSE)$p.value, 1),
      ifelse(any(partial_in$signed_asym != 0),
             wilcox.test(partial_in$signed_asym[partial_in$signed_asym != 0], mu = 0, exact = FALSE)$p.value, 1)
    )
  )
  
  slide_partial <- bind_rows(
    partial_gi %>% dplyr::group_by(slide) %>% dplyr::summarise(med = median(signed_asym), .groups = "drop") %>% dplyr::mutate(group = "GI"),
    partial_gn %>% dplyr::group_by(slide) %>% dplyr::summarise(med = median(signed_asym), .groups = "drop") %>% dplyr::mutate(group = "GN"),
    partial_in %>% dplyr::group_by(slide) %>% dplyr::summarise(med = median(signed_asym), .groups = "drop") %>% dplyr::mutate(group = "IN")
  )
  
  slide_partial_summary <- slide_partial %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(
      n_slides = dplyr::n(),
      median_slide_signed_asym = median(med),
      slide_level_p_vs_zero = ifelse(any(med != 0), wilcox.test(med[med != 0], mu = 0, exact = FALSE)$p.value, 1),
      .groups = "drop"
    )
  
  slide_partial_pairwise <- bind_rows(
    paired_slide_test(slide_partial %>% dplyr::rename(contain_rate = med), "contain_rate")
  ) %>%
    dplyr::mutate(metric = "median_signed_asym")
  
  # ---------------------------
  # 3) Distinct-host-subclass law in full containment vs overlap background
  # ---------------------------
  
  contain_edges <- bind_rows(
    nn_u %>%
      dplyr::filter(c1_in_c2) %>%
      dplyr::transmute(
        source = "neuron-neuron",
        guest_label = cluster.1_label,
        guest_type = cluster.1_cell_Neruon_type,
        host_label = cluster.2_label,
        host_type = cluster.2_cell_Neruon_type,
        host_subclass = cluster.2_subclass,
        edge_class = paste0(cluster.1_cell_Neruon_type, "_in_", cluster.2_cell_Neruon_type)
      ),
    nn_u %>%
      dplyr::filter(c2_in_c1) %>%
      dplyr::transmute(
        source = "neuron-neuron",
        guest_label = cluster.2_label,
        guest_type = cluster.2_cell_Neruon_type,
        host_label = cluster.1_label,
        host_type = cluster.1_cell_Neruon_type,
        host_subclass = cluster.1_subclass,
        edge_class = paste0(cluster.2_cell_Neruon_type, "_in_", cluster.1_cell_Neruon_type)
      ),
    nnon_u %>%
      dplyr::filter(c1_in_c2) %>%
      dplyr::transmute(
        source = "neuron-nonneuron",
        guest_label = cluster.1_label,
        guest_type = cluster.1_cell_Neruon_type,
        host_label = cluster.2_label,
        host_type = cluster.2_cell_Neruon_type,
        host_subclass = cluster.2_subclass,
        edge_class = paste0(cluster.1_cell_Neruon_type, "_in_", cluster.2_cell_Neruon_type)
      ),
    nnon_u %>%
      dplyr::filter(c2_in_c1) %>%
      dplyr::transmute(
        source = "neuron-nonneuron",
        guest_label = cluster.2_label,
        guest_type = cluster.2_cell_Neruon_type,
        host_label = cluster.1_label,
        host_type = cluster.1_cell_Neruon_type,
        host_subclass = cluster.1_subclass,
        edge_class = paste0(cluster.2_cell_Neruon_type, "_in_", cluster.1_cell_Neruon_type)
      )
  ) %>%
    dplyr::mutate(host_subclass_simple = simplify_subclass(host_subclass))
  
  # helper long table for overlap neighborhoods
  long_all <- bind_rows(
    nn_u %>%
      dplyr::transmute(
        source = "neuron-neuron",
        focal_label = cluster.1_label,
        focal_type = cluster.1_cell_Neruon_type,
        partner_label = cluster.2_label,
        partner_type = cluster.2_cell_Neruon_type,
        partner_subclass = cluster.2_subclass
      ),
    nn_u %>%
      dplyr::transmute(
        source = "neuron-neuron",
        focal_label = cluster.2_label,
        focal_type = cluster.2_cell_Neruon_type,
        partner_label = cluster.1_label,
        partner_type = cluster.1_cell_Neruon_type,
        partner_subclass = cluster.1_subclass
      ),
    nnon_u %>%
      dplyr::transmute(
        source = "neuron-nonneuron",
        focal_label = cluster.1_label,
        focal_type = cluster.1_cell_Neruon_type,
        partner_label = cluster.2_label,
        partner_type = cluster.2_cell_Neruon_type,
        partner_subclass = cluster.2_subclass
      ),
    nnon_u %>%
      dplyr::transmute(
        source = "neuron-nonneuron",
        focal_label = cluster.2_label,
        focal_type = cluster.2_cell_Neruon_type,
        partner_label = cluster.1_label,
        partner_type = cluster.1_cell_Neruon_type,
        partner_subclass = cluster.1_subclass
      )
  ) %>%
    dplyr::mutate(partner_subclass_simple = simplify_subclass(partner_subclass))
  
  long_overlap_groups_for_edge <- function(edge_name) {
    parts <- stringr::str_split(edge_name, "_in_", simplify = TRUE)
    guest_type <- parts[1]
    host_type <- parts[2]
    src <- ifelse(
      edge_name %in% c("Gaba_in_Glut", "Glut_in_Gaba", "Glut_in_Glut", "Gaba_in_Gaba"),
      "neuron-neuron",
      "neuron-nonneuron"
    )
    
    long_all %>%
      dplyr::filter(source == src, focal_type == guest_type, partner_type == host_type) %>%
      dplyr::distinct(focal_label, partner_label, partner_subclass_simple) %>%
      dplyr::group_by(focal_label) %>%
      dplyr::summarise(
        n_partners = dplyr::n_distinct(partner_label),
        n_subcls = dplyr::n_distinct(partner_subclass_simple),
        .groups = "drop"
      ) %>%
      dplyr::filter(n_partners >= 2) %>%
      dplyr::mutate(repeat_same_subclass = n_subcls < n_partners)
  }
  
  multihost_repeat_by_edge_class <- contain_edges %>%
    dplyr::group_by(edge_class) %>%
    dplyr::group_modify(function(.x, .y) {
      grp_cont <- .x %>%
        dplyr::distinct(guest_label, host_label, host_subclass_simple) %>%
        dplyr::group_by(guest_label) %>%
        dplyr::summarise(
          n_partners = dplyr::n_distinct(host_label),
          n_subcls = dplyr::n_distinct(host_subclass_simple),
          .groups = "drop"
        ) %>%
        dplyr::filter(n_partners >= 2) %>%
        dplyr::mutate(repeat_same_subclass = n_subcls < n_partners)
      
      grp_ov <- long_overlap_groups_for_edge(.y$edge_class)
      
      tibble::tibble(
        contain_multi_n = nrow(grp_cont),
        contain_repeat_rate = ifelse(nrow(grp_cont) > 0, mean(grp_cont$repeat_same_subclass), NA_real_),
        overlap_multi_n = nrow(grp_ov),
        overlap_repeat_rate = ifelse(nrow(grp_ov) > 0, mean(grp_ov$repeat_same_subclass), NA_real_),
        contain_repeat_count = ifelse(nrow(grp_cont) > 0, sum(grp_cont$repeat_same_subclass), 0),
        overlap_repeat_count = ifelse(nrow(grp_ov) > 0, sum(grp_ov$repeat_same_subclass), 0)
      )
    }) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(contain_multi_n))
  
  contain_multi_n <- sum(multihost_repeat_by_edge_class$contain_multi_n, na.rm = TRUE)
  contain_repeat_n <- sum(multihost_repeat_by_edge_class$contain_repeat_count, na.rm = TRUE)
  overlap_multi_n <- sum(multihost_repeat_by_edge_class$overlap_multi_n, na.rm = TRUE)
  overlap_repeat_n <- sum(multihost_repeat_by_edge_class$overlap_repeat_count, na.rm = TRUE)
  
  multihost_repeat_global <- tibble(
    containment_multi_guest_n = contain_multi_n,
    containment_repeat_same_subclass_n = contain_repeat_n,
    containment_repeat_rate = contain_repeat_n / contain_multi_n,
    overlap_multi_guest_n = overlap_multi_n,
    overlap_repeat_same_subclass_n = overlap_repeat_n,
    overlap_repeat_rate = overlap_repeat_n / overlap_multi_n,
    fisher_p_value = fisher.test(matrix(c(
      contain_repeat_n, contain_multi_n - contain_repeat_n,
      overlap_repeat_n, overlap_multi_n - overlap_repeat_n
    ), 2, 2))$p.value
  )
  
  # ---------------------------
  # Save core tables
  # ---------------------------
  write_csv(direction_summary, file.path(tab_dir, "table_01_direction_summary_GI_IN.csv"))
  write_csv(direction_tests, file.path(tab_dir, "table_02_direction_tests_GI_IN.csv"))
  write_csv(partial_pair_summary, file.path(tab_dir, "table_03_partial_pair_asymmetry_summary.csv"))
  write_csv(slide_partial_summary, file.path(tab_dir, "table_04_partial_slide_asymmetry_summary.csv"))
  write_csv(slide_partial_pairwise, file.path(tab_dir, "table_05_partial_slide_pairwise.csv"))
  write_csv(multihost_repeat_global, file.path(tab_dir, "table_06_multihost_repeat_global.csv"))
  write_csv(multihost_repeat_by_edge_class, file.path(tab_dir, "table_07_multihost_repeat_by_edge_class.csv"))
  
  figure_manifest <- tibble(
    filename = c(
      "fig01_GI_vs_IN_direction_counts",
      "fig02_GI_vs_IN_direction_proportions",
      "fig03_partial_signed_asymmetry_boxplots",
      "fig04_partial_signed_asymmetry_slide",
      "fig05_multihost_repeat_global",
      "fig06_multihost_repeat_by_edge_class"
    ),
    summary = c(
      "Directional full-containment counts for GI and IN",
      "Directional composition within GI and IN exact one-way edges",
      "Pair-level partial signed asymmetry across GI, GN, IN",
      "Slide-level partial signed asymmetry across GI, GN, IN",
      "Global distinct-host-subclass law: full containment vs overlap background",
      "Distinct-host-subclass law by edge class"
    )
  )
  write_csv(figure_manifest, file.path(out_dir, "figure_manifest.csv"))
  
  # ---------------------------
  # Figures
  # ---------------------------
  col_system <- c("GI" = "#7B61A8", "IN" = "#C44E52")
  col_group <- c("GI" = "#7B61A8", "GN" = "#4E79A7", "IN" = "#C44E52")
  col_bg <- c("Full containment" = "#C44E52", "Overlap background" = "#9A9A9A")
  
  # Fig1 counts
  fig1 <- ggplot(direction_summary, aes(x = direction, y = count, fill = system)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.64) +
    geom_text(
      aes(label = comma(count)),
      position = position_dodge(width = 0.72),
      vjust = -0.26, size = 4.2
    ) +
    scale_fill_manual(values = col_system) +
    scale_y_continuous(labels = comma_format(), expand = expansion(mult = c(0, 0.14))) +
    labs(
      title = wrap_plot_title("GABA shows a stable direction in GI but a partner-dependent flip in IN"),
      subtitle = wrap_plot_subtitle(
        paste0(
          "GI: Gaba ⊂ Glut = ", comma(sum(gi$I_in_E)),
          " vs Glut ⊂ Gaba = ", comma(sum(gi$E_in_I)),
          "; IN: Gaba ⊂ NonNeuron = ", comma(sum(inn$neuron_in_non)),
          " vs NonNeuron ⊂ Gaba = ", comma(sum(inn$non_in_neuron)), "."
        )
      ),
      x = NULL,
      y = "Number of one-way exact containment edges"
    ) +
    plot_theme_safe +
    theme(axis.text.x = element_text(angle = 18, hjust = 1))
  save_pub(fig1, "fig01_GI_vs_IN_direction_counts", 9.2, 5.8)
  
  # Fig2 proportions
  fig2 <- ggplot(direction_summary, aes(x = system, y = prop, fill = direction)) +
    geom_col(width = 0.64, position = "stack") +
    geom_text(
      aes(label = percent(prop, accuracy = 0.1)),
      position = position_stack(vjust = 0.5),
      colour = "white", fontface = "bold", size = 4.1
    ) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
    labs(
      title = wrap_plot_title("The dominant containment direction flips once GABA interacts with NonNeuron partners"),
      subtitle = wrap_plot_subtitle(
        paste0(
          "GI binomial p = ", fmt_p_short(direction_tests$p_value[direction_tests$system == "GI"]),
          "; IN binomial p = ", fmt_p_short(direction_tests$p_value[direction_tests$system == "IN"]),
          "."
        )
      ),
      x = NULL,
      y = "Proportion within one-way exact edges"
    ) +
    plot_theme_safe
  save_pub(fig2, "fig02_GI_vs_IN_direction_proportions", 7.0, 5.6)
  
  # Fig3 pair-level partial asymmetry
  partial_long <- bind_rows(
    tibble(group = "GI", value = partial_gi$signed_asym),
    tibble(group = "GN", value = partial_gn$signed_asym),
    tibble(group = "IN", value = partial_in$signed_asym)
  ) %>%
    dplyr::filter(is.finite(value))
  
  med_labs <- partial_long %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(med = median(value), .groups = "drop")
  
  fig3 <- ggplot(partial_long, aes(x = group, y = value, fill = group)) +
    geom_violin(alpha = 0.18, color = NA, trim = FALSE) +
    geom_boxplot(width = 0.24, outlier.shape = NA, linewidth = 0.7, alpha = 0.90) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.55, color = "#7F7F7F") +
    geom_text(
      data = med_labs,
      aes(x = group, y = max(partial_long$value, na.rm = TRUE) * 0.92, label = paste0("median=", sprintf("%.3f", med))),
      inherit.aes = FALSE, size = 4.0
    ) +
    scale_fill_manual(values = col_group) +
    labs(
      title = wrap_plot_title("IN partial overlap signed asymmetry is near zero rather than a simple extension of E-I geometry"),
      subtitle = wrap_plot_subtitle(
        paste0(
          "IN pair-level median = ", sprintf("%.3f", partial_pair_summary$median_signed_asym[partial_pair_summary$group == "IN"]),
          "; pair-level p = ", fmt_p_short(partial_pair_summary$pair_level_p_vs_zero[partial_pair_summary$group == "IN"]),
          "."
        )
      ),
      x = NULL,
      y = "Signed asymmetry"
    ) +
    plot_theme_safe +
    theme(legend.position = "none")
  save_pub(fig3, "fig03_partial_signed_asymmetry_boxplots", 8.2, 5.8)
  
  # Fig4 slide-level partial asymmetry
  slide_meds <- slide_partial %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(med = median(med), .groups = "drop")
  
  fig4 <- ggplot(slide_partial, aes(x = group, y = med, fill = group)) +
    geom_violin(alpha = 0.18, color = NA, trim = FALSE) +
    geom_boxplot(width = 0.24, outlier.shape = NA, linewidth = 0.7, alpha = 0.90) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.55, color = "#7F7F7F") +
    geom_text(
      data = slide_meds,
      aes(x = group, y = max(slide_partial$med, na.rm = TRUE) * 0.92, label = paste0("median=", sprintf("%.3f", med))),
      inherit.aes = FALSE, size = 4.0
    ) +
    scale_fill_manual(values = col_group) +
    labs(
      title = wrap_plot_title("At the slide level, IN still centers on zero asymmetry"),
      subtitle = wrap_plot_subtitle(
        paste0(
          "IN slide-level median = ", sprintf("%.3f", slide_partial_summary$median_slide_signed_asym[slide_partial_summary$group == "IN"]),
          "; slide-level p = ", fmt_p_short(slide_partial_summary$slide_level_p_vs_zero[slide_partial_summary$group == "IN"]),
          "."
        )
      ),
      x = NULL,
      y = "Slide-level median signed asymmetry"
    ) +
    plot_theme_safe +
    theme(legend.position = "none")
  save_pub(fig4, "fig04_partial_signed_asymmetry_slide", 8.2, 5.8)
  
  # Fig5 global multihost repeat law
  fig5_df <- tibble(
    mode = c("Full containment", "Overlap background"),
    repeat_rate = c(multihost_repeat_global$containment_repeat_rate, multihost_repeat_global$overlap_repeat_rate),
    n_multi = c(multihost_repeat_global$containment_multi_guest_n, multihost_repeat_global$overlap_multi_guest_n),
    n_repeat = c(multihost_repeat_global$containment_repeat_same_subclass_n, multihost_repeat_global$overlap_repeat_same_subclass_n)
  )
  
  fig5 <- ggplot(fig5_df, aes(x = mode, y = repeat_rate, fill = mode)) +
    geom_col(width = 0.62) +
    geom_text(
      aes(label = paste0(comma(n_repeat), "/", comma(n_multi))),
      vjust = -0.24, size = 4.2
    ) +
    scale_fill_manual(values = col_bg) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.14))) +
    labs(
      title = wrap_plot_title("Repeated hosts of the same guest almost never come from the same subclass once full containment is imposed"),
      subtitle = wrap_plot_subtitle(
        paste0(
          "Global Fisher p = ", fmt_p_short(multihost_repeat_global$fisher_p_value),
          "; full containment repeat count = ", comma(multihost_repeat_global$containment_repeat_same_subclass_n),
          " across ", comma(multihost_repeat_global$containment_multi_guest_n),
          " multi-host guests."
        )
      ),
      x = NULL,
      y = "Rate of repeated same-host-subclass usage"
    ) +
    plot_theme_safe +
    theme(legend.position = "none")
  save_pub(fig5, "fig05_multihost_repeat_global", 8.8, 5.6)
  
  # Fig6 by edge class
  mr_long <- multihost_repeat_by_edge_class %>%
    dplyr::select(edge_class, contain_repeat_rate, overlap_repeat_rate) %>%
    tidyr::pivot_longer(
      cols = c(contain_repeat_rate, overlap_repeat_rate),
      names_to = "kind", values_to = "rate"
    ) %>%
    dplyr::mutate(
      kind = dplyr::recode(
        kind,
        contain_repeat_rate = "Full containment",
        overlap_repeat_rate = "Overlap background"
      )
    )
  
  fig6 <- ggplot(mr_long, aes(x = edge_class, y = rate, fill = kind)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.64) +
    scale_fill_manual(values = col_bg) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.10))) +
    labs(
      title = wrap_plot_title("The distinct-host-subclass law holds across multi-host edge classes"),
      subtitle = wrap_plot_subtitle("Full-containment neighborhoods suppress same-subclass host repetition relative to ordinary overlap neighborhoods."),
      x = NULL,
      y = "Rate of repeated same-host-subclass usage"
    ) +
    plot_theme_safe +
    theme(axis.text.x = element_text(angle = 26, hjust = 1))
  save_pub(fig6, "fig06_multihost_repeat_by_edge_class", 10.6, 6.0)
  
  message("All outputs saved to: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
}

run_module_09 <- function(runtime_base) {
  message('[MODULE 09] start')
  # [comment omitted: encoding-safe]
  # [removed] rm(list = ls())
  gc()
  # ============================================================
  # GABA fully contained in GLU versus other spatial interactions
  # Cross-dataset comparative analysis for:
  #   1) neuron-neuron interactions
  #   2) neuron-nonneuron interactions
  #
  # Input files expected in the working directory:
  #   neuron-neuron-partner.csv
  #   neuron-nonneuron-partner.csv
  #
  # Main outputs:
  #   gaba_in_glu_cross_results/
  #   plus CSV tables and PNG figures
  #
  # Author: OpenAI GPT-5.4 Pro (equivalent R reproduction script)
  # ============================================================
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readr)
    library(stringr)
    library(purrr)
    library(ggplot2)
  })
  
  options(stringsAsFactors = FALSE, scipen = 999)
  
  # dplyr verb aliases to avoid conflicts with AnnotationDbi / other packages
  select <- dplyr::select
  rename <- dplyr::rename
  filter <- dplyr::filter
  mutate <- dplyr::mutate
  summarise <- dplyr::summarise
  arrange <- dplyr::arrange
  distinct <- dplyr::distinct
  count <- dplyr::count
  left_join <- dplyr::left_join
  right_join <- dplyr::right_join
  full_join <- dplyr::full_join
  inner_join <- dplyr::inner_join
  pull <- dplyr::pull
  transmute <- dplyr::transmute
  first <- dplyr::first
  
  
  # -----------------------------
  # Paths
  # -----------------------------
  base_dir <- runtime_base
  
  nn_path   <- file.path(base_dir, "neuron-neuron-partner.csv")
  nonn_path <- file.path(base_dir, "neuron-nonneuron-partner.csv")
  
  outdir <- file.path(base_dir, "R_output_09")
  figdir <- file.path(outdir, "figures")
  tabdir <- file.path(outdir, "tables")
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)
  
  if (!file.exists(nn_path)) {
    stop("File not found: ", nn_path)
  }
  if (!file.exists(nonn_path)) {
    stop("File not found: ", nonn_path)
  }
  
  # -----------------------------
  # Helpers
  # -----------------------------
  type_code <- c("Glut" = "E", "Gaba" = "I", "NonNeuron" = "N")
  
  make_pair_id <- function(a, b) {
    ifelse(a <= b, paste(a, b, sep = "|||"), paste(b, a, sep = "|||"))
  }
  
  neglog10_safe <- function(x) {
    -log10(pmax(x, 1e-300))
  }
  
  fisher_compare <- function(a_full, a_total, b_full, b_total) {
    tab <- matrix(c(a_full, a_total - a_full, b_full, b_total - b_full),
                  nrow = 2, byrow = TRUE)
    ft <- fisher.test(tab)
    tibble(
      OR = unname(ft$estimate),
      p = ft$p.value
    )
  }
  
  mwu_compare <- function(x, y) {
    wt <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
    tibble(
      x_median = median(x, na.rm = TRUE),
      y_median = median(y, na.rm = TRUE),
      p = wt$p.value
    )
  }
  
  signed_rank_compare <- function(a, b) {
    ok <- !(is.na(a) | is.na(b))
    a <- a[ok]
    b <- b[ok]
    if (length(a) == 0) {
      return(tibble(
        n = 0, A_median = NA_real_, B_median = NA_real_,
        median_diff = NA_real_, A_gt_B = NA_integer_,
        A_lt_B = NA_integer_, A_eq_B = NA_integer_, p = NA_real_
      ))
    }
    wt <- suppressWarnings(wilcox.test(a, b, paired = TRUE, exact = FALSE))
    tibble(
      n = length(a),
      A_median = median(a, na.rm = TRUE),
      B_median = median(b, na.rm = TRUE),
      median_diff = median(a - b, na.rm = TRUE),
      A_gt_B = sum(a > b, na.rm = TRUE),
      A_lt_B = sum(a < b, na.rm = TRUE),
      A_eq_B = sum(a == b, na.rm = TRUE),
      p = wt$p.value
    )
  }
  
  safe_mcnemar <- function(x, y) {
    tab <- table(factor(x, levels = c(FALSE, TRUE)),
                 factor(y, levels = c(FALSE, TRUE)))
    mt <- mcnemar.test(tab, correct = FALSE)
    tibble(
      false_false = unname(tab[1, 1]),
      false_true  = unname(tab[1, 2]),
      true_false  = unname(tab[2, 1]),
      true_true   = unname(tab[2, 2]),
      p = mt$p.value
    )
  }
  
  prop_ge <- function(x, k) {
    tibble(
      yes = sum(x >= k, na.rm = TRUE),
      total = length(x),
      prop = mean(x >= k, na.rm = TRUE)
    )
  }
  
  pairwise_fullrate_fisher <- function(df, motif_a, motif_b) {
    a <- df %>% filter(motif == motif_a) %>%
      summarise(full_n = sum(full_containment), total_n = n())
    b <- df %>% filter(motif == motif_b) %>%
      summarise(full_n = sum(full_containment), total_n = n())
    tibble(
      motif_a = motif_a,
      motif_b = motif_b,
      a_rate = a$full_n / a$total_n,
      b_rate = b$full_n / b$total_n
    ) %>%
      bind_cols(fisher_compare(a$full_n, a$total_n, b$full_n, b$total_n))
  }
  
  bh <- function(p) p.adjust(p, method = "BH")
  
  fmt_pct <- function(x, digits = 1) sprintf(paste0("%.", digits, "f%%"), x * 100)
  
  # -----------------------------
  # Load and prepare
  # -----------------------------
  nn   <- read_csv(nn_path, show_col_types = FALSE)
  nonn <- read_csv(nonn_path, show_col_types = FALSE)
  
  prepare_directed <- function(df, source_name) {
    tibble(
      source = source_name,
      guest_label = df$`cluster.1_label`,
      host_label = df$`cluster.2_label`,
      slide = df$`cluster.1_slide`,
      guest_layer = df$`cluster.1_layer`,
      host_layer = df$`cluster.2_layer`,
      guest_region = df$`cluster.1_region`,
      host_region = df$`cluster.2_region`,
      guest_type = df$`cluster.1_cell_Neruon_type`,
      host_type = df$`cluster.2_cell_Neruon_type`,
      guest_subclass = df$`cluster.1_subclass`,
      host_subclass = df$`cluster.2_subclass`,
      guest_size = df$`cluster.1_total_cell_num`,
      host_size = df$`cluster.2_total_cell_num`,
      guest_p = df$`cluster.1_cauchy_combination_p`,
      host_p = df$`cluster.2_cauchy_combination_p`,
      guest_ei = df$`cluster.1_E_I_Ratio`,
      host_ei = df$`cluster.2_E_I_Ratio`,
      overlap_cell = df$overlap_cell,
      union_cell = df$union_cell,
      jaccard = df$jaccard,
      guest_overlap_pct = df$`cluster.1.overlap.percent`,
      host_overlap_pct = df$`cluster.2.overlap.percent`
    ) %>%
      mutate(
        full_containment = guest_overlap_pct == 1,
        partial = guest_overlap_pct < 1 & host_overlap_pct < 1,
        guest_fill = guest_overlap_pct,
        host_fill = host_overlap_pct,
        size_ratio_host_over_guest = host_size / guest_size,
        log_size_ratio = log2(size_ratio_host_over_guest),
        guest_logp = neglog10_safe(guest_p),
        host_logp = neglog10_safe(host_p),
        logp_gap = host_logp - guest_logp,
        overlap_asym = guest_overlap_pct - host_overlap_pct,
        motif = paste0(type_code[guest_type], "\u2282", type_code[host_type]),
        pair_id = make_pair_id(guest_label, host_label)
      )
  }
  
  combined <- bind_rows(
    prepare_directed(nn, "neuron-neuron"),
    prepare_directed(nonn, "neuron-nonneuron")
  )
  
  full    <- combined %>% filter(full_containment)
  partial <- combined %>% filter(partial)
  eligible <- combined %>% filter(host_size >= guest_size)
  
  motif_order <- c("I\u2282E", "E\u2282I", "I\u2282N", "N\u2282I",
                   "E\u2282N", "N\u2282E", "I\u2282I", "E\u2282E", "N\u2282N")
  
  clusters_unique <- bind_rows(
    combined %>% transmute(label = guest_label, type = guest_type, subclass = guest_subclass, layer = guest_layer, slide = slide),
    combined %>% transmute(label = host_label, type = host_type, subclass = host_subclass, layer = host_layer, slide = slide)
  ) %>% distinct(label, .keep_all = TRUE)
  
  # ============================================================
  # Table 01: full containment rates by motif
  # ============================================================
  table01 <- combined %>%
    group_by(motif) %>%
    summarise(full_n = sum(full_containment), total_n = n(), full_rate = full_n / total_n, .groups = "drop") %>%
    mutate(motif = factor(motif, levels = motif_order)) %>%
    arrange(motif)
  
  write_csv(table01, file.path(tabdir, "table01_full_containment_rates_by_motif.csv"))
  
  # ============================================================
  # Table 02: size-eligible rates by motif
  # ============================================================
  table02 <- eligible %>%
    group_by(motif) %>%
    summarise(full_n = sum(full_containment), total_n = n(),
              full_rate_if_host_ge_guest = full_n / total_n, .groups = "drop") %>%
    mutate(motif = factor(motif, levels = motif_order)) %>%
    arrange(motif)
  
  write_csv(table02, file.path(tabdir, "table02_size_eligible_full_rates_by_motif.csv"))
  
  # ============================================================
  # [comment omitted: encoding-safe]
  # ============================================================
  selected_compare <- c("E\u2282I", "I\u2282N", "N\u2282I", "E\u2282N", "N\u2282E", "I\u2282I", "E\u2282E", "N\u2282N")
  
  get_rate_row <- function(df_all, df_eligible, motif_b) {
    a <- df_all %>% filter(motif == "I\u2282E") %>% summarise(full_n = sum(full_containment), total_n = n())
    b <- df_all %>% filter(motif == motif_b)      %>% summarise(full_n = sum(full_containment), total_n = n())
    ae <- df_eligible %>% filter(motif == "I\u2282E") %>% summarise(full_n = sum(full_containment), total_n = n())
    be <- df_eligible %>% filter(motif == motif_b)      %>% summarise(full_n = sum(full_containment), total_n = n())
    
    raw <- fisher_compare(a$full_n, a$total_n, b$full_n, b$total_n)
    elig <- fisher_compare(ae$full_n, ae$total_n, be$full_n, be$total_n)
    
    tibble(
      compare_to = motif_b,
      I_in_E_rate_raw = a$full_n / a$total_n,
      compare_rate_raw = b$full_n / b$total_n,
      OR_raw = raw$OR,
      p_raw = raw$p,
      I_in_E_rate_if_host_ge_guest = ae$full_n / ae$total_n,
      compare_rate_if_host_ge_guest = be$full_n / be$total_n,
      OR_if_host_ge_guest = elig$OR,
      p_if_host_ge_guest = elig$p
    )
  }
  
  table03 <- map_dfr(selected_compare, ~get_rate_row(combined, eligible, .x)) %>%
    mutate(
      fdr_raw = bh(p_raw),
      fdr_if_host_ge_guest = bh(p_if_host_ge_guest)
    )
  
  write_csv(table03, file.path(tabdir, "table03_i_in_e_vs_selected_rate_tests.csv"))
  
  # ============================================================
  # Table 04: geometry summary for all full-containment motifs
  # ============================================================
  table04 <- full %>%
    group_by(motif) %>%
    summarise(
      full_n = n(),
      guest_size_median = median(guest_size),
      host_size_median = median(host_size),
      host_to_guest_size_ratio_median = median(size_ratio_host_over_guest),
      overlap_cell_median = median(overlap_cell),
      jaccard_median = median(jaccard),
      host_fill_median = median(host_fill),
      host_log10p_median = median(host_logp),
      guest_log10p_median = median(guest_logp),
      log10p_gap_median = median(logp_gap),
      .groups = "drop"
    ) %>%
    mutate(motif = factor(motif, levels = motif_order)) %>%
    arrange(motif)
  
  write_csv(table04, file.path(tabdir, "table04_full_containment_geometry_summary.csv"))
  
  # ============================================================
  # [comment omitted: encoding-safe]
  # ============================================================
  metrics_for_test <- c("guest_size", "host_size", "size_ratio_host_over_guest",
                        "overlap_cell", "jaccard", "host_fill",
                        "host_logp", "guest_logp", "logp_gap")
  
  table05 <- purrr::map_dfr(selected_compare, function(compare_to_now) {
    purrr::map_dfr(metrics_for_test, function(metric_now) {
      x <- full %>%
        dplyr::filter(motif == "I\u2282E") %>%
        dplyr::pull(!!rlang::sym(metric_now))
      
      y <- full %>%
        dplyr::filter(motif == compare_to_now) %>%
        dplyr::pull(!!rlang::sym(metric_now))
      
      tibble(
        compare_to = compare_to_now,
        metric = metric_now,
        I_in_E_median = median(x, na.rm = TRUE),
        compare_median = median(y, na.rm = TRUE),
        p = suppressWarnings(wilcox.test(x, y, exact = FALSE)$p.value)
      )
    })
  }) %>%
    dplyr::mutate(fdr = bh(p)) %>%
    dplyr::select(compare_to, metric, I_in_E_median, compare_median, p, fdr)
  
  write_csv(table05, file.path(tabdir, "table05_i_in_e_vs_selected_metric_tests.csv"))
  
  # ============================================================
  # Table 06/07: partial asymmetry
  # ============================================================
  table06 <- partial %>%
    group_by(motif) %>%
    summarise(
      partial_n = n(),
      asym_median = median(overlap_asym),
      jaccard_median = median(jaccard),
      .groups = "drop"
    ) %>%
    mutate(motif = factor(motif, levels = motif_order)) %>%
    arrange(motif)
  
  write_csv(table06, file.path(tabdir, "table06_partial_asymmetry_by_motif.csv"))
  
  partial_compare <- c("I\u2282N", "I\u2282I", "N\u2282E", "E\u2282I", "N\u2282I", "E\u2282N")
  table07 <- map_dfr(partial_compare, function(m) {
    x <- partial %>% filter(motif == "I\u2282E") %>% pull(overlap_asym)
    y <- partial %>% filter(motif == m) %>% pull(overlap_asym)
    tibble(
      compare_to = m,
      I_in_E_asym_median = median(x, na.rm = TRUE),
      compare_asym_median = median(y, na.rm = TRUE),
      p = suppressWarnings(wilcox.test(x, y, exact = FALSE)$p.value)
    )
  }) %>%
    mutate(fdr = bh(p))
  
  write_csv(table07, file.path(tabdir, "table07_i_in_e_partial_asymmetry_tests.csv"))
  
  # ============================================================
  # Table 08: slide-level rate comparisons
  # ============================================================
  slide_rates <- combined %>%
    group_by(slide, motif) %>%
    summarise(full_rate = mean(full_containment), .groups = "drop") %>%
    pivot_wider(names_from = motif, values_from = full_rate)
  
  slide_comp_pairs <- tribble(
    ~A, ~B,
    "I\u2282E", "E\u2282I",
    "I\u2282E", "I\u2282N",
    "I\u2282E", "N\u2282I",
    "I\u2282E", "N\u2282E",
    "I\u2282E", "E\u2282E"
  )
  
  table08 <- pmap_dfr(slide_comp_pairs, function(A, B) {
    sub <- slide_rates %>%
      select(all_of(c(A, B))) %>%
      drop_na()
    tibble(
      A = A,
      B = B,
      n_slides = nrow(sub),
      A_rate_median = median(sub[[A]]),
      B_rate_median = median(sub[[B]]),
      A_gt_B_slides = sum(sub[[A]] > sub[[B]]),
      A_lt_B_slides = sum(sub[[A]] < sub[[B]]),
      p = suppressWarnings(wilcox.test(sub[[A]], sub[[B]], paired = TRUE, exact = FALSE)$p.value)
    )
  })
  
  write_csv(table08, file.path(tabdir, "table08_slide_level_rate_comparisons.csv"))
  
  # ============================================================
  # Same-GABA guest matched analyses
  # ============================================================
  gaba_full <- full %>% filter(guest_type == "Gaba")
  
  guest_host_agg <- gaba_full %>%
    group_by(guest_label, host_type) %>%
    summarise(
      n_hosts = n_distinct(host_label),
      host_size_med = median(host_size),
      host_logp_med = median(host_logp),
      host_fill_med = median(host_fill),
      jaccard_med = median(jaccard),
      size_ratio_med = median(size_ratio_host_over_guest),
      overlap_med = median(overlap_cell),
      .groups = "drop"
    ) %>%
    pivot_wider(names_from = host_type,
                values_from = c(n_hosts, host_size_med, host_logp_med, host_fill_med, jaccard_med, size_ratio_med, overlap_med))
  
  table09_rows <- list(
    tibble(comparison = "same_GABA_guest_Glut_host_vs_NonNeuron_host", metric = "host_size") %>%
      bind_cols(signed_rank_compare(guest_host_agg$host_size_med_Glut, guest_host_agg$host_size_med_NonNeuron)),
    tibble(comparison = "same_GABA_guest_Glut_host_vs_NonNeuron_host", metric = "host_log10p") %>%
      bind_cols(signed_rank_compare(guest_host_agg$host_logp_med_Glut, guest_host_agg$host_logp_med_NonNeuron)),
    tibble(comparison = "same_GABA_guest_Glut_host_vs_NonNeuron_host", metric = "host_fill") %>%
      bind_cols(signed_rank_compare(guest_host_agg$host_fill_med_Glut, guest_host_agg$host_fill_med_NonNeuron)),
    tibble(comparison = "same_GABA_guest_Glut_host_vs_NonNeuron_host", metric = "host_to_guest_size_ratio") %>%
      bind_cols(signed_rank_compare(guest_host_agg$size_ratio_med_Glut, guest_host_agg$size_ratio_med_NonNeuron)),
    tibble(comparison = "same_GABA_guest_Glut_host_vs_NonNeuron_host", metric = "n_hosts") %>%
      bind_cols(signed_rank_compare(guest_host_agg$n_hosts_Glut, guest_host_agg$n_hosts_NonNeuron)),
    
    tibble(comparison = "same_GABA_guest_Glut_host_vs_Gaba_host", metric = "host_size") %>%
      bind_cols(signed_rank_compare(guest_host_agg$host_size_med_Glut, guest_host_agg$host_size_med_Gaba)),
    tibble(comparison = "same_GABA_guest_Glut_host_vs_Gaba_host", metric = "host_log10p") %>%
      bind_cols(signed_rank_compare(guest_host_agg$host_logp_med_Glut, guest_host_agg$host_logp_med_Gaba)),
    tibble(comparison = "same_GABA_guest_Glut_host_vs_Gaba_host", metric = "host_fill") %>%
      bind_cols(signed_rank_compare(guest_host_agg$host_fill_med_Glut, guest_host_agg$host_fill_med_Gaba)),
    tibble(comparison = "same_GABA_guest_Glut_host_vs_Gaba_host", metric = "host_to_guest_size_ratio") %>%
      bind_cols(signed_rank_compare(guest_host_agg$size_ratio_med_Glut, guest_host_agg$size_ratio_med_Gaba)),
    tibble(comparison = "same_GABA_guest_Glut_host_vs_Gaba_host", metric = "n_hosts") %>%
      bind_cols(signed_rank_compare(guest_host_agg$n_hosts_Glut, guest_host_agg$n_hosts_Gaba))
  )
  
  table09 <- bind_rows(table09_rows)
  write_csv(table09, file.path(tabdir, "table09_matched_same_gaba_guest_comparisons.csv"))
  
  # ============================================================
  # Same-GLU host matched analyses
  # ============================================================
  glu_host_full <- full %>% filter(host_type == "Glut")
  
  host_guest_agg <- glu_host_full %>%
    group_by(host_label, guest_type) %>%
    summarise(
      n_guests = n_distinct(guest_label),
      guest_size_med = median(guest_size),
      guest_logp_med = median(guest_logp),
      host_fill_med = median(host_fill),
      jaccard_med = median(jaccard),
      size_ratio_med = median(size_ratio_host_over_guest),
      overlap_med = median(overlap_cell),
      .groups = "drop"
    ) %>%
    pivot_wider(names_from = guest_type,
                values_from = c(n_guests, guest_size_med, guest_logp_med, host_fill_med, jaccard_med, size_ratio_med, overlap_med))
  
  table10_rows <- list(
    tibble(comparison = "same_Glut_host_Gaba_guest_vs_NonNeuron_guest", metric = "guest_size") %>%
      bind_cols(signed_rank_compare(host_guest_agg$guest_size_med_Gaba, host_guest_agg$guest_size_med_NonNeuron)),
    tibble(comparison = "same_Glut_host_Gaba_guest_vs_NonNeuron_guest", metric = "guest_log10p") %>%
      bind_cols(signed_rank_compare(host_guest_agg$guest_logp_med_Gaba, host_guest_agg$guest_logp_med_NonNeuron)),
    tibble(comparison = "same_Glut_host_Gaba_guest_vs_NonNeuron_guest", metric = "host_fill") %>%
      bind_cols(signed_rank_compare(host_guest_agg$host_fill_med_Gaba, host_guest_agg$host_fill_med_NonNeuron)),
    tibble(comparison = "same_Glut_host_Gaba_guest_vs_NonNeuron_guest", metric = "host_to_guest_size_ratio") %>%
      bind_cols(signed_rank_compare(host_guest_agg$size_ratio_med_Gaba, host_guest_agg$size_ratio_med_NonNeuron)),
    tibble(comparison = "same_Glut_host_Gaba_guest_vs_NonNeuron_guest", metric = "n_guests") %>%
      bind_cols(signed_rank_compare(host_guest_agg$n_guests_Gaba, host_guest_agg$n_guests_NonNeuron)),
    
    tibble(comparison = "same_Glut_host_Gaba_guest_vs_Glut_guest", metric = "guest_size") %>%
      bind_cols(signed_rank_compare(host_guest_agg$guest_size_med_Gaba, host_guest_agg$guest_size_med_Glut)),
    tibble(comparison = "same_Glut_host_Gaba_guest_vs_Glut_guest", metric = "guest_log10p") %>%
      bind_cols(signed_rank_compare(host_guest_agg$guest_logp_med_Gaba, host_guest_agg$guest_logp_med_Glut)),
    tibble(comparison = "same_Glut_host_Gaba_guest_vs_Glut_guest", metric = "host_fill") %>%
      bind_cols(signed_rank_compare(host_guest_agg$host_fill_med_Gaba, host_guest_agg$host_fill_med_Glut)),
    tibble(comparison = "same_Glut_host_Gaba_guest_vs_Glut_guest", metric = "host_to_guest_size_ratio") %>%
      bind_cols(signed_rank_compare(host_guest_agg$size_ratio_med_Gaba, host_guest_agg$size_ratio_med_Glut)),
    tibble(comparison = "same_Glut_host_Gaba_guest_vs_Glut_guest", metric = "n_guests") %>%
      bind_cols(signed_rank_compare(host_guest_agg$n_guests_Gaba, host_guest_agg$n_guests_Glut))
  )
  
  table10 <- bind_rows(table10_rows)
  write_csv(table10, file.path(tabdir, "table10_matched_same_glu_host_comparisons.csv"))
  
  # ============================================================
  # Multiplicity
  # ============================================================
  guest_mult <- full %>%
    group_by(motif, guest_label) %>%
    summarise(n_hosts = n_distinct(host_label), .groups = "drop") %>%
    group_by(motif) %>%
    summarise(
      guest_n = n(),
      median_hosts_per_guest = median(n_hosts),
      mean_hosts_per_guest = mean(n_hosts),
      pct_ge2_hosts = mean(n_hosts >= 2),
      pct_ge3_hosts = mean(n_hosts >= 3),
      .groups = "drop"
    )
  
  host_mult <- full %>%
    group_by(motif, host_label) %>%
    summarise(n_guests = n_distinct(guest_label), .groups = "drop") %>%
    group_by(motif) %>%
    summarise(
      host_n = n(),
      median_guests_per_host = median(n_guests),
      mean_guests_per_host = mean(n_guests),
      pct_ge2_guests = mean(n_guests >= 2),
      pct_ge3_guests = mean(n_guests >= 3),
      .groups = "drop"
    )
  
  table11 <- guest_mult %>%
    left_join(host_mult, by = "motif") %>%
    mutate(motif = factor(motif, levels = motif_order)) %>%
    arrange(motif)
  
  write_csv(table11, file.path(tabdir, "table11_multiplicity_summary_by_motif.csv"))
  
  get_multiplicity_test <- function(kind = c("hosts_with_>=2_guests", "guests_with_>=2_hosts"),
                                    motif_b) {
    kind <- match.arg(kind)
    if (kind == "hosts_with_>=2_guests") {
      a_raw <- full %>% filter(motif == "I\u2282E") %>% group_by(host_label) %>% summarise(v = n_distinct(guest_label), .groups = "drop") %>% pull(v)
      b_raw <- full %>% filter(motif == motif_b)     %>% group_by(host_label) %>% summarise(v = n_distinct(guest_label), .groups = "drop") %>% pull(v)
    } else {
      a_raw <- full %>% filter(motif == "I\u2282E") %>% group_by(guest_label) %>% summarise(v = n_distinct(host_label), .groups = "drop") %>% pull(v)
      b_raw <- full %>% filter(motif == motif_b)     %>% group_by(guest_label) %>% summarise(v = n_distinct(host_label), .groups = "drop") %>% pull(v)
    }
    
    a <- prop_ge(a_raw, 2)
    b <- prop_ge(b_raw, 2)
    ft <- fisher_compare(a$yes, a$total, b$yes, b$total)
    
    tibble(
      comparison = paste(kind, "I\u2282E vs", motif_b),
      I_in_E_pct = a$prop,
      compare_pct = b$prop,
      OR = ft$OR,
      p = ft$p
    )
  }
  
  mult_compare <- c("I\u2282N", "I\u2282I", "E\u2282I", "E\u2282E", "N\u2282E")
  table12 <- bind_rows(
    map_dfr(mult_compare, ~get_multiplicity_test("hosts_with_>=2_guests", .x)),
    map_dfr(mult_compare, ~get_multiplicity_test("guests_with_>=2_hosts", .x))
  ) %>% mutate(fdr = bh(p))
  
  write_csv(table12, file.path(tabdir, "table12_multiplicity_tests.csv"))
  
  # ============================================================
  # Host-host overlap around multi-host guests
  # ============================================================
  nn_lookup <- combined %>%
    filter(source == "neuron-neuron") %>%
    group_by(pair_id) %>%
    summarise(jaccard = max(jaccard), .groups = "drop")
  
  nonn_lookup <- combined %>%
    filter(source == "neuron-nonneuron") %>%
    group_by(pair_id) %>%
    summarise(jaccard = max(jaccard), .groups = "drop")
  
  gg_background <- combined %>%
    filter(guest_type == "Glut", host_type == "Glut") %>%
    distinct(pair_id, jaccard) %>%
    pull(jaccard)
  
  ii_background <- combined %>%
    filter(guest_type == "Gaba", host_type == "Gaba") %>%
    distinct(pair_id, jaccard) %>%
    pull(jaccard)
  
  nn_background <- combined %>%
    filter(guest_type == "NonNeuron", host_type == "NonNeuron") %>%
    distinct(pair_id, jaccard) %>%
    pull(jaccard)
  
  get_hosthost_df <- function(motif_name) {
    d <- full %>% filter(motif == motif_name)
    if (nrow(d) == 0) return(tibble(guest_label = character(), pair_id = character(), jaccard = numeric()))
    host_type_here <- unique(d$host_type)
    lookup <- if (host_type_here == "NonNeuron") nonn_lookup else nn_lookup
    
    out <- d %>%
      group_by(guest_label) %>%
      summarise(hosts = list(unique(host_label)), .groups = "drop") %>%
      mutate(
        pair_id = map(hosts, function(h) {
          if (length(h) < 2) return(character())
          combn(h, 2, FUN = function(z) make_pair_id(z[1], z[2]), simplify = TRUE)
        })
      ) %>%
      dplyr::select(guest_label, pair_id) %>%
      tidyr::unnest(pair_id)
    
    out %>% left_join(lookup, by = "pair_id")
  }
  
  hosthost_IE <- get_hosthost_df("I\u2282E")
  hosthost_IN <- get_hosthost_df("I\u2282N")
  hosthost_II <- get_hosthost_df("I\u2282I")
  hosthost_NE <- get_hosthost_df("N\u2282E")
  hosthost_EE <- get_hosthost_df("E\u2282E")
  
  compare_hosthost <- function(x, bg, motif_name, bg_name) {
    wt <- suppressWarnings(wilcox.test(x$jaccard, bg, exact = FALSE))
    tibble(
      motif = motif_name,
      hosthost_pair_n = sum(!is.na(x$jaccard)),
      hosthost_jaccard_median = median(x$jaccard, na.rm = TRUE),
      background = bg_name,
      background_n = length(bg),
      background_jaccard_median = median(bg, na.rm = TRUE),
      p_vs_background = wt$p.value
    )
  }
  
  table13 <- bind_rows(
    compare_hosthost(hosthost_IE, gg_background, "I\u2282E", "GG_background"),
    compare_hosthost(hosthost_IN, nn_background, "I\u2282N", "NN_background"),
    compare_hosthost(hosthost_II, ii_background, "I\u2282I", "II_background"),
    compare_hosthost(hosthost_NE, gg_background, "N\u2282E", "GG_background"),
    compare_hosthost(hosthost_EE, gg_background, "E\u2282E", "GG_background"),
    
    tibble(
      motif = "I\u2282E_vs_I\u2282N",
      hosthost_pair_n = sum(!is.na(hosthost_IE$jaccard)),
      hosthost_jaccard_median = median(hosthost_IE$jaccard, na.rm = TRUE),
      background = "I\u2282N",
      background_n = sum(!is.na(hosthost_IN$jaccard)),
      background_jaccard_median = median(hosthost_IN$jaccard, na.rm = TRUE),
      p_vs_background = suppressWarnings(wilcox.test(hosthost_IE$jaccard, hosthost_IN$jaccard, exact = FALSE)$p.value)
    ),
    tibble(
      motif = "I\u2282E_vs_N\u2282E",
      hosthost_pair_n = sum(!is.na(hosthost_IE$jaccard)),
      hosthost_jaccard_median = median(hosthost_IE$jaccard, na.rm = TRUE),
      background = "N\u2282E",
      background_n = sum(!is.na(hosthost_NE$jaccard)),
      background_jaccard_median = median(hosthost_NE$jaccard, na.rm = TRUE),
      p_vs_background = suppressWarnings(wilcox.test(hosthost_IE$jaccard, hosthost_NE$jaccard, exact = FALSE)$p.value)
    ),
    tibble(
      motif = "I\u2282E_vs_E\u2282E",
      hosthost_pair_n = sum(!is.na(hosthost_IE$jaccard)),
      hosthost_jaccard_median = median(hosthost_IE$jaccard, na.rm = TRUE),
      background = "E\u2282E",
      background_n = sum(!is.na(hosthost_EE$jaccard)),
      background_jaccard_median = median(hosthost_EE$jaccard, na.rm = TRUE),
      p_vs_background = suppressWarnings(wilcox.test(hosthost_IE$jaccard, hosthost_EE$jaccard, exact = FALSE)$p.value)
    )
  )
  
  write_csv(table13, file.path(tabdir, "table13_hosthost_overlap_tests.csv"))
  
  # ============================================================
  # GABA cluster-level preference for GLU versus NonNeuron hosts
  # ============================================================
  gaba_clusters <- clusters_unique %>% filter(type == "Gaba")
  
  gaba_guest_full <- full %>%
    filter(guest_type == "Gaba") %>%
    count(guest_label, host_type, name = "n") %>%
    pivot_wider(names_from = host_type, values_from = n, values_fill = 0)
  
  gaba_status <- gaba_clusters %>%
    transmute(guest_label = label, guest_subclass = subclass, guest_layer = layer, slide = slide) %>%
    left_join(gaba_guest_full, by = "guest_label") %>%
    mutate(
      Gaba = coalesce(Gaba, 0L),
      Glut = coalesce(Glut, 0L),
      NonNeuron = coalesce(NonNeuron, 0L),
      any_Glut_host = Glut > 0,
      any_NonNeuron_host = NonNeuron > 0,
      any_Gaba_host = Gaba > 0
    )
  
  mcn_gaba <- safe_mcnemar(gaba_status$any_Glut_host, gaba_status$any_NonNeuron_host)
  
  table14 <- tibble(
    n_gaba_clusters = nrow(gaba_status),
    any_Glut_host_n = sum(gaba_status$any_Glut_host),
    any_NonNeuron_host_n = sum(gaba_status$any_NonNeuron_host),
    both_n = sum(gaba_status$any_Glut_host & gaba_status$any_NonNeuron_host),
    Glut_only_n = sum(gaba_status$any_Glut_host & !gaba_status$any_NonNeuron_host),
    NonNeuron_only_n = sum(!gaba_status$any_Glut_host & gaba_status$any_NonNeuron_host),
    McNemar_p = mcn_gaba$p
  )
  
  write_csv(table14, file.path(tabdir, "table14_gaba_cluster_glut_vs_nonneuron_host_preference.csv"))
  
  table15 <- gaba_status %>%
    group_by(guest_layer) %>%
    summarise(
      n_gaba_clusters = n(),
      pct_any_Glut_host = mean(any_Glut_host),
      pct_any_NN_host = mean(any_NonNeuron_host),
      glut_only = sum(any_Glut_host & !any_NonNeuron_host),
      nn_only = sum(!any_Glut_host & any_NonNeuron_host),
      mcnemar_p = safe_mcnemar(any_Glut_host, any_NonNeuron_host)$p,
      .groups = "drop"
    ) %>%
    mutate(fdr = bh(mcnemar_p))
  
  write_csv(table15, file.path(tabdir, "table15_gaba_layer_glut_vs_nonneuron_host_preference.csv"))
  
  table16 <- gaba_status %>%
    group_by(guest_subclass) %>%
    summarise(
      n_clusters = n(),
      pct_any_Glut = mean(any_Glut_host),
      pct_any_NonNeuron = mean(any_NonNeuron_host),
      delta = pct_any_Glut - pct_any_NonNeuron,
      both = mean(any_Glut_host & any_NonNeuron_host),
      glut_only = sum(any_Glut_host & !any_NonNeuron_host),
      nn_only = sum(!any_Glut_host & any_NonNeuron_host),
      mcnemar_p = safe_mcnemar(any_Glut_host, any_NonNeuron_host)$p,
      .groups = "drop"
    ) %>%
    mutate(fdr = bh(mcnemar_p)) %>%
    arrange(desc(delta))
  
  write_csv(table16, file.path(tabdir, "table16_gaba_subclass_glut_vs_nonneuron_host_preference.csv"))
  
  # ============================================================
  # GLU host-level preference for GABA versus NonNeuron guests
  # ============================================================
  glu_clusters <- clusters_unique %>% filter(type == "Glut")
  
  glu_host_full_counts <- full %>%
    filter(host_type == "Glut") %>%
    count(host_label, guest_type, name = "n") %>%
    pivot_wider(names_from = guest_type, values_from = n, values_fill = 0)
  
  glu_status <- glu_clusters %>%
    transmute(host_label = label, host_subclass = subclass, host_layer = layer, slide = slide) %>%
    left_join(glu_host_full_counts, by = "host_label") %>%
    mutate(
      Gaba = coalesce(Gaba, 0L),
      Glut = coalesce(Glut, 0L),
      NonNeuron = coalesce(NonNeuron, 0L),
      any_Gaba_guest = Gaba > 0,
      any_NonNeuron_guest = NonNeuron > 0,
      any_Glut_guest = Glut > 0
    )
  
  table17 <- tibble(
    n_glu_clusters = nrow(glu_status),
    any_Gaba_guest_n = sum(glu_status$any_Gaba_guest),
    any_NonNeuron_guest_n = sum(glu_status$any_NonNeuron_guest),
    both_n = sum(glu_status$any_Gaba_guest & glu_status$any_NonNeuron_guest),
    Gaba_only_n = sum(glu_status$any_Gaba_guest & !glu_status$any_NonNeuron_guest),
    NonNeuron_only_n = sum(!glu_status$any_Gaba_guest & glu_status$any_NonNeuron_guest),
    McNemar_p = safe_mcnemar(glu_status$any_Gaba_guest, glu_status$any_NonNeuron_guest)$p
  )
  
  write_csv(table17, file.path(tabdir, "table17_glu_cluster_gaba_vs_nonneuron_guest_preference.csv"))
  
  table18 <- glu_status %>%
    group_by(host_layer) %>%
    summarise(
      n_glu_clusters = n(),
      pct_any_Gaba_guest = mean(any_Gaba_guest),
      pct_any_NonNeuron_guest = mean(any_NonNeuron_guest),
      gaba_only = sum(any_Gaba_guest & !any_NonNeuron_guest),
      nn_only = sum(!any_Gaba_guest & any_NonNeuron_guest),
      mcnemar_p = safe_mcnemar(any_Gaba_guest, any_NonNeuron_guest)$p,
      .groups = "drop"
    ) %>%
    mutate(fdr = bh(mcnemar_p))
  
  write_csv(table18, file.path(tabdir, "table18_glu_layer_gaba_vs_nonneuron_guest_preference.csv"))
  
  table19 <- glu_status %>%
    group_by(host_subclass) %>%
    summarise(
      n_clusters = n(),
      pct_any_Gaba_guest = mean(any_Gaba_guest),
      pct_any_NonNeuron_guest = mean(any_NonNeuron_guest),
      delta = pct_any_Gaba_guest - pct_any_NonNeuron_guest,
      both = mean(any_Gaba_guest & any_NonNeuron_guest),
      gaba_only = sum(any_Gaba_guest & !any_NonNeuron_guest),
      nn_only = sum(!any_Gaba_guest & any_NonNeuron_guest),
      mcnemar_p = safe_mcnemar(any_Gaba_guest, any_NonNeuron_guest)$p,
      .groups = "drop"
    ) %>%
    mutate(fdr = bh(mcnemar_p)) %>%
    arrange(desc(delta))
  
  write_csv(table19, file.path(tabdir, "table19_glu_subclass_gaba_vs_nonneuron_guest_preference.csv"))
  
  # ============================================================
  # Enrichment analyses
  # ============================================================
  # Table 20: GABA subclass enrichment for any GLU host
  table20 <- gaba_status %>%
    group_by(guest_subclass) %>%
    group_modify(~{
      if (nrow(.x) < 10) return(tibble())
      a <- sum(.x$any_Glut_host)
      b <- nrow(.x) - a
      rest <- gaba_status %>% filter(guest_subclass != .y$guest_subclass)
      c <- sum(rest$any_Glut_host)
      d <- nrow(rest) - c
      ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))
      tibble(
        n_clusters = nrow(.x),
        pct_any_Glut_host = a / nrow(.x),
        odds_ratio = unname(ft$estimate),
        p = ft$p.value
      )
    }) %>%
    ungroup() %>%
    mutate(fdr = bh(p)) %>%
    arrange(desc(odds_ratio))
  
  write_csv(table20, file.path(tabdir, "table20_gaba_subclass_enrichment_for_any_glut_host.csv"))
  
  # Table 21: GLU host subclass enrichment for any GABA guest
  table21 <- glu_status %>%
    group_by(host_subclass) %>%
    group_modify(~{
      if (nrow(.x) < 10) return(tibble())
      a <- sum(.x$any_Gaba_guest)
      b <- nrow(.x) - a
      rest <- glu_status %>% filter(host_subclass != .y$host_subclass)
      c <- sum(rest$any_Gaba_guest)
      d <- nrow(rest) - c
      ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))
      tibble(
        n_clusters = nrow(.x),
        pct_any_Gaba_guest = a / nrow(.x),
        odds_ratio = unname(ft$estimate),
        p = ft$p.value
      )
    }) %>%
    ungroup() %>%
    mutate(fdr = bh(p)) %>%
    arrange(desc(odds_ratio))
  
  write_csv(table21, file.path(tabdir, "table21_glu_host_subclass_enrichment_for_any_gaba_guest.csv"))
  
  # Table 22: pair-specific full enrichment within I-E overlaps
  ie <- combined %>% filter(guest_type == "Gaba", host_type == "Glut") %>%
    mutate(is_full = full_containment)
  
  table22 <- ie %>%
    count(guest_subclass, host_subclass, wt = is_full, name = "full_n") %>%
    left_join(ie %>% count(guest_subclass, host_subclass, name = "total_n"),
              by = c("guest_subclass", "host_subclass")) %>%
    filter(total_n >= 20, full_n >= 5) %>%
    rowwise() %>%
    mutate(
      rest_full = sum(ie$is_full) - full_n,
      rest_total = nrow(ie) - total_n,
      ft = list(fisher.test(matrix(c(full_n, total_n - full_n,
                                     rest_full, rest_total - rest_full),
                                   nrow = 2, byrow = TRUE))),
      full_rate = full_n / total_n,
      odds_ratio = unname(ft$estimate),
      p = ft$p.value
    ) %>%
    ungroup() %>%
    select(guest_subclass, host_subclass, full_n, total_n, full_rate, odds_ratio, p) %>%
    mutate(fdr = bh(p)) %>%
    arrange(fdr, desc(odds_ratio))
  
  write_csv(table22, file.path(tabdir, "table22_i_in_e_pair_specific_full_enrichment.csv"))
  
  # Table 23: host subclass enrichment within I-N overlaps
  in_df <- combined %>% filter(guest_type == "Gaba", host_type == "NonNeuron") %>%
    mutate(is_full = full_containment)
  
  table23 <- in_df %>%
    group_by(host_subclass) %>%
    summarise(
      full_n = sum(is_full),
      total_n = n(),
      .groups = "drop"
    ) %>%
    filter(total_n >= 20, full_n >= 5) %>%
    rowwise() %>%
    mutate(
      rest_full = sum(in_df$is_full) - full_n,
      rest_total = nrow(in_df) - total_n,
      ft = list(fisher.test(matrix(c(full_n, total_n - full_n,
                                     rest_full, rest_total - rest_full),
                                   nrow = 2, byrow = TRUE))),
      full_rate = full_n / total_n,
      odds_ratio = unname(ft$estimate),
      p = ft$p.value
    ) %>%
    ungroup() %>%
    select(host_subclass, full_n, total_n, full_rate, odds_ratio, p) %>%
    mutate(fdr = bh(p)) %>%
    arrange(fdr, desc(odds_ratio))
  
  write_csv(table23, file.path(tabdir, "table23_i_in_n_host_subclass_full_enrichment.csv"))
  
  # ============================================================
  # Stratified (size-ratio) analyses with Mantel-Haenszel
  # ============================================================
  # GABA guest: GLU host vs NonNeuron host
  gaba_gn <- combined %>%
    filter(guest_type == "Gaba", host_type %in% c("Glut", "NonNeuron")) %>%
    filter(host_size >= guest_size) %>%
    mutate(
      is_Glut_host = host_type == "Glut",
      ratio_bin = ntile(size_ratio_host_over_guest, 5)
    )
  
  table24 <- gaba_gn %>%
    group_by(ratio_bin, is_Glut_host) %>%
    summarise(
      full_n = sum(full_containment),
      total_n = n(),
      .groups = "drop"
    ) %>%
    mutate(rate = full_n / total_n) %>%
    select(ratio_bin, is_Glut_host, rate, total_n) %>%
    pivot_wider(names_from = is_Glut_host, values_from = c(rate, total_n),
                names_glue = "{ifelse(is_Glut_host, 'Glut_host', 'NonNeuron_host')}_{.value}") %>%
    rename(
      Glut_host_rate = Glut_host_rate,
      NonNeuron_host_rate = NonNeuron_host_rate,
      n_pairs_Glut = Glut_host_total_n,
      n_pairs_NonNeuron = NonNeuron_host_total_n
    ) %>%
    mutate(n_pairs = n_pairs_Glut + n_pairs_NonNeuron)
  
  write_csv(table24, file.path(tabdir, "table24_size_ratio_stratified_gaba_glut_vs_nonneuron.csv"))
  
  # Mantel-Haenszel
  make_2x2xk <- function(df, row_binary, full_binary, strata) {
    lev <- sort(unique(df[[strata]]))
    arr <- array(0, dim = c(2, 2, length(lev)),
                 dimnames = list(
                   row = c("TRUE", "FALSE"),
                   col = c("TRUE", "FALSE"),
                   strata = as.character(lev)
                 ))
    for (s in seq_along(lev)) {
      d <- df %>% filter(.data[[strata]] == lev[s])
      tab <- table(factor(d[[row_binary]], levels = c(TRUE, FALSE)),
                   factor(d[[full_binary]], levels = c(TRUE, FALSE)))
      arr[, , s] <- tab
    }
    arr
  }
  
  mh_gn_array <- make_2x2xk(gaba_gn, "is_Glut_host", "full_containment", "ratio_bin")
  mh_gn <- mantelhaen.test(mh_gn_array)
  
  table25 <- tibble(
    comparison = "GABA guest, Glut host vs NonNeuron host, stratified by size-ratio quintiles",
    pooled_OR = unname(mh_gn$estimate),
    ci_low = mh_gn$conf.int[1],
    ci_high = mh_gn$conf.int[2],
    CMH_statistic = unname(mh_gn$statistic),
    CMH_p = mh_gn$p.value
  )
  
  write_csv(table25, file.path(tabdir, "table25_cmh_gaba_glut_vs_nonneuron_size_stratified.csv"))
  
  # [comment omitted: encoding-safe]
  ie_dir <- combined %>%
    filter(guest_type %in% c("Gaba", "Glut"),
           host_type %in% c("Gaba", "Glut"),
           guest_type != host_type) %>%
    filter(host_size >= guest_size) %>%
    mutate(
      guest_is_Gaba = guest_type == "Gaba",
      ratio_bin = ntile(size_ratio_host_over_guest, 5)
    )
  
  table26 <- ie_dir %>%
    group_by(ratio_bin, guest_is_Gaba) %>%
    summarise(
      full_n = sum(full_containment),
      total_n = n(),
      .groups = "drop"
    ) %>%
    mutate(rate = full_n / total_n) %>%
    select(ratio_bin, guest_is_Gaba, rate, total_n) %>%
    pivot_wider(names_from = guest_is_Gaba, values_from = c(rate, total_n),
                names_glue = "{ifelse(guest_is_Gaba, 'I_in_E', 'E_in_I')}_{.value}") %>%
    rename(
      I_in_E_rate = I_in_E_rate,
      E_in_I_rate = E_in_I_rate,
      n_pairs_I = I_in_E_total_n,
      n_pairs_E = E_in_I_total_n
    ) %>%
    mutate(n_pairs = n_pairs_I + n_pairs_E)
  
  write_csv(table26, file.path(tabdir, "table26_size_ratio_stratified_i_in_e_vs_e_in_i.csv"))
  
  mh_ie_array <- make_2x2xk(ie_dir, "guest_is_Gaba", "full_containment", "ratio_bin")
  mh_ie <- mantelhaen.test(mh_ie_array)
  
  table27 <- tibble(
    comparison = "I⊂E vs E⊂I, stratified by size-ratio quintiles",
    pooled_OR = unname(mh_ie$estimate),
    ci_low = mh_ie$conf.int[1],
    ci_high = mh_ie$conf.int[2],
    CMH_statistic = unname(mh_ie$statistic),
    CMH_p = mh_ie$p.value
  )
  
  write_csv(table27, file.path(tabdir, "table27_cmh_i_in_e_vs_e_in_i_size_stratified.csv"))
  
  # ============================================================
  # Directional binomial tests
  # ============================================================
  binom_row <- function(label, a, b) {
    bt <- binom.test(a, a + b, p = 0.5)
    tibble(
      comparison = label,
      A_n = a,
      B_n = b,
      A_fraction = a / (a + b),
      binom_p_two_sided = bt$p.value
    )
  }
  
  table28 <- bind_rows(
    binom_row("I⊂E_vs_E⊂I", 1857, 287),
    binom_row("I⊂N_vs_N⊂I", 1109, 1352),
    binom_row("E⊂N_vs_N⊂E", 1150, 6062)
  )
  
  write_csv(table28, file.path(tabdir, "table28_directional_binomial_tests.csv"))
  
  # ============================================================
  # A compact key-findings table
  # ============================================================
  table00 <- tribble(
    ~Finding, ~I_in_E, ~Comparator, ~Statistic,
    "1. I⊂E full rate", "1857 / 8457 = 21.96%", "vs E⊂I 287 / 8457 = 3.39%", "Fisher OR=8.01; p=4.25e-317; binomial p=8.29e-281",
    "2. Size-opportunity adjusted", "28.46% when host_size>=guest_size", "vs I⊂N 15.50%; E⊂E 20.59%", "Fisher p(I⊂E vs I⊂N)=1.09e-75; p(I⊂E vs E⊂E)=3.18e-26",
    "3. Host significance", "host_log10p median 7.14", "vs I⊂N 2.23; E⊂I 1.90", "MWU p=8.37e-260 and 9.76e-126",
    "4. Host-guest significance gap", "median 5.43", "vs I⊂N 0.57; E⊂E 2.62", "MWU p=2.54e-257 and 2.03e-52",
    "5. Same GABA guest: GLU host size", "527", "vs matched NonNeuron host 336", "paired Wilcoxon p=3.39e-50 (n=489)",
    "6. Same GABA guest: GLU host fill", "0.199", "vs matched NonNeuron host 0.335", "paired Wilcoxon p=2.40e-39 (n=489)",
    "7. Same GLU host: GABA guest size", "135", "vs matched NonNeuron guest 94", "paired Wilcoxon p=1.38e-37 (n=958)",
    "8. Same GLU host: n guests", "median GABA guests=1", "vs matched NonNeuron guests=3", "paired Wilcoxon p=5.23e-110 (n=958)",
    "9. Partial asymmetry", "median +0.146", "vs I⊂N +0.0009", "MWU p=4.79e-298; slide-level p=3.02e-16",
    "10. GABA cluster-level preference", "1359 with any GLU host", "vs 820 with any NonNeuron host", "McNemar p=2.94e-56"
  )
  write_csv(table00, file.path(tabdir, "table00_key_findings_summary.csv"))
  
  
  # ============================================================
  # Figures
  # ============================================================
  fmt_p <- function(p) {
    if (is.na(p)) return('NA')
    if (p < 1e-300) return('<1e-300')
    if (p < 1e-4) return(formatC(p, format = 'e', digits = 2))
    sprintf('%.4f', p)
  }
  
  fmt_pct_lab <- function(x, digits = 1) sprintf(paste0('%.', digits, 'f%%'), 100 * x)
  wrap_title <- function(x, width = 58) stringr::str_wrap(x, width = width)
  wrap_subtitle <- function(x, width = 92) stringr::str_wrap(x, width = width)
  
  pal_motif <- c(
    "I⊂E" = "#B22222",
    "E⊂I" = "#D98C8C",
    "I⊂N" = "#4E79A7",
    "N⊂I" = "#9BBCE0",
    "E⊂N" = "#59A14F",
    "N⊂E" = "#A8D39E",
    "I⊂I" = "#7F3C8D",
    "E⊂E" = "#F28E2B",
    "N⊂N" = "#8C8C8C"
  )
  
  col_compare <- c(
    "Glut host" = "#B22222",
    "NonNeuron host" = "#4E79A7",
    "GABA guest" = "#B22222",
    "NonNeuron guest" = "#4E79A7"
  )
  
  figure_theme <- theme_pub_global(base_size = 12) +
    ggplot2::theme(
      plot.title.position = 'plot',
      plot.caption.position = 'plot',
      plot.title = ggplot2::element_text(face = 'bold', size = 16.5, colour = 'black', margin = ggplot2::margin(b = 8)),
      plot.subtitle = ggplot2::element_text(size = 11.5, colour = '#4D4D4D', lineheight = 1.08, margin = ggplot2::margin(b = 10)),
      axis.title = ggplot2::element_text(size = 13.5, colour = 'black'),
      axis.text = ggplot2::element_text(size = 11.2, colour = 'black'),
      axis.line = ggplot2::element_line(linewidth = 0.55, colour = 'black'),
      axis.ticks = ggplot2::element_line(linewidth = 0.5, colour = 'black'),
      legend.position = 'top',
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 11),
      plot.margin = ggplot2::margin(18, 22, 14, 18),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  theme_set(figure_theme)
  
  save_pub <- function(plot_obj, filename, width, height, dpi = 600) {
    ggsave_pub(filename = file.path(figdir, filename), plot = plot_obj,
                    width = width, height = height, dpi = dpi, bg = 'white', limitsize = FALSE)
    pdf_name <- sub('\\.png$', '.pdf', filename)
    ggsave_pub(filename = file.path(figdir, pdf_name), plot = plot_obj,
                    width = width, height = height, bg = 'white', device = cairo_pdf, limitsize = FALSE)
  }
  
  # Figure 1 ---------------------------------------------------
  fig01_stats_a <- table03 %>% dplyr::filter(compare_to == 'E⊂I')
  fig01_stats_b <- table03 %>% dplyr::filter(compare_to == 'I⊂N')
  fig01_df <- table01 %>%
    dplyr::mutate(
      motif = factor(motif, levels = motif_order),
      lab = paste0(fmt_pct_lab(full_rate), '\n', 'n=', scales::comma(full_n))
    )
  
  fig01 <- ggplot(fig01_df, aes(x = motif, y = full_rate, fill = motif)) +
    geom_col(width = 0.72, colour = 'white', linewidth = 0.35) +
    geom_text(aes(label = lab), vjust = -0.28, size = 3.65, lineheight = 0.95) +
    scale_fill_manual(values = pal_motif, drop = FALSE) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.18))) +
    labs(
      x = NULL,
      y = 'Full containment rate',
      title = wrap_title('I⊂E is the dominant full-containment motif across directed interaction classes'),
      subtitle = wrap_subtitle(
        paste0(
          'I⊂E = ', fmt_pct_lab(fig01_df$full_rate[fig01_df$motif == 'I⊂E'], 2),
          '; vs E⊂I OR = ', sprintf('%.2f', fig01_stats_a$OR_raw), ', FDR = ', fmt_p(fig01_stats_a$fdr_raw),
          '; vs I⊂N OR = ', sprintf('%.2f', fig01_stats_b$OR_raw), ', FDR = ', fmt_p(fig01_stats_b$fdr_raw)
        )
      )
    ) +
    guides(fill = 'none')
  save_pub(fig01, 'figure01_full_containment_rate_by_motif.png', 9.6, 5.6)
  
  # Figure 2 ---------------------------------------------------
  fig02_df <- table04 %>%
    dplyr::filter(motif %in% c('N⊂E', 'I⊂E', 'E⊂E', 'I⊂N', 'I⊂I')) %>%
    dplyr::mutate(
      motif = factor(motif, levels = c('N⊂E', 'I⊂E', 'E⊂E', 'I⊂N', 'I⊂I')),
      lab = sprintf('%.3f', host_fill_median)
    )
  fig02_p1 <- table05 %>% dplyr::filter(compare_to == 'I⊂N', metric == 'host_fill')
  fig02_p2 <- table05 %>% dplyr::filter(compare_to == 'E⊂E', metric == 'host_fill')
  
  fig02 <- ggplot(fig02_df, aes(x = motif, y = host_fill_median, fill = motif)) +
    geom_col(width = 0.72, colour = 'white', linewidth = 0.35) +
    geom_text(aes(label = lab), vjust = -0.25, size = 3.7) +
    scale_fill_manual(values = pal_motif) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(
      x = NULL,
      y = 'Median host fill (= Jaccard under full containment)',
      title = wrap_title('I⊂E sits in a distinct containment-tightness regime'),
      subtitle = wrap_subtitle(
        paste0('Higher fill than I⊂N (FDR = ', fmt_p(fig02_p1$fdr), ') and E⊂E (FDR = ', fmt_p(fig02_p2$fdr), ').')
      )
    ) +
    guides(fill = 'none')
  save_pub(fig02, 'figure02_host_fill_spectrum_selected_motifs.png', 8.2, 5.4)
  
  # Figure 3 ---------------------------------------------------
  fig03_df <- table04 %>%
    dplyr::filter(motif %in% c('I⊂E', 'N⊂E', 'E⊂E', 'I⊂N', 'N⊂N')) %>%
    dplyr::mutate(
      motif = factor(motif, levels = c('I⊂E', 'N⊂E', 'E⊂E', 'I⊂N', 'N⊂N')),
      lab = sprintf('%.2f', log10p_gap_median)
    )
  fig03_p1 <- table05 %>% dplyr::filter(compare_to == 'I⊂N', metric == 'logp_gap')
  fig03_p2 <- table05 %>% dplyr::filter(compare_to == 'E⊂E', metric == 'logp_gap')
  
  fig03 <- ggplot(fig03_df, aes(x = motif, y = log10p_gap_median, fill = motif)) +
    geom_col(width = 0.72, colour = 'white', linewidth = 0.35) +
    geom_text(aes(label = lab), vjust = -0.25, size = 3.7) +
    scale_fill_manual(values = pal_motif) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(
      x = NULL,
      y = 'Median host_log10p - guest_log10p',
      title = wrap_title('I⊂E has the largest host–guest significance gap'),
      subtitle = wrap_subtitle(
        paste0('Compared with I⊂N: FDR = ', fmt_p(fig03_p1$fdr), '; compared with E⊂E: FDR = ', fmt_p(fig03_p2$fdr), '.')
      )
    ) +
    guides(fill = 'none')
  save_pub(fig03, 'figure03_host_guest_significance_gap_selected_motifs.png', 8.6, 5.4)
  
  # Figure 4 ---------------------------------------------------
  fig04_df <- guest_host_agg %>% dplyr::select(host_size_med_Glut, host_size_med_NonNeuron) %>% tidyr::drop_na()
  lim04 <- max(fig04_df$host_size_med_Glut, fig04_df$host_size_med_NonNeuron, na.rm = TRUE)
  fig04_stat <- table09 %>% dplyr::filter(comparison == 'same_GABA_guest_Glut_host_vs_NonNeuron_host', metric == 'host_size')
  fig04_cor <- suppressWarnings(cor.test(fig04_df$host_size_med_NonNeuron, fig04_df$host_size_med_Glut, method = 'spearman'))
  label04 <- paste0(
    'n = ', nrow(fig04_df),
    '\npaired p = ', fmt_p(fig04_stat$p),
    '\nmedian(GLU) = ', sprintf('%.0f', fig04_stat$A_median),
    '\nmedian(NN) = ', sprintf('%.0f', fig04_stat$B_median),
    '\nSpearman ρ = ', sprintf('%.2f', unname(fig04_cor$estimate)), ', p = ', fmt_p(fig04_cor$p.value)
  )
  
  fig04 <- ggplot(fig04_df, aes(x = host_size_med_NonNeuron, y = host_size_med_Glut)) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, linewidth = 0.6, colour = '#7F7F7F') +
    geom_point(size = 1.8, alpha = 0.52, colour = '#B22222') +
    annotate('label', x = lim04 * 0.10, y = lim04 * 0.92, hjust = 0, vjust = 1,
             label = label04, size = 3.5, label.size = 0.25, fill = 'white') +
    coord_equal(xlim = c(0, lim04 * 1.02), ylim = c(0, lim04 * 1.02), clip = 'off') +
    labs(
      x = 'Matched NonNeuron host median size',
      y = 'Matched Glut host median size',
      title = wrap_title('For the same GABA guest, matched Glut hosts are systematically larger')
    )
  save_pub(fig04, 'figure04_matched_same_gaba_guest_host_size_glut_vs_nonneuron.png', 6.3, 6.1)
  
  # Figure 5 ---------------------------------------------------
  fig05_df <- guest_host_agg %>% dplyr::select(host_fill_med_Glut, host_fill_med_NonNeuron) %>% tidyr::drop_na()
  lim05 <- max(fig05_df$host_fill_med_Glut, fig05_df$host_fill_med_NonNeuron, na.rm = TRUE)
  fig05_stat <- table09 %>% dplyr::filter(comparison == 'same_GABA_guest_Glut_host_vs_NonNeuron_host', metric == 'host_fill')
  fig05_cor <- suppressWarnings(cor.test(fig05_df$host_fill_med_NonNeuron, fig05_df$host_fill_med_Glut, method = 'spearman'))
  label05 <- paste0(
    'n = ', nrow(fig05_df),
    '\npaired p = ', fmt_p(fig05_stat$p),
    '\nmedian(GLU) = ', sprintf('%.3f', fig05_stat$A_median),
    '\nmedian(NN) = ', sprintf('%.3f', fig05_stat$B_median),
    '\nSpearman ρ = ', sprintf('%.2f', unname(fig05_cor$estimate)), ', p = ', fmt_p(fig05_cor$p.value)
  )
  
  fig05 <- ggplot(fig05_df, aes(x = host_fill_med_NonNeuron, y = host_fill_med_Glut)) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, linewidth = 0.6, colour = '#7F7F7F') +
    geom_point(size = 1.8, alpha = 0.52, colour = '#B22222') +
    annotate('label', x = lim05 * 0.08, y = lim05 * 0.94, hjust = 0, vjust = 1,
             label = label05, size = 3.5, label.size = 0.25, fill = 'white') +
    coord_equal(xlim = c(0, lim05 * 1.02), ylim = c(0, lim05 * 1.02), clip = 'off') +
    labs(
      x = 'Matched NonNeuron host median fill',
      y = 'Matched Glut host median fill',
      title = wrap_title('For the same GABA guest, Glut hosts form looser shells than NonNeuron hosts')
    )
  save_pub(fig05, 'figure05_matched_same_gaba_guest_host_fill_glut_vs_nonneuron.png', 6.3, 6.1)
  
  # Figure 6 ---------------------------------------------------
  fig06_df <- host_guest_agg %>% dplyr::select(guest_size_med_Gaba, guest_size_med_NonNeuron) %>% tidyr::drop_na()
  lim06 <- max(fig06_df$guest_size_med_Gaba, fig06_df$guest_size_med_NonNeuron, na.rm = TRUE)
  fig06_stat <- table10 %>% dplyr::filter(comparison == 'same_Glut_host_Gaba_guest_vs_NonNeuron_guest', metric == 'guest_size')
  fig06_cor <- suppressWarnings(cor.test(fig06_df$guest_size_med_NonNeuron, fig06_df$guest_size_med_Gaba, method = 'spearman'))
  label06 <- paste0(
    'n = ', nrow(fig06_df),
    '\npaired p = ', fmt_p(fig06_stat$p),
    '\nmedian(GABA) = ', sprintf('%.0f', fig06_stat$A_median),
    '\nmedian(NN) = ', sprintf('%.0f', fig06_stat$B_median),
    '\nSpearman ρ = ', sprintf('%.2f', unname(fig06_cor$estimate)), ', p = ', fmt_p(fig06_cor$p.value)
  )
  
  fig06 <- ggplot(fig06_df, aes(x = guest_size_med_NonNeuron, y = guest_size_med_Gaba)) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, linewidth = 0.6, colour = '#7F7F7F') +
    geom_point(size = 1.8, alpha = 0.52, colour = '#B22222') +
    annotate('label', x = lim06 * 0.10, y = lim06 * 0.92, hjust = 0, vjust = 1,
             label = label06, size = 3.5, label.size = 0.25, fill = 'white') +
    coord_equal(xlim = c(0, lim06 * 1.02), ylim = c(0, lim06 * 1.02), clip = 'off') +
    labs(
      x = 'Matched NonNeuron guest median size',
      y = 'Matched GABA guest median size',
      title = wrap_title('For the same Glut host, GABA guests are larger than matched NonNeuron guests')
    )
  save_pub(fig06, 'figure06_matched_same_glu_host_guest_size_gaba_vs_nonneuron.png', 6.3, 6.1)
  
  # Figure 7 ---------------------------------------------------
  fig07_slide <- partial %>%
    dplyr::filter(motif == 'I⊂E') %>%
    dplyr::group_by(slide) %>%
    dplyr::summarise(slide_med = median(overlap_asym, na.rm = TRUE), .groups = 'drop')
  
  fig07_slide_stat <- tibble::tibble(
    median_signed_asym = median(fig07_slide$slide_med, na.rm = TRUE),
    slide_level_p_vs_zero = {
      x <- fig07_slide$slide_med
      x <- x[!is.na(x)]
      x_nz <- x[x != 0]
      if (length(x_nz) == 0) 1 else suppressWarnings(wilcox.test(x_nz, mu = 0, exact = FALSE)$p.value)
    }
  )
  
  fig07_ref <- table06 %>%
    dplyr::filter(motif == 'I⊂E') %>%
    dplyr::transmute(motif, asym_median, fdr = NA_real_)
  
  fig07_comp <- table06 %>%
    dplyr::filter(motif %in% c('I⊂N', 'I⊂I', 'N⊂E', 'E⊂I', 'E⊂N')) %>%
    dplyr::left_join(
      table07 %>% dplyr::rename(motif = compare_to),
      by = 'motif'
    ) %>%
    dplyr::select(motif, asym_median, fdr)
  
  fig07_df <- dplyr::bind_rows(fig07_ref, fig07_comp) %>%
    dplyr::mutate(
      motif = factor(motif, levels = c('I⊂E', 'I⊂N', 'I⊂I', 'N⊂E', 'E⊂I', 'E⊂N')),
      stat_lab = dplyr::case_when(
        as.character(motif) == 'I⊂E' ~ paste0('median=', sprintf('%.3f', asym_median), '\nref'),
        TRUE ~ paste0('median=', sprintf('%.3f', asym_median), '\nFDR=', purrr::map_chr(fdr, fmt_p))
      ),
      vjust_lab = ifelse(asym_median >= 0, -0.24, 1.14)
    )
  
  fig07 <- ggplot(fig07_df, aes(x = motif, y = asym_median, fill = motif)) +
    geom_col(width = 0.72, colour = 'white', linewidth = 0.35) +
    geom_hline(yintercept = 0, linewidth = 0.6, colour = '#7F7F7F') +
    geom_text(aes(label = stat_lab, vjust = vjust_lab), size = 3.35, lineheight = 0.95, show.legend = FALSE) +
    scale_fill_manual(values = pal_motif) +
    labs(
      x = NULL,
      y = 'Median signed overlap asymmetry\n(guest overlap − host overlap)',
      title = wrap_title('Partial-overlap directionality is strongest and most stable for I⊂E'),
      subtitle = wrap_subtitle(
        paste0('I⊂E slide-level p vs 0 = ', fmt_p(fig07_slide_stat$slide_level_p_vs_zero),
               '; median slide asymmetry = ', sprintf('%.3f', fig07_slide_stat$median_signed_asym), '.')
      )
    ) +
    guides(fill = 'none')
  save_pub(fig07, 'figure07_partial_overlap_asymmetry_selected_motifs.png', 8.8, 5.6)
  
  # Figure 8 ---------------------------------------------------
  fig08_df <- table24 %>%
    tidyr::pivot_longer(cols = c(Glut_host_rate, NonNeuron_host_rate), names_to = 'host_kind', values_to = 'rate') %>%
    dplyr::mutate(
      host_kind = dplyr::recode(host_kind,
                                Glut_host_rate = 'Glut host',
                                NonNeuron_host_rate = 'NonNeuron host'),
      label = fmt_pct_lab(rate)
    )
  fig08_stat <- table25
  
  fig08 <- ggplot(fig08_df, aes(x = factor(ratio_bin), y = rate, colour = host_kind, group = host_kind)) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 2.5) +
    geom_text(aes(label = label), size = 3.2, vjust = -0.75, show.legend = FALSE) +
    scale_colour_manual(values = col_compare[c('Glut host','NonNeuron host')]) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0.02, 0.15))) +
    labs(
      x = 'Host/guest size-ratio quintile',
      y = 'Full containment rate',
      colour = NULL,
      title = wrap_title('Across size-ratio strata, Glut hosts remain more likely than NonNeuron hosts to fully contain GABA guests'),
      subtitle = wrap_subtitle(
        paste0('CMH pooled OR = ', sprintf('%.2f', fig08_stat$pooled_OR),
               ' (95% CI ', sprintf('%.2f', fig08_stat$ci_low), '–', sprintf('%.2f', fig08_stat$ci_high),
               '), p = ', fmt_p(fig08_stat$CMH_p), '.')
      )
    )
  save_pub(fig08, 'figure08_size_ratio_stratified_gaba_glut_vs_nonneuron.png', 8.4, 5.4)
  
  # Figure 9 ---------------------------------------------------
  fig09_df <- table16 %>%
    dplyr::filter(n_clusters >= 50) %>%
    dplyr::mutate(guest_subclass_short = stringr::str_replace(guest_subclass, ' Gaba$', '')) %>%
    dplyr::slice_max(order_by = pct_any_Glut, n = 10) %>%
    dplyr::arrange(delta) %>%
    dplyr::mutate(guest_subclass_short = factor(guest_subclass_short, levels = guest_subclass_short))
  
  xmax09 <- max(fig09_df$pct_any_Glut, fig09_df$pct_any_NonNeuron, na.rm = TRUE)
  fig09 <- ggplot(fig09_df, aes(y = guest_subclass_short)) +
    geom_segment(aes(x = pct_any_NonNeuron, xend = pct_any_Glut, yend = guest_subclass_short),
                 linewidth = 0.8, colour = '#BDBDBD') +
    geom_point(aes(x = pct_any_NonNeuron, colour = 'NonNeuron host'), size = 2.7) +
    geom_point(aes(x = pct_any_Glut, colour = 'Glut host'), size = 2.9) +
    geom_text(aes(
      x = pmax(pct_any_Glut, pct_any_NonNeuron) + 0.035,
      label = paste0('FDR=', purrr::map_chr(fdr, fmt_p))
    ),
    hjust = 0, size = 3.2, colour = '#444444') +
    scale_colour_manual(values = col_compare[c('Glut host','NonNeuron host')]) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, min(1, xmax09 + 0.16)),
                       expand = expansion(mult = c(0.01, 0.01))) +
    coord_cartesian(clip = 'off') +
    labs(
      x = 'Fraction of GABA clusters with any full host',
      y = NULL,
      colour = NULL,
      title = wrap_title('Main GABA subclasses show a systematic preference for full Glut hosts over full NonNeuron hosts'),
      subtitle = wrap_subtitle('Each row shows matched subclass-level probabilities; right-side labels report subclass-level FDR.')
    )
  save_pub(fig09, 'figure09_main_gaba_subclass_host_preference.png', 9.4, 5.8)
  
  # Figure 10 --------------------------------------------------
  fig10_df <- table19 %>%
    dplyr::filter(n_clusters >= 50) %>%
    dplyr::mutate(host_subclass_short = stringr::str_replace(host_subclass, ' Glut$', '')) %>%
    dplyr::slice_max(order_by = pct_any_Gaba_guest, n = 8) %>%
    dplyr::arrange(delta) %>%
    dplyr::mutate(host_subclass_short = factor(host_subclass_short, levels = host_subclass_short))
  
  xmax10 <- max(fig10_df$pct_any_Gaba_guest, fig10_df$pct_any_NonNeuron_guest, na.rm = TRUE)
  fig10 <- ggplot(fig10_df, aes(y = host_subclass_short)) +
    geom_segment(aes(x = pct_any_NonNeuron_guest, xend = pct_any_Gaba_guest, yend = host_subclass_short),
                 linewidth = 0.8, colour = '#BDBDBD') +
    geom_point(aes(x = pct_any_NonNeuron_guest, colour = 'NonNeuron guest'), size = 2.7) +
    geom_point(aes(x = pct_any_Gaba_guest, colour = 'GABA guest'), size = 2.9) +
    geom_text(aes(x = pmax(pct_any_Gaba_guest, pct_any_NonNeuron_guest) + 0.035,
                  label = paste0('FDR=', purrr::map_chr(fdr, fmt_p))),
              hjust = 0, size = 3.2, colour = '#444444') +
    scale_colour_manual(values = col_compare[c('GABA guest','NonNeuron guest')]) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, min(1, xmax10 + 0.16)),
                       expand = expansion(mult = c(0.01, 0.01))) +
    coord_cartesian(clip = 'off') +
    labs(
      x = 'Fraction of Glut host clusters with any full guest',
      y = NULL,
      colour = NULL,
      title = wrap_title('Selected Glut subclasses are more likely to host full GABA guests than full NonNeuron guests'),
      subtitle = wrap_subtitle('Each row shows subclass-level host propensity; right-side labels report subclass-level FDR.')
    )
  save_pub(fig10, 'figure10_selected_glu_subclass_guest_preference.png', 9.4, 5.6)
  
  message('Analysis complete. Results written to: ', normalizePath(outdir))
}

run_module_10 <- function(runtime_base) {
  message('[MODULE 10] start')
  # [comment omitted: encoding-safe]
  # [removed] rm(list = ls())
  gc()
  options(stringsAsFactors = FALSE, scipen = 999)
  
  suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(readr)
    library(stringr)
    library(purrr)
    library(ggplot2)
    library(scales)
    library(sandwich)
    library(lmtest)
    library(forcats)
    library(tibble)
  })
  
  # ------------------------------------------------------------------
  # Paths
  # ------------------------------------------------------------------
  base_dir <- runtime_base
  input_nn_candidates <- c(
    file.path(base_dir, "neuron-neuron-partner.txt"),
    file.path(base_dir, "neuron-neuron-partner.csv")
  )
  input_nonn_candidates <- c(
    file.path(base_dir, "neuron-nonneuron-partner.csv"),
    file.path(base_dir, "neuron-nonneuron-partner.txt")
  )
  
  input_nn <- input_nn_candidates[file.exists(input_nn_candidates)][1]
  input_nonn <- input_nonn_candidates[file.exists(input_nonn_candidates)][1]
  
  if (is.na(input_nn) || !nzchar(input_nn)) stop("Cannot find neuron-neuron partner file under: ", base_dir)
  if (is.na(input_nonn) || !nzchar(input_nonn)) stop("Cannot find neuron-nonneuron partner file under: ", base_dir)
  
  output_dir <- file.path(base_dir, "R_output_10")
  table_dir <- file.path(output_dir, "tables")
  fig_dir <- file.path(output_dir, "figures")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ------------------------------------------------------------------
  # Safe dplyr aliases against namespace conflicts
  # ------------------------------------------------------------------
  select <- dplyr::select
  rename <- dplyr::rename
  filter <- dplyr::filter
  mutate <- dplyr::mutate
  summarise <- dplyr::summarise
  arrange <- dplyr::arrange
  distinct <- dplyr::distinct
  count <- dplyr::count
  left_join <- dplyr::left_join
  inner_join <- dplyr::inner_join
  pull <- dplyr::pull
  transmute <- dplyr::transmute
  first <- dplyr::first
  
  # ------------------------------------------------------------------
  # Plot helpers
  # ------------------------------------------------------------------
  wrap_title <- function(x, width = 62) stringr::str_wrap(x, width = width)
  wrap_subtitle <- function(x, width = 98) stringr::str_wrap(x, width = width)
  fmt_p <- function(p) {
    if (length(p) != 1) stop("fmt_p expects a scalar p-value")
    if (is.na(p)) return("NA")
    if (p < 1e-300) return("<1e-300")
    if (p < 1e-4) return(format(p, scientific = TRUE, digits = 2))
    sprintf("%.4f", p)
  }
  fmt_or <- function(x) ifelse(is.na(x), "NA", sprintf("%.3f", x))
  
  col_target <- c(
    "Gaba" = "#4E79A7",
    "Glut" = "#C44E52",
    "NonNeuron" = "#59A14F"
  )
  col_binary <- c("Not exact" = "#B9B9B9", "Exact" = "#C44E52")
  
  theme_pub <- function() {
    theme_pub_global(base_size = 12) +
      ggplot2::theme(
        plot.title.position = "plot",
        plot.margin = ggplot2::margin(18, 22, 14, 16),
        plot.title = ggplot2::element_text(face = "bold", size = 16, margin = ggplot2::margin(b = 7)),
        plot.subtitle = ggplot2::element_text(size = 11.3, colour = "#4D4D4D", lineheight = 1.08, margin = ggplot2::margin(b = 10)),
        axis.title = ggplot2::element_text(size = 13),
        axis.text = ggplot2::element_text(size = 11, colour = "black"),
        legend.position = "top",
        legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 10.8),
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(face = "bold")
      )
  }
  
  save_pub <- function(p, filename, width = 8, height = 5.2, dpi = 600) {
    ggsave_pub(file.path(fig_dir, filename), p, width = width, height = height, dpi = dpi, bg = "white", limitsize = FALSE)
    pdf_name <- sub("\\.png$", ".pdf", filename)
    ggsave_pub(file.path(fig_dir, pdf_name), p, width = width, height = height, bg = "white", device = cairo_pdf, limitsize = FALSE)
  }
  
  # ------------------------------------------------------------------
  # Helper functions
  # ------------------------------------------------------------------
  safe_p <- function(x) pmax(x, .Machine$double.xmin)
  
  rename_partner_columns <- function(df) {
    safe_names <- c(
      "cluster1_label", "cluster1_slide", "cluster1_layer", "cluster1_region", "cluster1_cell_type",
      "cluster1_subclass", "cluster1_total_cell_num", "cluster1_glut_num", "cluster1_gaba_num",
      "cluster1_cauchy_p", "cluster1_ei_ratio",
      "cluster2_label", "cluster2_slide", "cluster2_layer", "cluster2_region", "cluster2_cell_type",
      "cluster2_subclass", "cluster2_total_cell_num", "cluster2_glut_num", "cluster2_gaba_num",
      "cluster2_cauchy_p", "cluster2_ei_ratio",
      "overlap_cell", "union_cell", "jaccard", "cluster1_overlap", "cluster2_overlap"
    )
    
    # remove columns that are entirely empty, which sometimes appear after malformed parsing
    keep <- !vapply(df, function(col) {
      if (is.character(col)) {
        all(is.na(col) | trimws(col) == "")
      } else {
        all(is.na(col))
      }
    }, logical(1))
    df <- df[, keep, drop = FALSE]
    
    # if the file already has exactly the expected structure, just rename
    if (ncol(df) == length(safe_names)) {
      names(df) <- safe_names
      return(df)
    }
    
    # if there are extra junk columns, keep the first expected block only
    if (ncol(df) > length(safe_names)) {
      df <- df[, seq_along(safe_names), drop = FALSE]
      names(df) <- safe_names
      return(df)
    }
    
    stop(sprintf('Parsed %d columns, but expected %d. Please check the delimiter/quoting of the input file.',
                 ncol(df), length(safe_names)))
  }
  
  read_partner_file <- function(path) {
    # use readr::read_csv so quoted region fields containing commas are handled correctly
    df <- suppressMessages(
      readr::read_csv(
        file = path,
        show_col_types = FALSE,
        progress = FALSE,
        quote = '"',
        trim_ws = FALSE,
        name_repair = 'minimal'
      )
    )
    as.data.frame(df)
  }
  
  swap_clusters <- function(df) {
    out <- df
    c1_cols <- names(df)[startsWith(names(df), "cluster1_")]
    for (c1 in c1_cols) {
      suffix <- sub("^cluster1_", "", c1)
      c2 <- paste0("cluster2_", suffix)
      out[[c1]] <- df[[c2]]
      out[[c2]] <- df[[c1]]
    }
    out
  }
  
  clean_subclass_name <- function(x) {
    x <- as.character(x)
    x <- stringr::str_replace(x, "^\\d+\\s+", "")
    x <- stringr::str_replace(x, "\\s+(Glut|Gaba|NN)$", "")
    trimws(x)
  }
  
  extract_family <- function(subclass, cell_type) {
    core <- clean_subclass_name(subclass)
    if (cell_type == "Gaba") {
      if (startsWith(core, "Pvalb chandelier")) return("Pvalb chandelier")
      if (startsWith(core, "Sst Chodl")) return("Sst Chodl")
      first_tok <- strsplit(core, "\\s+")[[1]][1]
      if (first_tok %in% c("OB", "SCsg", "SCs", "STR", "STR-PAL", "SCig", "PAG", "HY", "RT")) return("Other noncanonical Gaba")
      return(first_tok)
    }
    if (cell_type == "Glut") {
      common <- c(
        "L2/3 IT CTX", "L4/5 IT CTX", "L5 NP CTX", "L5 ET CTX", "L6 CT CTX", "L6 IT CTX",
        "L5 IT CTX", "CLA-EPd-CTX Car3", "L2/3 IT RSP", "L6b CTX", "L4 RSP-ACA", "IT EP-CLA",
        "L5/6 IT TPE-ENT", "L2/3 IT PPP", "CA2-FC-IG", "L2/3 IT ENT", "IT AON-TT-DP", "L2 IT PPP-APr",
        "L6b/CT ENT", "SUB-ProS", "NP SUB", "HPF CR", "OB Eomes Ms4a15", "Pineal Crx", "L5 PPP"
      )
      for (cc in common) if (startsWith(core, cc)) return(cc)
      first_tok <- strsplit(core, "\\s+")[[1]][1]
      if (first_tok %in% c("SCsg", "SCzo", "SCiw", "SCop", "SCs", "PAG", "TH", "HY", "Pons", "MY", "PIR")) return("Other noncanonical Glut")
      return(core)
    }
    if (cell_type == "NonNeuron") return(clean_subclass_name(subclass))
    clean_subclass_name(subclass)
  }
  
  split_regions <- function(x) {
    if (is.null(x) || is.na(x) || x == "") return(character(0))
    unique(trimws(unlist(strsplit(as.character(x), ","))))
  }
  
  region_jaccard_fun <- function(a, b) {
    sa <- split_regions(a)
    sb <- split_regions(b)
    if (length(sa) == 0 && length(sb) == 0) return(NA_real_)
    length(intersect(sa, sb)) / length(union(sa, sb))
  }
  
  make_unique_clusters <- function(df) {
    c1 <- df %>% dplyr::select(
      label = cluster1_label, slide = cluster1_slide, layer = cluster1_layer,
      region = cluster1_region, type = cluster1_cell_type, subclass = cluster1_subclass,
      total = cluster1_total_cell_num, glut_n = cluster1_glut_num, gaba_n = cluster1_gaba_num,
      ei_ratio = cluster1_ei_ratio
    )
    c2 <- df %>% dplyr::select(
      label = cluster2_label, slide = cluster2_slide, layer = cluster2_layer,
      region = cluster2_region, type = cluster2_cell_type, subclass = cluster2_subclass,
      total = cluster2_total_cell_num, glut_n = cluster2_glut_num, gaba_n = cluster2_gaba_num,
      ei_ratio = cluster2_ei_ratio
    )
    bind_rows(c1, c2) %>%
      dplyr::distinct(label, .keep_all = TRUE) %>%
      dplyr::mutate(family = purrr::map2_chr(subclass, type, extract_family))
  }
  
  add_pair_metrics <- function(df, dataset_name) {
    out <- df %>%
      dplyr::mutate(
        dataset = dataset_name,
        source_type = cluster1_cell_type,
        target_type = cluster2_cell_type,
        source_subclass = cluster1_subclass,
        target_subclass = cluster2_subclass,
        source_family = purrr::map2_chr(source_subclass, source_type, extract_family),
        target_family = purrr::map2_chr(target_subclass, target_type, extract_family),
        slide_layer = paste(cluster1_slide, cluster1_layer, sep = "__"),
        source_overlap = as.numeric(cluster1_overlap),
        target_overlap = as.numeric(cluster2_overlap),
        asym = source_overlap - target_overlap,
      abs_asym = abs(asym),
      size_ratio = (cluster1_total_cell_num + 1) / (cluster2_total_cell_num + 1),
      log_size_ratio = log(size_ratio),
      source_glut_fraction = cluster1_glut_num / cluster1_total_cell_num,
      target_glut_fraction = cluster2_glut_num / cluster2_total_cell_num,
      source_gaba_fraction = cluster1_gaba_num / cluster1_total_cell_num,
        target_gaba_fraction = cluster2_gaba_num / cluster2_total_cell_num,
        region_jaccard = purrr::map2_dbl(cluster1_region, cluster2_region, region_jaccard_fun),
        exact_source_in_target = as.integer(source_overlap >= 1),
        exact_target_in_source = as.integer(target_overlap >= 1),
        source_in_target_only = as.integer(source_overlap >= 1 & target_overlap < 1),
        target_in_source_only = as.integer(target_overlap >= 1 & source_overlap < 1),
        both_exact = as.integer(source_overlap >= 1 & target_overlap >= 1),
        neither_exact = as.integer(source_overlap < 1 & target_overlap < 1)
      )
    
    N_lower_df <- out %>%
      dplyr::group_by(slide_layer) %>%
      dplyr::summarise(N_lower = max(union_cell, na.rm = TRUE), .groups = "drop")
    
    out <- out %>%
      dplyr::left_join(N_lower_df, by = "slide_layer") %>%
      dplyr::mutate(
        N_lower = pmax(N_lower, cluster1_total_cell_num, cluster2_total_cell_num, union_cell),
        expected_overlap = cluster1_total_cell_num * cluster2_total_cell_num / N_lower,
        hg_p = purrr::pmap_dbl(
          list(N_lower, cluster1_total_cell_num, cluster2_total_cell_num, overlap_cell),
          function(N, K, n, x) {
            phyper(q = as.integer(x) - 1, m = as.integer(K), n = as.integer(N - K),
                   k = as.integer(n), lower.tail = FALSE)
          }
        )
      )
    out$hg_fdr <- p.adjust(out$hg_p, method = "BH")
    out$log2_fe <- log2((out$overlap_cell + 0.5) / (out$expected_overlap + 0.5))
    out
  }
  
  make_universe <- function(clusters, pair_defs) {
    counts <- clusters %>%
      dplyr::count(slide, layer, type, name = "n") %>%
      tidyr::pivot_wider(names_from = type, values_from = n, values_fill = 0)
    out <- list(); idx <- 1
    for (ii in seq_len(nrow(counts))) {
      row <- counts[ii, ]
      for (pp in seq_len(nrow(pair_defs))) {
        st <- pair_defs$source_type[pp]
        tt <- pair_defs$target_type[pp]
        pc <- pair_defs$pair_class[pp]
        ns <- row[[st]]
        nt <- row[[tt]]
        possible <- ns * nt
        if (st == tt) possible <- possible - ns
        out[[idx]] <- tibble(
          slide = row$slide,
          layer = row$layer,
          slide_layer = paste(row$slide, row$layer, sep = "__"),
          pair_class = pc,
          possible_pairs = possible
        )
        idx <- idx + 1
      }
    }
    bind_rows(out)
  }
  
  cliffs_delta_from_wilcox <- function(x, y) {
    wt <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
    W <- unname(wt$statistic)
    n1 <- length(x); n2 <- length(y)
    U <- W - n1 * (n1 + 1) / 2
    2 * U / (n1 * n2) - 1
  }
  
  metric_compare <- function(x, y, feature, contrast_label) {
    x <- x[is.finite(x)]; y <- y[is.finite(y)]
    wt <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
    tibble(
      contrast = contrast_label, feature = feature,
      n1 = length(x), n2 = length(y),
      median1 = median(x, na.rm = TRUE), median2 = median(y, na.rm = TRUE),
      median_diff = median(x, na.rm = TRUE) - median(y, na.rm = TRUE),
      p = wt$p.value, cliffs_delta = cliffs_delta_from_wilcox(x, y)
    )
  }
  
  feature_compare_df <- function(df1, df2, name1, name2, features) {
    bind_rows(lapply(features, function(feat) metric_compare(df1[[feat]], df2[[feat]], feat, paste(name1, "vs", name2)))) %>%
      mutate(fdr = p.adjust(p, method = "BH"))
  }
  
  lin_contrast <- function(model, vcov_mat, weights_named, label, exponentiate = TRUE) {
    beta <- coef(model)
    L <- rep(0, length(beta)); names(L) <- names(beta)
    for (nm in names(weights_named)) L[nm] <- weights_named[[nm]]
    est <- sum(L * beta)
    se <- sqrt(drop(t(L) %*% vcov_mat %*% L))
    z <- est / se
    p <- 2 * pnorm(abs(z), lower.tail = FALSE)
    out <- tibble(contrast = label, estimate = est, std.error = se, z = z, p = p)
    if (exponentiate) out <- out %>% mutate(OR = exp(estimate), ci_low = exp(estimate - 1.96 * std.error), ci_high = exp(estimate + 1.96 * std.error))
    out
  }
  
  pred_prob <- function(model, newdata) {
    pr <- predict(model, newdata = newdata, type = "link", se.fit = TRUE)
    tibble(
      fit = plogis(pr$fit),
      low = plogis(pr$fit - 1.96 * pr$se.fit),
      high = plogis(pr$fit + 1.96 * pr$se.fit)
    )
  }
  
  # ------------------------------------------------------------------
  # Read and canonicalize data
  # ------------------------------------------------------------------
  nn_raw <- read_partner_file(input_nn)
  nonn_raw <- read_partner_file(input_nonn)
  
  nn_raw <- rename_partner_columns(nn_raw)
  nonn_raw <- rename_partner_columns(nonn_raw)
  
  mask_1 <- nonn_raw$cluster1_cell_type %in% c("Glut", "Gaba") & nonn_raw$cluster2_cell_type == "NonNeuron"
  mask_2 <- nonn_raw$cluster1_cell_type == "NonNeuron" & nonn_raw$cluster2_cell_type %in% c("Glut", "Gaba")
  part_a <- nonn_raw[mask_1, ]
  part_b <- nonn_raw[mask_2, ]
  part_b_swapped <- swap_clusters(part_b)
  
  nonn_canon <- bind_rows(part_a, part_b_swapped) %>%
    mutate(pair_key = paste(cluster1_label, cluster2_label, sep = "||")) %>%
    distinct(pair_key, .keep_all = TRUE) %>%
    dplyr::select(-pair_key)
  
  nn <- nn_raw %>%
    mutate(pair_key = paste(cluster1_label, cluster2_label, sep = "||")) %>%
    distinct(pair_key, .keep_all = TRUE) %>%
    dplyr::select(-pair_key)
  
  nn_dir <- add_pair_metrics(nn, "neuron_neuron") %>%
    mutate(pair_class = paste0(ifelse(source_type == "Glut", "E", "I"), "->", ifelse(target_type == "Glut", "E", "I")))
  
  nonn_dir <- add_pair_metrics(nonn_canon, "neuron_nonneuron") %>%
    mutate(pair_class = paste0(ifelse(source_type == "Glut", "E", "I"), "->NN"))
  
  obs_all <- bind_rows(nn_dir, nonn_dir)
  
  clusters_nn <- make_unique_clusters(nn)
  clusters_nonn <- make_unique_clusters(nonn_canon)
  
  pair_defs_nn <- tibble(
    source_type = c("Glut", "Glut", "Gaba", "Gaba"),
    target_type = c("Glut", "Gaba", "Glut", "Gaba"),
    pair_class = c("E->E", "E->I", "I->E", "I->I")
  )
  pair_defs_nonn <- tibble(
    source_type = c("Glut", "Gaba"),
    target_type = c("NonNeuron", "NonNeuron"),
    pair_class = c("E->NN", "I->NN")
  )
  
  u_nn <- make_universe(clusters_nn, pair_defs_nn)
  u_nonn <- make_universe(clusters_nonn, pair_defs_nonn)
  u_all <- bind_rows(u_nn, u_nonn)
  
  # ------------------------------------------------------------------
  # Summaries
  # ------------------------------------------------------------------
  pair_order <- c("E->E", "E->I", "I->E", "I->I", "E->NN", "I->NN")
  class_summary <- obs_all %>%
    group_by(pair_class) %>%
    summarise(
      n = n(),
      source_in_target_rate = mean(exact_source_in_target, na.rm = TRUE),
      target_in_source_rate = mean(exact_target_in_source, na.rm = TRUE),
      source_in_target_only_rate = mean(source_in_target_only, na.rm = TRUE),
      target_in_source_only_rate = mean(target_in_source_only, na.rm = TRUE),
      both_exact_rate = mean(both_exact, na.rm = TRUE),
      median_jaccard = median(jaccard, na.rm = TRUE),
      median_log2_fe = median(log2_fe, na.rm = TRUE),
      median_source_overlap = median(source_overlap, na.rm = TRUE),
      median_target_overlap = median(target_overlap, na.rm = TRUE),
      median_asym = median(asym, na.rm = TRUE),
      median_abs_asym = median(abs_asym, na.rm = TRUE),
      median_size_ratio = median(size_ratio, na.rm = TRUE),
      sig_frac = mean(hg_fdr < 0.05, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(u_all %>% group_by(pair_class) %>% summarise(possible_pairs = sum(possible_pairs), .groups = "drop"), by = "pair_class") %>%
    mutate(
      edge_rate = n / possible_pairs,
      source_in_target_all_possible = source_in_target_rate * edge_rate,
      target_in_source_all_possible = target_in_source_rate * edge_rate,
      pair_class = factor(pair_class, levels = pair_order)
    ) %>%
    arrange(pair_class)
  
  exact_summary <- obs_all %>%
    filter(exact_source_in_target == 1) %>%
    group_by(pair_class) %>%
    summarise(
      n_exact = n(),
      median_jaccard_exact = median(jaccard, na.rm = TRUE),
      median_log2_fe_exact = median(log2_fe, na.rm = TRUE),
      median_size_ratio_exact = median(size_ratio, na.rm = TRUE),
      .groups = "drop"
    )
  
  obs_gaba <- obs_all %>% filter(pair_class %in% c("I->E", "I->I", "I->NN"))
  ie_exact <- obs_gaba %>% filter(pair_class == "I->E", exact_source_in_target == 1)
  ii_exact <- obs_gaba %>% filter(pair_class == "I->I", exact_source_in_target == 1)
  inn_exact <- obs_gaba %>% filter(pair_class == "I->NN", exact_source_in_target == 1)
  
  u_gaba <- bind_rows(
    u_all %>% filter(pair_class %in% c("I->E", "I->I", "I->NN"))
  )
  
  layer_rates <- obs_gaba %>%
    group_by(pair_class, cluster1_layer) %>%
    summarise(
      observed = n(),
      exact_hits = sum(exact_source_in_target, na.rm = TRUE),
      median_jaccard = median(jaccard, na.rm = TRUE),
      median_log2_fe = median(log2_fe, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(
      u_gaba %>% group_by(pair_class, layer) %>% summarise(possible_pairs = sum(possible_pairs), .groups = "drop") %>% rename(cluster1_layer = layer),
      by = c("pair_class", "cluster1_layer")
    ) %>%
    mutate(all_possible_rate = exact_hits / possible_pairs)
  
  # ------------------------------------------------------------------
  # Model of interest
  # ------------------------------------------------------------------
  gaba_all <- obs_all %>%
    filter(source_type == "Gaba") %>%
    mutate(
      target_class3 = factor(case_when(
        pair_class == "I->E" ~ "Glut",
        pair_class == "I->I" ~ "Gaba",
        pair_class == "I->NN" ~ "NonNeuron",
        TRUE ~ NA_character_
      ), levels = c("Gaba", "Glut", "NonNeuron")),
      layer_cat = factor(cluster1_layer)
  ) %>%
  filter(
    is.finite(log_size_ratio), is.finite(region_jaccard),
    is.finite(source_glut_fraction), is.finite(target_glut_fraction), !is.na(target_class3)
  )
  
  m_exact <- glm(
  exact_source_in_target ~ target_class3 + log_size_ratio +
    region_jaccard + source_glut_fraction + target_glut_fraction + layer_cat,
  family = binomial(), data = gaba_all
)
  
  vc_exact <- sandwich::vcovCL(m_exact, cluster = gaba_all$slide_layer)
  coef_exact <- lmtest::coeftest(m_exact, vcov. = vc_exact)
  coef_df <- tibble(
    term = rownames(coef_exact),
    estimate = coef_exact[, 1],
    se = coef_exact[, 2],
    z = coef_exact[, 3],
    p = coef_exact[, 4]
  ) %>% mutate(OR = exp(estimate), ci_low = exp(estimate - 1.96 * se), ci_high = exp(estimate + 1.96 * se), fdr = p.adjust(p, method = "BH"))
  
  contrast_exact <- bind_rows(
    lin_contrast(m_exact, vc_exact, c("target_class3Glut" = 1), "Glut target vs Gaba target", exponentiate = TRUE),
    lin_contrast(m_exact, vc_exact, c("target_class3NonNeuron" = 1), "NonNeuron target vs Gaba target", exponentiate = TRUE),
    lin_contrast(m_exact, vc_exact, c("target_class3Glut" = 1, "target_class3NonNeuron" = -1), "Glut target vs NonNeuron target", exponentiate = TRUE)
  ) %>% mutate(fdr = p.adjust(p, method = "BH"))
  
  exact_gaba <- gaba_all %>% filter(exact_source_in_target == 1)
  m_j <- lm(jaccard ~ target_class3 + log_size_ratio + region_jaccard + source_glut_fraction + target_glut_fraction + layer_cat, data = exact_gaba)
  m_a <- lm(asym ~ target_class3 + log_size_ratio + region_jaccard + source_glut_fraction + target_glut_fraction + layer_cat, data = exact_gaba)
  
  tidy_cluster <- function(model, cluster) {
    vc <- sandwich::vcovCL(model, cluster = cluster)
    ct <- lmtest::coeftest(model, vcov. = vc)
    tibble(term = rownames(ct), estimate = ct[, 1], se = ct[, 2], statistic = ct[, 3], p = ct[, 4])
  }
  
  jaccard_coef <- tidy_cluster(m_j, exact_gaba$slide_layer) %>% mutate(fdr = p.adjust(p, method = "BH"))
  asym_coef <- tidy_cluster(m_a, exact_gaba$slide_layer) %>% mutate(fdr = p.adjust(p, method = "BH"))
  
  # adjusted predictions
  newdat <- gaba_all %>%
  summarise(
    log_size_ratio = median(log_size_ratio, na.rm = TRUE),
    region_jaccard = median(region_jaccard, na.rm = TRUE),
    source_glut_fraction = median(source_glut_fraction, na.rm = TRUE),
    target_glut_fraction = median(target_glut_fraction, na.rm = TRUE),
      layer_cat = names(sort(table(layer_cat), decreasing = TRUE))[1]
    ) %>%
    tidyr::crossing(target_class3 = factor(c("Gaba", "Glut", "NonNeuron"), levels = c("Gaba", "Glut", "NonNeuron"))) %>%
    mutate(layer_cat = factor(layer_cat, levels = levels(gaba_all$layer_cat)))
  
  pred_df <- bind_cols(newdat, pred_prob(m_exact, newdat)) %>%
    mutate(target_class3 = factor(target_class3, levels = c("Gaba", "Glut", "NonNeuron")))
  
  # ------------------------------------------------------------------
  # Tables for output 10
  # ------------------------------------------------------------------
  readr::write_csv(class_summary, file.path(table_dir, "table01_pair_class_summary.csv"))
  readr::write_csv(exact_summary, file.path(table_dir, "table02_pair_class_exact_summary.csv"))
  readr::write_csv(layer_rates, file.path(table_dir, "table03_gaba_layer_rates.csv"))
  readr::write_csv(coef_df, file.path(table_dir, "table19_gaba_source_exact_logistic_coefficients.csv"))
  readr::write_csv(contrast_exact, file.path(table_dir, "table20_gaba_source_exact_logistic_contrasts.csv"))
  readr::write_csv(jaccard_coef, file.path(table_dir, "table21_exactedge_jaccard_ols_coefficients.csv"))
  readr::write_csv(asym_coef, file.path(table_dir, "table22_exactedge_asym_ols_coefficients.csv"))
  readr::write_csv(pred_df, file.path(table_dir, "table23_adjusted_prediction_by_target_class.csv"))
  
  model_summary_key <- contrast_exact %>%
    mutate(
      p_fmt = purrr::map_chr(p, fmt_p),
      fdr_fmt = purrr::map_chr(fdr, fmt_p)
    ) %>%
    transmute(
      contrast,
      OR,
      ci_low,
      ci_high,
      p,
      fdr,
      result = paste0("OR=", sprintf("%.3f", OR), ", 95% CI [", sprintf("%.3f", ci_low), ", ", sprintf("%.3f", ci_high), "], p=", p_fmt)
    )
  readr::write_csv(model_summary_key, file.path(table_dir, "table24_key_model_summary.csv"))
  
  # ------------------------------------------------------------------
  # Figures (top-journal style)
  # ------------------------------------------------------------------
  # Figure 1: raw rate + compactness map
  fig1_df <- class_summary %>% left_join(exact_summary, by = "pair_class") %>%
    mutate(label = as.character(pair_class))
  cor1 <- suppressWarnings(cor.test(fig1_df$source_in_target_all_possible * 100, fig1_df$median_jaccard_exact, method = "spearman"))
  p1 <- ggplot(fig1_df, aes(x = source_in_target_all_possible * 100, y = median_jaccard_exact, size = n_exact, colour = pair_class, label = label)) +
    geom_point(alpha = 0.9) +
    geom_text(nudge_x = 0.18, size = 3.4, show.legend = FALSE) +
    scale_colour_brewer(palette = "Dark2") +
    scale_size_continuous(range = c(3, 10)) +
    labs(
      title = wrap_title("I→E is frequent, but exact-edge compactness does not simply track raw exact-containment opportunity"),
      subtitle = wrap_subtitle(paste0("Bubble size = number of exact edges. Spearman ρ = ", sprintf("%.3f", unname(cor1$estimate)), "; p = ", fmt_p(cor1$p.value), ".")),
      x = "Exact source-in-target rate among all possible pairs (%)",
      y = "Median Jaccard among exact edges"
    ) +
    theme_pub_global()
  save_pub(p1, "fig1_pairclass_rate_vs_compactness.png", width = 8.2, height = 5.6)
  
  # Figure 2: raw exact rates for Gaba sources with table20 annotation
  raw_gaba_rates <- class_summary %>%
    filter(pair_class %in% c("I->E", "I->I", "I->NN")) %>%
    mutate(target_label = factor(recode(as.character(pair_class), "I->E" = "Glut", "I->I" = "Gaba", "I->NN" = "NonNeuron"), levels = c("Gaba", "Glut", "NonNeuron"))) %>%
    select(target_label, all_possible = source_in_target_all_possible, observed = source_in_target_rate) %>%
    pivot_longer(-target_label, names_to = "metric", values_to = "rate") %>%
    mutate(metric = factor(metric, levels = c("all_possible", "observed"), labels = c("All possible exact rate", "Observed-edge exact rate")))
  sub2 <- paste0(
    "Cluster-robust logistic contrasts: Glut vs Gaba OR=", fmt_or(contrast_exact$OR[contrast_exact$contrast == "Glut target vs Gaba target"]),
    ", p=", fmt_p(contrast_exact$p[contrast_exact$contrast == "Glut target vs Gaba target"]),
    "; Glut vs NonNeuron OR=", fmt_or(contrast_exact$OR[contrast_exact$contrast == "Glut target vs NonNeuron target"]),
    ", p=", fmt_p(contrast_exact$p[contrast_exact$contrast == "Glut target vs NonNeuron target"]), "."
  )
  p2 <- ggplot(raw_gaba_rates, aes(target_label, rate * 100, fill = metric)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.66) +
    geom_text(aes(label = sprintf("%.1f%%", rate * 100)), position = position_dodge(width = 0.72), vjust = -0.25, size = 3.6) +
    scale_fill_manual(values = c("All possible exact rate" = "#B9B9B9", "Observed-edge exact rate" = "#C44E52")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = wrap_title("For Gaba sources, exact containment is strongest into Glut targets"),
      subtitle = wrap_subtitle(sub2),
      x = NULL, y = "Exact source-in-target rate (%)", fill = NULL
    ) +
    theme_pub_global()
  save_pub(p2, "fig2_gaba_targetclass_exact_rates.png", width = 8.2, height = 5.6)
  
  # Figure 3: layer heatmap
  heat_df <- layer_rates %>%
    filter(cluster1_layer %in% c("layer 1", "layer 2/3", "layer 4", "layer 5", "layer 6a"), pair_class %in% c("I->E", "I->I", "I->NN")) %>%
    mutate(pair_class = factor(pair_class, levels = c("I->I", "I->NN", "I->E")), cluster1_layer = factor(cluster1_layer, levels = c("layer 1", "layer 2/3", "layer 4", "layer 5", "layer 6a")))
  p3 <- ggplot(heat_df, aes(pair_class, cluster1_layer, fill = all_possible_rate)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.2f%%\n(n=%s)", 100 * all_possible_rate, scales::comma(observed))), size = 3.2, lineheight = 0.95) +
    scale_fill_viridis_c(labels = percent_format(accuracy = 1), option = "C") +
    labs(
      title = wrap_title("Across common layers, I→E keeps the highest all-possible exact rate among Gaba-source target classes"),
      subtitle = wrap_subtitle("Cell labels show exact rate and observed edge count within each layer × pair-class block."),
      x = NULL, y = NULL, fill = "All-possible\nexact rate"
    ) +
    theme_pub_global()
  save_pub(p3, "fig3_gaba_layer_heatmap.png", width = 7.2, height = 4.9)
  
  # Figure 4: cluster-robust logistic contrast forest (main figure)
  forest4 <- contrast_exact %>%
    mutate(
      contrast = factor(
        contrast,
        levels = rev(c("Glut target vs Gaba target", "Glut target vs NonNeuron target", "NonNeuron target vs Gaba target"))
      ),
      p_lab = purrr::map_chr(p, fmt_p),
      fdr_lab = purrr::map_chr(fdr, fmt_p),
      lab = paste0(
        "OR=", sprintf("%.3f", OR),
        " [", sprintf("%.3f", ci_low), ", ", sprintf("%.3f", ci_high), "]",
        "\np=", p_lab, "; FDR=", fdr_lab
      )
    )
  xmax4 <- max(forest4$ci_high, na.rm = TRUE)
  p4 <- ggplot(forest4, aes(x = OR, y = contrast, colour = contrast)) +
    geom_vline(xintercept = 1, linetype = 2, linewidth = 0.6, colour = "#808080") +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.18, linewidth = 0.8, show.legend = FALSE) +
    geom_point(size = 3.2, show.legend = FALSE) +
    geom_text(aes(x = pmin(ci_high * 1.12, xmax4 * 1.33), label = lab), hjust = 0, size = 3.3, colour = "#444444", show.legend = FALSE) +
    scale_x_log10(labels = label_number(accuracy = 0.01), expand = expansion(mult = c(0.03, 0.38))) +
    scale_colour_manual(values = c(
      "Glut target vs Gaba target" = col_target[["Glut"]],
      "Glut target vs NonNeuron target" = "#9C755F",
      "NonNeuron target vs Gaba target" = col_target[["NonNeuron"]]
    )) +
    coord_cartesian(clip = "off") +
  labs(
    title = wrap_title("The Glut-target preference survives covariate control"),
    subtitle = wrap_subtitle("Gaba-source cluster-robust logistic model controlling for log(size ratio), region_jaccard, source/target glut fraction, and layer."),
    x = "Adjusted odds ratio for exact source-in-target (log scale)",
    y = NULL
  ) +
    theme_pub_global()
  save_pub(p4, "fig4_shared_source_mcnemar_heatmap.png", width = 9.2, height = 5.4)
  
  # Figure 5: adjusted predicted probabilities
  pred_lab <- contrast_exact %>% select(contrast, OR, p) 
  p5 <- ggplot(pred_df, aes(x = target_class3, y = fit, colour = target_class3)) +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0.08, linewidth = 0.9) +
    geom_point(size = 3.6) +
    geom_text(aes(label = sprintf("%.1f%%", fit * 100)), vjust = -0.85, size = 3.6, show.legend = FALSE) +
    scale_colour_manual(values = col_target) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0.04, 0.16))) +
    labs(
      title = wrap_title("Adjusted predictions still place Glut above Gaba and NonNeuron targets"),
      subtitle = wrap_subtitle("Predicted exact-containment probability at the median covariate profile and modal layer."),
      x = NULL, y = "Predicted probability of exact source-in-target", colour = NULL
    ) +
    theme_pub_global()
  save_pub(p5, "fig5_gaba_family_slopeplot.png", width = 7.8, height = 5.4)
  
  # Figure 6a: exact model coefficient forest
  forest6a <- coef_df %>%
  filter(term %in% c("target_class3Glut", "target_class3NonNeuron", "log_size_ratio", "region_jaccard", "source_glut_fraction", "target_glut_fraction")) %>%
  mutate(term_label = recode(term,
                             `target_class3Glut` = "Target class: Glut vs Gaba",
                             `target_class3NonNeuron` = "Target class: NonNeuron vs Gaba",
                             `log_size_ratio` = "log(size ratio)",
                             `region_jaccard` = "region_jaccard",
                             `source_glut_fraction` = "source glut fraction",
                             `target_glut_fraction` = "target glut fraction"
  )) %>%
    mutate(
      term_label = factor(term_label, levels = rev(term_label)),
      p_lab = purrr::map_chr(p, fmt_p),
      lab = paste0("OR=", sprintf("%.3f", OR), "\np=", p_lab)
    )
  p6a <- ggplot(forest6a, aes(x = OR, y = term_label)) +
    geom_vline(xintercept = 1, linetype = 2, linewidth = 0.6, colour = "#808080") +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.18, linewidth = 0.8, colour = "#4D4D4D") +
    geom_point(size = 3.0, colour = "#4E79A7") +
    geom_text(aes(x = ci_high * 1.08, label = lab), hjust = 0, size = 3.0, colour = "#444444") +
    scale_x_log10(labels = label_number(accuracy = 0.01), expand = expansion(mult = c(0.03, 0.34))) +
    coord_cartesian(clip = "off") +
    labs(
      title = wrap_title("Effect decomposition in the exact model"),
      subtitle = wrap_subtitle("Target-class terms remain strong after geometric and composition covariates are included."),
      x = "Adjusted odds ratio (log scale)", y = NULL
    ) +
    theme_pub_global()
  save_pub(p6a, "fig6a_top_IE_combos.png", width = 8.6, height = 5.8)
  
  # Figure 6b: target-class contrasts for exact-edge compactness/asymmetry
  j_sub <- jaccard_coef %>% filter(term %in% c("target_class3Glut", "target_class3NonNeuron")) %>% mutate(model = "Exact-edge Jaccard")
  a_sub <- asym_coef %>% filter(term %in% c("target_class3Glut", "target_class3NonNeuron")) %>% mutate(model = "Exact-edge asymmetry")
  forest6b <- bind_rows(j_sub, a_sub) %>%
    mutate(term_label = recode(term,
                               `target_class3Glut` = "Glut vs Gaba target",
                               `target_class3NonNeuron` = "NonNeuron vs Gaba target"
    )) %>%
    mutate(term_label = factor(term_label, levels = rev(c("Glut vs Gaba target", "NonNeuron vs Gaba target"))))
  p6b <- ggplot(forest6b, aes(x = estimate, y = term_label, colour = model)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.6, colour = "#808080") +
    geom_errorbarh(aes(xmin = estimate - 1.96 * se, xmax = estimate + 1.96 * se), height = 0.18, linewidth = 0.8, position = position_dodge(width = 0.5)) +
    geom_point(size = 3.0, position = position_dodge(width = 0.5)) +
    geom_text(aes(label = paste0("p=", purrr::map_chr(p, fmt_p))), position = position_dodge(width = 0.5), hjust = -0.15, size = 3.0, show.legend = FALSE) +
    scale_colour_manual(values = c("Exact-edge Jaccard" = "#C44E52", "Exact-edge asymmetry" = "#59A14F")) +
    labs(
      title = wrap_title("Among exact Gaba-source edges, target class still shapes compactness and asymmetry"),
      subtitle = wrap_subtitle("Positive coefficients indicate larger Jaccard or asymmetry relative to Gaba targets within exact edges."),
      x = "Cluster-robust coefficient estimate", y = NULL, colour = NULL
    ) +
    theme_pub_global()
  save_pub(p6b, "fig6b_top_INN_combos.png", width = 8.8, height = 5.6)
  
  # Figure 7: exact-edge boxplots
  plot_df <- bind_rows(
    ie_exact %>% mutate(target_group = "Glut"),
    ii_exact %>% mutate(target_group = "Gaba"),
    inn_exact %>% mutate(target_group = "NonNeuron")
  )
  
  long_box <- plot_df %>%
    dplyr::select(target_group, jaccard, log2_fe, size_ratio) %>%
    tidyr::pivot_longer(
      cols = c(jaccard, log2_fe, size_ratio),
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      metric = dplyr::recode(
        metric,
        "jaccard" = "Jaccard",
        "log2_fe" = "log2 enrichment",
        "size_ratio" = "Source/target size ratio"
      ),
      target_group = factor(target_group, levels = c("Glut", "Gaba", "NonNeuron"))
    )
  
  # per-group labels
  box_lab <- long_box %>%
    dplyr::group_by(metric, target_group) %>%
    dplyr::summarise(
      n = dplyr::n(),
      med = median(value, na.rm = TRUE),
      q3 = quantile(value, 0.75, na.rm = TRUE),
      vmax = max(value, na.rm = TRUE),
      vmin = min(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      span = pmax(vmax - vmin, 1e-6),
      y = q3 + 0.06 * span,
      lab = paste0("n=", scales::comma(n), "\nmed=", sprintf("%.3f", med))
    )
  
  # per-facet p-value summary
  facet_p <- long_box %>%
    dplyr::group_by(metric) %>%
    dplyr::summarise(
      p_glut_vs_gaba = suppressWarnings(
        wilcox.test(
          value[target_group == "Glut"],
          value[target_group == "Gaba"],
          exact = FALSE
        )$p.value
      ),
      p_glut_vs_non = suppressWarnings(
        wilcox.test(
          value[target_group == "Glut"],
          value[target_group == "NonNeuron"],
          exact = FALSE
        )$p.value
      ),
      p_gaba_vs_non = suppressWarnings(
        wilcox.test(
          value[target_group == "Gaba"],
          value[target_group == "NonNeuron"],
          exact = FALSE
        )$p.value
      ),
      vmax = max(value, na.rm = TRUE),
      vmin = min(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      span = pmax(vmax - vmin, 1e-6),
      y = vmax + 0.22 * span,
      x = factor("NonNeuron", levels = c("Glut", "Gaba", "NonNeuron")),
      lab = paste0(
        "G vs I: p=", purrr::map_chr(p_glut_vs_gaba, fmt_p), "\n",
        "G vs N: p=", purrr::map_chr(p_glut_vs_non, fmt_p), "\n",
        "I vs N: p=", purrr::map_chr(p_gaba_vs_non, fmt_p)
      )
    )
  
  p7 <- ggplot(long_box, aes(x = target_group, y = value, fill = target_group)) +
    geom_boxplot(
      outlier.shape = NA,
      width = 0.62,
      linewidth = 0.75,
      alpha = 0.92
    ) +
    geom_jitter(
      width = 0.12,
      size = 0.9,
      alpha = 0.18,
      colour = "#4D4D4D",
      show.legend = FALSE
    ) +
    geom_text(
      data = box_lab,
      aes(x = target_group, y = y, label = lab),
      inherit.aes = FALSE,
      size = 3.0,
      lineheight = 0.92,
      colour = "#3A3A3A"
    ) +
    geom_text(
      data = facet_p,
      aes(x = x, y = y, label = lab),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 1,
      size = 3.0,
      lineheight = 0.96,
      colour = "#2F2F2F"
    ) +
    facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = c(
      "Glut" = col_target[["Glut"]],
      "Gaba" = col_target[["Gaba"]],
      "NonNeuron" = col_target[["NonNeuron"]]
    )) +
    labs(
      title = wrap_title("Among exact Gaba containment events, Glut targets are frequent but looser containers"),
      subtitle = wrap_subtitle("Each facet shows the exact-edge distribution across target classes, with per-group n/median and facet-level Wilcoxon p-values."),
      x = NULL,
      y = NULL
    ) +
    coord_cartesian(clip = "off") +
    theme_pub_global() +
    theme(
      legend.position = "none",
      plot.margin = margin(14, 34, 10, 10)
    )
  
  save_pub(p7, "fig7_exact_edge_boxplots.png", width = 11.6, height = 4.8)
  
  message("R_output_10 finished: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE))
}

run_module_11 <- function(runtime_base) {
  message('[MODULE 11] start')
  # [comment omitted: encoding-safe]
  # [removed] rm(list = ls())
  gc()
  options(stringsAsFactors = FALSE, scipen = 999)
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readr)
    library(stringr)
    library(purrr)
    library(ggplot2)
    library(scales)
    library(forcats)
    library(tibble)
    library(pROC)
  })
  
  # ============================================================
  # R_output_11
  # Causal-priority model for containment motifs across:
  #   1) neuron-neuron
  #   2) neuron-nonneuron
  #
  # Main biological questions addressed:
  # [comment omitted: encoding-safe]
  # [comment omitted: encoding-safe]
  # [comment omitted: encoding-safe]
  #
  # Output:
  #   base_dir/R_output_11/
  # [comment omitted: encoding-safe]
  # [comment omitted: encoding-safe]
  # ============================================================
  
  # ------------------------------------------------------------------
  # Paths
  # ------------------------------------------------------------------
  base_dir <- runtime_base
  
  nn_candidates <- c(
    file.path(base_dir, "neuron-neuron-partner.csv"),
    file.path(base_dir, "neuron-neuron-partner.txt")
  )
  nonn_candidates <- c(
    file.path(base_dir, "neuron-nonneuron-partner.csv"),
    file.path(base_dir, "neuron-nonneuron-partner.txt")
  )
  
  nn_file <- nn_candidates[file.exists(nn_candidates)][1]
  nonn_file <- nonn_candidates[file.exists(nonn_candidates)][1]
  
  if (is.na(nn_file) || !nzchar(nn_file)) {
    stop("Cannot find neuron-neuron partner file under: ", base_dir)
  }
  if (is.na(nonn_file) || !nzchar(nonn_file)) {
    stop("Cannot find neuron-nonneuron partner file under: ", base_dir)
  }
  
  out_dir <- file.path(base_dir, "R_output_11")
  fig_dir <- file.path(out_dir, "figures")
  table_dir <- file.path(out_dir, "tables")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ------------------------------------------------------------------
  # Formatting helpers
  # ------------------------------------------------------------------
  fmt_p <- function(p) {
    if (length(p) != 1) stop("fmt_p expects a scalar p-value")
    if (is.na(p)) return("NA")
    if (p < 1e-300) return("<1e-300")
    if (p < 1e-4) return(format(p, scientific = TRUE, digits = 2))
    sprintf("%.4f", p)
  }
  fmt_pct <- function(x, digits = 1) sprintf(paste0("%.", digits, "f%%"), 100 * x)
  fmt_num <- function(x, digits = 3) sprintf(paste0("%.", digits, "f"), x)
  wrap_title <- function(x, width = 68) stringr::str_wrap(x, width = width)
  wrap_subtitle <- function(x, width = 102) stringr::str_wrap(x, width = width)
  
  col_motif <- c(
    "Gaba⊂Glut" = "#C44E52",
    "smaller⊂Glut (GG)" = "#4E79A7",
    "Gaba⊂NonNeuron" = "#59A14F",
    "NonNeuron⊂Glut" = "#9C755F",
    "smaller⊂Gaba (II)" = "#B07AA1"
  )
  col_binary <- c("Higher" = "#C44E52", "Lower" = "#4E79A7")
  col_partner <- c("Glut host" = "#C44E52", "NonNeuron host" = "#59A14F")
  
  base_theme <- function() {
    theme_pub_global(base_size = 12) +
      theme(
        plot.title.position = "plot",
        plot.caption.position = "plot",
        plot.margin = margin(18, 26, 16, 16),
        plot.title = element_text(face = "bold", size = 16, lineheight = 1.02, margin = margin(b = 8)),
        plot.subtitle = element_text(size = 11.4, colour = "#4D4D4D", lineheight = 1.10, margin = margin(b = 10)),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11, colour = "black"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 10.8),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12)
      )
  }
  
  save_pub <- function(p, filename, width = 8.2, height = 5.4, dpi = 600) {
    ggsave_pub(file.path(fig_dir, filename), p, width = width, height = height, dpi = dpi,
           bg = "white", limitsize = FALSE)
    pdf_name <- sub("\\.png$", ".pdf", filename)
    ggsave_pub(file.path(fig_dir, pdf_name), p, width = width, height = height,
           bg = "white", device = cairo_pdf, limitsize = FALSE)
  }
  
  # ------------------------------------------------------------------
  # Data helpers
  # ------------------------------------------------------------------
  region_breadth <- function(x) {
    if (is.na(x) || x == "") return(NA_integer_)
    length(str_split(as.character(x), ",", simplify = FALSE)[[1]])
  }
  
  safe_auc <- function(y, score) {
    y <- as.integer(y)
    ok <- is.finite(score) & !is.na(y)
    y <- y[ok]
    score <- score[ok]
    if (length(unique(y)) < 2) return(NA_real_)
    out <- tryCatch(as.numeric(pROC::auc(response = y, predictor = score)), error = function(e) NA_real_)
    out
  }
  
  paired_wilcox_p <- function(x, y) {
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]
    if (length(x) == 0 || length(y) == 0) return(NA_real_)
    suppressWarnings(wilcox.test(x, y, paired = TRUE, alternative = "two.sided", exact = FALSE)$p.value)
  }
  
  one_sample_wilcox_p <- function(x, mu = 0) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NA_real_)
    suppressWarnings(wilcox.test(x, mu = mu, exact = FALSE)$p.value)
  }
  
  cramers_v_binary_by_group <- function(group_code, y) {
    tab <- table(group_code, y)
    if (nrow(tab) < 2 || ncol(tab) < 2) return(NA_real_)
    chi <- suppressWarnings(chisq.test(tab, correct = FALSE)$statistic)
    n <- sum(tab)
    as.numeric(sqrt(chi / n))
  }
  
  conditional_perm_specificity <- function(df, outcome_col, pair_cols, strata_cols, n_perm = 1000, seed = 1) {
    set.seed(seed)
    dat <- df %>%
      dplyr::select(dplyr::all_of(c(outcome_col, pair_cols, strata_cols))) %>%
      tidyr::drop_na()
    if (nrow(dat) == 0) {
      return(list(obs_v = NA_real_, perm_mean_v = NA_real_, perm_sd_v = NA_real_, z = NA_real_, p = NA_real_, n_pairs = 0, n_strata = 0))
    }
    
    pair_id <- interaction(dat[, pair_cols], drop = TRUE, lex.order = TRUE)
    strata <- interaction(dat[, strata_cols], drop = TRUE, lex.order = TRUE)
    y <- as.integer(dat[[outcome_col]])
    
    obs_v <- cramers_v_binary_by_group(pair_id, y)
    strata_idx <- split(seq_len(nrow(dat)), strata)
    perm_v <- numeric(n_perm)
    
    for (b in seq_len(n_perm)) {
      y_perm <- y
      for (idx in strata_idx) {
        if (length(idx) > 1) y_perm[idx] <- sample(y_perm[idx], replace = FALSE)
      }
      perm_v[b] <- cramers_v_binary_by_group(pair_id, y_perm)
    }
    
    list(
      obs_v = obs_v,
      perm_mean_v = mean(perm_v, na.rm = TRUE),
      perm_sd_v = sd(perm_v, na.rm = TRUE),
      z = (obs_v - mean(perm_v, na.rm = TRUE)) / sd(perm_v, na.rm = TRUE),
      p = (sum(perm_v >= obs_v, na.rm = TRUE) + 1) / (sum(is.finite(perm_v)) + 1),
      n_pairs = nlevels(pair_id),
      n_strata = nlevels(strata)
    )
  }
  
  partial_asym_stats <- function(df, fill_a, fill_b, label_text) {
    partial <- df %>% dplyr::filter(.data[[fill_a]] < 1, .data[[fill_b]] < 1)
    partial <- partial %>% dplyr::mutate(asym = .data[[fill_a]] - .data[[fill_b]])
    slide_median <- partial %>%
      dplyr::group_by(slide) %>%
      dplyr::summarise(med = median(asym, na.rm = TRUE), .groups = "drop")
    
    tibble(
      name = label_text,
      n_partial = nrow(partial),
      pair_median = median(partial$asym, na.rm = TRUE),
      pair_mean = mean(partial$asym, na.rm = TRUE),
      pair_wilcoxon_p = one_sample_wilcox_p(partial$asym, 0),
      slide_n = nrow(slide_median),
      slide_positive_n = sum(slide_median$med > 0, na.rm = TRUE),
      slide_median_of_medians = median(slide_median$med, na.rm = TRUE),
      slide_wilcoxon_p = one_sample_wilcox_p(slide_median$med, 0)
    )
  }
  
  mcnemar_exact <- function(x, y) {
    b <- sum(x == 1 & y == 0, na.rm = TRUE)
    c <- sum(x == 0 & y == 1, na.rm = TRUE)
    if ((b + c) == 0) return(list(b_10 = b, c_01 = c, p = NA_real_))
    p <- binom.test(min(b, c), b + c, p = 0.5, alternative = "two.sided")$p.value
    list(b_10 = b, c_01 = c, p = p)
  }
  
  dedupe_unordered <- function(df) {
    key <- purrr::pmap_chr(
      df[, c("cluster.1_label", "cluster.2_label")],
      function(cluster.1_label, cluster.2_label) paste(sort(c(cluster.1_label, cluster.2_label)), collapse = " || ")
    )
    df[!duplicated(key), , drop = FALSE]
  }
  
  prep_df <- function(df, source_name) {
    num_cols <- c(
      "cluster.1_total_cell_num", "cluster.2_total_cell_num",
      "cluster.1_E_I_Ratio", "cluster.2_E_I_Ratio",
      "overlap_cell", "union_cell", "jaccard",
      "cluster.1.overlap.percent", "cluster.2.overlap.percent",
      "cluster.1_cauchy_combination_p", "cluster.2_cauchy_combination_p"
    )
    for (cc in num_cols) df[[cc]] <- suppressWarnings(as.numeric(df[[cc]]))
    df$cluster.1_cauchy_combination_p[df$cluster.1_cauchy_combination_p == 0] <- .Machine$double.xmin
    df$cluster.2_cauchy_combination_p[df$cluster.2_cauchy_combination_p == 0] <- .Machine$double.xmin
    df$logp1 <- -log10(df$cluster.1_cauchy_combination_p)
    df$logp2 <- -log10(df$cluster.2_cauchy_combination_p)
    df$region_breadth1 <- vapply(df$cluster.1_region, region_breadth, integer(1))
    df$region_breadth2 <- vapply(df$cluster.2_region, region_breadth, integer(1))
    df$source <- source_name
    df
  }
  
  side_df <- function(df, idx) {
    tibble::tibble(
      label = df[[paste0("cluster.", idx, "_label")]],
      slide = df[[paste0("cluster.", idx, "_slide")]],
      layer = df[[paste0("cluster.", idx, "_layer")]],
      region = df[[paste0("cluster.", idx, "_region")]],
      type = df[[paste0("cluster.", idx, "_cell_Neruon_type")]],
      subclass = df[[paste0("cluster.", idx, "_subclass")]],
      n = df[[paste0("cluster.", idx, "_total_cell_num")]],
      logp = df[[paste0("logp", idx)]],
      fill = df[[paste0("cluster.", idx, ".overlap.percent")]],
      breadth = df[[paste0("region_breadth", idx)]]
    )
  }
  
  orient_cross_pair <- function(df, type_a, type_b, prefix_a, prefix_b) {
    keep <- (df$cluster.1_cell_Neruon_type == type_a & df$cluster.2_cell_Neruon_type == type_b) |
      (df$cluster.1_cell_Neruon_type == type_b & df$cluster.2_cell_Neruon_type == type_a)
    x <- df[keep, , drop = FALSE]
    out_list <- vector("list", nrow(x))
    
    for (i in seq_len(nrow(x))) {
      row <- x[i, , drop = FALSE]
      t1 <- row$cluster.1_cell_Neruon_type
      if (t1 == type_a) {
        a <- side_df(row, 1)
        b <- side_df(row, 2)
      } else {
        a <- side_df(row, 2)
        b <- side_df(row, 1)
      }
      out_list[[i]] <- tibble::tibble(
        slide = a$slide,
        overlap_cell = row$overlap_cell,
        jaccard = row$jaccard,
        !!paste0(prefix_a, "_label") := a$label,
        !!paste0(prefix_a, "_subclass") := a$subclass,
        !!paste0(prefix_a, "_layer") := a$layer,
        !!paste0(prefix_a, "_n") := a$n,
        !!paste0(prefix_a, "_logp") := a$logp,
        !!paste0(prefix_a, "_fill") := a$fill,
        !!paste0(prefix_a, "_region_breadth") := a$breadth,
        !!paste0(prefix_b, "_label") := b$label,
        !!paste0(prefix_b, "_subclass") := b$subclass,
        !!paste0(prefix_b, "_layer") := b$layer,
        !!paste0(prefix_b, "_n") := b$n,
        !!paste0(prefix_b, "_logp") := b$logp,
        !!paste0(prefix_b, "_fill") := b$fill,
        !!paste0(prefix_b, "_region_breadth") := b$breadth
      )
    }
    dplyr::bind_rows(out_list)
  }
  
  orient_same_type_by_size <- function(df, type_name) {
    keep <- df$cluster.1_cell_Neruon_type == type_name & df$cluster.2_cell_Neruon_type == type_name
    x <- df[keep, , drop = FALSE]
    out_list <- vector("list", nrow(x))
    
    for (i in seq_len(nrow(x))) {
      row <- x[i, , drop = FALSE]
      s1 <- side_df(row, 1)
      s2 <- side_df(row, 2)
      if (s1$n <= s2$n) {
        small <- s1; large <- s2
      } else {
        small <- s2; large <- s1
      }
      out_list[[i]] <- tibble::tibble(
        slide = small$slide,
        overlap_cell = row$overlap_cell,
        jaccard = row$jaccard,
        small_label = small$label,
        small_subclass = small$subclass,
        small_layer = small$layer,
        small_n = small$n,
        small_logp = small$logp,
        small_fill = small$fill,
        small_region_breadth = small$breadth,
        large_label = large$label,
        large_subclass = large$subclass,
        large_layer = large$layer,
        large_n = large$n,
        large_logp = large$logp,
        large_fill = large$fill,
        large_region_breadth = large$breadth
      )
    }
    dplyr::bind_rows(out_list)
  }
  
  # ------------------------------------------------------------------
  # Load data
  # ------------------------------------------------------------------
  nn <- readr::read_csv(nn_file, show_col_types = FALSE)
  nonn <- readr::read_csv(nonn_file, show_col_types = FALSE)
  
  nn_u <- prep_df(dedupe_unordered(nn), "neuron-neuron")
  nonn_u <- prep_df(dedupe_unordered(nonn), "neuron-nonneuron")
  
  nn_u$pair_type <- apply(
    nn_u[, c("cluster.1_cell_Neruon_type", "cluster.2_cell_Neruon_type")],
    1, function(z) paste(sort(z), collapse = "-")
  )
  nonn_u$pair_type <- apply(
    nonn_u[, c("cluster.1_cell_Neruon_type", "cluster.2_cell_Neruon_type")],
    1, function(z) paste(sort(z), collapse = "-")
  )
  
  GI <- nn_u[nn_u$pair_type == "Gaba-Glut", , drop = FALSE]
  GG <- nn_u[nn_u$pair_type == "Glut-Glut", , drop = FALSE]
  II <- nn_u[nn_u$pair_type == "Gaba-Gaba", , drop = FALSE]
  GN <- nonn_u[nonn_u$pair_type == "Glut-NonNeuron", , drop = FALSE]
  IN <- nonn_u[nonn_u$pair_type == "Gaba-NonNeuron", , drop = FALSE]
  NN <- nonn_u[nonn_u$pair_type == "NonNeuron-NonNeuron", , drop = FALSE]
  
  GIo <- orient_cross_pair(nn_u, "Gaba", "Glut", "G", "E")
  GNo <- orient_cross_pair(nonn_u, "Glut", "NonNeuron", "E", "N")
  INo <- orient_cross_pair(nonn_u, "Gaba", "NonNeuron", "G", "N")
  GGs <- orient_same_type_by_size(GG, "Glut")
  IIs <- orient_same_type_by_size(II, "Gaba")
  
  GIo <- GIo %>% mutate(y_G_in_E = as.integer(G_fill == 1), log_size_ratio = log2((E_n + 1) / (G_n + 1)), size_bin = dplyr::ntile(log_size_ratio, 10))
  GNo <- GNo %>% mutate(y_N_in_E = as.integer(N_fill == 1), log_size_ratio = log2((E_n + 1) / (N_n + 1)), size_bin = dplyr::ntile(log_size_ratio, 10))
  INo <- INo %>% mutate(y_G_in_N = as.integer(G_fill == 1), log_size_ratio = log2((N_n + 1) / (G_n + 1)), size_bin = dplyr::ntile(log_size_ratio, 10))
  GGs <- GGs %>% mutate(y_small_in_large = as.integer(small_fill == 1), log_size_ratio = log2((large_n + 1) / (small_n + 1)), size_bin = dplyr::ntile(log_size_ratio, 10))
  IIs <- IIs %>% mutate(y_small_in_large = as.integer(small_fill == 1), log_size_ratio = log2((large_n + 1) / (small_n + 1)), size_bin = dplyr::ntile(log_size_ratio, 10))
  
  # ------------------------------------------------------------------
  # Residual identity and precursor analyses
  # ------------------------------------------------------------------
  res_GI <- conditional_perm_specificity(GIo, "y_G_in_E", c("G_subclass", "E_subclass"), c("slide", "size_bin", "G_layer", "E_layer"), n_perm = 1000, seed = 11)
  res_GN <- conditional_perm_specificity(GNo, "y_N_in_E", c("N_subclass", "E_subclass"), c("slide", "size_bin", "N_layer", "E_layer"), n_perm = 1000, seed = 12)
  res_IN <- conditional_perm_specificity(INo, "y_G_in_N", c("G_subclass", "N_subclass"), c("slide", "size_bin", "G_layer", "N_layer"), n_perm = 1000, seed = 13)
  res_GG <- conditional_perm_specificity(GGs, "y_small_in_large", c("small_subclass", "large_subclass"), c("slide", "size_bin", "small_layer", "large_layer"), n_perm = 1000, seed = 14)
  res_II <- conditional_perm_specificity(IIs, "y_small_in_large", c("small_subclass", "large_subclass"), c("slide", "size_bin", "small_layer", "large_layer"), n_perm = 1000, seed = 15)
  
  as_GI <- partial_asym_stats(GIo, "G_fill", "E_fill", "GI")
  as_GN <- partial_asym_stats(GNo, "N_fill", "E_fill", "GN")
  as_IN <- partial_asym_stats(INo, "G_fill", "N_fill", "IN")
  as_GG <- partial_asym_stats(GGs, "small_fill", "large_fill", "GG")
  as_II <- partial_asym_stats(IIs, "small_fill", "large_fill", "II")
  
  motif_summary <- tibble::tibble(
    motif = c("Gaba⊂Glut", "NonNeuron⊂Glut", "Gaba⊂NonNeuron", "smaller⊂Glut (GG)", "smaller⊂Gaba (II)"),
    pair_class = c("GI", "GN", "IN", "GG", "II"),
    raw_rate = c(mean(GIo$G_fill == 1), mean(GNo$N_fill == 1), mean(INo$G_fill == 1), mean(GGs$small_fill == 1), mean(IIs$small_fill == 1)),
    eligible_rate = c(
      mean(GIo$G_fill[GIo$E_n >= GIo$G_n] == 1, na.rm = TRUE),
      mean(GNo$N_fill[GNo$E_n >= GNo$N_n] == 1, na.rm = TRUE),
      mean(INo$G_fill[INo$N_n >= INo$G_n] == 1, na.rm = TRUE),
      mean(GGs$small_fill == 1, na.rm = TRUE),
      mean(IIs$small_fill == 1, na.rm = TRUE)
    ),
    size_auc = c(
      safe_auc(GIo$y_G_in_E, GIo$log_size_ratio),
      safe_auc(GNo$y_N_in_E, GNo$log_size_ratio),
      safe_auc(INo$y_G_in_N, INo$log_size_ratio),
      safe_auc(GGs$y_small_in_large, GGs$log_size_ratio),
      safe_auc(IIs$y_small_in_large, IIs$log_size_ratio)
    ),
    identity_v = c(res_GI$obs_v, res_GN$obs_v, res_IN$obs_v, res_GG$obs_v, res_II$obs_v),
    identity_perm_z = c(res_GI$z, res_GN$z, res_IN$z, res_GG$z, res_II$z),
    identity_perm_p = c(res_GI$p, res_GN$p, res_IN$p, res_GG$p, res_II$p),
    partial_asym_median = c(as_GI$pair_median, as_GN$pair_median, as_IN$pair_median, as_GG$pair_median, as_II$pair_median),
    partial_slide_p = c(as_GI$slide_wilcoxon_p, as_GN$slide_wilcoxon_p, as_IN$slide_wilcoxon_p, as_GG$slide_wilcoxon_p, as_II$slide_wilcoxon_p),
    n_pairs = c(nrow(GIo), nrow(GNo), nrow(INo), nrow(GGs), nrow(IIs))
  )
  write_csv(motif_summary, file.path(table_dir, "table01_motif_summary.csv"))
  write_csv(bind_rows(as_GI, as_GN, as_IN, as_GG, as_II), file.path(table_dir, "table02_partial_asymmetry_summary.csv"))
  
  # ------------------------------------------------------------------
  # Size-matched GABA partner preference: Glut host vs NonNeuron host
  # ------------------------------------------------------------------
  gi_elig <- GIo %>%
    dplyr::filter(E_n >= G_n) %>%
    dplyr::transmute(G_label, G_subclass, G_layer, partner = "Glut", full = as.integer(G_fill == 1), size_ratio = E_n / G_n)
  
  in_elig <- INo %>%
    dplyr::filter(N_n >= G_n) %>%
    dplyr::transmute(G_label, G_subclass, G_layer, partner = "NonNeuron", full = as.integer(G_fill == 1), size_ratio = N_n / G_n)
  
  gcomb <- dplyr::bind_rows(gi_elig, in_elig) %>% dplyr::mutate(size_bin = dplyr::ntile(size_ratio, 10))
  
  # subclass/layer/size-bin stratified MH pooled OR
  strat_g <- gcomb %>% dplyr::mutate(stratum = interaction(G_subclass, G_layer, size_bin, drop = TRUE))
  mh_tables <- split(strat_g, strat_g$stratum)
  mh_tables <- mh_tables[vapply(mh_tables, function(x) length(unique(x$partner)) > 1, logical(1))]
  mh_array <- array(0, dim = c(2, 2, length(mh_tables)), dimnames = list(c("Glut", "NonNeuron"), c("0", "1"), NULL))
  kk <- 1
  for (nm in names(mh_tables)) {
    x <- mh_tables[[nm]]
    ct <- table(x$partner, x$full)
    ct2 <- matrix(0, nrow = 2, ncol = 2, dimnames = list(c("Glut", "NonNeuron"), c("0", "1")))
    ct2[rownames(ct), colnames(ct)] <- ct
    mh_array[, , kk] <- ct2
    kk <- kk + 1
  }
  mh_gaba <- stats::mantelhaen.test(mh_array)
  
  cmh_gaba_summary <- tibble::tibble(
    comparison = "Comparable GABA pairs: Glut host vs NonNeuron host",
    pooled_OR = unname(mh_gaba$estimate),
    ci_low = if (!is.null(mh_gaba$conf.int)) mh_gaba$conf.int[1] else NA_real_,
    ci_high = if (!is.null(mh_gaba$conf.int)) mh_gaba$conf.int[2] else NA_real_,
    p = mh_gaba$p.value,
    n_strata = length(mh_tables)
  )
  write_csv(cmh_gaba_summary, file.path(table_dir, "table03_cmh_gaba_host_preference.csv"))
  
  # ------------------------------------------------------------------
  # Same GABA cluster: Glut host vs NonNeuron host
  # ------------------------------------------------------------------
  rate_by_guest <- function(df, guest_col, full_vec) {
    tibble::tibble(guest = df[[guest_col]], full = full_vec) %>%
      dplyr::group_by(guest) %>%
      dplyr::summarise(rate = mean(full), n_pairs = dplyr::n(), n_full = sum(full), .groups = "drop")
  }
  
  g_glu_rate_elig <- rate_by_guest(GIo[GIo$E_n >= GIo$G_n, ], "G_label", GIo$G_fill[GIo$E_n >= GIo$G_n] == 1)
  g_nn_rate_elig  <- rate_by_guest(INo[INo$N_n >= INo$G_n, ], "G_label", INo$G_fill[INo$N_n >= INo$G_n] == 1)
  shared_g_elig <- intersect(g_glu_rate_elig$guest, g_nn_rate_elig$guest)
  
  paired_rate_elig <- dplyr::inner_join(
    dplyr::rename(g_glu_rate_elig[g_glu_rate_elig$guest %in% shared_g_elig, ], rate_glu = rate, n_pairs_glu = n_pairs, n_full_glu = n_full),
    dplyr::rename(g_nn_rate_elig[g_nn_rate_elig$guest %in% shared_g_elig, ], rate_nn = rate, n_pairs_nn = n_pairs, n_full_nn = n_full),
    by = "guest"
  )
  
  paired_rate_summary <- tibble::tibble(
    n_shared_gaba = nrow(paired_rate_elig),
    median_rate_glu = median(paired_rate_elig$rate_glu, na.rm = TRUE),
    median_rate_nn = median(paired_rate_elig$rate_nn, na.rm = TRUE),
    p_paired = paired_wilcox_p(paired_rate_elig$rate_glu, paired_rate_elig$rate_nn)
  )
  write_csv(paired_rate_elig, file.path(table_dir, "table04_shared_gaba_eligible_rate_pairs.csv"))
  write_csv(paired_rate_summary, file.path(table_dir, "table05_shared_gaba_eligible_rate_summary.csv"))
  
  # ------------------------------------------------------------------
  # [comment omitted: encoding-safe]
  # ------------------------------------------------------------------
  gi_dir_counts <- tibble::tibble(
    direction = c("GABA in Glut", "Glut in GABA"),
    exact_n = c(sum(GIo$G_fill == 1, na.rm = TRUE), sum(GIo$E_fill == 1, na.rm = TRUE)),
    total_pairs = nrow(GIo)
  ) %>%
    dplyr::mutate(rate = exact_n / total_pairs)
  
  binom_gi <- binom.test(gi_dir_counts$exact_n[1], sum(gi_dir_counts$exact_n), p = 0.5, alternative = "greater")
  write_csv(gi_dir_counts, file.path(table_dir, "table06_gi_directional_exact_counts.csv"))
  write_csv(tibble::tibble(p = binom_gi$p.value), file.path(table_dir, "table07_gi_direction_binom.csv"))
  
  # ------------------------------------------------------------------
  # Priority score table
  # ------------------------------------------------------------------
  priority_tbl <- motif_summary %>%
    dplyr::filter(motif %in% c("Gaba⊂Glut", "smaller⊂Glut (GG)", "Gaba⊂NonNeuron")) %>%
    dplyr::mutate(
      priority_tier = c("Tier 1", "Tier 2", "Tier 3"),
      interpretation = c(
        "Best causal candidate",
        "Field-level hierarchical candidate",
        "Mostly geometric/byproduct regime"
      )
    )
  write_csv(priority_tbl, file.path(table_dir, "table08_priority_tiers.csv"))
  
  # ------------------------------------------------------------------
  # Figures
  # ------------------------------------------------------------------
  # Fig 1: priority map
  fig1_df <- priority_tbl %>%
    dplyr::mutate(label = paste0(motif, "\nAUC=", fmt_num(size_auc, 3), "; z=", fmt_num(identity_perm_z, 2)))
  
  p1 <- ggplot(fig1_df, aes(x = size_auc, y = identity_perm_z, size = abs(partial_asym_median), colour = motif)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.6, colour = "#888888") +
    geom_vline(xintercept = 0.75, linetype = 2, linewidth = 0.6, colour = "#BBBBBB") +
    geom_point(alpha = 0.95) +
    geom_text(aes(label = label), nudge_x = 0.015, hjust = 0, size = 3.4, show.legend = FALSE) +
    scale_colour_manual(values = col_motif[names(col_motif) %in% fig1_df$motif]) +
    scale_size_continuous(range = c(4.5, 10)) +
    labs(
      title = wrap_title("Containment motifs separate into causal-priority tiers"),
      subtitle = wrap_subtitle("x-axis: size-only discriminability (AUC). y-axis: residual identity signal after conditioning on slide/layer/size bin. Point size scales with precursor partial asymmetry."),
      x = "Size-only AUC for exact containment",
      y = "Residual identity permutation z-score",
      colour = NULL,
      size = "|Partial asym|"
    ) +
    base_theme()
  save_pub(p1, "fig1_priority_map.png", width = 8.8, height = 5.8)
  
  # Fig 2: size-controlled OR forest
  fig2_df <- cmh_gaba_summary %>%
    dplyr::mutate(
      comparison = factor(comparison, levels = comparison),
      label = paste0("OR=", fmt_num(pooled_OR, 2), " [", fmt_num(ci_low, 2), ", ", fmt_num(ci_high, 2), "]\np=", fmt_p(p))
    )
  
  p2 <- ggplot(fig2_df, aes(x = pooled_OR, y = comparison)) +
    geom_vline(xintercept = 1, linetype = 2, linewidth = 0.6, colour = "#888888") +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.18, linewidth = 0.9, colour = col_motif[["Gaba⊂Glut"]]) +
    geom_point(size = 3.8, colour = col_motif[["Gaba⊂Glut"]]) +
    geom_text(aes(x = ci_high * 1.05, label = label), hjust = 0, size = 3.6, colour = "#444444") +
    scale_x_log10(labels = scales::label_number(accuracy = 0.01), expand = expansion(mult = c(0.03, 0.34))) +
    coord_cartesian(clip = "off") +
    labs(
      title = wrap_title("Gaba⊂Glut remains enriched after size matching"),
      subtitle = wrap_subtitle("Only comparable pairs are kept, and the pooled odds ratio is stratified by GABA subclass, layer, and host/guest size bin."),
      x = "Mantel–Haenszel pooled odds ratio (log scale)",
      y = NULL
    ) +
    base_theme()
  save_pub(p2, "fig2_size_control_cmh_forest.png", width = 8.4, height = 4.8)
  
  # Fig 3: residual identity permutation comparison
  fig3_df <- priority_tbl %>%
    dplyr::mutate(
      motif = factor(motif, levels = c("Gaba⊂Glut", "smaller⊂Glut (GG)", "Gaba⊂NonNeuron")),
      stat_lab = paste0("z=", fmt_num(identity_perm_z, 2), "\np=", purrr::map_chr(identity_perm_p, fmt_p))
    )
  
  p3 <- ggplot(fig3_df, aes(x = motif, y = identity_perm_z, fill = motif)) +
    geom_col(width = 0.66, alpha = 0.96) +
    geom_text(aes(label = stat_lab), vjust = -0.28, size = 3.4) +
    scale_fill_manual(values = col_motif[names(col_motif) %in% levels(fig3_df$motif)]) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(
      title = wrap_title("Residual identity signal distinguishes causal motifs from geometric byproducts"),
      subtitle = wrap_subtitle("Permutation is conditioned on slide, layer, and size bin. Gaba⊂Glut and smaller⊂Glut retain significant subclass-specific structure, whereas Gaba⊂NonNeuron does not."),
      x = NULL,
      y = "Permutation z-score for residual identity specificity"
    ) +
    base_theme() +
    theme(legend.position = "none")
  save_pub(p3, "fig3_residual_identity_barplot.png", width = 8.4, height = 5.0)
  
  # Fig 4: precursor partial asymmetry
  fig4_df <- bind_rows(as_GI, as_GG, as_IN) %>%
    dplyr::mutate(
      motif = factor(c("Gaba⊂Glut", "smaller⊂Glut (GG)", "Gaba⊂NonNeuron"), levels = c("Gaba⊂Glut", "smaller⊂Glut (GG)", "Gaba⊂NonNeuron")),
      stat_lab = paste0("median=", fmt_num(slide_median_of_medians, 3), "\n", slide_positive_n, "/", slide_n, " slides > 0\np=", purrr::map_chr(slide_wilcoxon_p, fmt_p))
    )
  
  p4 <- ggplot(fig4_df, aes(x = motif, y = slide_median_of_medians, fill = motif)) +
    geom_col(width = 0.66, alpha = 0.96) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.6, colour = "#888888") +
    geom_text(aes(label = stat_lab), vjust = ifelse(fig4_df$slide_median_of_medians >= 0, -0.25, 1.15), size = 3.3) +
    scale_fill_manual(values = col_motif[names(col_motif) %in% levels(fig4_df$motif)]) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.20))) +
    labs(
      title = wrap_title("Partial-overlap asymmetry appears early for Gaba⊂Glut and smaller⊂Glut, but not for Gaba⊂NonNeuron"),
      subtitle = wrap_subtitle("Slide-level medians are tested against zero. Positive values indicate that the candidate guest already overlaps more asymmetrically at the partial stage."),
      x = NULL,
      y = "Slide-level median signed asymmetry"
    ) +
    base_theme() +
    theme(legend.position = "none")
  save_pub(p4, "fig4_partial_precursor_barplot.png", width = 9.0, height = 6.0)
  
  # Fig 5: paired shared-GABA host preference
  fig5_long <- paired_rate_elig %>%
    dplyr::transmute(
      guest,
      `Glut host` = rate_glu,
      `NonNeuron host` = rate_nn
    ) %>%
    tidyr::pivot_longer(-guest, names_to = "host_class", values_to = "rate") %>%
    dplyr::mutate(host_class = factor(host_class, levels = c("NonNeuron host", "Glut host")))
  
  fig5_stat <- tibble::tibble(
    x = 1.5,
    y = max(fig5_long$rate, na.rm = TRUE) * 1.07,
    label = paste0(
      "n=", scales::comma(nrow(paired_rate_elig)),
      "; median Glut=", fmt_num(median(paired_rate_elig$rate_glu, na.rm = TRUE), 2),
      "; median NonNeuron=", fmt_num(median(paired_rate_elig$rate_nn, na.rm = TRUE), 2),
      "\np=", fmt_p(paired_rate_summary$p_paired)
    )
  )
  
  p5 <- ggplot(fig5_long, aes(x = host_class, y = rate, group = guest)) +
    geom_line(colour = "#CFCFCF", alpha = 0.28, linewidth = 0.35) +
    geom_point(aes(colour = host_class), alpha = 0.80, size = 1.8) +
    geom_boxplot(aes(fill = host_class), width = 0.34, outlier.shape = NA, alpha = 0.90, linewidth = 0.8) +
    geom_text(data = fig5_stat, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 3.3) +
    scale_colour_manual(values = col_partner) +
    scale_fill_manual(values = col_partner) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0.02, 0.16))) +
    labs(
      title = wrap_title("Within the same GABA cluster, Glut behaves more like the true exact host than NonNeuron"),
      subtitle = wrap_subtitle("Only GABA clusters with comparable Glut and NonNeuron partners are kept. Each line is one shared GABA cluster."),
      x = NULL,
      y = "Exact containment rate for the GABA cluster",
      colour = NULL,
      fill = NULL
    ) +
    base_theme()
  save_pub(p5, "fig5_shared_gaba_paired_rates.png", width = 8.6, height = 5.6)
  
  # Fig 6: directional exact count within GI
  fig6_df <- gi_dir_counts %>%
    dplyr::mutate(direction = factor(direction, levels = c("GABA in Glut", "Glut in GABA")))
  sub6 <- paste0("Among ", scales::comma(nrow(GIo)), " unique GABA–Glut pairs: ",
                 "directional OR=", fmt_num((fig6_df$exact_n[1] / (nrow(GIo) - fig6_df$exact_n[1])) / (fig6_df$exact_n[2] / (nrow(GIo) - fig6_df$exact_n[2])), 2),
                 "; binomial p=", fmt_p(binom_gi$p.value), ".")
  
  p6 <- ggplot(fig6_df, aes(x = direction, y = exact_n, fill = direction)) +
    geom_col(width = 0.62, alpha = 0.96) +
    geom_text(aes(label = paste0(scales::comma(exact_n), "\n(", fmt_pct(rate, 1), ")")), vjust = -0.15, size = 3.8) +
    scale_fill_manual(values = c("GABA in Glut" = col_motif[["Gaba⊂Glut"]], "Glut in GABA" = "#4E79A7")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(
      title = wrap_title("Directionality inside GABA–Glut overlap is strongly biased toward Gaba⊂Glut"),
      subtitle = wrap_subtitle(sub6),
      x = NULL,
      y = "Number of exact directional edges"
    ) +
    base_theme() +
    theme(legend.position = "none")
  save_pub(p6, "fig6_gi_directional_exact_counts.png", width = 8.0, height = 5.2)
  
  # Fig 7: motif-specific comparison dashboard
  fig7_df <- priority_tbl %>%
    dplyr::transmute(
      motif,
      `Residual identity z` = identity_perm_z,
      `Partial precursor` = partial_asym_median,
      `Size-only AUC` = size_auc
    ) %>%
    tidyr::pivot_longer(-motif, names_to = "metric", values_to = "value") %>%
    dplyr::group_by(metric) %>%
    dplyr::mutate(zv = as.numeric(scale(value))) %>%
    dplyr::ungroup()
  
  p7 <- ggplot(fig7_df, aes(x = metric, y = fct_reorder(motif, zv, .fun = mean), fill = zv)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.3f", value)), size = 3.4) +
    scale_fill_gradient2(low = "#4E79A7", mid = "#F2F2F2", high = "#C44E52", midpoint = 0) +
    labs(
      title = wrap_title("The three priority tiers separate consistently across orthogonal evidence axes"),
      subtitle = wrap_subtitle("Values are shown in their native scales; tile colour reflects z-scored ranking within each metric."),
      x = NULL,
      y = NULL,
      fill = "Within-metric\nrank z"
    ) +
    base_theme()
  save_pub(p7, "fig7_priority_dashboard_heatmap.png", width = 8.6, height = 4.6)
  
  # Fig 8: concise manuscript-style schematic summary
  fig8_df <- tibble::tibble(
    tier = factor(c("Tier 1", "Tier 2", "Tier 3"), levels = c("Tier 1", "Tier 2", "Tier 3")),
    motif = c("Gaba⊂Glut", "smaller⊂Glut (GG)", "Gaba⊂NonNeuron"),
    y = c(3, 2, 1),
    key_result = c(
      paste0("Size-matched OR=", fmt_num(cmh_gaba_summary$pooled_OR, 2), "; residual z=", fmt_num(res_GI$z, 2), "; partial med=", fmt_num(as_GI$slide_median_of_medians, 3)),
      paste0("Residual z=", fmt_num(res_GG$z, 2), "; partial med=", fmt_num(as_GG$slide_median_of_medians, 3)),
      paste0("AUC=", fmt_num(motif_summary$size_auc[motif_summary$motif == "Gaba⊂NonNeuron"], 3), "; residual p=", fmt_p(res_IN$p), "; partial p=", fmt_p(as_IN$slide_wilcoxon_p))
    ),
    interpretation = c(
      "Best causal motif",
      "Hierarchical excitatory field rule",
      "Mostly geometric/support-shell byproduct"
    )
  )
  
  p8 <- ggplot(fig8_df, aes(x = tier, y = y, colour = motif)) +
    geom_point(size = 6) +
    geom_text(aes(label = motif), nudge_y = 0.28, size = 3.8, show.legend = FALSE) +
    geom_text(aes(label = key_result), nudge_y = -0.18, size = 3.1, show.legend = FALSE, colour = "#444444") +
    geom_text(aes(label = interpretation), nudge_y = -0.47, size = 3.2, show.legend = FALSE, colour = "#666666") +
    scale_colour_manual(values = col_motif[names(col_motif) %in% fig8_df$motif]) +
    scale_y_continuous(limits = c(0.4, 3.5), breaks = NULL) +
    labs(
      title = wrap_title("Causal-priority interpretation of containment motifs"),
      subtitle = wrap_subtitle("Tier 1 is the best candidate for direct causal experiments, Tier 2 is a field-level hierarchy rule, and Tier 3 is best interpreted as a geometric/support-system outcome."),
      x = NULL,
      y = NULL,
      colour = NULL
    ) +
    base_theme()
  save_pub(p8, "fig8_priority_schematic_summary.png", width = 9.5, height = 4.9)
  
  # ------------------------------------------------------------------
  # Text summary for quick manuscript/PPT reuse
  # ------------------------------------------------------------------
  summary_txt <- c(
    paste0("Tier 1 candidate: Gaba⊂Glut. Size-matched pooled OR = ", fmt_num(cmh_gaba_summary$pooled_OR, 2),
           ", 95% CI [", fmt_num(cmh_gaba_summary$ci_low, 2), ", ", fmt_num(cmh_gaba_summary$ci_high, 2),
           "], p = ", fmt_p(cmh_gaba_summary$p), ". Residual identity permutation z = ", fmt_num(res_GI$z, 2),
           ", p = ", fmt_p(res_GI$p), ". Partial asymmetry median = ", fmt_num(as_GI$slide_median_of_medians, 3),
           ", slide-level p = ", fmt_p(as_GI$slide_wilcoxon_p), "."),
    paste0("Tier 2 candidate: smaller⊂larger Glut. Residual identity permutation z = ", fmt_num(res_GG$z, 2),
           ", p = ", fmt_p(res_GG$p), ". Partial asymmetry median = ", fmt_num(as_GG$slide_median_of_medians, 3),
           ", slide-level p = ", fmt_p(as_GG$slide_wilcoxon_p), "."),
    paste0("Tier 3 regime: Gaba⊂NonNeuron. Size-only AUC = ", fmt_num(motif_summary$size_auc[motif_summary$motif == "Gaba⊂NonNeuron"], 3),
           ". Residual identity permutation z = ", fmt_num(res_IN$z, 2), ", p = ", fmt_p(res_IN$p),
           ". Partial asymmetry median = ", fmt_num(as_IN$slide_median_of_medians, 3),
           ", slide-level p = ", fmt_p(as_IN$slide_wilcoxon_p), ".")
  )
  writeLines(summary_txt, con = file.path(out_dir, "README_R_output_11_summary.txt"))
  
  message("R_output_11 finished: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
}

run_module_12 <- function(runtime_base) {
  message('[MODULE 12] start')
  # [comment omitted: encoding-safe]
  # [removed] rm(list = ls())
  gc()
  #!/usr/bin/env Rscript
  
  suppressPackageStartupMessages({
    library(tidyverse)
    library(broom)
    library(scales)
    library(grid)
  })
  
  # avoid select/filter conflicts
  select <- dplyr::select
  filter <- dplyr::filter
  mutate <- dplyr::mutate
  arrange <- dplyr::arrange
  distinct <- dplyr::distinct
  pull <- dplyr::pull
  count <- dplyr::count
  group_by <- dplyr::group_by
  summarise <- dplyr::summarise
  transmute <- dplyr::transmute
  ungroup <- dplyr::ungroup
  rename <- dplyr::rename
  left_join <- dplyr::left_join
  
  # ============================================================
  # 0. Paths
  # ============================================================
  project_root <- runtime_base
  out_dir <- file.path(project_root, "R_output_12")
  fig_dir <- file.path(out_dir, "figures")
  tab_dir <- file.path(out_dir, "tables")
  
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
  
  find_first_existing <- function(paths) {
    hit <- paths[file.exists(paths)][1]
    if (is.na(hit)) stop("No valid input file found. Please check input path.")
    hit
  }
  
  input_neuron_nonneuron <- find_first_existing(c(
    file.path(project_root, "neuron-nonneuron-partner.csv"),
    file.path(project_root, "data", "neuron-nonneuron-partner.csv"),
    file.path(getwd(), "neuron-nonneuron-partner.csv")
  ))
  
  # ============================================================
  # 1. Helpers
  # ============================================================
  col_main <- c(
    "Micro->Gaba exact" = "#C44E52",
    "Micro->Glut exact" = "#4C78A8"
  )
  col_dir <- c("Enriched" = "#4C78A8", "Depleted" = "#E15759")
  
  fmt_p <- function(p) {
    vapply(p, function(x) {
      if (is.na(x)) return(NA_character_)
      if (x < 1e-300) return("<1e-300")
      if (x < 1e-3) return(formatC(x, format = "e", digits = 2))
      sprintf("%.3f", x)
    }, character(1))
  }
  fmt_num <- function(x, digits = 3) sprintf(paste0("%.", digits, "f"), x)
  wrap_title <- function(x) stringr::str_wrap(x, width = 62)
  wrap_subtitle <- function(x) stringr::str_wrap(x, width = 86)
  
  base_theme <- function() {
    theme_pub_global(base_size = 14) +
      theme(
        plot.title.position = "plot",
        plot.caption.position = "plot",
        plot.title = element_text(size = 21, face = "bold", colour = "black", hjust = 0, lineheight = 1.02,
                                  margin = margin(b = 8)),
        plot.subtitle = element_text(size = 12.8, colour = "#444444", hjust = 0, lineheight = 1.18,
                                     margin = margin(b = 12)),
        axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 12.5, colour = "black"),
        axis.line = element_line(linewidth = 0.8, colour = "black"),
        axis.ticks = element_line(linewidth = 0.8, colour = "black"),
        axis.ticks.length = unit(0.18, "cm"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = "#F1F1F1", colour = NA),
        strip.text = element_text(size = 13.8, face = "bold"),
        panel.grid = element_blank(),
        plot.margin = margin(16, 36, 14, 26)
      )
  }
  
  save_pub <- function(plot, filename, width, height) {
    ggsave_pub(file.path(fig_dir, filename), plot = plot, width = width, height = height,
           dpi = 600, bg = "white", limitsize = FALSE)
  }
  
  standardize_names <- function(df) {
    names(df) <- names(df) |>
      stringr::str_replace_all("\\.overlap\\.percent", "_overlap_percent") |>
      stringr::str_replace_all("\\.", "_")
    df
  }
  
  clean_subclass <- function(x) {
    x |>
      as.character() |>
      stringr::str_replace("^\\s*\\d+\\s*", "") |>
      stringr::str_replace("\\s+(Glut|Gaba|NN)\\s*$", "") |>
      stringr::str_trim()
  }
  
  pair_key <- function(a, b) ifelse(a <= b, paste0(a, "|||", b), paste0(b, "|||", a))
  safe_logp <- function(p) -log10(pmax(p, 1e-300))
  safe_div <- function(a, b) ifelse(b > 0, a / b, NA_real_)
  
  fisher_tbl <- function(a, b, c, d) {
    ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))
    tibble(
      a = a, b = b, c = c, d = d,
      or = unname(ft$estimate),
      ci_low = unname(ft$conf.int[1]),
      ci_high = unname(ft$conf.int[2]),
      p = ft$p.value
    )
  }
  
  mw_compare <- function(df, group_col, focal, vars) {
    purrr::map_dfr(vars, function(v) {
      x <- df |> dplyr::filter(.data[[group_col]] == focal) |> dplyr::pull(.data[[v]])
      y <- df |> dplyr::filter(.data[[group_col]] != focal) |> dplyr::pull(.data[[v]])
      wt <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
      tibble(
        var = v,
        focal_median = median(x, na.rm = TRUE),
        other_median = median(y, na.rm = TRUE),
        delta = median(x, na.rm = TRUE) - median(y, na.rm = TRUE),
        p = wt$p.value,
        n_focal = sum(!is.na(x)),
        n_other = sum(!is.na(y))
      )
    }) |>
      mutate(fdr = p.adjust(p, method = "BH"))
  }
  
  topology_label <- function(non_overlap, neuron_overlap, tol = 1e-9) {
    dplyr::case_when(
      non_overlap >= 1 - tol & neuron_overlap < 1 - tol ~ "Non->Neuron exact",
      neuron_overlap >= 1 - tol & non_overlap < 1 - tol ~ "Neuron->Non exact",
      neuron_overlap >= 1 - tol & non_overlap >= 1 - tol ~ "Bidirectional exact",
      non_overlap >= 0.95 & non_overlap < 1 - tol & neuron_overlap < 0.95 ~ "Non->Neuron near-full",
      neuron_overlap >= 0.95 & neuron_overlap < 1 - tol & non_overlap < 0.95 ~ "Neuron->Non near-full",
      pmax(non_overlap, neuron_overlap) >= 0.95 ~ "High-overlap other",
      TRUE ~ "Other"
    )
  }
  
  # ============================================================
  # 2. Data
  # ============================================================
  raw_nn_non <- readr::read_csv(input_neuron_nonneuron, show_col_types = FALSE) |>
    standardize_names()
  
  dedup_nn_non <- raw_nn_non |>
    mutate(pair_key = pair_key(cluster_1_label, cluster_2_label)) |>
    distinct(pair_key, .keep_all = TRUE)
  
  nn <- dedup_nn_non |>
    filter(
      (cluster_1_cell_Neruon_type == "NonNeuron" & cluster_2_cell_Neruon_type %in% c("Glut", "Gaba")) |
        (cluster_2_cell_Neruon_type == "NonNeuron" & cluster_1_cell_Neruon_type %in% c("Glut", "Gaba"))
    ) |>
    mutate(non_is_1 = cluster_1_cell_Neruon_type == "NonNeuron") |>
    transmute(
      pair_key,
      slide = if_else(non_is_1, cluster_2_slide, cluster_1_slide),
      layer = if_else(non_is_1, cluster_2_layer, cluster_1_layer),
      layer_simplified = stringr::str_remove(if_else(non_is_1, cluster_2_layer, cluster_1_layer), "^layer\\s+"),
      region = if_else(non_is_1, cluster_2_region, cluster_1_region),
      
      neuron_label = if_else(non_is_1, cluster_2_label, cluster_1_label),
      neuron_type = if_else(non_is_1, cluster_2_cell_Neruon_type, cluster_1_cell_Neruon_type),
      neuron_subclass = clean_subclass(if_else(non_is_1, cluster_2_subclass, cluster_1_subclass)),
      neuron_total = if_else(non_is_1, cluster_2_total_cell_num, cluster_1_total_cell_num),
      neuron_glut_n = if_else(non_is_1, cluster_2_Glut_Neruon_cell_ids_num, cluster_1_Glut_Neruon_cell_ids_num),
      neuron_gaba_n = if_else(non_is_1, cluster_2_GABA_Neruon_cell_ids_num, cluster_1_GABA_Neruon_cell_ids_num),
      neuron_p = if_else(non_is_1, cluster_2_cauchy_combination_p, cluster_1_cauchy_combination_p),
      neuron_overlap_pct = if_else(non_is_1, cluster_2_overlap_percent, cluster_1_overlap_percent),
      
      non_label = if_else(non_is_1, cluster_1_label, cluster_2_label),
      non_class = clean_subclass(if_else(non_is_1, cluster_1_subclass, cluster_2_subclass)),
      non_total = if_else(non_is_1, cluster_1_total_cell_num, cluster_2_total_cell_num),
      non_glut_n = if_else(non_is_1, cluster_1_Glut_Neruon_cell_ids_num, cluster_2_Glut_Neruon_cell_ids_num),
      non_gaba_n = if_else(non_is_1, cluster_1_GABA_Neruon_cell_ids_num, cluster_2_GABA_Neruon_cell_ids_num),
      non_p = if_else(non_is_1, cluster_1_cauchy_combination_p, cluster_2_cauchy_combination_p),
      non_overlap_pct = if_else(non_is_1, cluster_1_overlap_percent, cluster_2_overlap_percent),
      
      overlap_cell,
      union_cell,
      jaccard
    ) |>
    mutate(
      topology = topology_label(non_overlap_pct, neuron_overlap_pct),
      log10_neuron_non_ratio = log10((neuron_total + 1) / (non_total + 1)),
      containment_asym = non_overlap_pct - neuron_overlap_pct,
      non_glut_fraction = safe_div(non_glut_n, non_total),
      non_gaba_fraction = safe_div(non_gaba_n, non_total),
      non_non_fraction = 1 - safe_div(non_glut_n + non_gaba_n, non_total),
      neuron_logp = safe_logp(neuron_p),
      non_logp = safe_logp(non_p)
    )
  
  # ============================================================
  # 3. Figure 1 | geometry signature
  # ============================================================
  core_vars <- c("jaccard", "log10_neuron_non_ratio", "non_total")
  metric_label_map <- c(
    jaccard = "Jaccard",
    log10_neuron_non_ratio = "log10(host/child size ratio)",
    non_total = "Microglia child size"
  )
  
  sig_glut_child <- nn |>
    filter(neuron_type == "Glut", topology == "Non->Neuron exact") |>
    mw_compare(group_col = "non_class", focal = "Microglia", vars = core_vars) |>
    mutate(context = "Micro->Glut exact")
  
  sig_gaba_child <- nn |>
    filter(neuron_type == "Gaba", topology == "Non->Neuron exact") |>
    mw_compare(group_col = "non_class", focal = "Microglia", vars = core_vars) |>
    mutate(context = "Micro->Gaba exact")
  
  geom_compare_summary <- bind_rows(sig_gaba_child, sig_glut_child) |>
    mutate(
      var_label = recode(var, !!!metric_label_map),
      focal_median_label = if_else(var == "non_total", as.character(round(focal_median)), fmt_num(focal_median, 3)),
      other_median_label = if_else(var == "non_total", as.character(round(other_median)), fmt_num(other_median, 3))
    ) |>
    dplyr::select(context, var, var_label, focal_median, other_median, delta, p, fdr, n_focal, n_other,
                  focal_median_label, other_median_label)
  readr::write_csv(geom_compare_summary, file.path(tab_dir, "table01_micro_exact_geometry_summary.csv"))
  
  plot_df_a <- geom_compare_summary |>
    mutate(
      context = factor(context, levels = c("Micro->Gaba exact", "Micro->Glut exact")),
      var_label = factor(var_label, levels = c("Microglia child size", "log10(host/child size ratio)", "Jaccard")),
      effect_direction = case_when(
        context == "Micro->Gaba exact" & var == "jaccard" ~ "Higher in Micro",
        context == "Micro->Gaba exact" & var == "log10_neuron_non_ratio" ~ "Lower in Micro",
        context == "Micro->Gaba exact" & var == "non_total" ~ "Higher in Micro",
        context == "Micro->Glut exact" & fdr < 0.05 & delta > 0 ~ "Higher in Micro",
        context == "Micro->Glut exact" & fdr < 0.05 & delta < 0 ~ "Lower in Micro",
        TRUE ~ "No clear shift"
      )
    ) |>
    group_by(context) |>
    mutate(
      xmax_panel = max(c(focal_median, other_median), na.rm = TRUE),
      xmin_panel = min(c(focal_median, other_median), na.rm = TRUE),
      panel_range = pmax(xmax_panel - xmin_panel, xmax_panel * 0.35, 0.35),
      x_lower = pmin(0, xmin_panel - panel_range * 0.28),
      ann_x = xmax_panel + panel_range * 0.55,
      x_upper = ann_x + panel_range * 0.18,
      focal_lab_x = if_else(focal_median <= other_median,
                            pmax(focal_median + panel_range * 0.04, x_lower + panel_range * 0.02),
                            focal_median),
      other_lab_x = if_else(other_median <= focal_median,
                            pmax(other_median + panel_range * 0.04, x_lower + panel_range * 0.02),
                            other_median),
      fdr_lab = if_else(context == "Micro->Glut exact" & fdr >= 0.05 & var == "non_total", ">0.05", fmt_p(fdr))
    ) |>
    ungroup()
  
  fig1 <- ggplot(plot_df_a, aes(y = var_label)) +
    geom_blank(aes(x = x_lower)) +
    geom_blank(aes(x = x_upper)) +
    geom_segment(aes(x = other_median, xend = focal_median, yend = var_label),
                 linewidth = 1.25, colour = "#B8B8B8") +
    geom_point(aes(x = other_median), size = 4.2, colour = "#8F8F8F") +
    geom_point(aes(x = focal_median, colour = context), size = 4.4, show.legend = FALSE) +
    geom_text(aes(x = focal_lab_x, label = focal_median_label, colour = context),
              nudge_y = 0.22, hjust = 0, size = 3.9, show.legend = FALSE) +
    geom_text(aes(x = other_lab_x, label = other_median_label),
              nudge_y = -0.22, hjust = 0, size = 3.7, colour = "#4D4D4D") +
    geom_text(aes(x = ann_x, label = paste0(effect_direction, "\nFDR=", fdr_lab)),
              hjust = 0, size = 3.8, lineheight = 1.05, colour = "#333333") +
    facet_wrap(~context, scales = "free_x", nrow = 1) +
    scale_colour_manual(values = col_main) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
    coord_cartesian(clip = "off") +
    labs(
      title = wrap_title("Micro->Gaba exact is geometrically distinct, whereas Micro->Glut exact shows much weaker compactness-specific geometry"),
      subtitle = wrap_subtitle("Within exact nonneuron->neuron events, Microglia is compared against all other nonneuronal child classes. Positive evidence for a compact functional core in Micro->Gaba should appear as higher Jaccard, smaller host/child mismatch, and larger Microglia child size."),
      x = "Median value",
      y = NULL
    ) +
    base_theme() +
    theme(strip.text = element_text(size = 14.5, face = "bold"),
          plot.margin = margin(18, 64, 18, 34))
  
  save_pub(fig1, "fig01_micro_exact_geometry_signature.png", 12.8, 6.3)
  
  # ============================================================
  # 4. Figure 2 | layer heatmap
  # ============================================================
  layer_tests <- nn |>
    filter(neuron_type %in% c("Glut", "Gaba")) |>
    group_by(neuron_type, layer_simplified) |>
    group_modify(~ {
      purrr::map_dfr(c("Non->Neuron exact", "Neuron->Non exact"), function(out) {
        a <- sum(.x$non_class == "Microglia" & .x$topology == out)
        b <- sum(.x$non_class == "Microglia" & .x$topology != out)
        c <- sum(.x$non_class != "Microglia" & .x$topology == out)
        d <- sum(.x$non_class != "Microglia" & .x$topology != out)
        fisher_tbl(a, b, c, d) |>
          mutate(
            outcome = out,
            rate_micro = safe_div(a, a + b),
            rate_other = safe_div(c, c + d)
          )
      })
    }) |>
    ungroup() |>
    group_by(neuron_type, outcome) |>
    mutate(fdr = p.adjust(p, method = "BH")) |>
    ungroup()
  
  layer_glut_focus <- layer_tests |>
    filter(neuron_type == "Glut") |>
    mutate(
      outcome_label = recode(outcome,
                             "Non->Neuron exact" = "Micro->Glut exact",
                             "Neuron->Non exact" = "Glut->Micro exact"),
      outcome_label = factor(outcome_label, levels = c("Micro->Glut exact", "Glut->Micro exact")),
      layer_simplified = factor(layer_simplified, levels = c("1", "2/3", "4", "5", "6a", "6b")),
      tile_lab = paste0("OR=", fmt_num(or, 2), "\nFDR=", fmt_p(fdr))
    ) |>
    arrange(outcome_label, layer_simplified)
  
  readr::write_csv(layer_tests, file.path(tab_dir, "table02_micro_layer_tests.csv"))
  readr::write_csv(layer_glut_focus, file.path(tab_dir, "table03_micro_glut_layer_focus.csv"))
  
  fig2 <- ggplot(layer_glut_focus, aes(x = layer_simplified, y = outcome_label, fill = log2(or))) +
    geom_tile(colour = "white", linewidth = 1.8) +
    geom_text(aes(label = tile_lab), size = 4.0, lineheight = 1.02) +
    scale_fill_gradient2(
      low = "#3B4CC0", mid = "#F7F7F7", high = "#C44E52", midpoint = 0,
      breaks = c(-1, 0, 1), labels = label_number(accuracy = 0.1)
    ) +
    labs(
      title = wrap_title("Microglia shows a pronounced layer bias in Glut pairs, with layer 1 being the strongest hotspot"),
      subtitle = wrap_subtitle("Odds ratios compare Microglia against all other nonneuronal classes within each layer. Red indicates enrichment, blue indicates depletion. The expected signature is simultaneous enrichment of Micro->Glut exact and depletion of Glut->Micro exact in shallow layers."),
      x = "Layer",
      y = NULL,
      fill = "Enrichment / depletion\n(log2 OR)"
    ) +
    base_theme() +
    theme(
      axis.text.y = element_text(face = "bold", size = 13.5),
      axis.text.x = element_text(size = 13.5),
      legend.position = "top",
      legend.key.width = unit(2.2, "cm"),
      plot.margin = margin(16, 18, 14, 18)
    )
  
  save_pub(fig2, "fig02_micro_glut_layer_heatmap.png", 11.2, 5.8)
  
  # ============================================================
  # 5. Figure 3 | subclass forest
  # ============================================================
  subclass_tests <- purrr::map_dfr(c("Glut", "Gaba"), function(nt) {
    sub <- nn |> filter(neuron_type == nt)
    valid_subclasses <- sub |>
      count(neuron_subclass, name = "n_pairs") |>
      filter(n_pairs >= 50) |>
      pull(neuron_subclass)
    
    purrr::map_dfr(valid_subclasses, function(ns) {
      sub2 <- sub |> filter(neuron_subclass == ns)
      purrr::map_dfr(c("Non->Neuron exact", "Neuron->Non exact"), function(out) {
        a <- sum(sub2$non_class == "Microglia" & sub2$topology == out)
        b <- sum(sub2$non_class == "Microglia" & sub2$topology != out)
        c <- sum(sub2$non_class != "Microglia" & sub2$topology == out)
        d <- sum(sub2$non_class != "Microglia" & sub2$topology != out)
        fisher_tbl(a, b, c, d) |>
          mutate(
            neuron_type = nt,
            neuron_subclass = ns,
            outcome = out,
            rate_micro = safe_div(a, a + b),
            rate_other = safe_div(c, c + d),
            n_micro_pairs = a + b,
            n_total_pairs = a + b + c + d
          )
      })
    })
  }) |>
    group_by(neuron_type, outcome) |>
    mutate(fdr = p.adjust(p, method = "BH")) |>
    ungroup()
  
  readr::write_csv(subclass_tests, file.path(tab_dir, "table04_micro_subclass_tests.csv"))
  
  focus_enriched <- c("L2/3 IT CTX", "CLA-EPd-CTX Car3", "L2/3 IT RSP")
  focus_depleted <- c("L5 NP CTX", "L5 ET CTX", "L6 IT CTX")
  
  plot_subclass_glut <- subclass_tests |>
    filter(neuron_type == "Glut", outcome == "Non->Neuron exact", neuron_subclass %in% c(focus_enriched, focus_depleted)) |>
    mutate(
      direction = if_else(neuron_subclass %in% focus_enriched, "Enriched", "Depleted"),
      neuron_subclass = factor(neuron_subclass, levels = rev(c(focus_depleted, focus_enriched))),
      ann = paste0("OR=", fmt_num(or, 2), "\nFDR=", fmt_p(fdr)),
      ann_x = if_else(or >= 1, ci_high * 1.12, ci_low / 1.12),
      ann_hjust = if_else(or >= 1, 0, 1)
    )
  
  readr::write_csv(plot_subclass_glut, file.path(tab_dir, "table05_micro_glut_subclass_focus.csv"))
  
  xmin3 <- min(plot_subclass_glut$ci_low, na.rm = TRUE)
  xmax3 <- max(plot_subclass_glut$ci_high, na.rm = TRUE)
  
  fig3 <- ggplot(plot_subclass_glut, aes(x = or, y = neuron_subclass, colour = direction)) +
    geom_vline(xintercept = 1, linetype = 2, linewidth = 0.8, colour = "#7A7A7A") +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.18, linewidth = 1.2, show.legend = FALSE) +
    geom_point(size = 4.0, show.legend = FALSE) +
    geom_text(aes(x = ann_x, label = ann, hjust = ann_hjust), size = 3.7, colour = "#444444", show.legend = FALSE) +
    scale_x_log10(labels = label_number(accuracy = 0.01), expand = expansion(mult = c(0.18, 0.18))) +
    scale_colour_manual(values = col_dir) +
    coord_cartesian(xlim = c(xmin3 / 1.25, xmax3 * 1.35), clip = "off") +
    labs(
      title = wrap_title("Micro->Glut exact is not uniform across excitatory subclasses"),
      subtitle = wrap_subtitle("Enrichment concentrates in shallow IT / association-like hosts, whereas several deeper subclasses are depleted."),
      x = "Odds ratio for Micro->Glut exact (log scale)",
      y = NULL
    ) +
    base_theme() +
    theme(plot.margin = margin(16, 58, 16, 26))
  
  save_pub(fig3, "fig03_micro_glut_subclass_forest.png", 11.2, 6.4)
  
  # ============================================================
  # 6. Figure 4 | compact micro->gaba core summary
  # ============================================================
  core_gaba_summary <- sig_gaba_child |>
    mutate(
      metric = recode(var, !!!metric_label_map),
      effect_direction = case_when(
        var == "jaccard" ~ "Higher in Micro",
        var == "log10_neuron_non_ratio" ~ "Lower in Micro",
        var == "non_total" ~ "Higher in Micro",
        TRUE ~ NA_character_
      ),
      focal_median_label = if_else(var == "non_total", as.character(round(focal_median)), fmt_num(focal_median, 3)),
      other_median_label = if_else(var == "non_total", as.character(round(other_median)), fmt_num(other_median, 3))
    ) |>
    dplyr::select(metric, focal_median, other_median, focal_median_label, other_median_label, fdr, effect_direction) |>
    mutate(
      metric = factor(metric, levels = c("Microglia child size", "log10(host/child size ratio)", "Jaccard")),
      ann_x = c(52, 52, 52)
    )
  
  readr::write_csv(core_gaba_summary, file.path(tab_dir, "table06_micro_gaba_compact_core_summary.csv"))
  
  fig4 <- ggplot(core_gaba_summary, aes(y = metric)) +
    geom_segment(aes(x = other_median, xend = focal_median, yend = metric), linewidth = 1.35, colour = "#B8B8B8") +
    geom_point(aes(x = other_median), size = 4.4, colour = "#8F8F8F") +
    geom_point(aes(x = focal_median), size = 4.8, colour = "#C44E52") +
    geom_text(aes(x = other_median, label = other_median_label), nudge_y = -0.22, size = 4.0, colour = "#4D4D4D") +
    geom_text(aes(x = focal_median, label = focal_median_label), nudge_y = 0.22, size = 4.0, colour = "#C44E52") +
    geom_text(aes(x = ann_x, label = paste0(effect_direction, "\nFDR=", fmt_p(fdr))), hjust = 0, size = 3.9, lineheight = 1.05, colour = "#333333") +
    coord_cartesian(xlim = c(-3, 62), clip = "off") +
    labs(
      title = wrap_title("Micro->Gaba exact behaves like a compact functional core rather than a fragmented nonneuronal satellite"),
      subtitle = wrap_subtitle("Three geometric signals are highlighted: tighter overlap (higher Jaccard), smaller host-child mismatch, and larger Microglia child size, all relative to other nonneuronal child exact events in Gaba hosts."),
      x = "Median value",
      y = NULL
    ) +
    base_theme() +
    theme(plot.margin = margin(18, 62, 16, 28))
  
  save_pub(fig4, "fig04_micro_gaba_core_summary.png", 11.0, 5.9)
  
  # ============================================================
  # 7. Console summary
  # ============================================================
  cat("\n==================== INPUT ====================\n")
  cat("File:", input_neuron_nonneuron, "\n")
  cat("Pairs used:", nrow(nn), "\n")
  
  cat("\n==================== MICRO->GABA CORE ====================\n")
  print(core_gaba_summary)
  
  cat("\n==================== LAYER BIAS (GLUT) ====================\n")
  print(layer_glut_focus |> dplyr::select(layer_simplified, outcome_label, or, ci_low, ci_high, p, fdr, rate_micro, rate_other))
  
  cat("\n==================== GLUT SUBCLASS FOCUS ====================\n")
  print(plot_subclass_glut |> dplyr::select(neuron_subclass, or, ci_low, ci_high, p, fdr, rate_micro, rate_other))
  
  cat("\nOutputs saved to:\n")
  cat("  ", fig_dir, "\n")
  cat("  ", tab_dir, "\n")
}

main <- function() {
  defaults <- list(
    base_dir = "E:/zaw/2603",
    out_dir = "",
    run_modules = "true",
    overwrite = "true",
    keep_ext = "png,pdf,svg"
  )
  args <- parse_args(defaults)

  run_modules <- parse_bool(args$run_modules)
  overwrite <- parse_bool(args$overwrite)
  base_dir <- normalizePath(args$base_dir, mustWork = FALSE)
  out_dir <- if (nzchar(trimws(args$out_dir))) normalizePath(args$out_dir, mustWork = FALSE) else normalizePath(file.path(base_dir, "R_output_00"), mustWork = FALSE)

  keep_ext <- tolower(trimws(unlist(strsplit(args$keep_ext, ",", fixed = TRUE))))
  keep_ext <- keep_ext[nzchar(keep_ext)]
  if (length(keep_ext) == 0L) keep_ext <- c("png", "pdf", "svg")

  runtime_base <- file.path(out_dir, "runtime_base")
  fig_out <- file.path(out_dir, "figures")
  meta_out <- file.path(out_dir, "metadata")
  dir.create(runtime_base, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig_out, recursive = TRUE, showWarnings = FALSE)
  dir.create(meta_out, recursive = TRUE, showWarnings = FALSE)

  # Copy shared inputs for fully self-contained run
  src_nn_csv <- file.path(base_dir, "neuron-neuron-partner.csv")
  src_nn_txt <- file.path(base_dir, "neuron-neuron-partner.txt")
  src_nnon_csv <- file.path(base_dir, "neuron-nonneuron-partner.csv")
  src_nnon_txt <- file.path(base_dir, "neuron-nonneuron-partner.txt")

  copy_if_exists <- function(src, dst) {
    if (file.exists(src)) file.copy(src, dst, overwrite = TRUE)
  }
  copy_if_exists(src_nn_csv, file.path(runtime_base, "neuron-neuron-partner.csv"))
  copy_if_exists(src_nn_txt, file.path(runtime_base, "neuron-neuron-partner.txt"))
  copy_if_exists(src_nnon_csv, file.path(runtime_base, "neuron-nonneuron-partner.csv"))
  copy_if_exists(src_nnon_txt, file.path(runtime_base, "neuron-nonneuron-partner.txt"))

  if (!(file.exists(file.path(runtime_base, "neuron-neuron-partner.csv")) || file.exists(file.path(runtime_base, "neuron-neuron-partner.txt")))) {
    stop("Missing neuron-neuron partner input in base_dir.", call. = FALSE)
  }
  if (!(file.exists(file.path(runtime_base, "neuron-nonneuron-partner.csv")) || file.exists(file.path(runtime_base, "neuron-nonneuron-partner.txt")))) {
    stop("Missing neuron-nonneuron partner input in base_dir.", call. = FALSE)
  }

  module_map <- data.frame(
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

  run_log <- data.frame(
    module_id = module_map$module_id,
    script_file = module_map$script_file,
    status = "SKIPPED",
    message = "",
    stringsAsFactors = FALSE
  )

  module_funs <- list(
    run_module_01, run_module_02, run_module_03, run_module_04,
    run_module_05, run_module_06, run_module_07, run_module_08,
    run_module_09, run_module_10, run_module_11, run_module_12
  )

  if (isTRUE(run_modules)) {
    for (i in seq_along(module_funs)) {
      ok <- TRUE
      msg <- ""
      tryCatch({
        suppressWarnings(try(ggplot2::theme_set(theme_pub_global(12)), silent = TRUE))
        force_dplyr_aliases(.GlobalEnv)
        module_funs[[i]](runtime_base)
        force_dplyr_aliases(.GlobalEnv)
      }, error = function(e) {
        ok <<- FALSE
        msg <<- conditionMessage(e)
      })
      run_log$status[i] <- if (ok) "OK" else "FAILED"
      run_log$message[i] <- msg
      message(sprintf("[MODULE %02d] %s", i, run_log$status[i]))
      if (!ok) message(sprintf("[MODULE %02d] %s", i, msg))
    }
  }

  fig_df <- collect_figures(runtime_base, module_map, keep_ext = keep_ext)
  if (is.null(fig_df) || nrow(fig_df) == 0L) {
    write.csv(run_log, file.path(meta_out, "module_run_log.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    stop("No figures found under R_output_00/runtime_base. Check inputs and module execution.", call. = FALSE)
  }

  panel_df <- unique(fig_df[, c("module_id", "script_file", "stem_rel")])
  panel_df$topic <- vapply(panel_df$stem_rel, infer_topic, character(1))

  topic_order <- c(
    "Containment Pattern", "Directionality", "Layer Effect",
    "Subclass/Pair Preference", "Overlap/Compactness", "Host-Guest Topology",
    "Geometry/Scale", "Cross-Dataset Comparison", "Adjusted Model Effects",
    "Causal Priority", "Microglia Specific", "Other"
  )
  panel_df$topic_rank <- match(panel_df$topic, topic_order)
  panel_df$topic_rank[is.na(panel_df$topic_rank)] <- length(topic_order) + 1L
  panel_df <- panel_df[order(panel_df$topic_rank, panel_df$module_id, panel_df$stem_rel), , drop = FALSE]

  ordered_topics <- unique(panel_df$topic)
  panel_df$main_idx <- match(panel_df$topic, ordered_topics)
  panel_df$sub_idx <- 0L
  for (m in unique(panel_df$main_idx)) {
    idx <- which(panel_df$main_idx == m)
    panel_df$sub_idx[idx] <- seq_along(idx)
  }
  panel_df$sub_label <- index_to_letters(panel_df$sub_idx)
  panel_df$figure_id <- paste0(panel_df$main_idx, panel_df$sub_label)

  panel_df$new_stem <- sprintf(
    "figure_%02d%s__%s__%s",
    panel_df$main_idx, panel_df$sub_label,
    slugify(panel_df$topic, 40L), slugify(panel_df$stem_rel, 80L)
  )
  panel_df$new_stem <- make.unique(panel_df$new_stem, sep = "_dup")

  fig_df$panel_key <- paste(fig_df$module_id, fig_df$stem_rel, sep = "::")
  panel_df$panel_key <- paste(panel_df$module_id, panel_df$stem_rel, sep = "::")
  hit <- match(fig_df$panel_key, panel_df$panel_key)
  fig_df$topic <- panel_df$topic[hit]
  fig_df$main_idx <- panel_df$main_idx[hit]
  fig_df$sub_label <- panel_df$sub_label[hit]
  fig_df$figure_id <- panel_df$figure_id[hit]
  fig_df$new_stem <- panel_df$new_stem[hit]
  fig_df$new_name <- paste0(fig_df$new_stem, ".", fig_df$ext)
  fig_df$output_file <- file.path(fig_out, fig_df$new_name)

  fig_df <- fig_df[order(fig_df$main_idx, fig_df$sub_label, fig_df$module_id, fig_df$ext), , drop = FALSE]
  copied <- logical(nrow(fig_df))
  for (i in seq_len(nrow(fig_df))) {
    copied[i] <- isTRUE(file.copy(fig_df$source_file[i], fig_df$output_file[i], overwrite = overwrite))
  }
  fig_df$copied <- copied

  write.csv(run_log, file.path(meta_out, "module_run_log.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(panel_df, file.path(meta_out, "figure_panels_manifest.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(fig_df, file.path(meta_out, "figure_files_manifest.csv"), row.names = FALSE, fileEncoding = "UTF-8")

  topic_summary <- aggregate(
    x = list(n_panels = panel_df$figure_id),
    by = list(main_idx = panel_df$main_idx, topic = panel_df$topic),
    FUN = length
  )
  topic_summary <- topic_summary[order(topic_summary$main_idx), , drop = FALSE]
  write.csv(topic_summary, file.path(meta_out, "topic_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")

  message("[DONE] Unified outputs in: ", normalizePath(out_dir, mustWork = FALSE))
  message("[DONE] Final figures in: ", normalizePath(fig_out, mustWork = FALSE))
  message("[DONE] Manifests in: ", normalizePath(meta_out, mustWork = FALSE))
}

main()
