############################################################
# gaba_glut_focal_exact_source_models.R
#
# Purpose
#   Main goal: model which factors influence exact containment
#   between neuron-neuron cluster pairs.
#
#   The primary analysis creates directed/oriented rows from each
#   neuron-neuron pair:
#     source cluster -> target cluster
#   and models whether the source cluster is exactly contained in the
#   target cluster.
#
#   Main models:
#     1) Total-effect model:
#        exact containment ~ source type * target type + source size +
#        target size + layer
#     2) Composition-adjusted model:
#        the same model plus source focal-cell fraction and target GABA
#        fraction.
#     3) Main-model stratified contrasts:
#        derive source-specific target effects and target-specific source
#        effects from the fitted total-effect model; these are forest
#        plots of model contrasts, not descriptive bar plots.
#
#   Supplementary focal models:
#     1) Gaba-focal model: whether a Gaba source is exactly contained
#        in its partner cluster.
#     2) Glut-focal model: whether a Glut source is exactly contained
#        in its partner cluster.
#
#   The model follows the exact-containment formula family in
#   neuron_ie_special_cross_dataset_analysis_full.R, but removes
#   class_1_region_and_class_2_region, and replaces raw size/E-I terms
#   with log source/target size, source focal purity, target GABA fraction,
#   and layer adjustment.
#   target_class3 is the type of the oriented target/partner cluster
#   after the focal Gaba or Glut cluster has been assigned to source.
#   For focal_position == cluster1 this is cluster2_cell_type; for
#   focal_position == cluster2 this is cluster1_cell_type.
#
# Important correction
#   The previous draft mixed neuron-neuron and neuron-nonneuron rows in
#   one model. NonNeuron targets come from a different pair table, so
#   target_class3NonNeuron was confounded with dataset/source-generation
#   logic. This version fits the primary interpretable model only inside
#   neuron-neuron pairs, comparing Gaba targets against Glut targets.
############################################################

options(stringsAsFactors = FALSE, scipen = 999)

# -----------------------------
# Package setup
# -----------------------------
required_pkgs <- c(
  "dplyr", "readr", "tibble", "sandwich", "lmtest", "ggplot2"
)
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing R packages: ", paste(missing_pkgs, collapse = ", "),
    ". Please install them before running this analysis.",
    call. = FALSE
  )
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

# -----------------------------
# User paths
# -----------------------------
project_root <- Sys.getenv("FOCAL_MODEL_PROJECT_ROOT", "E:/zaw/2603")
input_nn <- file.path(project_root, "neuron-neuron-partner.csv")
input_nonn <- file.path(project_root, "neuron-nonneuron-partner.csv")

output_dir <- file.path(project_root, "gaba_glut_focal_exact_source_model_outputs")
table_dir <- file.path(output_dir, "tables")
fig_dir <- file.path(output_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

make_top_or_diagnostic <- identical(Sys.getenv("FOCAL_MODEL_MAKE_TOP_OR", ""), "1")

# -----------------------------
# Helper functions
# -----------------------------
rename_partner_columns <- function(df) {
  name_map <- c(
    "cluster.1_label" = "cluster1_label",
    "cluster.1_slide" = "cluster1_slide",
    "cluster.1_layer" = "cluster1_layer",
    "cluster.1_region" = "cluster1_region",
    "cluster.1_cell_Neruon_type" = "cluster1_cell_type",
    "cluster.1_subclass" = "cluster1_subclass",
    "cluster.1_total_cell_num" = "cluster1_total_cell_num",
    "cluster.1_Glut_Neruon_cell_ids_num" = "cluster1_glut_num",
    "cluster.1_GABA_Neruon_cell_ids_num" = "cluster1_gaba_num",
    "cluster.1_cauchy_combination_p" = "cluster1_cauchy_p",
    "cluster.1_E_I_Ratio" = "cluster1_ei_ratio",
    "cluster.2_label" = "cluster2_label",
    "cluster.2_slide" = "cluster2_slide",
    "cluster.2_layer" = "cluster2_layer",
    "cluster.2_region" = "cluster2_region",
    "cluster.2_cell_Neruon_type" = "cluster2_cell_type",
    "cluster.2_subclass" = "cluster2_subclass",
    "cluster.2_total_cell_num" = "cluster2_total_cell_num",
    "cluster.2_Glut_Neruon_cell_ids_num" = "cluster2_glut_num",
    "cluster.2_GABA_Neruon_cell_ids_num" = "cluster2_gaba_num",
    "cluster.2_cauchy_combination_p" = "cluster2_cauchy_p",
    "cluster.2_E_I_Ratio" = "cluster2_ei_ratio",
    "overlap_cell" = "overlap_cell",
    "union_cell" = "union_cell",
    "jaccard" = "jaccard",
    "cluster.1.overlap.percent" = "cluster1_overlap",
    "cluster.2.overlap.percent" = "cluster2_overlap"
  )
  old_names <- names(name_map)
  present <- intersect(old_names, names(df))
  names(df)[match(present, names(df))] <- unname(name_map[present])
  required <- unname(name_map)
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("Input file missing required columns after name mapping: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  df <- df[, required, drop = FALSE]
  df
}

safe_fraction <- function(numerator, denominator) {
  ifelse(denominator > 0, numerator / denominator, NA_real_)
}

target_class_factor <- function(x) {
  factor(
    dplyr::case_when(
      x == "Gaba" ~ "Gaba",
      x == "Glut" ~ "Glut",
      x == "NonNeuron" ~ "NonNeuron",
      TRUE ~ NA_character_
    ),
    levels = c("Gaba", "Glut", "NonNeuron")
  )
}

clean_layer <- function(x) {
  y <- tolower(trimws(as.character(x)))
  y <- gsub("^layer\\s*", "", y)
  y <- gsub("\\s+", "", y)
  y <- gsub("^2[._-]?3$", "2/3", y)
  factor(y, levels = c("1", "2/3", "4", "5", "6a", "6b"))
}

make_focal_rows_from_position <- function(df, focal_type, focal_position) {
  if (focal_position == 1) {
    out <- df %>%
      filter(cluster1_cell_type == focal_type) %>%
      transmute(
        dataset = .data$dataset,
        focal_type = focal_type,
        focal_position = "cluster1",
        source_label = cluster1_label,
        target_label = cluster2_label,
        source_type = cluster1_cell_type,
        target_type = cluster2_cell_type,
        slide_layer = paste(cluster1_slide, cluster1_layer, sep = "__"),
        layer_f = clean_layer(cluster1_layer),
        source_overlap = as.numeric(cluster1_overlap),
        target_overlap = as.numeric(cluster2_overlap),
        source_size = as.numeric(cluster1_total_cell_num),
        target_size = as.numeric(cluster2_total_cell_num),
        log_source_size = log1p(as.numeric(cluster1_total_cell_num)),
        log_target_size = log1p(as.numeric(cluster2_total_cell_num)),
        source_glut_fraction = safe_fraction(cluster1_glut_num, cluster1_total_cell_num),
        target_glut_fraction = safe_fraction(cluster2_glut_num, cluster2_total_cell_num),
        source_gaba_fraction = safe_fraction(cluster1_gaba_num, cluster1_total_cell_num),
        target_gaba_fraction = safe_fraction(cluster2_gaba_num, cluster2_total_cell_num),
        cluster_1_E_I_Ratio = as.numeric(cluster1_ei_ratio),
        cluster_2_E_I_Ratio = as.numeric(cluster2_ei_ratio),
        exact_source_in_target = as.integer(as.numeric(cluster1_overlap) >= 1),
        target_class3 = target_class_factor(target_type)
      )
  } else if (focal_position == 2) {
    out <- df %>%
      filter(cluster2_cell_type == focal_type) %>%
      transmute(
        dataset = .data$dataset,
        focal_type = focal_type,
        focal_position = "cluster2",
        source_label = cluster2_label,
        target_label = cluster1_label,
        source_type = cluster2_cell_type,
        target_type = cluster1_cell_type,
        slide_layer = paste(cluster2_slide, cluster2_layer, sep = "__"),
        layer_f = clean_layer(cluster2_layer),
        source_overlap = as.numeric(cluster2_overlap),
        target_overlap = as.numeric(cluster1_overlap),
        source_size = as.numeric(cluster2_total_cell_num),
        target_size = as.numeric(cluster1_total_cell_num),
        log_source_size = log1p(as.numeric(cluster2_total_cell_num)),
        log_target_size = log1p(as.numeric(cluster1_total_cell_num)),
        source_glut_fraction = safe_fraction(cluster2_glut_num, cluster2_total_cell_num),
        target_glut_fraction = safe_fraction(cluster1_glut_num, cluster1_total_cell_num),
        source_gaba_fraction = safe_fraction(cluster2_gaba_num, cluster2_total_cell_num),
        target_gaba_fraction = safe_fraction(cluster1_gaba_num, cluster1_total_cell_num),
        cluster_1_E_I_Ratio = as.numeric(cluster2_ei_ratio),
        cluster_2_E_I_Ratio = as.numeric(cluster1_ei_ratio),
        exact_source_in_target = as.integer(as.numeric(cluster2_overlap) >= 1),
        target_class3 = target_class_factor(target_type)
      )
  } else {
    stop("focal_position must be 1 or 2")
  }

  out
}

make_focal_model_input <- function(pair_df, focal_type) {
  bind_rows(
    make_focal_rows_from_position(pair_df, focal_type, 1),
    make_focal_rows_from_position(pair_df, focal_type, 2)
  ) %>%
    distinct(dataset, focal_type, source_label, target_label, .keep_all = TRUE) %>%
    filter(
      dataset == "neuron_neuron",
      target_class3 %in% c("Gaba", "Glut")
    ) %>%
    mutate(
      source_focal_fraction = ifelse(
        focal_type == "Gaba",
        source_gaba_fraction,
        source_glut_fraction
      ),
      target_class3 = factor(as.character(target_class3), levels = c("Gaba", "Glut"))
    ) %>%
    filter(
      !is.na(target_class3),
      !is.na(layer_f),
      is.finite(source_focal_fraction), is.finite(target_gaba_fraction),
      is.finite(log_source_size), is.finite(log_target_size)
    ) %>%
    mutate(target_class3 = factor(target_class3, levels = c("Gaba", "Glut")))
}

write_model_outputs <- function(model, vcov_mat, output_prefix) {
  coef_test <- lmtest::coeftest(model, vcov. = vcov_mat)
  coef_df <- tibble(
    term = rownames(coef_test),
    estimate = coef_test[, 1],
    se = coef_test[, 2],
    z = coef_test[, 3],
    p = coef_test[, 4]
  ) %>%
    mutate(
      OR = exp(estimate),
      ci_low = exp(estimate - 1.96 * se),
      ci_high = exp(estimate + 1.96 * se)
    )

  readr::write_csv(coef_df, file.path(table_dir, paste0(output_prefix, "_coefficients.csv")))
  coef_df
}


format_p_label <- function(p) {
  ifelse(
    is.na(p),
    "p=NA",
    ifelse(p < 0.001, paste0("p=", formatC(p, format = "e", digits = 3)),
           paste0("p=", signif(p, 3)))
  )
}

wrap_label <- function(x, width = 88) {
  vapply(
    as.character(x),
    function(s) paste(strwrap(s, width = width), collapse = "\n"),
    character(1)
  )
}

format_or_p_label <- function(or, p) {
  paste0("OR=", sprintf("%.2f", or), "\n", format_p_label(p))
}

plot_top_or_terms <- function(coef_df, model_name, top_n = 12) {
  if (!make_top_or_diagnostic) {
    cat(sprintf("[Plot skipped] %s top-OR diagnostic disabled. Set FOCAL_MODEL_MAKE_TOP_OR=1 to enable.\n", model_name))
    return(invisible(NULL))
  }

  plot_df <- coef_df %>%
    filter(term != "(Intercept)", is.finite(OR), is.finite(z), is.finite(p)) %>%
    mutate(
      abs_z = abs(z),
      p_label = format_p_label(p),
      direction = ifelse(OR >= 1, "OR >= 1", "OR < 1")
    ) %>%
    arrange(desc(abs_z)) %>%
    slice_head(n = top_n) %>%
    arrange(OR) %>%
    mutate(term = factor(term, levels = term))

  if (nrow(plot_df) == 0) {
    cat(sprintf("[Plot skipped] %s has no finite coefficient rows to plot.\n", model_name))
    return(invisible(NULL))
  }

  max_or <- max(plot_df$OR, na.rm = TRUE)
  x_upper <- max(1.2, max_or * 1.12)
  plot_df <- plot_df %>%
    mutate(label_x = pmin(OR + 0.015 * x_upper, 0.98 * x_upper))

  p <- ggplot(plot_df, aes(x = OR, y = term, fill = direction)) +
    geom_col(width = 0.70) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray45") +
    geom_text(
      aes(x = label_x, label = p_label),
      hjust = 0,
      size = 3.2
    ) +
    scale_fill_manual(values = c("OR >= 1" = "#D95F02", "OR < 1" = "#0072B2")) +
    coord_cartesian(xlim = c(0, x_upper), clip = "off") +
    labs(
      title = paste0("gaba_glut_focal_exact_source_models.R :: ", model_name),
      subtitle = wrap_label("glm binomial logit model. Top terms are ranked by absolute robust z statistic; labels show robust p-values."),
      x = "OR",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 14, face = "plain"),
      plot.subtitle = element_text(size = 10),
      panel.grid.major.y = element_line(color = "gray90"),
      panel.grid.minor.y = element_blank(),
      plot.margin = margin(8, 35, 8, 8)
    )

  png_path <- file.path(fig_dir, paste0(model_name, "_top_or_terms.png"))
  pdf_path <- file.path(fig_dir, paste0(model_name, "_top_or_terms.pdf"))
  ggsave(png_path, plot = p, width = 10.5, height = 5.6, dpi = 220)
  ggsave(pdf_path, plot = p, width = 10.5, height = 5.6)

  cat(sprintf("[Plot saved] %s\n", png_path))
  cat(sprintf("[Plot saved] %s\n", pdf_path))
  invisible(p)
}

write_target_class_counts <- function(model_df, model_name) {
  target_levels <- c("Gaba", "Glut")
  observed_counts <- model_df %>%
    count(target_class3, name = "n_rows") %>%
    mutate(target_class3 = as.character(target_class3))

  target_class_counts <- tibble(target_class3 = target_levels) %>%
    left_join(observed_counts, by = "target_class3") %>%
    mutate(
      n_rows = ifelse(is.na(n_rows), 0L, as.integer(n_rows)),
      is_estimable_from_rows = n_rows > 0
    )

  readr::write_csv(target_class_counts, file.path(table_dir, paste0(model_name, "_target_class3_counts.csv")))

  missing_levels <- target_class_counts$target_class3[target_class_counts$n_rows == 0]
  if (length(missing_levels) > 0) {
    cat(sprintf(
      "[Target class check] %s has no rows for target_class3 level(s): %s. These levels cannot be estimated as model coefficients.\n",
      model_name,
      paste(missing_levels, collapse = ", ")
    ))
  }

  target_class_counts
}


lin_comb <- function(model, vcov_mat, weights_named) {
  beta <- coef(model)
  vc_terms <- intersect(names(beta), colnames(vcov_mat))
  vc_terms <- vc_terms[is.finite(beta[vc_terms])]
  beta <- beta[vc_terms]
  L <- rep(0, length(beta))
  names(L) <- names(beta)
  overlap_terms <- intersect(names(weights_named), names(L))
  L[overlap_terms] <- weights_named[overlap_terms]

  used_terms <- names(L)[L != 0]
  if (length(used_terms) == 0 || any(!is.finite(beta[used_terms]))) {
    return(tibble(estimate = NA_real_, se = NA_real_, z = NA_real_, p = NA_real_))
  }

  vc <- vcov_mat[names(beta), names(beta), drop = FALSE]
  var_est <- suppressWarnings(drop(t(L) %*% vc %*% L))
  se <- ifelse(is.finite(var_est) && var_est >= 0, sqrt(var_est), NA_real_)
  estimate <- sum(L * beta)
  z <- estimate / se
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)

  tibble(estimate = estimate, se = se, z = z, p = p)
}

write_target_class_effect_summaries <- function(model, vcov_mat, fit_df, model_name) {
  target_terms <- c(
    Gaba = "target_class3Gaba",
    Glut = "target_class3Glut"
  )
  covariate_terms <- c(
    "source_focal_fraction", "target_gaba_fraction",
    "log_source_size", "log_target_size"
  )
  covariate_means <- vapply(fit_df[covariate_terms], mean, numeric(1), na.rm = TRUE)

  # The three target-class coefficients are intercept-like absolute log-odds
  # at covariates equal to zero. They can all be negative if the modeled
  # containment probability is below 50% at that zero-covariate point. The
  # adjusted summaries below add typical covariate values, and the contrasts
  # below compare target classes directly.
  adjusted_probability_df <- bind_rows(lapply(names(target_terms), function(target_class) {
    weights <- c(setNames(1, target_terms[[target_class]]), covariate_means)
    lc <- lin_comb(model, vcov_mat, weights)
    lc %>%
      mutate(
        target_class3 = target_class,
        adjusted_logit = estimate,
        adjusted_probability = plogis(estimate),
        ci_low_probability = plogis(estimate - 1.96 * se),
        ci_high_probability = plogis(estimate + 1.96 * se)
      ) %>%
      select(target_class3, adjusted_logit, se, z, p, adjusted_probability,
             ci_low_probability, ci_high_probability)
  }))
  readr::write_csv(
    adjusted_probability_df,
    file.path(table_dir, paste0(model_name, "_target_class3_adjusted_probabilities.csv"))
  )

  contrast_defs <- list(
    "Glut_vs_Gaba" = c(target_class3Glut = 1, target_class3Gaba = -1)
  )
  contrast_df <- bind_rows(lapply(names(contrast_defs), function(contrast_name) {
    lc <- lin_comb(model, vcov_mat, contrast_defs[[contrast_name]])
    lc %>%
      mutate(
        contrast = contrast_name,
        OR = exp(estimate),
        ci_low_OR = exp(estimate - 1.96 * se),
        ci_high_OR = exp(estimate + 1.96 * se)
      ) %>%
      select(contrast, estimate, se, z, p, OR, ci_low_OR, ci_high_OR)
  }))
  readr::write_csv(
    contrast_df,
    file.path(table_dir, paste0(model_name, "_target_class3_pairwise_contrasts.csv"))
  )

  invisible(list(
    adjusted_probabilities = adjusted_probability_df,
    pairwise_contrasts = contrast_df
  ))
}

fit_exact_model <- function(model_df, model_name) {
  readr::write_csv(model_df, file.path(table_dir, paste0(model_name, "_model_input.csv")))

  model_input_summary <- model_df %>%
    group_by(focal_type, source_type, target_type, target_class3) %>%
    summarise(
      n_rows = n(),
      exact_source_in_target_n = sum(exact_source_in_target == 1, na.rm = TRUE),
      exact_source_in_target_rate = mean(exact_source_in_target == 1, na.rm = TRUE),
      .groups = "drop"
    )
  readr::write_csv(model_input_summary, file.path(table_dir, paste0(model_name, "_model_input_summary.csv")))
  target_class_counts <- write_target_class_counts(model_df, model_name)

  # Build explicit target-class dummy columns rather than relying on
  # factor contrasts. Primary models compare only neuron-neuron targets:
  # Gaba versus Glut.
  fit_df <- model_df %>%
    mutate(
      target_class3Gaba = as.integer(target_class3 == "Gaba"),
      target_class3Glut = as.integer(target_class3 == "Glut"),
      layer_f = factor(layer_f, levels = c("1", "2/3", "4", "5", "6a", "6b"))
    )

  cost_formula_exact <- exact_source_in_target ~ 0 +
    target_class3Gaba + target_class3Glut +
    source_focal_fraction + target_gaba_fraction +
    log_source_size + log_target_size +
    relevel(layer_f, ref = "1")

  mm_exact <- model.matrix(cost_formula_exact, data = fit_df)
  term_param_count_exact <- tibble(
    term = c(
      "target_class3",
      "source_focal_fraction",
      "target_gaba_fraction",
      "log_source_size",
      "log_target_size",
      "layer_f(ref=1)"
    ),
    n_parameters = c(
      sum(c("target_class3Gaba", "target_class3Glut") %in% colnames(mm_exact)),
      as.integer("source_focal_fraction" %in% colnames(mm_exact)),
      as.integer("target_gaba_fraction" %in% colnames(mm_exact)),
      as.integer("log_source_size" %in% colnames(mm_exact)),
      as.integer("log_target_size" %in% colnames(mm_exact)),
      sum(grepl("^relevel\\(layer_f, ref = \"1\"\\)", colnames(mm_exact)))
    )
  )
  readr::write_csv(term_param_count_exact, file.path(table_dir, paste0(model_name, "_parameter_count_preview.csv")))

  cat(sprintf("\n[Model cost preview | pre-fit] %s exact_source_in_target logistic regression\n", model_name))
  cat(sprintf("Rows used: %d\n", nrow(model_df)))
  cat(sprintf(
    "Outcome positives: %d / %d (%.2f%%)\n",
    sum(model_df$exact_source_in_target == 1, na.rm = TRUE),
    nrow(model_df),
    100 * mean(model_df$exact_source_in_target == 1, na.rm = TRUE)
  ))
  cat(sprintf("Design matrix: %d rows x %d columns\n", nrow(mm_exact), ncol(mm_exact)))
  cat("Note: primary model uses neuron-neuron rows only; target_class3 is encoded as target_class3Gaba and target_class3Glut.\n")
  print(term_param_count_exact, n = nrow(term_param_count_exact))

  fit_time_exact <- system.time({
    m_exact <- glm(
      cost_formula_exact,
      family = binomial(),
      data = fit_df
    )
  })
  cat(sprintf(
    "[Model fit timing] elapsed=%.3fs user=%.3fs sys=%.3fs\n",
    fit_time_exact["elapsed"], fit_time_exact["user.self"], fit_time_exact["sys.self"]
  ))
  print(summary(m_exact))

  vc_exact <- tryCatch(
    sandwich::vcovCL(
      m_exact,
      cluster = data.frame(
        slide_layer = fit_df$slide_layer,
        source_label = fit_df$source_label,
        target_label = fit_df$target_label
      )
    ),
    error = function(e) {
      message("Multi-way cluster-robust vcov failed; falling back to slide_layer clustering: ", conditionMessage(e))
      sandwich::vcovCL(m_exact, cluster = fit_df$slide_layer)
    }
  )
  coef_df <- write_model_outputs(m_exact, vc_exact, paste0(model_name, "_exact_logistic"))
  top_terms_plot <- plot_top_or_terms(coef_df, model_name)
  target_class_effect_summaries <- write_target_class_effect_summaries(m_exact, vc_exact, fit_df, model_name)

  invisible(list(
    model = m_exact,
    robust_vcov = vc_exact,
    coefficients = coef_df,
    top_terms_plot = top_terms_plot,
    parameter_count = term_param_count_exact,
    target_class_counts = target_class_counts,
    target_class_effect_summaries = target_class_effect_summaries
  ))
}

theme_focal_model <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 3),
      plot.subtitle = element_text(hjust = 0, color = "#444444", size = base_size, lineheight = 1.05),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "#111111"),
      axis.line = element_line(color = "#222222", linewidth = 0.45),
      axis.ticks = element_line(color = "#222222", linewidth = 0.35),
      legend.position = "top",
      legend.title = element_blank(),
      strip.background = element_rect(fill = "#F2F2F2", color = NA),
      strip.text = element_text(face = "bold"),
      panel.grid.major.y = element_line(color = "#E7E7E7", linewidth = 0.35),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(12, 28, 12, 14)
    )
}

save_pub_plot <- function(plot, stem, width, height) {
  ggsave(file.path(fig_dir, paste0(stem, ".png")), plot = plot,
         width = width, height = height, dpi = 600, bg = "white", limitsize = FALSE)
  ggsave(file.path(fig_dir, paste0(stem, ".pdf")), plot = plot,
         width = width, height = height, device = grDevices::cairo_pdf,
         bg = "white", limitsize = FALSE)
}

make_observed_rate_table <- function(gaba_df, glut_df) {
  bind_rows(
    gaba_df %>% mutate(model_name = "gaba_focal", focal_label = "Gaba source"),
    glut_df %>% mutate(model_name = "glut_focal", focal_label = "Glut source")
  ) %>%
    group_by(model_name, focal_label, target_class3) %>%
    summarise(
      n_rows = n(),
      exact_n = sum(exact_source_in_target == 1, na.rm = TRUE),
      exact_rate = exact_n / n_rows,
      approx_se = sqrt(pmax(exact_rate * (1 - exact_rate) / n_rows, 0)),
      ci_low = pmax(0, exact_rate - 1.96 * approx_se),
      ci_high = pmin(1, exact_rate + 1.96 * approx_se),
      .groups = "drop"
    ) %>%
    mutate(
      target_class3 = factor(as.character(target_class3), levels = c("Gaba", "Glut")),
      focal_label = factor(focal_label, levels = c("Gaba source", "Glut source")),
      target_label = paste0(as.character(target_class3), " target"),
      rate_label = sprintf("%.1f%%\n%d/%d", 100 * exact_rate, exact_n, n_rows)
    )
}

make_final_summary_plots <- function(gaba_df, glut_df, gaba_results, glut_results) {
  observed_df <- make_observed_rate_table(gaba_df, glut_df)
  readr::write_csv(observed_df, file.path(table_dir, "final_observed_exact_rates_by_target.csv"))

  adjusted_df <- bind_rows(
    gaba_results$target_class_effect_summaries$adjusted_probabilities %>%
      mutate(model_name = "gaba_focal", focal_label = "Gaba source"),
    glut_results$target_class_effect_summaries$adjusted_probabilities %>%
      mutate(model_name = "glut_focal", focal_label = "Glut source")
  ) %>%
    mutate(
      target_class3 = factor(as.character(target_class3), levels = c("Gaba", "Glut")),
      focal_label = factor(focal_label, levels = c("Gaba source", "Glut source")),
      target_label = paste0(as.character(target_class3), " target"),
      prob_label = sprintf("%.1f%%", 100 * adjusted_probability)
    )
  readr::write_csv(adjusted_df, file.path(table_dir, "final_adjusted_probabilities_by_target.csv"))

  contrast_df <- bind_rows(
    gaba_results$target_class_effect_summaries$pairwise_contrasts %>%
      mutate(model_name = "gaba_focal", focal_label = "Gaba source"),
    glut_results$target_class_effect_summaries$pairwise_contrasts %>%
      mutate(model_name = "glut_focal", focal_label = "Glut source")
  ) %>%
    mutate(
      focal_label = factor(focal_label, levels = c("Gaba source", "Glut source")),
      contrast_label = "Glut target vs Gaba target",
      OR_label = format_or_p_label(OR, p),
      label_x = ci_high_OR * 1.18
    )
  readr::write_csv(contrast_df, file.path(table_dir, "final_target_class_glut_vs_gaba_contrasts.csv"))

  target_cols <- c("Gaba" = "#4E79A7", "Glut" = "#D55E00")

  p_observed <- ggplot(observed_df, aes(x = target_class3, y = exact_rate, fill = target_class3)) +
    geom_col(width = 0.62, color = "black", linewidth = 0.25) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.16, linewidth = 0.35) +
    geom_text(aes(label = rate_label), vjust = -0.35, size = 3.3, fontface = "bold", lineheight = 0.9) +
    facet_wrap(~ focal_label, nrow = 1) +
    scale_fill_manual(values = target_cols) +
    scale_y_continuous(labels = function(x) sprintf("%.0f%%", 100 * x),
                       expand = expansion(mult = c(0, 0.18))) +
    labs(
      title = "Observed exact-containment rate by target class",
      subtitle = wrap_label("Rows are neuron-neuron pairs after orienting Gaba or Glut as the focal source. Bars show raw rates; error bars are approximate binomial 95% confidence intervals."),
      x = NULL,
      y = "Exact containment rate"
    ) +
    theme_focal_model(12)
  save_pub_plot(p_observed, "Fig_final_01_observed_exact_rate_by_target", 9.6, 4.8)

  p_adjusted <- ggplot(adjusted_df, aes(x = target_class3, y = adjusted_probability,
                                        color = target_class3)) +
    geom_errorbar(aes(ymin = ci_low_probability, ymax = ci_high_probability),
                  width = 0.12, linewidth = 0.65) +
    geom_point(size = 3.2) +
    geom_text(aes(label = prob_label), vjust = -1.1, size = 3.4, fontface = "bold") +
    facet_wrap(~ focal_label, nrow = 1) +
    scale_color_manual(values = target_cols) +
    scale_y_continuous(labels = function(x) sprintf("%.0f%%", 100 * x),
                       expand = expansion(mult = c(0.08, 0.2))) +
    labs(
      title = "Adjusted probability of exact containment",
      subtitle = wrap_label("Model-adjusted probabilities from focal models at mean source/target size and composition. Layer reference is layer 1."),
      x = NULL,
      y = "Adjusted probability"
    ) +
    theme_focal_model(12)
  save_pub_plot(p_adjusted, "Fig_final_02_adjusted_probability_by_target", 9.6, 4.8)

  p_or <- ggplot(contrast_df, aes(x = OR, y = focal_label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "#777777", linewidth = 0.55) +
    geom_errorbarh(aes(xmin = ci_low_OR, xmax = ci_high_OR),
                   height = 0.18, linewidth = 0.75, color = "#222222") +
    geom_point(size = 3.4, color = "#D55E00") +
    geom_text(aes(x = label_x, label = OR_label), hjust = 0, size = 3.3, fontface = "bold", lineheight = 0.9) +
    scale_x_log10(expand = expansion(mult = c(0.05, 0.38))) +
    coord_cartesian(clip = "off") +
    labs(
      title = "Target-class effect on exact containment",
      subtitle = wrap_label("Focal-model contrast. OR > 1 means the source is more likely to be exactly contained by a Glut target than by a Gaba target; labels show robust p-values."),
      x = "Odds ratio: Glut target vs Gaba target (log scale)",
      y = NULL
    ) +
    theme_focal_model(12)
  save_pub_plot(p_or, "Fig_final_03_glut_vs_gaba_target_OR_forest", 8.8, 4.0)

  invisible(list(
    observed_rates = observed_df,
    adjusted_probabilities = adjusted_df,
    contrasts = contrast_df,
    observed_plot = p_observed,
    adjusted_plot = p_adjusted,
    contrast_plot = p_or
  ))
}

build_main_oriented_input <- function(pair_df) {
  nn_pairs <- pair_df %>%
    filter(
      dataset == "neuron_neuron",
      cluster1_cell_type %in% c("Gaba", "Glut"),
      cluster2_cell_type %in% c("Gaba", "Glut")
    ) %>%
    mutate(
      pair_key = ifelse(
        cluster1_label <= cluster2_label,
        paste(dataset, cluster1_label, cluster2_label, sep = "||"),
        paste(dataset, cluster2_label, cluster1_label, sep = "||")
      )
    )

  oriented_1 <- nn_pairs %>%
    transmute(
      dataset,
      pair_key,
      direction = "cluster1_to_cluster2",
      source_label = cluster1_label,
      target_label = cluster2_label,
      source_type = cluster1_cell_type,
      target_type = cluster2_cell_type,
      source_subclass = cluster1_subclass,
      target_subclass = cluster2_subclass,
      slide_layer = paste(cluster1_slide, cluster1_layer, sep = "__"),
      layer_f = clean_layer(cluster1_layer),
      source_overlap = as.numeric(cluster1_overlap),
      target_overlap = as.numeric(cluster2_overlap),
      exact_source_in_target = as.integer(as.numeric(cluster1_overlap) >= 1),
      source_size = as.numeric(cluster1_total_cell_num),
      target_size = as.numeric(cluster2_total_cell_num),
      log_source_size = log1p(as.numeric(cluster1_total_cell_num)),
      log_target_size = log1p(as.numeric(cluster2_total_cell_num)),
      log_target_source_size_ratio = log1p(as.numeric(cluster2_total_cell_num)) -
        log1p(as.numeric(cluster1_total_cell_num)),
      source_glut_fraction = safe_fraction(cluster1_glut_num, cluster1_total_cell_num),
      source_gaba_fraction = safe_fraction(cluster1_gaba_num, cluster1_total_cell_num),
      target_glut_fraction = safe_fraction(cluster2_glut_num, cluster2_total_cell_num),
      target_gaba_fraction = safe_fraction(cluster2_gaba_num, cluster2_total_cell_num)
    )

  oriented_2 <- nn_pairs %>%
    transmute(
      dataset,
      pair_key,
      direction = "cluster2_to_cluster1",
      source_label = cluster2_label,
      target_label = cluster1_label,
      source_type = cluster2_cell_type,
      target_type = cluster1_cell_type,
      source_subclass = cluster2_subclass,
      target_subclass = cluster1_subclass,
      slide_layer = paste(cluster2_slide, cluster2_layer, sep = "__"),
      layer_f = clean_layer(cluster2_layer),
      source_overlap = as.numeric(cluster2_overlap),
      target_overlap = as.numeric(cluster1_overlap),
      exact_source_in_target = as.integer(as.numeric(cluster2_overlap) >= 1),
      source_size = as.numeric(cluster2_total_cell_num),
      target_size = as.numeric(cluster1_total_cell_num),
      log_source_size = log1p(as.numeric(cluster2_total_cell_num)),
      log_target_size = log1p(as.numeric(cluster1_total_cell_num)),
      log_target_source_size_ratio = log1p(as.numeric(cluster1_total_cell_num)) -
        log1p(as.numeric(cluster2_total_cell_num)),
      source_glut_fraction = safe_fraction(cluster2_glut_num, cluster2_total_cell_num),
      source_gaba_fraction = safe_fraction(cluster2_gaba_num, cluster2_total_cell_num),
      target_glut_fraction = safe_fraction(cluster1_glut_num, cluster1_total_cell_num),
      target_gaba_fraction = safe_fraction(cluster1_gaba_num, cluster1_total_cell_num)
    )

  bind_rows(oriented_1, oriented_2) %>%
    mutate(
      source_type = factor(source_type, levels = c("Gaba", "Glut")),
      target_type = factor(target_type, levels = c("Gaba", "Glut")),
      layer_f = factor(layer_f, levels = c("1", "2/3", "4", "5", "6a", "6b")),
      pair_type = paste0(as.character(source_type), "_to_", as.character(target_type)),
      same_type = source_type == target_type,
      same_subclass = source_subclass == target_subclass,
      source_focal_fraction = ifelse(
        source_type == "Gaba",
        source_gaba_fraction,
        source_glut_fraction
      )
    ) %>%
    filter(
      !is.na(source_type), !is.na(target_type), !is.na(layer_f),
      is.finite(exact_source_in_target),
      is.finite(log_source_size), is.finite(log_target_size),
      is.finite(log_target_source_size_ratio),
      is.finite(source_focal_fraction),
      is.finite(target_gaba_fraction)
    )
}

safe_cluster_vcov <- function(model, fit_df) {
  tryCatch(
    sandwich::vcovCL(
      model,
      cluster = data.frame(
        pair_key = fit_df$pair_key,
        slide_layer = fit_df$slide_layer,
        source_label = fit_df$source_label,
        target_label = fit_df$target_label
      )
    ),
    error = function(e) {
      message("Multi-way cluster-robust vcov failed; falling back to pair_key clustering: ", conditionMessage(e))
      sandwich::vcovCL(model, cluster = fit_df$pair_key)
    }
  )
}

tidy_robust_glm <- function(model, vcov_mat) {
  coef_test <- lmtest::coeftest(model, vcov. = vcov_mat)
  tibble(
    term = rownames(coef_test),
    estimate = coef_test[, 1],
    se = coef_test[, 2],
    z = coef_test[, 3],
    p = coef_test[, 4]
  ) %>%
    filter(is.finite(estimate), is.finite(se)) %>%
    mutate(
      OR = exp(estimate),
      ci_low = exp(estimate - 1.96 * se),
      ci_high = exp(estimate + 1.96 * se)
    )
}

predict_glm_robust <- function(model, vcov_mat, newdata) {
  tt <- delete.response(terms(model))
  mm <- model.matrix(tt, newdata)
  beta <- coef(model)
  keep_terms <- intersect(colnames(mm), names(beta))
  keep_terms <- intersect(keep_terms, colnames(vcov_mat))
  keep_terms <- keep_terms[is.finite(beta[keep_terms])]
  mm <- mm[, keep_terms, drop = FALSE]
  beta <- beta[keep_terms]
  vc <- vcov_mat[keep_terms, keep_terms, drop = FALSE]
  eta <- drop(mm %*% beta)
  se <- sqrt(pmax(diag(mm %*% vc %*% t(mm)), 0))
  newdata %>%
    mutate(
      adjusted_logit = eta,
      se = se,
      adjusted_probability = plogis(eta),
      ci_low_probability = plogis(eta - 1.96 * se),
      ci_high_probability = plogis(eta + 1.96 * se)
    )
}

contrast_glm_robust <- function(model, vcov_mat, newdata_a, newdata_b, contrast_name) {
  tt <- delete.response(terms(model))
  mm_a <- model.matrix(tt, newdata_a)
  mm_b <- model.matrix(tt, newdata_b)
  beta <- coef(model)
  keep_terms <- Reduce(intersect, list(colnames(mm_a), colnames(mm_b), names(beta), colnames(vcov_mat)))
  keep_terms <- keep_terms[is.finite(beta[keep_terms])]
  L <- mm_a[, keep_terms, drop = FALSE] - mm_b[, keep_terms, drop = FALSE]
  beta <- beta[keep_terms]
  vc <- vcov_mat[keep_terms, keep_terms, drop = FALSE]
  estimate <- drop(L %*% beta)
  se <- sqrt(pmax(diag(L %*% vc %*% t(L)), 0))
  tibble(
    contrast = contrast_name,
    estimate = estimate,
    se = se,
    z = estimate / se,
    p = 2 * pnorm(abs(estimate / se), lower.tail = FALSE),
    OR = exp(estimate),
    ci_low_OR = exp(estimate - 1.96 * se),
    ci_high_OR = exp(estimate + 1.96 * se)
  )
}

make_reference_grid <- function(model_df) {
  expand.grid(
    source_type = factor(c("Gaba", "Glut"), levels = c("Gaba", "Glut")),
    target_type = factor(c("Gaba", "Glut"), levels = c("Gaba", "Glut")),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    mutate(
      layer_f = factor("1", levels = c("1", "2/3", "4", "5", "6a", "6b")),
      log_source_size = mean(model_df$log_source_size, na.rm = TRUE),
      log_target_size = mean(model_df$log_target_size, na.rm = TRUE),
      log_target_source_size_ratio = mean(model_df$log_target_source_size_ratio, na.rm = TRUE),
      source_focal_fraction = mean(model_df$source_focal_fraction, na.rm = TRUE),
      target_gaba_fraction = mean(model_df$target_gaba_fraction, na.rm = TRUE)
    )
}

fit_main_factor_model <- function(model_df, model_name, formula_obj) {
  readr::write_csv(model_df, file.path(table_dir, paste0(model_name, "_model_input.csv")))
  input_summary <- model_df %>%
    group_by(source_type, target_type) %>%
    summarise(
      n_rows = n(),
      exact_n = sum(exact_source_in_target == 1, na.rm = TRUE),
      exact_rate = exact_n / n_rows,
      .groups = "drop"
    )
  readr::write_csv(input_summary, file.path(table_dir, paste0(model_name, "_input_summary.csv")))

  fit <- glm(formula_obj, family = binomial(), data = model_df)
  vc <- safe_cluster_vcov(fit, model_df)
  coef_df <- tidy_robust_glm(fit, vc)
  readr::write_csv(coef_df, file.path(table_dir, paste0(model_name, "_coefficients.csv")))

  grid <- make_reference_grid(model_df)
  pred <- predict_glm_robust(fit, vc, grid) %>%
    mutate(
      model_name = model_name,
      source_target = paste0(as.character(source_type), " -> ", as.character(target_type)),
      prob_label = sprintf("%.1f%%", 100 * adjusted_probability)
    )
  readr::write_csv(pred, file.path(table_dir, paste0(model_name, "_adjusted_source_target_probabilities.csv")))

  contrast_defs <- list(
    "Gaba_source_Glut_target_vs_Gaba_target" = list(
      a = grid %>% filter(source_type == "Gaba", target_type == "Glut"),
      b = grid %>% filter(source_type == "Gaba", target_type == "Gaba")
    ),
    "Glut_source_Glut_target_vs_Gaba_target" = list(
      a = grid %>% filter(source_type == "Glut", target_type == "Glut"),
      b = grid %>% filter(source_type == "Glut", target_type == "Gaba")
    ),
    "Cross_direction_Gaba_to_Glut_vs_Glut_to_Gaba" = list(
      a = grid %>% filter(source_type == "Gaba", target_type == "Glut"),
      b = grid %>% filter(source_type == "Glut", target_type == "Gaba")
    )
  )
  contrast_df <- bind_rows(lapply(names(contrast_defs), function(nm) {
    contrast_glm_robust(fit, vc, contrast_defs[[nm]]$a, contrast_defs[[nm]]$b, nm)
  }))
  readr::write_csv(contrast_df, file.path(table_dir, paste0(model_name, "_key_contrasts.csv")))

  invisible(list(
    model = fit,
    robust_vcov = vc,
    coefficients = coef_df,
    input_summary = input_summary,
    adjusted_probabilities = pred,
    key_contrasts = contrast_df
  ))
}

label_main_terms <- function(x) {
  dplyr::case_when(
    x == "source_typeGlut" ~ "Source type: Glut vs Gaba",
    x == "target_typeGlut" ~ "Target type: Glut vs Gaba",
    x == "source_typeGlut:target_typeGlut" ~ "Interaction: Glut source x Glut target",
    x == "log_source_size" ~ "log(source size + 1)",
    x == "log_target_size" ~ "log(target size + 1)",
    x == "log_target_source_size_ratio" ~ "log target/source size ratio",
    x == "source_focal_fraction" ~ "Source focal-cell fraction",
    x == "target_gaba_fraction" ~ "Target GABA-cell fraction",
    grepl("relevel\\(layer_f, ref = \"1\"\\)", x) ~
      paste0("Layer ", sub("relevel\\(layer_f, ref = \"1\"\\)", "", x), " vs layer 1"),
    TRUE ~ x
  )
}

plot_main_factor_outputs <- function(model_df, total_results, adjusted_results) {
  observed_df <- model_df %>%
    group_by(source_type, target_type) %>%
    summarise(
      n_rows = n(),
      exact_n = sum(exact_source_in_target == 1, na.rm = TRUE),
      exact_rate = exact_n / n_rows,
      .groups = "drop"
    ) %>%
    mutate(
      rate_label = sprintf("%.1f%%\n%d/%d", 100 * exact_rate, exact_n, n_rows),
      source_type = factor(source_type, levels = c("Gaba", "Glut")),
      target_type = factor(target_type, levels = c("Gaba", "Glut"))
    )
  readr::write_csv(observed_df, file.path(table_dir, "main_observed_source_target_exact_rates.csv"))

  p_observed_matrix <- ggplot(observed_df, aes(x = target_type, y = source_type, fill = exact_rate)) +
    geom_tile(color = "white", linewidth = 1.1) +
    geom_text(aes(label = rate_label), size = 4.1, fontface = "bold", lineheight = 0.9) +
    scale_fill_gradient(low = "#F6E5C8", high = "#D55E00", labels = function(x) sprintf("%.0f%%", 100 * x)) +
    labs(
      title = "Observed exact-containment rate by source and target type",
      subtitle = wrap_label("Each cell is a directed neuron-neuron pair. The source cluster is tested for complete containment in the target cluster; labels show raw rate and exact/total counts."),
      x = "Target / potential host",
      y = "Source / tested cluster",
      fill = "Exact rate"
    ) +
    theme_focal_model(12)
  save_pub_plot(p_observed_matrix, "Fig_main_00_observed_source_target_exact_rate_matrix", 8.2, 5.8)

  forest_df <- total_results$coefficients %>%
    filter(term != "(Intercept)", is.finite(OR), is.finite(ci_low), is.finite(ci_high)) %>%
    mutate(
      term_label = label_main_terms(term),
      term_label = factor(term_label, levels = rev(term_label)),
      direction = ifelse(OR >= 1, "Increases exact containment", "Decreases exact containment"),
      p_label = format_p_label(p),
      label_x = ci_high * 1.18
    )

  p_forest <- ggplot(forest_df, aes(x = OR, y = term_label, color = direction)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "#777777", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.18, linewidth = 0.7) +
    geom_point(size = 2.8) +
    geom_text(aes(x = label_x, label = p_label), hjust = 0, size = 3.0, color = "#222222") +
    scale_x_log10(expand = expansion(mult = c(0.05, 0.42))) +
    scale_color_manual(values = c("Increases exact containment" = "#D55E00", "Decreases exact containment" = "#0072B2")) +
    coord_cartesian(clip = "off") +
    labs(
      title = "Factors associated with exact containment",
      subtitle = wrap_label("Main total-effect logistic model. Points are odds ratios; bars are robust 95% confidence intervals clustered by pair, slide-layer, source cluster and target cluster. Labels show robust p-values."),
      x = "Odds ratio (log scale); label = robust p-value",
      y = NULL
    ) +
    theme_focal_model(12)
  save_pub_plot(p_forest, "Fig_main_01_total_effect_factor_OR_forest", 11.4, 6.6)

  pred_df <- total_results$adjusted_probabilities %>%
    mutate(
      source_type = factor(source_type, levels = c("Gaba", "Glut")),
      target_type = factor(target_type, levels = c("Gaba", "Glut"))
    )

  p_pred_matrix <- ggplot(pred_df, aes(x = target_type, y = source_type, fill = adjusted_probability)) +
    geom_tile(color = "white", linewidth = 1.1) +
    geom_text(aes(label = prob_label), size = 4.5, fontface = "bold") +
    scale_fill_gradient(low = "#F6E5C8", high = "#D55E00", labels = function(x) sprintf("%.0f%%", 100 * x)) +
    labs(
      title = "Adjusted probability of exact containment",
      subtitle = wrap_label("Predicted from the main total-effect model at mean source/target size covariates and reference layer 1. Values are adjusted probabilities, not raw rates."),
      x = "Target / potential host",
      y = "Source / tested cluster",
      fill = "Adjusted probability"
    ) +
    theme_focal_model(12)
  save_pub_plot(p_pred_matrix, "Fig_main_02_total_effect_adjusted_probability_matrix", 8.2, 5.8)

  contrast_df <- total_results$key_contrasts %>%
    mutate(
      contrast_label = dplyr::recode(
        contrast,
        "Gaba_source_Glut_target_vs_Gaba_target" = "Gaba source: Glut vs Gaba target",
        "Glut_source_Glut_target_vs_Gaba_target" = "Glut source: Glut vs Gaba target",
        "Cross_direction_Gaba_to_Glut_vs_Glut_to_Gaba" = "Cross direction: Gaba->Glut vs Glut->Gaba"
      ),
      contrast_label = factor(contrast_label, levels = rev(contrast_label)),
      OR_label = format_or_p_label(OR, p),
      label_x = ci_high_OR * 1.18
    )

  p_contrast <- ggplot(contrast_df, aes(x = OR, y = contrast_label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "#777777", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = ci_low_OR, xmax = ci_high_OR), height = 0.18, linewidth = 0.75) +
    geom_point(size = 3.2, color = "#D55E00") +
    geom_text(aes(x = label_x, label = OR_label), hjust = 0, size = 3.2, fontface = "bold", lineheight = 0.9) +
    scale_x_log10(expand = expansion(mult = c(0.05, 0.42))) +
    coord_cartesian(clip = "off") +
    labs(
      title = "Key source-target contrasts for exact containment",
      subtitle = wrap_label("Contrasts are computed from the main total-effect model. OR > 1 means the first direction/class in the label has higher odds of exact containment; labels show robust p-values."),
      x = "Odds ratio (log scale)",
      y = NULL
    ) +
    theme_focal_model(12)
  save_pub_plot(p_contrast, "Fig_main_03_total_effect_key_contrasts", 10.8, 4.9)

  adjusted_forest_df <- adjusted_results$coefficients %>%
    filter(term != "(Intercept)", is.finite(OR), is.finite(ci_low), is.finite(ci_high)) %>%
    mutate(
      term_label = label_main_terms(term),
      term_label = factor(term_label, levels = rev(term_label)),
      direction = ifelse(OR >= 1, "Increases exact containment", "Decreases exact containment"),
      p_label = format_p_label(p),
      label_x = ci_high * 1.18
    )

  p_adjusted_forest <- ggplot(adjusted_forest_df, aes(x = OR, y = term_label, color = direction)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "#777777", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.18, linewidth = 0.7) +
    geom_point(size = 2.8) +
    geom_text(aes(x = label_x, label = p_label), hjust = 0, size = 3.0, color = "#222222") +
    scale_x_log10(expand = expansion(mult = c(0.05, 0.42))) +
    scale_color_manual(values = c("Increases exact containment" = "#D55E00", "Decreases exact containment" = "#0072B2")) +
    coord_cartesian(clip = "off") +
    labs(
      title = "Composition-adjusted factor model",
      subtitle = wrap_label("Supplementary logistic model adding source focal-cell fraction and target GABA-cell fraction. Points are odds ratios; bars are robust 95% confidence intervals; labels show robust p-values."),
      x = "Odds ratio (log scale); label = robust p-value",
      y = NULL
    ) +
    theme_focal_model(12)
  save_pub_plot(p_adjusted_forest, "Fig_supp_01_composition_adjusted_factor_OR_forest", 11.4, 6.8)

  invisible(list(
    observed_matrix = p_observed_matrix,
    total_forest = p_forest,
    total_probability_matrix = p_pred_matrix,
    total_contrasts = p_contrast,
    composition_adjusted_forest = p_adjusted_forest
  ))
}

plot_stratified_source_target_outputs <- function(model_df, total_results) {
  observed_df <- model_df %>%
    group_by(source_type, target_type) %>%
    summarise(
      n_rows = n(),
      exact_n = sum(exact_source_in_target == 1, na.rm = TRUE),
      exact_rate = exact_n / n_rows,
      approx_se = sqrt(pmax(exact_rate * (1 - exact_rate) / n_rows, 0)),
      ci_low = pmax(0, exact_rate - 1.96 * approx_se),
      ci_high = pmin(1, exact_rate + 1.96 * approx_se),
      .groups = "drop"
    ) %>%
    mutate(
      source_type = factor(source_type, levels = c("Gaba", "Glut")),
      target_type = factor(target_type, levels = c("Gaba", "Glut")),
      rate_label = sprintf("%.1f%%\n%d/%d", 100 * exact_rate, exact_n, n_rows)
    )
  readr::write_csv(
    observed_df,
    file.path(table_dir, "extra_stratified_observed_source_target_rates.csv")
  )

  adjusted_df <- total_results$adjusted_probabilities %>%
    mutate(
      source_type = factor(source_type, levels = c("Gaba", "Glut")),
      target_type = factor(target_type, levels = c("Gaba", "Glut")),
      prob_label = sprintf("%.1f%%", 100 * adjusted_probability)
    )
  readr::write_csv(
    adjusted_df,
    file.path(table_dir, "extra_stratified_total_effect_adjusted_probabilities.csv")
  )

  source_contrast_df <- total_results$key_contrasts %>%
    filter(contrast %in% c(
      "Gaba_source_Glut_target_vs_Gaba_target",
      "Glut_source_Glut_target_vs_Gaba_target"
    )) %>%
    mutate(
      source_type = factor(
        ifelse(grepl("^Gaba", contrast), "Gaba", "Glut"),
        levels = c("Gaba", "Glut")
      ),
      target_type = factor("Glut", levels = c("Gaba", "Glut")),
      contrast_label = "Glut target vs Gaba target",
      label = format_or_p_label(OR, p)
    )
  readr::write_csv(
    source_contrast_df,
    file.path(table_dir, "extra_stratified_source_target_contrasts.csv")
  )

  observed_source_ann <- observed_df %>%
    group_by(source_type) %>%
    summarise(ann_y = min(1, max(ci_high, na.rm = TRUE) * 1.18 + 0.015), .groups = "drop") %>%
    left_join(source_contrast_df, by = "source_type")

  adjusted_source_ann <- adjusted_df %>%
    group_by(source_type) %>%
    summarise(
      ann_y = min(1, max(ci_high_probability, na.rm = TRUE) * 1.18 + 0.015),
      .groups = "drop"
    ) %>%
    left_join(source_contrast_df, by = "source_type")

  type_cols <- c("Gaba" = "#4E79A7", "Glut" = "#D55E00")

  p_observed_by_source <- ggplot(
    observed_df,
    aes(x = target_type, y = exact_rate, fill = target_type)
  ) +
    geom_col(width = 0.62, color = "black", linewidth = 0.25) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.16, linewidth = 0.35) +
    geom_text(aes(label = rate_label), vjust = -0.35, size = 3.35, fontface = "bold", lineheight = 0.9) +
    geom_label(
      data = observed_source_ann,
      aes(x = target_type, y = ann_y, label = label),
      inherit.aes = FALSE,
      hjust = 0.5,
      size = 3.1,
      label.size = 0.2,
      fill = "white",
      lineheight = 0.9
    ) +
    facet_wrap(~ source_type, nrow = 1, labeller = label_both) +
    scale_fill_manual(values = type_cols) +
    scale_y_continuous(labels = function(x) sprintf("%.0f%%", 100 * x),
                       expand = expansion(mult = c(0, 0.28))) +
    labs(
      title = "Observed exact-containment rate split by source type",
      subtitle = wrap_label("For each source class, Gaba and Glut targets are shown separately. Labels show raw exact-containment rates and counts; boxed labels show main-model target-class contrasts."),
      x = "Target / potential host",
      y = "Observed exact-containment rate"
    ) +
    theme_focal_model(12)
  save_pub_plot(
    p_observed_by_source,
    "Fig_extra_01_observed_exact_rate_split_by_source",
    9.6,
    4.9
  )

  p_adjusted_by_source <- ggplot(
    adjusted_df,
    aes(x = target_type, y = adjusted_probability, color = target_type)
  ) +
    geom_errorbar(aes(ymin = ci_low_probability, ymax = ci_high_probability),
                  width = 0.14, linewidth = 0.65) +
    geom_point(size = 3.3) +
    geom_text(aes(label = prob_label), vjust = -1.05, size = 3.4, fontface = "bold") +
    geom_label(
      data = adjusted_source_ann,
      aes(x = target_type, y = ann_y, label = label),
      inherit.aes = FALSE,
      hjust = 0.5,
      size = 3.1,
      label.size = 0.2,
      fill = "white",
      lineheight = 0.9
    ) +
    facet_wrap(~ source_type, nrow = 1, labeller = label_both) +
    scale_color_manual(values = type_cols) +
    scale_y_continuous(labels = function(x) sprintf("%.0f%%", 100 * x),
                       expand = expansion(mult = c(0.06, 0.3))) +
    labs(
      title = "Adjusted exact-containment probability split by source type",
      subtitle = wrap_label("Predictions from the main total-effect model at mean source/target size and reference layer 1. Boxed labels compare Glut target versus Gaba target within each source class."),
      x = "Target / potential host",
      y = "Adjusted probability"
    ) +
    theme_focal_model(12)
  save_pub_plot(
    p_adjusted_by_source,
    "Fig_extra_02_adjusted_probability_split_by_source",
    9.6,
    4.9
  )

  p_observed_by_target <- ggplot(
    observed_df,
    aes(x = source_type, y = exact_rate, fill = source_type)
  ) +
    geom_col(width = 0.62, color = "black", linewidth = 0.25) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.16, linewidth = 0.35) +
    geom_text(aes(label = rate_label), vjust = -0.35, size = 3.35, fontface = "bold", lineheight = 0.9) +
    facet_wrap(~ target_type, nrow = 1, labeller = label_both) +
    scale_fill_manual(values = type_cols) +
    scale_y_continuous(labels = function(x) sprintf("%.0f%%", 100 * x),
                       expand = expansion(mult = c(0, 0.18))) +
    labs(
      title = "Observed exact-containment rate split by target type",
      subtitle = wrap_label("For each target class, Gaba and Glut sources are shown separately. This view emphasizes how source identity behaves within the same potential host class."),
      x = "Source / tested cluster",
      y = "Observed exact-containment rate"
    ) +
    theme_focal_model(12)
  save_pub_plot(
    p_observed_by_target,
    "Fig_extra_03_observed_exact_rate_split_by_target",
    9.6,
    4.9
  )

  p_adjusted_by_target <- ggplot(
    adjusted_df,
    aes(x = source_type, y = adjusted_probability, color = source_type)
  ) +
    geom_errorbar(aes(ymin = ci_low_probability, ymax = ci_high_probability),
                  width = 0.14, linewidth = 0.65) +
    geom_point(size = 3.3) +
    geom_text(aes(label = prob_label), vjust = -1.05, size = 3.4, fontface = "bold") +
    facet_wrap(~ target_type, nrow = 1, labeller = label_both) +
    scale_color_manual(values = type_cols) +
    scale_y_continuous(labels = function(x) sprintf("%.0f%%", 100 * x),
                       expand = expansion(mult = c(0.06, 0.22))) +
    labs(
      title = "Adjusted exact-containment probability split by target type",
      subtitle = wrap_label("Predictions from the main total-effect model at mean source/target size and reference layer 1. This separates Gaba and Glut sources within each target class."),
      x = "Source / tested cluster",
      y = "Adjusted probability"
    ) +
    theme_focal_model(12)
  save_pub_plot(
    p_adjusted_by_target,
    "Fig_extra_04_adjusted_probability_split_by_target",
    9.6,
    4.9
  )

  invisible(list(
    observed_by_source = p_observed_by_source,
    adjusted_by_source = p_adjusted_by_source,
    observed_by_target = p_observed_by_target,
    adjusted_by_target = p_adjusted_by_target,
    observed_data = observed_df,
    adjusted_data = adjusted_df,
    source_target_contrasts = source_contrast_df
  ))
}

fit_source_specific_factor_model <- function(model_df, source_value, model_name, formula_obj) {
  source_df <- model_df %>%
    filter(source_type == source_value) %>%
    mutate(
      source_type = factor(source_type, levels = c("Gaba", "Glut")),
      target_type = factor(target_type, levels = c("Gaba", "Glut")),
      layer_f = factor(layer_f, levels = c("1", "2/3", "4", "5", "6a", "6b"))
    )

  readr::write_csv(source_df, file.path(table_dir, paste0(model_name, "_model_input.csv")))

  input_summary <- source_df %>%
    group_by(source_type, target_type) %>%
    summarise(
      n_rows = n(),
      exact_n = sum(exact_source_in_target == 1, na.rm = TRUE),
      exact_rate = exact_n / n_rows,
      .groups = "drop"
    )
  readr::write_csv(input_summary, file.path(table_dir, paste0(model_name, "_input_summary.csv")))

  fit <- glm(formula_obj, family = binomial(), data = source_df)
  vc <- safe_cluster_vcov(fit, source_df)
  coef_df <- tidy_robust_glm(fit, vc)
  readr::write_csv(coef_df, file.path(table_dir, paste0(model_name, "_coefficients.csv")))

  invisible(list(
    model = fit,
    robust_vcov = vc,
    coefficients = coef_df,
    input_summary = input_summary,
    source_df = source_df
  ))
}

plot_source_specific_factor_forest <- function(results, source_value, model_label, stem) {
  forest_df <- results$coefficients %>%
    filter(term != "(Intercept)", is.finite(OR), is.finite(ci_low), is.finite(ci_high)) %>%
    mutate(
      term_label = label_main_terms(term),
      term_label = factor(term_label, levels = rev(term_label)),
      direction = ifelse(OR >= 1, "Increases exact containment", "Decreases exact containment"),
      p_label = format_p_label(p),
      label_x = ci_high * 1.18
    )

  if (nrow(forest_df) == 0) {
    cat(sprintf("[Plot skipped] %s has no finite coefficient rows to plot.\n", stem))
    return(invisible(NULL))
  }

  n_rows <- nrow(results$source_df)
  n_positive <- sum(results$source_df$exact_source_in_target == 1, na.rm = TRUE)
  positive_rate <- mean(results$source_df$exact_source_in_target == 1, na.rm = TRUE)

  p <- ggplot(forest_df, aes(x = OR, y = term_label, color = direction)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "#777777", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.18, linewidth = 0.7) +
    geom_point(size = 2.8) +
    geom_text(aes(x = label_x, label = p_label), hjust = 0, size = 3.0, color = "#222222") +
    scale_x_log10(expand = expansion(mult = c(0.05, 0.42))) +
    scale_color_manual(values = c("Increases exact containment" = "#D55E00", "Decreases exact containment" = "#0072B2")) +
    coord_cartesian(clip = "off") +
    labs(
      title = paste0(source_value, " source: ", model_label),
      subtitle = wrap_label(sprintf(
        "Source-specific logistic model fitted only on %s-source directed neuron-neuron rows (n=%d; exact=%d, %.1f%%). Source type and source-target interaction are removed because source type is fixed. Points are odds ratios; bars are robust 95%% confidence intervals; labels show robust p-values.",
        source_value,
        n_rows,
        n_positive,
        100 * positive_rate
      )),
      x = "Odds ratio (log scale); label = robust p-value",
      y = NULL
    ) +
    theme_focal_model(12)

  save_pub_plot(p, stem, 11.4, 6.2)
  invisible(p)
}

plot_source_specific_factor_outputs <- function(model_df) {
  source_specific_total_formula <- exact_source_in_target ~
    target_type +
    log_source_size + log_target_size +
    relevel(layer_f, ref = "1")

  source_specific_adjusted_formula <- exact_source_in_target ~
    target_type +
    log_source_size + log_target_size +
    source_focal_fraction + target_gaba_fraction +
    relevel(layer_f, ref = "1")

  gaba_total_results <- fit_source_specific_factor_model(
    model_df,
    "Gaba",
    "extra_gaba_source_total_effect",
    source_specific_total_formula
  )
  glut_total_results <- fit_source_specific_factor_model(
    model_df,
    "Glut",
    "extra_glut_source_total_effect",
    source_specific_total_formula
  )
  gaba_adjusted_results <- fit_source_specific_factor_model(
    model_df,
    "Gaba",
    "extra_gaba_source_composition_adjusted",
    source_specific_adjusted_formula
  )
  glut_adjusted_results <- fit_source_specific_factor_model(
    model_df,
    "Glut",
    "extra_glut_source_composition_adjusted",
    source_specific_adjusted_formula
  )

  gaba_total_plot <- plot_source_specific_factor_forest(
    gaba_total_results,
    "Gaba",
    "total-effect factor model",
    "Fig_extra_05_gaba_source_total_effect_factor_OR_forest"
  )
  glut_total_plot <- plot_source_specific_factor_forest(
    glut_total_results,
    "Glut",
    "total-effect factor model",
    "Fig_extra_06_glut_source_total_effect_factor_OR_forest"
  )
  gaba_adjusted_plot <- plot_source_specific_factor_forest(
    gaba_adjusted_results,
    "Gaba",
    "composition-adjusted factor model",
    "Fig_extra_07_gaba_source_composition_adjusted_factor_OR_forest"
  )
  glut_adjusted_plot <- plot_source_specific_factor_forest(
    glut_adjusted_results,
    "Glut",
    "composition-adjusted factor model",
    "Fig_extra_08_glut_source_composition_adjusted_factor_OR_forest"
  )

  invisible(list(
    gaba_total_results = gaba_total_results,
    glut_total_results = glut_total_results,
    gaba_adjusted_results = gaba_adjusted_results,
    glut_adjusted_results = glut_adjusted_results,
    gaba_total_plot = gaba_total_plot,
    glut_total_plot = glut_total_plot,
    gaba_adjusted_plot = gaba_adjusted_plot,
    glut_adjusted_plot = glut_adjusted_plot
  ))
}

# -----------------------------
# Read all rows and canonicalize column names
# -----------------------------
if (!file.exists(input_nn)) {
  stop("Cannot find neuron-neuron input file: ", input_nn, call. = FALSE)
}
nn_raw <- readr::read_csv(input_nn, show_col_types = FALSE)

nn <- rename_partner_columns(nn_raw) %>%
  mutate(dataset = "neuron_neuron") %>%
  mutate(pair_key = paste(dataset, cluster1_label, cluster2_label, sep = "||")) %>%
  distinct(pair_key, .keep_all = TRUE) %>%
  select(-pair_key)

nonn <- tibble()
if (file.exists(input_nonn)) {
  nonn_raw <- readr::read_csv(input_nonn, show_col_types = FALSE)
  nonn <- rename_partner_columns(nonn_raw) %>%
    mutate(dataset = "neuron_nonneuron") %>%
    mutate(pair_key = paste(dataset, cluster1_label, cluster2_label, sep = "||")) %>%
    distinct(pair_key, .keep_all = TRUE) %>%
    select(-pair_key)
}

all_pairs <- bind_rows(nn, nonn)
readr::write_csv(all_pairs, file.path(table_dir, "all_input_rows_no_type_filter.csv"))
readr::write_csv(
  all_pairs %>% count(dataset, cluster1_cell_type, cluster2_cell_type, name = "n_rows"),
  file.path(table_dir, "all_input_rows_dataset_type_audit.csv")
)

# -----------------------------
# Build and fit the two focal models
# -----------------------------
main_model_df <- build_main_oriented_input(all_pairs)
readr::write_csv(main_model_df, file.path(table_dir, "main_oriented_factor_model_input.csv"))

main_model_manifest <- tibble(
  model_name = c("main_total_effect", "main_composition_adjusted"),
  purpose = c(
    "Total effects of source type, target type, source size, target size, and layer",
    "Composition-adjusted effects adding source focal fraction and target GABA fraction"
  ),
  n_rows = nrow(main_model_df),
  n_positive = sum(main_model_df$exact_source_in_target == 1, na.rm = TRUE),
  positive_rate = mean(main_model_df$exact_source_in_target == 1, na.rm = TRUE)
)
readr::write_csv(main_model_manifest, file.path(table_dir, "main_factor_model_manifest.csv"))

cat("\n[Fitting main model 1 of 2] main_total_effect\n")
main_total_formula <- exact_source_in_target ~
  source_type * target_type +
  log_source_size + log_target_size +
  relevel(layer_f, ref = "1")
main_total_results <- fit_main_factor_model(
  main_model_df,
  "main_total_effect",
  main_total_formula
)

cat("\n[Fitting main model 2 of 2] main_composition_adjusted\n")
main_adjusted_formula <- exact_source_in_target ~
  source_type * target_type +
  log_source_size + log_target_size +
  source_focal_fraction + target_gaba_fraction +
  relevel(layer_f, ref = "1")
main_adjusted_results <- fit_main_factor_model(
  main_model_df,
  "main_composition_adjusted",
  main_adjusted_formula
)

cat("\n[Making main factor-model figures]\n")
main_factor_summary <- plot_main_factor_outputs(main_model_df, main_total_results, main_adjusted_results)

cat("\n[Making source-specific factor-model figures]\n")
source_specific_factor_summary <- plot_source_specific_factor_outputs(main_model_df)

# -----------------------------
# Supplementary focal models
# -----------------------------
gaba_model_df <- make_focal_model_input(all_pairs, "Gaba")
glut_model_df <- make_focal_model_input(all_pairs, "Glut")

readr::write_csv(
  bind_rows(
    gaba_model_df %>% mutate(model_name = "gaba_focal"),
    glut_model_df %>% mutate(model_name = "glut_focal")
  ),
  file.path(table_dir, "combined_gaba_glut_focal_model_inputs.csv")
)

model_run_manifest <- tibble(
  model_name = c("gaba_focal", "glut_focal"),
  focal_type = c("Gaba", "Glut"),
  n_rows = c(nrow(gaba_model_df), nrow(glut_model_df))
)
readr::write_csv(model_run_manifest, file.path(table_dir, "model_run_manifest.csv"))

cat("\n[Fitting supplementary focal model 1 of 2] gaba_focal\n")
gaba_results <- fit_exact_model(gaba_model_df, "gaba_focal")
cat("\n[Fitting supplementary focal model 2 of 2] glut_focal\n")
glut_results <- fit_exact_model(glut_model_df, "glut_focal")

cat("\n[Making final presentation figures]\n")
final_summary <- make_final_summary_plots(gaba_model_df, glut_model_df, gaba_results, glut_results)

cat("\n[Output files]\n")
cat(sprintf("Main factor-model manifest: %s\n", file.path(table_dir, "main_factor_model_manifest.csv")))
cat(sprintf("Main oriented model input: %s\n", file.path(table_dir, "main_oriented_factor_model_input.csv")))
cat(sprintf("Main total-effect coefficients: %s\n", file.path(table_dir, "main_total_effect_coefficients.csv")))
cat(sprintf("Main composition-adjusted coefficients: %s\n", file.path(table_dir, "main_composition_adjusted_coefficients.csv")))
cat(sprintf("Main observed source-target matrix: %s\n", file.path(fig_dir, "Fig_main_00_observed_source_target_exact_rate_matrix.png")))
cat(sprintf("Main total-effect OR forest: %s\n", file.path(fig_dir, "Fig_main_01_total_effect_factor_OR_forest.png")))
cat(sprintf("Main adjusted probability matrix: %s\n", file.path(fig_dir, "Fig_main_02_total_effect_adjusted_probability_matrix.png")))
cat(sprintf("Main key contrasts: %s\n", file.path(fig_dir, "Fig_main_03_total_effect_key_contrasts.png")))
cat(sprintf("Supplement composition-adjusted OR forest: %s\n", file.path(fig_dir, "Fig_supp_01_composition_adjusted_factor_OR_forest.png")))
cat(sprintf("Gaba-source total-effect OR forest: %s\n", file.path(fig_dir, "Fig_extra_05_gaba_source_total_effect_factor_OR_forest.png")))
cat(sprintf("Glut-source total-effect OR forest: %s\n", file.path(fig_dir, "Fig_extra_06_glut_source_total_effect_factor_OR_forest.png")))
cat(sprintf("Gaba-source composition-adjusted OR forest: %s\n", file.path(fig_dir, "Fig_extra_07_gaba_source_composition_adjusted_factor_OR_forest.png")))
cat(sprintf("Glut-source composition-adjusted OR forest: %s\n", file.path(fig_dir, "Fig_extra_08_glut_source_composition_adjusted_factor_OR_forest.png")))
cat(sprintf("Model run manifest: %s\n", file.path(table_dir, "model_run_manifest.csv")))
cat(sprintf("All input rows audit: %s\n", file.path(table_dir, "all_input_rows_no_type_filter.csv")))
cat("Primary fitted models use directed neuron-neuron rows to estimate factor effects on exact containment.\n")
cat("Supplementary focal models use neuron-neuron rows only, with target_class3 restricted to Gaba or Glut.\n")
cat(sprintf("Gaba-focal model input: %s\n", file.path(table_dir, "gaba_focal_model_input.csv")))
cat(sprintf("Glut-focal model input: %s\n", file.path(table_dir, "glut_focal_model_input.csv")))
cat(sprintf("Gaba-focal target class counts: %s\n", file.path(table_dir, "gaba_focal_target_class3_counts.csv")))
cat(sprintf("Glut-focal target class counts: %s\n", file.path(table_dir, "glut_focal_target_class3_counts.csv")))
cat(sprintf("Gaba-focal adjusted target probabilities: %s\n", file.path(table_dir, "gaba_focal_target_class3_adjusted_probabilities.csv")))
cat(sprintf("Glut-focal adjusted target probabilities: %s\n", file.path(table_dir, "glut_focal_target_class3_adjusted_probabilities.csv")))
cat(sprintf("Gaba-focal target class pairwise contrasts: %s\n", file.path(table_dir, "gaba_focal_target_class3_pairwise_contrasts.csv")))
cat(sprintf("Glut-focal target class pairwise contrasts: %s\n", file.path(table_dir, "glut_focal_target_class3_pairwise_contrasts.csv")))
cat(sprintf("Final observed-rate figure: %s\n", file.path(fig_dir, "Fig_final_01_observed_exact_rate_by_target.png")))
cat(sprintf("Final adjusted-probability figure: %s\n", file.path(fig_dir, "Fig_final_02_adjusted_probability_by_target.png")))
cat(sprintf("Final target-class OR forest: %s\n", file.path(fig_dir, "Fig_final_03_glut_vs_gaba_target_OR_forest.png")))
cat(sprintf("Gaba-focal coefficient table: %s\n", file.path(table_dir, "gaba_focal_exact_logistic_coefficients.csv")))
cat(sprintf("Glut-focal coefficient table: %s\n", file.path(table_dir, "glut_focal_exact_logistic_coefficients.csv")))


