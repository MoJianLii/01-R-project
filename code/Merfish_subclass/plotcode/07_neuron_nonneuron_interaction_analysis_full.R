# 清空环境
rm(list = ls())
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
base_dir <- "E:/zaw/2603"
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
  ggsave(filename = file.path(fig_dir, paste0(filename, ".png")), plot = plot_obj,
         width = width, height = height, dpi = dpi, bg = "white", limitsize = FALSE)
  ggsave(filename = file.path(fig_dir, paste0(filename, ".pdf")), plot = plot_obj,
         width = width, height = height, bg = "white", device = cairo_pdf, limitsize = FALSE)
}

plot_theme_safe <- ggplot2::theme_classic(base_size = 12) +
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
