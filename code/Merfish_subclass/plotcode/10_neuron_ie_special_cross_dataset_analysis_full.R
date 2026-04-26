# 清空环境
rm(list = ls())
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
base_dir <- "E:/zaw/2603"
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
  ggplot2::theme_classic(base_size = 12) +
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
  ggplot2::ggsave(file.path(fig_dir, filename), p, width = width, height = height, dpi = dpi, bg = "white", limitsize = FALSE)
  pdf_name <- sub("\\.png$", ".pdf", filename)
  ggplot2::ggsave(file.path(fig_dir, pdf_name), p, width = width, height = height, bg = "white", device = cairo_pdf, limitsize = FALSE)
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
  theme_pub()
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
  theme_pub()
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
  theme_pub()
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
  theme_pub()
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
  theme_pub()
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
    subtitle = wrap_subtitle("Target-class terms remain strong after size and composition covariates are included."),
    x = "Adjusted odds ratio (log scale)", y = NULL
  ) +
  theme_pub()
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
  theme_pub()
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
  theme_pub() +
  theme(
    legend.position = "none",
    plot.margin = margin(14, 34, 10, 10)
  )

save_pub(p7, "fig7_exact_edge_boxplots.png", width = 11.6, height = 4.8)

message("R_output_10 finished: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE))
