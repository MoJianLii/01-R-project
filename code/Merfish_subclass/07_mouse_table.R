rm(list = ls()); gc()

suppressPackageStartupMessages({
  library(stringr)
  library(data.table)
  library(doParallel)
  library(foreach)
  library(parallel)
})

info <- function(...) cat(sprintf("[INFO %s] ", format(Sys.time(), "%F %T")),
                          sprintf(...), "\n")

Sys.setenv(
  OPENBLAS_NUM_THREADS = "1",
  OMP_NUM_THREADS      = "1",
  MKL_NUM_THREADS      = "1",
  GOTO_NUM_THREADS     = "1"
)
data.table::setDTthreads(1L)

setwd('E:/zaw/2511/mouseMerfish_zhuang_subclass')
inp_cluster_file <- "./mouse_subclass_cluster_total.txt"
loc_dir         <- "./merge_region"


force_id_cols_char <- function(df) {
  id_cols <- names(df)[grepl("cell_label$|cell_id$|cell_ids$", names(df))]
  for (cc in id_cols) {
    df[[cc]] <- as.character(df[[cc]])
  }
  df
}

stopifnot(file.exists(inp_cluster_file))
mouse_subclass_cluster_total <- fread(
  inp_cluster_file,
  sep         = "\t",
  header      = TRUE,
  data.table  = FALSE,
  check.names = FALSE
)

mouse_subclass_cluster_total <- force_id_cols_char(mouse_subclass_cluster_total)

num_cols <- c("total_cell_num",
              "Glut_Neruon_cell_ids_num",
              "GABA_Neruon_cell_ids_num")
for (nm in num_cols) {
  if (nm %in% names(mouse_subclass_cluster_total)) {
    mouse_subclass_cluster_total[[nm]] <- as.numeric(mouse_subclass_cluster_total[[nm]])
  }
}

info("read cluster_total: rows=%d, cols=%d",
     nrow(mouse_subclass_cluster_total), ncol(mouse_subclass_cluster_total))

need_cols_cluster <- c("slide","layer","region","cell_Neuron_type","subclass",
                       "total_cell_num","Glut_Neruon_cell_ids_num","GABA_Neruon_cell_ids_num",
                       "Glut_Neruon_cell_ids","GABA_Neruon_cell_ids","Non_Neruon_cell_ids",
                       "merge_regions")
miss_cluster <- setdiff(need_cols_cluster, names(mouse_subclass_cluster_total))
if (length(miss_cluster)) {
  stop(sprintf("missing columns in mouse_subclass_cluster_total: %s",
               paste(miss_cluster, collapse = ",")))
}

cauchy_combination_test <- function(pvals, weights = NULL) {
  pvals <- as.numeric(pvals)
  pvals[pvals < 1e-15]     <- 1e-15
  pvals[pvals > 1 - 1e-15] <- 1 - 1e-15
  if (is.null(weights)) weights <- rep(1 / length(pvals), length(pvals))
  t_stat <- sum(weights * tan((0.5 - pvals) * pi))
  0.5 - atan(t_stat) / pi
}

stopifnot(dir.exists(loc_dir))
file_names <- list.files(loc_dir,
                         pattern = "merged_rgions_loc\\.txt$",
                         full.names = FALSE)
info("found loc files: %d", length(file_names))

temp.total <- data.table()

if (length(file_names) > 0) {
  workers <- max(1L, min(length(file_names), detectCores() - 1L))
  info("ACAT workers: %d", workers)
  cl <- makeCluster(workers)
  registerDoParallel(cl)
  
  temp.list <- foreach(file_name = file_names,
                       .packages = c("data.table", "stringr"),
                       .export   = c("cauchy_combination_test")) %dopar% {
                         fp <- file.path(loc_dir, file_name)
                         dt <- tryCatch(
                           fread(fp, sep = "\t", header = TRUE,
                                 data.table = FALSE, check.names = FALSE),
                           error = function(e) NULL
                         )
                         if (is.null(dt) || !nrow(dt)) {
                           return(data.frame(
                             slide = character(0),
                             layer = character(0),
                             subclass = character(0),
                             merge_regions = character(0),
                             cauchy_combination_p = numeric(0),
                             stringsAsFactors = FALSE
                           ))
                         }
                         
                         req_cols <- c("merge_regions","layer","subclass","p")
                         if (!all(req_cols %in% names(dt))) {
                           return(data.frame(
                             slide = character(0),
                             layer = character(0),
                             subclass = character(0),
                             merge_regions = character(0),
                             cauchy_combination_p = numeric(0),
                             stringsAsFactors = FALSE
                           ))
                         }
                         
                         slide_id <- str_split(file_name, "[_]", simplify = TRUE)[,1]
                         lyr      <- unique(dt$layer)
                         cls      <- unique(dt$subclass)
                         if (length(lyr) != 1L) lyr <- lyr[1]
                         if (length(cls) != 1L) cls <- cls[1]
                         
                         regs <- sort(unique(dt$merge_regions))
                         if (!length(regs)) {
                           return(data.frame(
                             slide = character(0),
                             layer = character(0),
                             subclass = character(0),
                             merge_regions = character(0),
                             cauchy_combination_p = numeric(0),
                             stringsAsFactors = FALSE
                           ))
                         }
                         
                         out_list <- vector("list", length(regs))
                         k <- 1L
                         for (rg in regs) {
                           sub <- dt[dt$merge_regions %in% rg, , drop = FALSE]
                           if (!nrow(sub)) next
                           pc <- cauchy_combination_test(sub$p)
                           out_list[[k]] <- data.frame(
                             slide = slide_id,
                             layer = lyr,
                             subclass = cls,
                             merge_regions = rg,
                             cauchy_combination_p = pc,
                             stringsAsFactors = FALSE
                           )
                           k <- k + 1L
                         }
                         if (k == 1L) {
                           data.frame(
                             slide = character(0),
                             layer = character(0),
                             subclass = character(0),
                             merge_regions = character(0),
                             cauchy_combination_p = numeric(0),
                             stringsAsFactors = FALSE
                           )
                         } else {
                           do.call(rbind, out_list[1:(k-1)])
                         }
                       }
  
  stopCluster(cl)
  
  temp.total <- if (length(temp.list)) {
    rbindlist(temp.list, use.names = TRUE, fill = TRUE)
  } else data.table()
}

info("temp.total rows=%d", nrow(temp.total))

if (nrow(temp.total)) {
  temp.total$label <- paste(
    temp.total$slide,
    temp.total$layer,
    temp.total$subclass,
    temp.total$merge_regions, sep = "_"
  )
} else {
  temp.total$label <- character(0)
}

mouse_subclass_cluster_total$label <- paste(
  mouse_subclass_cluster_total$slide,
  mouse_subclass_cluster_total$layer,
  mouse_subclass_cluster_total$subclass,
  mouse_subclass_cluster_total$merge_regions, sep = "_"
)

lab_int <- length(intersect(unique(temp.total$label),
                            unique(mouse_subclass_cluster_total$label)))
info("labels: temp.total=%d, cluster_total=%d, intersect=%d",
     length(unique(temp.total$label)),
     length(unique(mouse_subclass_cluster_total$label)),
     lab_int)

keep_cols <- c("label","cauchy_combination_p")
mouse_subclass_cluster_total_cauchy_combination_test <-
  merge(mouse_subclass_cluster_total,
        temp.total[, ..keep_cols],
        by = "label", all.x = FALSE, all.y = FALSE)

info("merged rows=%d",
     nrow(mouse_subclass_cluster_total_cauchy_combination_test))

dir.create("./mouse_table", showWarnings = FALSE)
dir.create("./mouse_table/table1", showWarnings = FALSE)

tab1_cols <- c("label","slide","layer","region","cell_Neuron_type","subclass",
               "total_cell_num","Glut_Neruon_cell_ids_num",
               "GABA_Neruon_cell_ids_num","cauchy_combination_p")
Table_1_total_cell <- mouse_subclass_cluster_total_cauchy_combination_test[, tab1_cols, drop = FALSE]

Table_1_total_cell$E_I_Ratio <- with(
  Table_1_total_cell,
  Glut_Neruon_cell_ids_num / ifelse(GABA_Neruon_cell_ids_num == 0, NA, GABA_Neruon_cell_ids_num)
)

fwrite(Table_1_total_cell, "./mouse_table/table1/Table_1_total_cell.txt",
       sep = "\t", quote = FALSE, row.names = FALSE)
info("wrote Table_1_total_cell: rows=%d", nrow(Table_1_total_cell))

mouse_cluster_number <- as.data.frame(
  table(mouse_subclass_cluster_total$subclass),
  stringsAsFactors = FALSE
)
colnames(mouse_cluster_number) <- c("subclass","cluster_number")
mouse_cluster_number <- mouse_cluster_number[
  order(mouse_cluster_number$cluster_number, decreasing = TRUE), ]
top_subclasses <- head(as.character(mouse_cluster_number$subclass), 8)

for (cl in top_subclasses) {
  sub1 <- Table_1_total_cell[Table_1_total_cell$subclass %in% cl, , drop = FALSE]
  outp <- sprintf("./mouse_table/table1/Table_1_%s_cell.txt", cl)
  fwrite(sub1, outp, sep = "\t", quote = FALSE, row.names = FALSE)
  info("wrote %s (rows=%d)", outp, nrow(sub1))
}

pairship_file <- "./mouse_cluster_pairship_new_total.txt"
stopifnot(file.exists(pairship_file))
mouse_cluster_pairship_new_total <- fread(
  pairship_file,
  sep         = "\t",
  header      = TRUE,
  data.table  = FALSE,
  check.names = FALSE
)

mouse_cluster_pairship_new_total <- force_id_cols_char(mouse_cluster_pairship_new_total)

mouse_cluster_pairship_new_total <- mouse_cluster_pairship_new_total[
  mouse_cluster_pairship_new_total$nearby_dist == 0, , drop = FALSE
]
info("pairship rows (nearby_dist==0) = %d",
     nrow(mouse_cluster_pairship_new_total))

mouse_cluster_pairship_new_total$label.1 <- paste(
  mouse_cluster_pairship_new_total$slide,
  mouse_cluster_pairship_new_total$cluster.1, sep = "_"
)
mouse_cluster_pairship_new_total$label.2 <- paste(
  mouse_cluster_pairship_new_total$slide,
  mouse_cluster_pairship_new_total$cluster.2, sep = "_"
)
mouse_cluster_pairship_new_total <- mouse_cluster_pairship_new_total[, c("label.1","label.2"), drop = FALSE]

mouse_subclass_cluster_total_cauchy_combination_test$label <- paste(
  mouse_subclass_cluster_total_cauchy_combination_test$slide,
  mouse_subclass_cluster_total_cauchy_combination_test$layer,
  mouse_subclass_cluster_total_cauchy_combination_test$subclass,
  mouse_subclass_cluster_total_cauchy_combination_test$merge_regions, sep = "_"
)
mouse_subclass_cluster_total$label <- paste(
  mouse_subclass_cluster_total$slide,
  mouse_subclass_cluster_total$layer,
  mouse_subclass_cluster_total$subclass,
  mouse_subclass_cluster_total$merge_regions, sep = "_"
)

dir.create("./mouse_table/table2", showWarnings = FALSE)

n_pairs <- nrow(mouse_cluster_pairship_new_total)
info("start computing pair overlaps for %d pairs (parallel)", n_pairs)

if (n_pairs > 0) {
  workers2 <- max(1L, min(n_pairs, detectCores() - 1L))
  info("pair overlap workers: %d", workers2)
  cl2 <- makeCluster(workers2)
  registerDoParallel(cl2)
  
  pair_overlap_dt <- foreach(i = seq_len(n_pairs),
                             .packages = c("stringr","data.table"),
                             .export   = c("mouse_subclass_cluster_total")) %dopar% {
                               
                               pr <- mouse_cluster_pairship_new_total[i, ]
                               xlab <- pr$label.1
                               ylab <- pr$label.2
                               
                               xrow <- mouse_subclass_cluster_total[
                                 mouse_subclass_cluster_total$label %in% xlab, , drop = FALSE]
                               yrow <- mouse_subclass_cluster_total[
                                 mouse_subclass_cluster_total$label %in% ylab, , drop = FALSE]
                               
                               x_cells <- Reduce(
                                 union,
                                 list(
                                   unique(as.character(str_split(xrow$Glut_Neruon_cell_ids, "[,]", simplify = TRUE))),
                                   unique(as.character(str_split(xrow$GABA_Neruon_cell_ids,  "[,]", simplify = TRUE))),
                                   unique(as.character(str_split(xrow$Non_Neruon_cell_ids,   "[,]", simplify = TRUE)))
                                 )
                               )
                               x_cells <- x_cells[x_cells != "" & !is.na(x_cells)]
                               
                               y_cells <- Reduce(
                                 union,
                                 list(
                                   unique(as.character(str_split(yrow$Glut_Neruon_cell_ids, "[,]", simplify = TRUE))),
                                   unique(as.character(str_split(yrow$GABA_Neruon_cell_ids,  "[,]", simplify = TRUE))),
                                   unique(as.character(str_split(yrow$Non_Neruon_cell_ids,   "[,]", simplify = TRUE)))
                                 )
                               )
                               y_cells <- y_cells[y_cells != "" & !is.na(y_cells)]
                               
                               overlap <- length(intersect(x_cells, y_cells))
                               un      <- length(union(x_cells, y_cells))
                               
                               data.table(
                                 label.1      = xlab,
                                 label.2      = ylab,
                                 overlap_cell = overlap,
                                 union_cell   = un
                               )
                             }
  
  stopCluster(cl2)
  
  pair_overlap <- rbindlist(pair_overlap_dt, use.names = TRUE, fill = TRUE)
} else {
  pair_overlap <- data.table()
}

info("pair_overlap rows=%d", nrow(pair_overlap))

Table_1_min <- Table_1_total_cell

colnames(pair_overlap)[1] <- "label"
Table_2_total_cell <- merge(Table_1_min, pair_overlap, by = "label")
colnames(Table_2_total_cell)[1:11] <- paste0("cluster.1_", colnames(Table_2_total_cell)[1:11])

colnames(Table_2_total_cell)[12] <- "label"
Table_2_total_cell <- merge(Table_2_total_cell, Table_1_min, by = "label")
colnames(Table_2_total_cell)[c(1, 15:24)] <- paste0("cluster.2_", colnames(Table_2_total_cell)[c(1, 15:24)])

Table_2_total_cell <- Table_2_total_cell[, c(2:12, 1, 15:24, 13:14)]
Table_2_total_cell$jaccard <- Table_2_total_cell$overlap_cell /
  ifelse(Table_2_total_cell$union_cell == 0, NA, Table_2_total_cell$union_cell)
Table_2_total_cell$cluster.1.overlap.percent <-
  Table_2_total_cell$overlap_cell /
  ifelse(Table_2_total_cell$cluster.1_total_cell_num == 0, NA,
         Table_2_total_cell$cluster.1_total_cell_num)
Table_2_total_cell$cluster.2.overlap.percent <-
  Table_2_total_cell$overlap_cell /
  ifelse(Table_2_total_cell$cluster.2_total_cell_num == 0, NA,
         Table_2_total_cell$cluster.2_total_cell_num)

fwrite(Table_2_total_cell, "./mouse_table/table2/Table_2_total_cell.txt",
       sep = "\t", quote = FALSE, row.names = FALSE)
info("wrote Table_2_total_cell: rows=%d", nrow(Table_2_total_cell))

