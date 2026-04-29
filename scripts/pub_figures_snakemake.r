#!/usr/bin/env Rscript
# pub_figures_snakemake.r — Publication figures for the gpmapping Snakemake pipeline.
#
# Usage (from repo root):
#   Rscript scripts/pub_figures_snakemake.r \
#     --base-dir . \
#     --output-dir plots \
#     --yeast-seed 6174 \
#     --yeast-prefix yeast \
#     --yeast-plink-threshold 1.00e-05

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(stringr)
  library(optparse)
  library(data.table)
  library(patchwork)
  library(scales)
  library(cli)
  library(systemfonts)
  library(ragg)
  library(svglite)
  library(ggpubr)
  library(rstatix)
  library(ggstar)
  library(ggsurvfit)
})

# Try to load arcadiathemeR; fall back to theme_bw if unavailable or broken.
ARCADIA_AVAILABLE <- tryCatch(
  {
    suppressPackageStartupMessages(library(arcadiathemeR))
    extrafont::loadfonts(quiet = TRUE)
    TRUE
  },
  error = function(e) {
    cli_alert_warning(
      "arcadiathemeR unavailable ({conditionMessage(e)}); using theme_bw fallback."
    )
    FALSE
  }
)

if (ARCADIA_AVAILABLE) {
  arcadia_palette <- show_arcadia_palettes()
} else {
  gradient_fill_arcadia <- function(
    palette_name = "magma",
    reverse = FALSE,
    limits = NULL,
    trans = "identity",
    breaks = waiver(),
    labels = waiver(),
    ...
  ) {
    palette_map <- c(
      magma = "magma",
      purples = "magma",
      greens = "viridis",
      blues = "mako",
      yellows = "rocket"
    )
    option <- unname(palette_map[palette_name])
    if (is.na(option)) {
      option <- "magma"
    }
    scale_fill_viridis_c(
      option = option,
      direction = if (reverse) -1 else 1,
      limits = limits,
      trans = trans,
      breaks = breaks,
      labels = labels,
      ...
    )
  }
  arcadia_palette <- list(
    yellow_shades = c("#F5E4BE", "#FFD364", "#F7B846", "#D68D22", "#A85E28"),
    purple_shades = c("#DCDFEF", "#BABEE0", "#7A77AB", "#54448C", "#341E60"),
    blue_shades = c("#C6E7F4", "#73B5E3", "#5088C5", "#2B66A2", "#094468"),
    pink_shades = c("#FFE3D4", "#F8C5C1", "#F898AE", "#E2718F", "#C04C70"),
    primary_ordered = c(
      "#5088C5",
      "#F28360",
      "#F7B846",
      "#97CD78",
      "#7A77AB",
      "#F898AE",
      "#3B9886",
      "#C85152",
      "#73B5E3",
      "#FFB984",
      "#F5E4BE",
      "#BABEE0"
    ),
    teal_shades = c("#C3E2DB", "#6FBCAD", "#3B9886", "#2A6B5E", "#09473E"),
    orange_shades = c("#FFCFAF", "#FFB883", "#F28360", "#C85152", "#9E3F41"),
    warm_gray_shades = c("#EDE6DA", "#DBD1C3", "#B9AFA7", "#8F8885", "#635C5A")
  )
}

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
option_list <- list(
  make_option(
    "--base-dir",
    type = "character",
    default = ".",
    help = "Base directory containing all pipeline outputs [default: ..]"
  ),
  make_option(
    "--output-dir",
    type = "character",
    default = "plots",
    help = "Directory to save plots [default: plots]"
  ),
  make_option(
    "--intermediates-dir",
    type = "character",
    default = "figure_intermediates",
    help = "Subdirectory (under base-dir) for cached feathers [default: figure_intermediates]"
  ),
  make_option(
    "--yeast-seed",
    type = "character",
    default = "6174",
    help = "Seed string used in Snakemake for yeast [default: 6174]"
  ),
  make_option(
    "--yeast-prefix",
    type = "character",
    default = "yeast",
    help = "File prefix for yeast outputs [default: yeast]"
  ),
  make_option(
    "--yeast-plink-threshold",
    type = "character",
    default = "1.00e-05",
    help = "Plink p-value threshold string for yeast [default: 1.00e-05]"
  ),
  make_option(
    "--human-seed",
    type = "character",
    default = "",
    help = "Seed for human outputs (leave empty to skip human) [default: '
']"
  ),
  make_option(
    "--human-prefix",
    type = "character",
    default = "human",
    help = "File prefix for human outputs [default: human]"
  ),
  make_option(
    "--human-plink-threshold",
    type = "character",
    default = "4.40e-06",
    help = "Plink p-value threshold string for human [default: 4.40e-06]"
  ),
  make_option(
    "--lars-maxiter",
    type = "integer",
    default = 1000L,
    help = "LARS max iterations (used in directory/file names) [default: 1000]"
  ),
  make_option(
    "--recalc-cors",
    action = "store_true",
    default = FALSE,
    help = "Recalculate all per-trait correlations from raw files"
  ),
  make_option(
    "--collect-cors-params",
    action = "store_true",
    default = FALSE,
    help = "Re-collect per-SNP beta matrices from raw files"
  ),
  make_option(
    "--recalc-wider-longer",
    action = "store_true",
    default = FALSE,
    help = "Recompute wide/long beta pivot tables"
  ),
  make_option(
    "--remake-maf",
    action = "store_true",
    default = FALSE,
    help = "Reload minor allele frequencies"
  ),
  make_option(
    "--remake-cumulative",
    action = "store_true",
    default = FALSE,
    help = "Recompute cumulative TP/FP tables for ROC"
  ),
  make_option(
    "--remake-roc",
    action = "store_true",
    default = FALSE,
    help = "Recompute exact-match ROC curves"
  ),
  make_option(
    "--remake-roc-approx",
    action = "store_true",
    default = FALSE,
    help = "Recompute approximate ROC curves"
  ),
  make_option(
    "--redo-snp-pairs",
    action = "store_true",
    default = FALSE,
    help = "Recompute within-chromosome SNP pair distances"
  ),
  make_option(
    "--redo-distances",
    action = "store_true",
    default = FALSE,
    help = "Recompute nearest true-SNP distances per method"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

BASE_DIR <- opt$`base-dir`
OUTPUT_DIR <- file.path(opt$`base-dir`, opt$`output-dir`)
YEAST_SEED <- opt$`yeast-seed`
YEAST_PREFIX <- opt$`yeast-prefix`
YEAST_P_THR <- opt$`yeast-plink-threshold`
HUMAN_SEED <- opt$`human-seed`
HUMAN_PREFIX <- opt$`human-prefix`
HUMAN_P_THR <- opt$`human-plink-threshold`
LARS_MAXITER <- opt$`lars-maxiter`
LARS_TAG <- paste0("lars_maxiter", LARS_MAXITER)

# Intermediates are stored in a dataset-specific subdirectory so that runs
# with different prefixes or seeds never overwrite each other's cached
# files.
dataset_tag <- paste0(YEAST_PREFIX, "_", YEAST_SEED)
if (nchar(HUMAN_SEED) > 0) {
  dataset_tag <- paste0(dataset_tag, "_", HUMAN_PREFIX, "_", HUMAN_SEED)
}
INTERMEDIATES_DIR <- file.path(BASE_DIR, opt$`intermediates-dir`, dataset_tag)

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(INTERMEDIATES_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Intermediate file paths — all seeds embedded in filenames
# ---------------------------------------------------------------------------
ROC_APPROX_WINDOW <- 250L
int_path <- function(...) file.path(INTERMEDIATES_DIR, paste0(..., ".feather"))

INT_PRED_CORS <- int_path("all_methods_prediction_correlations_", dataset_tag)
INT_EFF_CORS <- int_path("all_methods_effects_correlations_", dataset_tag)
INT_MAX_CORS <- int_path(
  "all_methods_max_prediction_correlations_",
  dataset_tag
)
INT_PARAMS <- int_path("all_methods_parameters_", dataset_tag)
INT_BETAS_YEAST <- int_path("combined_all_betas_yeast_", YEAST_SEED)
INT_BETAS_HUMAN <- int_path("combined_all_betas_human_", HUMAN_SEED)
INT_LONGER_YEAST <- int_path("yeast_littlelonger_with_fullinfo_", YEAST_SEED)
INT_LONGER_HUMAN <- int_path("human_littlelonger_with_fullinfo_", HUMAN_SEED)
INT_CUMUL_YEAST <- int_path("yeast_cumulative_", YEAST_SEED)
INT_ROC_YEAST <- int_path("yeast_roc_", YEAST_SEED)
INT_ROC_APPROX_YEAST <- int_path(
  "yeast_roc_approx_",
  ROC_APPROX_WINDOW,
  "_",
  YEAST_SEED
)
# The following files are not provided due to potential privacy issues, but are
# included here as comments in case someone wants to repeat this analysis
# in the future.
# INT_CUMUL_HUMAN <- int_path("human_cumulative_", HUMAN_SEED)
# INT_ROC_HUMAN <- int_path("human_roc_", HUMAN_SEED)
# INT_ROC_APPROX_HUMAN <- int_path(
#   "human_roc_approx_",
#   ROC_APPROX_WINDOW,
#   "_",
#   HUMAN_SEED
# )

INT_SNP_PAIRS_YEAST <- int_path("yeast_snp_pairs_", YEAST_SEED)
INT_SNP_PAIRS_HUMAN <- int_path("human_snp_pairs_", HUMAN_SEED)
INT_DIST_YEAST <- int_path("yeast_true_avg_distance_", YEAST_SEED)
INT_DIST_HUMAN <- int_path("human_true_avg_distance_", HUMAN_SEED)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
SKLEARN_METHODS <- c("ridge", "elasticnet", "lars", "lasso")
PYTORCH_METHODS <- c("pytorch_ridge")
SKLEARN_SPLIT <- "test"
PYTORCH_SPLIT <- "train"
H2_prop <- c("_H2_0.25", "_H2_0.5", "_H2_0.75", "_H2_0.99", "_H2_1")


# theme_pub() wraps theme_arcadia when available, otherwise falls back to
# theme_bw. All plot code uses theme_pub() rather than calling theme_pub()
# directly, so the script works on machines without arcadiathemeR installed.
theme_pub <- function(x_axis_type = NULL) {
  if (ARCADIA_AVAILABLE) {
    arcadiathemeR::theme_arcadia(x_axis_type = x_axis_type)
  } else {
    theme_bw(base_size = 16)
  }
}

METHOD_LABELS <- c(
  plink.50 = "PLINK",
  plink_clumped.50 = "Clumped\nPLINK",
  ridge.50 = "Ridge",
  lasso.50 = "LASSO",
  elasticnet.50 = "Elastic\nNet",
  lars_maxiter1000.50 = "LARS",
  lars_maxiter10000.50 = "LARS (10000)",
  pytorch_ridge.50 = "Approx.\nRidge",
  plink.100 = "PLINK",
  plink_clumped.100 = "Clumped\nPLINK",
  ridge.100 = "Ridge",
  lasso.100 = "LASSO",
  elasticnet.100 = "Elastic\nNet",
  lars_maxiter1000.100 = "LARS",
  lars_maxiter10000.100 = "LARS (10000)",
  pytorch_ridge.100 = "Approx.\nRidge",
  plink = "PLINK",
  plink_clumped = "Clumped\nPLINK",
  ridge = "Ridge",
  lasso = "LASSO",
  elasticnet = "Elastic\nNet",
  lars_maxiter1000 = "LARS",
  lars_maxiter10000 = "LARS (10000)",
  pytorch_ridge = "Approx.\nRidge",
  plink_train = "PLINK",
  plink_clumped_train = "Clumped\nPLINK",
  plink_cov6_train = "PLINK",
  plink_clumped_cov6_train = "Clumped\nPLINK",
  ridge_test = "Ridge",
  lasso_test = "LASSO",
  elasticnet_test = "Elastic\nNet",
  lars_maxiter1000_test = "LARS",
  lars_maxiter10000_test = "LARS (10000)",
  pytorch_ridge_train = "Approx.\nRidge"
)


METHOD_COLOURS <- c(
  elasticnet = arcadia_palette$yellow_shades[3],
  elasticnet_test = arcadia_palette$yellow_shades[3],
  elasticnet.small = arcadia_palette$yellow_shades[2],
  elasticnet.big = arcadia_palette$yellow_shades[3],
  elasticnet_test.TRUE = arcadia_palette$yellow_shades[4],
  elasticnet_test.FALSE = arcadia_palette$yellow_shades[3],
  lasso = arcadia_palette$purple_shades[3],
  lasso_test = arcadia_palette$purple_shades[3],
  lasso.small = arcadia_palette$purple_shades[2],
  lasso.big = arcadia_palette$purple_shades[3],
  lasso_test.TRUE = arcadia_palette$purple_shades[4],
  lasso_test.FALSE = arcadia_palette$purple_shades[3],
  lars_maxiter1000 = arcadia_palette$blue_shades[2],
  lars_maxiter1000_test = arcadia_palette$blue_shades[2],
  lars_maxiter1000.small = arcadia_palette$blue_shades[1],
  lars_maxiter1000.big = arcadia_palette$blue_shades[2],
  lars_maxiter1000_test.TRUE = arcadia_palette$blue_shades[3],
  lars_maxiter1000_test.FALSE = arcadia_palette$blue_shades[2],
  pytorch_ridge = arcadia_palette$primary_ordered[4],
  pytorch_ridge_train = arcadia_palette$primary_ordered[4],
  pytorch_ridge.small = "#97CD78",
  pytorch_ridge.big = arcadia_palette$primary_ordered[4],
  pytorch_ridge_train.TRUE = "#72ad50",
  pytorch_ridge_train.FALSE = arcadia_palette$primary_ordered[4],
  ridge = arcadia_palette$teal_shades[3],
  ridge_test = arcadia_palette$teal_shades[3],
  ridge.small = arcadia_palette$teal_shades[2],
  ridge.big = arcadia_palette$teal_shades[3],
  ridge_test.TRUE = arcadia_palette$teal_shades[4],
  ridge_test.FALSE = arcadia_palette$teal_shades[3],
  plink_clumped = "#F8C5C1",
  plink_clumped_train = "#F8C5C1",
  plink_clumped_cov6_train = "#F8C5C1",
  plink_clumped_cov6.small = "#FFE3D4",
  plink_clumped.small = "#FFE3D4",
  plink_clumped_cov6.big = "#F8C5C1",
  plink_clumped.big = "#F8C5C1",
  plink_clumped_train.TRUE = "#FFE3D4",
  plink_clumped_train.FALSE = "#F8C5C1",
  plink_clumped_cov6_train.TRUE = "#FFE3D4",
  plink_clumped_cov6_train.FALSE = "#F8C5C1",
  plink = arcadia_palette$orange_shades[3],
  plink_train = arcadia_palette$orange_shades[3],
  plink_cov6_train = arcadia_palette$orange_shades[3],
  plink_cov6.small = arcadia_palette$orange_shades[2],
  plink.small = arcadia_palette$orange_shades[2],
  plink_cov6.big = arcadia_palette$orange_shades[3],
  plink.big = arcadia_palette$orange_shades[3],
  plink_train.TRUE = arcadia_palette$orange_shades[4],
  plink_train.FALSE = arcadia_palette$orange_shades[3],
  plink_cov6_train.TRUE = arcadia_palette$orange_shades[4],
  plink_cov6_train.FALSE = arcadia_palette$orange_shades[3],
  `50` = arcadia_palette$warm_gray_shades[5],
  `100` = arcadia_palette$warm_gray_shades[5],
  `500` = arcadia_palette$warm_gray_shades[5],
  small = arcadia_palette$warm_gray_shades[4],
  big = arcadia_palette$warm_gray_shades[5]
)

METHOD_COLOURS2 <- c(
  unlist(arcadia_palette$orange_shades),
  unlist(arcadia_palette$yellow_shades),
  unlist(arcadia_palette$teal_shades),
  c("#bedead", "#97CD78", "#6da64c", "#3f7322", "#2e5c14"),
  unlist(arcadia_palette$blue_shades),
  unlist(arcadia_palette$purple_shades),
  unlist(arcadia_palette$pink_shades),
  unlist(arcadia_palette$warm_gray_shades)
)

names(METHOD_COLOURS2) <- c(
  paste0("elasticnet_test.", 1:5),
  paste0("lasso_test.", 1:5),
  paste0("lars_maxiter1000_test.", 1:5),
  paste0("pytorch_ridge_train.", 1:5),
  paste0("ridge_test.", 1:5),
  paste0("plink_clumped_train.", 1:5),
  paste0("plink_train.", 1:5),
  c(
    '[0,0.0155]',
    '(0.0155,0.0351]',
    '(0.0351,0.0883]',
    '(0.0883,0.249]',
    '(0.249,0.5]'
  )
)


METHOD_COLOURS <- c(METHOD_COLOURS, METHOD_COLOURS2)

METHOD_ORDER_HUMAN <- c(
  "elasticnet",
  "lasso",
  LARS_TAG,
  "pytorch_ridge",
  "ridge",
  "plink_clumped",
  "plink"
)
METHOD_ORDER_YEAST <- setdiff(METHOD_ORDER_HUMAN, "elasticnet")


SPARSE_METHODS <- c("elasticnet_test", "lasso_test", paste0(LARS_TAG, "_test"))
PLINK_METHODS <- c(
  "plink_clumped_train",
  "plink_train",
  "plink_clumped_cov6_train",
  "plink_cov6_train"
)
RIDGE_METHODS <- c("pytorch_ridge_train", "ridge_test")
r2_group_vars <- c(
  "overlap",
  "target_vavg",
  "numQTL",
  "numChr",
  "target_H2",
  "method_append",
  "species",
  "target_h2"
)

# ---------------------------------------------------------------------------
# SEEDS configuration
# ---------------------------------------------------------------------------
SEEDS <- list()
SEEDS[[YEAST_SEED]] <- list(
  prefix = YEAST_PREFIX,
  label = "Yeast",
  plink_threshold = YEAST_P_THR,
  ID_col_select = 1
)
if (nchar(HUMAN_SEED) > 0) {
  SEEDS[[HUMAN_SEED]] <- list(
    prefix = HUMAN_PREFIX,
    label = "Human",
    plink_threshold = HUMAN_P_THR,
    ID_col_select = 2
  )
}

# ---------------------------------------------------------------------------
# SNP info loading
# ---------------------------------------------------------------------------
yeast_snp_file <- file.path(
  BASE_DIR,
  "input_data",
  "SNP_list_pos_corrected.txt"
)
yeast_snp_list <- fread(yeast_snp_file)
yeast_snp_list <- yeast_snp_list[yeast_snp_list$SNP != "END", ]
colnames(yeast_snp_list) <- c(
  "Index",
  "SNP",
  "SNP_0index",
  "Chr",
  "POS",
  "REF",
  "ALT"
)
yeast_numeric_part <- as.numeric(gsub("[^0-9]", "", yeast_snp_list$SNP))
yeast_snp_list_sorted <- yeast_snp_list[order(yeast_numeric_part)]
yeast_chr_along <- yeast_snp_list_sorted %>%
  group_by(Chr) %>%
  summarize(min_pos = min(POS), max_pos = max(POS)) %>%
  ungroup() %>%
  mutate(chrom_sum = cumsum(max_pos))
yeast_chr_along$chrom_sum <- c(
  0,
  yeast_chr_along$chrom_sum[1:(nrow(yeast_chr_along) - 1)]
)
yeast_chr_along <- yeast_chr_along %>%
  mutate(
    min_pos_new = chrom_sum + min_pos,
    max_pos_new = chrom_sum + max_pos,
    chrom_ave = chrom_sum + (min_pos + max_pos) / 2
  )
full_snp_info_yeast <- left_join(
  yeast_snp_list_sorted,
  yeast_chr_along,
  by = "Chr"
)

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
expand_trait_names <- function(df, keep = FALSE) {
  dt <- as.data.table(df)
  if (nrow(dt) == 0) {
    dt[, `:=`(
      overlap = character(),
      target_vavg = numeric(),
      numQTL = numeric(),
      numChr = numeric(),
      target_H2 = numeric(),
      target_h2 = numeric()
    )]
    if (!keep) {
      dt[, trait := NULL]
    }
    return(as.data.frame(dt))
  }
  parts <- dt[, tstrsplit(trait, "_", fixed = TRUE)]
  dt[, `:=`(
    overlap = parts[[2]],
    target_vavg = as.numeric(parts[[4]]),
    numQTL = as.numeric(parts[[6]]),
    numChr = as.numeric(parts[[8]]),
    target_H2 = as.numeric(parts[[10]]),
    target_h2 = as.numeric(parts[[4]]) * as.numeric(parts[[10]])
  )]
  if (!keep) {
    dt[, trait := NULL]
  }
  as.data.frame(dt)
}

collect_cors <- function(
  truth_mat,
  est_mat,
  id_col,
  method_name,
  method_task,
  seed,
  split,
  covariates = NA
) {
  est_mat <- est_mat %>% tibble::column_to_rownames(id_col)
  truth_mat <- truth_mat %>% tibble::column_to_rownames(id_col)
  shared_traits <- intersect(colnames(est_mat), colnames(truth_mat))
  if (length(shared_traits) == 0) {
    warning(sprintf("[%s seed %s] No shared traits.", method_name, seed))
    return(NULL)
  }
  shared_rows <- intersect(rownames(est_mat), rownames(truth_mat))
  if (length(shared_rows) == 0) {
    warning(sprintf("[%s seed %s] No shared rows.", method_name, seed))
    return(NULL)
  }
  est_mat <- est_mat[shared_rows, shared_traits, drop = FALSE]
  truth_mat <- truth_mat[shared_rows, shared_traits, drop = FALSE]
  cors <- sapply(shared_traits, function(t) {
    tryCatch(
      cor(
        as.numeric(est_mat[, t]),
        as.numeric(truth_mat[, t]),
        use = "pairwise.complete.obs"
      ),
      error = function(e) NA_real_
    )
  })
  data.frame(
    method = method_name,
    task = method_task,
    seed = as.character(seed),
    split = split,
    trait = shared_traits,
    cor = cors,
    r2 = cors^2,
    cov = covariates
  )
}

calculate_cors_all <- function(
  truth_phenos,
  truth_geneticvalues,
  truth_addvalues,
  id_col,
  eval,
  seed,
  split
) {
  gv <- truth_geneticvalues %>% tibble::column_to_rownames(id_col)
  add <- truth_addvalues %>% tibble::column_to_rownames(id_col)
  pheno <- truth_phenos %>% tibble::column_to_rownames(id_col)
  shared_traits <- intersect(
    intersect(colnames(gv), colnames(add)),
    colnames(pheno)
  )
  gv <- gv[rownames(gv) %in% eval, ]
  add <- add[rownames(add) %in% eval, ]
  pheno <- pheno[rownames(pheno) %in% eval, ]
  shared_ind <- intersect(
    intersect(rownames(gv), rownames(pheno)),
    rownames(add)
  )
  gv <- gv[shared_ind, shared_traits]
  add <- add[shared_ind, shared_traits]
  pheno <- pheno[shared_ind, shared_traits]
  true_H2_cor <- sapply(shared_traits, function(t) {
    cor(gv[[t]], pheno[[t]], use = "pairwise.complete.obs")
  })
  true_vavg <- sapply(shared_traits, function(t) var(add[[t]]) / var(gv[[t]]))
  true_h2_cor <- sapply(shared_traits, function(t) {
    cor(add[[t]], pheno[[t]], use = "pairwise.complete.obs")
  })
  expand_trait_names(data.frame(
    trait = shared_traits,
    seed = seed,
    split = split,
    true_H2 = true_H2_cor^2,
    true_h2 = true_h2_cor^2,
    true_vavg = true_vavg,
    row.names = NULL
  ))
}

get_indiv_sets <- function(dir, file_prefix, seed, biobank = FALSE) {
  prefix <- file.path(dir, paste0(file_prefix, "_seed_", seed))
  if (!biobank) {
    test_ids <- readLines(paste0(prefix, "_test_ids.txt"))[-1]
    train_ids <- readLines(paste0(prefix, "_train_ids.txt"))[-1]
    test_ids <- as.numeric(sub("^i", "", test_ids))
    train_ids <- as.numeric(sub("^i", "", train_ids))
  } else {
    test_ids <- read.table(paste0(prefix, "_test_ids.txt"))
    colnames(test_ids) <- c("FID", "IID")
    test_ids <- test_ids$IID
    train_ids <- read.table(paste0(prefix, "_train_ids.txt"))
    colnames(train_ids) <- c("FID", "IID")
    train_ids <- train_ids$IID
  }
  list(test_ids = test_ids, train_ids = train_ids)
}

# ---------------------------------------------------------------------------
# Recalculate correlations
# ---------------------------------------------------------------------------
if (opt$`recalc-cors`) {
  cat("Recalculating correlations...\n")
  eff_cors <- list()
  pred_cors <- list()
  true_max_cors <- list()
  method_params <- list()
  i <- 1
  for (seed_chr in names(SEEDS)) {
    cfg <- SEEDS[[seed_chr]]
    file_prefix <- cfg$prefix
    p_thr <- cfg$plink_threshold
    seed <- as.numeric(seed_chr)
    is_biobank <- seed_chr == HUMAN_SEED
    cli_alert_info("Starting seed: {seed_chr} ({cfg$label})")

    indiv_sets <- get_indiv_sets(
      file.path(BASE_DIR, paste0("test_train_seed_", seed_chr)),
      cfg$prefix,
      seed_chr,
      biobank = is_biobank
    )
    cfg$test_set <- indiv_sets$test_ids
    cfg$train_set <- indiv_sets$train_ids

    true_phenos <- fread(file.path(
      BASE_DIR,
      paste0("test_train_seed_", seed_chr),
      paste0(file_prefix, "_seed_", seed_chr, "_all_phenotypes_normalized.csv")
    ))
    if (!is_biobank) {
      true_phenos$IID <- as.numeric(sub("^i", "", true_phenos$IID))
    }

    seed_int_dir <- file.path(BASE_DIR, paste0("intermediates_seed_", seed_chr))
    true_gvs_pre <- read_feather(file.path(
      seed_int_dir,
      paste0(file_prefix, "_seed_", seed_chr, "_geneticvalues.feather")
    ))
    true_add_scores_pre <- read_feather(file.path(
      seed_int_dir,
      paste0(
        file_prefix,
        "_seed_",
        seed_chr,
        "_variancecomponents_additive_scaled.feather"
      )
    ))

    for_gv_big <- list()
    for_add_big <- list()
    j <- 1
    for (H2 in H2_prop) {
      gv2 <- true_gvs_pre
      colnames(gv2) <- paste0(colnames(gv2), H2)
      add2 <- true_add_scores_pre
      colnames(add2) <- paste0(colnames(add2), H2)
      for_gv_big[[j]] <- gv2
      for_add_big[[j]] <- add2
      j <- j + 1
    }
    true_gvs <- as.data.frame(do.call(cbind, for_gv_big))
    true_add_scores <- as.data.frame(do.call(cbind, for_add_big))
    true_gvs$IID <- true_phenos$IID
    true_add_scores$IID <- true_phenos$IID

    true_effects <- read_feather(file.path(
      seed_int_dir,
      paste0(file_prefix, "_seed_", seed_chr, "_additive_effects_wide.feather")
    )) %>%
      separate_wider_delim(
        SNP,
        names = c("SNP_clean", NA),
        delim = "_",
        too_few = "align_start"
      ) %>%
      rename(SNP = SNP_clean)

    # Plink methods
    for (data_subset in c("test", "train")) {
      evaluation_set <- if (data_subset == "train") {
        cfg$test_set
      } else {
        cfg$train_set
      }
      true_max_cors[[i]] <- calculate_cors_all(
        true_phenos,
        true_gvs,
        true_add_scores,
        "IID",
        evaluation_set,
        cfg$label,
        data_subset
      )
      i <- i + 1

      cov_tags_plink <- if (is_biobank) c("", "cov6_") else c("")
      for (cov_tag in cov_tags_plink) {
        plink_dir <- file.path(
          BASE_DIR,
          paste0("aggregate_outputs_", cov_tag, data_subset, "_seed_", seed_chr)
        )
        plink_suff <- paste0(cov_tag, data_subset, "_seed_", seed_chr)
        clump_suff <- if (is_biobank) "_clumped" else "_clumped_r2_0.5"

        if (!dir.exists(plink_dir)) {
          cat(
            "DEBUG recalc-cors plink_dir not found, skipping:",
            plink_dir,
            "\n"
          )
          next
        }

        scores_file <- file.path(
          plink_dir,
          paste0("polygenic_scores_", plink_suff, "_p", p_thr, ".feather")
        )

        estimated_phenos <- read_feather(scores_file)
        if (is_biobank) {
          estimated_phenos <- estimated_phenos %>% rename(IID = individual)
        } else {
          estimated_phenos$IID <- as.numeric(sub(
            "^i",
            "",
            estimated_phenos$individual
          ))
          estimated_phenos <- estimated_phenos %>% select(-individual)
        }
        estimated_phenos <- estimated_phenos[
          match(evaluation_set, estimated_phenos$IID),
        ]

        estimated_effects <- read_feather(file.path(
          plink_dir,
          paste0("beta_matrix_filtered_", plink_suff, "_p", p_thr, ".feather")
        ))
        estimated_effects_unfiltered <- read_feather(file.path(
          plink_dir,
          paste0("beta_matrix_", plink_suff, ".feather")
        ))
        estimated_effects <- estimated_effects[
          estimated_effects$SNP != "__bias__",
        ]
        estimated_effects <- estimated_effects[
          match(true_effects$SNP, estimated_effects$SNP),
        ]
        estimated_effects_unfiltered <- estimated_effects_unfiltered[
          estimated_effects_unfiltered$SNP != "__bias__",
        ]
        estimated_effects_unfiltered <- estimated_effects_unfiltered[
          match(true_effects$SNP, estimated_effects_unfiltered$SNP),
        ]

        eff_cors[[i]] <- collect_cors(
          true_effects,
          estimated_effects,
          "SNP",
          "plink",
          "effects",
          cfg$label,
          data_subset,
          cov_tag
        )
        pred_cors[[i]] <- collect_cors(
          true_phenos,
          estimated_phenos,
          "IID",
          "plink",
          "phenotypes",
          cfg$label,
          data_subset,
          cov_tag
        )
        i <- i + 1

        eff_cors[[i]] <- collect_cors(
          true_effects,
          estimated_effects_unfiltered,
          "SNP",
          "plink",
          "effects_unfiltered",
          cfg$label,
          data_subset,
          cov_tag
        )
        i <- i + 1

        # Clumped (optional)
        clumped_scores_file <- file.path(
          plink_dir,
          paste0(
            "polygenic_scores_",
            plink_suff,
            "_p",
            p_thr,
            clump_suff,
            ".feather"
          )
        )
        clumped_betas_file <- file.path(
          plink_dir,
          paste0(
            "beta_matrix_filtered_",
            plink_suff,
            "_p",
            p_thr,
            clump_suff,
            ".feather"
          )
        )
        if (
          file.exists(clumped_scores_file) && file.exists(clumped_betas_file)
        ) {
          estimated_phenos_clumped <- read_feather(clumped_scores_file)
          estimated_effects_clumped <- read_feather(clumped_betas_file)
          if (is_biobank) {
            estimated_phenos_clumped <- estimated_phenos_clumped %>%
              rename(IID = individual)
          } else {
            estimated_phenos_clumped$IID <- as.numeric(sub(
              "^i",
              "",
              estimated_phenos_clumped$individual
            ))
            estimated_phenos_clumped <- estimated_phenos_clumped %>%
              select(-individual)
          }
          estimated_phenos_clumped <- estimated_phenos_clumped[
            match(evaluation_set, estimated_phenos_clumped$IID),
          ]
          estimated_effects_clumped <- estimated_effects_clumped[
            match(true_effects$SNP, estimated_effects_clumped$SNP),
          ]
          pred_cors[[i]] <- collect_cors(
            true_phenos,
            estimated_phenos_clumped,
            "IID",
            "plink_clumped",
            "phenotypes",
            cfg$label,
            data_subset,
            cov_tag
          )
          eff_cors[[i]] <- collect_cors(
            true_effects,
            estimated_effects_clumped,
            "SNP",
            "plink_clumped",
            "effects",
            cfg$label,
            data_subset,
            cov_tag
          )
          i <- i + 1
        }
      } # end cov_tags_plink loop
    } # end data_subset loop

    # Sklearn + pytorch methods
    sklearn_methods_run <- c(
      "lasso",
      LARS_TAG,
      "ridge",
      "elasticnet",
      "pytorch_ridge"
    )
    cov_tags_sklearn <- if (is_biobank) c("", "cov6_") else c("")
    pb_sk <- cli_progress_bar(
      paste0("Correlations: ", seed_chr, " methods"),
      total = length(sklearn_methods_run) * length(cov_tags_sklearn)
    )
    for (method in sklearn_methods_run) {
      if (method == "pytorch_ridge") {
        data_subset <- "train"
        evaluation_set <- cfg$test_set
      } else {
        data_subset <- "test"
        evaluation_set <- cfg$train_set
      }
      for (cov_tag in cov_tags_sklearn) {
        method_dir <- file.path(
          BASE_DIR,
          paste0(
            "aggregate_outputs_",
            method,
            "_",
            cov_tag,
            data_subset,
            "_seed_",
            seed_chr
          )
        )
        suff_pre <- paste0(
          method,
          "_",
          cov_tag,
          data_subset,
          "_seed_",
          seed_chr
        )
        suff <- paste0(suff_pre, ".feather")

        pred_file <- file.path(method_dir, paste0("prediction_matrix_", suff))

        if (!file.exists(pred_file)) {
          cli_progress_update(id = pb_sk)
          next
        }

        params <- fread(file.path(
          method_dir,
          paste0("params_", suff_pre, ".csv")
        ))
        estimated_phenos <- read_feather(pred_file)
        estimated_phenos <- estimated_phenos[
          match(evaluation_set, estimated_phenos$IID),
        ]
        estimated_effects <- read_feather(file.path(
          method_dir,
          paste0("coeff_matrix_", suff)
        ))
        estimated_effects <- estimated_effects[
          estimated_effects$SNP != "__bias__",
        ] %>%
          separate_wider_delim(
            SNP,
            names = c("SNP_clean", NA),
            delim = "_",
            too_few = "align_start"
          ) %>%
          rename(SNP = SNP_clean)
        estimated_effects <- estimated_effects[
          match(true_effects$SNP, estimated_effects$SNP),
        ]

        method_params[[i]] <- params
        eff_cors[[i]] <- collect_cors(
          true_effects,
          estimated_effects,
          "SNP",
          method,
          "effects",
          cfg$label,
          data_subset
        )
        pred_cors[[i]] <- collect_cors(
          true_phenos,
          estimated_phenos,
          "IID",
          method,
          "phenotypes",
          cfg$label,
          data_subset
        )
        i <- i + 1
        cli_progress_update(id = pb_sk)
      } # end cov_tags_sklearn loop
    }
    #cli_progress_done(id = pb_sk)
  }

  all_pred_cors <- as.data.frame(do.call(rbind, pred_cors))
  all_eff_cors <- as.data.frame(do.call(rbind, eff_cors))
  all_true_max_cors <- as.data.frame(do.call(rbind, true_max_cors))
  all_method_params <- as.data.frame(do.call(rbind, method_params))

  all_pred_cors <- expand_trait_names(all_pred_cors)
  all_eff_cors <- expand_trait_names(all_eff_cors)
  all_method_params <- expand_trait_names(all_method_params)

  all_pred_cors_with_truth <- left_join(all_pred_cors, all_true_max_cors)
  all_eff_cors_with_truth <- left_join(all_eff_cors, all_true_max_cors)

  all_pred_cors$seed[all_pred_cors$seed == YEAST_SEED] <- "Yeast"
  all_eff_cors$seed[all_eff_cors$seed == YEAST_SEED] <- "Yeast"
  if (nchar(HUMAN_SEED) > 0) {
    all_pred_cors$seed[all_pred_cors$seed == HUMAN_SEED] <- "Human"
    all_eff_cors$seed[all_eff_cors$seed == HUMAN_SEED] <- "Human"
  }

  all_pred_cors_with_truth <- left_join(all_pred_cors, all_true_max_cors)
  all_eff_cors_with_truth <- left_join(all_eff_cors, all_true_max_cors)

  all_method_params$seed[all_method_params$seed == YEAST_SEED] <- "Yeast"
  if (nchar(HUMAN_SEED) > 0) {
    all_method_params$seed[all_method_params$seed == HUMAN_SEED] <- "Human"
  }

  all_method_params_with_truth <- left_join(
    all_method_params,
    all_true_max_cors
  )
  all_method_params_with_truth[
    all_method_params_with_truth$method == "lasso",
    "l1_ratio"
  ] <- 1
  if ("lars" %in% all_method_params_with_truth$method) {
    lars_rows <- all_method_params_with_truth$method == "lars"
    all_method_params_with_truth[lars_rows, "method"] <- paste0(
      "lars_maxiter",
      unlist(all_method_params_with_truth[lars_rows, "maxiter"])
    )
  }

  arrow::write_feather(all_pred_cors_with_truth, INT_PRED_CORS)
  arrow::write_feather(all_eff_cors_with_truth, INT_EFF_CORS)
  arrow::write_feather(all_true_max_cors, INT_MAX_CORS)
  arrow::write_feather(all_method_params_with_truth, INT_PARAMS)
} else {
  cat("Loading correlations...")
  all_pred_cors_with_truth <- read_feather(INT_PRED_CORS)
  all_eff_cors_with_truth <- read_feather(INT_EFF_CORS)
  all_true_max_cors <- read_feather(INT_MAX_CORS)
  all_method_params_with_truth <- read_feather(INT_PARAMS)
}

# Reconstruct trait column
reconstruct_trait <- function(df) {
  df %>%
    mutate(
      traitname = "trait",
      vavgname = "vavg",
      numQTLname = "numQTL",
      numChrname = "numChr",
      H2name = "H2"
    ) %>%
    unite(
      "trait",
      traitname,
      overlap,
      vavgname,
      target_vavg,
      numQTLname,
      numQTL,
      numChrname,
      numChr,
      H2name,
      target_H2,
      sep = "_",
      remove = FALSE
    )
}
all_pred_cors_with_truth <- reconstruct_trait(all_pred_cors_with_truth)
all_eff_cors_with_truth <- reconstruct_trait(all_eff_cors_with_truth)

# Trait filtering
traits_with_lots <- all_pred_cors_with_truth[
  !all_pred_cors_with_truth$method %in%
    c(
      "plink",
      "plink_clumped",
      "plink_cov6",
      "plink_clumped_cov6",
      "lars_maxiter10000"
    ),
] %>%
  group_by(trait, seed) %>%
  summarize(n = n(), .groups = "drop") %>%
  filter(n > 2)


compare_phenos_yeast <- traits_with_lots$trait[grepl(
  "numChr_16_",
  traits_with_lots$trait
)]
compare_phenos_human <- traits_with_lots$trait[grepl(
  "numChr_1_",
  traits_with_lots$trait
)]
#print(length(compare_phenos_yeast))
#print(length(compare_phenos_human))

traits_with_all <- all_pred_cors_with_truth[
  !all_pred_cors_with_truth$method %in% c("lars_maxiter10000"),
] %>%
  group_by(trait, seed) %>%
  summarize(n = n(), .groups = "drop") %>%
  arrange(desc(n)) %>%
  filter(n == 9 & seed == "Human" | n >= 8 & seed == "Yeast")


traits_with_two_numQTL <- all_pred_cors_with_truth[
  !all_pred_cors_with_truth$method %in% c("lars_maxiter10000"),
] %>%
  group_by(trait, numQTL, seed) %>%
  summarize(n = n()) %>%
  mutate(trait2 = gsub("_numQTL_[0-9]+_", "_numQTL_", trait)) %>%
  filter(n == 9 & seed == "Human" | n >= 8 & seed == "Yeast") %>%
  group_by(trait2) %>%
  mutate(n2 = n())


# Auto-select example traits for figures from available data.
# These replace hardcoded trait names from the original script that may not
# match the simulation parameters used in the current run.
yeast_qtl_vals <- sort(unique(as.numeric(
  regmatches(
    compare_phenos_yeast,
    regexpr("(?<=numQTL_)\\d+", compare_phenos_yeast, perl = TRUE)
  )
)))
yeast_qtl_high <- if (length(yeast_qtl_vals) > 0) max(yeast_qtl_vals) else NULL
yeast_qtl_low <- if (length(yeast_qtl_vals) > 1) {
  min(yeast_qtl_vals)
} else {
  yeast_qtl_high
}

yeast_h2_dummy_high <- "trait_overlap_vavg_0.5_numQTL_500_numChr_16_H2_0.25"
yeast_h2_dummy_low <- "trait_overlap_vavg_0.5_numQTL_100_numChr_16_H2_0.25"

if (nchar(HUMAN_SEED) > 0 && length(compare_phenos_human) > 0) {
  human_qtl_vals <- sort(unique(as.numeric(
    regmatches(
      compare_phenos_human,
      regexpr("(?<=numQTL_)\\d+", compare_phenos_human, perl = TRUE)
    )
  )))
  human_qtl_high <- if (length(human_qtl_vals) > 0) {
    max(human_qtl_vals)
  } else {
    NULL
  }
  human_qtl_low <- if (length(human_qtl_vals) > 1) {
    min(human_qtl_vals)
  } else {
    human_qtl_high
  }
  human_h2_dummy_high <- "trait_overlap_vavg_0.5_numQTL_100_numChr_1_H2_0.25"
  human_h2_dummy_low <- "trait_overlap_vavg_0.5_numQTL_50_numChr_1_H2_0.25"
} else {
  human_h2_dummy_high <- ""
  human_h2_dummy_low <- ""
}

# ---------------------------------------------------------------------------
# Collect per-SNP beta params
# ---------------------------------------------------------------------------
if (opt$`collect-cors-params`) {
  for (seed_chr in names(SEEDS)) {
    cat("Collecting betas...\n")
    all_betas <- list()
    i <- 1
    cfg <- SEEDS[[seed_chr]]
    file_prefix <- cfg$prefix
    p_thr <- cfg$plink_threshold
    is_biobank <- seed_chr == HUMAN_SEED
    phenos <- if (is_biobank) compare_phenos_human else compare_phenos_yeast
    cli_alert_info("Collecting betas for seed: {seed_chr}")

    seed_int_dir <- file.path(BASE_DIR, paste0("intermediates_seed_", seed_chr))
    true_effects <- read_feather(file.path(
      seed_int_dir,
      paste0(file_prefix, "_seed_", seed_chr, "_additive_effects_wide.feather")
    )) %>%
      separate_wider_delim(
        SNP,
        names = c("SNP_clean", NA),
        delim = "_",
        too_few = "align_start"
      ) %>%
      rename(SNP = SNP_clean) %>%
      select(SNP, any_of(phenos))
    true_effects$method <- "truth"
    true_effects$cov <- ""
    true_effects$split <- "none"
    all_betas[[i]] <- true_effects
    i <- i + 1

    for (data_subset in c("test", "train")) {
      cov_tags <- if (is_biobank) {
        c("", "cov0_", "cov6_", "cov10_", "cov15_", "cov40_")
      } else {
        c("")
      }

      for (cov_tag in cov_tags) {
        plink_dir <- file.path(
          BASE_DIR,
          paste0("aggregate_outputs_", cov_tag, data_subset, "_seed_", seed_chr)
        )
        plink_suff <- paste0(cov_tag, data_subset, "_seed_", seed_chr)
        clump_suff <- if (is_biobank) "_clumped" else "_clumped_r2_0.5"

        if (!dir.exists(plink_dir)) {
          cat("DEBUG plink_dir not found, skipping:", plink_dir, "\n")
          next
        }

        beta_file <- file.path(
          plink_dir,
          paste0("beta_matrix_filtered_", plink_suff, "_p", p_thr, ".feather")
        )
        if (!file.exists(beta_file)) {
          cat("DEBUG beta_file not found, skipping:", beta_file, "\n")
          next
        }

        estimated_effects <- read_feather(beta_file) %>%
          select(SNP, any_of(phenos))
        estimated_effects_unfiltered <- read_feather(file.path(
          plink_dir,
          paste0("beta_matrix_", plink_suff, ".feather")
        )) %>%
          select(SNP, any_of(phenos))

        pval_file <- file.path(
          plink_dir,
          paste0("pvalue_matrix_", plink_suff, ".feather")
        )
        pvals <- if (file.exists(pval_file)) {
          read_feather(pval_file) %>% select(SNP, any_of(phenos))
        } else {
          NULL
        }

        estimated_effects <- estimated_effects[
          estimated_effects$SNP != "__bias__",
        ]
        estimated_effects <- estimated_effects[
          match(true_effects$SNP, estimated_effects$SNP),
        ]
        estimated_effects_unfiltered <- estimated_effects_unfiltered[
          estimated_effects_unfiltered$SNP != "__bias__",
        ]
        estimated_effects_unfiltered <- estimated_effects_unfiltered[
          match(true_effects$SNP, estimated_effects_unfiltered$SNP),
        ]
        if (!is.null(pvals)) {
          pvals <- pvals[match(true_effects$SNP, pvals$SNP), ]
        }

        estimated_effects$method <- "plink"
        estimated_effects$cov <- cov_tag
        estimated_effects$split <- data_subset
        all_betas[[i]] <- estimated_effects
        i <- i + 1

        estimated_effects_unfiltered$method <- "plink_unfiltered"
        estimated_effects_unfiltered$cov <- cov_tag
        estimated_effects_unfiltered$split <- data_subset
        all_betas[[i]] <- estimated_effects_unfiltered
        i <- i + 1

        if (!is.null(pvals)) {
          pvals$method <- "plink_pval"
          pvals$cov <- cov_tag
          pvals$split <- data_subset
          all_betas[[i]] <- pvals
          i <- i + 1
        }

        clump_beta_file <- file.path(
          plink_dir,
          paste0(
            "beta_matrix_filtered_",
            plink_suff,
            "_p",
            p_thr,
            clump_suff,
            ".feather"
          )
        )
        clump_pval_file <- file.path(
          plink_dir,
          paste0("pvalue_matrix_", plink_suff, clump_suff, ".feather")
        )

        if (file.exists(clump_beta_file)) {
          ec_clumped <- read_feather(clump_beta_file) %>%
            select(SNP, any_of(phenos))
          ec_clumped <- ec_clumped[match(true_effects$SNP, ec_clumped$SNP), ]
          ec_clumped$method <- "plink_clumped"
          ec_clumped$cov <- cov_tag
          ec_clumped$split <- data_subset
          all_betas[[i]] <- ec_clumped
          i <- i + 1
        }
        if (file.exists(clump_pval_file)) {
          pv_clumped <- read_feather(clump_pval_file) %>%
            select(SNP, any_of(phenos))
          pv_clumped <- pv_clumped[match(true_effects$SNP, pv_clumped$SNP), ]
          pv_clumped$method <- "plink_pval_clumped"
          pv_clumped$cov <- cov_tag
          pv_clumped$split <- data_subset
          all_betas[[i]] <- pv_clumped
          i <- i + 1
        }
      } # end cov_tag loop
    } # end data_subset loop
    #cli_progress_done(id = pb_betas)
    betas_methods_run <- c(
      "lasso",
      LARS_TAG,
      "ridge",
      "elasticnet",
      "pytorch_ridge"
    )
    cov_tags_sklearn <- if (is_biobank) c("", "cov6_") else c("")

    for (method in betas_methods_run) {
      data_subset <- if (method == "pytorch_ridge") "train" else "test"
      for (cov_tag in cov_tags_sklearn) {
        method_dir <- file.path(
          BASE_DIR,
          paste0(
            "aggregate_outputs_",
            method,
            "_",
            cov_tag,
            data_subset,
            "_seed_",
            seed_chr
          )
        )
        suff <- paste0(
          method,
          "_",
          cov_tag,
          data_subset,
          "_seed_",
          seed_chr,
          ".feather"
        )
        coeff_file <- file.path(method_dir, paste0("coeff_matrix_", suff))
        cat(
          "DEBUG sklearn/pytorch method:",
          method,
          "| cov_tag:",
          cov_tag,
          "| coeff_file exists:",
          file.exists(coeff_file),
          "\n"
        )
        if (!file.exists(coeff_file)) {
          next
        }
        estimated_effects <- read_feather(coeff_file)
        estimated_effects <- estimated_effects[
          estimated_effects$SNP != "__bias__",
        ] %>%
          separate_wider_delim(
            SNP,
            names = c("SNP_clean", NA),
            delim = "_",
            too_few = "align_start"
          ) %>%
          rename(SNP = SNP_clean) %>%
          select(SNP, any_of(phenos))
        estimated_effects <- estimated_effects[
          match(true_effects$SNP, estimated_effects$SNP),
        ]
        estimated_effects$method <- method
        estimated_effects$cov <- cov_tag
        estimated_effects$split <- data_subset
        all_betas[[i]] <- estimated_effects
        i <- i + 1
      } # end cov_tags_sklearn loop
    } # end betas_methods_run loop

    combined <- do.call(bind_rows, all_betas)

    if (seed_chr == YEAST_SEED) {
      combined_all_betas_yeast <- combined
      arrow::write_feather(
        combined_all_betas_yeast,
        INT_BETAS_YEAST
      )
    } else if (seed_chr == HUMAN_SEED) {
      combined_all_betas_human <- combined
      arrow::write_feather(combined_all_betas_human, INT_BETAS_HUMAN)
    }
  }
} else {
  combined_all_betas_yeast <- read_feather(INT_BETAS_YEAST)

  if (nchar(HUMAN_SEED) > 0) {
    combined_all_betas_human <- read_feather(INT_BETAS_HUMAN)
  }
}

# ---------------------------------------------------------------------------
# Recalc wider/longer pivot tables
# ---------------------------------------------------------------------------
if (opt$`recalc-wider-longer`) {
  cat("Pivoting dataframes...\n")
  make_littlelonger <- function(combined, snp_info) {
    longer <- combined %>%
      pivot_longer(
        cols = -c("SNP", "method", "cov", "split"),
        names_to = "trait",
        values_to = "coeff"
      )
    littlewider <- longer %>%
      mutate(method_append = paste0(method, "_", cov, split)) %>%
      select(-method, -cov, -split) %>%
      pivot_wider(
        id_cols = c("SNP", "trait"),
        names_from = "method_append",
        values_from = "coeff"
      )
    littlelonger <- littlewider %>%
      pivot_longer(
        cols = -c("SNP", "trait", "truth_none"),
        names_to = "method_append",
        values_to = "coeff"
      )
    left_join(littlelonger, snp_info, by = "SNP") %>%
      mutate(pos_along = chrom_sum + POS)
  }
  yeast_littlelonger_with_fullinfo <- make_littlelonger(
    combined_all_betas_yeast,
    full_snp_info_yeast
  )
  arrow::write_feather(
    yeast_littlelonger_with_fullinfo,
    INT_LONGER_YEAST
  )
  if (nchar(HUMAN_SEED) > 0) {
    human_littlelonger_with_fullinfo <- make_littlelonger(
      combined_all_betas_human,
      full_snp_info_human
    )
    arrow::write_feather(
      human_littlelonger_with_fullinfo,
      INT_LONGER_HUMAN
    )
  }
} else {
  yeast_littlelonger_with_fullinfo <- read_feather(INT_LONGER_YEAST)
  if (nchar(HUMAN_SEED) > 0) {
    human_littlelonger_with_fullinfo <- read_feather(INT_LONGER_HUMAN)
  }
}

# ---------------------------------------------------------------------------
# Figure 2/3
# ---------------------------------------------------------------------------

# ── Shared setup ──────────────────────────────────────────────────────────────

methods_for_final_plot <- c(
  "elasticnet",
  "lasso",
  "lars_maxiter1000",
  "pytorch_ridge",
  "ridge"
)

line_segments <- data.frame(
  target_H2_2 = factor(
    c(0.25, 0.5, 0.75),
    levels = c(0.25, 0.5, 0.75, 0.99, 1)
  ),
  x = 0.6,
  y = c(0.25, 0.5, 0.75),
  xend = 1.4,
  yend = c(0.25, 0.5, 0.75)
)

h2_levels <- c(0.25, 0.5, 0.75, 0.99, 1)
method_levels <- c(methods_for_final_plot, "plink_clumped", "plink")

all_pred_cors_with_truth$target_H2_2 <- factor(
  all_pred_cors_with_truth$target_H2,
  levels = h2_levels
)
all_pred_cors_with_truth$method_2 <- factor(
  all_pred_cors_with_truth$method,
  levels = method_levels
)
all_eff_cors_with_truth$target_H2_2 <- factor(
  all_eff_cors_with_truth$target_H2,
  levels = h2_levels
)
all_eff_cors_with_truth$method_2 <- factor(
  all_eff_cors_with_truth$method,
  levels = method_levels
)

# ── Helper: shared row filter ─────────────────────────────────────────────────

filter_methods <- function(df, seed) {
  df[
    df$seed == seed &
      df$overlap != "singlesnp" &
      ((df$method == "plink" & df$split == "train") |
        (df$method == "plink_clumped" & df$split == "train") |
        df$method %in% methods_for_final_plot) &
      (df$cov %in% c("", "cov6_") | is.na(df$cov)),
  ]
}

# ── Main plotting function ────────────────────────────────────────────────────

make_predictability_plots <- function(
  seed,
  trait_subset_points,
  trait_subset_barplot,
  legend_on_boxplot = TRUE
) {
  pred <- filter_methods(all_pred_cors_with_truth, seed)
  eff <- filter_methods(all_eff_cors_with_truth, seed)

  colour_scale <- scale_color_manual(
    values = METHOD_COLOURS,
    name = "Method",
    labels = METHOD_LABELS,
    drop = FALSE
  )

  if (legend_on_boxplot) {
    guide_right <- guides(
      color = guide_legend(override.aes = list(size = 4, alpha = 1))
    )
    legend_place_right <- theme(legend.position = c(0.16, 0.67))
    guide_points <- guides(color = "none")
    legend_place_pts <- NULL
  } else {
    guide_right <- guides(color = "none")
    legend_place_right <- NULL
    guide_points <- guides(
      color = guide_legend(ncol = 3, override.aes = list(size = 4, alpha = 1))
    )
    legend_place_pts <- theme(legend.position = c(0.4, 0.8))
  }

  # A-left: phenotype predictability, narrow-sense traits
  pred_narrow <- ggplot(pred[pred$trait %in% trait_subset_points, ]) +
    geom_abline(slope = 1, intercept = 0) +
    geom_point(aes(x = true_h2, y = r2, color = method_2)) +
    xlim(0, 1) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      y = expression(italic(r)^2 ~ "(true phenotype, predicted phenotype)"),
      tag = "A)"
    ) +
    colour_scale +
    theme_pub(x_axis_type = "numerical") +
    theme(
      axis.title.x = element_blank(),
      legend.key.spacing.y = unit(0.3, "cm"),
      legend.background = element_rect(fill = "transparent", color = NA)
    ) +
    guide_points +
    legend_place_pts

  # A-right: phenotype predictability, broad-sense traits — mean ± 1 SE
  pred_broad <- ggplot(
    pred[
      pred$trait %in%
        trait_subset_barplot &
        pred$target_H2 %in% c(0.25, 0.5, 0.75),
    ] %>%
      {
        if (seed == "Yeast") filter(., method != "elasticnet") else .
      }
  ) +
    geom_segment(
      data = line_segments,
      aes(x = x, y = y, xend = xend, yend = yend)
    ) +
    stat_summary(
      aes(x = target_H2_2, y = r2, color = method_2),
      fun.data = mean_se,
      geom = "errorbar",
      width = 0.3,
      position = position_dodge(width = 0.7),
      show.legend = FALSE
    ) +
    stat_summary(
      aes(x = target_H2_2, y = r2, color = method_2),
      fun = mean,
      geom = "point",
      position = position_dodge(width = 0.7)
    ) +
    facet_wrap(~target_H2_2, scales = "free_x") +
    scale_x_discrete(
      drop = TRUE,
      expand = expansion(add = c(0.4, 0), mult = c(0, 0))
    ) +
    scale_y_continuous(limits = c(0, 0.8), breaks = c(0, 0.25, 0.5, 0.75)) +
    colour_scale +
    labs(tag = "B)") +
    theme_pub(x_axis_type = "numerical") +
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      legend.key.spacing.y = unit(0.3, "cm"),
      panel.grid.major.x = element_blank()
    ) +
    guide_right +
    legend_place_right

  # B-left: effect-size predictability, narrow-sense traits
  beta_narrow <- ggplot(eff[eff$trait %in% trait_subset_points, ]) +
    geom_point(aes(x = true_h2, y = r2, color = method_2)) + # size = 1.5
    guides(color = "none") +
    xlim(0, 1) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      y = expression(italic(r)^2 ~ "(true effect, predicted effect)"),
      tag = "C)"
    ) +
    xlab("Relative var(Additive Effects)") +
    colour_scale +
    theme_pub(x_axis_type = "numerical")

  # B-right: effect-size predictability, broad-sense traits — mean ± 1 SE
  beta_broad <- ggplot(
    eff[
      eff$trait %in%
        trait_subset_barplot &
        eff$target_H2 %in% c(0.25, 0.5, 0.75),
    ] %>%
      {
        if (seed == "Yeast") filter(., method != "elasticnet") else .
      }
  ) +
    stat_summary(
      aes(x = target_H2_2, y = r2, color = method_2),
      fun.data = mean_se,
      geom = "errorbar",
      width = 0.3,
      position = position_dodge(width = 0.7),
      show.legend = FALSE
    ) +
    stat_summary(
      aes(x = target_H2_2, y = r2, color = method_2),
      fun = mean,
      geom = "point",
      position = position_dodge(width = 0.7)
    ) +
    facet_wrap(~target_H2_2, scales = "free_x") +
    scale_x_discrete(
      drop = TRUE,
      expand = expansion(add = c(0.4, 0), mult = c(0, 0))
    ) +
    labs(x = "Broad-Sense Heritability", tag = "D)") +
    colour_scale +
    theme_pub(x_axis_type = "numerical") +
    theme(
      axis.title.y = element_blank(),
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      panel.grid.major.x = element_blank()
    ) +
    guides(color = "none")
  (pred_narrow + pred_broad) / (beta_narrow + beta_broad)
}

# ── Render plots ──────────────────────────────────────────────────────────────

if (nchar(HUMAN_SEED) > 0) {
  human_plot <- make_predictability_plots(
    "Human",
    traits_with_lots$trait,
    traits_with_all$trait,
    legend_on_boxplot = TRUE
  )
  ggsave(
    file.path(OUTPUT_DIR, "Human_r2.svg"),
    plot = human_plot,
    width = 12,
    height = 12,
    dpi = 100,
    units = "in"
  )
}

yeast_plot <- make_predictability_plots(
  "Yeast",
  traits_with_lots$trait,
  traits_with_all$trait,
  legend_on_boxplot = FALSE
)
ggsave(
  file.path(OUTPUT_DIR, "Yeast_r2.svg"),
  plot = yeast_plot,
  width = 12,
  height = 12,
  dpi = 100,
  units = "in"
)

# ---------------------------------------------------------------------------
# Figure 4 — by heritability and by method
# ---------------------------------------------------------------------------

# Plot grouped by method (x = method), with mean ± SE and stat tests
make_method_plot <- function(
  data,
  seed,
  type,
  dummy_trait_high,
  dummy_trait_low,
  low_val,
  ylim = NULL,
  ylab,
  ybreaks = NULL,
  hide_ylab = FALSE,
  hide_xlab = FALSE,
  tag,
  method_ordered,
  include_legend
) {
  df <- data[data$seed == seed, ] %>%
    filter(
      method != "lars_maxiter10000",
      (method %in% c("plink", "plink_clumped")) &
        split == "train" |
        !grepl("plink", method)
    ) %>%
    {
      if (seed == "Yeast") filter(., method != "elasticnet") else .
    } %>%
    mutate(mod = ifelse(numQTL == low_val, "small", "big"))
  df$mod <- factor(df$mod, levels = c("small", "big"))
  df$numQTL <- as.factor(df$numQTL)
  df$target_H2 <- as.factor(df$target_H2)
  df$method <- factor(df$method, levels = method_ordered)
  df <- df %>% mutate(inter = interaction(method, mod))

  stat_test <- df %>%
    group_by(method) %>%
    wilcox_test(r2 ~ numQTL) %>%
    add_xy_position(x = "method", fun = "mean", dodge = 0.8)

  first_method <- if (seed == "Human") "elasticnet" else "lasso"
  stat_test <- stat_test %>%
    mutate(
      label = if_else(
        method == first_method,
        paste0("p = ", signif(p, 2)),
        as.character(signif(p, 2))
      )
    )

  if (seed == "Human" && type == "pred") {
    stat_test$y.position[!grepl("plink", stat_test$method)] <-
      0.8 * stat_test$y.position[!grepl("plink", stat_test$method)]
  } else if (seed == "Human" && type == "effect") {
    stat_test$y.position[
      grepl("ridge", stat_test$method) & !grepl("plink", stat_test$method)
    ] <-
      0.8 *
      stat_test$y.position[
        grepl("ridge", stat_test$method) & !grepl("plink", stat_test$method)
      ]
    stat_test$y.position[
      !grepl("ridge", stat_test$method) & !grepl("plink", stat_test$method)
    ] <-
      0.7 *
      stat_test$y.position[
        !grepl("ridge", stat_test$method) & !grepl("plink", stat_test$method)
      ]
  } else if (seed == "Yeast" && type == "pred") {
    stat_test$y.position[grepl("plink", stat_test$method)] <- 0.9 *
      stat_test$y.position[grepl("plink", stat_test$method)]
    stat_test$y.position[!grepl("plink", stat_test$method)] <- 0.8 *
      stat_test$y.position[!grepl("plink", stat_test$method)]
  } else if (seed == "Yeast" && type == "effect") {
    stat_test$y.position <- 0.8 * stat_test$y.position
  }

  g <- ggplot(df, aes(y = r2, x = method, color = inter)) +
    scale_x_discrete(labels = METHOD_LABELS) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)), limits = ylim) +
    stat_summary(
      fun = mean,
      geom = "point",
      size = 3,
      shape = 19,
      position = position_dodge(width = 0.75),
      show.legend = FALSE
    ) +
    stat_summary(
      fun.data = mean_se,
      geom = "errorbar",
      width = 0.2,
      linewidth = 0.8,
      position = position_dodge(width = 0.75),
      show.legend = FALSE
    ) +
    labs(y = ylab, x = "Method", alpha = "# QTL", tag = tag) +
    stat_pvalue_manual(
      stat_test,
      label = "label",
      vjust = -0.5,
      tip.length = 0.01,
      bracket.nudge.y = 0.05
    ) +
    theme_pub(x_axis_type = "discrete") +
    scale_colour_manual(
      values = METHOD_COLOURS,
      name = "Method",
      labels = METHOD_LABELS
    ) +
    guides(
      color = "none",
      alpha = guide_legend(override.aes = list(size = 4))
    )
  if (hide_ylab) {
    g <- g +
      theme(axis.title.y = element_blank(), axis.text.y = element_blank())
  } else {
    g <- g +
      theme(axis.text.y = element_text(size = 10))
  }
  if (hide_xlab) {
    g <- g +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10)
      )
  }
  if (include_legend) {
    g <- g + theme(legend.position = c(0.9, 0.9))
  } else {
    g <- g + guides(alpha = "none")
  }
  if (!is.null(ybreaks)) {
    g <- g + scale_y_continuous(limits = ylim, breaks = ybreaks)
  }
  g
}

all_pred_cors_with_broad_trait <- left_join(
  all_pred_cors_with_truth,
  traits_with_two_numQTL
) %>%
  select(
    method,
    task,
    seed,
    split,
    r2,
    target_vavg,
    numQTL,
    target_H2,
    target_h2,
    true_H2,
    true_h2,
    true_vavg,
    trait,
    trait2,
    n2
  ) %>%
  filter(!is.na(trait2), n2 == 2) #%>%pivot_wider(id_cols = c("method","task","seed","split","target_vavg","target_H2","target_h2"),names_from = numQTL,values_from=r2)


all_eff_cors_with_broad_trait <- left_join(
  all_eff_cors_with_truth,
  traits_with_two_numQTL
) %>%
  select(
    method,
    task,
    seed,
    split,
    r2,
    target_vavg,
    numQTL,
    target_H2,
    target_h2,
    true_H2,
    true_h2,
    true_vavg,
    trait,
    trait2,
    n2
  ) %>%
  filter(!is.na(trait2), n2 == 2) #%>%pivot_wider(id_cols = c("method","task","seed","split","target_vavg","target_H2","target_h2"),names_from = numQTL,values_from=r2)


# Figure 4b — by method
yeast_pred_method <- make_method_plot(
  all_pred_cors_with_broad_trait,
  "Yeast",
  type = "pred",
  dummy_trait_high = yeast_h2_dummy_high,
  dummy_trait_low = yeast_h2_dummy_low,
  low_val = yeast_qtl_low,
  ylim = c(0, 0.5),
  hide_ylab = TRUE,
  hide_xlab = TRUE,
  ylab = expression(
    italic(r)^2 ~ "(true phenotype, predicted phenotype)"
  ),
  method_ordered = METHOD_ORDER_YEAST,
  tag = "B)",
  include_legend = TRUE
)

yeast_eff_method <- make_method_plot(
  all_eff_cors_with_broad_trait,
  "Yeast",
  type = "effect",
  dummy_trait_high = yeast_h2_dummy_high,
  dummy_trait_low = yeast_h2_dummy_low,
  low_val = yeast_qtl_low,
  ylim = c(0, 0.5),
  hide_ylab = TRUE,
  ylab = expression(italic(r)^2 ~ "(true effect, predicted effect)"),
  method_ordered = METHOD_ORDER_YEAST,
  tag = "D)",
  include_legend = FALSE
)

if (nchar(HUMAN_SEED) > 0) {
  human_pred_method <- make_method_plot(
    all_pred_cors_with_broad_trait,
    "Human",
    type = "pred",
    dummy_trait_high = human_h2_dummy_high,
    dummy_trait_low = human_h2_dummy_low,
    low_val = human_qtl_low,
    ylim = c(0, 0.5),
    hide_xlab = TRUE,
    ylab = expression(atop(
      italic(r)^2 ~ "(true phenotype,",
      "predicted phenotype)"
    )),
    method_ordered = METHOD_ORDER_HUMAN,
    tag = "A)",
    include_legend = TRUE
  )
  human_eff_method <- make_method_plot(
    all_eff_cors_with_broad_trait,
    "Human",
    type = "effect",
    dummy_trait_high = human_h2_dummy_high,
    dummy_trait_low = human_h2_dummy_low,
    low_val = human_qtl_low,
    ylim = c(0, 0.5),
    ylab = expression(atop(italic(r)^2 ~ "(true effect,", "predicted effect)")),
    method_ordered = METHOD_ORDER_HUMAN,
    tag = "C)",
    include_legend = FALSE
  )
  by_method <- (human_pred_method + yeast_pred_method) /
    (human_eff_method + yeast_eff_method)
} else {
  by_method <- (yeast_pred_method) / (yeast_eff_method)
}
ggsave(
  file.path(OUTPUT_DIR, "YeastHuman_numQTL.svg"),
  plot = by_method,
  width = 13,
  height = 8,
  dpi = 100,
  units = "in"
)

# ---------------------------------------------------------------------------
# Figures 6/7 — per-trait beta correlation scatter plots
# ---------------------------------------------------------------------------
build_explore_compare <- function(trait) {
  yeast <- yeast_littlelonger_with_fullinfo[
    yeast_littlelonger_with_fullinfo$trait == trait,
  ] %>%
    mutate(species = "Yeast") %>%
    select(SNP, trait, truth_none, method_append, coeff, species)

  if (nchar(HUMAN_SEED) > 0 && exists("human_littlelonger_with_fullinfo")) {
    human_trait_name <- gsub("_numChr_16_", "_numChr_1_", trait)
    human <- human_littlelonger_with_fullinfo[
      human_littlelonger_with_fullinfo$trait == human_trait_name,
    ] %>%
      mutate(species = "Human") %>%
      select(SNP, trait, truth_none, method_append, coeff, species)
    expand_trait_names(rbind(yeast, human))
  } else {
    expand_trait_names(yeast)
  }
}

compute_r2_labels <- function(df) {
  if (nrow(df) == 0) {
    cat("DEBUG compute_r2_labels: empty dataframe, returning empty labels\n")
    return(tibble(
      overlap = character(),
      target_vavg = numeric(),
      numQTL = numeric(),
      numChr = numeric(),
      target_H2 = numeric(),
      method_append = character(),
      species = character(),
      target_h2 = numeric(),
      r_squared = numeric(),
      label = character()
    ))
  }

  grp_summary <- df %>%
    group_by(across(all_of(r2_group_vars))) %>%
    summarise(
      n = n(),
      n_complete = sum(!is.na(truth_none) & !is.na(coeff)),
      n_nonzero_coeff = sum(coeff != 0, na.rm = TRUE),
      .groups = "drop"
    )

  df %>%
    group_by(across(all_of(r2_group_vars))) %>%
    summarise(
      r_squared = cor(truth_none, coeff, use = "complete.obs")^2,
      .groups = "drop"
    ) %>%
    mutate(label = paste0("R\u00b2 = ", signif(r_squared, 2)))
}

filter_plink_df <- function(df) {
  out <- df[complete.cases(df) & df$method_append %in% PLINK_METHODS, ]
  out
}

r2_text_layer <- function(labels, species) {
  geom_text(
    data = labels[labels$species == species, ],
    aes(label = label),
    x = Inf,
    y = -Inf,
    hjust = 1.1,
    vjust = -0.5,
    fontface = "bold",
    color = "#43413F"
  )
}

make_beta_cor_plot <- function(
  df,
  r2_labels,
  species,
  title = NULL,
  tag = NULL,
  hide_x_title = FALSE,
  hide_y_title = FALSE,
  y_formatted = TRUE,
  x_scale = scale_x_continuous(breaks = seq(-0.3, 0.3, by = 0.3))
) {
  p <- ggplot(
    df[df$species == species, ],
    aes(x = truth_none, y = coeff, color = method_append)
  ) +
    geom_point() +
    facet_wrap(
      ~method_append,
      labeller = labeller(method_append = METHOD_LABELS),
      ncol = 5
    ) +
    x_scale +
    xlab("True Effect") +
    ylab("Estimated Effect") +
    theme_pub(x_axis_type = "numerical") +
    scale_color_manual(values = METHOD_COLOURS) +
    theme(panel.grid.minor = element_blank()) +
    r2_text_layer(r2_labels, species) +
    guides(color = "none")
  if (y_formatted) {
    p <- p + scale_y_continuous(labels = scales::label_number(accuracy = 0.001))
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  if (!is.null(tag)) {
    p <- p + labs(tag = tag)
  }
  if (hide_x_title) {
    p <- p + theme(axis.title.x = element_blank())
  }
  if (hide_y_title) {
    p <- p + theme(axis.title.y = element_blank())
  }
  p
}

make_trait_plots <- function(trait, tag, hide_x = FALSE) {
  explore <- build_explore_compare(trait)

  sparse <- explore[
    complete.cases(explore) & explore$method_append %in% SPARSE_METHODS,
  ]
  sparse$method_append <- factor(sparse$method_append, levels = SPARSE_METHODS)
  r2_sparse <- compute_r2_labels(sparse)
  plink_data <- filter_plink_df(explore)
  r2_plink <- compute_r2_labels(plink_data)
  ridge_data <- explore[
    complete.cases(explore) & explore$method_append %in% RIDGE_METHODS,
  ]
  r2_ridge <- compute_r2_labels(ridge_data)
  yeast_x <- scale_x_continuous(limits = c(-0.3, 0.3))
  human_x <- scale_x_continuous(breaks = seq(-0.3, 0.3, by = 0.3))
  plots <- list(
    yeast_sparse = make_beta_cor_plot(
      sparse[!grepl("plink", sparse$method_append), ],
      r2_sparse,
      "Yeast",
      tag = tag,
      hide_x_title = TRUE,
      hide_y_title = FALSE,
      x_scale = yeast_x
    ),
    yeast_ridge = make_beta_cor_plot(
      ridge_data,
      r2_ridge,
      "Yeast",
      hide_x_title = hide_x,
      hide_y_title = TRUE,
      x_scale = yeast_x
    ),
    yeast_plink = make_beta_cor_plot(
      plink_data,
      r2_plink,
      "Yeast",
      hide_x_title = TRUE,
      hide_y_title = TRUE,
      y_formatted = FALSE,
      x_scale = yeast_x
    )
  )
  if (nchar(HUMAN_SEED) > 0 && exists("human_littlelonger_with_fullinfo")) {
    plots$human_sparse <- make_beta_cor_plot(
      sparse[!grepl("plink", sparse$method_append), ],
      r2_sparse,
      "Human",
      tag = tag,
      hide_x_title = TRUE,
      hide_y_title = FALSE,
      x_scale = human_x
    )
    plots$human_ridge <- make_beta_cor_plot(
      ridge_data,
      r2_ridge,
      "Human",
      hide_x_title = hide_x,
      hide_y_title = TRUE,
      x_scale = human_x
    )
    plots$human_plink <- make_beta_cor_plot(
      plink_data,
      r2_plink,
      "Human",
      hide_x_title = TRUE,
      hide_y_title = TRUE,
      y_formatted = FALSE,
      x_scale = human_x
    )
  }
  plots
}

p1 <- make_trait_plots(
  "trait_overlap_vavg_0.9_numQTL_100_numChr_16_H2_0.25",
  "A)",
  hide_x = TRUE
)
p2 <- make_trait_plots(
  "trait_overlap_vavg_0.98_numQTL_100_numChr_16_H2_0.5",
  "B)"
)

yeast_cor_plot <- (p1$yeast_sparse +
  p1$yeast_ridge +
  p1$yeast_plink +
  plot_layout(widths = c(30, 28.5, 29))) /
  (p2$yeast_sparse +
    p2$yeast_ridge +
    p2$yeast_plink +
    plot_layout(widths = c(30, 28.5, 29)))
ggsave(
  file.path(OUTPUT_DIR, "Yeast_betas.svg"),
  plot = yeast_cor_plot,
  width = 20,
  height = 8,
  dpi = 100,
  units = "in"
)

if (nchar(HUMAN_SEED) > 0 && !is.null(p1$human_sparse)) {
  human_cor_plot <- (p1$human_sparse +
    p1$human_ridge +
    p1$human_plink +
    plot_layout(widths = c(44, 27.5, 28))) /
    (p2$human_sparse +
      p2$human_ridge +
      p2$human_plink +
      plot_layout(widths = c(44, 27.5, 28)))
  ggsave(
    file.path(OUTPUT_DIR, "Human_betas.svg"),
    plot = human_cor_plot,
    width = 20,
    height = 8,
    dpi = 100,
    units = "in"
  )
}

# ---------------------------------------------------------------------------
# Figures 8/9 — genomic effect size plots
# ---------------------------------------------------------------------------
YEAST_METHOD_LIST <- c(
  elasticnet_test = "Elastic Net",
  lasso_test = "LASSO",
  lars_maxiter1000_test = "LARS",
  pytorch_ridge_train = "Approx. Ridge",
  ridge_test = "Ridge",
  plink_clumped_train = "Clumped PLINK",
  plink_train = "PLINK"
)
HUMAN_METHOD_LIST <- c(
  elasticnet_test = "Elastic Net",
  lasso_test = "LASSO",
  lars_maxiter1000_test = "LARS",
  pytorch_ridge_train = "Approx. Ridge",
  ridge_test = "Ridge",
  plink_clumped_cov6_train = "Clumped PLINK",
  plink_cov6_train = "PLINK"
)
YEAST_LABELLER <- gsub(
  "cutoff",
  YEAST_P_THR,
  gsub("train", "train_n", gsub("test", "test_n", YEAST_METHOD_LIST))
)
HUMAN_LABELLER <- if (nchar(HUMAN_SEED) > 0) {
  gsub(
    "cutoff",
    HUMAN_P_THR,
    gsub("train", "train_n", gsub("test", "test_n", HUMAN_METHOD_LIST))
  )
} else {
  YEAST_LABELLER
}

rescale_region <- function(df, chrom = NULL, lower, upper) {
  out <- df |> filter(POS >= lower, POS <= upper)
  if (!is.null(chrom)) {
    out <- out |> filter(Chr == chrom)
  }
  out |>
    group_by(method_append1, trait) |>
    mutate(
      POS_kb = POS / 1000,
      group_max_abs = max(abs(coeff), na.rm = TRUE),
      group_max_truth = max(abs(truth_none), na.rm = TRUE),
      coeff_rescaled = if_else(
        group_max_abs == 0,
        coeff,
        coeff / group_max_abs
      ),
      truth_none_rescaled = if_else(
        group_max_truth == 0,
        truth_none,
        truth_none / group_max_truth
      )
    ) |>
    ungroup()
}

make_genomic_plot <- function(
  df,
  labeller_vec,
  x_label,
  title = NULL,
  hide_y_axis = FALSE
) {
  p <- ggplot(
    df,
    aes(
      y = coeff_rescaled,
      x = POS_kb,
      color = inter,
      alpha = inter,
      size = inter2
    )
  ) +
    geom_abline(slope = 0, intercept = 0) +
    geom_point(data = df %>% filter(truth_none_rescaled == 0), shape = 16) +
    geom_star(
      data = df %>% filter(truth_none_rescaled != 0),
      aes(x = POS_kb, y = truth_none_rescaled),
      fill = "white",
      color = "black",
      size = 3,
      inherit.aes = FALSE
    ) +
    geom_point(data = df %>% filter(truth_none_rescaled != 0), shape = 16) +
    facet_wrap(
      ~method_append1,
      ncol = 1,
      scales = "free_y",
      labeller = labeller(method_append1 = labeller_vec)
    ) +
    scale_size_manual(values = c("TRUE" = 3, "FALSE" = 1)) +
    scale_alpha_manual(
      values = c(
        elasticnet_test.TRUE = 0.8,
        lasso_test.TRUE = 0.8,
        pytorch_ridge_train.TRUE = 0.8,
        ridge_test.TRUE = 0.8,
        plink_clumped_train.TRUE = 0.8,
        plink_train.TRUE = 0.8,
        lars_maxiter1000_test.TRUE = 0.8,
        elasticnet_test.FALSE = 0.5,
        lasso_test.FALSE = 0.5,
        pytorch_ridge_train.FALSE = 0.5,
        ridge_test.FALSE = 0.5,
        plink_clumped_train.FALSE = 0.5,
        plink_train.FALSE = 0.5,
        lars_maxiter1000_test.FALSE = 0.5
      )
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(limits = c(-1.1, 1.1)) +
    scale_colour_manual(
      values = METHOD_COLOURS,
      name = "Method",
      labels = METHOD_LABELS
    ) +
    guides(color = "none", alpha = "none", size = "none", fill = "none") +
    theme_pub(x_axis_type = "numerical") +
    theme(panel.grid.major.x = element_blank()) +
    ylab("Rescaled Estimated Coefficient") +
    xlab(x_label)
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  if (hide_y_axis) {
    p <- p +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()
      )
  }
  p
}

yeast_trait <- "trait_overlap_vavg_0.98_numQTL_100_numChr_16_H2_0.5"
yeast_littlelonger_with_fullinfo$method_append1 <- factor(
  yeast_littlelonger_with_fullinfo$method_append,
  levels = names(YEAST_METHOD_LIST)
)
yeast_data_subset1 <- yeast_littlelonger_with_fullinfo[
  yeast_littlelonger_with_fullinfo$trait == yeast_trait &
    yeast_littlelonger_with_fullinfo$method_append1 %in%
      names(YEAST_METHOD_LIST)[-1],
]
print(head(yeast_data_subset1))

make_genomic_df <- function(df, chrom = NULL, lower, upper) {
  rescale_region(df, chrom = chrom, lower = lower, upper = upper) %>%
    mutate(
      has_truth = truth_none != 0,
      has_estimate = coeff != 0 & !is.na(coeff),
      inter = interaction(method_append, has_truth),
      inter2 = has_truth & has_estimate
    )
}

yeast_chr_11_subset <- make_genomic_plot(
  make_genomic_df(yeast_data_subset1, 11, 160000, 260000),
  YEAST_LABELLER,
  "Yeast Chromosome 11 Position (kb)"
)
yeast_chr_12_subset <- make_genomic_plot(
  make_genomic_df(yeast_data_subset1, 12, 700000, 1060000),
  YEAST_LABELLER,
  "Yeast Chromosome 12 Position (kb)",
  hide_y_axis = TRUE
)
yeast_genome <- yeast_chr_11_subset + yeast_chr_12_subset
ggsave(
  file.path(OUTPUT_DIR, "Yeast_genomic.svg"),
  plot = yeast_genome,
  width = 10,
  height = 10,
  dpi = 100,
  units = "in"
)

if (nchar(HUMAN_SEED) > 0 && exists("human_littlelonger_with_fullinfo")) {
  human_trait <- "trait_overlap_vavg_0.98_numQTL_50_numChr_1_H2_0.5"
  human_littlelonger_with_fullinfo$method_append1 <- factor(
    human_littlelonger_with_fullinfo$method_append,
    levels = names(HUMAN_METHOD_LIST)
  )
  human_data_subset1 <- human_littlelonger_with_fullinfo[
    human_littlelonger_with_fullinfo$trait == human_trait &
      human_littlelonger_with_fullinfo$method_append1 %in%
        names(HUMAN_METHOD_LIST),
  ]
  human_1_subset <- make_genomic_plot(
    make_genomic_df(human_data_subset1, NULL, 3.26e7, 3.26e7 + 2e6),
    HUMAN_LABELLER,
    "Human Chromosome 21 Position (kb)"
  )
  human_2_subset <- make_genomic_plot(
    make_genomic_df(human_data_subset1, NULL, 4.21e7, 4.21e7 + 2e6),
    HUMAN_LABELLER,
    "Human Chromosome 21 Position (kb)",
    hide_y_axis = TRUE
  )
  human_genome <- human_1_subset + human_2_subset
  ggsave(
    file.path(OUTPUT_DIR, "Human_genomic.svg"),
    plot = human_genome,
    width = 10,
    height = 10,
    dpi = 100,
    units = "in"
  )
}
# ---------------------------------------------------------------------------
# ROC — approximate (window-based)
# ---------------------------------------------------------------------------
make_roc_approx <- function(df, total_snps, cutoffs, window = 100) {
  cutoff_vec <- sort(unique(cutoffs))

  add_interval <- function(ivs, new_start, new_end) {
    if (nrow(ivs) == 0) {
      return(matrix(c(new_start, new_end), nrow = 1))
    }
    overlaps <- ivs[, 2] >= new_start - 1L & ivs[, 1] <= new_end + 1L
    if (!any(overlaps)) {
      insert_at <- findInterval(new_start, ivs[, 1]) + 1L
      before <- if (insert_at > 1L) {
        ivs[seq_len(insert_at - 1L), , drop = FALSE]
      } else {
        matrix(integer(0), ncol = 2)
      }
      after <- if (insert_at <= nrow(ivs)) {
        ivs[seq(insert_at, nrow(ivs)), , drop = FALSE]
      } else {
        matrix(integer(0), ncol = 2)
      }
      return(rbind(before, c(new_start, new_end), after))
    }
    merged_start <- min(new_start, min(ivs[overlaps, 1]))
    merged_end <- max(new_end, max(ivs[overlaps, 2]))
    keep <- ivs[!overlaps, , drop = FALSE]
    merged <- rbind(keep, c(merged_start, merged_end))
    merged[order(merged[, 1]), , drop = FALSE]
  }

  in_any_interval <- function(pos, ivs) {
    if (nrow(ivs) == 0) {
      return(rep(FALSE, length(pos)))
    }
    idx <- findInterval(pos, ivs[, 1])
    inside <- rep(FALSE, length(pos))
    has_match <- idx >= 1L
    inside[has_match] <- pos[has_match] <= ivs[idx[has_match], 2]
    inside
  }

  all_pos <- df %>%
    select(trait, POS, truth_none) %>%
    distinct()

  groups <- df %>% filter(!is.na(coeff)) %>% distinct(trait, method_append)
  n_groups <- nrow(groups)
  group_counter <- 0L
  pb <- txtProgressBar(min = 0, max = n_groups, style = 3)

  result <- df %>%
    filter(!is.na(coeff)) %>%
    arrange(trait, method_append, snp_in_order) %>%
    group_by(trait, method_append) %>%
    group_modify(function(grp, keys) {
      group_counter <- group_counter + 1L
      setTxtProgressBar(pb, group_counter)

      trait_pos <- all_pos %>% filter(trait == keys$trait)
      pos_vec <- trait_pos$POS
      is_causal <- trait_pos$truth_none != 0

      max_order <- max(grp$snp_in_order)
      valid_cutoffs <- cutoff_vec[cutoff_vec <= max_order]

      if (length(valid_cutoffs) == 0) {
        return(data.frame(
          cutoff = integer(0),
          unique_tp = integer(0),
          unique_fp = integer(0)
        ))
      }

      called <- grp %>% filter(coeff != 0) %>% arrange(snp_in_order)
      results <- vector("list", length(valid_cutoffs))
      ivs <- matrix(integer(0), nrow = 0, ncol = 2)
      called_idx <- 1L

      for (i in seq_along(valid_cutoffs)) {
        k <- valid_cutoffs[i]

        while (
          called_idx <= nrow(called) && called$snp_in_order[called_idx] <= k
        ) {
          p <- called$POS[called_idx]
          ivs <- add_interval(ivs, p - window, p + window)
          called_idx <- called_idx + 1L
        }

        if (nrow(ivs) == 0L) {
          results[[i]] <- data.frame(cutoff = k, unique_tp = 0L, unique_fp = 0L)
        } else {
          hits <- in_any_interval(pos_vec, ivs)
          results[[i]] <- data.frame(
            cutoff = k,
            unique_tp = sum(hits & is_causal),
            unique_fp = sum(hits & !is_causal)
          )
        }
      }

      bind_rows(results)
    }) %>%
    ungroup() %>%
    left_join(
      df %>% select(trait, numQTL) %>% distinct(),
      by = "trait"
    ) %>%
    mutate(
      tpr = unique_tp / numQTL,
      fpr = unique_fp / (total_snps - numQTL)
    ) %>%
    select(trait, method_append, cutoff, tpr, fpr)

  close(pb)
  result
}

# ---------------------------------------------------------------------------
# True additive effects — yeast
# ---------------------------------------------------------------------------
yeast_true_effects <- read_feather(file.path(
  BASE_DIR,
  paste0("intermediates_seed_", YEAST_SEED),
  paste0(
    YEAST_PREFIX,
    "_seed_",
    YEAST_SEED,
    "_additive_effects_wide.feather"
  )
)) %>%
  separate_wider_delim(
    SNP,
    names = c("SNP_clean", NA),
    delim = "_",
    too_few = "align_start"
  ) %>%
  rename(SNP = SNP_clean)

yeast_true_effects_longer <- yeast_true_effects %>%
  pivot_longer(cols = -SNP, names_to = "trait", values_to = "addEff") %>%
  filter(addEff != 0) %>%
  left_join(yeast_snp_list_sorted, by = "SNP")

yeast_true_effects_longer_nofilter <- yeast_true_effects %>%
  pivot_longer(cols = -SNP, names_to = "trait", values_to = "addEff") %>%
  left_join(yeast_snp_list_sorted, by = "SNP")


# ---------------------------------------------------------------------------
# SNP pairs (within-chromosome distances)
# ---------------------------------------------------------------------------
if (opt$`redo-snp-pairs`) {
  all_snps_yeast <- yeast_snp_list_sorted |> distinct(SNP, Chr, POS)

  snp_pairs_yeast <- all_snps_yeast |>
    group_by(Chr) |>
    group_split() |>
    map(\(chunk) {
      chunk <- arrange(chunk, POS)
      n <- nrow(chunk)
      forward <- map(seq_len(n), \(i) {
        j_max <- findInterval(chunk$POS[i] + 1000000, chunk$POS)
        if (j_max <= i) {
          return(NULL)
        }
        j_idx <- (i + 1):j_max
        tibble(
          SNP_a = chunk$SNP[i],
          SNP_b = chunk$SNP[j_idx],
          distance = chunk$POS[j_idx] - chunk$POS[i]
        )
      }) |>
        list_rbind()
      self <- tibble(SNP_a = chunk$SNP, SNP_b = chunk$SNP, distance = 0L)
      bind_rows(forward, self)
    }) |>
    list_rbind()

  arrow::write_feather(snp_pairs_yeast, INT_SNP_PAIRS_YEAST)
} else {
  snp_pairs_yeast <- read_feather(INT_SNP_PAIRS_YEAST)
}

# ---------------------------------------------------------------------------
# Nearest true-SNP distances — yeast
# ---------------------------------------------------------------------------
if (opt$`redo-distances`) {
  yeast_dist_methods <- c(
    "plink_train",
    "plink_clumped_train",
    "lasso_test",
    "lars_maxiter1000_test",
    "ridge_test",
    "elasticnet_test",
    "pytorch_ridge_train"
  )

  yeast_chunks <- yeast_littlelonger_with_fullinfo |>
    filter(
      !is.na(coeff),
      coeff != 0,
      trait %in% traits,
      method_append %in% yeast_dist_methods
    ) |>
    group_by(trait, method_append) |>
    group_split()

  pb <- cli_progress_bar(
    "Yeast: nearest true-SNP distances",
    total = length(yeast_chunks)
  )

  yeast_true_avg_distance <- yeast_chunks |>
    map(\(chunk) {
      this_trait <- chunk$trait[[1]]
      this_method <- chunk$method_append[[1]]

      tr <- yeast_true_effects_longer |> filter(trait == this_trait)

      if (this_method %in% c("ridge_test", "pytorch_ridge_train")) {
        threshold <- quantile(abs(chunk$coeff), 0.99, na.rm = TRUE)
        chunk <- chunk |> filter(abs(coeff) >= threshold)
      }

      result <- chunk |>
        left_join(snp_pairs_yeast, by = c("SNP" = "SNP_a")) |>
        left_join(tr, by = c("SNP_b" = "SNP", "trait")) |>
        filter(!is.na(addEff)) |>
        group_by(SNP) |>
        slice_min(distance, n = 1, with_ties = FALSE) |>
        group_by(SNP_b, trait, method_append) |>
        summarise(
          mean_dist = mean(distance),
          min_dist = min(distance),
          n_mapped = n(),
          .groups = "drop"
        )

      missing_true_snps <- tr |>
        filter(!SNP %in% result$SNP_b) |>
        transmute(
          SNP_b = SNP,
          trait = this_trait,
          method_append = this_method,
          mean_dist = NA,
          min_dist = NA,
          n_mapped = 0L
        )

      cli_progress_update(id = pb)
      bind_rows(result, missing_true_snps)
    }) |>
    list_rbind()

  cli_progress_done(id = pb)
  arrow::write_feather(yeast_true_avg_distance, INT_DIST_YEAST)

  # ── Human distances ──────────────────────────────────────────────────────
  if (nchar(HUMAN_SEED) > 0) {
    human_dist_methods <- c(
      "plink_cov6_train",
      "plink_clumped_cov6_train",
      "lasso_test",
      "lars_maxiter1000_test",
      "ridge_test",
      "elasticnet_test",
      "pytorch_ridge_train"
    )

    human_chunks <- human_littlelonger_with_fullinfo |>
      filter(
        !is.na(coeff),
        coeff != 0,
        trait %in% traits,
        method_append %in% human_dist_methods
      ) |>
      group_by(trait, method_append) |>
      group_split()

    pb <- cli_progress_bar(
      "Human: nearest true-SNP distances",
      total = length(human_chunks)
    )

    human_true_avg_distance <- human_chunks |>
      map(\(chunk) {
        this_trait <- chunk$trait[[1]]
        this_method <- chunk$method_append[[1]]

        tr <- human_true_effects_longer |> filter(trait == this_trait)

        if (this_method %in% c("ridge_test", "pytorch_ridge_train")) {
          threshold <- quantile(abs(chunk$coeff), 0.99, na.rm = TRUE)
          chunk <- chunk |> filter(abs(coeff) >= threshold)
        }

        result <- chunk |>
          left_join(snp_pairs_human, by = c("SNP" = "SNP_a")) |>
          left_join(tr, by = c("SNP_b" = "SNP", "trait")) |>
          filter(!is.na(addEff)) |>
          group_by(SNP) |>
          slice_min(distance, n = 1, with_ties = FALSE) |>
          group_by(SNP_b, trait, method_append) |>
          summarise(
            mean_dist = mean(distance),
            min_dist = min(distance),
            n_mapped = n(),
            .groups = "drop"
          )

        missing_true_snps <- tr |>
          filter(!SNP %in% result$SNP_b) |>
          transmute(
            SNP_b = SNP,
            trait = this_trait,
            method_append = this_method,
            mean_dist = NA,
            min_dist = NA,
            n_mapped = 0L
          )

        cli_progress_update(id = pb)
        bind_rows(result, missing_true_snps)
      }) |>
      list_rbind()

    cli_progress_done(id = pb)
    arrow::write_feather(human_true_avg_distance, INT_DIST_HUMAN)
  }
} else {
  yeast_true_avg_distance <- read_feather(INT_DIST_YEAST)
  if (nchar(HUMAN_SEED) > 0) {
    human_true_avg_distance <- read_feather(INT_DIST_HUMAN)
  }
}

# ---------------------------------------------------------------------------
# Distance ECDF / survival plot helper
# ---------------------------------------------------------------------------
plot_ecdf_survfit <- function(
  data,
  x_col,
  color_col,
  linetype_col,
  legend_pos,
  dist_plot_legend_col,
  subtitle = NULL,
  facet_col = NULL,
  x_offset = 1,
  tag,
  dashed_methods = NULL,
  show_x_title = TRUE
) {
  data <- data %>%
    mutate(
      .x = .data[[x_col]] + x_offset,
      .color = .data[[color_col]],
      .event = as.integer(!is.na(.data[[x_col]]))
    ) %>%
    mutate(.x = ifelse(is.na(.x), Inf, .x))

  if (!is.null(dashed_methods)) {
    dashed_data <- data %>% filter(.color %in% dashed_methods)
    data <- data %>% filter(!(.color %in% dashed_methods))
  }

  sf <- survfit2(Surv(.x, .event) ~ .color, data = data)

  p <- sf %>%
    ggsurvfit(linetype_aes = FALSE, alpha = 1) +
    labs(
      x = "Avg. Dist. to Closest True SNP (bp)",
      y = "Cumulative Prop. of True SNPs",
      color = color_col,
      subtitle = subtitle
    ) +
    guides(color = guide_legend(ncol = dist_plot_legend_col)) +
    theme_pub(x_axis_type = "numerical") +
    scale_y_reverse(
      limits = c(1, 0),
      breaks = seq(1, 0, by = -0.25),
      labels = c("0.00", "0.25", "0.50", "0.75", "1.00")
    ) +
    scale_x_log10(
      limits = c(1, 1000000),
      breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000),
      labels = c(
        0,
        expression(10^1),
        expression(10^2),
        expression(10^3),
        expression(10^4),
        expression(10^5),
        expression(10^6)
      )
    ) +
    scale_linetype_manual(values = c("1" = "solid", "2" = "dashed")) +
    guides(color = "none", linetype = "none") +
    scale_colour_manual(
      values = METHOD_COLOURS,
      name = "Method",
      labels = METHOD_LABELS
    ) +
    theme(
      legend.key.spacing.y = unit(0.3, "cm"),
      legend.key.spacing.x = unit(0.4, "cm"),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      legend.position = legend_pos,
      plot.subtitle = element_text(size = 15),
      legend.background = element_rect(fill = "transparent", color = NA)
    ) +
    labs(tag = tag)

  if (!is.null(dashed_methods)) {
    for (m in dashed_methods) {
      sf_m <- survfit2(
        Surv(.x, .event) ~ 1,
        data = dashed_data %>% filter(.color == m)
      )
      sf_tidy <- data.frame(time = sf_m$time, estimate = sf_m$surv)
      p <- p +
        geom_step(
          data = sf_tidy,
          aes(x = time, y = estimate),
          linetype = "dashed",
          color = METHOD_COLOURS[[m]]
        )
    }
  }

  if (!show_x_title) {
    p <- p + theme(axis.title.x = element_blank())
  }
  p
}

# ---------------------------------------------------------------------------
# ROC plot helpers
# ---------------------------------------------------------------------------

# Base ROC plot
plot_roc <- function(
  data,
  trait_id,
  expand,
  show_y_title = TRUE,
  show_y_ticks = TRUE,
  title = NULL,
  subtitle = NULL,
  tag,
  dashed_methods = NULL,
  show_x_title = TRUE
) {
  p <- ggplot(data[data$trait == trait_id, ])

  if (!is.null(dashed_methods)) {
    p <- p +
      geom_line(
        aes(x = fpr, y = tpr, color = method_append, linetype = linetype),
        alpha = 1
      )
  } else {
    p <- p +
      geom_line(
        aes(x = fpr, y = tpr, color = method_append),
        alpha = 1
      )
  }

  p <- p +
    theme_pub(x_axis_type = "numerical") +
    labs(x = "False Positive Rate", y = "True Positive Rate", tag = tag) +
    scale_x_continuous(limits = c(0, 1), expand = expand) +
    scale_y_continuous(limits = c(0, 1), expand = expand) +
    geom_abline(intercept = 0, slope = 1) +
    guides(color = "none", linetype = "none") +
    scale_colour_manual(
      values = METHOD_COLOURS,
      name = "Method",
      labels = METHOD_LABELS
    ) +
    theme(
      strip.text = element_text(size = 10),
      plot.subtitle = element_text(size = 14)
    )

  if (!show_y_title) {
    p <- p + theme(axis.title.y = element_blank())
  }
  if (!show_x_title) {
    p <- p + theme(axis.title.x = element_blank())
  }
  if (!show_y_ticks) {
    p <- p +
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
  }
  p
}

# Inset ROC plot (zoomed, transparent background)
plot_roc_inset <- function(
  data,
  trait_id,
  xlim,
  ylim,
  expand,
  show_y_title = FALSE,
  breaks = NULL,
  dashed_methods,
  show_x_title = TRUE
) {
  p <- plot_roc(
    data,
    trait_id,
    expand,
    show_y_title = show_y_title,
    tag = NULL,
    dashed_methods = dashed_methods,
    show_x_title = show_x_title
  ) +
    theme(
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent"),
      legend.box.background = element_rect(fill = "transparent"),
      axis.title.x = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    scale_colour_manual(
      values = METHOD_COLOURS,
      name = "Method",
      labels = METHOD_LABELS
    )

  if (!is.null(breaks)) {
    p <- p + scale_x_continuous(breaks = breaks, labels = breaks)
  }
  p
}

# ---------------------------------------------------------------------------
# Cumulative ROC data — exact match
# ---------------------------------------------------------------------------
make_roc_fast <- function(df, total_snps, cutoffs) {
  df %>%
    filter(!is.na(coeff), snp_in_order %in% cutoffs) %>%
    group_by(trait, method_append, snp_in_order, has_coeff) %>%
    summarize(
      tpr = first(cum_truepos) / first(numQTL),
      fpr = first(cum_falsepos) / (total_snps - first(numQTL)),
      .groups = "drop"
    ) %>%
    rename(cutoff = snp_in_order)
}

yeast_cutoffs <- c(
  seq(1, 99),
  seq(100, 2000, by = 10),
  seq(2000, 41500, by = 100),
  41594
)
human_cutoffs <- c(
  seq(1, 99),
  seq(100, 2000, by = 10),
  seq(2000, 11300, by = 100),
  11342
)

if (opt$`remake-roc`) {
  yeast_for_roc_all <- make_roc_fast(
    yeast_with_cumulative_other_small,
    41594,
    yeast_cutoffs
  )
  arrow::write_feather(yeast_for_roc_all, INT_ROC_YEAST)
} else {
  yeast_for_roc_all <- read_feather(INT_ROC_YEAST)
}

# ---------------------------------------------------------------------------
# Cumulative ROC data — approximate (window-based)
# ---------------------------------------------------------------------------
if (opt$`remake-roc-approx`) {
  yeast_for_roc_approx <- make_roc_approx(
    yeast_with_cumulative_other_small[
      yeast_with_cumulative_other_small$trait %in% traits,
    ],
    41594,
    yeast_cutoffs,
    window = ROC_APPROX_WINDOW
  )
  arrow::write_feather(yeast_for_roc_approx, INT_ROC_APPROX_YEAST)
} else {
  yeast_for_roc_approx <- read_feather(INT_ROC_APPROX_YEAST)
}

# ---------------------------------------------------------------------------
# ROC panel builder
# ---------------------------------------------------------------------------
make_roc_panel <- function(
  roc_data,
  roc_approx_data,
  dist_data,
  trait_id,
  xlim1,
  ylim1,
  xlim2,
  ylim2,
  expand,
  expand_inset,
  tags,
  title = NULL,
  subtitle1 = NULL,
  subtitle2 = NULL,
  subtitle3 = NULL,
  show_y_title_inset = FALSE,
  breaks = NULL,
  dist_plot_legend_pos,
  dist_plot_legend_col,
  dashed_methods = NULL,
  show_x_title = TRUE
) {
  large <- plot_roc(
    roc_data,
    trait_id,
    expand,
    title = title,
    subtitle = subtitle1,
    tag = tags[1],
    dashed_methods = dashed_methods,
    show_x_title = show_x_title
  )
  inset <- plot_roc_inset(
    roc_data,
    trait_id,
    xlim1,
    ylim1,
    expand_inset,
    show_y_title = show_y_title_inset,
    breaks = breaks,
    dashed_methods = dashed_methods,
    show_x_title = show_x_title
  )
  large_approx <- plot_roc(
    roc_approx_data,
    trait_id,
    expand,
    show_y_title = FALSE,
    show_y_ticks = FALSE,
    subtitle = subtitle2,
    tag = tags[2],
    dashed_methods = dashed_methods,
    show_x_title = show_x_title
  )
  inset_approx <- plot_roc_inset(
    roc_approx_data,
    trait_id,
    xlim2,
    ylim2,
    expand_inset,
    show_y_title = show_y_title_inset,
    breaks = breaks,
    dashed_methods = dashed_methods,
    show_x_title = show_x_title
  )
  dists <- plot_ecdf_survfit(
    data = dist_data %>% filter(trait == trait_id),
    x_col = "min_dist",
    color_col = "method_append",
    linetype_col = "linetype",
    legend_pos = dist_plot_legend_pos,
    subtitle = subtitle3,
    dist_plot_legend_col = dist_plot_legend_col,
    tag = tags[3],
    dashed_methods = dashed_methods,
    show_x_title = show_x_title
  )

  list(
    large = large,
    inset = inset,
    large_approx = large_approx,
    inset_approx = inset_approx,
    dists = dists
  )
}

# ---------------------------------------------------------------------------
# Figures 10/11 — ROC curves
# ---------------------------------------------------------------------------
methods_for_roc <- c(
  "elasticnet_test",
  "lasso_test",
  "lars_maxiter1000_test",
  "pytorch_ridge_train",
  "ridge_test",
  "plink_clumped_train",
  "plink_train",
  "plink_clumped_cov6_train",
  "plink_cov6_train"
)

# ── Trait definitions ──────────────────────────────────────────────────────────
yeast_trait1 <- "trait_overlap_vavg_0.9_numQTL_100_numChr_16_H2_0.25"
yeast_trait2 <- "trait_overlap_vavg_0.98_numQTL_100_numChr_16_H2_0.5"

yeast_trait1_xlim1 <- c(0, 0.02)
yeast_trait1_ylim1 <- c(0, 0.25)
yeast_trait2_xlim1 <- c(0, 0.01)
yeast_trait2_ylim1 <- c(0, 0.22)
yeast_trait1_xlim2 <- c(0, 0.075)
yeast_trait1_ylim2 <- c(0, 0.25)
yeast_trait2_xlim2 <- c(0, 0.075)
yeast_trait2_ylim2 <- c(0, 0.22)

human_trait1 <- "trait_overlap_vavg_0.9_numQTL_100_numChr_1_H2_0.25"
human_trait2 <- "trait_overlap_vavg_0.98_numQTL_100_numChr_1_H2_0.5"

human_trait1_xlim1 <- c(0, 0.016)
human_trait1_ylim1 <- c(0, 0.7)
human_trait2_xlim1 <- c(0, 0.016)
human_trait2_ylim1 <- c(0, 0.8)
human_trait1_xlim2 <- c(0, 0.016)
human_trait1_ylim2 <- c(0, 0.7)
human_trait2_xlim2 <- c(0, 0.016)
human_trait2_ylim2 <- c(0, 0.8)

# ── Linetype annotations ───────────────────────────────────────────────────────
yeast_for_roc_all$linetype <- "1"
yeast_for_roc_all$linetype[
  yeast_for_roc_all$method_append == "lasso_test"
] <- "2"
yeast_for_roc_approx$linetype <- "1"
yeast_for_roc_approx$linetype[
  yeast_for_roc_approx$method_append == "lasso_test"
] <- "2"

yeast_true_avg_distance$method_append <- factor(
  yeast_true_avg_distance$method_append,
  levels = methods_for_roc
)
yeast_true_avg_distance$linetype <- "1"
yeast_true_avg_distance$linetype[
  yeast_true_avg_distance$method_append == "lasso_test"
] <- "2"

# ── Build yeast ROC panels ─────────────────────────────────────────────────────
yeast_p1 <- make_roc_panel(
  yeast_for_roc_all[
    yeast_for_roc_all$method_append %in%
      methods_for_roc &
      yeast_for_roc_all$has_coeff,
  ],
  yeast_for_roc_approx[
    yeast_for_roc_approx$method_append %in% methods_for_roc,
  ],
  yeast_true_avg_distance[
    yeast_true_avg_distance$method_append %in% methods_for_roc,
  ],
  yeast_trait1,
  yeast_trait1_xlim1,
  yeast_trait1_ylim1,
  yeast_trait1_xlim2,
  yeast_trait1_ylim2,
  expand = c(0, 0.01),
  expand_inset = expansion(mult = c(0, 0), add = c(0.001, 0)),
  show_y_title_inset = FALSE,
  dist_plot_legend_pos = c(0.3, 0.75),
  dist_plot_legend_col = 2,
  tags = c("A)", "B)", "C)"),
  show_x_title = FALSE
)

yeast_p2 <- make_roc_panel(
  yeast_for_roc_all[
    yeast_for_roc_all$method_append %in%
      methods_for_roc &
      yeast_for_roc_all$has_coeff,
  ],
  yeast_for_roc_approx[
    yeast_for_roc_approx$method_append %in% methods_for_roc,
  ],
  yeast_true_avg_distance[
    yeast_true_avg_distance$method_append %in% methods_for_roc,
  ],
  yeast_trait2,
  yeast_trait2_xlim1,
  yeast_trait2_ylim1,
  yeast_trait2_xlim2,
  yeast_trait2_ylim2,
  expand = c(0, 0.01),
  expand_inset = expansion(mult = c(0, 0), add = c(0.001, 0)),
  dist_plot_legend_pos = c(0.3, 0.75),
  dist_plot_legend_col = 2,
  tags = c("D)", "E)", "F)"),
  show_x_title = TRUE
)

# ── Compose and save Figure 10 (yeast) ────────────────────────────────────────
inset_pos <- list(left = 0.47, bottom = 0.01, right = 1.02, top = 0.51)

add_inset <- function(large, inset, inset_pos) {
  large +
    inset_element(
      inset,
      left = inset_pos$left,
      bottom = inset_pos$bottom,
      right = inset_pos$right,
      top = inset_pos$top
    )
}

yeast_roc <- (add_inset(yeast_p1$large, yeast_p1$inset, inset_pos) +
  add_inset(yeast_p1$large_approx, yeast_p1$inset_approx, inset_pos) +
  yeast_p1$dists) /
  (add_inset(yeast_p2$large, yeast_p2$inset, inset_pos) +
    add_inset(yeast_p2$large_approx, yeast_p2$inset_approx, inset_pos) +
    yeast_p2$dists)

ggsave(
  file.path(OUTPUT_DIR, "Yeast_ROC.svg"),
  plot = yeast_roc,
  width = 18,
  height = 10,
  dpi = 100,
  units = "in"
)

# ---------------------------------------------------------------------------
# Figure 12 — summary heatmaps
# ---------------------------------------------------------------------------

# ── Method ordering / labels (local to summary figure) ────────────────────────
summary_method_labels <- c(
  lars_maxiter1000_test = "LARS",
  lasso_test = "LASSO",
  pytorch_ridge_train = "Approx.\nRidge",
  ridge_test = "Ridge",
  elasticnet_test = "Elastic\nNet",
  plink_train = "PLINK",
  plink_clumped_train = "Clumped\nPLINK",
  plink_cov6_train = "PLINK",
  plink_clumped_cov6_train = "Clumped\nPLINK",
  lars_maxiter1000 = "LARS",
  lasso = "LASSO",
  pytorch_ridge = "Approx.\nRidge",
  ridge = "Ridge",
  elasticnet = "Elastic\nNet",
  plink = "PLINK",
  plink_clumped = "Clumped\nPLINK"
)

summary_method_order <- c(
  "elasticnet_test",
  "lasso_test",
  "pytorch_ridge_train",
  "ridge_test",
  "lars_maxiter1000_test",
  "plink_clumped_train",
  "plink_train",
  "plink_clumped_cov6_train",
  "plink_cov6_train",
  "elasticnet",
  "lasso",
  "pytorch_ridge",
  "ridge",
  "lars_maxiter1000",
  "plink_clumped",
  "plink"
)

summary_method_order_yeast <- c(
  "elasticnet_test",
  "lasso_test",
  "ridge_test",
  "pytorch_ridge_train",
  "plink_clumped_train",
  "plink_train",
  "plink_clumped_cov6_train",
  "plink_cov6_train",
  "lars_maxiter1000_test",
  "elasticnet",
  "lasso",
  "ridge",
  "pytorch_ridge",
  "plink_clumped",
  "plink",
  "lars_maxiter1000"
)

# ── Summary data preparation ───────────────────────────────────────────────────
yeast_mean_of_means <- expand_trait_names(
  yeast_true_avg_distance %>%
    group_by(trait, method_append) %>%
    replace(is.na(.), 10^6) %>%
    summarize(m = mean(mean_dist, na.rm = TRUE), .groups = "drop")
) %>%
  group_by(target_H2, method_append) %>%
  summarize(mean_by_tile = mean(m, na.rm = TRUE), .groups = "drop")

yeast_mean_of_means$method_append <- factor(
  yeast_mean_of_means$method_append,
  levels = summary_method_order_yeast
)

yeast_method_subset_pred <- all_pred_cors_with_truth[
  all_pred_cors_with_truth$seed == "Yeast" &
    all_pred_cors_with_truth$overlap != "singlesnp" &
    ((all_pred_cors_with_truth$method == "plink" &
      all_pred_cors_with_truth$split == "train") |
      (all_pred_cors_with_truth$method == "plink_clumped" &
        all_pred_cors_with_truth$split == "train") |
      all_pred_cors_with_truth$method %in%
        c("lars_maxiter1000", "lasso", "ridge", "pytorch_ridge")) &
    (all_pred_cors_with_truth$cov %in%
      c("", "cov6_") |
      is.na(all_pred_cors_with_truth$cov)),
]
yeast_method_subset_eff <- all_eff_cors_with_truth[
  all_eff_cors_with_truth$seed == "Yeast" &
    all_eff_cors_with_truth$overlap != "singlesnp" &
    ((all_eff_cors_with_truth$method == "plink" &
      all_eff_cors_with_truth$split == "train") |
      (all_eff_cors_with_truth$method == "plink_clumped" &
        all_eff_cors_with_truth$split == "train") |
      all_eff_cors_with_truth$method %in%
        c("lars_maxiter1000", "lasso", "ridge", "pytorch_ridge")) &
    (all_eff_cors_with_truth$cov %in%
      c("", "cov6_") |
      is.na(all_eff_cors_with_truth$cov)),
]

yeast_pred <- yeast_method_subset_pred %>%
  group_by(target_H2, method) %>%
  summarize(mean_by_tile = median(r2, na.rm = TRUE), .groups = "drop")

yeast_eff <- yeast_method_subset_eff %>%
  group_by(target_H2, method) %>%
  summarize(mean_by_tile = median(r2, na.rm = TRUE), .groups = "drop")

yeast_pred$method <- factor(
  yeast_pred$method,
  levels = summary_method_order_yeast
)
yeast_eff$method <- factor(
  yeast_eff$method,
  levels = summary_method_order_yeast
)

if (nchar(HUMAN_SEED) > 0) {
  human_mean_of_means <- expand_trait_names(
    human_true_avg_distance %>%
      group_by(trait, method_append) %>%
      replace(is.na(.), 10^6) %>%
      summarize(m = mean(mean_dist, na.rm = TRUE), .groups = "drop")
  ) %>%
    group_by(target_H2, method_append) %>%
    summarize(mean_by_tile = mean(m, na.rm = TRUE), .groups = "drop")

  human_mean_of_means$method_append <- factor(
    human_mean_of_means$method_append,
    levels = summary_method_order
  )

  human_method_subset_pred <- all_pred_cors_with_truth[
    all_pred_cors_with_truth$seed == "Human" &
      all_pred_cors_with_truth$overlap != "singlesnp" &
      ((all_pred_cors_with_truth$method == "plink" &
        all_pred_cors_with_truth$split == "train") |
        (all_pred_cors_with_truth$method == "plink_clumped" &
          all_pred_cors_with_truth$split == "train") |
        all_pred_cors_with_truth$method %in%
          c(
            "elasticnet",
            "lars_maxiter1000",
            "lasso",
            "ridge",
            "pytorch_ridge"
          )) &
      (all_pred_cors_with_truth$cov %in%
        c("", "cov6_") |
        is.na(all_pred_cors_with_truth$cov)),
  ]
  human_method_subset_eff <- all_eff_cors_with_truth[
    all_eff_cors_with_truth$seed == "Human" &
      all_eff_cors_with_truth$overlap != "singlesnp" &
      ((all_eff_cors_with_truth$method == "plink" &
        all_eff_cors_with_truth$split == "train") |
        (all_eff_cors_with_truth$method == "plink_clumped" &
          all_eff_cors_with_truth$split == "train") |
        all_eff_cors_with_truth$method %in%
          c(
            "elasticnet",
            "lars_maxiter1000",
            "lasso",
            "ridge",
            "pytorch_ridge"
          )) &
      (all_eff_cors_with_truth$cov %in%
        c("", "cov6_") |
        is.na(all_eff_cors_with_truth$cov)),
  ]

  human_pred <- human_method_subset_pred %>%
    group_by(target_H2, method) %>%
    summarize(mean_by_tile = median(r2, na.rm = TRUE), .groups = "drop")

  human_eff <- human_method_subset_eff %>%
    group_by(target_H2, method) %>%
    summarize(mean_by_tile = median(r2, na.rm = TRUE), .groups = "drop")

  human_pred$method <- factor(human_pred$method, levels = summary_method_order)
  human_eff$method <- factor(human_eff$method, levels = summary_method_order)
}

# ── Shared heatmap theme additions ────────────────────────────────────────────
theme_heatmap_shared <- function(hide_x = FALSE, hide_y_text = FALSE) {
  t <- theme_pub(x_axis_type = "discrete") +
    theme(
      plot.margin = margin(b = 0, r = 0, unit = "pt"),
      legend.direction = "horizontal",
      legend.title.position = "top"
    )
  if (hide_x) {
    t <- t +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank()
      )
  }
  if (hide_y_text) {
    t <- t + theme(axis.title.y = element_blank())
  }
  t
}

# ── Build heatmap panels ───────────────────────────────────────────────────────
make_heatmap <- function(
  data,
  x_col,
  fill_col,
  tag,
  fill_label,
  palette_name,
  limits,
  breaks = waiver(),
  fill_labels = waiver(),
  reverse = FALSE,
  hide_x = FALSE,
  hide_y_text = FALSE,
  hide_y_axis = FALSE,
  show_legend = TRUE
) {
  p <- ggplot(
    data,
    aes(
      x = .data[[x_col]],
      y = as.factor(target_H2),
      fill = .data[[fill_col]]
    )
  ) +
    geom_tile() +
    theme_pub(x_axis_type = "discrete") +
    gradient_fill_arcadia(
      palette_name = palette_name,
      reverse = reverse,
      limits = limits,
      breaks = breaks,
      labels = fill_labels
    ) +
    scale_x_discrete(labels = summary_method_labels, expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme(
      plot.margin = margin(b = 0, r = 0, unit = "pt"),
      legend.direction = "horizontal",
      legend.title.position = "top"
    ) +
    labs(fill = fill_label, tag = tag, y = "Broad-Sense Heritability")

  if (!show_legend) {
    p <- p + guides(fill = "none")
  }
  if (hide_x) {
    p <- p +
      theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  }
  if (hide_y_text) {
    p <- p + theme(axis.title.y = element_blank())
  }
  if (hide_y_axis) {
    p <- p +
      theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank()
      )
  }
  p
}

# ── Yeast panels ───────────────────────────────────────────────────────────────
yeast_predictive_plot <- make_heatmap(
  yeast_pred %>% filter(target_H2 != 1, method != "elasticnet"),
  x_col = "method",
  fill_col = "mean_by_tile",
  tag = "D",
  fill_label = expression(italic(r)^2),
  palette_name = "magma",
  limits = c(0, 1),
  reverse = TRUE,
  hide_x = TRUE,
  hide_y_axis = TRUE,
  show_legend = TRUE
) +
  theme(legend.text = element_text(angle = 45, vjust = 1, hjust = 1))

yeast_effect_plot <- make_heatmap(
  yeast_eff %>% filter(target_H2 != 1, method != "elasticnet"),
  x_col = "method",
  fill_col = "mean_by_tile",
  tag = "E",
  fill_label = expression(italic(r)^2),
  palette_name = "magma",
  limits = c(0, 1),
  reverse = TRUE,
  hide_x = TRUE,
  hide_y_axis = TRUE,
  show_legend = TRUE
) +
  theme(legend.text = element_text(angle = 45, vjust = 1, hjust = 1))

yeast_dist_cdf_plot <- make_heatmap(
  yeast_mean_of_means %>% filter(method_append != "elasticnet_test"),
  x_col = "method_append",
  fill_col = "mean_by_tile",
  tag = "F",
  fill_label = "Avg. Dist.",
  palette_name = "magma",
  limits = c(5.1, 6),
  breaks = log10(c(150000, 250000, 500000, 1000000)),
  fill_labels = c("150kb", "250kb", "500kb", "1Mb"),
  reverse = FALSE,
  hide_y_text = TRUE,
  hide_y_axis = TRUE,
  show_legend = TRUE
) +
  aes(fill = log10(mean_by_tile)) +
  theme(
    legend.box.spacing = unit(0.8, "cm"),
    legend.text = element_text(angle = 45, vjust = 1, hjust = 1)
  ) +
  labs(x = "Method")

if (nchar(HUMAN_SEED) > 0) {
  # ── Human panels ──────────────────────────────────────────────────────────
  human_predictive_plot <- make_heatmap(
    human_pred %>% filter(target_H2 != 1),
    x_col = "method",
    fill_col = "mean_by_tile",
    tag = "A",
    fill_label = expression(italic(r)^2),
    palette_name = "magma",
    limits = c(0, 1),
    reverse = TRUE,
    hide_x = TRUE,
    hide_y_text = TRUE,
    show_legend = FALSE
  )

  human_effect_plot <- make_heatmap(
    human_eff %>% filter(target_H2 != 1),
    x_col = "method",
    fill_col = "mean_by_tile",
    tag = "B",
    fill_label = expression(italic(r)^2),
    palette_name = "magma",
    limits = c(0, 1),
    reverse = TRUE,
    hide_x = TRUE,
    hide_y_text = FALSE,
    show_legend = FALSE
  )
  human_dist_cdf_plot <- make_heatmap(
    human_mean_of_means,
    x_col = "method_append",
    fill_col = "mean_by_tile",
    tag = "C",
    fill_label = "Avg. Dist.",
    palette_name = "magma",
    limits = c(5.1, 6),
    hide_y_text = TRUE,
    show_legend = FALSE
  ) +
    aes(fill = log10(mean_by_tile)) +
    labs(x = "Method")

  summary_fig <- (human_predictive_plot + yeast_predictive_plot) /
    (human_effect_plot + yeast_effect_plot) /
    (human_dist_cdf_plot + yeast_dist_cdf_plot)
} else {
  summary_fig <- yeast_predictive_plot /
    yeast_effect_plot /
    yeast_dist_cdf_plot
}

ggsave(
  file.path(OUTPUT_DIR, "YeastHuman_Summary.svg"),
  plot = summary_fig,
  width = 14,
  height = 10,
  dpi = 100,
  units = "in"
)

# ---------------------------------------------------------------------------
# Supplementary Figure 1 — LD decay
# ---------------------------------------------------------------------------
INT_LD_YEAST <- int_path("yeast_ld_subset_", YEAST_SEED)
INT_LD_HUMAN <- int_path("human_ld_subset_", HUMAN_SEED)

LD_SAMPLE_SIZE <- 1000000L
LD_MAX_DIST <- 750000L
LD_BIN_WIDTH <- c(10000, 0.04)

yeast_ld_file <- file.path(BASE_DIR, "yeast_plink_ld", "yeast_ld.vcor")
yeast_ld <- fread(yeast_ld_file) %>%
  mutate(dist = abs(POS_A - POS_B), species = "Yeast") %>%
  filter(dist <= LD_MAX_DIST)

set.seed(as.integer(YEAST_SEED))
yeast_ld_subset <- yeast_ld[sample(nrow(yeast_ld), size = LD_SAMPLE_SIZE), ]

yeast_ld_breakdown <- ggplot(
  yeast_ld_subset,
  aes(x = dist, y = UNPHASED_R2)
) +
  geom_hex(binwidth = LD_BIN_WIDTH) +
  theme_pub(x_axis_type = "numerical") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.01))) +
  gradient_fill_arcadia(
    palette_name = "purples",
    trans = "log10",
    breaks = c(1, 100, 10000),
    labels = c(expression(1), expression(10^2), expression(10^4))
  ) +
  labs(fill = "# of Pairs", tag = "A") +
  theme(
    legend.position = c(0.88, 0.72),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.direction = "horizontal",
    legend.title.position = "top"
  )

ld_smooth <- ggplot(
  yeast_ld_subset,
  aes(x = dist, y = UNPHASED_R2, color = species, linetype = species)
) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(
    breaks = seq(0, LD_MAX_DIST, by = 200000),
    labels = c(0, 200, 400, 600),
    expand = expansion(mult = c(0.02, 0.01))
  ) +
  geom_smooth(se = TRUE) +
  labs(
    x = "Distance (kb)",
    y = expression("               " ~ r^2 ~ "between pairs of variants"),
    color = "Species",
    linetype = "Species",
    tag = "C"
  ) +
  theme_pub(x_axis_type = "numerical") +
  scale_color_manual(values = c("#7A77AB", "#97CD78")) +
  theme(
    legend.position = c(0.88, 0.72) #,axis.title.y = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(fill = NA)))

ld_plot <- yeast_ld_breakdown / ld_smooth

ggsave(
  file.path(OUTPUT_DIR, "Yeast_LD.svg"),
  plot = ld_plot,
  width = 10,
  height = 6,
  dpi = 100,
  units = "in"
)
