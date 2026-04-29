library(tidyr, quietly = TRUE, warn.conflicts = FALSE)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(fst, quietly = TRUE, warn.conflicts = FALSE)
library(data.table, quietly = TRUE, warn.conflicts = FALSE)
library(argparser, quietly = TRUE, warn.conflicts = FALSE)
library(arrow, quietly = TRUE, warn.conflicts = FALSE)

# generate_phenotypes.r — Simulate quantitative traits from a genotype matrix.
#
# Run from the repository root. Example usage:
#
#   Yeast — simulate traits:
#     Rscript scripts/generate_phenotypes.r \
#       --snp_file input_data/SNP_list_pos_corrected.txt \
#       --geno_dir input_data --geno yeast_genotypes_binarized.feather \
#       --subset_geno 100 --ploidy 1 --file_save_prefix yeast_test \
#       --single_snp_trait FALSE --remake_traits TRUE --calculate_phenos TRUE \
#       --vavg_ratios 1 --QTL_numbers 100 --broad_sense 0.9 0.5 --seed 1237
#
#   Yeast — save train/test splits:
#     Rscript scripts/generate_phenotypes.r \
#       --snp_file input_data/SNP_list_pos_corrected.txt \
#       --geno_dir input_data --geno yeast_genotypes_binarized \
#       --ploidy 1 --file_save_prefix yeast_simulated_data \
#       --single_snp_trait TRUE --remake_traits FALSE --calculate_phenos FALSE \
#       --save_test_train TRUE --seed 1510 \
#       --allele_freqs yeast_simulated_data_allele_frequencies.Rda --normalize TRUE
#
#   UK Biobank — simulate traits:
#     Rscript scripts/generate_phenotypes.r \
#       --snp_file genotypes/ukb22418_c21_b0_v2.bim \
#       --geno_dir genotypes --geno ukb22418_c21_b0_v2.raw \
#       --subset_geno 100 --ploidy 2 --file_save_prefix ukbb_test \
#       --single_snp_trait FALSE --remake_traits TRUE --calculate_phenos TRUE \
#       --vavg_ratios 1 --QTL_numbers 100 --broad_sense 0.9 0.5 --seed 1237 \
#       --biobank TRUE --allele_freqs ukb22418_c21_b0_v2.afreq
#
#   UK Biobank — save train/test splits:
#     Rscript scripts/generate_phenotypes.r \
#       --snp_file genotypes/ukb22418_c21_b0_v2.bim \
#       --geno_dir genotypes --geno ukb22418_c21_b0_v2 \
#       --ploidy 2 --file_save_prefix ukbb_simulated_traits \
#       --single_snp_trait TRUE --remake_traits FALSE --calculate_phenos FALSE \
#       --save_test_train TRUE --seed 1105 \
#       --allele_freqs ukb22418_c21_b0_v2.afreq --normalize TRUE --biobank TRUE

p <- arg_parser("generate_phenotypes")

p <- add_argument(
  p,
  "--biobank",
  help = "Whether or not we are using UK Biobank (which has different file formats).",
  default = FALSE,
  nargs = 1,
  type = "logical"
)
p <- add_argument(
  p,
  "--geno",
  help = "Genotype file name. Provide either an uncentered genotype matrix file ending in '.feather' or a file prefix. If a prefix is provided, the script will use 'prefix_uncentered.feather','prefix_centered.feather' and 'prefix_ids.Rda'",
  default = "yeast_genotypes_binarized_subset_5000",
  nargs = 1,
  type = "character"
)
p <- add_argument(
  p,
  "--geno_dir",
  help = "Folder containing genotype files and allele frequency file.",
  default = "input_data/",
  nargs = 1,
  type = "character"
)
p <- add_argument(
  p,
  "--subset_geno",
  help = "Randomly subset the genotype file to this many individuals.",
  nargs = 1,
  type = "integer"
)
p <- add_argument(
  p,
  "--geno_seed",
  help = "Seed to use to subset the genotype file.",
  nargs = 1,
  type = "integer"
)
p <- add_argument(
  p,
  "--subset_snps",
  help = "Randomly subset the genotype matrix to this many SNPs (columns).",
  nargs = 1,
  type = "integer"
)
p <- add_argument(
  p,
  "--snp_file",
  help = "Path to SNP information file.",
  default = "input_data/SNP_list_pos_corrected.txt",
  nargs = 1,
  type = "character"
)
p <- add_argument(
  p,
  "--allele_freqs",
  help = "Allele frequency file. If not provided, allele frequencies will be calculated from genotype file.",
  nargs = 1,
  type = "character"
)
p <- add_argument(
  p,
  "--output_dir",
  help = "Output directory for intermediates created during this analysis.",
  nargs = 1,
  type = "character"
)
p <- add_argument(
  p,
  "--file_save_prefix",
  help = "File save prefix for intermediates created during this analysis",
  default = "all_simulated_traits",
  nargs = 1,
  type = "character"
)
p <- add_argument(
  p,
  "--test_train_output_dir",
  help = "Output directory for final genotypes and phenotypes created during this analysis.",
  nargs = 1,
  type = "character"
)
p <- add_argument(
  p,
  "--single_snp_trait",
  help = "Make a trait with a single causal snp.",
  default = TRUE,
  nargs = 1,
  type = "logical"
)
p <- add_argument(
  p,
  "--remake_single_snp_trait",
  help = "Whether to remake (TRUE) or just reload the single-snp trait (FALSE).",
  default = FALSE,
  nargs = 1,
  type = "logical"
)
p <- add_argument(
  p,
  "--remake_traits",
  help = "Whether to remake (TRUE) or just reload the traits (FALSE).",
  default = FALSE,
  nargs = 1,
  type = "logical"
)
p <- add_argument(
  p,
  "--calculate_phenos",
  help = "Use the generated trait architecture to calculate phenotypes for each individual.",
  default = TRUE,
  nargs = 1,
  type = "logical"
)
p <- add_argument(
  p,
  "--save_test_train",
  help = "Whether to save the test and train data as csv and feather files.",
  default = TRUE,
  nargs = 1,
  type = "logical"
)
p <- add_argument(
  p,
  "--seed",
  help = "Random seed for multi-snp trait generation.",
  default = 520,
  nargs = 1,
  type = "integer"
)
p <- add_argument(
  p,
  "--single_snp_seed",
  help = "Random seed for single-snp trait generation.",
  nargs = 1,
  type = "integer"
)
p <- add_argument(
  p,
  "--ploidy",
  help = "Ploidy of organism.",
  nargs = 1,
  type = "numeric"
)
p <- add_argument(
  p,
  "--vavg_ratios",
  help = "Comma-separated ratios of additive variance to total genetic variance (e.g. '0.5,0.9,1').",
  default = "1,0.98,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1",
  nargs = 1,
  type = "character"
)
p <- add_argument(
  p,
  "--QTL_numbers",
  help = "Comma-separated numbers of total QTLs to use (e.g. '50,100,500').",
  default = "50,100,500,1000,5000,10000,41594",
  nargs = 1,
  type = "character"
)
p <- add_argument(
  p,
  "--chromosome_numbers",
  help = "Numbers of chromosomes to spread QTLs over (list, default: all of them).",
  nargs = Inf,
  type = "numeric"
)
p <- add_argument(
  p,
  "--broad_sense",
  help = "Comma-separated broad-sense heritability values (e.g. '0.25,0.5,0.75').",
  default = "0.25,0.5,0.75,0.99,1",
  nargs = 1,
  type = "character"
)
p <- add_argument(
  p,
  "--train_prop",
  help = "Proportion of data to use in training split.",
  default = 0.85,
  nargs = 1,
  type = "numeric"
)
p <- add_argument(
  p,
  "--normalize",
  help = "Whether or not to normalize the phenotypes to have mean 0 and variance 1.",
  default = TRUE,
  nargs = 1,
  type = "logical"
)

args <- parse_args(p)

args$vavg_ratios <- as.numeric(strsplit(args$vavg_ratios, ",")[[1]])
args$QTL_numbers <- as.numeric(strsplit(args$QTL_numbers, ",")[[1]])
args$broad_sense  <- as.numeric(strsplit(args$broad_sense,  ",")[[1]])

file_prefix <- strsplit(args$geno, "\\.")[[1]][1]
if (is.na(args$ploidy)) {
  message("The organism's ploidy is required.")
  quit(save = "no")
}
if (is.na(args$snp_file)) {
  message("SNP information file is required.")
  quit(save = "no")
}
if (is.na(args$output_dir)) {
  args$output_dir <- paste0("intermediates_seed_", args$seed)
}
if (!endsWith(args$output_dir, "/")) {
  args$output_dir <- paste0(args$output_dir, "/")
}
if (!endsWith(args$geno_dir, "/")) {
  args$geno_dir <- paste0(args$geno_dir, "/")
}
if (is.na(args$test_train_output_dir)) {
  args$test_train_output_dir <- paste0("test_train_seed_", args$seed)
}
if (is.na(args$geno_seed)) {
  args$geno_seed <- args$seed + 1
}
if (is.na(args$single_snp_seed)) {
  args$single_snp_seed <- args$seed + 2
}
if (!dir.exists(args$output_dir)) {
  dir.create(args$output_dir, recursive = TRUE)
}
if (!dir.exists(args$test_train_output_dir)) {
  dir.create(args$test_train_output_dir, recursive = TRUE)
}

cat("==========================================================\n")
cat("               ANALYSIS CONFIGURATION SUMMARY            \n")
cat("==========================================================\n")
cat("\n--- Genotype Data ---\n")
cat(sprintf("  Genotype file/prefix:         %s\n", args$geno))
cat(sprintf("  Genotype directory:           %s\n", args$geno_dir))
cat(sprintf(
  "  Subset individuals:           %s\n",
  ifelse(is.na(args$subset_geno), "No subsetting", args$subset_geno)
))
cat(sprintf(
  "  Subset SNPs:                  %s\n",
  ifelse(is.na(args$subset_snps), "No subsetting", args$subset_snps)
))
cat("\n--- Input Files ---\n")
cat(sprintf("  SNP file:                     %s\n", args$snp_file))
cat(sprintf(
  "  Allele frequencies:           %s\n",
  ifelse(is.na(args$allele_freqs), "Not provided", args$allele_freqs)
))
cat("\n--- Output Files ---\n")
cat(sprintf("  Output directory:             %s\n", args$output_dir))
cat(sprintf("  File save prefix:             %s\n", args$file_save_prefix))
cat(sprintf("  Test/train output dir:        %s\n", args$test_train_output_dir))
cat("\n--- Actions ---\n")
cat(sprintf("  Single SNP trait:             %s\n", args$single_snp_trait))
cat(sprintf(
  "  Remake single SNP trait:      %s\n",
  args$remake_single_snp_trait
))
cat(sprintf("  Remake traits:                %s\n", args$remake_traits))
cat(sprintf("  Calculate phenotypes:         %s\n", args$calculate_phenos))
cat("\n--- Trait Generation Parameters ---\n")
cat(sprintf("  Random seed:                  %s\n", args$seed))
cat(sprintf("  Ploidy:                       %s\n", args$ploidy))
cat(sprintf(
  "  Va/Vg ratios:                 %s\n",
  paste(args$vavg_ratios, collapse = ", ")
))
cat(sprintf(
  "  QTL numbers:                  %s\n",
  paste(args$QTL_numbers, collapse = ", ")
))
cat(sprintf(
  "  Chromosome numbers:           %s\n",
  ifelse(
    is.null(args$chromosome_numbers) | all(is.na(args$chromosome_numbers)),
    "All",
    paste(args$chromosome_numbers, collapse = ", ")
  )
))
cat(sprintf(
  "  Broad-sense heritabilities:   %s\n",
  paste(args$broad_sense, collapse = ", ")
))
cat("\n--- Test/Train Split ---\n")
cat(sprintf("  Training proportion:      %s\n", args$train_prop))
cat("==========================================================\n\n")

# ==============================================================================
# FUNCTIONS
# ==============================================================================

# # Functions to distribute QTLs among chromosomes
# distribute_QTLs_evenly <- function(snps_per_chrom,mainChr,nQTLs) {
#   # Distributes QTLs evenly across specified chromosomes.
#   # Parameters:
#   # - nChr (integer): Total number of chromosomes in the genome
#   # - mainChr (vector): Chromosome numbers where QTLs should be placed
#   # - nQTLs (integer): Total number of QTLs to distribute
#   # Returns:
#   # A vector of length nChr with the number of QTLs assigned to each
#   # chromosome. Chromosomes not in mainChr receive 0 QTLs.

#   if (length(mainChr)>1) {
#     even_split<-table(cut(1:nQTLs, breaks = length(mainChr), labels = 1:nrow(snps_per_chrom)))
#   } else {
#     even_split <- nQTLs
#   }
#   even_split_randomized<-sample(even_split,length(even_split))
#   names(even_split_randomized)<-snps_per_chrom$Chromosome
#   return(even_split_randomized)
# }
# Functions to distribute QTLs among chromosomes
distribute_QTLs_evenly <- function(snps_per_chrom, mainChr, nQTLs) {
  # Distributes QTLs evenly across specified chromosomes.
  # Parameters:
  # - nChr (integer): Total number of chromosomes in the genome
  # - mainChr (vector): Chromosome numbers where QTLs should be placed
  # - nQTLs (integer): Total number of QTLs to distribute
  # Returns:
  # A vector of length nChr with the number of QTLs assigned to each
  # chromosome. Chromosomes not in mainChr receive 0 QTLs.

  if (length(mainChr) > 1) {
    even_split <- table(cut(
      1:nQTLs,
      breaks = length(mainChr),
      labels = 1:nrow(snps_per_chrom)
    ))
    even_split_randomized <- sample(even_split, length(even_split))
  } else {
    even_split_randomized <- nQTLs
  }

  names(even_split_randomized) <- snps_per_chrom$Chromosome
  return(even_split_randomized)
}

distribute_QTLs_proportionally <- function(snps_per_chrom, nQTLs) {
  # Distributes QTLs proportionally across all chromosomes based on
  # the real SNP density.
  # Parameters:
  # - snps_per_chrom (vector): Number of SNPs on each chromosome
  # - nQTLs (integer): Total number of QTLs to distribute
  # Returns:
  # A vector with the number of QTLs assigned to each chromosome,
  # proportional to SNP density. The function ensures the exact number
  # of QTLs is distributed through random adjustment if rounding causes
  # discrepancies.

  uneven_split <- round(snps_per_chrom / sum(snps_per_chrom) * nQTLs)
  while (sum(uneven_split) != nQTLs) {
    rand_chr <- sample(1:length(snps_per_chrom), 1)
    if ((sum(uneven_split) > nQTLs)) {
      uneven_split[rand_chr] <- uneven_split[rand_chr] - 1
    } else {
      uneven_split[rand_chr] <- uneven_split[rand_chr] + 1
    }
  }
  return(uneven_split)
}

choose_snps <- function(
  snp_map,
  loci_per_chr,
  relAA,
  trait_name,
  ploidy,
  mean = 0,
  var = 1,
  epistatic_overlap = TRUE,
  epistatic_loci_per_chr = NULL,
  choose_single_snp = FALSE,
  single_snp_effect = 10
) {
  # Selects SNPs to act as QTLs and assigns additive and epistatic effects.
  # Parameters:
  # - snp_map (data.frame): SNP information with columns including
  #   Chromosome and SNP
  # - loci_per_chr (vector or integer): Number of loci per chromosome.
  #   If a single value, applied to all chromosomes
  # - relAA (numeric): Relative variance of epistatic effects compared
  #   to additive effects
  # - trait_name (character): Name identifier for the trait
  # - mean (numeric): Mean of the effect size distribution (default: 0)
  # - var (numeric): Variance of the effect size distribution (default: 1)
  # - epistatic_overlap (logical): Whether epistatic pairs should be drawn
  #   from the previously chosen additive loci. Otherwise they are also
  #   chosen randomly. (default: TRUE)
  # - epistatic_loci_per_chr (vector or NULL): Number of epistatic loci per
  #   chromosome. If NULL, uses loci_per_chr
  # - choose_single_snp (logical): Modification of this function allowing
  #   selection of a single SNP with large effect rather than a large set
  #   (default: FALSE)
  # - single_snp_effect (numeric): Effect size when using single SNP mode
  #   (default: 10)
  # Returns:
  # A list containing:
  # - additive: Data frame of selected SNPs with additive effects (addEff
  #   column)
  # - epistatic: Data frame of epistatic pairs with interaction effects
  #   (epiEff column)
  # - mean: Effect size distribution mean
  # - var: Effect size distribution variance
  # Details:
  # - Additive effects are drawn from a normal distribution with specified
  #   mean and variance
  # - Epistatic effects are drawn from a normal distribution scaled by
  #   sqrt(relAA)
  # - When epistatic_overlap = TRUE, epistatic pairs are formed from the
  #   additive loci
  # - When epistatic_overlap = FALSE, epistatic loci are independently
  #   sampled
  # - The function pairs epistatic loci randomly and requires an even number
  #   of loci

  if (choose_single_snp) {
    additive_chosen <- snp_map %>% slice_sample(n = 1)
    additive_chosen$addEff <- single_snp_effect
    additive_chosen$name <- trait_name
    return(list(
      additive = as.data.frame(additive_chosen),
      snp_effect = single_snp_effect
    ))
  } else {
    if (length(loci_per_chr) == 1) {
      nChr <- length(unique(snp_map$Chromosome))
      loci_per_chr <- rep(loci_per_chr, nChr)
    }
    if (is.null(epistatic_loci_per_chr)) {
      epistatic_loci_per_chr <- loci_per_chr
    }
    additive_chosen <- snp_map %>%
      group_by(Chromosome) %>%
      slice_sample(prop = 1) %>%
      filter(row_number() <= loci_per_chr[cur_group_id()])
    additive_chosen$addEff <- rnorm(
      nrow(additive_chosen),
      mean = mean,
      sd = sqrt(var)
    )
    additive_chosen$name <- trait_name

    if (epistatic_overlap) {
      stopifnot(nrow(additive_chosen) %% 2 == 0)

      snps <- sample(additive_chosen$SNP, length(additive_chosen$SNP) / 2)
      epistatic_1 <- additive_chosen[additive_chosen$SNP %in% snps, ]
      epistatic_2 <- additive_chosen[!additive_chosen$SNP %in% snps, ]
    } else {
      snps_with_adds <- left_join(snp_map, additive_chosen)
      epistatic_chosen <- snps_with_adds %>%
        group_by(Chromosome) %>%
        slice_sample(prop = 1) %>%
        filter(row_number() <= epistatic_loci_per_chr[cur_group_id()])
      epistatic_chosen$addEff[is.na(epistatic_chosen$addEff)] <- 0

      snps <- sample(epistatic_chosen$SNP, length(epistatic_chosen$SNP) / 2)
      epistatic_1 <- epistatic_chosen[epistatic_chosen$SNP %in% snps, ]
      epistatic_2 <- epistatic_chosen[!epistatic_chosen$SNP %in% snps, ]
    }
    colnames(epistatic_2) <- paste0(colnames(epistatic_1), "_interactor")
    epistatic <- cbind(epistatic_1, epistatic_2)

    epistatic$epiEff <- rnorm(nrow(epistatic), mean = mean, sd = sqrt(var)) *
      sqrt(relAA * 2 * ploidy)
    epistatic$name <- trait_name
    return(list(
      additive = as.data.frame(additive_chosen),
      epistatic = as.data.frame(epistatic),
      mean = mean,
      var = var
    ))
  }
}

get_allele_freq <- function(centered_vector, ploidy) {
  # Calculates the frequency of the reference allele for a given locus.
  # Parameters:
  # - x (vector): Vector of raw genotype dosages for a single locus
  # - ploidy (integer): Ploidy level of the organism (e.g., 2 for diploid)
  # Returns:
  # A numeric value representing the reference allele frequency (p) at
  # the locus, estimated from genotype frequencies.
  # Details:
  # Raw genotypes are first centered and scaled to {-1, 0, 1} encoding.
  # Allele frequency is then estimated as: p = (freq_0 + ploidy * freq_1) / ploidy

  geno_freqs <- table(centered_vector) / length(centered_vector)

  freq_neg1 <- geno_freqs["-1"]
  freq_0 <- geno_freqs["0"]
  freq_1 <- geno_freqs["1"]

  if (is.na(freq_neg1)) {
    freq_neg1 <- 0
  }
  if (is.na(freq_0)) {
    freq_0 <- 0
  }
  if (is.na(freq_1)) {
    freq_1 <- 0
  }

  p <- (freq_0 + (ploidy) * freq_1) / ploidy

  return(p)
}

calculate_alphas <- function(
  first_name,
  second_name,
  g1,
  g2,
  first_snp_add,
  second_snp_add,
  epi_effect
) {
  # Calculates average allelic substitution effects (alphas) for an
  # epistatic SNP pair, accounting for the interaction between the two loci.
  # Parameters:
  # - first_name (character): Column name of the first SNP
  # - second_name (character): Column name of the second SNP
  # - g1 (vector): Centered genotype vector for the first SNP ({-1,0,1})
  # - g2 (vector): Centered genotype vector for the second SNP ({-1,0,1})
  # - first_snp_add (numeric): Additive effect of the first SNP
  # - second_snp_add (numeric): Additive effect of the second SNP
  # - epi_effect (numeric): Epistatic interaction effect
  # Returns:
  # A data frame with columns 'coeff' and 'SNP' containing the adjusted
  # alpha for each of the two SNPs.

  n <- length(g1)

  # Joint genotype frequencies — integer index avoids paste + table
  # Maps (-1,0,1) x (-1,0,1) → indices 1..9
  joint_idx <- (g1 + 2L) * 3L + (g2 + 2L) + 1L
  jf <- tabulate(joint_idx, nbins = 9L) / n

  # Genotype values laid out row-major (first SNP changes slowest):
  # idx: 1=-1/-1, 2=-1/0, 3=-1/1, 4=0/-1, 5=0/0, 6=0/1, 7=1/-1, 8=1/0, 9=1/1
  gv <- c(
    -first_snp_add - second_snp_add + epi_effect,
    -first_snp_add,
    -first_snp_add + second_snp_add - epi_effect,
    -second_snp_add,
    0,
    second_snp_add,
    first_snp_add - second_snp_add - epi_effect,
    first_snp_add,
    first_snp_add + second_snp_add + epi_effect
  )

  pop_mean <- sum(jf * gv)
  gvc <- gv - pop_mean

  # Marginal allele frequencies via tabulate
  mf1 <- tabulate(g1 + 2L, nbins = 3L) / n # [-1, 0, 1]
  mf2 <- tabulate(g2 + 2L, nbins = 3L) / n

  # Marginal means via matrix multiply — avoids 12 named scalar lookups
  gvc_mat <- matrix(gvc, nrow = 3, byrow = TRUE) # 3x3, rows=first, cols=second
  mm1 <- as.vector(gvc_mat %*% mf2) # marginal mean for first SNP
  mm2 <- as.vector(t(gvc_mat) %*% mf1) # marginal mean for second SNP

  add_exp <- outer(mm1, mm2, `+`) # 3x3 additive expectation
  deviations <- gvc_mat - add_exp # 3x3 interaction deviations

  # Analytic no-intercept OLS replaces lm() on a 9-row data frame
  first_vals <- rep(c(-1, 0, 1), each = 3)
  second_vals <- rep(c(-1, 0, 1), 3)
  dev_vec <- as.vector(t(deviations)) # row-major to match layout

  first_fixed <- all(g1 == -1)
  second_fixed <- all(g2 == -1)

  if (!first_fixed && !second_fixed) {
    ss1 <- sum(first_vals^2)
    ss2 <- sum(second_vals^2)
    ss12 <- sum(first_vals * second_vals)
    sy1 <- sum(first_vals * dev_vec)
    sy2 <- sum(second_vals * dev_vec)
    det <- ss1 * ss2 - ss12^2
    coeff <- c(
      (ss2 * sy1 - ss12 * sy2) / det,
      (ss1 * sy2 - ss12 * sy1) / det
    )
  } else if (first_fixed && !second_fixed) {
    coeff <- c(
      first_snp_add,
      sum(second_vals * dev_vec) / sum(second_vals^2)
    )
  } else if (second_fixed && !first_fixed) {
    coeff <- c(
      sum(first_vals * dev_vec) / sum(first_vals^2),
      second_snp_add
    )
  } else {
    coeff <- c(first_snp_add, second_snp_add)
  }

  names(coeff) <- c(first_name, second_name)
  data.frame(coeff = coeff, SNP = c(first_name, second_name))
}

breeding_vals_one_trait <- function(
  this_trait_epistatic_effects,
  this_trait_additive_effects,
  this_trait_additive_effects_matrix,
  centered_genotypes,
  snp_order,
  ploidy,
  allele_freq_matrix
) {
  # Calculates epistatic scores, alpha corrections, and breeding values
  # for a single trait.
  # Parameters:
  # - this_trait_epistatic_effects (data.frame): Epistatic SNP pairs and
  #   their interaction effects (columns: SNP, SNP_interactor, epiEff)
  # - this_trait_additive_effects (data.frame): Additive effects for SNPs
  #   (columns: SNP, addEff)
  # - this_trait_additive_effects_matrix (vector): Per-SNP additive effects
  #   aligned to snp_order
  # - centered_genotypes (matrix): Genotype matrix in {-1,0,1} encoding
  # - snp_order (vector): Ordered SNP names (columns of genotype matrix)
  # - ploidy (integer): Ploidy level
  # - allele_freq_matrix (matrix): Reference allele frequencies, individuals
  #   x SNPs
  # Returns:
  # A list with:
  # - epistatic_scores: Per-individual summed epistatic effects
  # - alphas: Alpha values corrected for epistatic contributions
  # - breeding_values: Per-individual breeding values

  n_pairs <- nrow(this_trait_epistatic_effects)
  first_names <- this_trait_epistatic_effects$SNP
  second_names <- this_trait_epistatic_effects$SNP_interactor
  epi_effs <- this_trait_epistatic_effects$epiEff

  # Named lookup avoids scanning the additive effects data frame each iteration
  add_lookup <- setNames(
    this_trait_additive_effects$addEff,
    this_trait_additive_effects$SNP
  )

  # Extract all needed genotype columns once up front
  g1_mat <- centered_genotypes[, first_names, drop = FALSE]
  g2_mat <- centered_genotypes[, second_names, drop = FALSE]

  # Epistatic scores: vectorized across all pairs at once
  epi_scores_mat <- sweep(g1_mat * g2_mat, 2, epi_effs, `*`)
  summed_epi <- rowSums(epi_scores_mat)
  rm(epi_scores_mat) # free immediately — summed_epi is all we need

  # Alpha corrections: lapply over pairs, passing pre-extracted vectors
  all_alpha_corrections <- lapply(seq_len(n_pairs), function(j) {
    calculate_alphas(
      first_names[j],
      second_names[j],
      g1_mat[, j],
      g2_mat[, j],
      add_lookup[first_names[j]],
      add_lookup[second_names[j]],
      epi_effs[j]
    )
  })

  # Free large column matrices as soon as calculate_alphas is done with them
  rm(g1_mat, g2_mat)
  gc()

  all_alphas_df <- do.call(rbind, all_alpha_corrections)
  alphas_ordered <- all_alphas_df[match(snp_order, all_alphas_df$SNP), ]
  alphas_ordered$coeff[is.na(alphas_ordered$coeff)] <- 0
  alphas_ordered$SNP[is.na(alphas_ordered$SNP)] <- snp_order[is.na(
    alphas_ordered$SNP
  )]

  alpha_vals <- this_trait_additive_effects_matrix - alphas_ordered$coeff

  alpha_matrix <- matrix(
    alpha_vals,
    nrow = nrow(centered_genotypes),
    ncol = ncol(centered_genotypes),
    byrow = TRUE
  )
  breeding_values <- as.data.frame(rowSums(
    (centered_genotypes + 1 - ploidy * allele_freq_matrix) * alpha_matrix
  ))
  rm(alpha_matrix)

  list(
    epistatic_scores = data.frame(name = summed_epi),
    alphas = alpha_vals,
    breeding_values = breeding_values
  )
}

calculate_breeding_vals <- function(
  do_epistatic,
  traits,
  epistatic,
  additive,
  additive_matrix,
  centered_geno,
  snp_order,
  ploidy,
  allele_freq_matrix,
  suffix,
  column_names,
  save_name
) {
  # Orchestrates calculation or loading of epistatic scores and breeding
  # values across one or multiple traits.
  # Parameters:
  # - do_epistatic (logical): If TRUE, calculates and saves epistatic scores
  #   and breeding values. If FALSE, loads previously saved results from disk
  # - traits (character vector): Trait name(s) to process
  # - epistatic (data.frame): Epistatic SNP pairs and interaction effects
  #   for all traits
  # - additive (data.frame): Additive SNP effects for all traits
  # - additive_matrix (matrix): Per-SNP additive effects as a matrix aligned
  #   to the SNP order of the genotype matrix
  # - centered_geno (matrix): Genotype matrix centered to {-1, 0, 1} encoding
  # - snp_order (vector): Ordered SNP names corresponding to genotype matrix
  #   columns
  # - ploidy (integer): Ploidy level of the organism
  # - allele_freq_matrix (matrix): Matrix of reference allele frequencies,
  #   with SNPs as columns duplicated down a number of rows equal to the number of
  #   individuals
  # - suffix (character): File name suffix indicating scaling status
  #   (e.g., "unscaled" or "scaled")
  # - column_names (character vector): Ordering of trait names to match output
  #   data frames to so that they can be added together later
  # - save_name (character): Base file path for saving or loading .fst files
  # Returns:
  # A list containing:
  # - epistatic_scores_ordered: Data frame of per-individual epistatic scores,
  #   columns ordered by column_names
  # - breeding_values_ordered: Data frame of per-individual breeding values,
  #   columns ordered by column_names
  # - alphas: Alpha correction values (single trait) or list of alpha vectors
  #   (multiple traits)

  if (do_epistatic) {
    if (is.vector(traits) & length(traits) > 1) {
      epistatic_all <- list()
      breeding_values_all <- list()
      alphas <- list()
      for (i in 1:length(traits)) {
        trait <- traits[i]
        this_epistatic <- epistatic[epistatic$name == trait, ]
        this_additive <- additive[additive$name == trait, ]
        this_additive_matrix <- additive_matrix[, trait]

        this_trait_results <- breeding_vals_one_trait(
          this_epistatic,
          this_additive,
          this_additive_matrix,
          centered_geno,
          snp_order,
          ploidy,
          allele_freq_matrix
        )
        epistatic_all[[i]] <- this_trait_results$epistatic_scores
        breeding_values_all[[i]] <- this_trait_results$breeding_values
        alphas[[i]] <- this_trait_results$alphas
      }
      epistatic_scores_all_traits <- as.data.frame(do.call(
        cbind,
        epistatic_all
      ))
      breeding_values_all_traits <- as.data.frame(do.call(
        cbind,
        breeding_values_all
      ))
    } else {
      this_epistatic <- epistatic[epistatic$name == traits, ]
      this_additive <- additive[additive$name == traits, ]
      this_additive_matrix <- additive_matrix[, traits]

      this_trait_results <- breeding_vals_one_trait(
        this_epistatic,
        this_additive,
        this_additive_matrix,
        centered_geno,
        snp_order,
        ploidy,
        allele_freq_matrix
      )

      epistatic_scores_all_traits <- data.frame(
        name = this_trait_results$epistatic_scores
      )
      breeding_values_all_traits <- data.frame(
        name = this_trait_results$breeding_values
      )
      alphas <- this_trait_results$alphas
    }
    stopifnot(length(alphas) == nrow(this_additive_matrix))
    colnames(epistatic_scores_all_traits) <- traits
    colnames(breeding_values_all_traits) <- traits
    epistatic_scores_ordered <- as.data.frame(epistatic_scores_all_traits) %>%
      select(all_of(column_names))
    breeding_values_ordered <- as.data.frame(breeding_values_all_traits) %>%
      select(all_of(column_names))
    write_fst(
      x = epistatic_scores_ordered,
      path = paste0(save_name, "_epistatic_", suffix, ".fst")
    )
    write_fst(
      x = breeding_values_ordered,
      path = paste0(save_name, "_breeding_values_", suffix, ".fst")
    )
  } else {
    epistatic_scores_ordered <- read_fst(paste0(
      save_name,
      "_epistatic_",
      suffix,
      ".fst"
    ))
    breeding_values_ordered <- read_fst(paste0(
      save_name,
      "_breeding_values_",
      suffix,
      ".fst"
    ))
  }
  result <- list(
    epistatic_scores_ordered = epistatic_scores_ordered,
    breeding_values_ordered = breeding_values_ordered,
    alphas = alphas
  )
  return(result)
}

calculate_traits <- function(
  snp_map,
  uncentered_geno_matrix,
  allele_freq_matrix,
  additive,
  epistatic = NULL,
  traits,
  ploidy,
  rescale_scores = TRUE,
  save_name = NULL,
  mean = 0,
  var = 1,
  do_epistatic = FALSE,
  centered_geno_precomputed = NULL
) {
  # Calculates genetic values for one or more traits by combining additive
  # and (optionally) epistatic effects, with optional rescaling to target
  # mean and variance.
  # Parameters:
  # - snp_map (data.frame): SNP information with at least columns SNP and
  #   Chromosome
  # - uncentered_geno_matrix (matrix): Raw genotype matrix (individuals x SNPs)
  #   in dosage encoding
  # - allele_freq_matrix (matrix): Matrix of reference allele frequencies,
  #   with SNPs as columns duplicated down a number of rows equal to the number of
  #   individuals
  # - additive (data.frame): Additive SNP effects with columns SNP, name
  #   (trait), and addEff
  # - epistatic (data.frame or NULL): Epistatic SNP pair effects with columns
  #   SNP, SNP_interactor, name (trait), and epiEff. Set to NULL to skip
  #   epistatic calculations (useful for single-snp case)
  # - traits (character vector): Trait name(s) to calculate
  # - ploidy (integer): Ploidy level of the organism
  # - rescale_scores (logical): If TRUE, rescales effects to achieve the
  #   target mean and variance (default: TRUE). The "FALSE" option should be used
  #   when the effect sizes given are already scaled, as is the case when we are
  #   trying to replicate a result from AlphaSimR
  # - save_name (character): Base file path for saving intermediate .fst files
  # - mean (numeric or vector): Target mean(s) of the genetic value distribution
  #   (default: 0)
  # - var (numeric or vector): Target variance(s) of the genetic value
  #   distribution (default: 1)
  # - do_epistatic (logical): If TRUE, calculates epistatic scores from
  #   scratch. If FALSE, loads previously saved results (default: FALSE)
  # - centered_geno_precomputed (matrix or NULL): Optional pre-centered genotype
  #   matrix. If provided, skips the centering step (saves ~31 GB allocation and
  #   ~14s per call when looping over many traits with the same genotype matrix).
  # Returns:
  # If epistatic effects are included, a list containing:
  # - additive: (Rescaled) additive effect data frame
  # - epistatic: (Rescaled) epistatic effect data frame
  # - genetic_values: Data frame of final per-individual genetic values
  #   for each trait
  # - add_scores: Matrix of per-individual additive scores for each trait
  # - epi_scores: Data frame of per-individual epistatic scores for each trait
  # - intercept: Data frame of per-trait intercept corrections with number of
  #   identical rows equal to the number of individuals
  # - breeding_vals: Data frame of per-individual breeding values for each trait
  # - alphas: Corrected effect sizes for each SNP (used to calculate breeding
  #   values)
  # - mean: Target mean value(s)
  # - var: Target variance value(s)
  # If no epistatic effects, a list with:
  # - additive: Additive effect data frame
  # - genetic_values: Matrix of additive scores

  if (rescale_scores) {
    suffix1 <- "unscaled"
    suffix2 <- "scaled"
  } else {
    suffix1 <- "scaled"
  }

  # center genotypes — skip if a pre-centered matrix was passed in
  if (!is.null(centered_geno_precomputed)) {
    centered_geno <- centered_geno_precomputed
  } else {
    centered_geno <- (uncentered_geno_matrix - (ploidy / 2)) * (2 / ploidy)
  }

  # Subset to traits of interest and format additive effects as matrix
  epistatic <- epistatic[epistatic$name %in% traits, ]

  additive_pivoted <- pivot_wider(
    additive[additive$name %in% traits, ],
    id_cols = "SNP",
    names_from = "name",
    values_from = "addEff",
    values_fill = 0
  ) %>%
    arrange(match(SNP, colnames(centered_geno)))

  additive_pivoted2 <- left_join(
    snp_map[, c("SNP", "Chromosome")],
    additive_pivoted,
    by = "SNP"
  ) %>%
    mutate_if(is.numeric, coalesce, 0)

  stopifnot(all(colnames(centered_geno) == additive_pivoted2$SNP))
  additive_matrix <- as.matrix(additive_pivoted2 %>% select(-SNP, -Chromosome))

  # Calculate additive value for each individual
  message("Calculating Additive Scores...")

  additive_scores <- centered_geno %*% additive_matrix
  additive_scores_df <- as.data.frame(additive_scores)
  write_fst(
    x = additive_scores_df,
    path = paste0(save_name, "_additive_", suffix1, ".fst")
  )

  if (!is.null(epistatic)) {
    # Calculate epistatic value for each individual
    message("Calculating Epistatic Scores...")

    bv_result <- calculate_breeding_vals(
      do_epistatic,
      traits,
      epistatic,
      additive,
      additive_matrix,
      centered_geno,
      additive_pivoted2$SNP,
      ploidy,
      allele_freq_matrix,
      suffix1,
      colnames(additive_scores),
      save_name
    )
    epistatic_scores_ordered <- bv_result$epistatic_scores_ordered
    breeding_values_ordered <- bv_result$breeding_values_ordered
    alphas <- bv_result$alphas

    # Calculate genetic values for rescaling
    stopifnot(all(
      colnames(additive_scores) == colnames(epistatic_scores_ordered)
    ))
    genetic_values <- additive_scores + epistatic_scores_ordered

    intercept <- mean - sapply(genetic_values, mean)
    intercept_df <- as.data.frame(matrix(
      ncol = length(intercept),
      nrow = nrow(centered_geno),
      data = as.numeric(intercept),
      byrow = TRUE
    ))
    colnames(intercept_df) <- traits

    genetic_values_new <- genetic_values + intercept_df

    if (rescale_scores) {
      # calculate scale factors
      stopifnot(
        (length(mean) == 1 | length(mean) == length(traits)),
        (length(var) == 1 | length(var) == length(traits))
      )
      scale_vals <- sqrt(var) /
        sapply(breeding_values_ordered, function(x) {
          sqrt(var(x) * ((length(x) - 1) / length(x)))
        })
      scales_df <- as.data.frame(matrix(
        ncol = length(scale_vals),
        nrow = nrow(centered_geno),
        data = as.numeric(scale_vals),
        byrow = TRUE
      ))

      scales <- data.frame(
        scale = scale_vals,
        name = colnames(additive_scores),
        intercept = intercept
      )

      # Rescale the effect sizes
      epistatic_merged <- left_join(epistatic, scales, by = "name") %>%
        dplyr::mutate(scaled_epiEff = epiEff * scale) %>%
        select(-epiEff, -scale) %>%
        rename(epiEff = scaled_epiEff)

      additive_merged <- left_join(
        additive[additive$name %in% traits, ],
        scales,
        by = "name"
      ) %>%
        dplyr::mutate(scaled_addEff = addEff * scale) %>%
        select(-addEff, -scale) %>%
        rename(addEff = scaled_addEff)

      additive_pivoted_neweff <- pivot_wider(
        additive_merged,
        id_cols = "SNP",
        names_from = "name",
        values_from = "addEff",
        values_fill = 0
      ) %>%
        arrange(match(SNP, colnames(centered_geno)))

      additive_pivoted_neweff2 <- left_join(
        snp_map[, c("SNP", "Chromosome")],
        additive_pivoted_neweff,
        by = "SNP"
      ) %>%
        mutate_if(is.numeric, coalesce, 0)

      # Calculate additive and epistatic scores for each individual using rescaled values
      message("Calculating Rescaled Additive Scores...")

      stopifnot(all(colnames(centered_geno) == additive_pivoted_neweff2$SNP))
      additive_matrix_neweff <- as.matrix(
        additive_pivoted_neweff2 %>% select(-SNP, -Chromosome)
      )

      additive_scores_new <- centered_geno %*% additive_matrix_neweff
      additive_scores_new_df <- as.data.frame(additive_scores_new)
      write_fst(
        x = additive_scores_new_df,
        path = paste0(save_name, "_additive_", suffix2, ".fst")
      )

      message("Calculating Rescaled Epistatic Scores...")

      bv_result_rescaled <- calculate_breeding_vals(
        do_epistatic,
        traits,
        epistatic_merged,
        additive_merged,
        additive_matrix_neweff,
        centered_geno,
        additive_pivoted2$SNP,
        ploidy,
        allele_freq_matrix,
        suffix2,
        colnames(additive_scores_new),
        save_name
      )
      epistatic_scores_rescaled_ordered <- bv_result_rescaled$epistatic_scores_ordered
      stopifnot(all(
        colnames(additive_scores_new) ==
          colnames(epistatic_scores_rescaled_ordered)
      ))
      genetic_values_rescaled <- additive_scores_new +
        epistatic_scores_rescaled_ordered

      intercept_new <- mean - sapply(genetic_values_rescaled, mean)
      intercept_new_df <- as.data.frame(matrix(
        ncol = length(intercept_new),
        nrow = nrow(centered_geno),
        data = as.numeric(intercept_new),
        byrow = TRUE
      ))
      colnames(intercept_new_df) <- traits

      genetic_values_rescaled_new <- genetic_values_rescaled + intercept_new_df

      return(list(
        additive = additive_merged,
        epistatic = epistatic_merged,
        genetic_values = genetic_values_rescaled_new,
        add_scores = additive_scores_new,
        epi_scores = epistatic_scores_rescaled_ordered,
        intercept = intercept_new_df,
        breeding_vals = breeding_values_ordered,
        alphas = alphas,
        mean = mean,
        var = var
      ))
    } else {
      return(list(
        additive = additive,
        epistatic = epistatic,
        genetic_values = genetic_values_new,
        add_scores = additive_scores,
        epi_scores = epistatic_scores_ordered,
        intercept = intercept_df,
        breeding_vals = breeding_values_ordered,
        alphas = alphas,
        mean = mean,
        var = var
      ))
    }
  } else {
    return(list(additive = additive, genetic_values = additive_scores))
  }
}


add_environment <- function(genetic_values, broad_sense, varG) {
  # Adds environmental noise to genetic values to achieve a target
  # broad-sense heritability.
  # Parameters:
  # - genetic_values (matrix): Matrix of genetic values (individuals x traits)
  # - broad_sense (numeric): Target broad-sense heritability (H2)
  # - varG (numeric vector): Genetic variance for each trait
  # Returns:
  # - A matrix of phenotypic values with environmental noise added. Column
  #   names include the heritability suffix (e.g., _H2_0.5).

  nTraits <- ncol(genetic_values)
  nInd <- nrow(genetic_values)

  stopifnot(length(broad_sense) == 1)
  name_suffix <- paste0("_H2_", broad_sense)
  broad_sense <- rep(broad_sense, nTraits)

  stopifnot(length(broad_sense) == nTraits)
  varE <- numeric(nTraits)
  for (i in seq_len(nTraits)) {
    tmp <- varG[i] / broad_sense[i] - varG[i]
    varE[i] <- tmp
  }

  # Create phenotypes
  stopifnot(length(varE) == nTraits)

  error <- lapply(varE, function(x) {
    if (is.na(x)) {
      return(rep(NA_real_, nInd))
    } else {
      return(rnorm(nInd, sd = sqrt(x)))
    }
  })
  error <- do.call("cbind", error)

  pheno <- genetic_values + error
  colnames(pheno) <- paste0(colnames(genetic_values), name_suffix)

  return(pheno)
}
combine_files <- function(
  file_path,
  file_name,
  file_type,
  mode,
  save_path,
  save_name = NULL,
  save_name_gv = NULL,
  save_name_others = NULL,
  write_csv = FALSE
) {
  # Merges intermediate output files produced by the trait loop.
  #
  # Two modes:
  #   "fst"          - reads a set of .fst files and cbinds or rbinds them
  #                    into one combined .fst (and optionally a .csv).
  #   "trait_results"- loads a set of per-trait .Rda files, extracts the
  #                    standard saved objects, concatenates them, and writes
  #                    a combined genetic-values .fst and a scaled-effects .Rda.
  #
  # Parameters:
  # - file_path (character): Directory to search for input files
  # - file_name (character): Regex pattern passed to list.files()
  # - file_type (character): "fst" or "trait_results"
  # - mode (character): "cbind" or "rbind" (only used for file_type = "fst")
  # - save_path (character): Directory for output files
  # - save_name (character): Output filename stem for "fst" mode
  # - save_name_gv (character): Output filename stem for genetic values in
  #   "trait_results" mode
  # - save_name_others (character): Output filename stem for the .Rda of scaled
  #   effects in "trait_results" mode
  # - write_csv (logical): If TRUE, also writes a .csv in "fst" mode
  input_files <- list.files(file_path, pattern = file_name, full.names = TRUE)
  if (length(input_files) == 0) {
    stop(sprintf(
      "combine_files: no files matched pattern '%s' in '%s'",
      file_name,
      file_path
    ))
  }

  if (file_type == "fst") {
    # ---- load ----------------------------------------------------------------
    data_list <- vector("list", length(input_files))
    for (i in seq_along(input_files)) {
      message(sprintf(
        "  Loading file %d / %d: %s",
        i,
        length(input_files),
        basename(input_files[i])
      ))
      data_list[[i]] <- read_fst(input_files[i])
    }

    # ---- combine -------------------------------------------------------------
    combined_df <- switch(
      mode,
      cbind = as.data.frame(do.call(cbind, data_list)),
      rbind = as.data.frame(do.call(rbind, data_list)),
      stop(sprintf("combine_files: unsupported mode '%s'", mode))
    )
    message(sprintf(
      "  Combined dimensions: %d rows x %d columns",
      nrow(combined_df),
      ncol(combined_df)
    ))

    # ---- save ----------------------------------------------------------------
    write_fst(combined_df, path = paste0(save_path, save_name, ".fst"))
    message(sprintf("  Saved: %s", paste0(save_path, save_name, ".fst")))

    write_feather(combined_df, paste0(save_path, save_name, ".feather"))
    message(sprintf(
      "  Saved as feather: %s",
      paste0(save_path, save_name, ".feather")
    ))

    if (write_csv) {
      out_csv <- paste0(save_path, save_name, ".csv")
      write.table(
        combined_df,
        out_csv,
        append = FALSE,
        sep = ",",
        dec = ".",
        row.names = FALSE,
        col.names = TRUE
      )
      message(sprintf("  Saved: %s", out_csv))
    }
  } else if (file_type == "trait_results") {
    # ---- accumulators --------------------------------------------------------
    n <- length(input_files)
    all_trait_result_list <- vector("list", n)
    additive_scaled_list <- vector("list", n)
    epistatic_scaled_list <- vector("list", n)
    add_scores_list <- vector("list", n)
    epi_scores_list <- vector("list", n)
    breeding_vals_list <- vector("list", n)
    intercept_list <- vector("list", n)
    alphas_list <- vector("list", n)
    combined_mean_result <- list()
    combined_var_result <- list()
    all_these_traits <- list()

    # ---- load ----------------------------------------------------------------
    for (i in seq_along(input_files)) {
      message(sprintf(
        "  Loading file %d / %d: %s",
        i,
        length(input_files),
        basename(input_files[i])
      ))
      load(input_files[i]) # provides: these_traits, all_trait_result,
      # additive_scaled, epistatic_scaled,
      # add_scores, epi_scores, intercept,
      # breeding_vals, alphas, mean, var
      all_trait_result_list[[i]] <- all_trait_result
      additive_scaled_list[[i]] <- additive_scaled
      epistatic_scaled_list[[i]] <- epistatic_scaled
      add_scores_list[[i]] <- add_scores
      epi_scores_list[[i]] <- epi_scores
      breeding_vals_list[[i]] <- breeding_vals
      intercept_list[[i]] <- intercept
      alphas_list[[i]] <- alphas
      combined_mean_result <- c(combined_mean_result, mean)
      combined_var_result <- c(combined_var_result, var)
      all_these_traits <- c(all_these_traits, these_traits)
    }

    # ---- combine -------------------------------------------------------------
    combined_traits_result <- as.data.frame(do.call(
      cbind,
      all_trait_result_list
    ))
    combined_additive_result <- do.call(rbind, additive_scaled_list)
    combined_epistatic_result <- do.call(rbind, epistatic_scaled_list)
    combined_add_scores <- do.call(cbind, add_scores_list)
    combined_epi_scores <- do.call(cbind, epi_scores_list)
    combined_breeding_vals <- do.call(cbind, breeding_vals_list)

    message(sprintf(
      "  Combined dimensions: %d rows x %d columns",
      nrow(combined_traits_result),
      ncol(combined_traits_result)
    ))

    # ---- save ----------------------------------------------------------------
    write_fst(
      combined_traits_result,
      path = paste0(save_path, save_name_gv, ".fst")
    )
    message(sprintf("  Saved: %s", paste0(save_path, save_name_gv, ".fst")))

    write_feather(
      combined_traits_result,
      paste0(save_path, save_name_gv, ".feather")
    )
    message(sprintf(
      "  Saved as feather: %s",
      paste0(save_path, save_name_gv, ".feather")
    ))

    out_rda <- paste0(save_path, save_name_others, ".Rda")
    save(
      combined_additive_result,
      combined_epistatic_result,
      all_these_traits,
      combined_var_result,
      combined_mean_result,
      combined_add_scores,
      combined_epi_scores,
      combined_breeding_vals,
      file = out_rda
    )
    message(sprintf("  Saved: %s", out_rda))
  } else {
    stop(sprintf("combine_files: unknown file_type '%s'", file_type))
  }
}


center_geno_matrix <- function(uncentered_geno_matrix, ploidy) {
  # Centers a raw genotype matrix to {-1, 0, 1} encoding using (X - ploidy/2) * (2/ploidy).
  # Parameters:
  # - uncentered_geno_matrix (matrix): Raw genotype dosage matrix (individuals x SNPs)
  # - ploidy (integer): Ploidy level of the organism (e.g., 1 for haploid, 2 for diploid)
  # Returns:
  # The centered genotype matrix with the same dimensions as the input.
  message("Pre-computing centered genotype matrix...")

  centered_geno_precomputed <- (uncentered_geno_matrix - (ploidy / 2)) *
    (2 / ploidy)
  message("Centered genotype matrix ready!")

  return(centered_geno_precomputed)
}

split_and_write_data <- function(
  geno_matrix_uncentered,
  geno_matrix_centered,
  pheno_matrix,
  seed,
  suffix,
  output_dir = "."
) {
  # Splits genotype and phenotype matrices into train/test sets and writes
  # feather, CSV, and ID files to disk.
  # Parameters:
  # - geno_matrix_uncentered (matrix or data.frame): Raw genotype matrix
  #   (individuals x SNPs)
  # - geno_matrix_centered (matrix or data.frame): Centered genotype matrix
  #   (individuals x SNPs)
  # - pheno_matrix (matrix or data.frame): Phenotype matrix
  #   (individuals x traits)
  # - seed (integer): Random seed for the train/test split
  # - suffix (character): Suffix appended to phenotype output filenames
  #   (e.g., "_normalized")
  # - output_dir (character): Directory to write all output files (default: ".")
  # Outputs written to output_dir:
  # - {prefix}_seed_{seed}_train/test_genotypes_centered.feather
  # - {prefix}_seed_{seed}_train/test_genotypes_uncentered.feather
  # - {prefix}_seed_{seed}_train/test_phenotypes{suffix}.feather
  # - {prefix}_seed_{seed}_all/train/test_phenotypes{suffix}.csv
  # - {prefix}_seed_{seed}_train/test_genotypes_uncentered/centered.csv
  # - {prefix}_seed_{seed}_train/test_ids.txt
  set.seed(seed)

  # --- Validate inputs ---
  if (
    !is.data.frame(geno_matrix_uncentered) && !is.matrix(geno_matrix_uncentered)
  ) {
    stop("geno_matrix_uncentered must be a matrix or data frame.")
  }
  if (
    !is.data.frame(geno_matrix_centered) && !is.matrix(geno_matrix_centered)
  ) {
    stop("geno_matrix_centered must be a matrix or data frame.")
  }
  if (!is.data.frame(pheno_matrix) && !is.matrix(pheno_matrix)) {
    stop("pheno_matrix must be a matrix or data frame.")
  }

  geno_matrix_uncentered <- as.data.frame(geno_matrix_uncentered)
  geno_matrix_centered <- as.data.frame(geno_matrix_centered)

  geno_matrix_uncentered$IID <- ids$ID
  geno_matrix_centered$IID <- ids$ID
  pheno_matrix <- as.data.frame(pheno_matrix)
  pheno_matrix$IID <- ids$ID
  stopifnot(nrow(geno_matrix_uncentered) == nrow(pheno_matrix))

  n_samples <- nrow(geno_matrix_uncentered)
  id_col <- "IID"

  # --- Create train/test split ---
  n_train <- round(n_samples * 0.85)
  train_idx <- sort(sample(ids$ID, size = n_train, replace = FALSE))
  test_idx <- sort(setdiff(ids$ID, train_idx))

  message("Splitting data into test and train sets...")
  geno_train_uncentered <- geno_matrix_uncentered[
    match(train_idx, geno_matrix_uncentered[[id_col]]),
  ] %>%
    relocate("IID")
  geno_test_uncentered <- geno_matrix_uncentered[
    match(test_idx, geno_matrix_uncentered[[id_col]]),
  ] %>%
    relocate("IID")
  geno_train_centered <- geno_matrix_centered[
    match(train_idx, geno_matrix_centered[[id_col]]),
  ] %>%
    relocate("IID")
  geno_test_centered <- geno_matrix_centered[
    match(test_idx, geno_matrix_centered[[id_col]]),
  ] %>%
    relocate("IID")

  # Re-order phenotype rows to match genotype row order
  pheno_train <- pheno_matrix[match(train_idx, pheno_matrix[[id_col]]), ] %>%
    relocate("IID")
  pheno_test <- pheno_matrix[match(test_idx, pheno_matrix[[id_col]]), ] %>%
    relocate("IID")
  pheno_matrix <- pheno_matrix %>% relocate("IID")
  if (args$biobank) {
    pheno_matrix$FID <- ids$FID
    pheno_matrix <- pheno_matrix %>% relocate("FID")
  } else {
    pheno_matrix$IID <- paste0("i", pheno_matrix$IID)
  }

  # --- Helper: ensure output directory exists ---
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  path <- function(filename) file.path(output_dir, filename)

  message("Writing test and train IDs...")
  # --- Write ID files ---
  write_ids <- function(id_df, filename) {
    write.table(
      id_df,
      file = path(filename),
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE
    )
  }

  if (args$biobank) {
    id_train_df <- data.frame(
      `#FID` = ids$FID[ids$ID %in% train_idx],
      `IID` = ids$ID[ids$ID %in% train_idx],
      check.names = FALSE
    )
    id_test_df <- data.frame(
      `#FID` = ids$FID[ids$ID %in% test_idx],
      `IID` = ids$ID[ids$ID %in% test_idx],
      check.names = FALSE
    )
  } else {
    corrected_train_ids <- paste0("i", train_idx)
    id_train_df <- data.frame(`#IID` = corrected_train_ids, check.names = FALSE)
    corrected_test_ids <- paste0("i", test_idx)
    id_test_df <- data.frame(`#IID` = corrected_test_ids, check.names = FALSE)
  }
  write_ids(
    id_train_df,
    paste0(args$file_save_prefix, "_seed_", args$seed, "_train_ids.txt")
  )
  write_ids(
    id_test_df,
    paste0(args$file_save_prefix, "_seed_", args$seed, "_test_ids.txt")
  )

  message("Writing genotype feather files...")
  # --- Write Feather files ---
  write_feather(
    geno_train_uncentered,
    path(paste0(
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_train_genotypes_uncentered.feather"
    ))
  )
  write_feather(
    geno_test_uncentered,
    path(paste0(
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_test_genotypes_uncentered.feather"
    ))
  )
  write_feather(
    geno_train_centered,
    path(paste0(
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_train_genotypes_centered.feather"
    ))
  )
  write_feather(
    geno_test_centered,
    path(paste0(
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_test_genotypes_centered.feather"
    ))
  )
  write_feather(
    pheno_train,
    path(paste0(
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_train_phenotypes",
      suffix,
      ".feather"
    ))
  )
  write_feather(
    pheno_test,
    path(paste0(
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_test_phenotypes",
      suffix,
      ".feather"
    ))
  )

  message("Writing phenotype csv files...")
  # --- Write CSV files ---
  write.table(
    pheno_matrix,
    path(paste0(
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_all_phenotypes",
      suffix,
      ".csv"
    )),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
  write.table(
    pheno_train,
    path(paste0(
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_train_phenotypes",
      suffix,
      ".csv"
    )),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
  write.table(
    pheno_test,
    path(paste0(
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_test_phenotypes",
      suffix,
      ".csv"
    )),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )

  message("Writing genotype csv files...")
  # --- Write CSV files ---
  write.table(
    geno_train_uncentered,
    path(paste0(
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_train_genotypes_uncentered.csv"
    )),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
  write.table(
    geno_test_uncentered,
    path(paste0(
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_test_genotypes_uncentered.csv"
    )),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
  write.table(
    geno_train_centered,
    path(paste0(
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_train_genotypes_centered.csv"
    )),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
  write.table(
    geno_test_centered,
    path(paste0(
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_test_genotypes_centered.csv"
    )),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )

  # --- Summary message ---
  message(sprintf(
    "Split complete.\n  Training samples : %d (%.1f%%)\n  Testing samples  : %d (%.1f%%)\n  Output directory : %s",
    nrow(geno_train_uncentered),
    100 * nrow(geno_train_uncentered) / n_samples,
    nrow(geno_test_uncentered),
    100 * nrow(geno_test_uncentered) / n_samples,
    normalizePath(output_dir)
  ))
}
remove_relevant_files <- function(output_dir, file_pattern) {
  # Deletes all files in output_dir whose names match file_pattern, exiting
  # with an error if any file cannot be removed.
  # Parameters:
  # - output_dir (character): Directory to search for files to delete
  # - file_pattern (character): Regex pattern passed to list.files()
  files_to_remove <- list.files(
    path = gsub("/$", "", output_dir),
    pattern = file_pattern,
    full.names = TRUE # Get full file paths
  )
  file.remove(files_to_remove)
  if (any(file.exists(files_to_remove))) {
    cat("Failed to remove files.")
    quit(save = "no")
  }
}


## Read in the SNP information list from original data
if (args$biobank) {
  if (endsWith(args$snp_file, ".bim")) {
    snp_map <- fread(args$snp_file, header = FALSE, sep = "\t")
    colnames(snp_map) <- c("Chromosome", "SNP_pre", "cM", "pos", "alt", "ref")
    snp_map <- snp_map %>% mutate(SNP = paste0(SNP_pre, "_", ref))

    write_feather(
      snp_map,
      paste0(args$geno_dir, args$file_save_prefix, "_snpmap.feather")
    )
  } else {
    snp_map <- read_feather(paste0(
      args$geno_dir,
      args$file_save_prefix,
      "_snpmap.feather"
    ))
  }
  snps_per_chrom <- snp_map %>%
    group_by(Chromosome) %>%
    summarize(n_snp = n(), min_snp = min(pos), max_snp = max(pos))
} else {
  snplist_full_pre <- read.csv(paste0(args$snp_file), sep = "\t")
  snp_map <- snplist_full_pre[!grepl("END", snplist_full_pre$SNP), ] %>%
    dplyr::mutate(SNP_1index = SNP_0index + 1)

  snps_per_chrom <- snp_map %>%
    group_by(Chromosome) %>%
    summarize(
      n_snp = n(),
      min_snp0 = min(SNP_0index),
      max_snp0 = max(SNP_0index),
      min_snp1 = min(SNP_1index),
      max_snp1 = max(SNP_1index)
    )
}
if (is.na(args$chromosome_numbers)) {
  args$chromosome_numbers <- nrow(snps_per_chrom)
}


## Load in nontransposed binarized yeast genotype data and (maybe) save it to R data
if (endsWith(args$geno, ".feather") | endsWith(args$geno, ".raw")) {
  message("Loading uncentered genotype matrix...")
  if (endsWith(args$geno, ".feather")) {
    genotypes_binarized <- read_feather(paste0(args$geno_dir, args$geno))
  } else if (endsWith(args$geno, ".raw") & args$biobank) {
    raw_genotypes <- fread(
      paste0(args$geno_dir, args$geno),
      header = TRUE,
      sep = "\t"
    )
    raw_genotypes_nona <- raw_genotypes %>%
      mutate(across(everything(), ~ replace_na(.x, 0)))
    genotypes_binarized <- raw_genotypes_nona %>%
      select(-PAT, -MAT, -SEX, -PHENOTYPE) %>%
      rename(ID = IID)
  }

  if (!is.na(args$subset_geno)) {
    if (args$subset_geno < nrow(genotypes_binarized)) {
      message(paste("Subsetting to", args$subset_geno, "individuals..."))
      set.seed(args$geno_seed)
      sample_ids <- sort(sample(
        1:nrow(genotypes_binarized),
        size = args$subset_geno,
        replace = FALSE
      ))
      genotypes_binarized <- genotypes_binarized[sample_ids, ]
      suffix <- paste0("_subset_", args$subset_geno)
    } else {
      message(paste(
        "Data contains",
        nrow(genotypes_binarized),
        "individuals. Not subsetting."
      ))
      suffix <- ""
    }
  } else {
    suffix <- ""
  }

  if (!is.na(args$subset_snps)) {
    id_cols <- intersect(c("FID", "ID"), colnames(genotypes_binarized))
    snp_cols <- setdiff(colnames(genotypes_binarized), id_cols)
    if (args$subset_snps < length(snp_cols)) {
      message(paste("Subsetting to", args$subset_snps, "SNPs..."))
      set.seed(args$geno_seed)
      sampled_snps <- sample(snp_cols, size = args$subset_snps, replace = FALSE)
      snp_map <- snp_map[snp_map$SNP %in% sampled_snps, ]
      # Order columns to match snp_map so downstream assertions on column order hold
      genotypes_binarized <- genotypes_binarized[, c(id_cols, snp_map$SNP)]
      snps_per_chrom <- snp_map %>%
        group_by(Chromosome) %>%
        summarize(n_snp = n())
      suffix <- paste0(suffix, "_snpsubset_", args$subset_snps)
    } else {
      message(paste(
        "Data contains",
        length(snp_cols),
        "SNPs. Not subsetting."
      ))
    }
  }

  if (args$biobank) {
    ids <- genotypes_binarized[, c("FID", "ID")]
    geno_matrix_uncentered_pre <- genotypes_binarized %>% select(-FID, -ID)
  } else {
    ids <- genotypes_binarized[, "ID"]
    geno_matrix_uncentered_pre <- genotypes_binarized %>% select(-ID)
  }
  geno_matrix_uncentered <- as.matrix(geno_matrix_uncentered_pre)

  geno_matrix_centered_pre <- center_geno_matrix(
    geno_matrix_uncentered,
    ploidy = args$ploidy
  )
  geno_matrix_centered_pre <- as.data.frame(geno_matrix_centered_pre)
  geno_matrix_centered_pre$ID <- ids$ID
  geno_matrix_centered_pre <- geno_matrix_centered_pre %>% relocate(ID)

  geno_matrix_uncentered_pre$ID <- ids$ID
  geno_matrix_uncentered_pre <- geno_matrix_uncentered_pre %>% relocate(ID)

  output_filepath_uncentered <- paste0(
    args$geno_dir,
    file_prefix,
    suffix,
    "_uncentered.feather"
  )
  output_filepath_centered <- paste0(
    args$geno_dir,
    file_prefix,
    suffix,
    "_centered.feather"
  )
  ids_filepath <- paste0(args$geno_dir, file_prefix, suffix, "_ids.feather")

  message("Saving centered and uncentered genotypes as feather...")
  write_feather(geno_matrix_centered_pre, output_filepath_centered)
  write_feather(geno_matrix_uncentered_pre, output_filepath_uncentered)
  write_feather(ids, ids_filepath)
} else if (
  file.exists(paste0(args$geno_dir, args$geno, "_uncentered.feather")) &
    file.exists(paste0(args$geno_dir, args$geno, "_centered.feather")) &
    file.exists(paste0(args$geno_dir, args$geno, "_ids.feather"))
) {
  message("Loading genotype matrix...")
  geno_matrix_centered_pre <- read_feather(paste0(
    args$geno_dir,
    args$geno,
    "_centered.feather"
  ))
  geno_matrix_uncentered_pre <- read_feather(paste0(
    args$geno_dir,
    args$geno,
    "_uncentered.feather"
  ))
  ids <- read_feather(paste0(
    args$geno_dir,
    args$geno,
    "_ids.feather"
  ))
} else {
  message(
    "Please provide a .feather genotype file to be centered or the
    file prefix for centered genotypes, uncentered genotypes, and ids."
  )
  quit(save = "no")
}
geno_matrix_centered <- as.matrix(geno_matrix_centered_pre[,
  2:ncol(geno_matrix_centered_pre)
])
geno_matrix_uncentered <- as.matrix(geno_matrix_uncentered_pre[,
  2:ncol(geno_matrix_uncentered_pre)
])

cat("Genotypes ready!\n\n")
message(paste(
  "Genotype data dimensions:",
  nrow(geno_matrix_centered_pre),
  "individuals x ",
  ncol(geno_matrix_centered_pre) - 1,
  "variants"
))

if (is.na(args$allele_freqs)) {
  if (!endsWith(args$geno, ".feather")) {
    suffix <- ""
  }
  message("Calculating allele frequencies...")
  allele_freqs <- apply(
    geno_matrix_centered,
    2,
    FUN = get_allele_freq,
    ploidy = args$ploidy
  )
  message("Writing allele frequencies...")
  save(
    x = allele_freqs,
    file = paste0(
      args$geno_dir,
      args$file_save_prefix,
      suffix,
      "_allele_frequencies.Rda"
    )
  )
} else if (
  !is.na(args$allele_freqs) &
    file.exists(paste0(args$geno_dir, args$allele_freqs))
) {
  message("Loading allele frequencies...")
  if (args$biobank & endsWith(args$allele_freqs, ".afreq")) {
    allele_freqs_raw <- fread(
      paste0(args$geno_dir, args$allele_freqs),
      header = TRUE,
      sep = "\t"
    )
    allele_freq_raw_df <- as.data.frame(allele_freqs_raw) %>%
      unite(col = "new_ID", ID, REF, sep = "_", remove = FALSE)

    allele_freqs <- as.numeric(unlist(allele_freqs_raw$ALT_FREQS))
    names(allele_freqs) <- allele_freq_raw_df$new_ID
  } else {
    load(paste0(args$geno_dir, args$allele_freqs))
  }
} else {
  message(
    "Please specify an allele frequency file or omit --allele_freqs
    to recalculate allele frequencies."
  )
  quit(save = "no")
}
allele_freq_matrix <- matrix(
  ncol = ncol(geno_matrix_centered),
  nrow = nrow(geno_matrix_centered),
  data = allele_freqs,
  byrow = TRUE
)
colnames(allele_freq_matrix) <- names(allele_freqs)
cat("Allele frequencies ready!\n\n")

## Initialize "simulation" and assign true traits to the individuals
if (args$single_snp_trait & args$remake_single_snp_trait) {
  remove_relevant_files(
    args$output_dir,
    paste0(
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_phenotypes_singlesnp_H2_"
    )
  )
  remove_relevant_files(
    args$output_dir,
    paste0(
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_phenotypes_singlesnp.csv"
    )
  )
  remove_relevant_files(
    args$output_dir,
    paste0(
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_phenotypes_singlesnp.f"
    )
  )
  set.seed(args$single_snp_seed)

  name <- paste(
    "trait",
    "singlesnp",
    "vavg",
    1,
    "numQTL",
    1,
    "numChr",
    1,
    sep = "_"
  )
  chosen <- choose_snps(
    snp_map,
    loci_per_chr = 1,
    relAA = 0,
    trait_name = name,
    choose_single_snp = TRUE
  )
  additive_all <- chosen$additive
  snp_effect_all <- chosen$snp_effect

  save(
    additive_all,
    snp_effect_all,
    file = paste0(
      args$output_dir,
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_chosensnps_singlesnp_unscaled.Rda"
    )
  )

  message("Calculating Single Snp Traits...")

  trait_result <- calculate_traits(
    snp_map,
    geno_matrix_uncentered,
    additive = additive_all,
    allele_freq_matrix = allele_freq_matrix,
    traits = name,
    ploidy = args$ploidy,
    rescale_scores = TRUE,
    do_epistatic = FALSE,
    centered_geno_precomputed = geno_matrix_centered
  )

  message("Finished Calculating Single Snp.")

  all_trait_result <- trait_result$genetic_values
  additive_scaled <- trait_result$additive
  additive_scaled$intercept <- 0
  save(
    additive_scaled,
    name,
    file = paste0(
      args$output_dir,
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_chosensnps_singlesnp_scaled.Rda"
    )
  )
  combined_traits_result <- as.data.frame(all_trait_result)
  write_feather(
    combined_traits_result,
    paste0(
      args$output_dir,
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_geneticvalues_singlesnp.feather"
    )
  )

  varG <- var(all_trait_result)
  for (H2 in args$broad_sense) {
    this_pheno <- add_environment(as.data.frame(all_trait_result), H2, varG)
    write_fst(
      this_pheno,
      path = paste0(
        args$output_dir,
        args$file_save_prefix,
        "_seed_",
        args$seed,
        "_phenotypes_singlesnp_H2_",
        H2,
        ".fst"
      )
    )
    message(paste0("Done with H2=", H2, "."))
  }

  combine_files(
    args$output_dir,
    paste0(
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_phenotypes_singlesnp_H2_.*\\.fst"
    ),
    mode = "cbind",
    "fst",
    args$output_dir,
    save_name = paste0(
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_phenotypes_singlesnp"
    ),
    write_csv = TRUE
  )
  message("Done with Single Snp!")
  cat("\n")
} else if (
  args$single_snp_trait &
    file.exists(paste0(
      args$output_dir,
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_phenotypes_singlesnp.feather"
    ))
) {
  combined_pheno_result_singlesnp <- read_feather(
    paste0(
      args$output_dir,
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_phenotypes_singlesnp.feather"
    )
  )
} else if (args$single_snp_trait) {
  message(
    "Please make the single-snp trait by setting --remake_single_snp_trait
    to TRUE or set --single_snp_trait to FALSE."
  )
  quit(save = "no")
}

## Initialize "simulation" and assign true traits to the individuals (1.5m)
if (args$remake_traits) {
  message("Calculating Multi-Snp Traits...")
  remove_relevant_files(
    args$output_dir,
    paste0(
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_chosensnps_unscaled.Rda"
    )
  )

  trait_names <- c()
  vars_all <- c()
  means_all <- c()
  additive_all <- matrix(ncol = 9, nrow = 0)
  epistatic_all <- matrix(ncol = 19, nrow = 0)

  set.seed(args$seed)
  for (additive_ratio in args$vavg_ratios) {
    relative_additive <- (1 - additive_ratio) / additive_ratio
    message(paste0("Working on Va/Vg=", additive_ratio, "..."))
    for (QTL_number in args$QTL_numbers) {
      message(paste0("Working on # QTLs=", QTL_number, "..."))

      for (num_chromosomes_affected in args$chromosome_numbers) {
        chromosomes <- sample(snps_per_chrom$Chromosome[
          snps_per_chrom$n_snp > (QTL_number / num_chromosomes_affected)
        ])
        if (length(chromosomes) >= num_chromosomes_affected) {
          loci_per_chr <- distribute_QTLs_evenly(
            snps_per_chrom,
            chromosomes[1:num_chromosomes_affected],
            QTL_number
          )
          stopifnot(length(loci_per_chr) == nrow(snps_per_chrom))
          name <- paste(
            "trait",
            "overlap",
            "vavg",
            additive_ratio,
            "numQTL",
            QTL_number,
            "numChr",
            num_chromosomes_affected,
            sep = "_"
          )

          chosen <- choose_snps(
            snp_map,
            loci_per_chr,
            relative_additive,
            name,
            ploidy = args$ploidy,
            epistatic_overlap = TRUE
          )
          additive_all <- rbind(additive_all, chosen$additive)
          epistatic_all <- rbind(epistatic_all, chosen$epistatic)
          chosen_var <- chosen$var
          chosen_mean <- chosen$mean
          trait_names <- c(trait_names, name)
          vars_all <- c(vars_all, chosen_var)
          means_all <- c(means_all, chosen_mean)
        } else if (num_chromosomes_affected == nrow(snps_per_chrom)) {
          loci_per_chr <- distribute_QTLs_proportionally(
            snps_per_chrom$n_snp,
            QTL_number
          )
          stopifnot(length(loci_per_chr) == nrow(snps_per_chrom))

          name <- paste(
            "trait",
            "overlap",
            "vavg",
            additive_ratio,
            "numQTL",
            QTL_number,
            "numChr",
            num_chromosomes_affected,
            sep = "_"
          )
          chosen <- choose_snps(
            snp_map,
            loci_per_chr,
            relative_additive,
            name,
            ploidy = args$ploidy,
            epistatic_overlap = TRUE
          )
          additive_all <- rbind(additive_all, chosen$additive)
          epistatic_all <- rbind(epistatic_all, chosen$epistatic)
          chosen_var <- chosen$var
          chosen_mean <- chosen$mean
          trait_names <- c(trait_names, name)
          vars_all <- c(vars_all, chosen_var)
          means_all <- c(means_all, chosen_mean)
        } else {
          message(paste(
            "Too many SNPs per chromosome, skipping",
            QTL_number,
            "QTLs for",
            num_chromosomes_affected,
            "chromosomes."
          ))
        }
      }
    }
  }
  message(paste(length(trait_names), "traits created."))
  message("Saving Traits...")
  save(
    additive_all,
    epistatic_all,
    trait_names,
    vars_all,
    means_all,
    file = paste0(
      args$output_dir,
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_chosensnps_unscaled.Rda"
    )
  )
} else if (
  !args$remake_traits &
    file.exists(paste0(
      args$output_dir,
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_chosensnps_unscaled.Rda"
    ))
) {
  load(paste0(
    args$output_dir,
    args$file_save_prefix,
    "_seed_",
    args$seed,
    "_chosensnps_unscaled.Rda"
  ))
  message(paste(length(trait_names), "traits loaded."))
} else {
  message("Please make the traits by setting --remake_traits to TRUE.")
  quit(save = "no")
}
cat("Traits ready!\n\n")

# ==============================================================================
# PHENOTYPE CALCULATION
# ==============================================================================

if (args$calculate_phenos) {
  remove_relevant_files(
    args$output_dir,
    paste0(
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_variancecomponents_trait_"
    )
  )
  remove_relevant_files(
    args$output_dir,
    paste0(
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_intermediates_trait_"
    )
  )
  remove_relevant_files(
    args$output_dir,
    paste0(
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_chosensnps_scaled.Rda"
    )
  )
  remove_relevant_files(
    args$output_dir,
    paste0(args$file_save_prefix, "_seed_", args$seed, "_geneticvalues.f")
  )

  message("Calculating Phenotypes...")
  num_traits <- length(trait_names)
  for (i in 1:num_traits) {
    message(paste("Starting trait:", i))
    these_traits <- trait_names[i]
    trait_result <- calculate_traits(
      snp_map,
      geno_matrix_uncentered,
      additive = additive_all,
      epistatic = epistatic_all,
      traits = these_traits,
      ploidy = args$ploidy,
      rescale_scores = TRUE,
      allele_freq_matrix = allele_freq_matrix,
      save_name = paste0(
        args$output_dir,
        args$file_save_prefix,
        "_seed_",
        args$seed,
        "_variancecomponents_trait_",
        i,
        "_of_",
        num_traits
      ),
      mean = means_all[i],
      var = vars_all[i],
      do_epistatic = TRUE,
      centered_geno_precomputed = geno_matrix_centered
    )
    all_trait_result <- trait_result$genetic_values
    additive_scaled <- trait_result$additive
    epistatic_scaled <- trait_result$epistatic
    add_scores <- trait_result$add_scores
    epi_scores <- trait_result$epi_scores
    intercept <- trait_result$intercept
    breeding_vals <- trait_result$breeding_vals
    alphas <- trait_result$alphas
    mean <- trait_result$mean
    var <- trait_result$var

    message(paste("Saving trait:", i))
    save(
      these_traits,
      all_trait_result,
      additive_scaled,
      epistatic_scaled,
      add_scores,
      epi_scores,
      intercept,
      breeding_vals,
      alphas,
      mean,
      var,
      file = paste0(
        args$output_dir,
        args$file_save_prefix,
        "_seed_",
        args$seed,
        "_intermediates_trait_",
        i,
        "_of_",
        num_traits,
        ".Rda"
      )
    )
  }

  cat("\n")
  message("Combining Files...")

  combine_files(
    args$output_dir,
    paste0(
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_intermediates_trait_.*\\_of_",
      num_traits,
      ".Rda"
    ),
    mode = NA,
    "trait_results",
    args$output_dir,
    save_name_gv = paste0(
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_geneticvalues"
    ),
    save_name_others = paste0(
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_chosensnps_scaled"
    )
  )

  for (vc_suffix in c(
    "additive_scaled",
    "epistatic_scaled",
    "additive_unscaled",
    "epistatic_unscaled"
  )) {
    combine_files(
      args$output_dir,
      paste0(
        args$file_save_prefix,
        "_seed_",
        args$seed,
        "_variancecomponents_trait_.*\\_of_",
        num_traits,
        "_",
        vc_suffix,
        ".fst"
      ),
      mode = "cbind",
      "fst",
      args$output_dir,
      save_name = paste0(
        args$file_save_prefix,
        "_seed_",
        args$seed,
        "_variancecomponents_",
        vc_suffix
      )
    )
  }

  # ----------------------------------------------------------------------------
  # Wide additive effects matrix
  # Rows = all SNPs (snp1 ... snp41594), columns = one per (trait x H2 level).
  # The scaled additive effects are the same regardless of environmental variance,
  # so each trait column is duplicated once per broad_sense value and named
  # <trait_name>_H2_<h2>, matching the phenotype column naming convention.
  # ----------------------------------------------------------------------------
  message("Building wide additive effects matrix...")
  load(paste0(
    args$output_dir,
    args$file_save_prefix,
    "_seed_",
    args$seed,
    "_chosensnps_scaled.Rda"
  ))
  # combined_additive_result has columns: SNP, name (trait), addEff, intercept

  additive_wide <- pivot_wider(
    combined_additive_result[, c("SNP", "name", "addEff")],
    id_cols = "SNP",
    names_from = "name",
    values_from = "addEff",
    values_fill = 0
  )
  # Ensure every SNP in the full map is present (zeroes for non-QTL SNPs)
  additive_wide <- left_join(
    data.frame(SNP = snp_map$SNP),
    additive_wide,
    by = "SNP"
  ) %>%
    mutate(across(where(is.numeric), ~ replace_na(.x, 0)))

  # Duplicate each trait column for each heritability level
  base_trait_cols <- setdiff(colnames(additive_wide), "SNP")
  h2_copies <- lapply(args$broad_sense, function(h2) {
    df <- additive_wide[, base_trait_cols, drop = FALSE]
    colnames(df) <- paste0(base_trait_cols, "_H2_", h2)
    df
  })
  additive_wide_full <- cbind(additive_wide["SNP"], do.call(cbind, h2_copies))

  out_add_wide <- paste0(
    args$output_dir,
    args$file_save_prefix,
    "_seed_",
    args$seed,
    "_additive_effects_wide.feather"
  )
  write_feather(additive_wide_full, out_add_wide)
  message(sprintf(
    "  Saved wide additive effects (%d SNPs x %d trait-H2 columns): %s",
    nrow(additive_wide_full),
    ncol(additive_wide_full) - 1L,
    out_add_wide
  ))

  message("Adding Environment...")
  combined_traits_result <- read_feather(
    paste0(
      args$output_dir,
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_geneticvalues.feather"
    )
  )
  varG <- sapply(combined_traits_result, var)

  remove_relevant_files(
    args$output_dir,
    paste0(args$file_save_prefix, "_seed_", args$seed, "_phenotypes_H2_")
  )
  remove_relevant_files(
    args$output_dir,
    paste0(args$file_save_prefix, "_seed_", args$seed, "_phenotypes.csv")
  )
  remove_relevant_files(
    args$output_dir,
    paste0(args$file_save_prefix, "_seed_", args$seed, "_phenotypes.f")
  )

  for (H2 in args$broad_sense) {
    this_pheno <- add_environment(combined_traits_result, H2, varG)
    write_fst(
      x = this_pheno,
      path = paste0(
        args$output_dir,
        args$file_save_prefix,
        "_seed_",
        args$seed,
        "_phenotypes_H2_",
        H2,
        ".fst"
      )
    )
    message(paste0("Done with H2=", H2, "."))
  }

  message("Saving Files...")
  combine_files(
    args$output_dir,
    paste0(
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_phenotypes_H2_.*\\.fst"
    ),
    mode = "cbind",
    "fst",
    args$output_dir,
    save_name = paste0(
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_phenotypes"
    ),
    write_csv = TRUE
  )

  combined_pheno_result <- read_feather(
    paste0(
      args$output_dir,
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_phenotypes.feather"
    )
  )
} else if (
  !args$calculate_phenos &
    file.exists(paste0(
      args$output_dir,
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_phenotypes.feather"
    ))
) {
  combined_pheno_result <- read_feather(
    paste0(
      args$output_dir,
      args$file_save_prefix,
      "_seed_",
      args$seed,
      "_phenotypes.feather"
    )
  )
} else {
  message(
    "Please calculate the phenotype values by setting --calculate_phenos to TRUE."
  )
  quit(save = "no")
}
cat("Phenotypes ready!\n\n")

if (args$save_test_train) {
  if (exists("combined_pheno_result_singlesnp")) {
    all_phenos <- cbind(combined_pheno_result, combined_pheno_result_singlesnp)
  } else {
    all_phenos <- combined_pheno_result
  }
  if (args$normalize) {
    message("Normalizing phenotypes...")
    all_phenos[, -1] <- scale(all_phenos[, -1])
    norm_suffix <- "_normalized"
  } else {
    norm_suffix <- ""
  }

  split_and_write_data(
    geno_matrix_uncentered = geno_matrix_uncentered,
    geno_matrix_centered = geno_matrix_centered,
    pheno_matrix = all_phenos,
    seed = args$seed,
    suffix = norm_suffix,
    output_dir = args$test_train_output_dir
  )
}

cat("\n")
message("Done!")
