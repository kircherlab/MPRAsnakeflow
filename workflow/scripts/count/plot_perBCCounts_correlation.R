# This script generates correlation plots and statistics for DNA and RNA counts per barcode across multiple replicates.
# It uses the following libraries: optparse, cowplot, and tidyverse.
#
# Usage:
# Rscript plot_perBCCounts_correlation.R -c <condition> -f <files> -r
#     <replicates> --mindnacounts <mindnacounts> --minrnacounts <minrnacounts> -o <outdir>
#
# Arguments:
# -c, --condition: Condition name (character)
# -f, --files: Comma-separated input files of assigned counts (character)
# -r, --replicates: Comma-separated names of the replicates (same order as files) (character)
# --mindnacounts: Minimum DNA counts required (integer)
# --minrnacounts: Minimum RNA counts required (integer)
# -o, --outdir: Output directory for the plots and table (character)
#
# The script performs the following steps:
# 1. Parses command-line arguments.
# 2. Checks for required arguments and stops execution if any are missing.
# 3. Reads input files and replicates.
# 4. For each pair of replicates, reads the data, filters based on
#    minimum counts, normalizes the counts, and calculates log2 values.
# 5. Generates pairwise correlation plots for DNA, RNA, and RNA/DNA ratio counts.
# 6. Calculates correlation statistics (Pearson and Spearman) for the normalized and log2-transformed counts.
# 7. Writes the correlation plots and statistics to the specified output directory.
#
# Functions:
# - plot_correlations_dna: Generates a scatter plot for DNA counts correlation between two replicates.
# - plot_correlations_rna: Generates a scatter plot for RNA counts correlation between two replicates.
# - plot_correlations_ratio: Generates a scatter plot for RNA/DNA ratio counts correlation between two replicates.
# - correlate: Calculates the correlation coefficient between two vectors.
# - get_correlation_stats: Computes correlation statistics and other metrics for a pair of replicates.
# - write_correlation_plots: Saves the correlation plots to a file.
# - write_correlation: Writes the correlation statistics to a file.
# - read_data: Reads and processes the input data file, filtering and normalizing the counts.
library(optparse)
library(cowplot)
library(tidyverse)


option_list <- list(
  make_option(c("-c", "--condition"),
    type = "character",
    help = "Condition name",
    metavar = "character"
  ),
  make_option(c("-f", "--files"),
    type = "character",
    help = "Comma separated input files of assigned counts",
    metavar = "character"
  ),
  make_option(c("-r", "--replicates"),
    type = "character",
    help = "Comma separated name of the replicates (same order than files)",
    metavar = "character"
  ),
  make_option(c("--mindnacounts"),
    type = "integer",
    help = "Minimum DNA counts required",
    metavar = "integer"
  ),
  make_option(c("--minrnacounts"),
    type = "integer",
    help = "Minimum RNA counts required",
    metavar = "integer"
  ),
  make_option(c("-o", "--outdir"),
    type = "character",
    help = "Output directory for the plots and table",
    metavar = "character"
  )
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

required_args <- c("condition", "files", "replicates", "mindnacounts", "minrnacounts")
missing_args <- setdiff(required_args, names(opt))

if (length(missing_args) > 0) {
  stop(sprintf("Missing required arguments: %s. See script usage (--help)", paste(missing_args, collapse = ", ")))
}

outdir <- ifelse("outdir" %in% names(opt), opt$outdir, "./unknown")

# global scaling factor
scaling <- 10**6
plot_sampling <- 10**5

# condition
cond <- opt$condition

# replicates and count files
files <- strsplit(opt$files, ",")[[1]]
replicates <- strsplit(opt$replicates, ",")[[1]]
if (length(files) != length(replicates)) {
  stop("Number of input files must be euqal to number of replicates")
}

num_replicates <- length(replicates)

data <- data.frame(File = files, Replicate = replicates)
data$Condition <- cond

# pairwise comparison only if more than one replicate

#' Plot Correlations of DNA Counts
#'
#' This function generates a scatter plot to visualize the correlation between
#' log2 normalized DNA counts per barcode for two replicates. It includes Pearson
#' and Spearman correlation coefficients on the plot.
#'
#' @param data A data frame containing the DNA normalized log2 counts for two replicates.
#' @param plot_data A data frame to be used for plotting.
#' @param condition A character string specifying the condition (not used in the function).
#' @param r1 A character string specifying the first replicate.
#' @param r2 A character string specifying the second replicate.
#' @param name A character string specifying the name (not used in the function).
#'
#' @return A ggplot object representing the scatter plot with correlation coefficients.
#'
#' @import ggplot2
#' @export
#'
#' @examples
#' # Example usage:
#' data <- data.frame(DNA_normalized_log2.x = rnorm(100), DNA_normalized_log2.y = rnorm(100))
#' plot_data <- data
#' plot_correlations_dna(data, plot_data, "condition", "Rep1", "Rep2", "name")
plot_correlations_dna <- function(data, plot_data, condition, r1, r2, name) {
  max <- max(data$`DNA_normalized_log2.y`)
  min <- min(data$`DNA_normalized_log2.x`)

  dna_p <- ggplot(plot_data, aes(DNA_normalized_log2.x, DNA_normalized_log2.y)) +
    geom_point() +
    ggtitle("Log2 normalized DNA counts per barcode", subtitle = sprintf("Replicate %s and %s", r1, r2)) +
    xlab(sprintf("Replicate %s", r1)) +
    ylab(sprintf("Replicate %s", r2)) +
    geom_text(
      x = min + 0.5, y = max - 0.5,
      label = sprintf("   r = %.2f", cor(data$DNA_normalized_log2.x, data$DNA_normalized_log2.y, method = "pearson")),
      size = 10
    ) +
    geom_text(
      x = min + 0.5, y = max - 1.0,
      label = sprintf("rho = %.2f", cor(data$DNA_normalized_log2.x, data$DNA_normalized_log2.y, method = "spearman")),
      size = 10
    ) +
    geom_abline(intercept = 0, slope = 1) +
    theme_classic(base_size = 30)

  return(dna_p)
}

#' Plot Correlations of RNA Counts
#'
#' This function generates a scatter plot to visualize the correlation between
#' RNA counts from two replicates. It includes Pearson and Spearman correlation
#' coefficients on the plot.
#'
#' @param data A data frame containing RNA counts with columns `RNA_normalized_log2.x`
#' and `RNA_normalized_log2.y` for log2 normalized counts, and `RNA_normalized.x`
#' and `RNA_normalized.y` for raw counts.
#' @param plot_data A data frame to be used for plotting, typically the same as `data`.
#' @param condition A character string specifying the condition or experiment name.
#' @param r1 A character string specifying the name of the first replicate.
#' @param r2 A character string specifying the name of the second replicate.
#' @param name A character string specifying the name for the plot.
#'
#' @return A ggplot object representing the scatter plot with correlation coefficients.
#'
#' @import ggplot2
#' @export
plot_correlations_rna <- function(data, plot_data, condition, r1, r2, name) {
  max <- max(data$`RNA_normalized_log2.y`)
  min <- min(data$`RNA_normalized_log2.x`)
  rna_p <- ggplot(plot_data, aes(RNA_normalized_log2.x, RNA_normalized_log2.y)) +
    geom_point() +
    ggtitle("Log2 normalized RNA counts per barcode", subtitle = sprintf("Replicate %s and %s", r1, r2)) +
    xlab(sprintf("Replicate %s", r1)) +
    ylab(sprintf("Replicate %s", r2)) +
    geom_text(
      x = min + 0.5, y = max - 0.5,
      label = sprintf("   r = %.2f", cor(data$RNA_normalized_log2.x, data$RNA_normalized_log2.y, method = "pearson")),
      size = 10
    ) +
    geom_text(
      x = min + 0.5, y = max - 1.0,
      label = sprintf("rho = %.2f", cor(data$RNA_normalized.x, data$RNA_normalized.y, method = "spearman")), size = 10
    ) +
    geom_abline(intercept = 0, slope = 1) +
    theme_classic(base_size = 30)
  return(rna_p)
}

#' Plot Correlations Ratio
#'
#' This function creates a scatter plot to visualize the correlation between
#' log2 RNA/DNA count ratios per barcode for two replicates. It calculates and
#' displays both Pearson and Spearman correlation coefficients on the plot.
#'
#' @param data A data frame containing the columns `Ratio_log2.x` and `Ratio_log2.y`
#'   which represent the log2 RNA/DNA count ratios for two replicates.
#' @param plot_data A data frame used for plotting, containing the same columns
#'   as `data`.
#' @param condition A character string representing the condition being analyzed.
#' @param r1 A character string representing the first replicate.
#' @param r2 A character string representing the second replicate.
#' @param name A character string representing the name of the plot.
#'
#' @return A ggplot object representing the scatter plot with correlation coefficients.
#'
#' @import ggplot2
#' @export
plot_correlations_ratio <- function(data, plot_data, condition, r1, r2, name) {
  max <- max(data$`Ratio_log2.y`)
  min <- min(data$`Ratio_log2.x`)
  ratio_p <- ggplot(plot_data, aes(Ratio_log2.x, Ratio_log2.y)) +
    geom_point() +
    ggtitle("Log2 RNA/DNA count ratio per barcode", subtitle = sprintf("Replicate %s and %s", r1, r2)) +
    xlab(sprintf("Replicate %s", r1)) +
    ylab(sprintf("Replicate %s", r2)) +
    geom_text(
      x = min + 0.5, y = max - 0.5,
      label = sprintf("   r = %.2f", cor(data$Ratio_log2.x, res$Ratio_log2.y, method = "pearson")), size = 10
    ) +
    geom_text(
      x = min + 0.5, y = max - 1.0,
      label = sprintf("rho = %.2f", cor(data$Ratio.x, data$Ratio.y, method = "spearman")),
      size = 10
    ) +
    geom_abline(intercept = 0, slope = 1) +
    theme_classic(base_size = 30)
  return(ratio_p)
}

#' Correlate Two Vectors
#'
#' This function calculates the correlation between two numeric vectors using the specified method.
#'
#' @param x A numeric vector.
#' @param y A numeric vector.
#' @param method A character string indicating which correlation coefficient is to be computed.
#'        One of "pearson" (default), "kendall", or "spearman".
#' @return A character string representing the correlation coefficient rounded to five decimal places.
#' @examples
#' x <- c(1, 2, 3, 4, 5)
#' y <- c(2, 4, 6, 8, 10)
#' correlate(x, y, method = "pearson")
#' @export
correlate <- function(x, y, method) {
  return(sprintf("%.5f", cor(x, y, method = method)))
}

#' Calculate Correlation Statistics for Barcode Counts
#'
#' This function computes various correlation statistics between two replicates of barcode counts data.
#'
#' @param data A data frame containing the barcode counts and normalized values for two replicates.
#' @param n_bc_r1 An integer representing the number of barcodes in replicate A.
#' @param n_bc_r2 An integer representing the number of barcodes in replicate B.
#' @param condition A character string specifying the condition under which the data was collected.
#' @param r1 A character string specifying the name of replicate A.
#' @param r2 A character string specifying the name of replicate B.
#' @param name A character string specifying the name of the dataset.
#'
#' @return A data frame containing the correlation statistics, including:
#' \item{Condition}{The condition under which the data was collected.}
#' \item{ReplicateA}{The name of replicate A.}
#' \item{ReplicateB}{The name of replicate B.}
#' \item{number_BC_ReplicateA}{The number of barcodes in replicate A.}
#' \item{number_BC_ReplicateB}{The number of barcodes in replicate B.}
#' \item{number_BC_Joined}{The number of barcodes present in both replicates.}
#' \item{fraction_BC_ReplicateA}{The fraction of barcodes in replicate A that are present in both replicates.}
#' \item{fraction_BC_ReplicateB}{The fraction of barcodes in replicate B that are present in both replicates.}
#' \item{DNA_spearman}{The Spearman correlation coefficient for DNA normalized values between the two replicates.}
#' \item{RNA_spearman}{The Spearman correlation coefficient for RNA normalized values between the two replicates.}
#' \item{Ratio_spearman}{The Spearman correlation coefficient for ratio values between the two replicates.}
#' \item{DNA_pearson}{The Pearson correlation coefficient for DNA normalized values between the two replicates.}
#' \item{RNA_pearson}{The Pearson correlation coefficient for RNA normalized values between the two replicates.}
#' \item{Ratio_pearson}{The Pearson correlation coefficient for ratio values between the two replicates.}
#' \item{DNA_log2_pearson}{The Pearson cor coefficient for log2-transformed DNA normalized values between the two replicates.}
#' \item{RNA_log2_pearson}{The Pearson cor coefficient for log2-transformed RNA normalized values between the two replicates.}
#' \item{Ratio_log2_pearson}{The Pearson correlation coefficient for log2-transformed ratio values between the two replicates.}
#' \item{NormSymmetry}{A measure of the symmetry of the normalized ratio differences between the two replicates.}
#'
#' @examples
#' \dontrun{
#' data <- data.frame(
#'   Ratio.x = rnorm(100),
#'   Ratio.y = rnorm(100),
#'   DNA_normalized.x = rnorm(100),
#'   DNA_normalized.y = rnorm(100),
#'   RNA_normalized.x = rnorm(100),
#'   RNA_normalized.y = rnorm(100),
#'   DNA_normalized_log2.x = rnorm(100),
#'   DNA_normalized_log2.y = rnorm(100),
#'   RNA_normalized_log2.x = rnorm(100),
#'   RNA_normalized_log2.y = rnorm(100),
#'   Ratio_log2.x = rnorm(100),
#'   Ratio_log2.y = rnorm(100)
#' )
#' get_correlation_stats(data, 100, 100, "Condition1", "Rep1", "Rep2", "Dataset1")
#' }
get_correlation_stats <- function(data, n_bc_r1, n_bc_r2, condition, r1, r2, name) {
  norm <- abs(length(which((data$Ratio.x - data$Ratio.y) > 0)) - length(which((data$Ratio.x - data$Ratio.y) < 0)))
  +abs(length(which((data$Ratio.x - data$Ratio.y) > 0)) - length(which((data$Ratio.x - data$Ratio.y) < 0)))
  +abs(length(which((data$Ratio.x - data$Ratio.y) > 0)) - length(which((data$Ratio.x - data$Ratio.y) < 0)))
  outs <- data.frame(
    Condition = condition,
    ReplicateA = r1,
    ReplicateB = r2,
    number_BC_ReplicateA = n_bc_r1,
    number_BC_ReplicateB = n_bc_r2,
    number_BC_Joined = data %>% nrow(),
    fraction_BC_ReplicateA = (data %>% nrow() / n_bc_r1),
    fraction_BC_ReplicateB = (data %>% nrow() / n_bc_r2),
    DNA_spearman = correlate(data$DNA_normalized.x, data$DNA_normalized.y, "spearman"),
    RNA_spearman = correlate(data$RNA_normalized.x, data$RNA_normalized.y, "spearman"),
    Ratio_spearman = correlate(data$Ratio.x, data$Ratio.y, "spearman"),
    DNA_pearson = correlate(data$DNA_normalized.x, data$DNA_normalized.y, "pearson"),
    RNA_pearson = correlate(data$RNA_normalized.x, data$RNA_normalized.y, "pearson"),
    Ratio_pearson = correlate(data$Ratio.x, data$Ratio.y, "pearson"),
    DNA_log2_pearson = correlate(data$DNA_normalized_log2.x, data$DNA_normalized_log2.y, "pearson"),
    RNA_log2_pearson = correlate(data$RNA_normalized_log2.x, data$RNA_normalized_log2.y, "pearson"),
    Ratio_log2_pearson = correlate(data$Ratio_log2.x, data$Ratio_log2.y, "pearson"),
    NormSymmetry = norm, stringsAsFactors = FALSE
  )
  return(outs)
}


#' Write Correlation Plots
#'
#' This function takes a list of plots and a file name, combines the plots into a single plot grid,
#' and saves the resulting plot grid to a file.
#'
#' @param plots A list of ggplot objects to be combined into a single plot grid.
#' @param name A string specifying the file name to save the plot grid.
#'
#' @details
#' The function uses `cowplot::plot_grid` to combine the list of plots into a single plot grid with one column.
#' The combined plot grid is then saved to a file using `ggsave` with specified width, height, dpi, and type.
#'
#' @return None. The function is called for its side effect of saving the plot grid to a file.
#'
#' @examples
#' \dontrun{
#' plots <- list(plot1, plot2, plot3)
#' write_correlation_plots(plots, "correlation_plots.png")
#' }
write_correlation_plots <- function(plots, name) {
  correlation_plots <- cowplot::plot_grid(plotlist = plots, ncol = 1)
  # correlation_plots <- do.call("grid.arrange", c(plots))

  ggsave(name, correlation_plots, width = 15, height = 10 * length(plots), dpi = 96, type = "cairo", limitsize = FALSE)
}

#' Write Correlation Data to File
#'
#' This function writes a data frame of correlation values to a specified file.
#'
#' @param correlations A data frame containing the correlation values to be written to the file.
#' @param name A string specifying the name of the file to which the correlation values will be written.
#'
#' @return None. This function writes the data frame to a file and does not return a value.
#'
#' @examples
#' correlations <- data.frame(sample1 = c(0.1, 0.2, 0.3), sample2 = c(0.4, 0.5, 0.6))
#' write_correlation(correlations, "correlations.txt")
#'
#' @export
write_correlation <- function(correlations, name) {
  write.table(correlations, file = name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

#' Read and process data from a file
#'
#' This function reads data from a specified file, filters the data based on
#' minimum DNA and RNA counts, and normalizes the counts. It also calculates
#' the ratio of DNA to RNA and their log2 transformations.
#'
#' @param file A string specifying the path to the input file. The file should
#'   be a tab-separated text file with three columns: Barcode, DNA, and RNA.
#' @param mindnacounts An integer specifying the minimum DNA counts required
#'   for a barcode to be included in the analysis.
#' @param minrnacounts An integer specifying the minimum RNA counts required
#'   for a barcode to be included in the analysis.
#' @param scaling A numeric value used to scale the normalized counts.
#'
#' @return A data frame with the following columns:
#'   \item{Barcode}{The barcode identifier.}
#'   \item{DNA}{The original DNA counts.}
#'   \item{RNA}{The original RNA counts.}
#'   \item{DNA_normalized}{The normalized DNA counts.}
#'   \item{RNA_normalized}{The normalized RNA counts.}
#'   \item{Ratio}{The ratio of DNA_normalized to RNA_normalized.}
#'   \item{DNA_normalized_log2}{The log2 transformation of DNA_normalized.}
#'   \item{RNA_normalized_log2}{The log2 transformation of RNA_normalized.}
#'   \item{Ratio_log2}{The log2 transformation of Ratio.}
#'
#' @import dplyr
#' @importFrom magrittr %>%
#' @export
read_data <- function(file, mindnacounts, minrnacounts, scaling) {
  data <- read.table(file,
    as.is = TRUE,
    sep = "\t", header = FALSE, comment.char = "", stringsAsFactors = FALSE
  )
  colnames(data) <- c("Barcode", "DNA", "RNA")

  pseudocountdna <- if (mindnacounts == 0) 1 else 0
  pseudocountrna <- if (minrnacounts == 0) 1 else 0
  data <- data %>%
    filter(DNA >= mindnacounts, RNA >= minrnacounts) %>%
    mutate(
      DNA_normalized = (DNA + pseudocountdna) / sum(DNA + pseudocountdna) * scaling,
      RNA_normalized = (RNA + pseudocountrna) / sum(RNA + pseudocountrna) * scaling,
      Ratio = DNA_normalized / RNA_normalized,
      DNA_normalized_log2 = log2(DNA_normalized),
      RNA_normalized_log2 = log2(RNA_normalized),
      Ratio_log2 = log2(Ratio)
    )
  return(data)
}

# This script performs pairwise correlation analysis on barcode count data from multiple replicates.
# It generates correlation plots and statistics for DNA, RNA, and ratio data.
#
# Steps:
# 1. Check if the input data has more than one row.
# 2. Create pairwise combinations of replicates.
# 3. Initialize lists to store correlation plots and a data frame for correlation statistics.
# 4. Loop through each pairwise combination of replicates:
#    a. Read data for each replicate.
#    b. Perform an inner join on the barcode column to get common barcodes.
#    c. Sample the data if it exceeds a specified threshold.
#    d. Generate correlation plots for DNA, RNA, and ratio data.
#    e. Calculate and store correlation statistics.
# 5. Write the correlation statistics and plots to output files.
#
# Functions used:
# - read_data(file, mindnacounts, minrnacounts, scaling): Reads and processes the data file.
# - plot_correlations_dna(res, res_plot, cond, r1, r2, type): Generates DNA correlation plot.
# - plot_correlations_rna(res, res_plot, cond, r1, r2, type): Generates RNA correlation plot.
# - plot_correlations_ratio(res, res_plot, cond, r1, r2, type): Generates ratio correlation plot.
# - get_correlation_stats(res, n_bc_r1, n_bc_r2, cond, r1, r2, type): Calculates correlation statistics.
# - write_correlation(stats, filepath): Writes correlation statistics to a file.
# - write_correlation_plots(plots, filepath): Writes correlation plots to a file.
if (data %>% nrow() > 1) {
  # make pairwise combinations
  selected <- combn(data$Replicate, 2)
  print("sel")
  print(selected)

  plots_correlations_rna <- list()
  plots_correlations_dna <- list()
  plots_correlations_ratio <- list()
  stats_correlations <- data.frame()
  print("reps")
  for (i in seq(1, dim(selected)[2])) {
    print(selected[, i])
    r1 <- selected[1, i]
    r2 <- selected[2, i]
    data1 <- read_data(
      as.character((data %>% filter(Replicate == r1))$File),
      opt$mindnacounts, opt$minrnacounts, scaling
    )
    data2 <- read_data(
      as.character((data %>% filter(Replicate == r2))$File),
      opt$mindnacounts, opt$minrnacounts, scaling
    )

    n_bc_r1 <- data1 %>% nrow()
    n_bc_r2 <- data2 %>% nrow()

    res <- data1 %>% inner_join(data2, by = c("Barcode"))
    if (res %>% nrow() > plot_sampling) {
      res_plot <- sample_n(res, plot_sampling)
    } else {
      res_plot <- res
    }
    plots_correlations_dna[[i]] <- plot_correlations_dna(res, res_plot, cond, r1, r2, "pairwise")
    plots_correlations_rna[[i]] <- plot_correlations_rna(res, res_plot, cond, r1, r2, "pairwise")
    plots_correlations_ratio[[i]] <- plot_correlations_ratio(res, res_plot, cond, r1, r2, "pairwise")

    stats_correlations <- stats_correlations %>% bind_rows(
      get_correlation_stats(res, n_bc_r1, n_bc_r2, cond, r1, r2, "correlation")
    )
  }

  print("write correlation stats")
  write_correlation(stats_correlations, sprintf("%s_barcode_correlation.tsv", outdir))
  print("write correlation plots")
  write_correlation_plots(plots_correlations_dna, sprintf("%s_barcode_DNA_pairwise.png", outdir))
  write_correlation_plots(plots_correlations_rna, sprintf("%s_barcode_RNA_pairwise.png", outdir))
  write_correlation_plots(plots_correlations_ratio, sprintf("%s_barcode_Ratio_pairwise.png", outdir))
}
