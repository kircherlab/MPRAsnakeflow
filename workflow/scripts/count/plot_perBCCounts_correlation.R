library(optparse)
library(cowplot)
library(tidyverse)


option_list <- list(
  make_option(c("-c", "--condition"),
    type = "character",
    help = "Condition name"
  ),
  make_option(c("-f", "--files"),
    type = "character",
    help = "Comma separated input files of assigned counts"
  ),
  make_option(c("-r", "--replicates"),
    type = "character",
    help = "Comma separated name of the replicates (same order than files)"
  ),
  make_option(c("--mindnacounts"),
    type = "integer",
    help = "minimum DNA counts required"
  ),
  make_option(c("--minrnacounts"),
    type = "integer",
    help = "minimum RNA counts required"
  ),
  make_option(c("-o", "--outdir"),
    type = "character",
    help = "Outdir of the plots and table."
  )
)

parser <- OptionParser(option_list = option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

if (!"condition" %in% names(opt)) {
  stop("--condition parameter must be provided. See script usage (--help)")
}
if (!"files" %in% names(opt)) {
  stop("--files parameter must be provided. See script usage (--help)")
}
if (!"replicates" %in% names(opt)) {
  stop("--replicates parameter must be provided. See script usage (--help)")
}
if (!"minrnacounts" %in% names(opt)) {
  stop("--minrnacounts parameter must be provided. See script usage (--help)")
}
if (!"mindnacounts" %in% names(opt)) {
  stop("--mindnacounts parameter must be provided. See script usage (--help)")
}
if (!"outdir" %in% names(opt)) {
  outdir <- "./unknown"
} else {
  outdir <- opt$outdir
}

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

correlate <- function(x, y, method) {
  return(sprintf("%.5f", cor(x, y, method = method)))
}

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

write_correlation_plots <- function(plots, name) {
  correlation_plots <- cowplot::plot_grid(plotlist = plots, ncol = 1)
  # correlation_plots <- do.call("grid.arrange", c(plots))

  ggsave(name, correlation_plots, width = 15, height = 10 * length(plots), dpi = 96, type = "cairo")
}

write_correlation <- function(correlations, name) {
  write.table(correlations, file = name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

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
