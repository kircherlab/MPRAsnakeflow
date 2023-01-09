
library(optparse)
library(cowplot)
library(tidyverse)
library(Cairo)


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

readData <- function(file, mindnacounts, minrnacounts) {
  data <- read.table(file, as.is = T,
    sep = "\t", header = F, stringsAsFactors = F
  )
  colnames(data) <- c("Barcode", "DNA", "RNA")

  data <- data %>% filter(DNA >= mindnacounts, RNA >= minrnacounts)
  return(data)
}

print("hist")

plots_dna <- list()
plots_rna <- list()

for (n in 1:(data %>% nrow())) {
  counts <- readData(as.character(data[n, ]$File), opt$mindnacounts, opt$minrnacounts)
  intercept_median <- median(counts$DNA)
  intercept_mean <- mean(counts$DNA)
  plots_dna[[n]] <- ggplot(counts, aes(x = DNA)) +
    geom_histogram(bins=100) +
    geom_vline(xintercept = intercept_median, colour = "red") +
    geom_vline(xintercept = intercept_mean, colour = "blue") +
    xlim(0,100) +
    ggtitle(paste("replicate", data[n, ]$Replicate, sep = " "))
  intercept_median <- median(counts$RNA)
  intercept_mean <- mean(counts$RNA)
  plots_rna[[n]] <- ggplot(counts, aes(x = RNA)) +
    geom_histogram(bins=100) +
    geom_vline(xintercept = intercept_median, colour = "red") +
    geom_vline(xintercept = intercept_mean, colour = "blue") +
    xlim(0,100) +
    ggtitle(paste("replicate", data[n, ]$Replicate, sep = " ")) 
}

hist_plot <- do.call("plot_grid", c(plots_rna))
ggsave(sprintf("%s_RNA_perBarcode.png", outdir), hist_plot, dpi=300, type="cairo")

hist_plot <- do.call("plot_grid", c(plots_dna))
ggsave(sprintf("%s_DNA_perBarcode.png", outdir), hist_plot, dpi=300, type="cairo")
