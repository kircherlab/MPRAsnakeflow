# adapted from Vikram Agarwal by Gracie Gordon

library(ggplot2)
library(optparse)
library(cowplot)
library(tidyverse)



option_list <- list(
  make_option(c("-c", "--condition"),
    type = "character",
    help = "Condition name"
  ),
  make_option(c("-l", "--label"),
    type = "character",
    help = "Label file. (optional)"
  ),
  make_option(c("-f", "--files"),
    type = "character",
    help = "Comma separated input files of assigned counts"
  ),
  make_option(c("-r", "--replicates"),
    type = "character",
    help = "Comma separated name of the replicates (same order than files)"
  ),
  make_option(
    c("-t", "--threshold"),
    type = "integer",
    default = 10,
    help = "Number of required barcodes (default 10)"
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
if (!"outdir" %in% names(opt)) {
  outdir <- "./unknown"
} else {
  outdir <- opt$outdir
}

# condition
cond <- opt$condition
# labels
if ("label" %in% names(opt)) {
  label_f <- as.data.frame(read.table(
    opt$label,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
  ))
  colnames(label_f) <- c("oligo_name", "label")
  use_labels <- TRUE
} else {
  use_labels <- FALSE
}


# replicates and count files
files <- strsplit(opt$files, ",")[[1]]
replicates <- strsplit(opt$replicates, ",")[[1]]
if (length(files) != length(replicates)) {
  stop("Number of input files must be euqal to number of replicates")
}

data <- data.frame(File = files, Replicate = replicates)
data["Condition"] <- cond

print(data)

# pairwise comparison only if more than one replicate
thresh <- opt$threshold

read_data <- function(file) {
  data <- read.table(
    file,
    as.is = TRUE,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
  ) %>%
    filter(oligo_name != "no_BC") %>%
    mutate(
      dna_normalized_log2 = log2(dna_normalized),
      rna_normalized_log2 = log2(rna_normalized),
    )
  return(data)
}


print("Read data")
all <- data.frame()

for (n in 1:(data %>% nrow())) {
  print(data[n, ]$File)
  assigned_counts <- read_data(as.character(data[n, ]$File))
  assigned_counts["replicate"] <- toString(data[n, ]$Replicate)
  all <- all %>% bind_rows(assigned_counts)
}

if (use_labels) {
  all <- all %>%
    left_join(label_f, by = c("oligo_name")) %>%
    mutate(label = replace_na(label, "NA"))
} else {
  all$label <- "NA"
}

print("Histogram plots, RNA/DNA correlation plots, Violin plots")


plot_median_dna_rna_cor <- function(data) {
  data <- data %>%
    group_by(oligo_name) %>%
    summarise(dna_normalized = median(log10(dna_normalized)), rna_normalized = median(log10(rna_normalized)), n = n())
  data <- data %>% filter(n == length(replicates))
  p <- ggplot(data, aes(x = dna_normalized, y = rna_normalized)) +
    geom_point() +
    ggtitle("Median normalized counts across replicates") +
    xlab("DNA [log10]") +
    ylab("RNA [log10]") +
    geom_abline(intercept = 0, slope = 1) +
    theme_classic(base_size = 30) +
    theme(
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      legend.text = element_text(size = 15),
      plot.title = element_text(size = 20),
    ) +
    guides(fill = "none")
  return(p)
}

plot_group_bc_per_insert <- function(data) {
  bp <- ggplot(data, aes(x = label, y = log2FoldChange, fill = label)) +
    geom_violin() +
    geom_boxplot(width = 0.1, fill = "white") +
    xlab("insert") +
    ylab("log2 fold change") +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 15
      ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      legend.text = element_text(size = 15)
    ) +
    guides(fill = "none")
  return(bp)
}

ggsave(sprintf("%s_dna_vs_rna.png", outdir),
  plot_median_dna_rna_cor(all),
  width = 10, height = 10
)
ggsave(sprintf("%s_dna_vs_rna_minThreshold.png", outdir),
  plot_median_dna_rna_cor(all %>% filter(n_bc >= thresh)),
  width = 10, height = 10
)


hist_plot_list <- list()
box_plot_list <- list()
box_plot_thresh_list <- list()
box_plot_insert_list <- list()
box_plot_insert_thresh_list <- list()

for (n in 1:(data %>% nrow())) {
  rep <- toString(data[n, ]$Replicate)
  assigned_counts <- all %>% filter(replicate == rep)

  # Histograms
  x_lim_n_bc <- min(max(assigned_counts$n_bc), 300)
  intercept_median <- median(assigned_counts$n_bc)
  intercept_mean <- mean(assigned_counts$n_bc)
  hist_plot_list[[n]] <-
    ggplot(assigned_counts, aes(x = n_bc)) +
    geom_histogram(bins = x_lim_n_bc) +
    geom_vline(xintercept = intercept_median, colour = "red") +
    geom_vline(xintercept = intercept_mean, colour = "blue") +
    xlim(0, x_lim_n_bc) +
    ylab("Frequency") +
    xlab("Barcodes per oligo") +
    ggtitle(paste("replicate", rep, sep = " "))

  # Boxplots
  box_plot_insert_list[[n]] <-
    plot_group_bc_per_insert(assigned_counts) +
    ggtitle(paste("replicate", rep, sep = " "))

  box_plot_insert_thresh_list[[n]] <-
    plot_group_bc_per_insert(assigned_counts %>% filter(n_bc >= thresh)) +
    ggtitle(paste("replicate", rep, sep = " "))
}

hist_plot <- do.call("plot_grid", c(hist_plot_list))
ggsave(sprintf("%s_barcodesPerInsert.png", outdir), hist_plot, dpi = 300, type = "cairo")



box_plot_insert <- do.call("plot_grid", c(box_plot_insert_list))
ggsave(sprintf("%s_group_barcodesPerInsert_box.png", outdir),
  box_plot_insert,
  dpi = 300, type = "cairo"
)

box_plot_insert_thresh <-
  do.call("plot_grid", c(box_plot_insert_thresh_list))
ggsave(
  sprintf("%s_group_barcodesPerInsert_box_minThreshold.png", outdir),
  box_plot_insert_thresh,
  dpi = 300, type = "cairo"
)

print("Script done")
