# adapted from Vikram Agarwal by Gracie Gordon

library(ggplot2)
library(optparse)
library(cowplot)
library(dplyr)
library(tidyr)



option_list <- list(
    make_option(c("-c", "--condition"),
                type = "character",
                help = "Condition name"),
    make_option(c("-l", "--label"),
                type = "character",
                help = "Label file. (optional)"),
    make_option(c("-f", "--files"),
                type = "character",
                help = "Comma separated input files of assigned counts"),
    make_option(c("-r", "--replicates"),
                type = "character",
                help = "Comma separated name of the replicates (same order than files)"),
    make_option(
        c("-t", "--threshold"),
        type = "integer",
        default = 10,
        help = "Number of required barcodes (default 10)"
    ),
    make_option(c("-o", "--outdir"),
                type = "character",
                help = "Outdir of the plots and table.")
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
        header = T,
        stringsAsFactors = F
    ))
    colnames(label_f) <- c("name", "label")
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
        as.is = T,
        sep = "\t",
        header = T,
        stringsAsFactors = F
    ) %>%
        filter(name != "no_BC") %>%
        mutate(
            dna_normalized_log2 = log2(dna_normalized),
            rna_normalized_log2 = log2(rna_normalized),
            ratio_log2 = log2(ratio)
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
        left_join(label_f, by = c("name")) %>%
        mutate(label = replace_na(label, "NA"))
} else {
    all$label <- "NA"
}


print("Histogram plots, Boxplots, Violinplots")

plot_all_bc_per_insert <- function(data) {
    data$name <- factor(data$name)
    data$label <- as.factor(data$label)

    data <- data[order(data$log2), ]

    bymedian <-
        with(data, reorder(name, -log2, median, order = TRUE)) # nolint
    data$name <- factor(data$name, levels = levels(bymedian))
    bp <- ggplot(data, aes(x = name, y = log2, color = label)) +
        geom_boxplot() +
        xlab("insert") +
        ylab("log2 fold change") +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.title.x = element_text(size = 15),
            axis.title.y = element_text(size = 15),
            axis.text.y = element_text(size = 15),
            legend.text = element_text(size = 15)
        )
    return(bp)
}

plot_group_bc_per_insert <- function(data) {
    bp <- ggplot(data, aes(x = label, y = log2, fill = label)) +
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

hist_plot_list <- list()
box_plot_list <- list()
box_plot_thresh_list <- list()
box_plot_insert_list <- list()
box_plot_insert_thresh_list <- list()

for (n in 1:(data %>% nrow())) {
    rep <- toString(data[n, ]$Replicate)
    assigned_counts <- all %>% filter(replicate == rep)

    # Histograms
    intercept <- median(assigned_counts$n_obs_bc)
    hist_plot_list[[n]] <-
        ggplot(assigned_counts, aes(x = n_obs_bc)) +
        geom_histogram(bins = 300) +
        geom_vline(xintercept = intercept, colour = "red") +
        xlim(0, 300) +
        ggtitle(paste("replicate", rep, sep = " "))

    # Boxplots
    assigned_counts_subsample <- assigned_counts %>%
        sample_n(min(10000, assigned_counts %>% nrow()))

    box_plot_list[[n]] <-
        plot_all_bc_per_insert(assigned_counts_subsample) +
        ggtitle(paste("replicate", rep, sep = " "))

    box_plot_thresh_list[[n]] <-
        plot_all_bc_per_insert(assigned_counts_subsample %>%
                                   filter(n_obs_bc >= thresh)) +
        ggtitle(paste("replicate", rep, sep =
                          " "))

    box_plot_insert_list[[n]] <-
        plot_group_bc_per_insert(assigned_counts) +
        ggtitle(paste("replicate", rep, sep = " "))

    box_plot_insert_thresh_list[[n]] <-
        plot_group_bc_per_insert(assigned_counts %>%
                                     filter(n_obs_bc >= thresh)) +
        ggtitle(paste("replicate", rep, sep = " "))
}

hist_plot <- do.call("plot_grid", c(hist_plot_list))
ggsave(sprintf("%s_barcodesPerInsert.png", outdir), hist_plot)

box_plot <- do.call("plot_grid", c(box_plot_list))
ggsave(sprintf("%s_all_barcodesPerInsert_box.png", outdir),
       box_plot)

box_plot_thresh <- do.call("plot_grid", c(box_plot_thresh_list))
ggsave(
    sprintf("%s_all_barcodesPerInsert_box_minThreshold.png", outdir),
    box_plot_thresh
)

box_plot_insert <- do.call("plot_grid", c(box_plot_insert_list))
ggsave(sprintf("%s_group_barcodesPerInsert_box.png", outdir),
       box_plot_insert)

box_plot_insert_thresh <-
    do.call("plot_grid", c(box_plot_insert_thresh_list))
ggsave(
    sprintf("%s_group_barcodesPerInsert_box_minThreshold.png", outdir),
    box_plot_insert_thresh
)

print("Script done")
