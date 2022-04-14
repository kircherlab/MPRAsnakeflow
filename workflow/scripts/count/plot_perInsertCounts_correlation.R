


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

num_replicates <- length(replicates)

data <- data.frame(File = files, Replicate = replicates)
data["Condition"] <- cond

print(data)

# pairwise comparison only if more than one replicate
thresh <- opt$threshold

plot_correlations_dna <- function(data, condition, r1, r2, name) {
    dna_p <-
        ggplot(data, aes(dna_normalized_log2.x, dna_normalized_log2.y)) +
        geom_point(aes(colour = label.x), show.legend = TRUE) +
        xlim(-5, 5) +
        ylim(-5, 5) +
        xlab(sprintf(
            paste("log2 Normalized DNA count per insert,\n replicate",
                  r1)
        )) +
        ylab(sprintf(
            paste("log2 Normalized DNA count per insert,\n replicate",
                  r2)
        )) +
        geom_text(
            x = 0,
            y = 4.5,
            label = sprintf(
                "   r = %.2f",
                cor(
                    data$dna_normalized_log2.x,
                    data$dna_normalized_log2.y,
                    method = "pearson"
                )
            ),
            size = 10
        ) +
        geom_text(
            x = 0,
            y = 4,
            label = sprintf(
                "rho = %.2f",
                cor(data$dna_normalized.x, data$dna_normalized.y,
                    method = "spearman")
            ),
            size = 10
        ) +
        geom_abline(intercept = 0, slope = 1) +
        theme_classic(base_size = 30)
    return(dna_p)
}
plot_correlations_rna <- function(data, condition, r1, r2, name) {
    rna_p <-
        ggplot(data, aes(rna_normalized_log2.x, rna_normalized_log2.y)) +
        geom_point(aes(colour = label.x), show.legend = TRUE) +
        xlim(-5, 5) +
        ylim(-5, 5) +
        xlab(sprintf(
            paste("log2 Normalized RNA count per insert,\n replicate", r1)
        )) +
        ylab(sprintf(
            paste("log2 Normalized RNA count per insert,\n replicate", r2)
        )) +
        geom_text(
            x = 0,
            y = 4.5,
            label = sprintf(
                "   r = %.2f",
                cor(
                    data$rna_normalized_log2.x,
                    data$rna_normalized_log2.y,
                    method = "pearson"
                )
            ),
            size = 10
        ) +
        geom_text(
            x = 0,
            y = 4,
            label = sprintf(
                "rho = %.2f",
                cor(data$rna_normalized.x, data$rna_normalized.y,
                    method = "spearman")
            ),
            size = 10
        ) +
        geom_abline(intercept = 0, slope = 1) +
        theme_classic(base_size = 30)
    return(rna_p)
}
plot_correlations_ratio <- function(data, condition, r1, r2, name) {
    ratio_p <- ggplot(data, aes(ratio_log2.x, ratio_log2.y)) +
        geom_point(aes(colour = label.x), show.legend = TRUE) +
        xlim(-5, 5) +
        ylim(-5, 5) +
        xlab(sprintf(paste(
            "log2 RNA/DNA per insert,\n replicate", r1
        ))) +
        ylab(sprintf(paste(
            "log2 RNA/DNA per insert,\n replicate", r2
        ))) +
        geom_text(
            x = 0,
            y = 4.5,
            label = sprintf(
                "   r = %.2f",
                cor(data$ratio_log2.x, res$ratio_log2.y, method = "pearson")
            ),
            size = 10
        ) +
        geom_text(
            x = 0,
            y = 4,
            label = sprintf(
                "rho = %.2f",
                cor(data$ratio.x, data$ratio.y, method = "spearman")
            ),
            size = 10
        ) +
        geom_abline(intercept = 0, slope = 1) +
        theme_classic(base_size = 30)
    return(ratio_p)
}

correlate <- function(x, y, method) {
    return(sprintf("%.5f", cor(x, y, method = method)))
}

get_correlation_stats <-
    function(data,
             n_oligos_r1,
             n_oligos_r2,
             condition,
             r1,
             r2,
             name) {
        norm <-
            abs(length(which((
                data$ratio.x - data$ratio.y
            ) > 0)) - length(which((
                data$ratio.x - data$ratio.y
            ) < 0)))
        + abs(length(which((
            data$ratio.x - data$ratio.y
        ) > 0)) - length(which((
            data$ratio.x - data$ratio.y
        ) < 0)))
        + abs(length(which((
            data$ratio.x - data$ratio.y
        ) > 0)) - length(which((
            data$ratio.x - data$ratio.y
        ) < 0)))
        outs <- data.frame(
            Condition = condition,
            ReplicateA = r1,
            ReplicateB = r2,
            number_Oligos_ReplicateA = n_oligos_r1,
            number_Oligos_ReplicateB = n_oligos_r2,
            number_Oligos_Joined = data %>% nrow(),
            fraction_Oligos_ReplicateA = (data %>% nrow() / n_oligos_r1),
            fraction_Oligos_ReplicateB = (data %>% nrow() / n_oligos_r2),
            DNA_spearman = correlate(data$dna_normalized.x, data$dna_normalized.y, "spearman"),
            RNA_spearman = correlate(data$rna_normalized.x, data$rna_normalized.y, "spearman"),
            Ratio_spearman = correlate(data$ratio.x, data$ratio.y, "spearman"),
            DNA_pearson = correlate(data$dna_normalized.x, data$dna_normalized.y, "pearson"),
            RNA_pearson = correlate(data$rna_normalized.x, data$rna_normalized.y, "pearson"),
            Ratio_pearson = correlate(data$ratio.x, data$ratio.y, "pearson"),
            DNA_log2_pearson = correlate(
                data$dna_normalized_log2.x,
                data$dna_normalized_log2.y,
                "pearson"
            ),
            RNA_log2_pearson = correlate(
                data$rna_normalized_log2.x,
                data$rna_normalized_log2.y,
                "pearson"
            ),
            Ratio_log2_pearson = correlate(data$ratio_log2.x, data$ratio_log2.y, "pearson"),
            NormSymmetry = norm,
            stringsAsFactors = FALSE
        )
        return(outs)
    }

write_correlation_plots <- function(plots, name) {
    correlation_plots <- cowplot::plot_grid(plotlist = plots, ncol = 1)
    # correlation_plots <- do.call("grid.arrange", c(plots))

    ggplot2::ggsave(name,
                    correlation_plots,
                    width = 15,
                    height = 10 * length(plots))
}

write_correlation <- function(correlations, name) {
    write.table(
        correlations,
        file = name,
        quote = FALSE,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE
    )
}

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

if (data %>% nrow() > 1) {
    print("Pairwise comparisons")
    # make pairwise combinations
    selected <- combn(data$Replicate, 2)

    plots_correlations_rna <- list()
    plots_correlations_dna <- list()
    plots_correlations_ratio <- list()
    plots_cor_min_thresh_rna <- list()
    plots_cor_min_thresh_dna <- list()
    plots_cor_min_thresh_ratio <- list()
    stats_correlations <- data.frame()
    stats_cor_min_thresh <- data.frame()

    for (i in seq(1, dim(selected)[2])) {
        print(selected[, i])
        r1 <- selected[1, i]
        r2 <- selected[2, i]
        data1 <- all %>% filter(replicate == r1)
        data2 <- all %>% filter(replicate == r2)

        n_oligos_r1 <- data1 %>% nrow()
        n_oligos_r2 <- data2 %>% nrow()

        n_oligos_r1_thres <- data1 %>%
            filter(n_obs_bc >= thresh) %>%
            nrow()
        n_oligos_r2_thres <- data2 %>%
            filter(n_obs_bc >= thresh) %>%
            nrow()

        res <- data1 %>% inner_join(data2, by = c("name"))

        plots_correlations_dna[[i]] <-
            plot_correlations_dna(res, cond, r1, r2, "pairwise")
        plots_correlations_rna[[i]] <-
            plot_correlations_rna(res, cond, r1, r2, "pairwise")
        plots_correlations_ratio[[i]] <-
            plot_correlations_ratio(res, cond, r1, r2, "pairwise")

        stats_correlations <- stats_correlations %>%
            bind_rows(
                get_correlation_stats(
                    res,
                    n_oligos_r1,
                    n_oligos_r2,
                    cond,
                    r1,
                    r2,
                    "correlation"
                )
            )

        # Min Threshold
        res <-
            res %>% filter(n_obs_bc.x >= thresh, n_obs_bc.y >= thresh)
        plots_cor_min_thresh_dna[[i]] <-
            plot_correlations_dna(res, cond, r1, r2, "pairwise_minThreshold")
        plots_cor_min_thresh_rna[[i]] <-
            plot_correlations_rna(res, cond, r1, r2, "pairwise_minThreshold")
        plots_cor_min_thresh_ratio[[i]] <-
            plot_correlations_ratio(res, cond, r1, r2, "pairwise_minThreshold")

        stats_cor_min_thresh <- stats_cor_min_thresh %>%
            bind_rows(
                get_correlation_stats(
                    res,
                    n_oligos_r1_thres,
                    n_oligos_r2_thres,
                    cond,
                    r1,
                    r2,
                    "correlation_minThreshold"
                )
            )
    }

    write_correlation_plots(plots_correlations_dna,
                            sprintf("%s_DNA_pairwise.png", outdir))
    write_correlation_plots(plots_correlations_rna,
                            sprintf("%s_RNA_pairwise.png", outdir))
    write_correlation_plots(plots_correlations_ratio,
                            sprintf("%s_Ratio_pairwise.png", outdir))
    write_correlation_plots(
        plots_cor_min_thresh_dna,
        sprintf("%s_DNA_pairwise_minThreshold.png", outdir)
    )
    write_correlation_plots(
        plots_cor_min_thresh_rna,
        sprintf("%s_RNA_pairwise_minThreshold.png", outdir)
    )
    write_correlation_plots(
        plots_cor_min_thresh_ratio,
        sprintf("%s_Ratio_pairwise_minThreshold.png", outdir)
    )

    write_correlation(stats_correlations,
                      sprintf("%s_correlation.tsv", outdir))
    write_correlation(stats_cor_min_thresh,
                      sprintf("%s_correlation_minThreshold.tsv", outdir))
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
