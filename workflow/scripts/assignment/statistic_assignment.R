library(tidyverse)
library(optparse)


option_list <- list(
    make_option(c("-i", "--input"),
        type = "character",
        help = "Assignment parcode file"
    ),
    make_option(c("-p", "--plot"),
        type = "character",
        help = "histogram plot"
    ),
    make_option(c("-s", "--statistic"),
        type = "character",
        help = "Statistics"
    )
)

arguments <- parse_args(OptionParser(option_list = option_list), positional_arguments = TRUE)
opt <- arguments$options

if (!"input" %in% names(opt)) {
    stop("--input parameter must be provided. See script usage (--help)")
}
if (!"plot" %in% names(opt)) {
    stop("--plot parameter must be provided. See script usage (--help)")
}
if (!"statistic" %in% names(opt)) {
    stop("--statistic parameter must be provided. See script usage (--help)")
}

bcs <- read_delim(opt$input, delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

bcs <- bcs %>%
    filter(!(X2 %in% c("other", "ambiguous"))) %>%
    separate(X4, c("support", "all"), sep = "/", convert = TRUE)



avg_support <- mean(bcs$support / bcs$all)
median_support <- median(bcs$support / bcs$all)

bcs <- bcs %>%
    group_by(X2) %>%
    count()

oligos_support <- bcs %>%
    filter(n >= 15) %>%
    nrow()

output <- data.frame(Value = c(avg_support, median_support, oligos_support))

rownames(output) <- c("Average support of BC for Oligo:", "Median support of BC for Oligo:", "Oligos with >= 15 BCs:")

write.table(output, file = opt$statistic, quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)

p <- ggplot(bcs, aes(n)) +
    geom_histogram(binwidth = 1) +
    geom_vline(aes(xintercept = mean(n)), col = "red", size = 2) +
    geom_text(aes(
        y = 500, label = sprintf("\nmean = %f", mean(n)),
        x = max(n) / 2
    ), colour = "red") +
    geom_vline(aes(xintercept = median(n)), col = "blue", size = 2) +
    geom_text(aes(
        y = 500, label = sprintf("median = %f\n", median(n)),
        x = max(n) / 2
    ), colour = "blue") +
    theme_bw()

ggsave(opt$plot, p, width = 15, height = 12)
