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
  count() %>%
  ungroup()

oligos_support <- bcs %>%
  filter(n >= 15) %>%
  nrow()

output <- data.frame(Value = c(avg_support, median_support, oligos_support))

rownames(output) <- c("Average support of BC for Oligo:", "Median support of BC for Oligo:", "Oligos with >= 15 BCs:")

write.table(output, file = opt$statistic, quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)

p <- ggplot(bcs, aes(n)) +
  geom_histogram(bins = 150) +
  xlab("Number of BCs") +
  ylab("Number of Oligos") +
  geom_vline(aes(xintercept = mean(bcs$n)), col = "red", linewidth = 2) +
  geom_text(
    aes(
      label = sprintf("\nmean = %.2f", mean(bcs$n)),
      x = Inf, y = Inf
    ),
    hjust = 1.1, vjust = 1.3, colour = "red"
  ) +
  geom_vline(aes(xintercept = median(bcs$n)), col = "blue", linewidth = 2) +
  geom_text(
    aes(
      label = sprintf("median = %.2f\n", median(bcs$n)),
      x = Inf, y = Inf
    ),
    hjust = 1.1, vjust = 1.3, colour = "blue"
  ) +
  theme_bw()

ggsave(opt$plot, p, type = "cairo", dpi = 300)
