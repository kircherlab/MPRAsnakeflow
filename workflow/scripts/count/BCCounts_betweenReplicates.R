# This script calculates the overlap statistics between replicates of barcode counts.
# It uses the optparse library to handle command-line arguments and the tidyverse library for data manipulation.
#
# Command-line arguments:
# -c, --condition: Condition name (character)
# -f, --files: Comma-separated input files of assigned counts (character)
# -r, --replicates: Comma-separated names of the replicates (same order as files) (character)
# -o, --outfile: Output file of the correlation table (character)
#
# The script performs the following steps:
# 1. Parse command-line arguments and check for required arguments.
# 2. Read the input files and replicates.
# 3. If the number of files does not match the number of replicates, stop with an error.
# 4. Create a data frame with the files, replicates, and condition.
# 5. Define functions to:
#    - Calculate overlap statistics between two data sets (get_overlap_stats).
#    - Read data from a file (read_data).
#    - Write the output to a file (write_output).
# 6. If there is more than one replicate, perform pairwise comparisons:
#    - Generate all pairwise combinations of replicates.
#    - For each pair, read the data, calculate overlap statistics, and store the results.
#    - Calculate the mean Lincoln-Peterson estimator and write the results to the output file.
# 7. If there is only one replicate, create a single row with NA values for overlap statistics and write it to the output file.
library(optparse)
library(tidyverse)
# library(viridis)


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
  make_option(c("-o", "--outfile"),
    type = "character",
    help = "Output file of the correlation table.",
    metavar = "character"
  )
)

parser <- OptionParser(option_list = option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

required_args <- c("condition", "files", "replicates", "outfile")
missing_args <- setdiff(required_args, names(opt))

if (length(missing_args) > 0) {
  stop(paste("Missing required arguments:", paste(missing_args, collapse = ", "), ". See script usage (--help)"))
}

# condition
cond <- opt$condition
# outfile
outfile <- opt$outfile
# replicates and count files
files <- strsplit(opt$files, ",")[[1]]
replicates <- strsplit(opt$replicates, ",")[[1]]
if (length(files) != length(replicates)) {
  stop("Number of input files must be euqal to number of replicates")
}

num_replicates <- length(replicates)

data <- data.frame(File = files, Replicate = replicates)
data$Condition <- cond


get_overlap_stats <- function(data1, data2, condition, r1, r2) {
  s_data1 <- data1 %>% summarize(size = n(), count = sum(Counts))
  print(s_data1)
  s_data2 <- data2 %>% summarize(size = n(), count = sum(Counts))
  print(s_data2)
  data <- data1 %>% inner_join(data2, by = c("Barcode"))
  s_data <- data %>% summarize(size = n(), count1 = sum(Counts.x), count2 = sum(Counts.y))
  print(s_data)
  print(s_data$count1)
  print(s_data$count1[1])
  outs <- data.frame(
    Condition = condition,
    ReplicateA = r1,
    ReplicateB = r2,
    BCs_ReplicateA = s_data1$size,
    BCs_ReplicateB = s_data2$size,
    BCs_Overlap = s_data$size,
    Counts_ReplicateA = s_data1$count,
    Counts_ReplicateB = s_data2$count,
    Counts_OverlapA = s_data$count1,
    Counts_OverlapB = s_data$count2,
    Overlap_BC_ReplicateA = (s_data$size / s_data1$size),
    Overlap_BC_ReplicateB = (s_data$size / s_data2$size),
    Overlap_Counts_ReplicateA = (s_data$count1 / s_data1$count),
    Overlap_Counts_ReplicateB = (s_data$count2 / s_data2$count),
    Lincoln_Peterson_estimator = round(s_data1$size * (s_data2$size / s_data$size)),
    stringsAsFactors = FALSE
  )
  return(outs)
}

read_data <- function(file) {
  data <- read.table(file, sep = "\t", as.is = T, header = F, stringsAsFactors = F)
  colnames(data) <- c("Barcode", "Counts")
  return(data)
}

write_output <- function(correlations, name) {
  write.table(correlations, file = name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

# pairwise comparison only if more than one replicate
if (data %>% nrow() > 1) {
  # make pairwise combinations
  selected <- combn(data$Replicate, 2)
  print("sel")
  print(selected)

  overlap_stats <- data.frame()
  print("reps")
  for (i in seq(1, dim(selected)[2])) {
    print(selected[, i])
    r1 <- selected[1, i]
    r2 <- selected[2, i]
    data1 <- read_data(as.character((data %>% filter(Replicate == r1))$File))
    print("Read data 1")
    data2 <- read_data(as.character((data %>% filter(Replicate == r2))$File))
    print("Read data 2")

    overlap_stats_new <- get_overlap_stats(data1, data2, cond, r1, r2)
    overlap_stats <- overlap_stats %>% bind_rows(overlap_stats_new)
  }
  overlap_stats$Mean_Lincoln_Peterson_estimator <- mean(overlap_stats$Lincoln_Peterson_estimator)
  write_output(overlap_stats, outfile)
} else {
  # If only one replicate, create a single row with NA values for overlap statistics
  single_replicate_stats <- data.frame(
    Condition = cond,
    ReplicateA = replicates[1],
    ReplicateB = NA,
    BCs_ReplicateA = NA,
    BCs_ReplicateB = NA,
    BCs_Overlap = NA,
    Counts_ReplicateA = NA,
    Counts_ReplicateB = NA,
    Counts_OverlapA = NA,
    Counts_OverlapB = NA,
    Overlap_BC_ReplicateA = NA,
    Overlap_BC_ReplicateB = NA,
    Overlap_Counts_ReplicateA = NA,
    Overlap_Counts_ReplicateB = NA,
    Lincoln_Peterson_estimator = NA,
    Mean_Lincoln_Peterson_estimator = NA,
    stringsAsFactors = FALSE
  )
  write_output(single_replicate_stats, outfile)
}
