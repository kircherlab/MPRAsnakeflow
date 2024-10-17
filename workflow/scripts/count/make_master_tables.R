# Adapted from Vikram Agarwal by Gracie Gordon and Max Schubach


library(dplyr)
library(optparse)


option_list <- list(
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
  make_option(c("-a", "--output-all"),
    type = "character",
    help = "Output file of master table. No BC threshold filter (optional)."
  ),
  make_option(c("-o", "--output"),
    type = "character",
    help = "Output file of master table filtered by --threshold"
  ),
  make_option(c("-s", "--statistic"),
    type = "character",
    help = "Statistic of master table and filtered master table"
  ),
  make_option(c("-t", "--threshold"),
    type = "integer", default = 10,
    help = "Number of required barcodes (default 10)"
  )
)

arguments <- parse_args(OptionParser(option_list = option_list), positional_arguments = TRUE)
opt <- arguments$options

if (!"files" %in% names(opt)) {
  stop("--files parameter must be provided. See script usage (--help)")
}
if (!"replicates" %in% names(opt)) {
  stop("--replicates parameter must be provided. See script usage (--help)")
}
if (!"output" %in% names(opt)) {
  stop("--output parameter must be provided. See script usage (--help)")
}
if (!"statistic" %in% names(opt)) {
  stop("--statistic parameter must be provided. See script usage (--help)")
}


# exp=args[1]
cond <- opt$cond
thresh <- opt$threshold
files <- strsplit(opt$files, ",")[[1]]
replicates <- strsplit(opt$replicates, ",")[[1]]
if (length(files) != length(replicates)) {
  stop("Number of input files must be euqal to number of replicates")
}

outfile <- opt$output
avg_outfile <- opt$statistic

precision <- 4

## MAKE MASTER TABLE
master_table <- data.frame()
for (i in seq_len(files)) {
  file <- files[i]
  rep <- replicates[i]

  table <- as.data.frame(read.table(file, header = TRUE), stringsAsFactors = FALSE)
  table$replicate <- rep

  master_table <- master_table %>% bind_rows(table)
}

master_table <- master_table %>%
  group_by(replicate) %>%
  filter(!oligo_name %in% c("no_BC")) %>%
  select(replicate, oligo_name, dna_counts, rna_counts, dna_normalized, rna_normalized, log2FoldChange, n_bc) %>%
  mutate(
    dna_normalized = round(dna_normalized, precision),
    rna_normalized = round(rna_normalized, precision),
    log2FoldChange = round(log2FoldChange, precision)
  )

# TODO use this? Maybe not because it normalizes across all replicates
# mutate(ratio = ratio / median(ratio))

master_table_filtered <- master_table %>%
  filter(n_bc >= thresh) %>%
  mutate(ratio = ratio / median(ratio), log2FoldChange = round(log2(ratio), 8))
master_table <- master_table %>% mutate(ratio = ratio / median(ratio), log2FoldChange = round(log2(ratio), 8))

write_file <- function(file, table) {
  gz <- gzfile(file, "w")
  write.table(table, file = gz, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  close(gz)
}

write_file(outfile, master_table_filtered)

if ("output-all" %in% names(opt)) {
  writeFile(opt$`output-all`, master_table)
}

## MAKE AVERAGED ACROSS REPLICATES
make_average_across_replicates <- function(table, name) {
  avg <- table %>% summarize(
    mean_ratio = mean(ratio),
    mean_log2FoldChange = log2(mean(ratio)),
    mean_n_bc = mean(n_bc)
  )
  avg$BC_filter <- name

  return(avg)
}

all_avg <- make_average_across_replicates(master_table, "None") %>%
  bind_rows(make_average_across_replicates(master_table_filtered, paste0("n_bc >= ", thresh)))

writeFile(avg_outfile, all_avg)


print("done")
