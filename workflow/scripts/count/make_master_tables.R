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
for (i in seq_len(length(files))) {
  file <- files[i]
  rep <- replicates[i]

  table <- as.data.frame(read.table(file, header = TRUE), stringsAsFactors = FALSE)
  table$replicate <- rep

  master_table <- master_table %>% bind_rows(table)
}

master_table <- master_table %>%
  group_by(replicate) %>%
  filter(!oligo_name %in% c("no_BC")) %>%
  mutate(
    dna_normalized = round(dna_normalized, precision),
    rna_normalized = round(rna_normalized, precision),
    ratio = round(ratio, precision),
    log2FoldChange = round(log2FoldChange, precision)
  )

master_table_filtered <- master_table %>%
  filter(n_bc >= thresh)
# TODO use this? Maybe not because it normalizes across all replicates
# master_table_filtered %>% mutate(ratio = ratio / median(ratio))

# TODO use this? Maybe not because it normalizes across all replicates
# master_table %>% mutate(ratio = ratio / median(ratio))

write_file <- function(file, table) {
  gz <- gzfile(file, "w")
  write.table(
    table,
    file = gz, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
  )
  close(gz)
}

write_file(
  outfile,
  master_table_filtered %>%
    select(replicate, oligo_name, dna_counts, rna_counts, dna_normalized, rna_normalized, log2FoldChange, n_bc)
)

if ("output-all" %in% names(opt)) {
  write_file(
    opt$`output-all`,
    master_table %>%
      select(replicate, oligo_name, dna_counts, rna_counts, dna_normalized, rna_normalized, log2FoldChange, n_bc)
  )
}

## MAKE AVERAGED ACROSS REPLICATES
make_average_across_replicates <- function(table, name) {
  avg <- table %>% summarize(
    mean_ratio = mean(ratio),
    mean_log2FoldChange = mean(log2FoldChange),
    mean_n_bc = mean(n_bc)
  )
  avg$BC_filter <- name

  return(avg)
}

all_avg <- make_average_across_replicates(master_table, "None") %>%
  bind_rows(make_average_across_replicates(master_table_filtered, paste0("n_bc >= ", thresh)))

write_file(avg_outfile, all_avg)

print("done")
