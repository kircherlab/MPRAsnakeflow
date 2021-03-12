
library(optparse)
library(dplyr)
library(readr)


option_list <- list(
    make_option(c("-c", "--count"), type="character",
        help="Statistic file of Barcode and UMI counts"),
    make_option(c("-s", "--shared"), type="character",
        help="Statistic file of shared RNA and DNA barcodes"),
    make_option(c("-o", "--output"), type="character",
        help="Output file of final count stats")
)

parser <- OptionParser(option_list=option_list)
arguments <- parse_args(parser, positional_arguments=TRUE)
opt <- arguments$options

if (!"count" %in% names(opt)) {
  stop("--count parameter must be provided. See script usage (--help)")
}
if (!"shared" %in% names(opt)) {
  stop("--shared parameter must be provided. See script usage (--help)")
}
if (!"output" %in% names(opt)) {
  stop("--output parameter must be provided. See script usage (--help)")
}



stats <- read_delim(opt$count, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

stats_shared <- read_delim(opt$shared, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

umiLength=16

colnames(stats) <- c("Condition","Replicate", "Type", "Reads","Barcodes x UMI", "Barcodes", "Unique UMIs")
colnames(stats_shared) <- c("Condition","Replicate", "Barcodes shared RNA&DNA")


stats <- stats %>% full_join(stats_shared)

stats <- stats %>%
  mutate(`Read Duplicates`=1-(`Barcodes x UMI`/Reads),`UMIs/Barcode`=`Barcodes x UMI`/Barcodes, RatioUMIPossible=`Unique UMIs`/pmin(`Barcodes x UMI`,4**umiLength), `Fraction barcodes shared`=`Barcodes shared RNA&DNA`/Barcodes) %>%
  select(1,2,3,4,5,9,6,10,7,11,8,12)

write.table(stats, file=opt$output, quote=FALSE, sep='\t',row.names = FALSE)
