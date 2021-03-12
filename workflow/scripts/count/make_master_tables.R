#Adapted from Vikram Agarwal by Gracie Gordon

cpath <- grep('conda', .libPaths(), value=TRUE, ignore.case=TRUE)
.libPaths(cpath)

library(dplyr)
library(optparse)


option_list <- list(
    make_option(c("-c", "--condition"), type="character",
        help="Condition name"),
    make_option(c("-l", "--label"), type="character",
        help="Label file. (optional)"),
    make_option(c("-f", "--files"), type="character",
        help="Comma separated input files of assigned counts"),
    make_option(c("-r", "--replicates"), type="character",
        help="Comma separated name of the replicates (same order than files)"),
    make_option(c("-a", "--output-all"), type="character",
        help="Output file of master table. No BC threshold filter (optional)."),
    make_option(c("-o", "--output"), type="character",
        help="Output file of master table filtered by --threshold"),
    make_option(c("-s", "--statistic"), type="character",
        help="Statistic of master table and filtered master table"),
    make_option(c("-t", "--threshold"), type="integer", default=10,
        help="Number of required barcodes (default 10)")
)

arguments <- parse_args(OptionParser(option_list=option_list), positional_arguments=TRUE)
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
if (!"output" %in% names(opt)) {
  stop("--output parameter must be provided. See script usage (--help)")
}
if (!"statistic" %in% names(opt)) {
  stop("--statistic parameter must be provided. See script usage (--help)")
}



#exp=args[1]
cond=opt$cond
thresh=opt$threshold
files <- strsplit(opt$files,",")[[1]]
replicates=strsplit(opt$replicates,",")[[1]]
if (length(files) != length(replicates)) {
    stop("Number of input files must be euqal to number of replicates")
}

outfile=opt$output
avg_outfile=opt$statistic
#out=args[3]

##MAKE MASTER TABLE
masterTable = data.frame()
for (i in 1:length(files)){
   file=files[i]
   rep=replicates[i]

   table = as.data.frame(read.table(file,header=TRUE),stringsAsFactors=FALSE)
   table$condition = cond
   table$replicate = rep

   masterTable <- masterTable %>% bind_rows(table)
}

masterTable <- masterTable %>% group_by(condition,replicate) %>%
                select(condition, replicate, name, dna_counts, rna_counts, dna_normalized, rna_normalized, ratio, log2, n_obs_bc)

masterTableFiltered <- masterTable %>%
                      filter(n_obs_bc >= thresh) %>%
                      mutate(ratio = ratio/median(ratio), log2 = round(log2(ratio),8))
masterTable <- masterTable %>% mutate(ratio = ratio/median(ratio), log2 = round(log2(ratio),8))

writeFile <- function(file,table) {
  gz <- gzfile(file, "w")
  write.table(table,file=gz,quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE )
  close(gz)
}

writeFile(outfile,masterTableFiltered)

if ("output-all" %in% names(opt)) {
  writeFile(opt$`output-all`, masterTable)
}

##MAKE AVERAGED ACROSS REPLICATES
makeAverageAcrossReplicates <- function(table, name) {
  avg <- table %>% summarize(
                      mean_ratio=mean(ratio),
                      mean_log2=log2(mean(ratio)),
                      mean_n_obs_bc=mean(n_obs_bc)
                    )
  avg$BC_filter = name

  return(avg)
}

all_avg <- makeAverageAcrossReplicates(masterTable, 'None') %>%
            bind_rows(makeAverageAcrossReplicates(masterTableFiltered, paste0("n_obs_bc >= ", thresh)))

writeFile(avg_outfile, all_avg)


print('done')
