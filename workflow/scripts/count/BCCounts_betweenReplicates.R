
cpath <- grep('conda', .libPaths(), value=TRUE, ignore.case=TRUE)
.libPaths(cpath)

library(optparse)
library(tidyverse)
#library(viridis)


option_list <- list(
    make_option(c("-c", "--condition"), type="character",
        help="Condition name"),
    make_option(c("-f", "--files"), type="character",
        help="Comma separated input files of assigned counts"),
    make_option(c("-r", "--replicates"), type="character",
        help="Comma separated name of the replicates (same order than files)"),
    make_option(c("-o", "--outfile"), type="character",
        help="Output file of the correlation table.")
)

parser <- OptionParser(option_list=option_list)
arguments <- parse_args(parser, positional_arguments=TRUE)
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
if (!"outfile" %in% names(opt)) {
  stop("--outfile parameter must be provided. See script usage (--help)")
}

# condition
cond=opt$condition
# outfile
outfile=opt$outfile
# replicates and count files
files <- strsplit(opt$files,",")[[1]]
replicates=strsplit(opt$replicates,",")[[1]]
if (length(files) != length(replicates)) {
    stop("Number of input files must be euqal to number of replicates")
}

num_replicates=length(replicates)

data <- data.frame(File=files,Replicate=replicates)
data$Condition <- cond

print(data)





getOverlapStats <- function(data1,data2,condition,r1,r2){
  sData1 <- data1 %>% summarize(size=n(), count=sum(Counts))
  print(sData1)
  sData2 <- data2 %>% summarize(size=n(), count=sum(Counts))
  print(sData2)
  data <- data1 %>% inner_join(data2,by=c('Barcode'))
  sData <- data %>% summarize(size=n(), count1=sum(Counts.x),count2=sum(Counts.y))
  print(sData)
  print(sData$count1)
  print(sData$count1[1])
  outs <- data.frame(
                Condition = condition,
                ReplicateA = r1,
                ReplicateB = r2,
                BCs_ReplicateA = sData1$size,
                BCs_ReplicateB = sData2$size,
                BCs_Overlap = sData$size,
                Counts_ReplicateA = sData1$count,
                Counts_ReplicateB = sData2$count,
                Counts_OverlapA = sData$count1,
                Counts_OverlapB = sData$count2,
                Overlap_BC_ReplicateA = (sData$size / sData1$size),
                Overlap_BC_ReplicateB = (sData$size / sData2$size),
                Overlap_Counts_ReplicateA = (sData$count1 / sData1$count),
                Overlap_Counts_ReplicateB = (sData$count2 / sData2$count),
                Lincoln_Peterson_estimator = round(sData1$size * (sData2$size/sData$size)),
                stringsAsFactors=FALSE)
  return(outs)

}

readData <- function(file) {
  data <- read.table(file,as.is=T,header=F,stringsAsFactors = F)
  colnames(data) <- c("Barcode","Counts")
  return(data)
}

writeOutput <- function(correlations, name){
    write.table(correlations,file=name,quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE )
}

# pairwise comparison only if more than one replicate
if(data %>% nrow >1){

  # make pairwise combinations
  selected <- combn(data$Replicate,2)
  print('sel')
  print(selected)

  overlapStats = data.frame()
  print('reps')
  for(i in seq(1,dim(selected)[2])){
    print(selected[,i])
    r1=selected[1,i]
    r2=selected[2,i]
    data1 <- readData(as.character((data %>% filter(Replicate == r1))$File))
    print("Read data 1")
    data2 <- readData(as.character((data %>% filter(Replicate == r2))$File))
    print("Read data 2")
    
    overlapStats_new <- getOverlapStats(data1, data2, cond,r1,r2)
    overlapStats <- overlapStats %>% bind_rows(overlapStats_new)
  }
  overlapStats$Mean_Lincoln_Peterson_estimator <- mean(overlapStats$Lincoln_Peterson_estimator)
  writeOutput(overlapStats, outfile)
}
