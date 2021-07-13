library(getopt)
library(uuid)

opt = getopt(matrix(c(
  'inputMatrix', 'i', 1, "character",
  'inputStat', 's', 1, "character",
  'outputFile', 'o', 1, "character"
),ncol=4,byrow=TRUE));

data <- read.csv(file = opt$inputMatrix, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
stats <- read.csv(file = opt$inputStat, header = FALSE)

sample_base <- unlist(strsplit(opt$inputMatrix, split="/"))

v <- c()

plot_title <- paste(sample_base[length(sample_base)], "\n", stats[1,1], "%", " of reads have mate mapped to a different chromosome", sep="")

for(i in 1:nrow(data)){
  v[i] <- sum(data[i,])
}

# fix for https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows/issues/63
oldWD = getwd()
setwd(dirname(opt$outputFile))
tmpFileName=paste0(UUIDgenerate(),".tmp")

pdf(tmpFileName, height = 10, width = 10)
  barplot(v, names.arg = rownames(data), col = rainbow(nrow(data)),ylab = "# reads", xlab = "Chromosome", main = plot_title)
dev.off()

file.rename(tmpFileName, opt$outputFile)
setwd(oldWD)