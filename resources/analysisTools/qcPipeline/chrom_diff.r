#!/usr/bin/env Rscript
#
# This file is part of the AlignmentAndQCWorkflow plugin.
#
# This script is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 or 3 of the License.
#
# This script is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this script.  If not, see <https://www.gnu.org/licenses/>.
#

library(getopt)

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

pdf(opt$outputFile, height = 10, width = 10)

barplot(v, names.arg = rownames(data), col = rainbow(nrow(data)),ylab = "# reads", xlab = "Chromosome", main = plot_title)

dev.off()
