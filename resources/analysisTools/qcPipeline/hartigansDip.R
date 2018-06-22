#!/usr/bin/env Rscript
#
# This file is part of the AlignmentAndQCWorkflow plugin.
#
# This script is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License.
#
# This script is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this script.  If not, see <https://www.gnu.org/licenses/>.
#
# install.packages("diptest")
# Run 
# Rscript DipStatistcs.R   5e4bbb6b-66b2-4787-b8ce-70d17bc80ba8.MultipleMetrics.insert_size_metrics
require(diptest)
require(ggplot2)
print("The dip test measures multimodality in a sample by the maximum difference, over all sample points, between the empirical distribution function, and the unimodal distribution function that minimizes that maximum difference")
args = commandArgs(trailingOnly=T)
mydat = read.delim(args[1], skip=11)
datDir = dirname(args[1])
head(mydat)
sampleName = args[2]
dat = sample(rep(mydat$insert_size, mydat$All_Reads.fr_count), 10000)
dipTest = (d.t <- dip.test(dat))
dipStat = round(dipTest$statistic,3)
dipPval = round(dipTest$p.value,3)

qplot(dat, geom="density", main=paste(sampleName,"\nInsert size density distr. Hartigans' dip statistic, Dn"), xlab="insert size")  + 
  geom_rug() + 
  annotate("text", label = paste("Test for unimodality\nDn =", dipStat, "\npval =",dipPval), x = mean(dat)*1.5, y = max(density(dat)$y)*0.8, size = 6, colour = " dodgerblue4", alpha=0.8)


fileOut = paste(datDir, "/",sampleName, "_HartigansDip.txt", sep="")

ggsave(paste(datDir, "/", sampleName, "_HartigansDip_densityPlot.png", sep=""))

datOut = paste("Name\tDn\t\tpval\n",sampleName,"\t", dipStat, "\t", dipPval, sep="")
write.table(datOut, fileOut , quote=F, col.names = F, row.names=F)
