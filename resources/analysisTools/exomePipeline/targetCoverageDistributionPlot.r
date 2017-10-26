cmdArgs = commandArgs(TRUE)
cmd <- paste("awk '{print $NF}'", cmdArgs[1], sep = " ")
covs <- read.table(pipe(cmd))
title <- cmdArgs[3]
library(Hmisc)
xmax <- quantile(covs$V1, c(0.9))
png(cmdArgs[2])
Ecdf(covs$V1, what="1-F", xlab = 'Coverage', q=c(.25,.5,.75), xlim=c(0,xmax), main = title)
dev.off()
