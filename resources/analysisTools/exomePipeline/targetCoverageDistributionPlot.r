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
cmdArgs = commandArgs(TRUE)
cmd <- paste("awk '{print $NF}'", cmdArgs[1], sep = " ")
covs <- read.table(pipe(cmd))
title <- cmdArgs[3]
library(Hmisc)
xmax <- quantile(covs$V1, c(0.9))
png(cmdArgs[2])
Ecdf(covs$V1, what="1-F", xlab = 'Coverage', q=c(.25,.5,.75), xlim=c(0,xmax), main = title)
dev.off()
