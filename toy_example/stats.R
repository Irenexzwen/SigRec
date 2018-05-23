# This script is used to calculate the correlation of the result from smallCS.m 

#!/usr/bin env Rscript
args <- commandArgs(TRUE)
input_mat <- args[1]
library(parallel)
library(R.matlab)

ten <- readMat(input_mat)
RX<-ten$recoverX
RX<-round(RX)
RX[RX<0] <-0
OX<-read.table("pan_expression.txt")
OX<-as.matrix(OX,nrow=nrow(OX))

ccorr <- list()
for (i in 1:ncol(RX))
{
	  ccorr <- c(ccorr,cor(OX[,i],RX[,i]))
}

rcorr <- list()
	for (i in 1:nrow(RX))
{
	  rcorr <- c(rcorr,cor(OX[i,],RX[i,]))
}


ccorr <- sapply(ccorr,'[')
rcorr <- sapply(rcorr,'[')

ccorr<-ccorr[!is.na(ccorr)]
rcorr<-rcorr[!is.na(rcorr)]

rcor_median <- median(rcorr)
rcor_mean <- mean(rcorr)
ccor_median <- median(ccorr)
ccor_mean <- mean(ccorr)

name <- paste(gsub(".mat","",input_mat),"_rowmedian:",rcor_median,"_rowmean:",rcor_mean,"_colmedian:",ccor_median,"_colmean:",ccor_mean,"_.txt",sep="")

write.table(RX,name)
