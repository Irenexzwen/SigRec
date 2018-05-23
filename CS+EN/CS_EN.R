#	This script applies elastic net as regression method for genes with poor performance in compressed sensing
#
#	input: Output of star.m 
#	output: recovered expression profile 


#!/usr/bin env Rscript
args <- commandArgs(TRUE)
input_mat <- args[1]

library(parallel)
library(R.matlab)
library(glmnet)
library(caret)

#read in .mat file which is the output of compressed sensing
cat(input_mat,'\n')
m <- readMat(input_mat)
XX0 = m$X # list of sub original X
M = m$M # list of pool design matrix for each sub block 
R = m$R # list of sub recovered X
rm(m)

len<-length(XX0)            #equals the number of sub blocks 
cat("load finished",'\n')

# pick out the "dense" genes from compressed sensing reults and use elastic net for regression again
RX = matrix(nrow=0,ncol=ncol(R[[1]][[1]])) 
OX = matrix(nrow=0,ncol=ncol(R[[1]][[1]])) 

#parallel 
coef_list <- mclapply(1:len,FUN=function(i)
  {
  XX <- R[[i]][[1]] #block of Recover X
  
  XX[XX<0] <- 0     # discard negative recovered expression level 
  XX <- round(XX)    

  r=apply(XX, 2,function(x){sum(x>0)})  # find dense genes
  bb <- nrow(M[[1]][[1]])               # choose threshold 
  dense=which(r>=bb)
  X = XX0[[i]][[1]][,dense]

  Y = M[[i]][[1]] %*% X 
  ncolM <- ncol(M[[i]][[1]])
  
  coeff <- matrix(nrow=ncolM,ncol=0)
  alpha_grid <- seq(0,1,0.2)
  trnCtrl <- trainControl(method = "repeatedCV",number = 5,repeats = 3)
  namex <- paste('c',1:ncol(M[[i]][[1]]),sep="")
  colnames(M[[i]][[1]]) <- namex

  for (j in 1:length(dense)){
   if(any(Y[,j]==0)){coeff <- cbind(coeff,X[,j])} 
   else{
    g1 <- cv.glmnet(M[[i]][[1]],Y[,j],nfolds = 3,alpha=0.5,lower.limits = 0)
	srchGrid <- expand.grid(.alpha = alpha_grid, .lambda = g1$lambda.min)
    my_train <- train(M[[i]][[1]], Y[,j],method = "glmnet",tuneGrid = srchGrid,trControl = trnCtrl) 
	alpha_caret <- as.numeric(my_train$bestTune[1])
    nn <- glmnet(M[[i]][[1]],Y[,j],family = "gaussian",alpha = alpha_caret,lower.limits = 0,lambda = g1$lambda.min) 
    coeff <- cbind(coeff,coef(nn,s=g1$lambda.min)[2:(ncolM+1),])
    rm(g1)
    rm(nn)}
  }
  
  XX[,dense]=coeff
  return(XX)
},mc.cores = len)


# row bind orginal and recovered X
for (i in 1:len){
	  RX <- rbind(RX,coef_list[[i]])
	  OX <- rbind(OX,XX0[[i]][[1]])
}


rownames(OX)<-NULL
rownames(RX)<-NULL
colnames(OX)<-NULL
colnames(RX)<-NULL

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

name <- paste("rowmedian_",rcor_median,"_rowmean_",rcor_mean,"_colmedian_",ccor_median,"_colmean_",ccor_mean,"_.txt",sep="")


RX <- t(RX)

drop <- as.matrix(read.table("drop.txt"))
drop <- drop[3:length(drop)]
colnames(RX) <- drop


genenames <- as.matrix(read.table("genename.txt"))
nalinenu <- as.vector(is.na(genenames))
genenames <- genenames[!nalinenu]
rownames(RX) <- genenames

write.table(RX,name)


