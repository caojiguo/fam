# set working directory
# setwd("/Users/cao/Dropbox/Research/Peijun/FAM/nfam2017Feb19/r code/pm 2.5 analysis")
pm25 <- read.csv("pm25_2000.csv")
#windows()
#pdf("profiles of pm2.5 of 101 cities.pdf")
par(pch="o", font.axis=2, font.lab=2, cex.axis=1.2, cex.lab=1.2)
matplot(x=1:366, y=t(pm25[,3:368]),type="l", lty=1, col=1, xlab="Day", ylab="PM 2.5")
# dev.off()
#======================================================

#data analysis
source("smoothing.R")

library(fda)
library(MASS)
library(cosso)
library(caTools)
#library(bisoreg)
library(mgcv)
#library(mda)
library(earth)
#library(fdapace)
#library(PACE)
library(R.matlab)
library(gglasso)
library(zoo)
library(splines)

# FPCA using PACE in Matlab
# The PACE code is provided in the folder "matlab code for PACE"
# Use the PACE code in Matlab to obtain results for the train data set
# The PACE results on the train data set are stored in train.mat

Res <- readMat("train.mat")
#PACE part
trainRes <- Res$trainRes
#pc score

ntrain <- 80 #80 cities in training data

trainPCscore <- matrix(unlist(trainRes[[6]]),nrow=ntrain)

#eigenvalue
lamhat <- unlist(trainRes[[3]])
#length(lamhat)


# Use the PACE code in Matlab to obtain results for the test data set
# The PACE results on the test data set are stored in test.mat
pred.test <- readMat("test.mat")
testPCscore <- pred.test$testPCscore
ntest <- nrow(testPCscore)

zeta.train <- pnorm(c(trainPCscore),sd=sqrt(rep(lamhat, each=ntrain)))
zeta.test <- pnorm(c(testPCscore),sd=sqrt(rep(lamhat, each=ntest)))
zeta.train <- matrix(zeta.train, nrow=ntrain)
zeta.test <- matrix(zeta.test, nrow=ntest)
P <- ncol(zeta.train)

head(zeta.train);
head(zeta.test)

ID <- pm25$ID  #indicate whether is training (ID=0) and test data (ID=1)
#resp date rate for 2000
train_res <- pm25[ID==0,2]
test_res <- pm25[ID==1,2]
hist(train_res)


#=====================================================
cat("MARS\n")
#mars with 18 FPCs
marsModel <- earth(y=train_res, x= zeta.train, degree=1) # build model
# plot(evimp(marsModel))
mars.pred <- predict(marsModel, newdata=zeta.test)
MSPE.mars <- mean((mars.pred - test_res)^2)
#1.676306e-08

#==================================
cat("saturated FAM\n")  #saturated FAM
full.data.train <- data.frame(cbind(train_res, zeta.train))
full.data.test <- data.frame(cbind(test_res, zeta.test))
colnames(full.data.train) <- c("resp",paste("x", 1:P, sep=""))
colnames(full.data.test) <- c("resp",paste("x", 1:P, sep=""))
#fam with 16 FPCs
predic16 <- paste( "s(", paste(paste("x", 1:P, sep=""), paste("k", "=4", sep=""), sep=","), ")", sep="",collapse=" + ")
# predic16 <- paste( "s(", paste("x", 1:P, sep=""),")", sep="",collapse=" + ")
fam <- gam(as.formula(paste("resp ~", predic16)),data=full.data.train, method="REML")
ypred <- predict(fam, full.data.test)
MSPE.fam <- mean((ypred - test_res)^2)
#1.766624e-08






#=================================
set.seed(123456)

cat("CSEFAM\n") #proposed by Hongxiao Zhu et al.
#using 20 FPC
G.Obj=cosso(x=zeta.train,y=train_res,family="Gaussian")
#windows()
tun1 <- tune.cosso(G.Obj,5,TRUE)
#the optimal one tun1$OptM = 11.45
#fit.sig <- predict.cosso(G.Obj,M=M0,type="nonzero")
ypred <- predict.cosso(G.Obj,xnew=zeta.test,M=tun1$OptM,type="fit")
# predict.cosso(G.Obj, M=tun1$OptM, type="nonzero")
MSPE.cosso <- mean((test_res-ypred)^2)
#1.607642e-08



#===========================================

#our proposed approach
cat("SSFAM\n")
#step 1 group lasso using cubic B-spline
M <- 15
Bbasis <- create.bspline.basis(c(0, 1), norder=4, nbasis=M)
#plot(Bbasis)
Z.train <- matrix(0, nrow=ntrain, ncol=P*M)
Z.test <- matrix(0, nrow=ntest, ncol=P*M)
for(j in 1:P)
{ #j = 1
  index <- ((j-1)*M+1):(j*M)
  #a <- eval.basis(zeta.est[,j], Bbasis); dim(a); dim(Z[,index])
  Z.train[,index] <- eval.basis(zeta.train[,j], Bbasis)
  Z.test[,index] <- eval.basis(zeta.test[,j], Bbasis)
}

mean.vec <- apply(Z.train, 2, mean); length(mean.vec)
Z.train <- Z.train - matrix(rep(mean.vec, times=ntrain), nrow=ntrain, byrow=T)
mean.vec <- apply(Z.test, 2, mean); length(mean.vec)
Z.test <- Z.test - matrix(rep(mean.vec, times=ntest), nrow=ntest, byrow=T)

G <- rep(1:P, each=M)
cv <- cv.gglasso(x=Z.train, y=train_res, group=G, loss="ls", pred.loss="L2")
#windows()
plot(cv)
fit <- gglasso(x=Z.train, y=train_res, group=G, lambda=3*1e-06)
#==================
#step 2 adaptive group lasso
w <- numeric(P)
for(j in 1:P)
{
  index <- (M*(j-1)+1):(M*j)
  w[j] <- abs(1/sqrt(sum((fit$beta[index])^2)))
}
w[w==Inf] <- 1e8

cv <- cv.gglasso(x=Z.train, y=train_res, group=G, loss="ls", pred.loss="L2",pf=w)
#windows()
plot(cv)

fit <- gglasso(x=Z.train, y=train_res, group=G, lambda=cv$lambda.min, pf=w)
#which(abs(fit$beta)>0) 1.503872e-10
ypred <- predict(fit,newx=Z.test,type="link")
MSPE.grp <- mean((ypred - test_res)^2)
# 1.677994e-08





#============================
#define weight vector and active sets after adaptive group lasso
w <- numeric(P)
act.set <- numeric(P*M)
numberBasis <- numeric(P)
for(j in 1:P)
{
  index <- (M*(j-1)+1):(M*j)
  w[j] <- abs(1/sqrt(sum((fit$beta[index])^2)))
  act.set[index] <- as.numeric(sum((fit$beta[index])^2) > 0)
  numberBasis[j] <- sum(fit$beta[index]!= 0)
}
w <- w[w!=Inf]
numberBasis <- numberBasis[numberBasis > 0]
w <- w*numberBasis/10
#length(w); table(act.set); length(numberBasis)

#penalize roughness


# Use quadrature to get integral - Composite Simpson's Rule
quadpts <- c(0, Bbasis$params, 1)
nquadpts <- length(quadpts)
quadwts <- as.vector(c(1,rep(c(4,2),(nquadpts-2)/2),4,1),mode="any")
quadwts <- c(1,rep(c(4,2),(nquadpts-1)/2))
quadwts[nquadpts] <- 1
delta <- diff(quadpts)[1]
quadwts <- quadwts*delta/3

# Second derivative of basis functions at quadrature points
Q2basismat   = eval.basis(quadpts,Bbasis,2);
# estimates for basis coefficients
R = t(Q2basismat)%*%(Q2basismat*(quadwts%*%t(rep(1,M))))


lam <- 10^(seq(-7,-4,by=0.01))
GCV.value <- GCV.fit(lam, y=train_res-mean(train_res), x=Z.train[,act.set==1],w = w, R0=R, M=M)
#windows()
plot(x=lam, y=GCV.value,type="l")

opt.lam <- lam[which.min(GCV.value)]
# 1.174898e-05
chat <- est.smooth(opt.lam, y=train_res-mean(train_res), x=Z.train[,act.set==1], w=w, R0=R, M=M)
ypred <- Z.test[,act.set==1]%*%chat + mean(train_res)
MSPE.ssfam <- mean((ypred-test_res)^2)
#1.487567e-08



set <- which(abs(fit$beta)>0)
#figure 3
# windows()
# pdf("estimated components of pm2.5 analysis.pdf", width=9, height=9)
par(mfrow=c(2,2), font.axis=2, oma = c(5,4,0,0) + 2, mar = c(0,0,1,1) + 2)
for(j in 1:P)
{ #j = 11
  varin.mv <- ((j-1)*M+1):(j*M)
  st <- sum(act.set[1:((j-1)*M+1)]); st
  en <- sum(act.set[1:(j*M)]); en
  index <- intersect(varin.mv, set) 
  #len = len + length(index)
  if(length(index)==0){
    next
  }
  
  plot(x=zeta.train[,j][order(zeta.train[,j])], y = (Z.train[,varin.mv]%*%chat[st:en])[order(zeta.train[,j])], lwd=2, type="l",
       xlab=" ", ylab=" ")
 # points(x=zeta.train[,j][order(zeta.train[,j])], y=(Z.train[,varin.mv]%*%fit$beta[varin.mv])[order(zeta.train[,j])],lwd=2, type="l", lty=1)
  mtext(bquote(bold(hat(f))[.(j)]), side=3, line=.5, cex.lab=1.3,cex.axis=1)
  
}
# dev.off()



#functional linear regression with lasso
#dim(trainPCscore); dim(testPCscore)
flr <- cv.glmnet(x=trainPCscore, y=train_res)
#plot(flr)
fit <- glmnet(x=trainPCscore, y=train_res, lambda=flr$lambda.min)
# fit$beta
ypred <- predict(fit, testPCscore)
MSPE.flr <- mean((ypred - test_res)^2)
#  1.676306e-08


