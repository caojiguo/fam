
source("smoothing.R")

#you need to install some R packages first before running the following code
#an example:
#install.packages(c("fda", "MASS", "cosso", "caTools", "mgcv", "earth", 
#"R.matlab","gglasso","zoo","splines"))

#==================================
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
#library(R.matlab)
library(gglasso)
library(zoo)
library(splines)
#=================================
# set working directory
# setwd("/Users/cao/Dropbox/Research/Peijun/FAM/nfam2017Feb19/r code/tecator analysis")
tecator <- read.table("tecator.txt", sep="")
#head(tecator)
tecator <- unlist(t(tecator))
all_samples <- matrix(NA, nrow=240, ncol=125)
for (i in 1:240)
  all_samples[i,] <- tecator[((i-1)*125+1):(i*125)]
allPCs <- all_samples[,101:122]
allcontent <- all_samples[,123:125]
train <- all_samples[1:129,]
monit <- all_samples[130:172,]
test <- all_samples[173:215,]
E <- all_samples[216:240,]


set.seed(18265)
ID <- sample(x=rep(0:1, times=c(15,10)),size=25)
#table(ID)

#figure 2
#windows()
#pdf("Tecator profile.pdf")
par(pch="o", font.axis=2, font.lab=2, cex.axis=1.2, cex.lab=1.2)
matplot(x=seq(851,1050,by=2), y=t(all_samples[,1:100]),type="l", lty=1, col=1,
        xlab="Wavelengths (nm)", ylab=" ")
#dev.off()


#as suggested in the Tecator data set.
train_spec = rbind(all_samples[1:172,1:100], E[ID==0,1:100])
test_spec = rbind(all_samples[173:215,1:100], E[ID==1, 1:100])
#train_protein = c(all_samples[1:172,125], E[ID==0,125])
#test_protein = c(all_samples[173:215,125], E[ID==1, 125])


#use protein as the response
#===============================================
#FPCA using PACE
ntrain = dim(train_spec)[1]
#187
#numgrid = dim(train_spec)[2]
ntest = dim(test_spec)[1]
#53
nm <- seq(851,1050,by=2)


# FPCA using PACE in Matlab
# The PACE code is provided in the folder "matlab code for PACE"
# Use the PACE code in Matlab to obtain results for the train data set
# The PACE results on the train data set are stored in traintec.csv
# use FPCA function

train <- read.csv("traintec.csv",header=T)
#response in training data
train_res <-train[, 21]
#scaled FPC score for training data
zeta.train <- train[,1:20]
P <- ncol(zeta.train)


# Use the PACE code in Matlab to obtain results for the test data set
# The PACE results on the test data set are stored in test.mat
test <- read.csv("testtec.csv",header=T)
#response in test data
test_res <- test[, 21]
#scaled FPC score for test data
zeta.test <- test[,1:P]


set.seed(123456)
#=====================================================
cat("MARS\n")
#mars with 20 FPCs
marsModel <- earth(y=train_res, x= zeta.train) # build model
mars.pred <- predict(marsModel, newdata=zeta.test)
MSPE.mars <- mean((mars.pred - test_res)^2)
#0.99

#=================================================
cat("LAF\n") #functional linear regression with lasso
trainPCscore <- apply(zeta.train, 2, qnorm)
#all(is.finite(trainPCscore))
testPCscore <- apply(zeta.test, 2, qnorm)
#all(is.finite(testPCscore))
flr <- cv.glmnet(x=trainPCscore, y=train_res)
#plot(flr)
fit <- glmnet(x=trainPCscore, y=train_res, lambda=flr$lambda.min)
ypred <- predict(fit, testPCscore)
MSPE.laf <- mean((ypred - test_res)^2)
#0.66

#==================================
cat("saturated FAM\n")  #saturated FAM with 20 FPCs
full.data.train <- data.frame(cbind(train_res, zeta.train))
full.data.test <- data.frame(cbind(test_res, zeta.test))
colnames(full.data.train) <- c("protein",paste("x", 1:P, sep=""))
colnames(full.data.test) <- c("protein",paste("x", 1:P, sep=""))
#fam with 20 FPCs
#predic20 <- paste( "s(", paste(paste("x", 1:P, sep=""), paste("k", "=9", sep=""), sep=","), ")", sep="",collapse=" + ")
predic <- paste( "s(", paste("x", 1:P, sep=""),")", sep="",collapse=" + ")
fam <- gam(as.formula(paste("protein ~", predic)),data=full.data.train, method="REML")
ypred <- predict(fam, full.data.test)
MSPE.fam <- mean((ypred - test_res)^2)
# 0.56


#=================================
cat("CSE-FAM\n") #proposed by Fang
#using 20 FPC
G.Obj=cosso(x=zeta.train,y=train_res,family="Gaussian")
#windows()
tun1 <- tune.cosso(G.Obj,5,F)
#the optimal one tun1$OptM = 11.45
#fit.sig <- predict.cosso(G.Obj,M=M0,type="nonzero")
#M0 = 12
ypred <- predict.cosso(G.Obj,xnew=zeta.test,M=tun1$OptM,type="fit")
MSPE.cosso <- mean((test_res-ypred)^2)
#0.55


#===========================================
#our proposed approach without smoothing, only the first two steps
cat("CSS-FAM\n")
#step 1 group lasso using cubic B-spline
M <- 10
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

mean.vec <- apply(Z.train, 2, mean)
Z.train <- Z.train - matrix(rep(mean.vec, times=ntrain), nrow=ntrain, byrow=T)
mean.vec <- apply(Z.test, 2, mean)
Z.test <- Z.test - matrix(rep(mean.vec, times=ntest), nrow=ntest, byrow=T)

G <- rep(1:P, each=M)
cv <- cv.gglasso(x=Z.train, y=train_res, group=G, lambda=exp(seq(-8,-5,by=0.05)),
                 loss="ls",  pred.loss="L2")
#windows()
plot(cv)
fit <- gglasso(x=Z.train, y=train_res, group=G, lambda=cv$lambda.min)
#step 2 adaptive group lasso
w <- numeric(P)
for(j in 1:P)
{
  index <- (M*(j-1)+1):(M*j)
  w[j] <- abs(1/sqrt(sum((fit$beta[index])^2)))
}
w[w==Inf] <- 1e8

cv <- cv.gglasso(x=Z.train, y=train_res, group=G, lambda=exp(seq(-6,-2.5,by=0.05)),
                   loss="ls", pred.loss="L2", pf=w)
#windows()
#plot(cv)

fit <- gglasso(x=Z.train, y=train_res, group=G, lambda=cv$lambda.min, pf=w)
#which(abs(fit$beta)>0)
ypred <- predict(fit,newx=Z.test,type="link")
MSPE.grp <- mean((ypred - test_res)^2)
#0.92

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

lam <- 10^(seq(-5,-2,by=0.01))
GCV.value <- GCV.fit(lam, y=train_res-mean(train_res), x=Z.train[,act.set==1],w = w, R0=R, M=M)
#windows()
plot(x=lam, y=GCV.value,type="l")


opt.lam <- lam[which.min(GCV.value)]
chat <- est.smooth(opt.lam, y=train_res-mean(train_res), x=Z.train[,act.set==1], w=w, R0=R, M=M)
ypred <- Z.test[,act.set==1]%*%chat + mean(train_res)
MSPE.ssfam <- mean((ypred-test_res)^2)
#0.51



set <- which(abs(fit$beta)>0)
#figure 3
# windows()
# pdf("estimated components of tecator analysis.pdf", width=9, height=9)
par(mfrow=c(4,4), font.axis=2, oma = c(5,4,0,0) + 2, mar = c(0,0,1,1) + 2)
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
  
# plot(x=zeta.train[,j][order(zeta.train[,j])], y = (Z.train[,varin.mv]%*%chat[st:en])[order(zeta.train[,j])], lwd=2, type="l",
#      xlab=" ", ylab=" ")
plot(x=zeta.train[,j][order(zeta.train[,j])], y=(Z.train[,varin.mv]%*%fit$beta[varin.mv])[order(zeta.train[,j])],lwd=2, type="l", lty=1)
  mtext(bquote(bold(hat(f))[.(j)]), side=3, line=.5, cex.lab=1.3,cex.axis=1)
  
}
# dev.off()















































