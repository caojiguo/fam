
source("smoothing.R")


#to implement the method, you need to install the following R packages, an example:
#install.packages(c("fda", "MASS", "cosso", "caTools", "mgcv", y


library(fda)
library(fdapace)
library(caTools)
library(earth)
library(mgcv)
library(cosso)
library(gglasso)
library(zoo)
library(splines)
library(MASS)




cat("generate functional data and response\n")
#200 equally-spaced points between [0,10]
m <- 200
t <- seq(0,10,length.out=m)

#training size
n.train <- 500
#n.train <- 100
#test size
n.test <- 200
#total sample size
N <- n.train+n.test

mean.fun <- t + sin(t) #mean function
a <- 31.5
b <- 0.6
K <- 20  #20 FPCs
eigenv <- a*(b^c(1:K)) #eigen values
fourb <- create.fourier.basis(rangeval=range(t), nbasis=K+1) #the first fourier basis is constant
basismat <- eval.basis(t, fourb)
B = eval.basis(t, fourb)[,-1] #eigen functions
sim.dat <- matrix(0, nrow=N, ncol=m)
eigens <- zeta <- matrix(0,nrow=N,ncol=K) #FPC score and uniform transformation of FPC score
Run <- 1000 #100 simulations

M <- 20 #choose 20  cubic  B-spline basis functions in our method
Bbasis <- create.bspline.basis(c(0, 1), norder=4, nbasis=M)
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

MSPE.mars <- MSPE.fam1 <- MSPE.fams <- MSPE.cosso <- MSPE.grp <- MSPE.ssfam <- numeric(Run)
time.mars <- time.fam1 <- time.fams <- time.cosso <- time.grp <- time.ssfam <- numeric(Run)
#model size for fams, cosso, ssfam
TP <- matrix(NA, nrow=Run, ncol=5) #True positive
FP <- matrix(NA, nrow=Run, ncol=5)

#fam1, fams, cosso, grp and ssfam for f_1, f_2 and f_4
ISE.fam1 <- ISE.fams <- ISE.cosso <- ISE.grp <- ISE.ssfam <- matrix(NA,ncol=3, nrow=Run)

#L <- 5
set.seed(20161104)
for(r in 1:Run)
{ #r = 1
  cat("This is Iteration: ", r, "\n")
  for(k in 1:K)
  {
    eigens[,k] <- rnorm(N,sd=sqrt(eigenv[k]))
    zeta[,k] <- pnorm(eigens[,k], sd=sqrt(eigenv[k]))
  }
  
  #observation of functional predictor: X(t_{ij}) + \eps_{ij}, eps: measurement errors
  for(i in 1:N)
    sim.dat[i,] <- mean.fun + B %*% eigens[i,] + rnorm(m, sd=sqrt(0.1))
  y0 <- 1.2 + zeta[,1]*exp(zeta[,1])-1+cos(2*pi*zeta[,2])+3*(zeta[,4]-1/4)^2-7/16+rnorm(N,0,sd=0.82)  
  
  # compute signal-to noise ratio
  #f1 <- function(x) (x*exp(x)-1)^2
  #f2 <- function(x) (cos(2*pi*x))^2
  #f3 <- function(x) (3*(x-1/4)^2-7/16)^2
  #integrate(f1,0,1)$value+integrate(f2,0,1)$value+integrate(f3,0,1)$value
  
  #SNR: 2
  #sqrt(1.3348*0.5)
  #0.82
  #index <- sample(x=rep(0:1, times=c(n.test, n.train)),size=N)
  #repeat the prodecure for N = 400, index <- sample(x=rep(0:1, times=c(0.5*N, 0.5*N)),size=N)
  #to get the result for n = 200. 
 
  # train <- which( index == 1 | index == 2 | index == 3)
  # train <- which(index ==1)
  # test <- (1:N)[-train]
  # train.sim.dat <- sim.dat[1:n.train,] #training data
  # test.sim.dat <- sim.dat[1:n.train+n.test,]  #test data
  #===========================================
  #PACE on training data
  cat("FPCA on simulated data\n")
  sim.yList <- list()
  sim.tList <- list()
  for(i in 1:N)
  {
    sim.yList[[i]] <-  sim.dat[i,]
    sim.tList[[i]] <-  t
  }
  
  #PACE from Fang Yao
  res <- FPCA(sim.yList[1:n.train], sim.tList[1:n.train], list(dataType='Dense',kernel='epan', verbose=TRUE, maxK=20))
  
  P <- min(which(res$cumFVE > 99.9)) #choose # of FPCs, explain over 99.9% variability
  #P <- min(P, 16)
  eigens.train <- eigens[1:n.train,] #true FPC scores for training data
  eigv.est <- res$lambda[1:P] #estimated eigenvalues
  score.train <- res$xiEst[,1:P] #estimated FPC scores, pay attention to sign
  
  zeta.train <- matrix(0, nrow=n.train, ncol=P)
  flag <- numeric(P)
  for(j in 1:P)
  { #j = 2
    index <- which.max(abs(eigens.train[,j]))
    flag[j] <- sign(eigens.train[index,j])*sign(score.train[index,j])
    score.train[,j] <- flag[j]*score.train[,j] #update the sign of estimated FPC scores
    zeta.train[,j] <- pnorm(score.train[,j], sd=sqrt(eigv.est[j])) #estimated transformation for traning data
  }
  
  #head(score.train[,1:4]); head(eigens.train[,1:4])
  
  # mean.train <- res$mu
  # cen.test <- test.sim.dat - matrix(rep(mean.train,each=N-n),nrow=N-n)
  # score.test <- matrix(0,nrow=N-n,ncol=P)
  # #estimate FPC scores for test data using information from training data
  # for(j in 1:P)
  # {
  #   index <- which.max(abs(B[,j]))
  #   flag <- sign(B[index,j])*sign(res$phi[index,j])
  #   c.phi <- flag*res$phi[,j]
  #   for(i in 1:N-n)
  #     #j = 1; i = 1
  #     score.test[i,j] <- trapz(t,cen.test[i,]*(c.phi))
  # }
  
  score.test= predict(res, newLy=sim.yList[(n.train+1):N], newLt=sim.tList[(n.train+1):N], sigma2=res$sigma2, K=P)  
  score.test = score.test %*% diag(flag)
  
  #head(score.test[,1:4]); head(eigens[(n.train+1):N,1:4])
  
  
  zeta.test <- matrix(0,nrow=nrow(score.test),ncol=P)
  for(k in 1:P)
    zeta.test[,k] <- pnorm(score.test[,k], sd=sqrt(eigv.est[k])) #estimated transformation for test data
  
  #head(zeta[1:n.train,1:4]);head(zeta.train[,1:4]); head(zeta[1:n.test+n.train,1:4]); head(zeta.test[,1:4])
  
  train.vec <- 1:n.train
  test.vec <- (n.train+1):N
  
  #==============================================
  #compare different fitting algorithms, as shown in Fang yao, including 
  #method proposed by Fangyao and our method (SSFAM)
  cat("MARS\n")
  ptm <- proc.time()
  marsModel <- earth(y=y0[train.vec], x= zeta.train) # build model
  time.mars[r] <- (proc.time()-ptm)[1]
  ev <- evimp (marsModel)[,1] # estimate variable importance
  index <- as.numeric(ev)
  TP[r, 1] <- length(intersect(index, c(1,2,4)))/3
  FP[r, 1] <- length(setdiff(index, c(1,2,4)))/(P-3)
  ypred <- predict(marsModel, newdata=zeta.test)
  MSPE.mars[r] <- mean((ypred - y0[test.vec])^2)
  
  
  #========================================
  cat("FAM_1 with known compoent structure and estimated zeta\n")
  est.predictor.train <- matrix(0, nrow=n.train,ncol=3)
  est.predictor.train[,1] <- zeta.train[,1]*exp(zeta.train[,1])- 1
  est.predictor.train[,2] <- cos(2*pi*zeta.train[,2])
  est.predictor.train[,3] <- 3*(zeta.train[,4]-1/4)^2-7/16
  est.predictor.test <- matrix(0, nrow=n.test,ncol=3) #three true functions
  est.predictor.test[,1] <- zeta.test[,1]*exp(zeta.test[,1])- 1
  est.predictor.test[,2] <- cos(2*pi*zeta.test[,2])
  est.predictor.test[,3] <- 3*(zeta.test[,4]-1/4)^2-7/16
  train.data <- data.frame(cbind(y0[train.vec], est.predictor.train))
  colnames(train.data) <- c("y","pred1", "pred2", "pred4")
  test.data <- data.frame(cbind(y0[test.vec], est.predictor.test))
  colnames(test.data) <- c("y","pred1", "pred2", "pred4")
  
  ptm <- proc.time()
  fit.fam1 <- lm(y ~ ., data=train.data)
  time.fam1[r] <- (proc.time()-ptm)[1]
  
  ypred <- predict(fit.fam1, test.data)
  MSPE.fam1[r] <- mean((ypred - y0[test.vec])^2)
  ISE.fam1[r,] <- as.numeric((coef(fit.fam1)[2:4]-1)^2)
  #0.6775216
  TP[r, 2] <- length((summary(fit.fam1)$coefficients[-1,4] < 0.05))/3
  FP[r, 2] <- 0
  #=============================================
  cat("saturated FAM\n") #fit an additive model using the first P FPC scores
  full.data.train <- data.frame(cbind(y0[train.vec], zeta.train))
  full.data.test <- data.frame(cbind(y0[test.vec], zeta.test))
  colnames(full.data.train) <- c("y",paste("x", 1:P, sep=""))
  colnames(full.data.test) <- c("y",paste("x", 1:P, sep=""))
  predic1 <- paste( "s(", paste("x", 1:P, sep=""),")", sep="",collapse=" + ")
 #for small training set, choose the following predict1.
 # predic1 <- paste( "s(", paste(paste("x", 1:P, sep=""), paste("k", "=4", sep=""), sep=","), ")", sep="",collapse=" + ")
  
  #as.formula(paste("y ~", predic2))
  ptm <- proc.time()
  fit.fams <- gam(as.formula(paste("y ~", predic1)), data=full.data.train)
  time.fams[r] <- (proc.time()-ptm)[1]
  #freq.fams[which(summary(fit.fams)$s.pv < 0.05)] <- freq.fams[which(summary(fit.fams)$s.pv < 0.05)] + 1
  TP[r,3] <- length(intersect(which(summary(fit.fams)$s.pv < 0.05), c(1,2,4)))/3
  FP[r,3] <- length(setdiff(which(summary(fit.fams)$s.pv < 0.05), c(1,2,4)))/(P-3)
  ypred <- predict(fit.fams, full.data.test)
  MSPE.fams[r] <- mean((ypred - y0[test.vec])^2)
  #0.6788667
  
  results.fams <- predict(fit.fams, type="terms")
  arg1 <- zeta.train[,1][order(zeta.train[,1])]
  predictor1.est <- results.fams[,1][order(zeta.train[,1])]
  predictor1.true <- arg1*exp(arg1)-1
  ise.1 <- trapz(x=arg1, y=(predictor1.est-predictor1.true)^2)
  #0.078
  arg2 <- zeta.train[,2][order(zeta.train[,2])]
  predictor2.est <- results.fams[,2][order(zeta.train[,2])]
  predictor2.true <- cos(2*pi*arg2)
  ise.2 <- trapz(x=arg2, y=(predictor2.est-predictor2.true)^2)
  arg4 <- zeta.train[,4][order(zeta.train[,4])]
  predictor4.est <- results.fams[,4][order(zeta.train[,4])]
  predictor4.true <- 3*(arg4-1/4)^2-7/16
  ise.4 <- trapz(x=arg4, y=(predictor4.est-predictor4.true)^2)
  ISE.fams[r,] <- c(ise.1,ise.2,ise.4)
  
  #==============================================
  cat("CSEFAM\n") #from Hongxiao Zhu et al. JRSSB
  ptm <- proc.time()
  fit.cosso <- cosso(x=zeta.train, y=y0[train.vec],family="Gaussian", nbasis=n.train)
  tune <- tune.cosso(fit.cosso, plot.it=F)
  time.cosso[r] <- (proc.time()-ptm)[1]
  
  ypred <- predict.cosso(fit.cosso, xnew=zeta.test, M=tune$OptM, type="fit")
  MSPE.cosso[r] <- mean((ypred - y0[test.vec])^2)
  
  #comp.pred <- predict.cosso(fit.cosso, xnew=zeta.test, M=tune$OptM, type="nonzero")
  fit.sig <- predict.cosso(fit.cosso,M=tune$OptM,type="nonzero")
  TP[r,4] <- length(intersect(fit.sig, c(1,2,4)))/3
  FP[r,4] <- length(setdiff(fit.sig, c(1,2,4)))/(P-3)
  #freq.cosso[fit.sig] <- freq.cosso[fit.sig] + 1
  fit.coef <- predict.cosso(fit.cosso,M=tune$OptM,type="coefficients")
  theta1 <- c(fit.coef$theta[1],rep(0,P-1))
  theta2 <- c(0, fit.coef$theta[2], rep(0, P-2))
  theta4 <- c(rep(0,3),fit.coef$theta[4],rep(0,P-4))
  fitObj=twostep.Gaussian(fit.cosso$Kmat[,fit.cosso$basis.id,],fit.cosso$Kmat[fit.cosso$basis.id,fit.cosso$basis.id,],
                          fit.cosso$y, fit.cosso$wt, fit.cosso$tune$OptLam, tune$OptM)
  
  predictor.1 <- as.numeric(wsGram(bigGram(zeta.train, fit.cosso$x[fit.cosso$basis.id,]), theta1/(fit.cosso$wt^2))%*%fitObj$coefs)
  predictor.2 <- as.numeric(wsGram(bigGram(zeta.train, fit.cosso$x[fit.cosso$basis.id,]), theta2/(fit.cosso$wt^2))%*%fitObj$coefs)
  predictor.4 <- as.numeric(wsGram(bigGram(zeta.train, fit.cosso$x[fit.cosso$basis.id,]), theta4/(fit.cosso$wt^2))%*%fitObj$coefs)
# arg1 <- zeta.train[,1][order(zeta.train[,1])]
  predictor1.est <- predictor.1[order(zeta.train[,1])]
 # predictor1.true <- arg1*exp(arg1)-1
  ise.1 <- trapz(x=arg1, y=(predictor1.est-predictor1.true)^2)
 #  arg2 <- zeta.train[,2][order(zeta.train[,2])]
  predictor2.est <- predictor.2[order(zeta.train[,2])]
 # predictor2.true <- cos(2*pi*arg2)
  ise.2 <- trapz(x=arg2, y=(predictor2.est-predictor2.true)^2)
 # arg4 <-  zeta.train[,4][order(zeta.train[,4])]
  predictor4.est <- predictor.4[order(zeta.train[,4])]
 # predictor4.true <- 3*(arg4-1/4)^2-7/16
  ise.4 <- trapz(x=arg4, y=(predictor4.est-predictor4.true)^2)
  
  ISE.cosso[r,] <- c(ise.1, ise.2, ise.4)
  
  
  #===========================================
  cat("SSFAM\n") #our proposed method
  #step 1 group lasso using cubic B-spline
  Z.train <- matrix(0, nrow=n.train, ncol=P*M) #define Z matrix for training data
  Z.test <- matrix(0, nrow=n.test, ncol=P*M) #for test data
  for(j in 1:P)
  { #j = 1
    index <- ((j-1)*M+1):(j*M)
    #a <- eval.basis(zeta.est[,j], Bbasis); dim(a); dim(Z[,index])
    Z.train[,index] <- eval.basis(zeta.train[,j], Bbasis)
    Z.test[,index] <- eval.basis(zeta.test[,j], Bbasis)
  }
  
  mean.vec <- apply(Z.train, 2, mean); #for centralization
  Z.train <- Z.train - matrix(rep(mean.vec, times=n.train), nrow=n.train, byrow=T)
  mean.vec <- apply(Z.test, 2, mean); #for centralization
  Z.test <- Z.test - matrix(rep(mean.vec, times=n.test), nrow=n.test, byrow=T)
  
  #fit a group lasso
  G <- rep(1:P, each=M)
  ptm <- proc.time()
  cv1 <- cv.gglasso(x=Z.train, y=y0[train.vec], group=G, loss="ls", pred.loss="L2") #choosing tuning parameter lambda_1
  #quartz()
  #plot(cv1)
  # cv$lambda.1se
  fit <- gglasso(x=Z.train, y=y0[train.vec], group=G, lambda=(cv1$lambda.min + cv1$lambda.1se)/2)
  # fit$b0
  # mean(y0[train])
  #==================
  #step 2 adaptive group lasso
  w <- numeric(P) #define weight for each group
  for(j in 1:P)
  {
    index <- (M*(j-1)+1):(M*j)
    w[j] <- abs(1/sqrt(sum((fit$beta[index])^2)))
  }
  w[w==Inf] <- 1e8
  
  cv2 <- cv.gglasso(x=Z.train, y=y0[train.vec], group=G, loss="ls", pred.loss="L2", pf=w) #choosing tuning parameter for adaptive group lasso, lambda_2
  #windows()
  #plot(cv2)
  # cv$lambda.1se
  #fit an adaptive group lasso with specific weight
  fit <- gglasso(x=Z.train, y=y0[train.vec], group=G, lambda= (cv2$lambda.min + cv2$lambda.1se)/2, pf=w)
  time.grp[r] <- (proc.time()-ptm)[1]
  # fit$b0
  # mean(y0[train])
  # which(abs(fit$beta)>0)
  ypred <- predict(fit,newx=Z.test,type="link")
  MSPE.grp[r] <- mean((ypred - y0[test.vec])^2)
  
  
  
  # arg1 <- zeta.train[,1][order(zeta.train[,1])]
  predictor.1 <- Z.train[,1:M]%*%fit$beta[1:M]
  predictor1.est <- predictor.1[order(zeta.train[,1])]
  # predictor1.true <- arg1*exp(arg1)-1
  ise.1 <- trapz(x=arg1, y=(predictor1.est-predictor1.true)^2)
  #  arg2 <- zeta.train[,2][order(zeta.train[,2])]
  predictor.2 <- Z.train[,(M+1):(2*M)]%*%fit$beta[(M+1):(2*M)]
  predictor2.est <- predictor.2[order(zeta.train[,2])]
  # predictor2.true <- cos(2*pi*arg2)
  ise.2 <- trapz(x=arg2, y=(predictor2.est-predictor2.true)^2)
  # arg4 <-  zeta.train[,4][order(zeta.train[,4])]
  predictor.4 <- Z.train[,(3*M+1):(4*M)]%*%fit$beta[(3*M+1):(4*M)]
  predictor4.est <- predictor.4[order(zeta.train[,4])]
  # predictor4.true <- 3*(arg4-1/4)^2-7/16
  ise.4 <- trapz(x=arg4, y=(predictor4.est-predictor4.true)^2)
  # 0.8006597
  ISE.grp[r,] <- c(ise.1, ise.2, ise.4)
  #============================
  #step3: penalize roughness, smoothing spline: second derivative matrix or difference matrix
  act.set <- numeric(P)
  w <- numeric(P) #define weight for each group
  for(j in 1:P)
  {
    index <- (M*(j-1)+1):(M*j)
    w[j] <- abs(1/sqrt(sum((fit$beta[index])^2)))
    act.set[j] <- sum((fit$beta[index])^2) > 0
  }
  w <- w[w!=Inf]
  TP[r,5] <- length(intersect(which(act.set==1), c(1,2,4)))/3
  FP[r,5] <- length(setdiff(which(act.set==1), c(1,2,4)))/(P-3)
  
  lam <- seq(0.001,0.3,0.001) #candidate of lambda_3, smoothing parameter
  ptm <- proc.time()
  GCV.value <- GCV.fit(lam, y=y0[train.vec]-mean(y0[train.vec]), x=Z.train[,which(abs(fit$beta)>0)], w=w, R0=R, M=M)
  #windows()
 #plot(x=lam, y=GCV.value,type="l")
  # lam[which.min(GCV.value)]
  opt.lam <- lam[which.min(GCV.value)] #choose the value based on GCV 
  chat <- est.smooth(opt.lam, y=y0[train.vec]-mean(y0[train.vec]), x=Z.train[,which(abs(fit$beta)>0)],w=w, R0=R, M=M)
  time.ssfam[r] <- (proc.time()-ptm)[1] + time.grp[r]
  ypred <- Z.test[,which(abs(fit$beta)>0)]%*%chat + mean(y0[train.vec])
  MSPE.ssfam[r] <- mean((ypred-y0[test.vec])^2)
  # 0.6806995
  #length(which(abs(fit$beta)>0))
 # size[r,3] <- length(which(abs(fit$beta)>0))/M
  #fit$beta[seq(1, P*M, by=M)]
 # freq.ssfam[which(fit$beta[seq(1, P*M, by=M)] != 0)] <- freq.ssfam[which(fit$beta[seq(1, P*M, by=M)] != 0)] + 1
  
  predictor.1 <- Z.train[,1:M]%*%chat[1:M]
  predictor1.est <- predictor.1[order(zeta.train[,1])]
  # predictor1.true <- arg1*exp(arg1)-1
  ise.1 <- trapz(x=arg1, y=(predictor1.est-predictor1.true)^2)
  #  arg2 <- zeta.train[,2][order(zeta.train[,2])]
  predictor.2 <- Z.train[,(M+1):(2*M)]%*%chat[(M+1):(2*M)]
  predictor2.est <- predictor.2[order(zeta.train[,2])]
  # predictor2.true <- cos(2*pi*arg2)
  ise.2 <- trapz(x=arg2, y=(predictor2.est-predictor2.true)^2)
  # arg4 <-  zeta.train[,4][order(zeta.train[,4])]
  predictor.4 <- Z.train[,(3*M+1):(4*M)]%*%chat[(2*M+1):(3*M)]
  predictor4.est <- predictor.4[order(zeta.train[,4])]
  # predictor4.true <- 3*(arg4-1/4)^2-7/16
  ise.4 <- trapz(x=arg4, y=(predictor4.est-predictor4.true)^2)
  # 0.8006597
  ISE.ssfam[r,] <- c(ise.1, ise.2, ise.4)
  
}

save(MSPE.mars,MSPE.fam1,MSPE.fams,MSPE.cosso,MSPE.grp,MSPE.ssfam,
     time.mars,time.fam1,time.fams,time.cosso,time.grp,time.ssfam,
     ISE.fam1,ISE.fams,ISE.cosso,ISE.grp,ISE.ssfam,TP,FP, file="result.RData")


#repeat the above procedure for n.train = 100






MSPE <- cbind(MSPE.mars, MSPE.fam1, MSPE.fams, MSPE.cosso, MSPE.grp, MSPE.ssfam)
colnames(MSPE) <- c("MARS","FAM", "S-FAM", "CSE-FAM", "AGL-FAM", "CSS-FAM")
index <- which(is.na(ISE.ssfam[,3]))
#MSPE.ssfam[index]
windows()
#pdf("prediction of each model.pdf")
boxplot(MSPE,use.cols=T)
#dev.off()
apply(MSPE,2,summary)
apply(MSPE, 2, sd)




#estimate of each component 
ISE1 <- cbind(ISE.fam1[,1], ISE.fams[,1], ISE.cosso[,1], ISE.grp[,1], ISE.ssfam[,1])
ISE2 <- cbind(ISE.fam1[,2], ISE.fams[,2], ISE.cosso[,2], ISE.grp[,2], ISE.ssfam[,2])
ISE4 <- cbind(ISE.fam1[,3], ISE.fams[,3], ISE.cosso[,3], ISE.grp[,3], ISE.ssfam[,3])

colnames(ISE1) <- colnames(ISE2) <- colnames(ISE4) <- c("FAM", "S-FAM", "CSE-FAM", "AGL-FAM", "CSS-FAM")
#windows()
#pdf("estimation error for each component.pdf", width=10)
par(mfrow=c(2,2))
boxplot(ISE1, use.cols=T)
boxplot(ISE2, use.cols=T)
boxplot(ISE4, use.cols=T)
#dev.off()
apply(ISE1, 2, summary); apply(ISE1, 2, sd)
apply(ISE2, 2, summary); apply(ISE2, 2, sd)
apply(ISE4, 2, summary); apply(ISE4, 2, sd)


#TP and FP

apply(TP, 2, mean)
apply(TP, 2, sd)

apply(FP, 2, mean)
apply(FP,2,sd)






plot1.cv.gglasso <- function(x, sign.lambda = 1, ...) {
  cvobj <- x
  xlab <- "log(Lambda)"
  if (sign.lambda < 0) 
    xlab <- paste("-", xlab, sep = "")
  plot.args <- list(x = sign.lambda * log(cvobj$lambda), y = cvobj$cvm, ylim = range(cvobj$cvupper, 
                cvobj$cvlo), xlab = xlab, ylab = cvobj$name, type = "n")
  new.args <- list(...)
  if (length(new.args)) 
    plot.args[names(new.args)] <- new.args
  #do.call("plot", plot.args)
  Hmisc::errbar(sign.lambda * log(cvobj$lambda), cvobj$cvm, cvobj$cvupper, cvobj$cvlo, cap = 0.01, 
             col = "darkgrey", xlab='', ylab='')
  points(sign.lambda * log(cvobj$lambda), cvobj$cvm, pch = 20, col = "red")
  #axis(side = 3, at = sign.lambda * log(cvobj$lambda), tick = FALSE, line = 0)
  abline(v = sign.lambda * log(cvobj$lambda.min), lty = 3)
  #abline(v = sign.lambda * log(cvobj$lambda.1se), lty = 3)
  invisible()
} 




#the following figure is for n = 300, sample size of training data


#windows()
#pdf("estimation from 1 simulation.pdf", width=9)
par(mfrow=c(2,3),mar=c(5,6, 6 ,2) +.1, cex.axis=1, cex.lab=1.6, font.lab=2,lwd=2,adj=0.5) 
plot1.cv.gglasso(cv1)
title(xlab=expression(paste("log(", lambda[1], ")")), ylab="Cross-validation Errors")
plot1.cv.gglasso(cv2)
title(xlab=expression(paste("log(", lambda[2], ")")), ylab="Cross-validation Errors")
#plot(cv2,xlab=expression(paste("log(", lambda[2], ")")), ylab="Cross-validation Errors")
plot(x=lam, y=GCV.value,type="l", xlab=expression(lambda[3]), ylab="GCV")
abline(v=lam[which.min(GCV.value)],lty=3)

curve(x*exp(x)-1, 0, 1, xlab=" ",ylab=expression(bold(hat(f)[1])),lwd=2)
lines(x=zeta.train[,1][order(zeta.train[,1])], y=(Z.train[,1:M]%*%chat[1:M])[order(zeta.train[,1])],lty=2)
lines(x=zeta.train[,1][order(zeta.train[,1])], y=(Z.train[,1:M]%*%fit$beta[1:M])[order(zeta.train[,1])],lty=3,,col="red")


curve(cos(2*pi*x), 0, 1, xlab=" ",ylab=expression(bold(hat(f)[2])),lwd=2,ylim=c(-1, 1.2))
lines(x=zeta.train[,2][order(zeta.train[,2])], y=(Z.train[,(M+1):(2*M)]%*%chat[(M+1):(2*M)])[order(zeta.train[,2])],lty=2)
lines(x=zeta.train[,2][order(zeta.train[,2])], y=(Z.train[,(M+1):(2*M)]%*%fit$beta[(M+1):(2*M)])[order(zeta.train[,2])],lty=3,col="red")


curve(3*(x-1/4)^2-7/16,xlab=" ", ylab=expression(bold(hat(f)[4])),lwd=2, ylim=c(-0.4,1))
lines(x=zeta.train[,4][order(zeta.train[,4])], y=(Z.train[,(3*M+1):(4*M)]%*%chat[(2*M+1):(3*M)])[order(zeta.train[,4])],lty=2)
lines(x=zeta.train[,4][order(zeta.train[,4])], y=(Z.train[,(3*M+1):(4*M)]%*%fit$beta[(3*M+1):(4*M)])[order(zeta.train[,4])],lty=3,col="red")
dev.off()


