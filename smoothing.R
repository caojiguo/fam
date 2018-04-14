#function using GCV to choose tuning parameter in smoothing spline method
GCV.fit <- function(lambda, y, x, w, R0, M)
{ #lambda= 0.1; y=y0; x=Z.train[,which(abs(fit$beta)>0)]
  GCV <- numeric(length(lambda))
  basismat2 = t(x)%*%x;
  n <- length(y)
  Rmat <- matrix(0, nrow=nrow(basismat2), ncol=ncol(basismat2))
  for(j in 1:(ncol(x)/M))
  {
    index <- ((j-1)*M+1):(j*M)
    Rmat[index,index] <- R*w[j]
  }
  for(i in 1:length(lambda))
  {
    Bmat  = basismat2 + lambda[i]*Rmat;
    chat = ginv(Bmat)%*%t(x)%*%y;
    Mmat = ginv(Bmat)%*%t(x)
    Smat = x%*%Mmat
    df = sum(diag(Smat))
    yhat = x%*%chat;
    SSE = t(yhat-y)%*%(yhat-y)
    GCV[i] <- n*SSE/(n-df)^2
  }
  GCV
}

#function used to estimate the coefficient in spline estimation 
est.smooth <- function(lambda, y, x, w, R0, M)
{ #lambda=opt.lam; y=y0[train]-mean(y0[train]); x= x=Z.test[,which(abs(fit$beta)>0)]
  basismat2 = t(x)%*%x;
  n <- length(y)
  Rmat <- matrix(0, nrow=nrow(basismat2), ncol=ncol(basismat2))
  for(j in 1:(ncol(x)/M))
  {
    index <- ((j-1)*M+1):(j*M)
    Rmat[index,index] <- R*w[j]
  }
  Bmat  = basismat2 + lambda*Rmat;
  chat = ginv(Bmat)%*%t(x)%*%y;
  chat
}
