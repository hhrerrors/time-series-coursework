largep<-function(Y,d,k){
  p<-nrow(Y)
  n<-ncol(Y) # p,n分别为Y的行和列
  sig1<-Y[,1:(n-1)]%*%t(Y[,2:n])/n     # Sigma_1
  sig0<-Y[,1:(n-1)]%*%t(Y[,1:(n-1)])/n # Sigma_0
  delta<-numeric(p)
  Ahhat<-matrix(0,nrow = p,ncol = p) #p阶零矩阵,A的估计A_hat
  Bhhat<-matrix(0,nrow = p,ncol = p)
  for(j in 1:p){
    z <- sig1[,j] #z_j=(Sigma_1)'*e_i
    if(k==0){
      tau1 <- 0
    }
    else{
      if(j==1){
        tau1 <- 2:min(1+k,p)
      }
      else if(j==p){
        tau1 <- max(p-k,1):(p-1)
      }
      else{
        tau1 <- c(max(j-k,1):(j-1),(j+1):min(j+k,p))
      }
    }
    tau2 <- c(max(j-k,1):min(j+k,p))
    tau <- length(tau1)+length(tau2)  # 式(2.6)
    u <- matrix(nrow = tau, ncol = n-1)
    for(t in 2:n){
      u[,(t-1)] <- c(Y[tau1,t],Y[tau2,t-1]) # 得出u_{t,j}
    }
    for(l in 1:p){
      de <- vector(length=n-1)
      for(t in 2:n){
        de[t-1] <- abs(Y[l,t-1])*sum(abs(Y[tau1,t]))+abs(Y[l,t-1])*sum(abs(Y[tau2,t-1]))
      }
      delta[l] <- sum(de)/n  #得到delta_l^j
    }
    del<-order(delta,decreasing = T)[1:d] # 得到前d_j个最大的delta_l^j的角标构成的向量
    w<-matrix(0,nrow = length(del),ncol = n)
    for(t in 1:(n-1)){
      w[,t] <- Y[del,t]  # 得到w_{t}^j
    }
    W<-matrix(0,nrow = d,ncol = tau)
    z<-matrix(0,nrow=d,ncol=1)
    for(t in 2:n){
      W = W + w[,t-1] %*% t(u[,(t-1)])
      z = z + w[,t-1] * Y[j,t]
    }
    W = (1/n)*W   # W_j的估计
    z = (1/n)*z   # z_j的估计
    bbb <- t(W)%*%W
    bb <- solve(bbb,tol=1e-21)
    beta<-bb %*% t(W)%*%z
    Ahhat[tau1,j]<-beta[1:length(tau1)]    # A[,j]的估计
    Bhhat[tau2,j]<-beta[(length(tau1)+1):length(beta)]
  }# B[,j]的估计
  return(list(Ahhat,Bhhat))
}

