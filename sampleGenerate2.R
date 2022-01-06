SampleGenerate2<-function(k0,p,n,r=1){
  A <- matrix(0,nrow=p,ncol=p)
  B <- matrix(0,nrow=p,ncol=p)
  epsilon <- matrix(0,nrow=p,ncol=n)
  for(i in 1:p){
    for(j in 1:p){
      if(abs(i-j)==k0){
        a<-runif(1,-2.5,2.5)
        b<-runif(1,-2.5,2.5)
        while(a<1.5 & a>-1.5){
          a<-runif(1,-2.5,2.5)
        }
        A[i,j] =a
        while(b<1.5 & b>-1.5){
          b<-runif(1,-2.5,2.5)
        }
        B[i,j] =b  # uniform distribution on [-2.5,-1.5]∪[1.5,2.5]
      }
      if(abs(i-j)>0 && abs(i-j)<k0){
        A[i,j] = runif(1,-1,1)
        B[i,j] = runif(1,-1,1)
      }
      if(abs(i-j)==0){
        B[i,j] = runif(1,-1,1)
      }
    }
  }
  n0 = runif(2,0.4,0.8)
  n1 = n0[1]
  n2 = n0[2]
  A = n1*A/norm(A,type="2")
  B = n2*B/norm(B,type="2")
  Y0 = matrix(0,nrow=p,ncol=(n+1))
  for(t in 2:(n+1)){
    epsilon[,t-1] = rnorm(p,0,1)
    Y0[,t] = solve(diag(p)-A) %*%B%*% Y0[,t-1] + solve(diag(p)-A) %*% epsilon[,t-1]
  }
  Y = Y0[,2:(n+1)]
  return(list(Y,A,B,epsilon))
}
```


```{r}
parameterEstimate<-function(Y,k){
  p<-nrow(Y)
  n<-ncol(Y) # p,n分别为Y的行和列
  sig1<-Y[,1:(n-1)]%*%t(Y[,2:n])/n      # Sigma_1
  sig0<-Y[,1:(n-1)]%*%t(Y[,1:(n-1)])/n # Sigma_0
  rss<-numeric(p) #p维零向量,RS                                    S
  nss<-numeric(p)
  Ahat<-matrix(0,nrow = p,ncol = p) #p阶零矩阵,A的估计A_hat
  Bhat<-matrix(0,nrow = p,ncol = p) #B的估计B_hat
  for(j in 1:p){
    z<-sig1[,j] #z_j=(Sigma_1)'*e_j
    if(k==0){
      tau1<-0
    }
    else{
      if(j==1){
        tau1 <- c(2:(min(j+k,p)))
      }
      else if(j==p){
        tau1 <- c(max(j-k,1):(j-1))
      }
      else{
        tau1<-c(max(j-k,1):(j-1),(j+1):min(j+k,p))
      }
    }
    tau2<-c(max(j-k,1):min(j+k,p))
    tau<-length(tau1)+length(tau2)  # 式(2.6)
    v<-cbind(sig1[,tau1],sig0[,tau2])
    bb <- solve(t(v)%*%v,tol=1e-21)
    beta<-bb%*%t(v)%*%z #式(2.9),最小二乘法算beta,beta是tau维向量
    y<-rbind(Y[tau1,2:n],Y[tau2,1:(n-1)]) #y为tau*(n-1)矩阵
    resd=Y[j,2:n]-t(beta)%*%y #计算残差向量
    rss[j]<-sum(resd^2) #残差平方和
    nss[j]<-length(beta) #即tau
    if(k==0){
      Bhat[tau2,j]<-beta #若k=0,即带宽为0，又因为A的对角元为0,A为零矩阵
    }
    else{
      Ahat[tau1,j]=beta[1:length(tau1)]
      Bhat[tau2,j]=beta[(length(tau1)+1):nss[j]]
    }
  }
  return(list(rss[1:p],nss[1:p],Ahat[1:p,1:p],Bhat[1:p,1:p]))
}