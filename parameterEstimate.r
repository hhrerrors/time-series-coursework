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