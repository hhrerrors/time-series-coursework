sig <- function(r,Y){
  p <- nrow(Y)
  n <- ncol(Y)
  x <- matrix(0,nrow=r*p,ncol=p)
  for (i in 1:r){
    n1 <- (i-1)*p+1
    n2 <- i*p
    x[n1:n2,] <- t((Y[,(i+1):n]) %*% t(Y[,1:(n-i)]))/n
  }
  return(x)
}

mult <- function(r,Y,k){
  p <- nrow(Y)
  n <- ncol(Y)
  Ahhat<-matrix(0,nrow = p,ncol = p) #p阶零矩阵,A的估计A_hat
  Bhhat<-matrix(0,nrow = p,ncol = p)
  for(j in 1:p){
    xj <- sig(r,Y)[,j]
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
    tau2 <- c(max(j-k,1):min(j+k,p))
    tau <- length(tau1)+length(tau2)  # 式(2.6)
    u <- matrix(0,nrow = tau, ncol = n-1)
    for(m in 2:n){
      u[,(m-1)]<-c(Y[tau1,m],Y[tau2,(m-1)]) # 得出u_{t,j}
    }
    Gj <- matrix(0,nrow=r*p,ncol=tau)
    for (i in 1:r){
      n1 <- (i-1)*p+1
      n2 <- i*p
      n3 <- n-i
      n4 <- i
      Gj[n1:n2,] <- Y[,1:n3]%*%t(u[,n4:(n-1)])/n
    }
    if(k==0){
      Bhhat[tau2,j]<-beta #若k=0,即带宽为0，又因为A的对角元为0,A为零矩阵
    }
    else{
    beta <- solve(t(Gj)%*%Gj,tol=1e-21) %*% t(Gj) %*% xj
    Ahhat[tau1,j]<-beta[1:length(tau1)]    #A[,j]的估计
    Bhhat[tau2,j]<-beta[(length(tau1)+1):tau]  #B[,j]的估计
    }
  }
  return(list(Ahhat,Bhhat))
}

