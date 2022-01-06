SampleGenerate1 <-function(k0,p,n,r=1){
  A <- matrix(0,nrow=p,ncol=p)
  B <- matrix(0,nrow=p,ncol=p)
  epsilon <- matrix(0,nrow=p,ncol=n)
  for(i in 1:p){
    for(j in 1:p){
      if(abs(i-j)==k0){
        A[i,j] = sample(c(-2,2),1,prob = c(0.5,0.5)) ## uniform distribution on {-2,2}
        B[i,j] = sample(c(-2,2),1,prob = c(0.5,0.5)) 
      }
      if((abs(i-j)>0) && (abs(i-j)<k0)){
        A[i,j] = rnorm(1)*(1-sample(c(0,1),1,prob=c(0.6,0.4)))
        B[i,j] = rnorm(1)*(1-sample(c(0,1),1,prob=c(0.6,0.4)))
      }
      if(abs(i-j)==0){
        w = sample(c(0,1),1,prob=c(0.6,0.4))
        A[i,j] = 0
        B[i,j] = w+rnorm(1)*(1-w)
      }
      if(abs(i-j)>k0){
        A[i,j] = 0
        B[i,j] = 0
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
    Y0[,t] = solve(diag(p)-A) %*% B %*% Y0[,t-1] + solve(diag(p)-A) %*% epsilon[,t-1]
  }
  Y <- Y0[,2:(n+1)]
  return(list(Y,A,B,epsilon))
}
```