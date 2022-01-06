Bandwidth <- function(K,Y,p,n,c){
  ## calculate k0
  ## Input : 
  ## K : integer,upper bound of k0
  ## Y ：p*p matrix,sample
  ## p : dimension of sample
  ## n : numbers of sample
  ## c : constant
  ## Output:
  ## k0 : integer,bandwidth
  rss <- matrix(nrow = p,ncol = K)
  for (k in 1:K){
    rss[,k] <- parameterEstimate(Y,k)[[1]]
  }
  w <- c/n
  result <- apply(rss,1,function(x){
    y <- vector(length=(K-1))
    for (i in 2:K){
      y[i-1] <- (x[i-1]+w)/(x[i]+w)
    }
    which.max(y)
  })
  return(max(result))
}