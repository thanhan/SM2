library(MASS)

sim_mvn <- function(n, mu, Sigma){
  return (mvrnorm(n = n, mu = mu, Sigma = Sigma))
}

mu <- c(0,0)
Sigma <- matrix(c(10,3,3,2),2,2)

mle_mvn <- function(X){
  mu = colMeans(X)
  n = nrow(X)
  Sigma = cov(X)*(n-1)/n
  return(list(mu, Sigma))
}

boot_mle <- function(bn, X){
  m = ncol(X)
  boot_mu = matrix(0, ncol = m, nrow = bn)
  boot_cov = matrix(0, ncol = m*m, nrow = bn)
  for(i in 1:bn)
  {
    idx = sample(1:nrow(X), nrow(X), replace=TRUE)
    bx = X[idx,]
    temp = mle_mvn(bx)
    boot_mu[i,] = temp[[1]]
    boot_cov[i,] = temp[[2]]
  }
  return (list(boot_mu, boot_cov))
}