library(MASS)

# Matern(5/2)
cov_m52 = function(b, tau1_sq, tau2_sq, d){
  f1 = tau1_sq * (1 + sqrt(5) * d / b + 5 * d*d / (3*b*b)) * exp(-sqrt(5) * d / b)
  f2 = tau2_sq * (d == 0)
  
  return(f1 + f2)
}

#squared exponential function
cov_se = function(b, tau1_sq, tau2_sq, d){
  f1 = tau1_sq * exp( - 0.5 * d*d / (b*b) )
  f2 = tau2_sq * (d == 0)
  
  return(f1 + f2)
}

gp_sample = function(X, b, tau1_sq, tau2_sq, cov_fun = 'm52'){
  n = length(X)
  
  mu = vector("numeric", length = n)
  C = matrix(0, nrow = n, ncol = n)
  for (i in 1:n){
    for (j in 1:n){
      if (cov_fun == 'm52')
        C[i, j] = cov_m52(b, tau1_sq, tau2_sq, abs(X[i] - X[j]))
      else if (cov_fun == 'se')
        C[i, j] = cov_se(b, tau1_sq, tau2_sq, abs(X[i] - X[j]))
    }
  }

  y = mvrnorm(mu = mu, Sigma = C)
  return(y)
}

run = function(){
  X = seq(0, 1, 0.01)
  
  
  y = gp_sample(X, b = 1, tau1_sq = 1, tau2_sq = 1e-6, cov_fun = 'm52'); 
  plot(X, y, ylim = c(-2, 2), "l")
  plot(X, y, "l")
}