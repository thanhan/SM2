library(MASS)
library(plotrix)



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

compute_c = function(X, b, tau1_sq, tau2_sq, cov_fun){
  n = length(X)
  
  C = matrix(0, nrow = n, ncol = n)
  for (i in 1:n){
    for (j in 1:n){
      if (cov_fun == 'm52')
        C[i, j] = cov_m52(b, tau1_sq, tau2_sq, abs(X[i] - X[j]))
      else if (cov_fun == 'se')
        C[i, j] = cov_se(b, tau1_sq, tau2_sq, abs(X[i] - X[j]))
    }
  }
  
  return(C)
}

gp_sample = function(X, b, tau1_sq, tau2_sq, cov_fun = 'm52'){
  n = length(X)
  
  mu = vector("numeric", length = n)
  C = compute_c(X, b, tau1_sq, tau2_sq, cov_fun)

  y = mvrnorm(mu = mu, Sigma = C)
  return(y)
}

# sample GP function and plot
run1 = function(){
  X = seq(0, 1, 0.01)
  
  y = gp_sample(X, b = 1, tau1_sq = 1, tau2_sq = 1e-6, cov_fun = 'm52'); 
  plot(X, y, ylim = c(-2, 2), "l")
  plot(X, y, "l")
}


gp_predict = function(X, y, x_star, sigma_sq = 1, b = 0.001, tau1_sq = 2, tau2_sq = 1e-6, cov_fun = 'm52'){
  n1 = length(x_star)
  n2 = length(X)
  C = compute_c(c(x_star, X), b, tau1_sq, tau2_sq, cov_fun)
  C_ss = C[1:n1,1:n1]
  C_sx = C[1:n1, (n1+1):(n1+n2)]
  C_xx = C[(n1+1):(n1+n2), (n1+1):(n1+n2)]
  
  A = solve( C_xx + sigma_sq * diag(n))
  
  m = C_sx %*% A %*% y
  v = C_ss - C_sx * A %*% t(C_sx)
  
  return(list(m, v))
}

# apply to utilities data
run2 = function(){
  data = read.csv('utilities.csv') 
  daily_gasbill = log(data$gasbill / data$billingdays)
  
  # want to predict daily_gasbill from temperature
  X = data$temp
  y = daily_gasbill
  res = gp_predict(X, y, X, sigma_sq = 0.01)
  m = res[[1]]
  v = res[[2]]
  
  plot(data$temp, daily_gasbill, xlab = 'Temperature', ylab = 'Log Daily gasbill', ylim = c(-3, 3))
  
  sd_e = diag(v) ^(0.5)
  plotCI(data$temp, p, sd_e*1.96, sd_e * 1.96,  col = 'red', add = TRUE)  
}