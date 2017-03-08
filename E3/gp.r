library(MASS)
library(plotrix)
library(stats)


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
  X = as.matrix(X)
  n = nrow(X)
  d = as.matrix(dist(X))
  
  C = matrix(0, nrow = n, ncol = n)
  for (i in 1:n){
    for (j in 1:n){
      if (cov_fun == 'm52')
        C[i, j] = cov_m52(b, tau1_sq, tau2_sq, d[i,j])
      else if (cov_fun == 'se')
        C[i, j] = cov_se(b, tau1_sq, tau2_sq, d[i,j])
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


gp_predict = function(X, y, x_star, sigma_sq = 1, b = 10, tau1_sq = 5, tau2_sq = 0, cov_fun = 'm52'){
  n1 = nrow(x_star)
  n2 = nrow(X)
  C = compute_c(rbind(x_star, X), b, tau1_sq, tau2_sq, cov_fun)
  C_ss = C[1:n1,1:n1]
  C_sx = C[1:n1, (n1+1):(n1+n2)]
  C_xx = C[(n1+1):(n1+n2), (n1+1):(n1+n2)]
  
  A = solve( C_xx + sigma_sq * diag(n2))
  
  m = C_sx %*% A %*% y
  
  v = C_ss - t(C_sx) %*% A  %*% C_sx
  
  #return (x_star)
  #return(C_xx)
  #return( compute_c(x_star, b, tau1_sq, tau2_sq, cov_fun))
  return(list(m, v))
}

mll = function(X, y, sigma_sq = 1, b = 10, tau1_sq = 5, tau2_sq = 0, cov_fun = 'm52'){
  C = compute_c(X, b, tau1_sq, tau2_sq, cov_fun)
  n = length(X)
  S = C + sigma_sq * diag(n)
  return ( dmvnorm(y, mean = rep(0, n), sigma = S) )
  
}

# apply to utilities data
run2 = function(){
  data = read.csv('utilities.csv') 
  daily_gasbill = log(data$gasbill / data$billingdays)
  
  # want to predict daily_gasbill from temperature
  X = as.matrix(data$temp)
  y = daily_gasbill
  res = gp_predict(X, y, X,tau1_sq = 1, b = 30, cov_fun = 'se')
  m = res[[1]]
  v = res[[2]]
  
  plot(data$temp, daily_gasbill, xlab = 'Temperature', ylab = 'Log Daily gasbill', ylim = c(-3, 3))
  
  sd_e = diag(v) ^(0.5)
  plotCI(data$temp, m, sd_e*1.96, sd_e * 1.96,  col = 'red', add = TRUE)  
}

# optimize parameter using marginal likelihood
run3 = function() {
  data = read.csv('utilities.csv') 
  daily_gasbill = log(data$gasbill / data$billingdays)
  
  # want to predict daily_gasbill from temperature
  X = data$temp
  y = daily_gasbill
  
  for (tau1_sq in c(0.5, 1, 1.5, 2, 1.5, 3))
    for (b in c(5, 10, 15, 20, 25, 30))
      print (paste(tau1_sq, b, mll(X, y, tau1_sq = tau1_sq, b = b)))
  
}

# weather dataset
run4 = function() {
  data = read.csv('weather.csv')
  X = data[, 3:4]
  y = data[, 2]
  
  bin_x = seq(min(X[, 1]), max(X[, 1])); len_x = length(bin_x)
  bin_y = seq(min(X[, 2]), max(X[, 2])); len_y = length(bin_y)
  f = matrix(0, len_x, len_y)
  
  for (i in 1:len_x)
    for (j in 1:len_y){
      x_star = data.frame(lon = bin_x[i], lat = bin_y[j])
      res = gp_predict(X, y, x_star ,tau1_sq = 1, b = 30, cov_fun = 'se')
      m = res[[1]]
      v = res[[2]]    
      f[i, j] = m
    }
  
}