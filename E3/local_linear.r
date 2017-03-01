library(plotrix)

data = read.csv('utilities.csv')

#prediction with local linear Kernel
# if return_i > 0 then return the row the hat matrix row
predict_llk = function(X, y, x_star, h, kernel_fun, D = 1, return_i = -1){
  n = length(y)
  
  #create the matrix R
  R = matrix(nrow = n, ncol = D + 1)
  for (i in 1:n){
    R[i, 1] = 1
    for (j in 2:(D+1)){
        R[i, j] = (X[i] - x_star)^(j-1)
    }
  }
  
  
  #create diagonal matrix K
  K = vector(length = n)
  for (i in 1:n) K[i] = 1/h * kernel_fun( (X[i] - x_star) / h )
  K = diag(K)
  
  # the first row of this matrix is the row in the Hat matrix
  H = solve(t(R) %*% K %*% R) %*% t(R) %*% K 
  
  #calculate polynomial coefficient a
  a = H %*% y
  
  if (return_i <= 0)
    return (a[1, 1])
  else
    return (list(a[1, 1], H[1, ]))
}

# calculate the loocv error
loocv = function(X, y, h, kernel_fun, D = 1){
  n = length(X)
  error = 0
  for (i in 1:n){
    x_star = X[i]
    res = predict_llk(X, y, x_star, h, kernel_fun, D, return_i = 1)
    error = error + ( (y[i] - res[[1]]) / (1 - res[[2]][i]) )^2
    
  }
  return (error)
}

#predict y_star for each x
predict_llk_all = function(X, y, h, kernel_fun, D = 1){
  n = length(X)
  error = 0
  y_star = vector(length = n)
  
  for (i in 1:n){
    x_star = X[i]
    res = predict_llk(X, y, x_star, h, kernel_fun, D)
    y_star[i] = res
  }
  
  return (y_star)
}

get_hat_matrix = function(X, y, h, kernel_fun, D = 1){
  n = length(X)
  error = 0
  y_star = vector(length = n)
  H = matrix(nrow = n, ncol = n)
  
  for (i in 1:n){
    x_star = X[i]
    res = predict_llk(X, y, x_star, h, kernel_fun, D, return_i = 1)
    H[i,] = res[[2]]
  }
  
  return (H)
}


run = function(){
  source('kernel.r')
  temp = simulate_data(100, 0.1, -5, 5, sin); X = temp[[1]]; y = temp[[2]]
  n = length(data$temp)
  
  #calculate the daily gasbill and take log
  daily_gasbill = log(data$gasbill / data$billingdays)
  
  data$temp = data$temp
  
  # print h with the loocv
  list_h = seq(1, 20, 1)
  for (h in list_h)
    print(paste(h, loocv(data$temp, daily_gasbill, h, k_gaussian)))
  # h = 7 is the best
  plot(data$temp, daily_gasbill, xlab = 'Temperature', ylab = 'Log Daily gasbill')
  p = predict_llk_all(data$temp, daily_gasbill, h = 7, kernel_fun = k_gaussian)
  #points(data$temp, p, col = 'red')
  legend('topright', c("data", "prediction"), pch = c(1,1), col = c("black", "red"))
  
  H = get_hat_matrix(data$temp, daily_gasbill, h = 7, kernel_fun = k_gaussian)
  r = abs(p - daily_gasbill)
  
  sigma_sq = sum(r*r) / ( n - 2 * sum(diag(H)) + sum(diag(t(H) %*% H)))
  v = sigma_sq * sum(diag(H) * diag(H))
  sd_e = v ^(0.5)
  plotCI(data$temp, p, sd_e*1.96, sd_e * 1.96,  col = 'red', add = TRUE)
  
}
