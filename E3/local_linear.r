data = read.csv('utilities.csv')

#prediction with local linear Kernel
# if return_i >=0 then return the (return_i) value in the hat matrix row
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
  
  if (return_i < 0)
    return (a[1, 1])
  else
    return (data.frame(y_star = a[1, 1], hii = H[1, return_i]))
}

# calculate the loocv error
loocv = function(X, y, h, kernel_fun, D = 1){
  n = length(X)
  error = 0
  for (i in 1:n){
    x_star = X[i]
    res = predict_llk(X, y, x_star, h, kernel_fun, D, return_i = i)
    error = error + ( (y[i] - res$y_star) / (1 - res$hii) )^2
    
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


run = function(){
  source('kernel.r')
  #temp = simulate_data(100, 0.1, -5, 5, sin); X = temp[[1]]; y = temp[[2]]
  
  daily_gasbill = data$gasbill / data$billingdays
  
  # print h with the loocv
  list_h = seq(1, 20, 1)
  for (h in list_h)
    print(paste(h, loocv(data$temp, daily_gasbill, h, k_gaussian)))
  # h = 7 is the best
  plot(data$temp, daily_gasbill, xlab = 'Temperature', ylab = 'Daily gasbill')
  points(data$temp, p, col = 'red')
  legend('topright', c("data", "prediction"), pch = c(1,1), col = c("black", "red"))
}