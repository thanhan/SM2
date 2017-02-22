data = read.csv('utilities.csv')

#prediction with local linear Kernel
predict_llk = function(X, y, x_star, h, kernel_fun, D = 1){
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
  
  #return (K)
  #calculate a
  a = solve(t(R) %*% K %*% R) %*% t(R) %*% K %*% y
  
  return (a[1, 1])
}


run = function(){
  source('kernel.r')
  temp = simulate_data(100, 0.1, -5, 5, sin); X = temp[[1]]; y = temp[[2]]
  
}