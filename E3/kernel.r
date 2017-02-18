k_gaussian = function (x){
  return (1.0/sqrt(2*pi) * exp(-x*x/2))
}

k_uniform = function(x){
  if (x >= -1 && x <= 1) return (1.0/2)
  return (0)
  
}

# kernel function with bandwith h and distance d
k = function(h, d, type = "Gaussian"){
  if (type == 'Gaussian'){
    return (1.0/h * k_gaussian(d/h))
  }
  else if (type == 'uniform'){
    return (1.0/h * k_uniform(d/h))
  }
  
}

#return a vector of weight
weight = function(X, x_star, h, kernel_type = "Gaussian"){
  n = length(X)
  w = vector(length = n)
  
  for (i in 1:n){
    d = X[i] - x_star
    w[i] = k(h, d, kernel_type)
  }
  return (w)
}

# make a prediction at x_star
predict = function(X, y, x_star, h, kernel_type = "Gaussian"){
  w = weight(X, x_star, h, kernel_type)
  s = sum(w)
  res = t(w) %*% y / s
  return (res[1,1])
}



fit_kernel_smoother = function(X, Y, l, r, h, kernel_type = "Gaussian"){
  a = seq(l, r, 0.01)
  n = length(a)
  b = vector(length = n)
  
  for (i in 1:n){
    x_star = a[i]
    b[i] = predict(X, Y, x_star, h, kernel_type)
  }
  return (list(a,b))
}


# sample X_i from uniform
# eval f(X_i) then add Normal noise
simulate_data = function(n = 100, sd = 0.1, l = -5, r = 5, f){
  X = runif(n, l, r)
  y = vector(length = n)
  for (i in 1:n){
    y[i] = f(X[i]) + rnorm(1, mean = 0, sd = sd)
  }
  return (list(c(X), c(y)))
}


# simulate data and run the kernel smoother
run = function(h = 0.01, l = -5, r = 5){
  
  temp = simulate_data(100, 0.5, l, r, sin)
  X = temp[[1]]
  Y = temp[[2]]
  
  plot(X, Y)
  
  temp = fit_kernel_smoother(X, Y, l, r, h)
  
  lines(temp[[1]], temp[[2]], col = 'blue')
}

#generate train/test data
gen_train_test = function(fun, l = -5, r = 5, sd, train_size = 100, test_size = 100){
  temp = simulate_data(train_size, sd, l, r, fun)
  X_train = temp[[1]]
  Y_train = temp[[2]]
  
  temp = simulate_data(test_size, sd, l, r, fun)
  X_test = temp[[1]]
  Y_test = temp[[2]]
  return (data.frame(X_train, Y_train, X_test, Y_test))
}

# Do cross validation:
cv = function(data, h = 0.1){
  
  test_size = length(df$X_test)
  sq_error = 0
  for (i in 1:test_size){
    x_star = data$X_test[i]
    y_star = p = predict(data$X_train, data$Y_train, x_star, h)
    sq_error = sq_error + (y_star - data$Y_test[i])^2
  }
  
  sq_error = sq_error / test_size
  return (sqrt(sq_error))
}

run_cv = function(fun, l = 0, r = 1, sd = 0.5){
  #list_h = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10)
  #list_h = c(0.1, 0.09, 0.08, 0.07, 0.06, 0.05)
  list_h = seq(0.01, 0.20, 0.01)
  n = length(list_h)
  rmse   = vector(length = n)
  
  
  data = gen_train_test(fun, l, r, sd)
  for (i in 1:n){
    h = list_h[i]
    rmse[i] = cv(data, h)
  }
  res = data.frame(list_h, rmse)
  return(res)
}

#wiggly function
w_fun = function(x){
  return (sin(x*10))
}

#smooth function
s_fun = function(x){
  return (x*x)
}

# For Cross Validation part (B)
run_cv_b = function(){
  print ("smooth function, not so noisy")
  res = run_cv(s_fun, l = 0, r = 1, sd = 0.05)
  plot(res, col = 'black', type = 'l', ylim = c(0,1), xlab = "h")
  
  print("smooth function, noisy")
  res = run_cv(s_fun, l = 0, r = 1, sd = 0.25)
  lines(res, col = 'blue')
  
  print("wiggly function, not so noisy")
  res = run_cv(w_fun, l = 0, r = 1, sd = 0.05)
  lines(res, col = 'red')
  
  print("wiggly function, noisy")
  res = run_cv(w_fun, l = 0, r = 1, sd = 0.25)
  lines(res, col = 'yellow')
}
  
