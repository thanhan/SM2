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



fit_kernel_smoother = function(X, Y, r, h, kernel_type = "Gaussian"){
  a = seq(-r, r, 0.01)
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
simulate_data = function(n = 100, sd = 0.1, range = 10, f){
  X = runif(n, -range, range)
  y = vector(length = n)
  for (i in 1:n){
    y[i] = f(X[i]) + rnorm(1, mean = 0, sd = sd)
  }
  return (list(c(X), c(y)))
}



run = function(h = 0.01){
  r = 5 # range of data
  
  temp = simulate_data(100, 0.5, r, sin)
  X = temp[[1]]
  Y = temp[[2]]
  
  plot(X, Y)
  
  temp = fit_kernel_smoother(X, Y, r, h)
  lines(temp[[1]], temp[[2]], col = 'blue')
  
}

