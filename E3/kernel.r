k_gaussian = function (x){
  return (1.0/sqrt(2*pi) * exp(-x*x/2))
}

k_uniform = function(x){
  if (x >= -1 && x <= 1) return (1.0/2)
  return (0)
  
}

k = function(h, d, type = "Gaussian"){
  if (type == 'Gaussian'){
    return (1.0/h * k_gaussian(d/h))
  }
  else if (type == 'uniform'){
    return (1.0/h * k_uniform(d/h))
  }
  
}

weight = function(X, x_star, h, kernel_type = "Gaussian"){
  n = length(X)
  w = vector(length = n)
  
  for (i in 1:n){
    d = X[i] - x_star
    w[i] = k(h, d, kernel_type)
  }
  return (w)
}