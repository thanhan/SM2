library(mosaic)

data = read.csv('mathtest.csv')
y = data$mathscore
sch = data$school
# number of scores
n = nrow(data)
# average score for each school
y_bar = mean(mathscore ~ school, data = data)
# number of schools
m = length(y_bar)
# number of scores for each school
ni = tally( ~ school, data = data)



NMC = 2000

theta = vector("numeric", m)
mu = mean(y)
t = 1.0      # 1 / sigma^2
lambda = 1.0 # 1 / tau^2

save.theta = matrix(nrow = NMC, ncol = m)
save.mu     = vector("numeric", NMC)
save.t      = vector("numeric", NMC)
save.lambda = vector("numeric", NMC)

p = vector("numeric", m)
theta_hat = vector("numeric", m)

for (it in 1:NMC){
  if (it %% 100 == 0) print(it)
  for (i in 1:m){
    p[i] = t * lambda + t * ni[i]
    theta_hat[i] = (mu * t * lambda + y_bar[i] * ni[i] * t) / p[i]
    theta[i] = rnorm(1, theta_hat[i], sqrt(1 / p[i]))
  }
  
  mu = rnorm(1, mean(theta), sqrt(m * t * lambda))
  
  s = sum((theta - mu)^2)
  for (j in 1:n)
    s = s + (y[j] - theta[sch[j]])^2
  t = rgamma(1, (m+n)/2, s/2)
  
  s = sum((theta - mu)^2)
  lambda = rgamma(1, 0.1 + m/2, 0.1 + t*s/2)
  
  save.theta[it, ] = theta
  save.mu[it] = mu
  save.t[it] = t
  save.lambda[it] = lambda
}


sc = vector("numeric", m)
for (i in 1:m){
  sc[i] = (y_bar[i] - mean(save.theta[,i]))/ y_bar[i]
}

