library(MASS)

data = read.csv('gdpgrowth.csv')
y = data$GR6096
x = cbind(1, data$DEF60)

n = nrow(x)
p = ncol(x)

# prior
Lambda = diag(n)
K = diag(c(0.01, 0.01))
m = rep(0, p)
d = 2
eta = 2
h = 2

plot(data$DEF60, y, xlab = 'Defense spending', ylab = 'GDP growth rate')

# Gibbs sampling
ng = 1000
save_beta = matrix(nrow = ng, ncol = p)

for (it in 1:ng){
  # posterior params:
  d_star = d + n
  K_star = t(x) %*% Lambda %*% x + K
  m_star = solve(K_star) %*% (t(x) %*% Lambda %*% y + K %*% m)
  eta_star = eta + t(y) %*% Lambda %*% y +  t(m) %*% K %*% m - t(m_star) %*% K_star %*% m_star
  eta_star = eta_star[1,1]
  
  omega = rgamma(1, d_star/2, eta_star/2)
  S = solve(omega * K_star)
  beta = mvrnorm(n = 1, mu = m_star, Sigma = S)
  for (i in 1:n){
    Lambda[i, i] = rgamma(1, (h+1)/2, (h + omega * (y[i] - t(x[i,]) %*% beta)^2)/2 )
  }
  #if (it > 9990){
  #  abline(beta, col = 'blue')
  #}
  save_beta[it,] = beta
}

beta = apply(save_beta, c(2), mean)

abline(beta, col = 'blue')