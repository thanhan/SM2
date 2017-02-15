data = read.csv('gdpgrowth.csv')

y = data$GR6096
x = cbind(1, data$DEF60)

#fit a frequentest linear model and plot
lm1 = lm(y~x-1)
plot(data$DEF60, y, xlab = 'Defense spending', ylab = 'GDP growth rate')
abline(lm1$coefficients)

n = nrow(data)
p = 2

# prior
Lambda = diag(n)
K = diag(c(0.01, 0.01))
m = rep(0, p)


#fit the Bayesian model and plot
beta = m - solve(K + t(x) %*% Lambda %*% x ) %*% (t(x) %*% Lambda) %*% (x %*% m - y)
abline(beta, col = 'blue')
