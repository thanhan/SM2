data = read.csv('gdpgrowth.csv')

y = data$GR6096
x = cbind(1, data$DEF60)

lm1 = lm(y~x-1)
plot(data$DEF60, y, xlab = 'Defense spending', ylab = 'GDP growth rate')
abline(lm1$coefficients)

n = nrow(data)
p = 2

Lambda = diag(n)
K = diag(0.1, 0.1)
m = rep(0, p)

beta = m - solve(K + T(x) %*% Lambda %*% x ) (T(x) %*% Lambda) (X %*% m - y)