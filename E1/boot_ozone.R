# Load the library
# you might have to install this the first time
nboots = 10000
library(mlbench)

# Load the data
ozone = data(Ozone, package='mlbench')

# Look at the help file for details
?Ozone

# Scrub the missing values
# Extract the relevant columns 
ozone = na.omit(Ozone)[,4:13]

y = ozone[,1]
x = as.matrix(ozone[,2:10])

# add an intercept
x = cbind(1,x)

boot = matrix(0, ncol = 10, nrow = nboots)
for(i in 1:nboots)
{
  idx = sample(1:nrow(x), nrow(x), replace=TRUE)
  bx = x[idx,]
  by = y[idx]
  betahat = solve(t(bx) %*% bx) %*% t(bx) %*% by
  boot[i, ] = betahat
}

bcov = cov(boot)
sqrt(diag(bcov))

# compute the estimator


# Fill in the blank
# betacov = ?

# Now compare to lm
# the 'minus 1' notation says not to fit an intercept (we've already hard-coded it as an extra column)
lm1 = lm(y~x-1)

summary(lm1)
betacovlm = vcov(lm1)
sqrt(diag(betacovlm))

round(betacovlm, 3) # theory estimate
round(bcov, 3) # bootstrap estimate

