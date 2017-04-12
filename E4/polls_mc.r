library(MASS)
library(truncnorm)

polls = read.csv("polls.csv")

#levels(polls$edu) = c("NoHS", "HS", "SomeColl", "Bacc")


#polls$HS = as.numeric(polls$edu == 'HS')
#polls$BA = as.numeric(polls$edu == 'Bacc')
#polls$NH = as.numeric(polls$edu == 'NoHS')
#polls$SC = as.numeric(polls$edu == 'SomeColl')

#polls$a18 = as.numeric(polls$age == '18to29')
#polls$a30 = as.numeric(polls$age == '30to44')
#polls$a45 = as.numeric(polls$age == '45to64')
#polls$a65 = as.numeric(polls$age == '65plus')

#X = polls[ , c('black', 'female', 'NH', 'HS', 'SC', 'BA', 'a18', 'a30', 'a45', 'a65' )]
#X = as.matrix(X)

X = model.matrix( ~ edu + age + female + black, polls)
y = polls$bush

# add state id
polls = transform(polls, sid=as.numeric(factor(state)))

n = nrow(X)
m = ncol(X)
ni = xtabs(~sid, polls)

mu = vector("numeric", 49)
beta = vector("numeric", m)


t = vector("numeric", n)

NMC = 2000
save_beta = matrix(nrow = NMC, ncol = m)
save_mu = matrix(nrow = NMC, ncol = 49)
save_t = matrix(nrow = NMC, ncol = n)


V = solve( t(X) %*% X )

for (it in 1:NMC){
  tmm = t - mu[polls$sid]
  beta_hat = V %*% t(X) %*% (tmm)
  
  beta = mvrnorm(mu = beta_hat, Sigma = V)
  
  sum_s = vector("numeric", 49)
  xb = X %*% beta
  for (i in 1:n){
    s = polls$sid[i]
    sum_s[s] = sum_s[s] + t[i] - xb[i]
  }
  for (s in 1:49) {
    sum_s[s] = sum_s[s] / ni[s]
    mu[s] = rnorm(1, mean = sum_s[s], sd = sqrt(1/ni[s]))
  }
  
  for (i in 1:n){
    if (!is.na(y[i]))
    if (y[i] == 0){
      t[i] = rtruncnorm(1, b = 0, mean = xb[i] + mu[polls$sid[i]], sd = 1)
    }
    else{
      t[i] = rtruncnorm(1, a = 0, mean = xb[i] + mu[polls$sid[i]], sd = 1)
    }
  }
  
  
  save_t[it, ] = t
  save_beta[it, ] = beta
  save_mu[it, ] = mu
  
}




