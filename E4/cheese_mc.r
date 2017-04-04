
cheese = read.csv('cheese.csv')

cheese$lvol = log(cheese$vol)
cheese$lprice = log(cheese$price)

#add store id
cheese = transform(cheese,sid=as.numeric(factor(store)))
n = nrow(cheese)
m = max(cheese$sid) # number of stores

#model: lvol = alpha[sid] + lprice * beta[sid] + disp * gamma[sid] + epsilon

# init parameters:
alpha = vector("numeric", m)
beta  = vector("numeric", m)
gamma = vector("numeric", m)

mua =  0
mub =  0
mug =  0

tauv = 10
taua = 10
taub = 10
taug = 10

# setup
NMC = 20
save_alpha = matrix(nrow = NMC, ncol = m)
save_beta  = matrix(nrow = NMC, ncol = m)
save_gamma = matrix(nrow = NMC, ncol = m)

save_mua = vector("numeric", NMC)
save_mub = vector("numeric", NMC)
save_mug = vector("numeric", NMC)

save_tauv = vector("numeric", NMC)
save_taua = vector("numeric", NMC)
save_taub = vector("numeric", NMC)
save_taug = vector("numeric", NMC)


# sampling
for (it in 1:NMC){
  if (it %% 100 == 0) print(it)
  # sample alpha
  sum_a = rep(taua, m)
  sum_b = rep(taua*mua, m)
  for (i in 1:n){
    j = cheese$sid[i]
    sum_a[j] = sum_a[j] + tauv
    sum_b[j] = sum_b[j] + cheese$lvol[i] - cheese$lprice[i] * beta[j] - cheese$disp[i] * gamma[j]
  }
  
  for (j in 1:m){
    # posterior mean and precision
    pm = sum_b[j] / sum_a[j]
    pp = sum_a[j]
    alpha[j] = rnorm(1, pm, sqrt(1/pp))
  }
   
  # sample beta 
  sum_a = rep(taub, m)
  sum_b = rep(taub*mub, m)
  for (i in 1:n){
    j = cheese$sid[i]
    sum_a[j] = sum_a[j] + tauv * cheese$lprice[i]^2
    sum_b[j] = sum_b[j] + cheese$lvol[i] - alpha[j] - cheese$disp[i] * gamma[j]
  }
  
  for (j in 1:m){
    # posterior mean and precision
    pm = sum_b[j] / sum_a[j]
    pp = sum_a[j]
    beta[j] = rnorm(1, pm, sqrt(1/pp))
  }
  
  # sample gamma
  sum_a = rep(taug, m)
  sum_b = rep(taug*mug, m)
  for (i in 1:n){
    j = cheese$sid[i]
    sum_a[j] = sum_a[j] + tauv * cheese$disp[i]^2
    sum_b[j] = sum_b[j] + cheese$lvol[i] - alpha[j] - cheese$lprice[i] * beta[j]
  }
  
  for (j in 1:m){
    # posterior mean and precision
    pm = sum_b[j] / sum_a[j]
    pp = sum_a[j]
    gamma[j] = rnorm(1, pm, sqrt(1/pp))
  }
  
  # sample mu
  mua = rnorm(1, mean(alpha), sqrt(1/ (taua * m)) )
  mub = rnorm(1, mean(beta), sqrt(1/ (taub * m)) )
  mug = rnorm(1, mean(gamma), sqrt(1/ (taug * m)) )
  
  # sample tau
  taua = rgamma(1, m/2 + 1, sum((alpha - mua)^2)/2 )
  taub = rgamma(1, m/2 + 1, sum((beta - mub)^2)/2 )
  taug = rgamma(1, m/2 + 1, sum((gamma - mug)^2)/2 )
  ss = 0
  for (i in 1:n){
    j = cheese$sid[i]
    ss = ss + (alpha[j] + cheese$lprice[i]* beta[j] + cheese$disp[i] * gamma[j] - cheese$lvol[i])^2
  } 
  tauv = rgamma(1, n/2 + 1, ss/2 )
  
  # save samples
  save_alpha[it,] = alpha
  save_beta[it,] = beta
  save_gamma[it,] = gamma
  
  save_mua[it] = mua
  save_mub[it] = mub
  save_mug[it] = mug
  
  save_taua[it] = taua
  save_taub[it] = taub
  save_taug[it] = taug
}
