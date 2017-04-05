
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

tauv = 1
taua = 1
taub = 1
taug = 1

# setup
NMC = 2000
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
    sum_b[j] = sum_b[j] + tauv * (cheese$lvol[i] - cheese$lprice[i] * beta[j] - cheese$disp[i] * gamma[j])
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
    sum_b[j] = sum_b[j] + tauv * cheese$lprice[i] * (cheese$lvol[i] - alpha[j] - cheese$disp[i] * gamma[j])
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
    sum_b[j] = sum_b[j] + tauv * cheese$disp[i] * (cheese$lvol[i] - alpha[j] - cheese$lprice[i] * beta[j])
  }
  
  for (j in 1:m){
    # posterior mean and precision
    pm = sum_b[j] / sum_a[j]
    pp = sum_a[j]
    gamma[j] = rnorm(1, pm, sqrt(1/pp))
  }
  
  # sample mu
  pa = taua * m + 1
  mua = rnorm(1, taua * m * mean(alpha) / pa, sqrt(1/ pa) )
  pb = taub * m + 1 
  mub = rnorm(1, taub * m * mean(beta) / pb, sqrt(1/ pb) )
  pg = taug * m + 1
  mug = rnorm(1, taug * m * mean(gamma) / pg, sqrt(1/ pg) )
  
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

store_names = rep("abc", m)
for (i in 1:n){
  j = cheese$sid[i]
  store_names[j] = as.character(cheese$store[i])
}

disp_ad = rep(0, m)
for (i in 1:n){
  j = cheese$sid[i]
  disp_ad[j] = disp_ad[j] + cheese$disp[i]
}

# plot
par(mfrow = c(4, 4))
for (i in 1:16)
  hist(save_beta[1000:2000,i], xlab = "beta", ylab = 'Freq', main = store_names[i], xlim = c(-5, 0), ylim = c(0, 400))

par(mfrow = c(4, 4))
for (i in 1:16)
  hist(save_gamma[1000:2000,i], xlab = cat("gamma", disp_ad[i]), ylab = 'Freq', main = store_names[i], xlim = c(-1, 1), ylim = c(0, 400))

par(mfrow = c(1, 1))
hist(save_gamma[1000:2000,])

# compare to lmer
