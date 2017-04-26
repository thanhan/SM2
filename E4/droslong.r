library(lattice)
library(mvtnorm)
library(MCMCpack)


d = read.csv('droslong.csv')

xyplot(log2exp~time | gene, data=d)
xyplot(log2exp~time | group, data=d)

# model: log2exp = a[group] + b1[gene] + b2[gene] * time + b3[gene] * time2 + e
# e ~ N(0, pre = p)
# a ~ N(ma, pre = pa)
# (b1, b2, b3) ~ N(mb, var = sb)
# (mb, sb) ~ NIW(mu0, k0, L0, nu0)


d = transform(d,geid=as.numeric(factor(gene)))
d = transform(d,grid=as.numeric(factor(group)))

ni = c(180, 180, 144)
n = nrow(d)
n_gr = 3
n_ge = 14

a = vector("numeric", n_gr)

b = matrix(nrow = 3, ncol = n_ge, 0)

p = 1
mb = c(0, 0, 0)
sb = diag(3)

# prior for NIW
mu0 = c(0, 0, 0)
k0 = 1
L0 = diag(3)
nu0 = 4

# prior for a
ma = 0
pa = 1

sum_xx = array(dim =  c(n_ge, 3, 3), 0)
for (i in 1:n){
  ge = d$geid[i]
  t = d$time[i]
  sum_xx[ge,,] = sum_xx[ge,,] + c(1, t, t*t) %*% t(c(1,2,3))
}

NMC = 1000

for (it in 1:NMC){
  if (it %% 100 == 0) print(it)
  # sample a
  
  s = vector("numeric", n_gr)
  for (i in 1:n){
    ge = d$geid[i]
    t = d$time[i]
    s = s + d$log2exp[i] - b[1, ge] - b[2, ge] * t - b[3, ge] * t*t
  }
  
  for (i in 1:n_gr){
    new_p = ni[i] * p + pa
    new_mu = (p * s[i] + ma * pa) / new_p
    a[i] = rnorm(1, new_mu, sqrt(1/new_p))
  }
  
  # sample b
  sum_xlma = matrix(nrow = n_ge, ncol = 3, 0)
  for (i in 1:n){
    ge = d$geid[i]
    gr = d$grid[i]
    t = d$time[i]
    sum_xlma[ge,] = sum_xlma[ge,] + c(1, t, t*t) * (d$log2exp[i] - a[gr])
  }
  sb_inv = solve(sb)
  for (i in 1:n_ge){ # i = the gene
    b_pre = sb_inv + p * sum_xx[i]
    b_var = solve(b_pre)
    b_m   = b_var %*% (sb_inv %*% mb + p * sum_xlma[i])
    b[,i] = rmvnorm(1, mean = b_m, sigma = b_var)
  }
  
  # sample mb, sb
  b_bar = rowMeans(b)
  b_S = matrix(nrow = 3, ncol = 3, 0)
  for (i in 1:n_ge){
    b_S = b_S + (b[,i] - b_bar) %*% t(b[,i] - b_bar)
  }
  
  k1 = k0 + n_ge
  mu1 = (k0 / k1) * mu0 +  (n_ge / k1) * b_bar
  L1 = L0 + b_S + k0 * n_ge / (k0 + n_ge) * (b_bar - mu0) %*% t(b_bar - mu0)
  nu1 = nu0 + n_ge
  
  sb = riwish(nu1, L1)
  
  mb = t(rmvnorm(1, mean = mu1, sb / k1))
  
  
  # sample p
  sum_lmp = 0
  for (i in 1:n){
    ge = d$geid[i]
    gr = d$grid[i]
    sum_lmp = sum_lmp + (d$log2exp - a[gr] - b[1, ge] - b[2, ge] * t - b[3, ge] * t*t)^2
  }
  p = rgamma(1, shape = 1.5, rate = sum_lmp/2)

  # sample ma and pa
  sum_ama = (a[1] - ma)^2 + (a[2] - ma)^2 + (a[3] - ma)^2
  pa = rgamma(1, shape = 1.5, rate = sum_ama / 2)
  ma= rnorm(1, mean(a), 1/ sqrt(3 * pa))
  
}

