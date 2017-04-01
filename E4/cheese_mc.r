
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

taua = 1
taub = 1
taug = 1

# setup
NMC = 1000
save_alpha = matrix(nrow = NMC, ncol = m)
save_beta  = matrix(nrow = NMC, ncol = m)
save_gamma = matrix(nrow = NMC, ncol = m)

save_mua = vector("numeric", NMC)
save_mub = vector("numeric", NMC)
save_mug = vector("numeric", NMC)

save_taua = vector("numeric", NMC)
save_taub = vector("numeric", NMC)
save_taug = vector("numeric", NMC)

# sampling
for (it in 1:NMC){
  
  
}
