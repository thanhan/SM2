library(mosaic)
library(lme4)

cheese = read.csv('cheese.csv')

cheese$lvol = log(cheese$vol)
cheese$lprice = log(cheese$price)

# split train/test
n = nrow(cheese)
ind <- sample(seq_len(n), size = floor(n * 0.7))
train = cheese[ind, ]
test = cheese[-ind, ]

m0 = lm(lvol ~ lprice, data = train)
p0 = predict(m0, test)
sum((p0 - test$lvol)^2)
#[1] 941.5814

m1 = lm(lvol ~ disp + lprice, data = train)
p1 = predict(m1, test)
sum((p1 - test$lvol)^2)
#[1] 920.1723

m2 = lmer(lvol ~  (lprice | store), data = train)
p2 = predict(m2, test)
sum((p2 - test$lvol)^2)
#[1] 115.4952


m3 = lmer(lvol ~ (1 + lprice| disp) + (1 + lprice | store), data = train)
p3 = predict(m3, test)
sum((p3 - test$lvol)^2)
#[1] 109.2668

m4 = lmer(lvol ~  disp + (lprice | store), data = train)
p4 = predict(m4, test)
sum((p4 - test$lvol)^2)
#[1] 111.3955

m = lmer(lvol ~ (1 + lprice| disp) + (1 + lprice | store), data = cheese)

store_intercept = coef(m2)$store$`(Intercept)`
store_slope     = coef(m2)$store$lprice

plot(1, type = 'n', xlim = c(0, 2), ylim = c(5, 12))
for (i in 1:88)
  abline(store_intercept[i], store_slope[i])
  