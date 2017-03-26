library(mosaic)
library(lme4)

cheese = read.csv('cheese.csv')

cheese$lvol = log(cheese$vol)
cheese$lprice = log(cheese$price)

m1 = lmer(lvol ~ (1 + lprice| disp) + (1 + lprice | store), data = cheese)
