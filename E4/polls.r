library(mosaic)
library(lme4)

polls = read.csv("polls.csv")

#gm1 = glm(bush ~ factor(state) + factor(edu) + factor(age) + female + black, 
#          family=binomial(link="probit"), data = polls)


gm1 = glm(bush ~ state + edu + age + female + black, 
                    family=binomial(link="probit"), data = polls)


xtabs(~state, data=polls)

gm2 = glmer(bush ~ (1|state) + edu + age + female + black, 
          family=binomial(link="probit"), data = polls)
