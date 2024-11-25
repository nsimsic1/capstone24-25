# load necessary libraries
library(e1071)
library(ppcc)
library(fitdistrplus)
library(FAdist)

# set working directory
setwd("D:/GoogleDrive/CE5500/Rcode")
source("utils.R")

# load data of annual maxima at Azibe Soltane on the Sebou River in Morocco
maxQ = read.csv("../Data/Problem6.17.csv")

# test if skewness is 0
n = length(maxQ$Flow)
z = skewness(maxQ$Flow, type=2) / sqrt(6*(n-2) / ((n+1)*(n+3)))
p_one_sided = 1 - pnorm(z)
p_two_sided = 2*p_one_sided

# test if skewness of log(data) is 0
z_log = skewness(log(maxQ$Flow), type=2) / sqrt(6*(n-2) / ((n+1)*(n+3)))
p_one_sided_log = 1 - pnorm(z_log)
p_two_sided_log = 2*p_one_sided_log

# test if data is normally distributed with K-S test
norm.fit.mle = fitdist(maxQ$Flow, "norm", method="mle")
ks.test(maxQ$Flow,"pnorm", norm.fit.mle$estimate[1],
        norm.fit.mle$estimate[2], alternative="two.sided")

# test if data is log-normally distributed with K-S test
LN2.fit.mle = fitdist(maxQ$Flow, "lnorm", method="mle")
ks.test(maxQ$Flow,"plnorm", LN2.fit.mle$estimate[1],
        LN2.fit.mle$estimate[2], alternative="two.sided")

# test if data is LN3-distributed with K-S test
LN3.params = lnorm3MOM(maxQ$Flow)
LN3.fit.mle = fitdist(maxQ$Flow, "lnorm3", 
                      start=list(shape=LN3.params$sigma, 
                                 scale=LN3.params$mu, 
                                 thres=LN3.params$tau), 
                      method="mle")
ks.test(maxQ$Flow,"plnorm3", LN3.fit.mle$estimate[1],
        LN3.fit.mle$estimate[2], LN3.fit.mle$estimate[3],
        alternative="two.sided")

# test if data is Gamma-distributed with K-S test
# first fit Gamma distribution
gamma2.fit.mle = fitdist(maxQ$Flow, "gamma", method="mle")
ks.test(maxQ$Flow, "pgamma", gamma2.fit.mle$estimate[1], 
        gamma2.fit.mle$estimate[2], alternative="two.sided")

# try 3-parameter Gamma
gamma3.params = gamma3MOM(maxQ$Flow)
gamma3.fit.mom = 
  fitdist(maxQ$Flow, "gamma3", 
          order=c(1,2,3),
          memp=memp.centered,
          start=list(shape=gamma3.params$alpha, 
                     scale=1/gamma3.params$beta, 
                     thres=gamma3.params$xi),
          method="mme")
ks.test(maxQ$Flow, "pgamma3", gamma3.fit.mom$estimate[1], 
        gamma3.fit.mom$estimate[2], gamma3.fit.mom$estimate[3],
        alternative="two.sided")

# try LP3 distribution
LP3.params = gamma3MOM(log(maxQ$Flow))
LP3.fit.mom = 
  fitdist(log(maxQ$Flow), "gamma3", 
          order=c(1,2,3),
          memp=memp.centered,
          start=list(shape=LP3.params$alpha, 
                     scale=1/LP3.params$beta, 
                     thres=LP3.params$xi),
          method="mme")
ks.test(log(maxQ$Flow), "pgamma3", LP3.fit.mom$estimate[1], 
        LP3.fit.mom$estimate[2], LP3.fit.mom$estimate[3],
        alternative="two.sided")

# try Gumbel distribution
gumbel.params = gumbelMOM(maxQ$Flow)
gumbel.fit.mom = 
  fitdist(maxQ$Flow, "gumbel", 
          order=c(1,2), memp=memp.centered,
          start=list(scale = gumbel.params$alpha, 
                     location = gumbel.params$xi), 
          method="mme")
ks.test(maxQ$Flow, "pgumbel", gumbel.fit.mom$estimate[1], 
        gumbel.fit.mom$estimate[2], alternative="two.sided")

# try GEV distribution
gev.params = gevMOM(maxQ$Flow)
gev.fit.mom = fitdist(maxQ$Flow, "gev", 
                      order=c(1,2,3), memp=memp.centered,
                      start=list(shape=gev.params$kappa, 
                                 scale=gev.params$alpha, 
                                 location=gev.params$xi), 
                      method="mme")
ks.test(maxQ$Flow, "pgev", gev.fit.mom$estimate[1], 
        gev.fit.mom$estimate[2], gev.fit.mom$estimate[3],
        alternative="two.sided")


# perform probability plot correlation test with normal distribution
p = ppPositions(n,"Blom")
cor(qnorm(p, norm.fit.mle$estimate[1], norm.fit.mle$estimate[2]),
    sort(maxQ$Flow))
ppccTest(maxQ$Flow, "qnorm")

# perform probability plot correlation test with log-normal
cor(qlnorm(p, LN2.fit.mle$estimate[1], LN2.fit.mle$estimate[2]), sort(maxQ$Flow))
ppccTest(maxQ$Flow, "qlnorm")
# write your own PPCC function
lnormPPCC(maxQ$Flow, LN2.fit.mle$estimate[1], LN2.fit.mle$estimate[2])

# perform probability plot correlation test with LN3
cor(qlnorm3(p, LN3.fit.mle$estimate[1], LN3.fit.mle$estimate[2],
            LN3.fit.mle$estimate[3]), sort(maxQ$Flow))
# write your own PPCC function
lnorm3PPCC(maxQ$Flow, LN3.fit.mle$estimate[1], LN3.fit.mle$estimate[2],
           LN3.fit.mle$estimate[3])

# perform probability plot correlation test with Gamma
cor(qgamma(p, gamma2.fit.mle$estimate[1], gamma2.fit.mle$estimate[2]), 
    sort(maxQ$Flow))
gamma2PPCC(maxQ$Flow, gamma2.fit.mle$estimate[1], gamma2.fit.mle$estimate[2])

# perform probability plot correlation test with 3 parameter Gamma
cor(qgamma3(p, gamma3.fit.mom$estimate[1], gamma3.fit.mom$estimate[2],
            gamma3.fit.mom$estimate[3]), sort(maxQ$Flow))
ppccTest(maxQ$Flow, "qpearson3", shape=gamma3.fit.mom$estimate[1])
# write your own PPCC function
gamma3PPCC(maxQ$Flow, gamma3.fit.mom$estimate[1], gamma3.fit.mom$estimate[2],
           gamma3.fit.mom$estimate[3])

# perform probability plot correlation test with LP3
cor(qgamma3(p, LP3.fit.mom$estimate[1], LP3.fit.mom$estimate[2],
            LP3.fit.mom$estimate[3]), sort(log(maxQ$Flow)))
ppccTest(log(maxQ$Flow), "qpearson3", shape=LP3.fit.mom$estimate[1])
# write your own PPCC function
gamma3PPCC(log(maxQ$Flow), LP3.fit.mom$estimate[1], LP3.fit.mom$estimate[2],
           LP3.fit.mom$estimate[3])

cor(exp(qgamma3(p, LP3.fit.mom$estimate[1], LP3.fit.mom$estimate[2],
                LP3.fit.mom$estimate[3])), sort(maxQ$Flow))
# write your own PPCC function
LP3PPCC(maxQ$Flow, LP3.fit.mom$estimate[1], LP3.fit.mom$estimate[2],
        LP3.fit.mom$estimate[3])

# perform probability plot correlation test with Gumbel
p = ppPositions(n,"Gringorton")
cor(qgumbel(p, gumbel.fit.mom$estimate[1], gumbel.fit.mom$estimate[2]), 
    sort(maxQ$Flow))
ppccTest(maxQ$Flow, "qgumbel")

# perform probability plot correlation test with GEV
p = ppPositions(n,"Cunane")
cor(qgev(p, gev.fit.mom$estimate[1], gev.fit.mom$estimate[2],
         gev.fit.mom$estimate[3]), sort(maxQ$Flow))
ppccTest(maxQ$Flow, "qgev", shape=gev.fit.mom$estimate[1])
# write your own PPCC function
gevPPCC(maxQ$Flow, gev.fit.mom$estimate[1], gev.fit.mom$estimate[2],
        gev.fit.mom$estimate[3])
