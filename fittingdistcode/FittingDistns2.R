# load libraries used last class for calculating the skew and fitting distributions
library(e1071)
library(fitdistrplus)

# set working directory
setwd("D:/GoogleDrive/CE5500/Rcode")

# load data of annual maxima at Azibe Soltane on the Sebou River in Morocco
maxQ = read.csv("../Data/Problem6.17.csv")


# load library for non-standard distributions
library(FAdist)

# calculate parameters from own MOM function - source from utils.R
source("utils.R")
LN3.params = lnorm3MOM(maxQ$Flow)

# use as estimates to fitdist for MLE estimation
LN3.fit.mle = fitdist(maxQ$Flow, "lnorm3", 
                      start=list(shape=LN3.params$sigma, 
                                 scale=LN3.params$mu, 
                                 thres=LN3.params$tau), 
                      method="mle")
LN3.fit.mle$estimate

plot(LN3.fit.mle)


# repeat with method of moments
LN3.fit.mom = fitdist(maxQ$Flow, "lnorm3", 
                      order=c(1,2,3), 
                      memp = memp.centered,
                      start=list(shape=LN3.params$sigma, 
                                 scale=LN3.params$mu, 
                                 thres=LN3.params$tau), 
                      method="mme")

LN3.fit.mom$estimate

plot(LN3.fit.mom)


# Fit Gamma distribution
gamma2.fit.mle = fitdist(maxQ$Flow, "gamma", method="mle")
gamma2.fit.mle$estimate
plot(gamma2.fit.mle)

gamma2.fit.mom = fitdist(maxQ$Flow, "gamma", method="mme")
gamma2.fit.mom$estimate
plot(gamma2.fit.mom)

# compare with theoretical moments
gamma2.params = gamma2MOM(maxQ$Flow)
gamma2.params


# Fit 3-parameter Gamma distribution
gamma3.params = gamma3MOM(maxQ$Flow)
gamma3.fit.mom = 
  fitdist(maxQ$Flow, "gamma3", 
          order=c(1,2,3),
          memp=memp.centered,
          start=list(shape=gamma3.params$alpha, 
                     scale=1/gamma3.params$beta, 
                     thres=gamma3.params$xi),
          method="mme")

gamma3.fit.mom$estimate
plot(gamma3.fit.mom)

gamma3.fit.mle = 
  fitdist(maxQ$Flow, "gamma3",
          start=list(shape=gamma3.params$alpha,
                     scale=1/gamma3.params$beta, 
                     thres=gamma3.params$xi), 
          method="mle")

gamma3.fit.mle$estimate
plot(gamma3.fit.mle)


# Fit LP3 distribution
LP3.params = gamma3MOM(log(maxQ$Flow))
LP3.fit.mom = 
  fitdist(log(maxQ$Flow), "gamma3", 
          order=c(1,2,3),
          memp=memp.centered,
          start=list(shape=LP3.params$alpha, 
                     scale=1/LP3.params$beta, 
                     thres=LP3.params$xi),
          method="mme")

LP3.fit.mom$estimate
plot(LP3.fit.mom)

LP3.fit.mle = 
  fitdist(log(maxQ$Flow), "gamma3", 
          start=list(shape=LP3.params$alpha, 
                     scale=1/LP3.params$beta, 
                     thres=LP3.params$xi), 
          method="mle")

LP3.fit.mle$estimate
plot(LP3.fit.mle)

if(skewness(log(maxQ$Flow), type=2)<0){
  LP3.params = gamma3MOM(-log(maxQ$Flow))
  LP3.fit.mom = fitdist(-log(maxQ$Flow), "gamma3", order=c(1,2,3), memp=memp.centered,
                        start=list(shape=LP3.params$alpha, scale=1/LP3.params$beta, 
                                   thres=LP3.params$xi), 
                        lower=c(-Inf,0,-Inf), upper=c(Inf,Inf,Inf), method="mme")
  LP3.fit.mle = fitdist(-log(maxQ$Flow), "gamma3",
                        start=list(shape=LP3.params$alpha, scale=1/LP3.params$beta, 
                                   thres=LP3.params$xi), 
                        method="mle")
}


# Fit Gumbel Distribution
gumbel.params = gumbelMOM(maxQ$Flow)
gumbel.fit.mom = 
  fitdist(maxQ$Flow, "gumbel", 
          order=c(1,2), memp=memp.centered,
          start=list(scale = gumbel.params$alpha, 
                     location = gumbel.params$xi), 
          method="mme")

gumbel.fit.mom$estimate

plot(gumbel.fit.mom)

gumbel.fit.mle = 
  fitdist(maxQ$Flow, "gumbel", 
          start=list(scale = gumbel.params$alpha, 
                     location = gumbel.params$xi), 
          method="mle")

gumbel.fit.mle

plot(gumbel.fit.mle)


# Fit GEV distribution
gev.params = gevMOM(maxQ$Flow)
gev.fit.mom = fitdist(maxQ$Flow, "gev", 
                      order=c(1,2,3), memp=memp.centered,
                      start=list(shape=gev.params$kappa, 
                                 scale=gev.params$alpha, 
                                 location=gev.params$xi), 
                      method="mme")

gev.fit.mom$estimate
plot(gev.fit.mom)

gev.fit.mle = fitdist(maxQ$Flow, "gev", 
                      start=list(shape=gev.params$kappa, 
                                 scale=gev.params$alpha, 
                                 location=gev.params$xi), 
                      method="mle")

gev.fit.mle$estimate
plot(gev.fit.mle)