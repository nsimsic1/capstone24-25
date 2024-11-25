# load necessary libraries
library(fitdistrplus)
library(FAdist)
library(FlowScreen)

# set working directory
setwd("D:/GoogleDrive/CE5500/Rcode")

# load Colorado River daily flow data
allQ = read.flows("../Data/ColoradoStateLine.csv")

# convert it to a time series for FlowScreen
allQ.ts = create.ts(allQ)

# find minimum 7-day flows each year
sevenQ = MAMn(allQ.ts)

# fit 2-parameter Weibull distribution to 7-day minima
source("utils.R")
weibull2.params = Weibull2MOM(as.vector(sevenQ))
weibull2.fit.mle = fitdist(as.vector(sevenQ), "weibull", method="mle")
weibull2.fit.mom = fitdist(as.vector(sevenQ), "weibull", order=c(1,2),
                           memp=memp.centered,
                           start = list(shape=weibull2.params$kappa,
                                        scale=weibull2.params$alpha),
                           method="mme")

plot(weibull2.fit.mle)
plot(weibull2.fit.mom)

# find 7Q10
sevenQ10_weibull2.mle = qweibull(0.1, weibull2.fit.mle$estimate[1], 
                                 weibull2.fit.mle$estimate[2])
sevenQ10_weibull2.mom = qweibull(0.1, weibull2.fit.mom$estimate[1], 
                                 weibull2.fit.mom$estimate[2])

# fit 3-parameter Weibull distribution to 7-day minima
weibull3.params = Weibull3MOM(as.vector(sevenQ))
weibull3.fit.mle = fitdist(as.vector(sevenQ), "weibull3", 
                           start = list(shape=weibull3.params$kappa,
                                        scale=weibull3.params$alpha,
                                        thres=weibull3.params$xi),
                           method="mle")
weibull3.fit.mom = fitdist(as.vector(sevenQ), "weibull3", order=c(1,2,3),
                           memp=memp.centered,
                           start = list(shape=weibull3.params$kappa,
                                        scale=weibull3.params$alpha,
                                        thres=weibull3.params$xi),
                           method="mme")

plot(weibull3.fit.mle)
plot(weibull3.fit.mom)

# find 7Q10
sevenQ10_weibull3.mle = qweibull3(0.1, weibull3.fit.mle$estimate[1], 
                                 weibull3.fit.mle$estimate[2],
                                 weibull3.fit.mle$estimate[3])
sevenQ10_weibull3.mom = qweibull3(0.1, weibull3.fit.mom$estimate[1], 
                                 weibull3.fit.mom$estimate[2],
                                 weibull3.fit.mom$estimate[3])
