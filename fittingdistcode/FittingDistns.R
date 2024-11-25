# set working directory
setwd("D:/GoogleDrive/CE5500/Rcode")

# load data of annual maxima at Azibe Soltane on the Sebou River in Morocco
maxQ = read.csv("../Data/Problem6.17.csv")

# look at the data; construct a histogram
hist(maxQ$Flow, freq=FALSE, xlab="Max Annual Flow", ylab="Density", main="")

library(ggplot2)
ggplot() + aes(x=maxQ$Flow) + 
  geom_histogram(aes(y=..density..), fill=NA, colour="black", 
                 bins = nclass.Sturges(maxQ$Flow)) + 
  labs(x="Max Annual Flow", y="Density") 

# calculate the unbiased sample mean, variance and skew
mu = mean(maxQ$Flow)
sigma = sd(maxQ$Flow)
sigmaSq = var(maxQ$Flow)

library(e1071)
?skewness
gamma = skewness(maxQ$Flow, type=2)

# fit a normal distribution to the data using MLE
library(fitdistrplus)
?fitdist
norm.fit.mle = fitdist(maxQ$Flow, "norm", method="mle")

# find the parameter estimates
norm.fit.mle$estimate

# plot the fit
plot(norm.fit.mle)



# clearly we need a skewed distribution
# fit two-parameter log-normal distribution with MLE
LN2.fit.mle = fitdist(maxQ$Flow, "lnorm", method="mle")
LN2.fit.mle$estimate

# compare with estimates calculated in class
LN2.mu.mle = mean(log(maxQ$Flow))
LN2.sigma.mle = sqrt( mean( (log(maxQ$Flow) - LN2.mu.mle)^2 ) )

LN2.mu.mle
LN2.sigma.mle

# plot fit
plot(LN2.fit.mle)

# estimate 100-yr and 500-yr floods
q100.LN2.mle = qlnorm(1-1/100, LN2.fit.mle$estimate[1], LN2.fit.mle$estimate[2])
q500.LN2.mle = qlnorm(1-1/500, LN2.fit.mle$estimate[1], LN2.fit.mle$estimate[2])


# repeat with MOM
LN2.fit.mom = fitdist(maxQ$Flow, "lnorm", method="mme") # methods matching estimation
LN2.fit.mom$estimate

# compare with estimates calculated in class
LN2.sigma.mom = sqrt(log(sigma^2/mu^2+1))
LN2.mu.mom = log(mu) - LN2.sigma.mom^2/2

LN2.mu.mom
LN2.sigma.mom

# plot fit
plot(LN2.fit.mom)

# estimate 100-yr and 500-yr floods
q100.LN2.mom = qlnorm(1-1/100, LN2.fit.mom$estimate[1], LN2.fit.mom$estimate[2])
q500.LN2.mom = qlnorm(1-1/500, LN2.fit.mom$estimate[1], LN2.fit.mom$estimate[2])
