# load necessary libraries
library(fitdistrplus)
library(FAdist)
library(ismev)
library(POT)
library(imputeTS)
library(dplyr)

# set working directory
setwd("D:/GoogleDrive/CE5500/Rcode")

# load data of annual maxima at Azibe Soltane on the Sebou River in Morocco
allQ = read.csv("../Data/Potomac_DC.csv")

# get maxima each water year
allQ$Month = as.numeric(format(as.Date(allQ$Date, format="%m/%d/%Y"),"%m"))
allQ$Year = as.numeric(format(as.Date(allQ$Date, format="%m/%d/%Y"),"%Y"))
allQ$WY = allQ$Year
allQ$WY[which(allQ$Month>=10)] = allQ$WY[which(allQ$Month>=10)] + 1

maxQ = aggregate(allQ$Flow, by=list(allQ$WY), max)
maxQ = maxQ %>% rename(WY=Group.1, Flow=x) # in dplyr package

# make a mean residual life plot of the Potomac flows at DC
mrl.plot(allQ$Flow) # in ismev package

# impute missing data
allQ$Flow = na_interpolation(allQ$Flow) # in imputeTS package

mrl.plot(allQ$Flow)

# choose a threshold for POT
x0 = 1.25e5

# find clusters of storms above threshold so we choose one exceedance per cluster
?clust # in POT package

# rename dataframe columns to format required by clust
newQ = allQ %>% rename(time=Date, obs=Flow)
newQ$time = c(seq(1,length(newQ$time)))

# find clusters and their peaks
clusters = clust(newQ, x0, time.cond=5, clust.max=TRUE)
peaks = clusters[,2]


# fit GPD to peaks using fitdistr package (and FAdist)
source("utils.R")
gp.Params = gpMOM(peaks-x0)
gp.fit.mle = fitdist(peaks-x0, "gp", 
                       start=list(shape=gp.Params$kappa,
                                  scale=gp.Params$alpha),
                       method="mle")
gp.fit.mom = fitdist(peaks-x0, "gp", order = c(1,2), memp=memp.centered, 
                      start=list(shape=gp.Params$kappa,
                                 scale=gp.Params$alpha),
                      method="mme")

gp.fit.mle
plot(gp.fit.mle)

gp.fit.mom
plot(gp.fit.mom)

# find arrival rate of floods above the threshold each year
lambda = length(peaks) / nrow(maxQ)

# estimate 100-yr flood 
# find equivalent parameters of GEV distribution
gev.Params.POT.mle = GPDtoGEV(gp.fit.mle, x0, lambda)
gev.Params.POT.mom = GPDtoGEV(gp.fit.mom, x0, lambda)

q100.POT.mle = qgev(1-1/100, shape=gev.Params.POT.mle$kappa, 
                    scale=gev.Params.POT.mle$alpha, 
                    location=gev.Params.POT.mle$xi)
q100.POT.mom = qgev(1-1/100, shape=gev.Params.POT.mom$kappa, 
                    scale=gev.Params.POT.mom$alpha, 
                    location=gev.Params.POT.mom$xi)

# estimate 100-yr flood by directly fitting GEV to annual maxima
gev.Params = gevMOM(maxQ$Flow)
gev.fit.mle = fitdist(maxQ$Flow, "gev", 
                     start=list(shape=gev.Params$kappa,
                                scale=gev.Params$alpha,
                                location=gev.Params$xi),
                     method="mle")
gev.fit.mom = fitdist(maxQ$Flow, "gev", order = c(1,2,3), memp=memp.centered,
                      start=list(shape=gev.Params$kappa,
                                 scale=gev.Params$alpha,
                                 location=gev.Params$xi),
                     method="mme")

q100.AMS.mle = qgev(1-1/100, shape=gev.fit.mle$estimate[1], 
                    scale=gev.fit.mle$estimate[2], 
                    location=gev.fit.mle$estimate[3])
q100.AMS.mom = qgev(1-1/100, shape=gev.fit.mom$estimate[1], 
                    scale=gev.fit.mom$estimate[2], 
                    location=gev.fit.mom$estimate[3])

