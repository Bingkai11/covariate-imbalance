library(dplyr)
source("adjust_estimator.R")
load("ACDS.rdata")
d <- d[!d$arm == "Vitamin E",]
y <- d$Y18
delta <- is.na(y)
w <- select(d,female,cototscr, mmscore, adtotscr, gdstot, age, Y0)
a <- d$arm == "Donepezil"

# calculate adjusted estimator with missing outcome model
adjust_estimator(y, a, w, delta, method = "DR-WLS-U")

# remove missing outcomes
indi <- !is.na(y)
yr <- y[indi]
wr <- w[indi,]
ar <- a[indi]
w1_bar <- colMeans(wr[ar == 1,])
w0_bar <- colMeans(wr[ar == 0,])
y1_bar <- mean(yr[ar == 1])
y0_bar <- mean(yr[ar == 0])
y1 <- yr[ar == 1]
y0 <- yr[ar == 0]
w1 <- wr[ar == 1,]
w0 <- wr[ar == 0,]
n1 <- length(y1)
n0 <- length(y0)

#caculate unadjusted estimator
y1_bar - y0_bar
adjust_estimator(yr, ar, wr, method = "DR-WLS")

#calcuate conditional bias reduction
ppi <- n1/(n1+n0)
sigma211 <- 0
sigma221 <- 0
for(i in 1: n1){
  bv <- as.numeric(w1[i,])
  sigma211 <- sigma211 + (y1[i] - y1_bar) * (bv - w1_bar)
  sigma221 <- sigma221 + (bv - w1_bar) %*% t((bv - w1_bar))
}
sigma211 <- sigma211/n1
sigma221 <- sigma221/n1

sigma210 <- 0
sigma220 <- 0
for(i in 1: n0){
  bv <- as.numeric(w0[i,])
  sigma210 <- sigma210 + (y0[i] - y0_bar) * (bv - w0_bar)
  sigma220 <- sigma220 + (bv - w0_bar) %*% t((bv - w0_bar))
}
sigma210 <- sigma210/n0
sigma220 <- sigma220/n0
sigma21 <- sigma211/ppi + sigma210/(1-ppi)
sigma22 <- sigma221/ppi + sigma220/(1-ppi)
imbalance <- w1_bar - w0_bar
t(sigma21) %*% solve(sigma22) %*% imbalance
t(sigma21) %*% solve(sigma22) %*% sigma21/(var(y1)/ppi + var(y0)/(1-ppi))
