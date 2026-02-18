


library(MASS) # For mvrnorm to generate multivariate normal distributions
library(GeneralisedCovarianceMeasure) # Generalised Covariance Measure (GCM)
library(CondIndTests) # kernel-based CI test (KCI)
library(scales)
library(Matrix)
library(MASS)
library(lme4)
library(mgcv)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(dplyr)
library(gamm4)



generate_two_different_functions <- function() {
  functions <- list(
    linear = function(z) z,
    square = function(z) z^2,
    cubic = function(z) z^3,
    tanh = function(z) tanh(z),
    negexp = function(z) exp(-sqrt(z^2))
    # ,
    # cosin = function(z) cos(z)
  )
  
  indices <- sample(length(functions), 2)
  f <- functions[[indices[1]]]
  g <- functions[[indices[2]]]
  
  return(list(f = f, g = g))
}
functions_x_y <- generate_two_different_functions()
f_x <- functions_x_y$f
g_y <- functions_x_y$g


ni = 13                                                 # number of subjects
RE = mvrnorm(ni, mu=c(0,0), Sigma=rbind(c(1.0, 0.3),
                                        c(0.3, 1.0) ))
alpha_vector <- rnorm(ni, mean = 0, sd = 5)
beta_vector <- rnorm(ni, mean = 0, sd = 5)
alpha_vector
colnames(RE) = c("ints","slopes");  t(round(RE,2))



nj   = 10                              # number of timepoints
data = data.frame(ID   = rep(1:ni,   each=nj), 
                  time = rep(1:nj,   times=ni),
                  RE.i = rep(RE[,1], each=nj),
                  RE.s = rep(RE[,2], each=nj),
                  z=NA,
                  x = NA,
                  y    = NA                    )
head(data, 14)
z <- with(data, rnorm(n=ni*nj, mean=0, sd=1))
data$z  = z
x       = with(data, (0 + RE.i) + f_x(z) + (.3 + RE.s)*time + rnorm(n=ni*nj, mean=0, sd=1))
y       = with(data, (0 + RE.i) + g_y(z) + (.3 + RE.s)*time + rnorm(n=ni*nj, mean=0, sd=1))
m1       = rbinom(n=ni*nj, size=1, prob=.1)  
m2       = rbinom(n=ni*nj, size=1, prob=.1) 
m3 = rbinom(n=ni*nj, size=1, prob=.1)
y[m1==1] = NA
x[m2==1] = NA
z[m3==1] = NA
data$y  = y
data$x  = x
data$z  = z
head(data, 14)


res_x <- comp.resids(data$x, data$z, regr.pars = list(), regr.method = "gam")
res_y <- comp.resids(data$y, data$z, regr.pars = list(), regr.method = "gam")

res_x
res_y
plot(res_y)
plot(res_x)


pval_egcm <- gcm.test(resid.XonZ = res_x, resid.YonZ = res_y)$p.value
pval_egcm
gamm_re_y <- gamm4(y~ s(z), data = data, random = ~(time|ID))$mer
gamm_re_x <- gamm4(x~ s(z), data = data, random = ~(time|ID))$mer
res_re_y <- resid(gamm_re_y)
res_re_x <- resid(gamm_re_x)
plot(res_re_y)
plot(res_re_x)

pval_egcm <- gcm.test(resid.XonZ = res_re_x, resid.YonZ = res_re_y)$p.value
pval_egcm


#####################

