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


n_subjects <- 500

b <- 0


generate_two_different_functions <- function() {
  functions <- list(
    linear = function(z) z,
    square = function(z) z^2,
    cubic = function(z) z^3
   ,
   tanh = function(z) tanh(z),
     negexp = function(z) exp(-sqrt(z^2)),
     cosin = function(z) cos(z)
  )
  
  indices <- sample(length(functions), 2)
  f <- functions[[indices[1]]]
  g <- functions[[indices[2]]]
  
  return(list(f = f, g = g))
}


data <- data.frame(subject = integer(0),
                   measurement = integer(0),
                   Z1 = numeric(0),
                   alpha=numeric(0),
                   beta=numeric(0))


alpha_vector <- rnorm(n_subjects, mean = 0, sd = 5)
beta_vector <- rnorm(n_subjects, mean = 0, sd = 5)

for (i in 1:n_subjects) {
  n_measurements <- round(rnorm(1, mean = 5, sd = 1))
  n_measurements <- max(min(n_measurements, 10), 1)
  cor_matrix_i <- matrix(0.5, nrow = n_measurements, ncol = n_measurements)
  diag(cor_matrix_i) <- 1
  Z1 <- mvrnorm(n = 1, mu = rep(2, n_measurements), Sigma = cor_matrix_i)
  # Z2 <- mvrnorm(n = 1, mu = rep(0, n_measurements), Sigma = cor_matrix_i)
  # Z3 <- mvrnorm(n = 1, mu = rep(0, n_measurements), Sigma = cor_matrix_i)
  # Z <- cbind(Z1,Z2,Z3)
  
  alpha <- rep(alpha_vector[i], n_measurements)
  beta <- rep(beta_vector[i], n_measurements)
  
  
  data <- rbind(data, data.frame(subject = rep(i, n_measurements),
                                 measurement = 1:n_measurements,
                                 Z1 = Z1,
                                 alpha=alpha,
                                 beta=beta))
  
}


functions_x_y <- generate_two_different_functions()
f_x <- functions_x_y$f
g_y <- functions_x_y$g

data$X <- NA
data$Y <- NA

for(i in 1:nrow(data)) {
  
  data$X[i] <- f_x(data$Z1[i])  + data$alpha[i] + 0.3*rnorm(1, mean = 0, sd = 1)
  data$Y[i] <- g_y(data$Z1[i])  + data$beta[i] + 0.3*rnorm(1, mean = 0, sd = 1) + b*data$X[i]
}


head(data)
#Z1 <- scale(data$Z1)
# Z2 <- scale(data$Z2)
# Z3 <- scale(data$Z3)

# Z <- cbind(Z1,Z2,Z3)
#X <- scale(data$X)
#Y <- scale(data$Y)

#data_scaled <- data.frame(X, Y, Z1)



gcm <- gcm.test(data$X, data$Y, data$Z1, regr.method = "gam", plot.residuals = F, regr.pars = list())
gcm

# res_gam_Y <- comp.resids(data$Y, data$Z1, regr.pars = list(), regr.method = "gam")
# res_gam_X <- comp.resids(data$X, data$Z1, regr.pars = list(), regr.method = "gam")
# plot(res_gam_Y)
# plot(res_gam_X)
#  

# data_scaled <- data.frame(X, Y, Z1, subject=data$subject)

gamm_Y <- gamm4(Y~ s(Z1), data = data, random = ~(1|subject))
gamm_X <- gamm4(X~ s(Z1), data = data, random = ~(1|subject))
res_gamm_Y <- resid(gamm_Y$mer)
res_gamm_X <- resid(gamm_X$mer)
plot(res_gamm_Y)
plot(res_gamm_X)
gcm_gamm <- gcm.test(resid.XonZ = res_gamm_X, resid.YonZ = res_gamm_Y)
gcm_gamm

# cor(res_gamm_X, res_gamm_Y)
# cor(res_gam_X, res_gam_Y)
# f_x
# g_y
# plot(Z1,X)
# plot(Z1,Y)

