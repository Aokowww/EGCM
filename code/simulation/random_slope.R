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
# Define vectors for subject numbers and time points
ni_vector <- c(50,100,200,300,400)

nj <- 10  # Number of time points

# Initialize array to store p-values
pvals <- array(NA, dim = c(length(ni_vector), 100, 2))  # 2 is the length of ni_vector

# Initialize rejection counts
rejection_count_resid <- rep(0, length(ni_vector))
rejection_count_gamm <- rep(0, length(ni_vector))

# Loop over different numbers of subjects
for (idx in 1:length(ni_vector)) {
  ni <- ni_vector[idx]
  
  for (i in 1:100) {
    # Generate random effects
    alpha_vector <- rnorm(ni, mean = 0, sd = 5)
    beta_vector <- rnorm(ni, mean = 0, sd = 5)
   # RE_x <- mvrnorm(ni, mu = c(0,0), Sigma = rbind(c(1.0, 0.8), c(0.8, 1.0)))
   # RE_y <- mvrnorm(ni, mu = c(0,0), Sigma = rbind(c(1.0, 0.8), c(0.8, 1.0)))
    RE_x <- rnorm(ni, mean = 0, sd = 0.5)
    RE_y <- rnorm(ni, mean = 0, sd = 0.5)
    # colnames(RE_x) <- c("ints", "slopes")
    # colnames(RE_y) <- c("ints", "slopes")
    
    functions_x_y <- generate_two_different_functions()
    f_x <- functions_x_y$f
    g_y <- functions_x_y$g
    
    # Create data frame
    data <- data.frame(ID = rep(1:ni, each = nj), 
                       time = rep(1:nj, times = ni),
                       RE.ix = rep(alpha_vector, each = nj),
                       RE.sx = rep(RE_x, each = nj),
                       RE.iy = rep(beta_vector, each = nj),
                       RE.sy = rep(RE_y, each = nj),
                       z = NA,
                       x = NA,
                       y = NA)
    
    # Generate random numbers
    z <- with(data, rnorm(n = ni * nj, mean = 0, sd = 1))
    data$z <- z
    
    # Generate X and Y variables
    x <- with(data, (3*RE.ix) + f_x(z) + (3 + RE.sx) * z + rnorm(n = ni * nj, mean = 0, sd = 1))
    
    y <- with(data, (3*RE.iy) + g_y(z) + (3 + RE.sy) * z + 0.2*f_x(z) + rnorm(n = ni * nj, mean = 0, sd = 1))
    
    # Introduce missing values
    m1 <- rbinom(n = ni * nj, size = 1, prob = 0.2)  
    y[m1 == 1] <- NA
    x[m1 == 1] <- NA
    z[m1 == 1] <- NA
    
    data$x <- x
    data$y <- y
    data$z <- z
    
    # Calculate residuals
    res_x <- comp.resids(data$x, data$z, regr.pars = list(), regr.method = "gam")
    res_y <- comp.resids(data$y, data$z, regr.pars = list(), regr.method = "gam")
    
    # Perform GCM test
    pval_gcm <- gcm.test(resid.XonZ = res_x, resid.YonZ = res_y)$p.value
    cat("Iteration", i, "GCM test p-value for ni =", ni, ":", pval_gcm, "\n")
    
    # Print p-value and update rejection count
    if (pval_gcm < 0.05) {
      cat("Iteration", i, "GCM test rejected for ni =", ni, "\n")
      rejection_count_resid[idx] <- rejection_count_resid[idx] + 1
      cat("Iteration", i, "GCM test rejected for ni =", ni,":", rejection_count_resid[idx], "\n")
    }
    
    # Save p-value for residuals
    pvals[idx, i, 1] <- pval_gcm
    
    # Perform GAMM fitting
    gamm_re_y <- gamm4(y ~ s(z), data = data, random = ~(1| ID))$mer
    gamm_re_x <- gamm4(x ~ s(z), data = data, random = ~(1| ID))$mer
    
    # Calculate residuals for random effects
    res_re_y <- resid(gamm_re_y)
    res_re_x <- resid(gamm_re_x)
    
    # Perform GCM test for GAMM residuals
    pval_egcm <- gcm.test(resid.XonZ = res_re_x, resid.YonZ = res_re_y)$p.value
    cat("Iteration", i, "EGCM test p-value for ni =", ni, ":", pval_egcm, "\n")
    
    # Print p-value and update rejection count
    if (pval_egcm < 0.05) {
      rejection_count_gamm[idx] <- rejection_count_gamm[idx] + 1
      cat("Iteration", i, "EGCM test rejected for ni =", ni,":", rejection_count_gamm[idx], "\n")
    }
    
    # Save p-value for GAMM residuals
    pvals[idx, i, 2] <- pval_egcm
  }
}

#save(pvals, file = "pvals_random_slope_alpha_beta_missp_sd1.5_power.RData")
load ("pvals_random_slope_alpha_beta_power_missp.RData")

head(pvals,10)
# # Calculate rejection rates
# rejection_rate_resid <- rejection_count_resid / (100 * length(ni_vector))
# rejection_rate_gamm <- rejection_count_gamm / (100 * length(ni_vector))
# 
# # Print rejection rates
# cat("GCM rejection rate:", rejection_rate_resid, "\n")
# cat("EGCM rejection rate:", rejection_rate_gamm, "\n")

# Calculate the proportion of p-values less than 0.05 for each combination of ni_vector
proportion_pvals <- array(NA, dim = c(length(ni_vector), 2), 
                          dimnames = list(ni_vector, c("GCM", "EGCM")))  # Initialize array to store proportions

# Loop through each combination
for (i in 1:length(ni_vector)) {
  for (j in 1:2) {
    pvals_subset <- pvals[i, , j]  # Extract subset for current combination
    proportion_pvals[i, j] <- 1-mean(pvals_subset < 0.05)  # Calculate proportion
  }
}

# Output the result
proportion_pvals





