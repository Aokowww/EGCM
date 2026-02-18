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


methodsvec <- c("gcm", "gcm+gamm")
nrep <- 100
nsubvec <- c(300, 600, 1000, 1200, 1500)
pvals <- array(NA, dim = c(length(methodsvec), nrep, length(nsubvec))) 
tstat <- array(NA, dim = c(length(methodsvec), nrep, length(nsubvec))) 
b <- 0.3

k <- 1; j <- 1

for(k in 1:length(nsubvec)){
  show(paste(k,':',nsubvec[k]))
    for(j in 1:nrep){
      show(j)
      set.seed(100*k + 10*j)
      n_subjects <- nsubvec[k]
      data <- data.frame(subject = integer(0),
                         measurement = integer(0),
                         Z1 = numeric(0),
                         Z2 = numeric(0),
                      #   Z3 = numeric(0),
                         alpha=numeric(0),
                         beta=numeric(0))
      
      
      alpha_vector <- rnorm(n_subjects, mean = 0, sd = 5)
      beta_vector <- rnorm(n_subjects, mean = 0, sd = 5)
      
      for (i in 1:n_subjects) {
        n_measurements <- round(rnorm(1, mean = 6.730864, sd = 4.746604))
        n_measurements <- max(min(n_measurements, 24), 1)
        cor_matrix_i <- matrix(0.5, nrow = n_measurements, ncol = n_measurements)
        diag(cor_matrix_i) <- 1
        Z1 <- mvrnorm(n = 1, mu = rep(0, n_measurements), Sigma = cor_matrix_i)
        Z2 <- mvrnorm(n = 1, mu = rep(0, n_measurements), Sigma = cor_matrix_i)
       # Z3 <- mvrnorm(n = 1, mu = rep(0, n_measurements), Sigma = cor_matrix_i)
        #Z <- cbind(Z1,Z2,Z3)
        
        alpha <- rep(alpha_vector[i], n_measurements)
        beta <- rep(beta_vector[i], n_measurements)
        
        
        data <- rbind(data, data.frame(subject = rep(i, n_measurements),
                                       measurement = 1:n_measurements,
                                       Z1 = Z1,
                                       Z2 = Z2,
                                    #   Z3 = Z3,
                                       alpha=alpha,
                                       beta=beta))
        
      }
      
      
      functions_x_y <- generate_two_different_functions()
      f_x <- functions_x_y$f
      g_y <- functions_x_y$g
      
      data$X <- NA
      data$Y <- NA
      
      for(i in 1:nrow(data)) {
        data$X[i] <- f_x(data$Z1[i]) + f_x(data$Z2[i]) + data$alpha[i] + 0.3*rnorm(1, mean = 0, sd = 1)
        data$Y[i] <- g_y(data$Z1[i]) + g_y(data$Z2[i]) + data$beta[i] + 0.3*rnorm(1, mean = 0, sd = 1) + b*data$X[i]
        #data$X[i] <- f_x(data$Z1[i]) + f_x(data$Z2[i]) + f_x(data$Z3[i]) + data$alpha[i] + 0.3*rnorm(1, mean = 0, sd = 1)
        #data$Y[i] <- g_y(data$Z1[i]) + g_y(data$Z2[i]) + g_y(data$Z3[i]) + data$beta[i] + 0.3*rnorm(1, mean = 0, sd = 1) + b*data$X[i]
      }
      
      # Z1 <- scale(data$Z1)
      # Z2 <- scale(data$Z2)
      # Z3 <- scale(data$Z3)
      # 
      # Z <- cbind(Z1,Z2,Z3)
      # X <- scale(data$X)
      # Y <- scale(data$Y)
      
      # data_scaled <- data.frame(X, Y, Z1, Z2, Z3, subject=data$subject)
      
      for(l in 1:length(methodsvec)){
        switch(methodsvec[l],
               "gcm" = {
                 gcm_original <- gcm.test(data$X, data$Y, cbind(data$Z1,data$Z2,data$Z3), regr.method = "gam", plot.residuals = F, regr.pars = list())
                 tstat[l,j,k] <- gcm_original$test.statistic
                 pvals[l,j,k] <- gcm_original$p.value
               },
               "gcm+gamm" = {
                 gamm_Y <- gamm4(Y~ s(Z1)+ s(Z2), data = data, random = ~(1|subject))
                 gamm_X <- gamm4(X~ s(Z1)+ s(Z2), data = data, random = ~(1|subject))
                # gamm_Y <- gamm4(Y~ s(Z1)+ s(Z2) + s(Z3), data = data, random = ~(1|subject))
                # gamm_X <- gamm4(X~ s(Z1)+ s(Z2) + s(Z3), data = data, random = ~(1|subject))
                 res_gamm_Y <- resid(gamm_Y$mer)
                 res_gamm_X <- resid(gamm_X$mer)
                 gcm_gamm <- gcm.test(resid.XonZ = res_gamm_X, resid.YonZ = res_gamm_Y)
                 tstat[l,j,k] <- gcm_gamm$test.statistic
                 pvals[l,j,k] <- gcm_gamm$p.value
                 
               })
      }
      show(paste(nsubvec[k],': ', sum(pvals[1,1:j,k] < 0.05)))
      show(paste(nsubvec[k],': ', sum(pvals[2,1:j,k] < 0.05)))
      
    }
}

save(pvals, file = "pvals_mimic_adni_real_mean_sd_zbiv_power.RData")
save(tstat, file = "tstat_mimic_adni_real_mean_sd_zbiv_power.RData")
#load("pvals_mimic_adni_real_mean_sd_zmulti.RData")

significance_level <- 0.05


rejectionRates <- array(NA, dim = c(length(methodsvec), length(nsubvec)))

for (m in 1:length(methodsvec)) {
    for (k in 1:length(nsubvec)) {
      
      pvals_subset <- pvals[m, , k]
      rejection_rate <- mean(pvals_subset < significance_level, na.rm = TRUE)
      rejectionRates[m, k] <- rejection_rate
    }
  }



print(rejectionRates)


rejectionRates_df <- data.frame(
  Method = character(),
  NSub = integer(),
  RejectionRate = numeric()
)


for (m in 1:length(methodsvec)) {
    for (k in 1:length(nsubvec)) {
      rejectionRates_df <- rbind(rejectionRates_df, data.frame(
        Method = methodsvec[m],
        NSub = nsubvec[k],
        RejectionRate = rejectionRates[m, k]
      ))
    }
  }



print(rejectionRates_df)

rejectionRates_df$NSub <- factor(rejectionRates_df$NSub)

# Create the line plot
ggplot(data = rejectionRates_df, aes(x = NSub, y = RejectionRate, group = Method, color = Method)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +  
  labs(title = "Type II Error for Mimic ADNI with multi-variate Z",
       x = "Number of Subject",
       y = "Type II Error") +
  theme_minimal() 


