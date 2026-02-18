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
library(ADNIMERGE)
# 
data(adnimerge)
# # Load necessary libraries
# library(ggplot2)  # For plotting
# 
# # Assuming adnimerge is your data frame and RID is the column of interest
# # Let's calculate the statistics first
# 
# # Count occurrences of each RID
# RID_counts <- table(adnimerge$RID)
# 
# # Calculate statistics
# max_count <- max(RID_counts)
# min_count <- min(RID_counts)
# avg_count <- mean(RID_counts)
# var_count <- var(RID_counts)
# 
# # Print the statistics
# cat("Maximum count of RID:", max_count, "\n")
# cat("Minimum count of RID:", min_count, "\n")
# cat("Average count of RID:", avg_count, "\n")
# cat("Variance of count of RID:", var_count, "\n")
# 
# # Now let's plot the frequency distribution
# # Convert RID_counts to a data frame for plotting
# RID_counts_df <- data.frame(RID = names(RID_counts), Count = as.numeric(RID_counts))
# 
# # Plot
# ggplot(RID_counts_df, aes(x = Count)) +
#   geom_histogram(binwidth = 1, fill = "skyblue", color = "black", aes(y = ..count..)) +
#   labs(title = "Frequency of RID Occurrences", x = "Number of Occurrences", y = "Frequency") +
#   theme_minimal()
# mean_count <- avg_count
# sd_count <- sqrt(var_count)

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


methodsvec <- c("GCM", "RGCM")
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
                       alpha=numeric(0),
                       beta=numeric(0))
    
    
    alpha_vector <- rnorm(n_subjects, mean = 0, sd = 5)
    beta_vector <- rnorm(n_subjects, mean = 0, sd = 5)
    RE_x <- mvrnorm(ni, mu = c(0,0), Sigma = rbind(c(1.0, 0.5), c(0.5, 1.0)))
    RE_y <- mvrnorm(ni, mu = c(0,0), Sigma = rbind(c(1.0, 0.5), c(0.5, 1.0)))
    colnames(RE_x) <- c("ints", "slopes")
    colnames(RE_y) <- c("ints", "slopes")
    
    for (i in 1:n_subjects) {
      n_measurements <- round(rnorm(1, mean =  6.730864, sd = 4.746604))
      n_measurements <- max(min(n_measurements, 24), 1)
      cor_matrix_i <- matrix(0.5, nrow = n_measurements, ncol = n_measurements)
      diag(cor_matrix_i) <- 1
      Z1 <- mvrnorm(n = 1, mu = rep(0, n_measurements), Sigma = cor_matrix_i)
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
      
      data$X[i] <- f_x(data$Z1[i]) + data$alpha[i] + 0.3*rnorm(1, mean = 0, sd = 1)
      data$Y[i] <- g_y(data$Z1[i]) + data$beta[i] + 0.3*rnorm(1, mean = 0, sd = 1) + b*data$X[i]
    }
    
    

    
   # X <- scale(data$X)
    #Y <- scale(data$Y)
    
 #   data_scaled <- data.frame(X, Y, Z1, subject=data$subject)
    
    for(l in 1:length(methodsvec)){
      switch(methodsvec[l],
             "gcm" = {
               
               pvals[l,j,k] <- gcm.test(data$X, data$Y, data$Z1, regr.method = "gam", plot.residuals = F, regr.pars = list())$p.value
             },
             "gcm+gamm" = {
               gamm_Y <- gamm4(Y~ s(Z1), data = data, random = ~(1|subject))
               gamm_X <- gamm4(X~ s(Z1), data = data, random = ~(1|subject))
               res_gamm_Y <- resid(gamm_Y$mer)
               res_gamm_X <- resid(gamm_X$mer)
               
               pvals[l,j,k] <- gcm.test(resid.XonZ = res_gamm_X, resid.YonZ = res_gamm_Y)$p.value
             })
    }
    show(paste(nsubvec[k],': ', sum(pvals[1,1:j,k] < 0.05)))
    show(paste(nsubvec[k],': ', sum(pvals[2,1:j,k] < 0.05)))
    
  }
}

save(pvals, file = "pvals_mimic_adni_zuni_4_real_mean_sd_power.RData")

load(file = "pvals_mimic_adni_zuni_4_real_mean_sd_power.RData")
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
  labs(title = "Type II Error for Mimic ADNI with univariate Z",
       x = "Number of Subject",
       y = "1-Type II Error") +
  theme_minimal() 


