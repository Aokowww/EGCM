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


methodsvec <- c("GCM", "EGCM")
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

#save(pvals, file = "pvals_mimic_adni_zuni_4_real_mean_sd_power.RData")

#load(file = "pvals_mimic_adni_real_mean_sd_zbiv.RData")
significance_level <- 0.05

load(file = "pvals_mimic_adni_zuni_4_real_mean_sd.RData")

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
# 
# # Create the line plot
# ggplot(data = rejectionRates_df, aes(x = NSub, y = RejectionRate, group = Method, color = Method)) +
#   geom_line() +
#   geom_point() +
#   geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +  
#   labs(title = "Type II Error for Mimic ADNI with univariate Z",
#        x = "Number of Subject",
#        y = "1-Type II Error") +
#   theme_minimal() 



base_size <- 12
title_size <- 16
subtitle_size <- 14
axis_title_size <- 14
axis_text_size <- 12
legend_text_size <- 12
strip_text_size <- 14


p <- ggplot(data = rejectionRates_df, aes(x = NSub, y = RejectionRate, group = Method, color = Method)) +
  geom_line(aes(linetype = Method), size = 1, alpha = 0.6) +
  geom_point(aes(shape = Method), size = 4, alpha = 0.6) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 0.8, alpha = 0.8) +
  scale_color_manual(values = c("GCM" = "#377eb8", "EGCM" = "#e7882c")) +
  scale_shape_manual(values = c("GCM" = 16, "EGCM" = 17)) +
  scale_y_continuous(
    breaks = c(0, 0.05, seq(0.1, 1, by = 0.1)), 
    labels = c("0", "0.05", seq(0.1, 1, by = 0.1)),
    limits = c(0, 1)
  ) +
  labs(
    title = "Type I Error in Mimicking ADNI with bivariate Z",
    x = "Number of Subjects",
    y = "Type I Error"
  ) +
  theme_minimal(base_size = base_size) +
  theme(
    plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
    axis.title = element_text(size = axis_title_size),
    axis.text = element_text(size = axis_text_size),
    legend.text = element_text(size = legend_text_size),
    strip.text = element_text(size = strip_text_size, face = "bold"),
    legend.position = "bottom",
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "mm"),
    panel.grid.major.y = element_line(color = "grey80", linetype = "solid"),
    panel.grid.minor.y = element_blank()
  )

print(p)


base_size <- 12
title_size <- 16
subtitle_size <- 14
axis_title_size <- 14
axis_text_size <- 12
legend_text_size <- 12
strip_text_size <- 14


x_breaks <- pretty(rejectionRates_df$NMeas, n = 5)

p2 <- ggplot(rejectionRates_df, aes(x = NMeas, y = RejectionRate, group = interaction(Method, NSub))) +
  geom_line(aes(color = Method, linetype = Method), size = 1, alpha = 0.6) +  
  geom_point(aes(color = Method, shape = Method), size = 4, alpha = 0.6) +  
  geom_hline(yintercept = 0.95, linetype = "dashed", color ="red", size = 0.8, alpha = 0.8) + 
  scale_color_manual(values = c("GCM" = "#377eb8", "EGCM" = "#e7882c")) +  
  facet_wrap(~ NSub, scales = "free", ncol = 2, labeller = as_labeller(nsub_labels)) +  
  labs(title = "Power results of the synthetic simulation b)",
       subtitle = "by different number of measurements across different numbers of subjects",
       x = "Number of Measurements", y = "Power") +
  theme_minimal(base_size = 12) +  
  theme(
    legend.title = element_blank(),
    plot.title = element_text(size = title_size, face = "bold", hjust = 0.5), 
    plot.subtitle = element_text(size = subtitle_size, hjust = 0.5), 
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "mm")
  ) +  
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(breaks = x_breaks) 

p2

