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
    negexp = function(z) exp(-sqrt(z^2)),
    cosin = function(z) cos(z)
  )
  
  indices <- sample(length(functions), 2)
  f <- functions[[indices[1]]]
  g <- functions[[indices[2]]]
  
  return(list(f = f, g = g))
}


methodsvec <- c("GCM", "EGCM")
nrep <- 100
nsubvec <- c(50, 100, 200, 400, 600) 
nmeasvec <- c(3, 5, 10, 20, 25) 
pvals <- array(NA, dim = c(length(methodsvec), nrep, length(nmeasvec), length(nsubvec))) 
tstat <- array(NA, dim = c(length(methodsvec), nrep, length(nmeasvec), length(nsubvec))) 
b <- 0

t <- 1; k <- 1; j <- 1
for(t in 1:length(nmeasvec)){
  show(paste(t,':',nmeasvec[t]))
  for(k in 1:length(nsubvec)){
    show(paste(k,':',nsubvec[k]))
    for(j in 1:nrep){
      show(j)
      set.seed(100*k + 10*j + t)
      
      n_subjects <- nsubvec[k]
      n_measurements <- nmeasvec[t]
      
      # cor_matrix <- matrix(0.5, nrow = n_measurements, ncol = n_measurements)
      # diag(cor_matrix) <- 1
      unit_matrix <- diag(n_measurements)
      
      data <- data.frame(subject = rep(1:n_subjects, each = n_measurements),
                         measurement = rep(1:n_measurements, times = n_subjects),
                         Z = NA)
      
      # Generate Z values
      for (i in 1:n_subjects) {
        data$Z[data$subject == i] <- mvrnorm(n = 1, mu = rep(0, n_measurements), Sigma = unit_matrix)
      #  data$Z[data$subject == i] <- mvrnorm(n = 1, mu = rep(0, n_measurements), Sigma = cor_matrix)
      }
      
      data$X <- NA
      data$Y <- NA
    
      alpha <- rnorm(n_subjects,sd=5)
      alpha <- rep(alpha, each = n_measurements)
      
      beta <- rnorm(n_subjects,sd=5)
      beta <- rep(beta, each = n_measurements)
      
      functions_x_y <- generate_two_different_functions()
      f_x <- functions_x_y$f
      g_y <- functions_x_y$g
      
      for(i in 1:nrow(data)) {
        data$X[i] <- f_x(data$Z[i]) + alpha[i] + rnorm(1, mean = 0, sd = 1)  
        data$Y[i] <- g_y(data$Z[i]) + beta[i] + rnorm(1, mean = 0, sd = 1)  
      }
      
      
      for(l in 1:length(methodsvec)){
        switch(methodsvec[l],
               "GCM" = {
                 gcm_original <- gcm.test(data$X, data$Y, data$Z, regr.method = "gam", plot.residuals = FALSE, regr.pars = list())
                 pvals[l,j,t,k] <-gcm_original$p.value
                 tstat[l,j,t,k] <- gcm_original$test.statistic
               },
               "EGCM" = {
                 model_X_gamm <- gamm4(X ~ s(Z), data = data, random = ~(1|subject))
                 model_Y_gamm <- gamm4(Y ~ s(Z), data = data, random = ~(1|subject))
                 residuals_X_gamm <- resid(model_X_gamm$mer)
                 residuals_Y_gamm <- resid(model_Y_gamm$mer)
                 egcm <- gcm.test(resid.XonZ = residuals_X_gamm, resid.YonZ = residuals_Y_gamm)
                 pvals[l,j,t,k] <- egcm$p.value
                 tstat[l,j,t,k] <- egcm$test.statistic
               })
      }
      show(paste(nmeasvec[t],',',nsubvec[k],': ', sum(pvals[1,1:j,t,k] < 0.05)))
      show(paste(nmeasvec[t],',',nsubvec[k],': ', sum(pvals[2,1:j,t,k] < 0.05)))
      
    }
  }
}


# save(pvals, file = "pvals_Exp1_TpyI.RData")
# save(tstat, file = "tstat_Exp1_TpyI.RData")

load("pvals_post_nonlinear_nmeas_nsubj.RData")
# 
# save(pvals, file = "pvals_post_nonlinear_nmeas_nsubj.RData")
# 
# 
# load("pvals_post_nonlinear_nmeas_nsubj.RData")

significance_level <- 0.05


rejectionRates <- array(NA, dim = c(length(methodsvec), length(nmeasvec), length(nsubvec)))

for (m in 1:length(methodsvec)) {
  for (j in 1:length(nmeasvec)) {
    for (k in 1:length(nsubvec)) {

      pvals_subset <- pvals[m, , j, k]
      rejection_rate <- mean(pvals_subset < significance_level, na.rm = TRUE)
      rejectionRates[m, j, k] <- rejection_rate
    }
  }
}


print(rejectionRates)


rejectionRates_df <- data.frame(
  Method = character(),
  NSub = integer(),
  NMeas = integer(),
  RejectionRate = numeric()
)


for (m in 1:length(methodsvec)) {
  for (j in 1:length(nmeasvec)) {
    for (k in 1:length(nsubvec)) {
      rejectionRates_df <- rbind(rejectionRates_df, data.frame(
        Method = methodsvec[m],
        NSub = nsubvec[k],
        NMeas = nmeasvec[j],
        RejectionRate = rejectionRates[m, j, k]
      ))
    }
  }
}


print(rejectionRates_df)

# greens <- c("green", "darkgreen", "springgreen", "palegreen", "forestgreen", "olivedrab", "darkolivegreen", "limegreen", "seagreen", "mediumseagreen", "lightgreen", "darkseagreen", "yellowgreen", "lawngreen", "chartreuse", "lime", "greenyellow")
# yellows <- c("yellow", "gold", "goldenrod", "darkgoldenrod", "lightyellow", "lemonchiffon", "lightgoldenrodyellow", "papayawhip", "moccasin", "peachpuff", "palegoldenrod", "khaki", "darkkhaki")
# oranges = c("darkorange", "orange", "orangered", "tomato", "coral")
# reds = c("darkred", "red", "orangered", "firebrick")

my_colors <- c("darkgreen","forestgreen", "green3", "yellow","gold", "orange","darkorange", "red","darkred")
values <- c(0, 0.05, 0.1, 0.2, 0.3, 0.5 ,0.6 , 0.7, 0.8, 1) 


color_gradient <- colorRampPalette(my_colors)


p <- ggplot(rejectionRates_df, aes(x = factor(NSub), y = factor(NMeas), fill = RejectionRate)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colours = color_gradient(100), 
                       values = scales::rescale(values), 
                       breaks = c(0, 0.05, 0.3, 0.7, 1),
                       labels = c("0", "0.05", "0.3", "0.7", "1"),
                       guide = guide_colourbar(ticks = TRUE, nbin = 500)) +
  facet_wrap(~Method, ncol = 1) +
  labs(x = "Sample Size", y = "Number of Measurements", fill = "Rejection\nRate") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 10, face = "bold"),
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))


ggsave("Heatmaps_Exp1_Tpy1.pdf", plot = p, width = 12, height = 10)





nsub_labels <- function(n) {
  paste(n, "subjects")
}


p2 <- ggplot(rejectionRates_df, aes(x = NMeas, y = RejectionRate, color = Method, group = interaction(Method, NSub))) +
  geom_line(aes(linetype = Method), size = 1.2, alpha = 0.8) +  
  geom_point(aes(shape = Method), size = 4, alpha = 0.8) +  
  scale_color_manual(values = c("GCM" = "#377eb8", "EGCM" = "#e41a1c")) +  
  facet_wrap(~ NSub, scales = "free", ncol = 2, labeller = as_labeller(nsub_labels)) +  
  labs(title = "Type I Error, Exp. 1",
       subtitle = "by Number of Measurements across Different Numbers of Subjects",
       x = "Number of Measurements", y = "Type I Error Rate") +
  theme_minimal(base_size = 12) +  
  theme(legend.title = element_blank(),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "bottom",  
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "mm")) +  
  scale_x_continuous(breaks = c(3, 5, 10, 20, 50, 100)) 

ggsave("Linegraph_Exp1_Tpy1.pdf", plot = p2, width = 12, height = 10)



