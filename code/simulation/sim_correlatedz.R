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


methodsvec <- c("GCM","EGCM")
nrep <- 100
nsubvec <- c(20, 50, 100, 200) 
nmeasvec <- c(3, 5, 10, 20, 50, 100) 
pvals <- array(NA, dim = c(length(methodsvec), nrep, length(nmeasvec), length(nsubvec))) 
b <- 0.3

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
      
      cor_matrix <- matrix(0.5, nrow = n_measurements, ncol = n_measurements)
      diag(cor_matrix) <- 1
      
      data <- data.frame(subject = rep(1:n_subjects, each = n_measurements),
                         measurement = rep(1:n_measurements, times = n_subjects),
                         Z = NA)
      
      # Generate Z values
      for (i in 1:n_subjects) {
        data$Z[data$subject == i] <- mvrnorm(n = 1, mu = rep(0, n_measurements), Sigma = cor_matrix)
      }
      
      data$X <- NA
      data$Y <- NA
      
      # alpha <- rnorm(n_subjects,sd=5)
      # alpha <- rep(alpha, each = n_measurements)
      # 
      # beta <- rnorm(n_subjects,sd=5)
      # beta <- rep(beta, each = n_measurements)
      # 
      functions_x_y <- generate_two_different_functions()
      f_x <- functions_x_y$f
      g_y <- functions_x_y$g
      
      for(i in 1:nrow(data)) {
        data$X[i] <- f_x(data$Z[i]) + rnorm(1, mean = 0, sd = 1)  
        data$Y[i] <- g_y(data$Z[i]) + rnorm(1, mean = 0, sd = 1)  
        # data$X[i] <- f_x(data$Z[i]) + alpha[i] + rnorm(1, mean = 0, sd = 1)  
        # data$Y[i] <- g_y(data$Z[i]) + beta[i] + rnorm(1, mean = 0, sd = 1)    
      }
      #pvals[1,j,t,k] <- gcm.test(data$X, data$Y, data$Z, regr.method = "gam", plot.residuals = FALSE, regr.pars = list())$p.value
      # 
      for(l in 1:length(methodsvec)){
        switch(methodsvec[l],
               "gcm" = {
                 pvals[l,j,t,k] <- gcm.test(data$X, data$Y, data$Z, regr.method = "gam", plot.residuals = FALSE, regr.pars = list())$p.value
               },
               "gcm+gamm" = {
                 model_X_gamm <- gamm4(X ~ s(Z), data = data, random = ~(1|subject))
                 model_Y_gamm <- gamm4(Y ~ s(Z), data = data, random = ~(1|subject))
                 residuals_X_gamm <- resid(model_X_gamm$mer)
                 residuals_Y_gamm <- resid(model_Y_gamm$mer)

                 pvals[l,j,t,k] <- gcm.test(resid.XonZ = residuals_X_gamm, resid.YonZ = residuals_Y_gamm)$p.value
               })
      }
      show(paste(nmeasvec[t],',',nsubvec[k],': ', sum(pvals[1,1:j,t,k] < 0.05)))
      show(paste(nmeasvec[t],',',nsubvec[k],': ', sum(pvals[2,1:j,t,k] < 0.05)))
      
    }
  }
}


#save(pvals, file = "pvals_gcm_works_power.RData")
load("pvals_gcm_works_power.RData")

#save(pvals, file = "pvals_gcm_works.RData")


#load("pvals_gcm_works.RData")
# load("pvals_gcm_works.RData")
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
        RejectionRate = 1-rejectionRates[m, j, k]
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
  geom_text(aes(label = sprintf("%.2f", RejectionRate)), vjust = 1) +
  scale_fill_gradientn(colours = color_gradient(100), 
                       values = scales::rescale(values), 
                       breaks = c(0, 0.05, seq(0.25, 1, by = 0.25)),  # Include special mark at 0.05
                       labels = c("", "0.05", seq(0.25, 1, by = 0.25)),  # Adjust labels here
                       # breaks = c(0, 0.05, 0.1, 0.2, 0.3, 0.5 ,0.6 , 0.7, 0.8, 1),
                       # labels = c("0", "0.05", "0.1", "0.2", "0.3", "0.5" ,"0.6" , "0.7", "0.8", "1"),
                       guide = guide_colourbar(ticks = TRUE, nbin = 1000),
                       limits = c(0, 1)) +  # Set limits from 0 to 1
  facet_wrap(~Method, ncol = 1) +
  labs(x = "Number of subjects", y = "Number of measurements", fill = "Power") +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
    legend.box.margin = margin(t = 0, r = -20, b = 0, l = 0, unit = "pt"),
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(0.5, "cm"),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"), 
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    axis.title = element_text(face = "bold"),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )


p




# ggsave("Heatmaps_correlatedz.pdf", plot = p, width = 12, height = 10)


# nsub_labels <- function(n) {
#   paste(n, "subjects")
# }
# 
# 
# p2 <- ggplot(rejectionRates_df, aes(x = NMeas, y = RejectionRate, color = Method, group = interaction(Method, NSub))) +
#   geom_line(aes(linetype = Method), size = 1.2, alpha = 0.8) +  
#   geom_point(aes(shape = Method), size = 4, alpha = 0.8) +  
#   scale_color_manual(values = c("gcm" = "#377eb8", "gcm+gamm" = "#e41a1c")) +  
#   facet_wrap(~ NSub, scales = "free", ncol = 2, labeller = as_labeller(nsub_labels)) +  
#   labs(title = "Type II Error, Case I",
#        subtitle = "by Number of Measurements across Different Numbers of Subjects",
#        x = "Number of Measurements", y = "1-Type II Error Rate") +
#   theme_minimal(base_size = 12) +  
#   theme(legend.title = element_blank(),
#         plot.title = element_text(size = 16, face = "bold"),
#         plot.subtitle = element_text(size = 14),
#         axis.title = element_text(size = 14),
#         axis.text = element_text(size = 12),
#         legend.text = element_text(size = 12),
#         strip.text = element_text(size = 14, face = "bold"),
#         legend.position = "bottom",  
#         plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "mm")) +  
#   scale_x_continuous(breaks = c(3, 5, 10, 20, 50, 100)) 
# p2
# 
# ggsave("Linegraph_correlatedz_power.pdf", plot = p2, width = 12, height = 10)
# 
# library(ggplot2)
# library(scales)
# 
# 

#
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






