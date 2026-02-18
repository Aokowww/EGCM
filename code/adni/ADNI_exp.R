# install.packages("Hmisc")
#install.packages("ADNIMERGE_0.0.1.tar.gz", repos = NULL, type = "source")
# help(package = "ADNIMERGE")
# install.packages("summarytools")
library(summarytools)
library(ADNIMERGE)
library(dplyr)
library(GeneralisedCovarianceMeasure) # Generalised Covariance Measure (GCM)
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
library(lmtest)
library(ggrepel)
library(patchwork)
library(RColorBrewer)
library(GGally)
library(gridExtra)

data(adnimerge)
selected_data <- adnimerge[, c("RID", "EXAMDATE", "DX", "Ventricles", "Hippocampus", 
                               "Entorhinal", "MOCA", "MMSE", "CDRSB", "AV45")]

adnidd <- selected_data %>%
  na.omit()%>%
  mutate(DX_binary = ifelse(DX %in% c("CN", "MCI"), 0, 1))

adnidd <- adnidd %>%
  group_by(RID) %>%
  mutate(MEASURE = row_number())


head(adnidd)
# write.csv(adnidd, "adnidd_data.csv", row.names = FALSE)


#RID_tab <- table(adnidd$RID)
#id_sel <- names(RID_tab[RID_tab > 5])
#adnidd_sel <- adnidd %>% filter(RID %in% id_sel)
#adnidd <- adnidd_sel

X <- cbind(adnidd$Ventricles, adnidd$Hippocampus, adnidd$Entorhinal)
X <- scale(X)
X1 <- X[,1]
X2 <- X[,2]
X3 <- X[,3]

DX <- adnidd$DX
Y <- scale(adnidd$AV45)


cor_matrix_X
#cor_matrix_Z <- cor(adnidd[, c("MOCA", "MMSE", "CDRSB")])
#cor_matrix_Z

# cor_matrix_X <- cor(adnidd[, c('Ventricles', 'Hippocampus','Entorhinal')])
# cor_matrix_X

MOCA <- adnidd$MOCA+rnorm(n = nrow(adnidd), mean = 0, sd = 1)
MMSE <- adnidd$MMSE+rnorm(n = nrow(adnidd), mean = 0, sd = 1)
CDRSB <- adnidd$CDRSB+rnorm(n = nrow(adnidd), mean = 0, sd = 1)

correlation <- cor(dummy_DX, Y_df)

X_df <- as.data.frame(cbind(adnidd$Ventricles, adnidd$Hippocampus, adnidd$Entorhinal))
Z_df <- as.data.frame(cbind(adnidd$MOCA, adnidd$MMSE, adnidd$CDRSB))
Y_df <- as.data.frame(adnidd$AV45)

DX_df <- as.data.frame(DX)
dummy_DX <- model.matrix(~DX - 1, data = DX_df) # -1 removes intercept
cor_DX_Y <- cor(DX_df, Y_df, use = "pairwise.complete.obs")

colnames(X_df) <- c("Ventricles", "Hippocampus", "Entorhinal")
colnames(Z_df) <- c("MOCA", "MMSE", "CDRSB")
colnames(Y_df) <- c("AV45")
# cor_XZ <- cor(X_df, Z_df, use = "pairwise.complete.obs")
# labels_XZ <- c(colnames(X_df), colnames(Z_df))
# cor_YZ <- cor(Y_df, Z_df, use = "pairwise.complete.obs")
# labels_YZ <- c(colnames(Y_df), colnames(Z_df))
# cor_YX <- cor(Y_df, X_df, use = "pairwise.complete.obs")
# labels_YX <- c(colnames(Y_df), colnames(X_df))
# 
# pdf("Correlation_XYZ.pdf", width = 8, height = 5)
# layout_matrix <- matrix(c(1, 2, 
#                           1, 3), 
#                         byrow = TRUE, 
#                         nrow = 2)
# layout(layout_matrix)
# corrplot(cor_XZ, method = "color", addCoef.col = NA, tl.col = "black", 
#          tl.srt = 45, tl.cex = 0.8, main = "Correlation between X and Z", 
#          cl.pos = "n")
# 
# corrplot(cor_YZ, method = "color", addCoef.col = NA, tl.col = "black", 
#          tl.srt = 45, tl.cex = 0.8, main = "Correlation between Y and Z", 
#          cl.pos = "n")
# corrplot(cor_YX, method = "color", addCoef.col = NA, tl.col = "black", 
#          tl.srt = 45, tl.cex = 0.8, main = "Correlation between Y and X", 
#          cl.pos = "n")
# dev.off()
# 
# 

# theme_set(theme_bw(base_size = 14) + theme(legend.position = "bottom"))
# 
# 
# py <- ggplot(adnidd, aes(x = DX, y = AV45, fill = DX)) + 
#   geom_boxplot(outlier.shape = NA) +  
#   geom_jitter(width = 0.2, alpha = 0.1, size = 1.2, color = "black") +
#   scale_fill_brewer(palette = "Paired") + 
#   labs(x = "Diagnosis", y = "AV45 ratio (cortical grey matter/whole cerebellum)", title = "AV45 Across Diagnoses") +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     axis.title.x = element_text(face = "bold"),
#     axis.title.y = element_text(face = "bold"),
#     legend.title = element_blank(),
#     legend.background = element_rect(fill = "white"),
#     legend.key = element_rect(fill = "white", colour = "white")
#   )
# print(py)
# 
# ggsave("plot_DX_AV45.pdf", py, width = 14, height = 7)
# 
# px1 <- ggplot(adnidd, aes(x = DX, y = Ventricles, fill = DX)) + 
#   geom_boxplot(outlier.shape = NA) +  
#   geom_jitter(width = 0.2, alpha = 0.1, size = 1, color = "black") +
#   scale_fill_brewer(palette = "Paired") + 
#   labs(x = "Diagnosis", y = "Ventricles", title = "Ventricles Across Diagnoses") +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     axis.title.x = element_text(face = "bold"),
#     axis.title.y = element_text(face = "bold"),
#     legend.title = element_blank(),
#     legend.background = element_rect(fill = "white"),
#     legend.key = element_rect(fill = "white", colour = "white")
#   )
# 
# px2 <- ggplot(adnidd, aes(x = DX, y = Hippocampus, fill = DX)) + 
#   geom_boxplot(outlier.shape = NA) +  
#   geom_jitter(width = 0.2, alpha = 0.1, size = 1, color = "black") +
#   scale_fill_brewer(palette = "Paired") + 
#   labs(x = "Diagnosis", y = "Hippocampus", title = "Hippocampus Across Diagnoses") +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     axis.title.x = element_text(face = "bold"),
#     axis.title.y = element_text(face = "bold"),
#     legend.title = element_blank(),
#     legend.background = element_rect(fill = "white"),
#     legend.key = element_rect(fill = "white", colour = "white")
#   )
# px3 <- ggplot(adnidd, aes(x = DX, y = Entorhinal, fill = DX)) + 
#   geom_boxplot(outlier.shape = NA) +  
#   geom_jitter(width = 0.2, alpha = 0.1, size = 1, color = "black") +
#   scale_fill_brewer(palette = "Paired") + 
#   labs(x = "Diagnosis", y = "Entorhinal", title = "Entorhinal Across Diagnoses") +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     axis.title.x = element_text(face = "bold"),
#     axis.title.y = element_text(face = "bold"),
#     legend.title = element_blank(),
#     legend.background = element_rect(fill = "white"),
#     legend.key = element_rect(fill = "white", colour = "white")
#   )
# 
# # Combine plots
# plot_X <- px1 + px2 +px3 + plot_layout(ncol = 3)
# 
# # Display combined plot
# plot_X
# 
# ggsave("plot_DX_X.pdf", plot_X, width = 14, height = 7)
# 
# 
# 
# pz1 <- ggplot(adnidd, aes(x = DX, y = adnidd$MOCA, fill = DX)) + 
#   geom_boxplot(outlier.shape = NA) +  
#   geom_jitter(width = 0.2, alpha = 0.1, size = 1, color = "black") +
#   scale_fill_brewer(palette = "Paired") + 
#   labs(x = "Diagnosis", y = "MOCA", title = "MOCA Across Diagnoses") +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     axis.title.x = element_text(face = "bold"),
#     axis.title.y = element_text(face = "bold"),
#     legend.title = element_blank(),
#     legend.background = element_rect(fill = "white"),
#     legend.key = element_rect(fill = "white", colour = "white")
#   )
# 
# pz2 <- ggplot(adnidd, aes(x = DX, y = adnidd$MMSE, fill = DX)) + 
#   geom_boxplot(outlier.shape = NA) +  
#   geom_jitter(width = 0.2, alpha = 0.1, size = 1, color = "black") +
#   scale_fill_brewer(palette = "Paired") + 
#   labs(x = "Diagnosis", y = "MMSE", title = "MMSE Across Diagnoses") +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     axis.title.x = element_text(face = "bold"),
#     axis.title.y = element_text(face = "bold"),
#     legend.title = element_blank(),
#     legend.background = element_rect(fill = "white"),
#     legend.key = element_rect(fill = "white", colour = "white")
#   )
# pz3 <- ggplot(adnidd, aes(x = DX, y = adnidd$CDRSB, fill = DX)) + 
#   geom_boxplot(outlier.shape = NA) +  
#   geom_jitter(width = 0.2, alpha = 0.1, size = 1, color = "black") +
#   scale_fill_brewer(palette = "Paired") + 
#   labs(x = "Diagnosis", y = "CDRSB", title = "CDRSB Across Diagnoses") +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     axis.title.x = element_text(face = "bold"),
#     axis.title.y = element_text(face = "bold"),
#     legend.title = element_blank(),
#     legend.background = element_rect(fill = "white"),
#     legend.key = element_rect(fill = "white", colour = "white")
#   )
# 
# # Combine plots
# plot_Z <- pz1 + pz2 +pz3 + plot_layout(ncol = 3)
# 
# # Display combined plot
# plot_Z
# 
# # ggsave("plot_DX_Z.pdf", plot_Z, width = 14, height = 7)
# 
# 
res_X1Z <- comp.resids(X1, Z, regr.pars = list(), regr.method = "gam")
#
res_X2Z <- comp.resids(X2, Z, regr.pars = list(), regr.method = "gam") 
res_X3Z <- comp.resids(X3, Z, regr.pars = list(), regr.method = "gam")
#
res_XallZ <- comp.resids(X, Z, regr.pars = list(), regr.method = "gam") #
res_X_gam <- cbind(res_X1Z,res_X2Z,res_X3Z) # dim(res_X_gam) # plot(res_X_gam)
# # plot(res_XallZ)
#
#
res_YallZ <- comp.resids(Y, Z, regr.pars = list(), regr.method = "gam") #
# plot(res_YallZ)


# together
gcm <- gcm.test(X, Y, Z, regr.method = "gam", plot.residuals = F, regr.pars = list())
gcm
# gcm
# $p.value
# [1] 0.002
# 
# $test.statistic
# [1] 5.996177
# 
# $reject
# [1] TRUE


# pvalue1 <- gcm.test(resid.XonZ = res_X1Z, resid.YonZ = res_YallZ)
# pvalue2 <- gcm.test(resid.XonZ = res_X2Z, resid.YonZ = res_YallZ)
# pvalue3 <- gcm.test(resid.XonZ = res_X3Z, resid.YonZ = res_YallZ)

# FDR control


# seperate
gcm1 <- gcm.test(resid.XonZ = res_X1Z, resid.YonZ = res_YallZ)
gcm2 <- gcm.test(resid.XonZ = res_X2Z, resid.YonZ = res_YallZ)
gcm3 <- gcm.test(resid.XonZ = res_X3Z, resid.YonZ = res_YallZ)
gcm1
gcm2
gcm3

p_gcm_fdr <- p.adjust(c(0.01759646, 2.020165e-09, 0.0008367561))
p_gcm_fdr
# > gcm1
# $p.value
# [1] 0.01759646
# 
# $test.statistic
# [1] 2.374002
# 
# $reject
# [1] TRUE
# 
# > gcm2
# $p.value
# [1] 2.020165e-09
# 
# $test.statistic
# [1] -5.996177
# 
# $reject
# [1] TRUE
# 
# > gcm3
# $p.value
# [1] 0.0008367561
# 
# $test.statistic
# [1] -3.340341
# 
# $reject
# [1] TRUE

# > p_gcm_fdr
# [1] 1.759646e-02 6.060495e-09 1.673512e-03

# GAMM

adnidd_scale1 <- scale(adnidd[,c("MOCA","MMSE", "CDRSB", "Ventricles", "Hippocampus",
                                 "Entorhinal", "AV45")]) %>% as.data.frame()
adnidd_scale1$RID <- adnidd$RID
adnidd_scale1$EXAMDATE <- adnidd$EXAMDATE
gamm_X1Z <- gamm4(Ventricles ~ s(MMSE)+ s(CDRSB) + s(MOCA), data = adnidd_scale1, random = ~(1|RID))
gamm_X2Z <- gamm4(Hippocampus ~ s(MMSE)+ s(CDRSB) + s(MOCA), data = adnidd_scale1, random = ~(1|RID))
gamm_X3Z <- gamm4(Entorhinal ~ s(MMSE)+s(CDRSB) + s(MOCA), data = adnidd_scale1, random = ~(1|RID))
head(adnidd_scale1, 5)

adnidd_scale2 <- adnidd_scale1 %>% 
  arrange(RID, EXAMDATE)
adnidd_scale2 <- adnidd_scale2 %>% 
  group_by(RID) %>% 
  mutate(t = row_number())


# head(adnidd_scale2, 5)
# 
# gamm_X1Z_re <- gamm4(Ventricles ~ s(MOCA)+s(MMSE)+s(CDRSB), data = adnidd_scale2, random = ~ (1|RID))
# gamm_X2Z_re <- gamm4(Hippocampus ~ s(MOCA)+s(MMSE)+s(CDRSB), data = adnidd_scale2, random = ~(1|RID))
# gamm_X3Z_re <- gamm4(Entorhinal ~ s(MOCA)+s(MMSE)+s(CDRSB), data = adnidd_scale2, random = ~(1|RID))


res_X_gamm <- cbind(resid(gamm_X1Z$mer),resid(gamm_X2Z$mer),resid(gamm_X3Z$mer))

# res_X_gamm_re <- cbind(resid(gamm_X1Z_re$mer),resid(gamm_X2Z_re$mer),resid(gamm_X3Z_re$mer))
# plot(res_X_gamm_re[,1])

# correlation 
cor(res_X_gamm[,1], res_Y_gamm)
cor.test(res_X_gamm[,1], res_Y_gamm)
# total randomly, gcm ok
gcm.test(resid.XonZ = rnorm(2063), resid.YonZ = rnorm(2063))
# $p.value
# [1] 0.8853513
# 
# $test.statistic
# [1] 0.144189
# 
# $reject
# [1] FALSE

# with even small value of cor, gcm reject h0 because sample is too large
matrix(c(1,0.1,0.1,1),nrow=2)
dat <- mvrnorm(n = 2063, mu = c(0, 1), Sigma = matrix(c(1, 0.1, 0.1, 1), nrow = 2))
gcm.test(resid.XonZ = dat[,1], resid.YonZ =  dat[,2])

# $p.value
# [1] 0.001123913
# 
# $test.statistic
# [1] 3.257518
# 
# $reject
# [1] TRUE


# plot(res_X_gamm[,1])

gamm_YZ <- gamm4(AV45 ~ s(MMSE)+ s(CDRSB) + s(MOCA), data = adnidd_scale1, random = ~(1|RID))
res_Y_gamm <- resid(gamm_YZ$mer)
# plot(res_Y_gamm)


# gamm_YZ_re <- gamm4(AV45 ~ s(MOCA)+s(MMSE)+s(CDRSB), data = adnidd_scale1, random = ~(EXAMDATE|RID))
# res_Y_gamm_re <- resid(gamm_YZ_re$mer)


gcm_gamm <- gcm.test(resid.XonZ = res_X_gamm, resid.YonZ = res_Y_gamm)
gcm_gamm

# > gcm_gamm
# $p.value
# [1] 0.002
# 
# $test.statistic
# [1] 12.15996
# 
# $reject
# [1] TRUE

# gcm_gamm_re <- gcm.test(resid.XonZ = res_X_gamm_re, resid.YonZ = res_Y_gamm_re)
# gcm_gamm_re 


gcm_gamm1 <- gcm.test(resid.XonZ = res_X_gamm[,1], resid.YonZ = res_Y_gamm)
gcm_gamm1
gcm_gamm2 <- gcm.test(resid.XonZ = res_X_gamm[,2], resid.YonZ = res_Y_gamm)
gcm_gamm2
gcm_gamm3 <- gcm.test(resid.XonZ = res_X_gamm[,3], resid.YonZ = res_Y_gamm)
gcm_gamm3

p_egcm_fdr <- p.adjust(c(5.078436e-34, 6.831176e-27, 0.8987019))
p_egcm_fdr
# > p_egcm_fdr
# [1] 1.523531e-33 1.366235e-26 8.987019e-01
# > gcm_gamm1
# $p.value
# [1] 5.078436e-34
# 
# $test.statistic
# [1] 12.15996
# 
# $reject
# [1] TRUE
# 
# > gcm_gamm2
# $p.value
# [1] 6.831176e-27
# 
# $test.statistic
# [1] -10.73688
# 
# $reject
# [1] TRUE
# 
# > gcm_gamm3
# $p.value
# [1] 0.8987019
# 
# $test.statistic
# [1] 0.1273014
# 
# $reject
# [1] FALSE


# analyse
# clustered outliers obviously
plot(res_X1Z)

# Entorhinal with some correlation
plot(res_X3Z)
# residuals is more iid like normal case residuals
plot(res_X_gamm[,3])

# X1 reject with relative smaller test statistic, think conditional dependence, with slightly dependence
plot(res_X1Z, res_YallZ)
# X1 with GAMM, dependence structure is more clear
plot(res_X_gamm[,1], res_Y_gamm)
# X2 reject with relative obviously larger test statistic
plot(res_X2Z, res_YallZ)

# X3 reject with relative smaller test statistic, think conditional dependence, with slightly dependence

# X3 with GAMM, dependence structure is more clear
plot(res_X_gamm[,3], res_Y_gamm)

# pvalue1 <- gcm.test(resid.XonZ = res_X1Z, resid.YonZ = res_YallZ)
# pvalue2 <- gcm.test(resid.XonZ = res_X2Z, resid.YonZ = res_YallZ)
# pvalue3 <- gcm.test(resid.XonZ = res_X3Z, resid.YonZ = res_YallZ)
# 
# # FDR control
# p <- p.adjust(c(pvalue1$p.value, pvalue2$p.value, pvalue3$p.value))


library(ggplot2)
library(ggrepel)
library(dplyr)
library(RColorBrewer)




# Plot for res_X1Z from GAM
df_X1 <- data.frame(index = 1:nrow(adnidd), ID = unlist(adnidd$RID), res = res_X1Z)

highlighted_ids <- c("4424", "4714", "2045", "6321", "2363", "4334", "337", "4736", "5234", "4331", "751", "2045", "1346", "4094", "2155", "4415")
colors <- rainbow(length(highlighted_ids))

p1 <- ggplot(df_X1, aes(x = index, y = res)) + 
  geom_point(alpha = 0.1) + 
  geom_point(data = df_X1 %>% filter(ID %in% highlighted_ids), 
             aes(x = index, y = res, color = factor(ID)), alpha = 0.8) + 
  scale_color_manual(values = setNames(colors, highlighted_ids)) + 
  geom_text_repel(data = df_X1 %>% filter(ID %in% highlighted_ids), 
                  aes(x = index, y = res, label = ID)) + 
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.7) +
  labs(title = "Residuals Plot for X1 utilizing the Original GCM",
       x = "Index",
       y = "Residual of X1 from Z based on GAM") +  
  theme(legend.position = "none", 
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))
p1

# Plot for res_X_gamm[,1] from GAMM
df_gamm_X1 <- data.frame(index = 1:nrow(adnidd), ID = unlist(adnidd$RID), res = res_X_gamm[,1])

p2 <- ggplot(df_gamm_X1, aes(x = index, y = res)) + 
  geom_point(alpha = 0.1) + 
  geom_point(data = df_gamm_X1 %>% filter(ID %in% highlighted_ids), 
             aes(x = index, y = res, color = factor(ID)), alpha = 0.8) + 
  scale_color_manual(values = setNames(colors, highlighted_ids)) + 
  geom_text_repel(data = df_gamm_X1 %>% filter(ID %in% highlighted_ids), 
                  aes(x = index, y = res, label = ID)) + 
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.7) +
  labs(title = "Residuals Plot for X1 utilizing the extended GCM with GAMM",
       x = "Index",
       y = "Residual of X1 from Z based on GAMM") +
  theme(legend.position = "none", 
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))


# Combine plots
res_plot_X1 <- p1 + p2 + plot_layout(ncol = 2)

# Display combined plot
res_plot_X1

ggsave("res_plot_X1.pdf", res_plot_X1, width = 14, height = 7)


# plot for res_X3Z from GAM
df_X3 <- data.frame(index = 1:nrow(adnidd),ID = unlist(adnidd$RID), res = res_X3Z)

df_summary_X3 <- df_X3 %>% group_by(ID) %>% summarise(mean = mean(res), sd = sd(res), n = n())

highlighted_ids <- c("4222", "4277","4636", "4173", "4654", "4713","2047","6038","4356","6073","6810","4974")
colors <- rainbow(length(highlighted_ids))

p1 <- ggplot(df_X3, aes(x = index, y = res)) + 
  geom_point(alpha = 0.1) + 
  geom_point(data = df_X3 %>% filter(ID %in% highlighted_ids), 
             aes(x = index, y = res, color = factor(ID)), alpha = 0.8) + 
  scale_color_manual(values = setNames(colors, highlighted_ids)) + 
  geom_text_repel(data = df_X3 %>% filter(ID %in% highlighted_ids), 
                  aes(x = index, y = res, label = ID)) + 
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.7) +
  labs(title = "Residuals Plot for X3 utilizing the Original GCM",
       x = "Index",
       y = "Residual of X3 from Z based on GAM") +
  theme(legend.position = "none", 
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))

p1

############

library(ggplot2)
library(dplyr)
library(ggrepel)
library(viridis)



colors <- viridis::viridis(length(highlighted_ids), option = "D")

p1 <- ggplot(df_X3, aes(x = index, y = res)) + 
  geom_point(alpha = 0.1) + 
  geom_point(data = df_X3 %>% filter(ID %in% highlighted_ids), 
             aes(x = index, y = res, color = factor(ID)), alpha = 0.8, size = 5) + 
  scale_color_manual(values = setNames(colors, highlighted_ids)) + 
  geom_text_repel(data = df_X3 %>% filter(ID %in% highlighted_ids), 
                  aes(x = index, y = res, label = ID), 
                  size = 6, 
                  box.padding = 0.35, 
                  point.padding = 0.5,
                  max.overlaps = Inf,
                  nudge_x = 0.5,
                  nudge_y = 0.5) +
  geom_hline(yintercept = 0, color = "#D03060", linetype = "dashed", size = 0.7) +
  labs(title = "Residuals Plot for X3 from the GCM Analysis",
       x = "Index",
       y = "Residual of X3 from Z based on GAM") +
  theme_minimal() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))


print(p1)


#############
# Plot for res_X_gamm[,3] from GAMM
df_gamm_X3 <- data.frame(index = 1:nrow(adnidd), ID = unlist(adnidd$RID), res = res_X_gamm[,3])

library(ggplot2)
library(dplyr)
library(ggrepel)
library(viridis)



p2 <- ggplot(df_gamm_X3, aes(x = index, y = res)) + 
  geom_point(alpha = 0.1) + 
  geom_point(data = df_gamm_X3 %>% filter(ID %in% highlighted_ids), 
             aes(x = index, y = res, color = factor(ID)), alpha = 0.8, size = 5) +
  scale_color_manual(values = setNames(colors, highlighted_ids)) + 
  geom_text_repel(data = df_gamm_X3 %>% filter(ID %in% highlighted_ids), 
                  aes(x = index, y = res, label = ID), 
                  size = 6,
                  max.overlaps = Inf,
                  nudge_x = 0, nudge_y = 0, 
                  point.padding = 0.3, # 
                  box.padding = 0.35) +
  geom_hline(yintercept = 0, color = "#D03060", linetype = "dashed", size = 0.7) +
  labs(title = "Residuals Plot for X3 from EGCM Analysis",
       x = "Index",
       y = "Residual of X3 from Z based on GAMM") +
  theme_minimal() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

print(p2)



p2

p2 <- ggplot(df_gamm_X3, aes(x = index, y = res)) + 
  geom_point(alpha = 0.1) + 
  geom_point(data = df_gamm_X3 %>% filter(ID %in% highlighted_ids), 
             aes(x = index, y = res, color = factor(ID)), alpha = 0.8) + 
  scale_color_manual(values = setNames(colors, highlighted_ids)) + 
  geom_text_repel(data = df_gamm_X3 %>% filter(ID %in% highlighted_ids), 
                  aes(x = index, y = res, label = ID)) + 
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.7) +
  labs(title = "Residuals Plot for X3 utilizing the extended GCM with GAMM",
       x = "Index",
       y = "Residual of X3 from Z based on GAMM") +
  theme(legend.position = "none", 
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))

# Combine plots
res_plot_X3 <- p1 + p2 + plot_layout(ncol = 2)

# Display combined plot
res_plot_X3

ggsave("res_plot_X3.pdf", res_plot_X3, width = 14, height = 7)



# plot for res_YallZ from GAM
df_YallZ <- data.frame(index = 1:nrow(adnidd),ID = unlist(adnidd$RID), res = res_YallZ)

df_summary_YallZ <- df_YallZ %>% group_by(ID) %>% summarise(mean = mean(res), sd = sd(res), n = n())

highlighted_ids <- c("731", "4386", "4431","5290","5277", "5127", "4777", "5166","4271","4007","4172")
colors <- rainbow(length(highlighted_ids))

p1 <- ggplot(df_YallZ, aes(x = index, y = res)) + 
  geom_point(alpha = 0.1) + 
  geom_point(data = df_YallZ %>% filter(ID %in% highlighted_ids), 
             aes(x = index, y = res, color = factor(ID)), alpha = 0.8) + 
  scale_color_manual(values = setNames(colors, highlighted_ids)) + 
  geom_text_repel(data = df_YallZ %>% filter(ID %in% highlighted_ids), 
                  aes(x = index, y = res, label = ID)) + 
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.7) +
  labs(title = "Residuals Plot for Y utilizing the Original GCM",
       x = "Index",
       y = "Residual of Y from Z based on GAM") +
  theme(legend.position = "none", 
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))

# plot for res_Y_gamm from GAMM (compared with plot for GAM for res_YallZ)
df_Y_gamm <- data.frame(index = 1:nrow(adnidd),ID = unlist(adnidd$RID), res = res_Y_gamm)

p2 <- ggplot(df_Y_gamm, aes(x = index, y = res)) + 
  geom_point(alpha = 0.1) + 
  geom_point(data = df_Y_gamm %>% filter(ID %in% highlighted_ids), 
             aes(x = index, y = res, color = factor(ID)), alpha = 0.9) + 
  scale_color_manual(values = setNames(colors, highlighted_ids)) + 
  geom_text_repel(data = df_Y_gamm %>% filter(ID %in% highlighted_ids), 
                  aes(x = index, y = res, label = ID)) + 
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.7) +
  labs(title = "Residuals Plot for Y utilizing the extended GCM with GAMM",
       x = "Index",
       y = "Residual of Y from Z based on GAMM") +
  theme(legend.position = "none", 
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))

# Combine plots
res_plot_Y <- p1 + p2 + plot_layout(ncol = 2)

# Display combined plot
res_plot_Y

ggsave("res_plot_Y.pdf", res_plot_Y, width = 14, height = 7)



