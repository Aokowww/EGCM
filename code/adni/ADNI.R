# install.packages("Hmisc")
install.packages("ADNIMERGE_0.0.1.tar.gz", repos = NULL, type = "source")
# help(package = "ADNIMERGE")
# install.packages("summarytools")
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

data(adnimerge)
selected_data <- adnimerge[, c("RID", "EXAMDATE", "DX", "Ventricles", "Hippocampus", 
                               "Entorhinal", "MOCA", "MMSE", "CDRSB")]

adnidd <- selected_data %>%
  # select('RID', 'EXAMDATE', 'DX', 'Ventricles', 'Hippocampus',
  #        'Entorhinal', 'MOCA', 'MMSE', 'CDRSB') %>%
  filter(!is.na(DX) & !is.na(Ventricles) & !is.na(Hippocampus) & 
           !is.na(Entorhinal) & !is.na(MOCA) & !is.na(MMSE) & !is.na(CDRSB))%>%
  mutate(DX_binary = ifelse(DX %in% c("CN", "MCI"), 0, 1))

adnidd <- adnidd %>%
  group_by(RID) %>%
  mutate(MEASURE = row_number())
write.csv(adnidd, "adnidd_data.csv", row.names = FALSE)
head(adnidd)

X <- cbind(adnidd$Ventricles, adnidd$Hippocampus, adnidd$Entorhinal)
X1 <- adnidd$Ventricles
X2 <- adnidd$Hippocampus
X3 <- adnidd$Entorhinal

           
Y <- cbind(adnidd$DX_binary)

MOCA <- adnidd$MOCA+rnorm(n = nrow(adnidd), mean = 0, sd = 1)
MMSE <- adnidd$MMSE+rnorm(n = nrow(adnidd), mean = 0, sd = 1)
CDRSB <- adnidd$CDRSB+rnorm(n = nrow(adnidd), mean = 0, sd = 1)

cor_matrix_Z <- cor(adnidd[, c("MOCA", "MMSE", "CDRSB")])
cor_matrix_Z

cor_matrix_X <- cor(adnidd[, c('Ventricles', 'Hippocampus','Entorhinal')])
cor_matrix_X

Z <- cbind(MOCA,MMSE,CDRSB)
Z1 <- MOCA


res_X1Z <- comp.resids(X1, Z, regr.pars = list(), regr.method = "gam")
plot(res_X1Z)
res_X2Z <- comp.resids(X2, Z, regr.pars = list(), regr.method = "gam")
res_X3Z <- comp.resids(X3, Z, regr.pars = list(), regr.method = "gam")
res_XallZ <- comp.resids(X, Z, regr.pars = list(), regr.method = "gam")
res_X_gam <- cbind(res_X1Z,res_X2Z,res_X3Z)
dim(res_X_gam)
plot(res_X_gam)
plot(res_XallZ)

res_YallZ <- comp.resids(Y, Z, regr.pars = list(), regr.method = "gam")
plot(res_YallZ)
gam_Y_Z <- gam(Y~s(Z), family = binomial)
res_YbinZ <- residuals(gam_Y_Z)
plot(res_YbinZ)

data_X_gam <- data.frame(
  index = 1:length(res_X_gam[,1]), # only the first col for plotting
  residuals_x = res_X_gam[,1], # only the first col for plotting
  subject = adnidd$RID
)


data_Y_gam_nobin <- data.frame(
  index = 1:length(res_YallZ),
  residuals_y = res_YallZ,
  subject = adnidd$RID
)

data_Y_gam_bin <- data.frame(
  index = 1:length(res_YbinZ),
  residuals_y = res_YbinZ,
  subject = adnidd$RID
)

p1 <- ggplot(data_X_gam, aes(x = index, y = residuals_x, colour = factor(subject))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ggtitle("Residuals for X1 by Sample Index") +
  theme_minimal() +
  labs(colour = "Subject")+
  theme(legend.position = "none") 
ggsave("adni_gam_res_x1.pdf", plot = p3, width = 12, height = 10)

p2_1 <- ggplot(data_Y_gam_nobin, aes(x = index, y = residuals_y, colour = factor(subject))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ggtitle("Residuals for Y by Sample Index from GAM") +
  theme_minimal() +
  labs(colour = "Subject")+
  theme(legend.position = "none") 

ggsave("adni_gam_res_y_nobin.pdf", plot = p4, width = 12, height = 10)

p2_2 <- ggplot(data_Y_gam_bin, aes(x = index, y = residuals_y, colour = factor(subject))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ggtitle("Residuals for Y by Sample Index from GAM with binomial family") +
  theme_minimal() +
  labs(colour = "Subject")+
  theme(legend.position = "none") 

ggsave("adni_gam_res_y_bin.pdf", plot = p4, width = 12, height = 10)


gcm_sep <- gcm.test(resid.XonZ = cbind(res_X1Z,res_X2Z,res_X3Z), resid.YonZ = res_YbinZ)
gcm_sep
# > gcm_sep
# $p.value
# [1] 0.002
# 
# $test.statistic
# [1] 5.264625
# 
# $reject
# [1] TRUE

gcm_tog <- gcm.test(X, Y, Z, regr.method = "gam", plot.residuals = F, regr.pars = list())
gcm_tog
# > gcm_tog
# $p.value
# [1] 0.002
# 
# $test.statistic
# [1] 5.622249
# 
# $reject
# [1] TRUE


gamm_X1Z <- gamm4(Ventricles ~ s(MOCA)+s(MMSE)+s(CDRSB), data = adnidd, random = ~(1|RID))
gamm_X2Z <- gamm4(Hippocampus ~ s(MOCA)+s(MMSE)+s(CDRSB), data = adnidd, random = ~(1|RID))
gamm_X3Z <- gamm4(Entorhinal ~ s(MOCA)+s(MMSE)+s(CDRSB), data = adnidd, random = ~(1|RID))

gamm_X1Z_re <- gamm4(Ventricles ~ s(MOCA)+s(MMSE)+s(CDRSB), data = adnidd, random = ~(EXAMDATE|RID))
gamm_X2Z_re <- gamm4(Hippocampus ~ s(MOCA)+s(MMSE)+s(CDRSB), data = adnidd, random = ~(EXAMDATE|RID))
gamm_X3Z_re <- gamm4(Entorhinal ~ s(MOCA)+s(MMSE)+s(CDRSB), data = adnidd, random = ~(EXAMDATE|RID))

res_X_gamm <- cbind(resid(gamm_X1Z$mer),resid(gamm_X2Z$mer),resid(gamm_X3Z$mer))
plot(res_X_gamm[,1])

res_X_gamm_re <- cbind(resid(gamm_X1Z_re$mer),resid(gamm_X2Z_re$mer),resid(gamm_X3Z_re$mer))
plot(res_X_gamm_re[,1])

gamm_YZ <- gamm4(DX_binary ~ s(MOCA)+s(MMSE)+s(CDRSB), data = adnidd, random = ~(1|RID))
res_Y_gamm <- resid(gamm_YZ$mer)
plot(res_Y_gamm)



gcm_gamm <- gcm.test(resid.XonZ = res_X_gamm, resid.YonZ = res_Y_gamm)
gcm_gamm

# > gcm_gamm
# $p.value
# [1] 0.002
# 
# $test.statistic
# [1] 5.47234
# 
# $reject
# [1] TRUE


data_X_gamm <- data.frame(
  index = 1:length(res_X_gamm[,1]), # only the first col for plotting
  residuals_x = res_X_gamm[,1], # only the first col for plotting
  subject = adnidd$RID
)


data_Y_gamm <- data.frame(
  index = 1:length(res_Y_gamm),
  residuals_y = res_Y_gamm,
  subject = adnidd$RID
)

p3 <- ggplot(data_X_gamm, aes(x = index, y = residuals_x, colour = factor(subject))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ggtitle("Residuals for X by Sample Index") +
  theme_minimal() +
  labs(colour = "Subject")+
  theme(legend.position = "none") 
ggsave("adni_gamm_res_x1.pdf", plot = p3, width = 12, height = 10)

p4 <- ggplot(data_Y_gamm, aes(x = index, y = residuals_y, colour = factor(subject))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ggtitle("Residuals for Y by Sample Index") +
  theme_minimal() +
  labs(colour = "Subject")+
  theme(legend.position = "none") 

ggsave("adni_gamm_res_y.pdf", plot = p4, width = 12, height = 10)








# tstat_gcm_gam1 <- gcm.test(resid.XonZ = resid(model_X1_gam$mer), resid.YonZ = residuals_Y_gam)$test.statistic
# tstat_gcm_gam1
# tstat_gcm_gam2 <- gcm.test(resid.XonZ = resid(model_X2_gam$mer), resid.YonZ = residuals_Y_gam)$test.statistic
# tstat_gcm_gam2
# tstat_gcm_gam3 <- gcm.test(resid.XonZ = resid(model_X3_gam$mer), resid.YonZ = residuals_Y_gam)$test.statistic
# tstat_gcm_gam3
# 
# test_statistic <- max(tstat_gcm_gam1,tstat_gcm_gam2,tstat_gcm_gam2)
# p_value <- 1 - pnorm(test_statistic)
# print(p_value)
# tstat_gcm_gam <- gcm.test(resid.XonZ = residuals_X_gam, resid.YonZ = residuals_Y_gam)$test.statistic
# tstat_gcm_gam
# 
# plot(residuals_Y_gam)
# gcm.test(resid.XonZ = residuals_X_gam, resid.YonZ = residuals_Y_gam)$p.value
# p.value <- 2 * pnorm(abs(tstat_gcm_gam), lower.tail = FALSE)
# p.value
# 
# pvals_gcm_gam1 <- gcm.test(resid.XonZ = resid(model_X1_gam$mer), resid.YonZ = residuals_Y_gam)$p.value
# 
# residuals_X_gam_1 <- resid(model_X_gam_1$mer)
# residuals_Y_gam <- resid(model_Y_gam$mer)
# plot(residuals_X_gam[,1])
# summary(adnidd$MEASURE)
# plot(adnidd$MEASURE)
# dfSummary(adnidd)
# 
# 
# 
# # tstat_gcm <- gcm.test(X, Y, Z, regr.method = "gam", plot.residuals = F, regr.pars = list())$test.statistic
# # tstat_gcm
# # 
# # pvals_gcm <- gcm.test(X, Y, Z, regr.method = "gam", plot.residuals = F, regr.pars = list())$p.value
# # pvals_gcm
# # pvals_gcm_Z1 <- gcm.test(X1, Y, Z1, regr.method = "gam", plot.residuals = F, regr.pars = list())$p.value
# # pvals_gcm_Z1
# # pvals_gcm_Z2 <- gcm.test(X, Y, Z2, regr.method = "gam", plot.residuals = F, regr.pars = list())$p.value
# # pvals_gcm_Z2
# # 
# # res <- comp.resids(X, Z3, regr.pars = list(), regr.method = "gam")
# # res <- comp.resids(X, Z, regr.pars = list(), regr.method = "gam")
# # res2 <- comp.resids(Y, Z, regr.pars = list(), regr.method = "gam")
# # res_Z1 <- comp.resids(X1, Z1, regr.pars = list(), regr.method = "gam")
# # res_Z2 <- comp.resids(X, Z2, regr.pars = list(), regr.method = "gam")
# 
# # 
# # plot(res_Z1)
# # plot(res_Z2)
# # length(unique(adnidd$DX_binary))
# # 
# # length(unique(adnidd$MOCA))
# # length(unique(adnidd$MMSE))
# # length(unique(adnidd$CDRSB))
# 
# # library(gamm4)
# # dat <- gamSim(1,n=600,dist="binary",scale=.33)
# # lr.fit0 <- gam(y~s(x0)+s(x1)+s(x2), family = binomial, data = dat)
# # residuals(lr.fit0)
# # plot(residuals(lr.fit0))
# # lr.fit <- gamm4(y~s(x0)+s(x1)+s(x2), family = binomial, data = dat)
# # resid(lr.fit$mer)
# # plot(resid(lr.fit$mer))
# 
# 
# 
