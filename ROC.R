library(tidyverse)
library(ggplot2)
library(timeROC)

data <- dat
data <- read.table("~/file.txt", header = T)

fit_timeROC <-
  timeROC(`T` = data$time,
          delta = data$event,
          marker = data$group,
          cause = 1,  weighting="marginal",
          times = c(1*365, 3*365, 5*365), iid = T) 
fit_timeROC
# Time-dependent-Roc curve estimated using IPCW  (n=418, without competing risks). 
#        Cases Survivors Censored AUC (%)   se
# t=365     30       388        0   82.30 3.64
# t=1095    81       312       25   83.42 2.48
# t=1825   115       197      106   86.23 2.15
# 
# Method used for estimating IPCW:marginal 
# 
# Total computation time : 0.24  secs.

confint(fit_timeROC)
# $CI_AUC
#         2.5% 97.5%
# t=365  75.16 89.44
# t=1095 78.55 88.28
# t=1825 82.02 90.44
# 
# $CB_AUC
#         2.5% 97.5%
# t=365  73.97 90.63
# t=1095 77.74 89.09
# t=1825 81.32 91.14
# 
# $C.alpha
#      95% 
# 2.286114 


cut_value = sort(unique(data$group))
tp = rev(fit_timeROC$TP[,1])
fp = rev(fit_timeROC$FP[,1])
idx <- which.max(tp - fp)[1]
SeSpPPVNPV(cutpoint=cut_value[idx - 1],
           `T` = data$time,
           delta = data$event,
           marker = data$group,
           cause= 1,  
           weighting="marginal",
           times = c(1*365))
# Predictive accuracy measures at cutpoint c=2.3 estimated using IPCW (n=418, with competing risks). 
# No. of positive (X>c) =143, No. of negative (X<=c) =275. 
# 
#       Cases Survivors Censored Se (%) Sp (%) PPV (%) NPV (%)
# t=0       0       418        0     NA  65.79      NA  100.00
# t=365    30       388        0  86.67  69.85   18.18   98.55
# 
# Method used for estimating IPCW:marginal 
# 
# Total computation time : 0  secs.


data2 <- data.frame(group = "1-Year",
                    x = fit_timeROC$FP[,1],
                    y = fit_timeROC$TP[,1])
data2 <- data2[order(data2$x, data2$y),]
tmp <- data.frame(group = "3-Year",
                  x = fit_timeROC$FP[,2],
                  y = fit_timeROC$TP[,2])
tmp <- tmp[order(tmp$x, tmp$y),]
data2 <- rbind(data2, tmp)

ggplot() +
  geom_line(data = data2, aes(x = x, y = y, colour = group)) +
  labs(x = "1-Specificity (FPR)", y = "Sensitivity (TPR)")
