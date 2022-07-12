if (!requireNamespace("survminer", quietly = TRUE))
  install.packages("survminer")
library(survival)
library(survminer)

# data <- lung
# colnames(data)[5] <- "variable"


fit <- survfit(Surv(time, status) ~ variable, data = data)
print(fit)
# Call: survfit(formula = Surv(time, status) ~ variable, data = data)
# 
# n events median 0.95LCL 0.95UCL
# variable=1 138    112    270     212     310
# variable=2  90     53    426     348     550


survdiff(Surv(time, status) ~ variable, data = data)
# survdiff(formula = Surv(time, status) ~ variable, data = data)
# 
# N Observed Expected (O-E)^2/E (O-E)^2/V
# variable=1 138      112     91.6      4.55      10.3
# variable=2  90       53     73.4      5.68      10.3
# 
# Chisq= 10.3  on 1 degrees of freedom, p= 0.001 


fit2 <- coxph(Surv(time, status) ~ variable, data = data)
summary(fit2)
# Call:
#   coxph(formula = Surv(time, status) ~ variable, data = data)
# 
# n= 228, number of events= 165 
# 
# coef exp(coef) se(coef)      z Pr(>|z|)   
# variable -0.5310    0.5880   0.1672 -3.176  0.00149 **
#   ---
#   Signif. codes:  0 鈥�***鈥� 0.001 鈥�**鈥� 0.01 鈥�*鈥� 0.05 鈥�.鈥� 0.1 鈥� 鈥� 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# variable     0.588      1.701    0.4237     0.816
# 
# Concordance= 0.579  (se = 0.021 )
# Likelihood ratio test= 10.63  on 1 df,   p=0.001
# Wald test            = 10.09  on 1 df,   p=0.001
# Score (logrank) test = 10.33  on 1 df,   p=0.001


# plot
ggsurvplot(fit = fit, data = data, pval = T)