library(ggplot2)
set.seed(100)
data <- data.frame(x = rnorm(100, 2, 1), y = rnorm(100, 1, 1))
data2 <- melt(data)
data3 <- lapply(data, function(x) get_summary_stats(data.frame(x)))
data3
# $x
# # A tibble: 1 x 13
# variable     n    min   max median    q1    q3   iqr   mad  mean    sd    se    ci
# <chr>    <dbl>  <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#   1 x          100 -0.272  4.58   1.94  1.39  2.66  1.26 0.974  2.00  1.02 0.102 0.203
# 
# $y
# # A tibble: 1 x 13
# variable     n   min   max median    q1    q3   iqr   mad  mean    sd    se    ci
# <chr>    <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#   1 x          100 -1.14  3.17  0.927 0.568  1.45 0.878 0.648  1.01 0.796  0.08 0.158


cor.test(data[,1], data[,2], method = "pearson")
# Pearson's product-moment correlation
# 
# data:  data[, 1] and data[, 2]
# t = -1.1205, df = 98, p-value = 0.2652
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.30221314  0.08584318
# sample estimates:
#        cor 
# -0.1124713 
cor.test(data[,1], data[,2], method = "spearman")
# Spearman's rank correlation rho
# 
# data:  data[, 1] and data[, 2]
# S = 192064, p-value = 0.1297
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.1524992 


ggplot(data, aes(x = x, y = y)) + 
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm") +
  theme_bw()