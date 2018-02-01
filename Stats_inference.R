# 1. A random sample of size 6 from the exp(λ) 
# distribution results in observations: 1.636, 0.374, 0.534, 3.015, 0.932, 0.179. 
# Find the MLE on this data set in two ways: by numerical optimization of the likelihood and by the analytic formula.

data1 <- c(1.636,0.374,0.534,3.015,0.932,0.179)
# Analytic
nlik <- function(lam) -prod(lam*exp(-lam*data1))
optim(par = 1, nlik)$par
# Numerical 
nloglikEXP1 <- function(lam) -sum(log2(lam*exp(-lam*data1)))
optim(par = 1, nloglikEXP1)$par


# 2. 