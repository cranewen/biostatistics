#Question 1
#(b)
data1 <- c(1.636,0.374,0.534,3.015,0.932,0.179)
#analytic
lik <- function(lam) prod(lam*exp(-lam*data1))
optim(par = 1, lik)$par
#numerical
nloglikEXP <- function(lam) -(6*log(lam)-lam*sum(c(1.636,0.374,0.534,3.015,0.932,0.179)))
optim(par = 1, nloglikEXP)$par


#Question 2
#(a)
# m = Xbar = 100.8

#(b)
100.8 + qt(c(0.1), df = 53-1)*12.4/sqrt(53)
# So one-sided 90% lower CI of m is (98.58908, inf)

#Question 3
data(golub, package = "multtest")
Zyxin.row <- grep("Zyxin",golub.gnames[,2])
gol.fac <- factor(golub.cl, levels=0:1, labels = c("ALL","AML"))
Zyxin.ALL <- golub[Zyxin.row, gol.fac=="ALL"]
Zyxin.AML <- golub[Zyxin.row, gol.fac=="AML"]

#(a)
nboot <- 1000
nALL <- length(Zyxin.ALL)
nAML <- length(Zyxin.AML)
boot.zyxinALLmean <- rep(NA, nboot)
boot.zyxinALLvar <- rep(NA, nboot)
boot.zyxinAMLmean <- rep(NA, nboot)
boot.zyxinAMLvar <- rep(NA, nboot)
for(i in 1:nboot){
  data.starALL <- sample(Zyxin.ALL, replace = TRUE)
  boot.zyxinALLmean[i] <- mean(data.starALL)
  boot.zyxinALLvar[i] <- var(data.starALL)
  data.starAML <- sample(Zyxin.AML, replace = TRUE)
  boot.zyxinAMLmean[i] <- mean(data.starAML)
  boot.zyxinAMLvar[i] <- var(data.starAML)
}

#Zyxin.All 95% CI mean is:
quantile(boot.zyxinALLmean,c(0.025,0.975))
#Zyxin.All 95% CI variance is:
quantile(boot.zyxinALLvar,c(0.025,0.975))


#Zyxin.AML 95% CI mean is:
quantile(boot.zyxinAMLmean,c(0.025,0.975))
#Zyxin.AMl 95% CI variance is:
quantile(boot.zyxinAMLvar,c(0.025,0.975))


#(b)
#Zyxin.All 95% CI mean is:
mean(Zyxin.ALL)+qt(c(0.025,0.975),df=nALL-1)*sd(Zyxin.ALL)/sqrt(nALL)
#Zyxin.All 95% CI variance is:
sqrt(nALL-1)*sd(Zyxin.ALL)/sqrt(qchisq(c(0.975,0.025),df=nALL-1))
#Zyxin.AML 95% CI mean is:
mean(Zyxin.AML)+qt(c(0.025,0.975),df=nAML-1)*sd(Zyxin.AML)/sqrt(nAML)
#Zyxin.AMl 95% CI variance is:
sqrt(nAML-1)*sd(Zyxin.AML)/sqrt(qchisq(c(0.975,0.025),df=nAML-1))

#(c)
boot.zyxinALLmedian <- rep(NA, nboot)
boot.zyxinAMLmedian <- rep(NA, nboot)
for(i in 1:nboot){
  data.starALL <- sample(Zyxin.ALL, replace = TRUE)
  boot.zyxinALLmedian[i] <- median(data.starALL)
  data.starAML <- sample(Zyxin.AML, replace = TRUE)
  boot.zyxinAMLmedian[i] <- median(data.starAML)
}
#Zyxin ALL 95% CI is:
quantile(boot.zyxinALLmedian, c(0.025,0.975))
#Zyxin AML 95% CI is:
quantile(boot.zyxinAMLmedian, c(0.025,0.975))

#(d)
#Yes, gene expressions in ALL and AML patients are different.

#Question 4
#(a)
nsim <- 1000
data4mean <- rep(NA, nsim)
data4var <- rep(NA, nsim)
fx <- function(lam){ 
for(i in 1:nsim){
  data <- rpois(50, lambda = lam)
  data4mean[i] <<- mean(data)
  data4var[i] <<- var(data)
  }
}

#(b)
fx(0.1)
quantile(data4mean, c(0.05,0.95))
quantile(data4var, c(0.05,0.95))
fx(1)
quantile(data4mean, c(0.05,0.95))
quantile(data4var, c(0.05,0.95))
fx(10)
quantile(data4mean, c(0.05,0.95))
quantile(data4var, c(0.05,0.95))

#(c)
fx(10)
mean(data4mean)+qt(c(0.05,0.95),49)*sqrt(mean(data4mean)/50)
mean(data4var)*49/qchisq(c(0.95,0.05),49)
