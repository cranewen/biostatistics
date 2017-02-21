#Question 1
library(multtest)
data("golub")
gol.fac <- factor(golub.cl, levels=0:1, labels = c("ALL","AML"))
dataH4j <- golub[2972,]
dataAPS <- golub[2989,]
dataH4j_ALL <- golub[2972, gol.fac=="ALL"]
dataH4j_AML <- golub[2972, gol.fac=="AML"]
dataAPS_ALL <- golub[2989, gol.fac=="ALL"]
dataAPS_AML <- golub[2989, gol.fac=="AML"]

#(a)
t.test(dataH4j_ALL, alternative = "greater", mu = -0.9)
#reject null hypothesis

#(b)
t.test(dataH4j_ALL, dataAPS_AML, var.equal = FALSE)
#reject null hypothesis

#(c)
t.test(dataH4j_ALL-dataAPS_ALL, alternative = "less", mu = 0)
#reject null hypothesis

#(d)
xd <- sum((dataH4j_ALL-dataAPS_ALL)<0)
binom.test(x=xd,n=27,p=0.5,alternative = "greater")
#Accept null hypothesis, Plow = 0.6296296, so we agree with the conclusion in part (c)

#(e)
xe <- sum(dataH4j_ALL>-0.6)
binom.test(x=xe,n=27,p=0.5,alternative = "less")
#Accept null hypothesis, equals to 0.5

#(f)
x_all <- sum(dataH4j_ALL>-0.6)
x_aml <- sum(dataH4j_AML>-0.6)
prop.test(x=c(x_all,x_aml),n=c(27,11),alternative = "two.sided",correct = F)
#Reject null hypothesis, they are not equal.

#Question 2
#(a)
ex_reject <- 2000*0.05
#(b)
binom.test(x=ex_reject,n=2000,p=90/2000,alternative = "less")
#The probability is 5%

#Question 3
#(a)
nsim3 <- 10000
p.value3 <- rep(NA, nsim3)
p.bon3 <- rep(NA, nsim3)
p.fdr3 <- rep(NA, nsim3)
for(i in 1:nsim3){
  data3 <- rnorm(20,3,4)
  p.value3[i] <- t.test(data3, alternative = "greater", mu = 3, conf.level = 0.9)$p.value
  p.bon3[i] <- p.adjust(p=p.value3[i], method = "bonferroni")
  p.fdr3[i] <- p.adjust(p=p.value3[i], method = "fdr")
}
errorRate3 <- sum(p.value3<0.1)/nsim3
#The error rate is closed to a=0.1 level, which is 10000*0.1=1000
#(b)
#No, we shouldn't. Because the high rate of error.

#Question 4
#(a)
p.values4 <- rep(NA, 3051)
for(i in 1:3051){
  p.values4[i] <- t.test(golub[i,]~gol.fac)$p.value
}
p.bon4 <- p.adjust(p=p.values4, method = "bonferroni")
p.fdr4 <- p.adjust(p=p.values4, method = "fdr")
cat("From Bonferron, there are: ", sum(p.bon4<0.05))
cat("From FDR, there are: ", sum(p.fdr4<0.05))
#(b)
index4b <- which(p.fdr4<0.05)
sortedData4b <- order(p.values4[index4b], decreasing = FALSE)
cat("The top three strongeet differentially expressed genes are: ")
for(i in 1:3){
  cat(golub.gnames[sortedData4b[1:3],2])
}

    
