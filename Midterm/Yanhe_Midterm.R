#Question 1
#(a)
fx1 <- function(x){2.469862*x*exp(-x^2)}
fy1 <- function(y){2*y*exp(-y^2)}
x1 <- c(1,2,3)
Ex1 <- sum(x1*fx1(x1))
Ey1 <- integrate(function(y) y*fy1(y),lower = 0, upper = Inf)$value
Sdx1 <- sqrt(sum((x1-Ex1)^2*fx1(x1)))
Sdy1 <- sqrt(integrate(function(y) (y-Ey1)^2*fy1(y) , lower = 0, upper = Inf)$value)
#(b)


#Question 2
X <- rnorm(1000, mean = 0, sd = 1)
Y <- rchisq(1000, df = 4)
mean(X^2/(X^2+Y))


#Question 3
nsim <- 1000
data3 <- rep(NA, nsim)
for(i in 1:nsim){
  data3[i] <- rnorm(1000, mean = 0, sd = 1)
}


#Question 4
y <- as.numeric(t(read.table(file = "normalData.txt", header = T)))


#Question 5
library(multtest)
data("golub")

#(a)
p.value5 <- t.test(golub[201,], alternative = "greater", mu = 0.6, conf.level = 0.9)$p.value
p.values <- apply(golub[1:3025,], 1, function(x) t.test(x, alternative = "greater", mu = 0.6, conf.level = 0.9)$p.value)
p.fdrs <- p.adjust(p.values, method = "fdr")
sum(p.fdrs>0.1)
#(b)
cat("Top five genes with mean expression values greater than 0.6: \n")
for(i in 1:5){
  cat(golub.gnames[order(p.fdrs,decreasing = FALSE)[i],2],"\n")
}

#Question 6
GRO3 <- golub[2715,]
MYC <- golub[2302,]
gol.fac <- factor(golub.cl, levels=0:1, labels = c("ALL","AML"))
#(a)
hist(GRO3)
#(b)
GRO3_ALL <- golub[2715,gol.fac=="ALL"]
GRO3_AML <- golub[2715,gol.fac=="AML"]
MYC_ALL <- golub[2302,gol.fac=="ALL"]
MYC_AML <- golub[2302,gol.fac=="AML"]
plot(GRO3_ALL, MYC_ALL, col = "green", pch = "+", xlab = "GRO3", ylab = "MYC", ylim = range(-1.6,1.6), xlim = range(-2,1))
points(GRO3_AML, MYC_AML, col = "blue")
legend("topright", c("ALL","AML"), col = c("green","blue"), pch = c(3,1))
#(c)
t.test(GRO3-MYC,alternative = "less")
#From GRO3 - MYC, we get p vlaue is 0.03718, so reject null hypothesis, GRO3 < MYC
#Here t.test(GRO3,MYC,alternative = "less", paired = TRUE) is equivalent to t.test(GRO3-MYC,alternative = "less")
#(d)
diff_GRO3_MYC <- GRO3 - MYC
t_statistic <- (mean(diff_GRO3_MYC)-0)/(sd(diff_GRO3_MYC)/sqrt(37))
1-pnorm(t_statistic)
#matches the t-test

#(e)
wilcox.test(GRO3, MYC, paired = T, alternative = "less")
#p-value is greater than 0.01, so we reject null hypothesis

#(f)






#Question 7
#(a)
HPCA_row <- grep("HPCA Hippocalcin", golub.gnames[,2])
cat("The row number of HPCA Hippocalcin is: ", HPCA_row)
#(b)
HPCA <- golub[HPCA_row,]
HPCA_ALL <- golub[HPCA_row,gol.fac=="ALL"]
HPCA_AML <- golub[HPCA_row,gol.fac=="AML"]
cat("The proportion is:",sum(HPCA_ALL<0)/length(HPCA_ALL))
#(c)
#H0: 0.5 population of ALL patients's gene is negatively expressed; H1: >0.5 ...
binom.test(sum(HPCA_ALL<0), n=27, p=0.5, alternative = "greater")
#p-value = 0.221, so we reject null hypothesis
#(d)




 