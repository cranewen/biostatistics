#Question 2
fx <- function(x) 2^x/factorial(x)*exp(-2)
#Find P(X = 1)
px1 <- fx(1)
#Find P(-2<X<4)
px2 <- fx(0) + fx(1) + fx(2) + fx(3)

#Question 3
# n = 3, p = 0.25

#Question 4
#P(X<=2)
vec4.range <- c(0,1,2)
sum(dbinom(vec4.range, 3, 0.25))
Ex4 <- sum(vec4.range*dbinom(vec4.range, 3, 0.25))
Var4 <- sum((vec4.range-Ex4)^2*dbinom(vec4.range, 3, 0.25))

#Question 5
fx5 <- function(x) dchisq(x,3)
#P(1 < X < 4)
px5 <- integrate(fx5, lower = 1, upper = 4)$value
Ex5 <- integrate(function(x) x*fx5(x), lower = 1, upper = 4)$value
Var5 <- integrate(function(x) (x-Ex5)^2*fx5(x), lower = 1, upper = 4)$value
#Monte Carlo Simulation
mcData5 <- rchisq(100000, 3)
sum <- 0
for (i in 1:100000){
  if(mcData5[i]<4 && mcData5[i]>1) {
    sum = sum + 1
  }
}
pm5 <- sum/100000 #yes

#Question 6
#E(aX+b) = aE(X) + b    and   Var(aX+b) = a^2*Var(X)
#E(Y) = 4*E(X) + 10 = 30
#Var(Y) = 4^2*10 = 160
#It doesn't follow m = 10

#Question 7
#(a)
p7a <- pnorm(1.6, 1.6, 0.4) - pnorm(1, 1.6, 0.4)
#(b)
nsim7 <- 500000
mcData7 <- rnorm(nsim7, 1.6, 0.4)
sum7 <- 0
for (i in 1:500000) {
  if(mcData7[i]<=1.6 && mcData7[i]>=1){
    sum7 = sum7 + 1
  }
}
p7b <- sum7/nsim7

#(c)
p7c <- dbinom(2, 5, p7a)

#Question 8
#(a)
X <- rf(10000000, 2, 5)
meanX <- mean(X)
varX <- var(X)
Y <- rf(10000000, 10, 5)
meanY <- mean(Y)
varY <- var(Y)
#(b)
meanXformula <- 5/(5-2)
varXformula <- 2*5^2*(5+2-2)/(2*(5-2)^2*(5-4))
meanYformula <- 5/(5-2)
varYformula <- 2*5^2*(5+10-2)/(10*(5-2)^2*(5-4))
#Yes
