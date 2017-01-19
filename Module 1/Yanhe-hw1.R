#Module 1 Homework
#Question 1

#(a)
vec <- c(5,TRUE)
class(vec)

#(b)
x <- 1:4
y <- 1:2
x+y

#(c)
fsin<-function(x) sin(pi*x)
fsin(1)

#(d)
c(1,2) %*% t(c(1,2))

#(e)
f <- function(x) {
  g <- function(y) {
    y+z
  }
  z <- 4
  x + g(x)
}
z <- 15
f(3)


#Question 2
fs <- function(x) x^2
sum(fs(1:1000))


#Question 3
#(a)
k <- c(1:20)
X <- 2*k
print(X)

#(b)
Y <- rep(0,20)
print(Y)

#(c)
integrand <- function(x) sqrt(x)
for (i in 1:20) {
  if (i < 12) {
    Y[i] = 3*k[i]
  }else {
    Y[i] = integrate(integrand, lower = 0, upper = k[i])$value
  }
}
print(Y)
