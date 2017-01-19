#Question 1
source("http://www.bioconductor.org/biocLite.R")
biocLite()
library(multtest)
data("golub")

#(a)
gol.fac <- factor(golub.cl, levels=0:1, labels = c("ALL","AML"))
meanALL <- apply(golub[,gol.fac=="ALL"], 1, mean)

#(b)
meanAML <- apply(golub[,gol.fac=="AML"], 1, mean)

#(c)
sortedALL <- order(abs(meanALL), decreasing = TRUE)
print(golub.gnames[sortedALL[1:3],2])

#(d)
sortedAML <- order(abs(meanAML), decreasing = TRUE)
print(golub.gnames[sortedAML[1:3],2])

#Question 2
#(a)
AML5 <- golub[1:5, gol.fac=="AML"]
write.csv(AML5, file = "AML5.csv")

#(b)
ALL5 <- golub[1:5, gol.fac=="ALL"]
write.csv(ALL5, file = "ALL5.txt")

#(c)
dataC <- golub[100:200, 1]
print(sd(dataC))

#(d)
sdAllPatient <- apply(golub, 1, sd)
print(sum(sdAllPatient>1))

#(e)
g101 <- golub[101,]
g102 <- golub[102,]
gname101 <- golub.gnames[101,2]
gname102 <- golub.gnames[102,2]
plot(g101, g102, xlab = gname101, ylab = gname102)

#Question 3
library(ALL)
data("ALL")
str(ALL)
openVignette("ALL")

#(a)
ALLB1 <- exprs(ALL[,ALL$BT=="B1"])
hist(ALLB1)

#(b)
meanB1 <- apply(ALLB1, 1, mean)

#(c)
sortedMeanB1 <- order(abs(meanB1), decreasing = TRUE)
print(meanB1[sortedMeanB1[1:3]])


#Question 4
#(a)
data("trees")
print(typeof(trees))

#(b)
plot(trees$Girth,trees$Height, col = "blue", type = "o", pch = "+", xlab = "Girth", ylab= "Height&Valume", ylim = range(1:100))
lines(trees$Girth, trees$Volume, col = "red", type = "o")
legend("bottomright", c("Height", "Valume"), col = c("blue", "red"), lty = c(1,1))
