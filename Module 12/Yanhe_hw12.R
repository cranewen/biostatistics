#Question 1
#(a)
library(ALL)
data(ALL)
ALL.fac <- factor(ALL$BT, levels = c('B','B1','B2','B3','B4','T','T1','T2','T3','T4'), 
                  labels = c(1,1,1,1,1,2,2,2,2,2))

#(b)
par(mfrow = c(1,3))
for(i in 1:3){
  hist(exprs(ALL[i,]))
}
#(c)
ALL1_5 <- exprs(ALL[1:5,])
pairs(t(ALL1_5))

#(d)
require(scatterplot3d)
ALL_39317_at <- data.frame(exprs(ALL["39317_at",]))
ALL_32649_at <- data.frame(exprs(ALL["32649_at",]))
ALL_481_at <- data.frame(exprs(ALL["481_at",]))
group <- rbind(ALL_32649_at,ALL_39317_at,ALL_481_at)
scatterplot3d(t(group), color = ALL.fac)

#(e)
cl.2means <- kmeans(t(group), centers = 2, nstart = 10)
table(cl.2means$cluster)
cl.3means <- kmeans(t(group), centers = 3, nstart = 10)
table(cl.3means$cluster)

#(f)
#alldata <- as.matrix(sapply(data.frame(ALL), as.numeric))
### Firstly, convert the ALL data into dataframe, then use t() to transpose the table.
ALL.data <- t(data.matrix(data.frame(ALL)[,1:nrow(ALL)]))
#alldata <- matrix(as.numeric(unlist(data.frame(ALL))), nrow = nrow(data.frame(ALL)))
# Z <- eigen(cor(ALL.data))
pca <- prcomp(ALL.data, scale = TRUE)
summary(pca)

#(g)
biplot(prcomp(ALL.data, scale = TRUE))

#(h)
print(sort(pca$x[,2], decreasing = TRUE)[1:3])
print(sort(pca$x[,2], decreasing = FALSE)[1:3])

#(i)
library(annotate)
ALL.annotation <- annotation(ALL)
library(hgu95av2.db)
o.decreasing <- order(pca$x[,2], decreasing = TRUE)
o.increasing <- order(pca$x[,2], decreasing = FALSE)
cat("The biggest PC2 value gene name and chromosome:")
mget(c(rownames(ALL[o.decreasing[1]])),hgu95av2GENENAME)
mget(c(rownames(ALL[o.decreasing[1]])),hgu95av2CHRLOC)
cat("The smallest PC2 value gene name and chromosome:")
mget(c(rownames(ALL[o.increasing[1]])),hgu95av2GENENAME)
mget(c(rownames(ALL[o.increasing[1]])),hgu95av2CHRLOC)

#Question 2
#(a)
data("iris")
iris$Species <- NULL
pca.2 <- prcomp(iris, scale = TRUE, center = TRUE)
#(b)
eigen(cor(iris))
prcomp(iris, scale = TRUE)
#(c)
dist(1-cor(pca.2$x))
#(d)
summary(prcomp(pca.2$x, scale = FALSE)) #scaled
summary(prcomp(iris, scale = FALSE)) #unscaled

#(e)
#(f)
nboot <- 1000
bootData <- pca.2$x
p <- ncol(bootData)
n <- nrow(bootData)
sdevs <- array(dim = c(nboot, p))
varsProp <- rep(NA, nboot)
for(i in 1:nboot){
  dat.star <- bootData[sample(1:n,replace = TRUE),]
  sdevs[i,] <- prcomp(dat.star, scale = TRUE)$sdev
  varsProp[i] <- sdevs[i,2]^2/sum(sdevs[i,1:p]^2)
}
quantile(varsProp, c(0.05,0.95))

