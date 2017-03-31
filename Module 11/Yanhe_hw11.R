#Question 1
#(a)
library(multtest)
data("golub")
gol.fac <- factor(golub.cl, levels=0:1, labels = c("ALL","AML"))
ccnd3 <- golub[grep("CCND3 Cyclin D3", golub.gnames[,2]),]
clust.data1a <- data.frame(ccnd3)
hc.single <- hclust(dist(clust.data1a, method = "euclidean"), method = "single")
hc.ward <- hclust(dist(clust.data1a, method = "euclidean"), method = "ward.D2")
plot(hc.single, labels=gol.fac)
plot(hc.ward, labels=gol.fac)
table(cutree(hc.single, k = 2), gol.fac)
table(cutree(hc.ward, k = 2), gol.fac)
#ward is better
#(b)
cl.2means1b <- kmeans(clust.data1a, centers = 2, nstart = 10)
table(cl.2means1b$cluster, labels=gol.fac)
#(c)
#Hierarchical is the best
#(d)
initial <- cl.2means1b$centers
n <- dim(clust.data1a)[1]
nboot <- 2000
boot.cl <- matrix(NA, nrow = nboot, ncol = 4) # column 4 to store CI for each mean, here is 2 means
for(i in 1:nboot){
  dat.star <- clust.data1a[sample(1:n, replace = TRUE),]
  cl <- kmeans(dat.star, initial, nstart = 10)
  boot.cl[i,] <- c(cl$centers[1], cl$centers[2])
}
apply(boot.cl,2,mean)
quantile(boot.cl[,1],c(0.025,0.975))
quantile(boot.cl[,2],c(0.025,0.975))
quantile(boot.cl[,3],c(0.025,0.975))
quantile(boot.cl[,4],c(0.025,0.975))
#(e)
K <- c(1:30)
sse <- rep(NA,length(K))
for(k in K){
  sse[k] <- kmeans(clust.data1a, centers = k, nstart = 10)$tot.withinss
}
plot(K, sse, type = 'o', xaxt = 'n')
axis(1, at = K, las = 2)
#From the plot, it suggests 5 clusters.

#Problem 2
#(a)
oncogenes <- golub[grep("oncogene",golub.gnames[,2]),]
antigens <- golub[grep("antigen",golub.gnames[,2]),]
#(b)
#Using t(data frame object or matrix) to tranpose the matrix
cl.2meansOncogenes <- kmeans(t(oncogenes), centers = 2, nstart = 10)
cl.2meansAntigens <- kmeans(t(antigens), centers = 2, nstart = 10)
