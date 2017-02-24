#Question 1
#(a)
library(multtest)
data("golub")
gol.fac <- factor(golub.cl, levels=0:1, labels = c("ALL","AML"))
p.values1 <- rep(NA, 3051)
for(i in 1:3051){
  p.values1[i] <- wilcox.test(golub[i,gol.fac=="ALL"],golub[i,gol.fac=="AML"], paired = F, alternative = "greater")$p.value
}
p.fdr1 <- p.adjust(p=p.values1, method = "fdr")
sum(p.fdr1<0.05)
#(b)
pIndex1 <- order(p.fdr1, decreasing = FALSE)
for(i in 1:3){
  cat(golub.gnames[pIndex1[i],2],"\n")
}

#Question 2
p.values2 <- rep(NA, 3051)
for(i in 1:3051){
  p.values2[i] <- shapiro.test(golub[i,gol.fac=="AML"])$p.value
}
p.fdr2 <- p.adjust(p=p.values2, method = "fdr")
sum(p.fdr2<0.05)

#Question 3
gene1.num <- grep("HOXA9 Homeo box A9", golub.gnames[,2])
gene2.num <- grep("CD33", golub.gnames[,2])
wilcox.test(golub[gene1.num,gol.fac=="ALL"],golub[gene2.num,gol.fac=="ALL"], paired = T, alternative = "two.sided")
#They do not express at the same level

#Question 4
fisher.test(margin.table(UCBAdmissions,c(2,1)))
#Reject null hypothesis of independence between the admission decision and gender

#Question 5
gol.cd33 <- golub[808,]
n <- length(gol.cd33)
T.obs <- abs(mean(gol.cd33[gol.fac=="ALL"])- mean(gol.cd33[gol.fac=="AML"]))
n.perm <- 2000
T.perm <- rep(NA, n.perm)
for(i in 1:n.perm){
  gol.cd33.perm <- sample(gol.cd33, n, replace = F)
  T.perm[i] <- abs(mean(gol.cd33.perm[gol.fac=="ALL"]) - mean(gol.cd33.perm[gol.fac=="AML"]))
}
mean(T.perm>=T.obs)




