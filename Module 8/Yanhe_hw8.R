#Question 1
#(a)
library(ALL)
data(ALL)
library(lmtest)
ALL_BB1234 <- ALL[,ALL$BT %in% c("B","B1","B2","B3","B4")]
y_109at <- exprs(ALL_BB1234)["109_at",]
anova(lm(y_109at~ALL_BB1234$BT))
#p-value is less than 0.05, we reject null hypothesis, there is difference among different 
#stages of patients, so the stages do affect the mean

#(b)
mean(exprs(ALL[,ALL$BT %in% c("B3")])["109_at",])

#(c)
pairwise.t.test(y_109at, ALL_BB1234$BT)
#No group's mean is different from group B.

#(d)
pairwise.t.test(y_109at, ALL_BB1234$BT, p.adjust.method = "fdr")
#Still showing all the p-values are greater than 0.05, so their means are no different from group B.

#(e)
shapiro.test(residuals(lm(y_109at~ALL_BB1234$BT)))
bptest(lm(y_109at~ALL_BB1234$BT))
kruskal.test(y_109at~ALL_BB1234$BT)

#Question 2
#(a)
y_2 <- exprs(ALL_BB1234)
# test for apply():
# count <- 0
# for(i in 1:nrow(y_2)){
#   if(kruskal.test(y_2[i,]~ALL_BB1234$BT)$p.value<0.05){
#     count = count + 1
#   }
# }
p_values <- apply(y_2, 1, function(x) kruskal.test(x~ALL_BB1234$BT)$p.value)
p_fdrs <- p.adjust(p_values, method = "fdr")
cat("There are",sum(p_fdrs<0.05),"genes expressed differently.")
#(b)
for(i in 1:5){
cat(rownames(y_2)[order(p_fdrs, decreasing = FALSE)][i],"\n")
}

#Question 3
#(a)
ALLBg <- ALL[,which(ALL$BT %in% c("B1","B2","B3","B4") & ALL$sex %in% c("M","F"))]
y_3 <- exprs(ALLBg)["38555_at",]
Bcell<-ALLBg$BT
gender <- ALLBg$sex
anova(lm(y_3~Bcell*gender))
anova(lm(y_3~Bcell+gender))
summary(lm(y_3~Bcell+gender))
#(b)
shapiro.test(residuals(lm(y_3~Bcell*gender)))
bptest(lm(y_3~Bcell*gender))
kruskal.test(y_3~Bcell*gender)

#Question 4
ALLB123 <- ALL[,ALL$BT %in% c("B1","B2","B3")]
y_4 <- exprs(ALLB123)["1242_at",]
group <- ALLB123$BT[,drop=T]
n <- length(y_4)
#T.obs <- anova(lm(y_4~group))$F[1]
T.obs <- max(by(y_4,group,mean)) - min(by(y_4,group,mean))
n.perm <- 2000
T.perm <- rep(NA, n.perm)
for(i in 1:n.perm) {
  y_4.perm = sample(y_4, n, replace=F) 
  #T.perm[i] = anova(lm(y_4.perm~group))$F[1]
  T.perm[i] = max(by(y_4.perm,group,mean))-min(by(y_4.perm,group,mean))
}
mean(T.perm>=T.obs)
