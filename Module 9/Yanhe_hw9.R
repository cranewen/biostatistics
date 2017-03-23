#Question 1
library(multtest)
data("golub")
#(a)
GRO2 <- golub[grep("GRO2 GRO2",golub.gnames[,2]),]
GRO3 <- golub[grep("GRO3 GRO3",golub.gnames[,2]),]
cor(GRO2,GRO3)
#(b)
cor.test(GRO2,GRO3, conf.level = 0.9)
#(c)
nboot1 <- 2000
boot1.cor <- matrix(0, nrow = nboot1, ncol = 1)
data1 <- cbind(GRO2,GRO3)
for(i in 1:nboot1){
  dat.star1 <- data1[sample(1:nrow(data1),replace = TRUE),]
  boot1.cor[i,] <- cor(dat.star1[,1],dat.star1[,2])
}
quantile(boot1.cor[,1],c(0.05,0.95))
#(d)
t.test(boot1.cor, mu = 0.64, alternative = "greater", conf.level = 0.95)

#Question 2
#(a)
zyxin_rowNum <- grep("Zyxin",golub.gnames[,2])
zyxin <- golub[zyxin_rowNum,]
# count <- 0
# for(i in 1:nrow(golub)){
#   if(cor(golub[i,],zyxin) < -0.5&i!=zyxin_rowNum){
#     count <- count + 1
#   }
# }
cor2 <- apply(golub,1,function(x) cor(x,zyxin))
sum(cor2 < -0.5)
#(b)
for(i in 1:5){
  cat(golub.gnames[order(cor2,decreasing = FALSE)[i],2],"\n")
}
#(c)
p_values2 <- apply(golub, 1, function(x) cor.test(x,zyxin,alternative = "less")$p.value)
p_fdrs2 <- p.adjust(p_values2, method = "fdr")
sum(p_fdrs2<0.05)

#Question 3
#(a)
reg3.fit <- lm(formula = GRO3~GRO2)
summary(reg3.fit)
#(b)
confint(reg3.fit, level = 0.95)
#(c)
#predict(reg3.fit, newdata = data.frame(GRO2=0), interval = "confidence", level = 0.8)
predict(reg3.fit, newdata = data.frame(GRO2=0), interval = "prediction", level = 0.8)
#(d)
qqnorm(resid(reg3.fit))
qqline(resid(reg3.fit))
shapiro.test(resid(reg3.fit))

#Quesiton 4
#(a)
data("stackloss")
attach(stackloss)
lin_reg4 <- lm(stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.)
summary(lin_reg4)
#stack = 0.7156Air.Flow + 1.2953Water.Temp -0.1521Acid.Conc
#(b)
#airFlow and waterTemp are significant effect on stackloss, acidConc is not. 91.36% of variation 
#(c)
predict(lin_reg4, data.frame(Air.Flow=60, Water.Temp=20, Acid.Conc.=90), interval = "confidence", level = 0.9)
predict(lin_reg4, data.frame(Air.Flow=60, Water.Temp=20, Acid.Conc.=90), interval = "prediction", level = 0.9)

