#Question 1
#(a)
y<-as.numeric(t(read.table(file = "DataPois.txt", header = TRUE)))
size_y <- length(y)
Ey <- mean(y)
#(b)
nloglik <- function(x) -sum(log2(dpois(y, lambda = exp(x))))
optim(par = 1, nloglik)$par

#(c)
nboot <- 1000
boot.data <- rep(NA, nboot)

for(i in 1:nboot){
  data.star <- y[sample(1:size_y,replace = TRUE)]
  nloglik1c <- function(x) -sum(log2(dpois(data.star, lambda = exp(x))))
  boot.data[i] <- optim(par = 1, nloglik1c)$par
}
t.test(boot.data, alternative = "two.sided", mu = 1, conf.level = 0.95)

#Question 2
#(a)
library(ISLR)
ncidata <- NCI60$data
ncilabs <- NCI60$labs
#ncidata <- data.frame(t(ncidata))

a <- which(ncilabs %in% names(which(table(ncilabs)<3))) # an array with indexes which only 1 or 2 repeated
ncidata <- ncidata[-a,] # delete the data by adding - to the indexes
#c<- c(21,24,30,35,36,49,51)
ncilabs <- ncilabs[-a]
#(b)
p <- pairwise.t.test(ncidata[,1], ncilabs, p.adjust.method = "fdr")
#all the p-values are greater than 0.05, so there are no differences in groups 
#(c)
shapiro.test(residuals(lm(ncidata[,1]~ncilabs)))
#p-value is 0.4414 proves the pairwise test is appropriate
#(d)
p_values <- apply(ncidata, 2, function(x) anova(lm(x~ncilabs))[["Pr(>F)"]][1])
p_fdrs <- p.adjust(p_values, method = "fdr")
sum(p_fdrs < 0.05)

#Question 3
#(a)
#x77.data <- as.data.frame(state.x77)
pairs(state.x77)
#From the plot, Illiteracy and Murder variables have linear corraltion with it.
#(b)
statex77.data <- as.data.frame(state.x77)
lin.reg <- lm(statex77.data$`Life Exp`~statex77.data$Income+statex77.data$Illiteracy+statex77.data$Frost)
summary(lin.reg)
# y = 72 + 0.00018*Income - 1.56*Illiteracy - 0.006*Frost
#(c)
require(VGAM)
statex77.lgr <- vglm()


#Question 4
library(ALL)
data(ALL)
#(a)
ALLB <- ALL[, ALL$BT %in% c("B","B1","B2","B3","B4")]
#ALLB <- data.frame(ALLB)
#(b)
library("genefilter")
f_cv <- function(x) {(sd(x)/mean(x))>0.2 }
#patientB <- factor(ALL$BT %in% c("B","B1","B2","B3","B4"))
sel1 <- genefilter(exprs(ALLB), filterfun(f_cv))
sel_ALLB <- data.frame(ALLB)[,sel1]
#sel_ALLB_lables <- factor(sel_ALLB$BT %in% c("B","B1","B2","B3","B4"))
sum(sel1)
#(c)
# We use principle components analysis

#(d)
hcluster.data <- sel_ALLB
#hc.single <- hclust(dist(hcluster.data, method = "euclidean"), method = "single")
hc.ward <- hclust(dist(hcluster.data, method = "euclidean"), method = "ward.D2")
# table(cutree(hc.single, k = 4), ALLB$mol.biol)
# table(cutree(hc.single, k = 4), ALLB$BT)
table(cutree(hc.ward, k = 4), ALLB$mol.biol)
table(cutree(hc.ward, k = 4), ALLB$BT)

#(e)
library(gplots)
#heatmap for B stages
heatmap.2(as.matrix(sel_ALLB), 
          hclustfun = function(d) hclust(dist(d, method = "euclidean"), method = "ward.D2"),
          col=topo.colors(75))
#heatmap for bio names
sel_ALLB.biol <- sel_ALLB[, ALLB$mol.biol]
heatmap.2(as.matrix(sel_ALLB.biol), 
        hclustfun = function(d) hclust(dist(d, method = "euclidean"), method = "ward.D2"),
        col=topo.colors(75))


#(f)
library(limma)
library(plyr)

ALLB1234 <- ALL[, ALL$BT %in% c("B1","B2","B3","B4")]
ALLB1234.factor <- factor(ALLB1234$BT)
ALLB1234.factor <- revalue(ALLB1234.factor, c("B3"="B34","B4"="B34"))

design.ma <- model.matrix(~ 0 + ALLB1234.factor)
colnames(design.ma) <- c("B1","B2","B34")
fit1234 <- lmFit(ALLB1234, design.ma)
fit1234 <- eBayes(fit1234)
cont.ma <- makeContrasts(B1-B2,B2-B34, levels=ALLB1234.factor)
fit_f <- contrasts.fit(fit1234, cont.ma)
fit_f <- eBayes(fit_f)
p_values_f <- fit_f$F.p.value
p_fdrs_f <- p.adjust(p_values_f, method = "fdr")
sum(p_fdrs_f<0.05)

#(g)
library(e1071)
sel_ALLB1234 <- ALLB1234[which(p_fdrs_f<0.05),]
data.svm <- data.frame(sel_ALLB1234, ALLB1234.factor)
mcr.svm.M<-matrix(NA, nrow=dim(sel_ALLB1234)[1], ncol=dim(sel_ALLB1234)[2])
for(i in 1:dim(sel_ALLB1234)[2]){
  data.tmp<-data.svm[-i,]
  p.t <-function(x) t.test(x~data.tmp$ALLB1234.factor)$p.value 
  p.value<-apply(data.tmp[,1:dim(sel_ALLB1234)[1]], 2, p.t)
  for(k in 1:dim(sel_ALLB1234)[1]){
    topknames<-names(data.tmp[,1:dim(sel_ALLB1234)[1]]) 
    topknames<-topknames[order(p.value)] 
    topknames<-topknames[1:k] 
    
    fml<-as.formula(paste("y~",paste(topknames,collapse="+")))
    svm.fit<-svm(fml, data=data.tmp, type = "C-classification", kernel = "linear") 
    svm.pred<-predict(svm.fit,data.svm[i,]) 
    mcr.svm.M[k,i]<- (svm.pred !=data.svm$ALLB1234.factor[i]) 
  }
}

#(h)

#Question 5
#(a)
my.dat <- read.table(file = "DataPoisReg.txt", header = TRUE)
my.dat_x <- my.dat[,1]
my.dat_y <- my.dat[,2]
nloglik5 <- function(x) -sum(log2(dpois(my.dat_y, lambda = exp(x))))
theta5 <- optim(par = 1, nloglik5)$par
lm(my.dat_y~my.dat_x)


#(b)
lgr5 <- lm(dpois(my.dat_y, lambda = theta5)~my.dat_x)
mean(predict(lgr5))
# The slope is -0.5296, not 2.