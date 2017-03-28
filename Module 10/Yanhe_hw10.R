#Question 1
#(a)
library(ArrayExpress)
library(affy)
#yeast.raw <- ArrayExpress('E-MEXP-1551')
getAE('E-MEXP-1551', path = 'yeast',  type = "full")
yeast.raw <-  ReadAffy(celfile.path= 'yeast' )
eset <- expresso(yeast.raw, bgcorrect.method = "mas", normalize.method = "quantiles", 
                 pmcorrect.method = "pmonly", summary.method = "medianpolish")
#(b)
eset_df <- data.frame(exprs(eset))
firstFiveMean <- apply(eset_df[1:5,], 1, mean)
#(c)
numSamples <- ncol(eset_df)
numGenes <- nrow(eset_df)
#or just type eset to see the information

#Question 2
source("https://bioconductor.org/biocLite.R")
#(a)
annotation(eset)
biocLite("yeast2.db")

#(b)
library(yeast2.db)
library(annotate)
go1769308 <- get("1769308_at", env = yeast2GO)
go1769308MF <- getOntology(go1769308, "MF")
length(go1769308MF)

#(c)
library("GO.db")
go1769308MF_parents <- getGOParents(go1769308MF)
go1769308MF_parents_names <- sapply(go1769308MF_parents, function(x) {x$Parents})
cat("There are:",length(unique(unlist(go1769308MF_parents_names))),"parents")
#(d)
go1769308MF_children <- getGOChildren(go1769308MF)
go1769308MF_children_names <- sapply(go1769308MF_children, function(x) {x$Children})
cat("There are",length(unique(unlist(go1769308MF_children_names))),"children")

#Question 3
#(a)
library(ALL)
data(ALL)
library("genefilter")
patients_B23 <- factor(ALL$BT %in% c("B2","B3"))
f1 <- function(x) (wilcox.test(x, exact=FALSE)$p.value<0.001)
f2 <- function(x) {t.test(x~patients_B23)$p.value<0.001}

#(b)
library(limma)
wilcoxon <- genefilter(exprs(ALL[,patients_B23==TRUE]),filterfun(f1))
welcht <- genefilter(exprs(ALL),filterfun(f2))
x <- apply(cbind(wilcoxon,welcht), 2, as.integer)
vc <- vennCounts(x, include="both")
vennDiagram(vc)

#(c)
cat("11876 genes pass wilcoxon test, 749 genes pass both filters")

#(d)
annotation(ALL)
biocLite("hgu95av2.db")
library(hgu95av2.db)
go3d <- get("oncogene", env = hgu95av2GO)

#Question 4
#(a)
allB123 <- ALL[,which(ALL$BT %in% c("B1","B2","B3"))]
#(b)
design.ma <- model.matrix(~ 0 + factor(allB123$BT))
colnames(design.ma) <- c("B1","B2","B3")
fit123 <- lmFit(allB123, design.ma)
fit123 <- eBayes(fit123)
print( topTable(fit123, coef=3, number=5, adjust.method="fdr"), digits=4)
#(c)
cont.ma <- makeContrasts(B1-B2,B2-B3, levels=factor(allB123$BT))
fit_c <- contrasts.fit(fit123, cont.ma)
fit_c <- eBayes(fit_c)
p_values <- fit_c$F.p.value
p_fdrs <- p.adjust(p_values, method = "fdr")
sum(p_fdrs<0.01)
print( topTable(fit_c, number=5, p.value=0.01, adjust.method="fdr"), digits=4)
