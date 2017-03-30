library(limma)
library(affy)

Data <- ReadAffy()
eset <- rma(Data)
pData(eset)
strain <- c("lrp+","lrp+","lrp+","lrp+","lrp+","lrp-","lrp-","lrp-","lrp-")
design <- model.matrix(~factor(strain))
colnames(design) <- c("lrp+","lrp-vs+")
design
fit <- lmFit(eset, design)
fit <- eBayes(fit)
options(digits=2)
topTable(fit, coef=2, n=40, adjust="BH")


res<-topTable(fit, number=Inf, adjust.method="none", coef=1)
write.table(res,"oct4_dif_exp.txt",sep="\t")
