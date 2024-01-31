## pull our GC content percent as a named vector
GC<-setNames(bacteria.data[,"GC_Content_Percent"],rownames(bacteria.data))
head(GC)

## set up for side-by-side plots
par(mfrow=c(1,2),mar=c(6.1,4.1,2.1,1.1))
## histogram of GC content percents on original scale
hist(GC,main="",las=2,xlab="",cex.axis=0.7,cex.lab=0.9,breaks=seq(min(GC),max(GC),length.out=12))
mtext("(a)",adj=0,line=1)
mtext("rate",side=1,line=4,cex=0.9)
## histogram of GC content percents on log scale
ln_GC<-log(GC)
hist(ln_GC,main="",las=2,xlab="",cex.axis=0.7,cex.lab=0.9,breaks=seq(min(ln_GC),max(ln_GC),length.out=12))
mtext("(b)",adj=0,line=1)
mtext("ln(rate)",side=1,line=4,cex=0.9)

## fit BM model
fitBM_GC<-fitContinuous(bacteria.tree,GC,model="BM")
##fit lambda model
fitlambda_GC<-fitContinuous(bacteria.tree,GC,model="lambda",bounds=list(lambda=c(0,1.005)))
## fit OU model
fitOU_GC<-fitContinuous(bacteria.tree,GC,model="OU",bounds=list(alpha=c(0,100)))
## accumulate AIC scores in a vector
aic_GC<-setNames(c(AIC(fitBM_GC),AIC(fitlambda_GC),AIC(fitOU_GC)),c("BM","lambda", "OU"))
## compute and print Akaike weights
aic.w(aic_GC)

#Best model
fitlambda_GC

##phylogenetic signal
fitBMlambda_GC<-phylosig(bacteria.tree,GC,method="lambda", test=T)
fitBMlambda_GC
