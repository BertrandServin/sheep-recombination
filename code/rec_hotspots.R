## read in LD map
rhomap = read.table('results/population/HD_SNP_array_map51.txt', head=T)
## model the distribution of recombination intensities lambda (log scale)
## as a mixture of normal (note the prior distribution is p(log(lambda)~N(.)))
library(mclust)
mc=Mclust(log10(rhomap$lambda),G=2)
## get parameters of the two components
pars = mc$parameters
xx=seq(-3,3,0.01)
d1.pars=c(pars$mean[1],sqrt(pars$variance$scale[1]))
d2.pars=c(pars$mean[2],sqrt(pars$variance$scale[2]))
d1=dnorm(xx, mean = d1.pars[1], sd = d1.pars[2])
d2=dnorm(xx, mean = d2.pars[1], sd = d2.pars[2])
## compute p-values corresponding to H0: interval in background distribution
## vs. H1: hotspot
if (pars$pro[1]>pars$pro[2]) {
    null.pars=d1.pars
} else {
    null.pars=d2.pars
}
pval=pnorm(log10(rhomap$lambda), mean=null.pars[1],
    sd=null.pars[2], lower.tail=F)
## Graphical representation of the fit
par(mfrow=c(1,2))
hist(log10(rhomap$lambda), n=100, freq=F,
     xlab=expression(log10(lambda[i])), main='')
lines(xx, pars$pro[1]*d1, col=2, lwd=2)
lines(xx, pars$pro[2]*d2, col=4, lwd=2)
hist(pval,main='P-Value distribution',xlab='',freq=F,n=100)
abline(h=1,lwd=2,lty=2,col='gray')

library(qvalue)
qval = qvalue(pval)
length(pval)
## Estimate of the number of hotspots in the sheep genome
nhs=length(pval)*(1-qval$pi0)
print(paste('Estimated # hotspots intervals:',nhs))
## call hotspots 
rhomap$hotspot = qval$qvalues < 0.05
rhomap$qvalue = qval$qvalues
hsmap=subset(rhomap, select = c('chr', 'left', 'right', 'qvalue', 'hotspot' ))
write.table(hsmap,row.names=F,quote=F,file='results/population/hotspots51.txt')

require(MESS)
## physical distance for each chromosome
physd = as.numeric(by(rhomap$right, rhomap$chr, max))
## genetic distance for each chromosome
gend = as.numeric(by(rhomap$delta, rhomap$chr, sum))
## Genomic Length (in bp)
genome.L = sum(physd)
## Genetic Length (in rho scale)
genome.G = sum(gend)

rhomap$len=rhomap$right-rhomap$left

## order intervals based on decreasing LD genetic distance
o = sort(rhomap$delta, decreasing=T, index.return=T)$ix
## cumulative proportion of physical distance covered
x = cumsum(as.numeric(rhomap$len[o]))/genome.L
## corresponding proportion of genetic distance covered
y = cumsum(rhomap$delta[o])/genome.G
## Gini coefficient
gini=1-2*auc(y,x)
par(mar=c(5,5,1,1))
plot.new()
plot.window(xlim=c(0,1),ylim=c(0,1))
title( xlab='Cumulated Genetic Distance', ylab='Cumulated Physical Distance',cex.lab=1.5)
xpol=c(x,rev(x))
ypol=c(y,rev(x))
polygon(ypol,xpol,col="Sienna",border=NA)
lines(y, x, lwd=4,axes=F)
axis(1,cex.axis=1.5)
axis(2,cex.axis=1.5)
segments(x0=0,y0=0,x1=1,y1=1,lwd=4,lty=2)
text(0,0.8,labels=paste("Gini coef. =",round(gini,digits=2)),adj=0,cex=1.5,col='Sienna')
