library(lme4)
library(reshape2)
## read in data and combine samples
fam=read.table('results/combined/comb_1Mb_map_family.txt',head=T)
pop=read.table('results/combined/comb_1Mb_map_pop51.txt',head=T)
tot=rbind(fam,pop)
tot$seg=paste(tot$chr,tot$left,tot$right,sep='-')
tot$seg=as.factor(tot$seg)
chrsize=do.call(rbind,list(by(tot$right,tot$chr,max)))
## remove extremities
for(i in 1:26) {
    end=chrsize[i]
    torm=(tot$chr==i)&(tot$right<=4e6)
    tot=tot[!torm,]
    torm=(tot$chr==i)&(tot$right>=(end-4e6))
    tot=tot[!torm,]
}
## The response considered is the log(10) of recombination rate
tot$logc=log10(tot$value)

## Fit a linear mixed model to estimate local Ne
mymod=lmer(logc~method+(1|chr)+(1|seg/method),data=tot)

library(utils)
library(ascii)
options(asciiType="org")

ss=capture.output(summary(mymod))
ss=as.data.frame(ss)
## rownames and colnames argument to avoid output bug from CRAN version 
ascii(ss,include.rownames=F,include.colnames=F,rownames=1,colnames=2)

require(lme4)
## Genome wide scaling factor for LD-based estimates
log.Ne=fixef(mymod)[2]
## Get predicted values in each segment
pred=tot[tot$rep==0,]
fit.val=predict(mymod,pred)
pred$fit=fit.val
pred.wide=dcast(pred,chr + left + right + seg ~ method,value.var='fit')
linmod=lm(I(pop-log.Ne)~family,data=pred.wide)


fam.c=pred.wide$family
## rescale LD-based estimates
pop.c=pred.wide$pop-log.Ne

## correlation between LD and family based estimates
cor.c=cor(pred.wide$family,pred.wide$pop)

par(mar=c(5,6,1,1))

plot(fam.c,pop.c,
     pch=19,col=rgb(0.63,0.32,0.17,0.2),axes=F,
     xlab='Meiotic recombination rate (cM/Mb)',
     ylab='',
     xlim=c(-1,1),
     ylim=c(-1.5,1),
     cex.lab=1.5
     )
lab.x=c(0.1,0.25,0.5,1,2,4,10)
tks.x=log10(lab.x)
lab.y=c(0.05,0.25,0.1,0.5,1,2,4,10)
tks.y=log10(lab.y)
axis(1,at=tks.x,lab=lab.x,cex.axis=1.3)
axis(2,at=tks.y,lab=lab.y,las=2,cex.axis=1.3)
title(ylab='Scaled historical recombination rate (cM/Mb)',line=4,cex.lab=1.5)
abline(coef=coef(linmod),lwd=2,lty=2)
text(log10(4),log10(0.1),labels=paste('Correlation = ',round(cor.c,digits=2)),cex=1.5)

## Get HD SNP map
snps=read.table('data/population/HDpanel/Lacaune.bim')
colnames(snps)=c('chr','name','gen','pos','A1','A2')

## function to compute number of SNP within a window
nsnp.seg=function(coord,pmap) {
    sub.snp=(pmap$chr==as.integer(coord[1]))&(pmap$pos>=coord[2])&(pmap$pos<=coord[3])
    return(sum(sub.snp))
}

pred.wide$nsnps=apply(pred.wide[,c(1,2,3)],1,nsnp.seg,pmap=snps)

pred.wide$len=pred.wide$right-pred.wide$left
pred.wide$snpdens=pred.wide$nsnps/pred.wide$len

## fit a linear model of historical rate with meiotic rate,
## adjusting for snp density (on a log scale)

linmod=lm(I(pop-log.Ne)~I(log10(snpdens))+family,data=pred.wide)

## We want to look at regions where pop and family are significantly different
require(MASS)
str.resid.pop=studres(linmod)
pred.wide$pval.hi=as.vector(pnorm(str.resid.pop,lower.tail=F))
pred.wide$pval.lo=as.vector(pnorm(str.resid.pop))

hh=hist(str.resid.pop,n=100,
    main='Historical Rates Residuals',
    xlab='Studentized Residuals',
    freq=FALSE)
xx=seq(-10,10,0.01)
lines(xx,dnorm(xx),lwd=2)

require(qvalue)
require(xtable)

pred.wide$resid=as.double(str.resid.pop)
pred.wide$regresid=as.double(residuals(linmod))

pred.wide$quantile=pnorm(qqnorm(pred.wide$family,plot.it=F)$x)
## compute values from 
pred.wide$pval.hi=as.vector(pnorm(pred.wide$resid,lower.tail=F))
pred.wide$pval.lo=as.vector(pnorm(pred.wide$resid))

pval=2*pnorm(abs(pred.wide$resid),lower.tail=F)
pred.wide$pval2side=as.double(pval)
qval=qvalue(pval,pi0.method='bootstrap')

fdr.th=max(qval$pvalues[qval$qvalues<0.02])

zones=subset(pred.wide[pred.wide$pval2side<fdr.th,],
    select=c(chr,left,right,quantile,regresid,pval2side))

zones$chr=as.factor(zones$chr)
zones$left=zones$left*1e-6
zones$right=zones$right*1e-6
zones$quantile=round(zones$quantile,digits=3)
zones$regresid=round(10^zones$regresid,digits=2)
zones$pval2side=format.pval(zones$pval2side,digits=2)
colnames(zones)=c('Chromosome','Left (Mbp)',
            'Right (Mbp)','Meiotic rec. rank','Ratio','p-value')

print(xtable(zones,
  caption='Genome regions where meiotic and historical %
recombination rates differ significantly.'),include.rownames=FALSE,
  type='latex',sanitize.text.function=function(x){x})

write.table(zones,quote=F,row.names=F,file='results/combined/outliers_regions.txt')

## local.Ne estimates log10(4Ne)
  local.Ne=pred.wide$pop-pred.wide$family

  ## get a global estimate for sub-telomeric regions
  genome.Ne=median(local.Ne)

  ## gather local.Ne estimates
  ## For sub-telomeric regions, this is set to genome.Ne
  mypop=pop[pop$rep==0,c(1,2,3)]
  mypop$seg=paste(mypop$chr,mypop$left,mypop$right,sep='-')
  mypop$scale=genome.Ne
  mypop$scale[match(pred.wide$seg,mypop$seg)]=local.Ne

  colnames(mypop)[5]='scale (log10(4Ne))'
  write.table(mypop[,c(1,2,3,5)],row.names=F,quote=F,
              file='results/combined/pop_scale.txt')

  ## Scale our map
  popscale=read.table('results/combined/pop_scale.txt',skip=1)

  pop.hd=read.table('results/population/HD_SNP_array_map.txt',h=T)

  scale.delta=function(seg,scale=popscale) {
      idx=NULL
      seg=unlist(seg)
      idx=which((popscale[,1]==seg[1])&(popscale[,2]<=seg[2])&(popscale[,3]>=seg[3]))
      if (length(idx)<1) {
          idx=max(which((popscale[,1]==seg[1])&(popscale[,2]<=seg[2])))
          return(popscale[idx,4])
      } else {
          return(popscale[idx,4])
      }
  }

  tt=unlist(apply(pop.hd,1,scale.delta))

  pop.hd$scale=tt
  pop.hd$d=100*pop.hd$delta*10^-pop.hd$scale ## in cM
  pop.hd$c=1e6*pop.hd$d/(pop.hd$right-pop.hd$left) ## in cM/Mb
  pop.hd$rho=1e3*pop.hd$delta/(pop.hd$right-pop.hd$left) ## per Kb

write.table(pop.hd,quote=F,row.names=F,
            file='results/combined/HD_SNP_array_map_scaled.txt')

## Now create a bim file with genetic distances, 
## using linear approximation for markers
## snps object contains original bim file data

  for (chrom in 1:26) {
      mysnps=snps$chr==chrom
      mymap=pop.hd$chr==chrom
      myd=c(0,cumsum(pop.hd$d[mymap]))
      myp=c(pop.hd$left[mymap][1],pop.hd$right[mymap])
      ff=approxfun(myp,myd,rule=2)
      snps$gen[mysnps]=round(ff(snps$pos[mysnps]),digits=3)
  }

write.table(snps,quote=F,row.names=F,
            file='results/combined/Illumina_Ovine_HD.bim')
