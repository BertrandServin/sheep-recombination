dat.1m=read.table('results/family/1Mb_map.txt',h=T,fill=NA)
dat.snp=read.table('results/family/SNP_array_map.txt',h=T,fill=NA)
pdf(file='figures/recombination_map.pdf',w=14,h=12)
for (ichr in 1:26) {

    thechr=dat.snp[ dat.snp$chr==ichr,]
    thechr$mids=apply(thechr[,c(2,3)],1,mean)
    thechr.1m=dat.1m[ dat.1m$chr==ichr,]
    thechr.1m$mids=apply(thechr.1m[,c(2,3)],1,mean)

    for (xl in seq(0,max(thechr$left),100e6) ) {
        xr=xl+100e6

        par(mfrow=c(2,1))
        ## 1Mb map
        mychr=thechr.1m[thechr.1m$left>xl & thechr.1m$left<xr,]
        if (dim(mychr)[1]==0) {
            next
        }
        xx=mychr$left
        xx=rep(xx,each=2)
        xx=xx[-1]
        yy.q5=mychr$q5_cj
        yy.q5=rep(yy.q5,each=2)
        yy.q5=yy.q5[-length(yy.q5)]
        yy.q95=mychr$q95_cj
        yy.q95=rep(yy.q95,each=2)
        yy.q95=yy.q95[-length(yy.q95)]

        plot(mychr$mids,mychr$m_cj,pch=16,type='n',ylim=c(0,max(dat.1m$q95_cj)),
             main=paste('Chromosome',ichr,'1Mb Windows'),axes=F,xlim=c(xl,xr),
             xlab='Position (Mbp)',ylab='Recombination Rate (cM/Mb)')
        axis(2)
        axis(1,at=seq(0,max(mychr$right),10e6),
             label=round(1e-6*seq(0,max(mychr$right),10e6),digit=0))
        polygon(c(xx,rev(xx)),c(yy.q5,rev(yy.q95)),col='gray',border='gray')
        points(mychr$mids,mychr$m_cj,pch=16,type='b',lty=3)

        ## SNP array map
        xr=xl+100e6
        mychr=thechr[thechr$left>xl & thechr$left<xr,]
        xx=mychr$left
        xx=rep(xx,each=2)
        xx=xx[-1]
        yy.q5=mychr$q5_cj
        yy.q5=rep(yy.q5,each=2)
        yy.q5=yy.q5[-length(yy.q5)]
        yy.q95=mychr$q95_cj
        yy.q95=rep(yy.q95,each=2)
        yy.q95=yy.q95[-length(yy.q95)]

        plot(mychr$mids,mychr$m_cj,pch=16,type='n',ylim=c(0,max(dat.snp$q95_cj)),
             main=paste('Chromosome',ichr,'SNP array'),axes=F,xlim=c(xl,xr),
             xlab='Position (Mbp)',ylab='Recombination Rate (cM/Mb)')
        axis(2)
        axis(1,at=seq(0,max(mychr$right),10e6),
             label=round(1e-6*seq(0,max(mychr$right),10e6),digit=0))
        polygon(c(xx,rev(xx)),c(yy.q5,rev(yy.q95)),col='gray',border='gray')
        points(mychr$mids,mychr$m_cj,pch=16,type='b',lty=3)

    }

}
dev.off()

ichr=24
thechr=dat.snp[ dat.snp$chr==ichr,]
thechr$mids=apply(thechr[,c(2,3)],1,mean)
thechr.1m=dat.1m[ dat.1m$chr==ichr,]
thechr.1m$mids=apply(thechr.1m[,c(2,3)],1,mean)
mychr=thechr.1m
xx=mychr$left
xx=rep(xx,each=2)
xx=xx[-1]
yy.q5=mychr$q5_cj
yy.q5=rep(yy.q5,each=2)
yy.q5=yy.q5[-length(yy.q5)]
yy.q95=mychr$q95_cj
yy.q95=rep(yy.q95,each=2)
yy.q95=yy.q95[-length(yy.q95)]
par(mfrow=c(2,1))
plot(mychr$mids,mychr$m_cj,pch=16,type='n',ylim=c(0,max(dat.1m$q95_cj)),
     main=paste('Chromosome',ichr,'1Mb Windows'),axes=F,
     xlab='Position (Mbp)',ylab='Recombination Rate (cM/Mb)')
axis(2)
axis(1,at=seq(0,max(mychr$right),10e6),
     label=round(1e-6*seq(0,max(mychr$right),10e6),digit=0))
polygon(c(xx,rev(xx)),c(yy.q5,rev(yy.q95)),col='gray',border='gray')
points(mychr$mids,mychr$m_cj,pch=16,type='b',lty=3)

mychr=thechr
xx=mychr$left
xx=rep(xx,each=2)
xx=xx[-1]
yy.q5=mychr$q5_cj
yy.q5=rep(yy.q5,each=2)
yy.q5=yy.q5[-length(yy.q5)]
yy.q95=mychr$q95_cj
yy.q95=rep(yy.q95,each=2)
yy.q95=yy.q95[-length(yy.q95)]

plot(mychr$mids,mychr$m_cj,pch=16,type='n',ylim=c(0,max(dat.snp$q95_cj)),
     main=paste('Chromosome',ichr,'SNP array'),axes=F,
     xlab='Position (Mbp)',ylab='Recombination Rate (cM/Mb)')
axis(2)
axis(1,at=seq(0,max(mychr$right),10e6),
     label=round(1e-6*seq(0,max(mychr$right),10e6),digit=0))
polygon(c(xx,rev(xx)),c(yy.q5,rev(yy.q95)),col='gray',border='gray')
points(mychr$mids,mychr$m_cj,pch=16,type='b',lty=3)

dat.1m$x=0
dat.1m$dtelo=0
dat.1m$mids=apply(dat.1m[,c(2,3)],1,mean)
for (ichr in 1:26) {
    subset=dat.1m$chr==ichr
    chr.size=max(dat.1m$right[subset])
    dat.1m$x[subset]=dat.1m$mids[subset]/chr.size
    dpossible=cbind(dat.1m$mids[subset],chr.size-dat.1m$mids[subset])
    dat.1m$dtelo[subset]=apply(dpossible,1,min)
}

par(mfrow=c(1,3),mar=c(6,6,4,1))
## metacentric
plot.new()
plot.window(xlim=c(0,1),ylim=c(0,7))
title(xlab='Relative Position on chromosome',ylab='Recombination Rate (cM/Mb)',main='Metacentric autosomes',cex.lab=1.5)
axis(1)
axis(2)
metacent=dat.1m[as.integer(dat.1m$chr)<4,]
ss.m=smooth.spline(metacent$x,metacent$m_cj,df=20)
points(metacent$x,metacent$m_cj,col='gray',pch=16)
lines(ss.m,col=2,lwd=3)
## acrocentric
plot.new()
plot.window(xlim=c(0,1),ylim=c(0,8))
title(xlab='Relative Position on chromosome',ylab='Recombination Rate (cM/Mb)',main='Acrocentric autosomes',cex.lab=1.5)
axis(1)
axis(2)
acrocent=dat.1m[as.integer(dat.1m$chr)>3,]
ss.a=smooth.spline(acrocent$x,acrocent$m_cj,df=5)
points(acrocent$x,acrocent$m_cj,col='gray',pch=16)
lines(ss.a,col=2,lwd=3)

## distance to telomere
plot(dat.1m$dtelo,dat.1m$m_cj,pch=16,col=rgb(0,0,0,0.2),xlim=c(0,6e7),
     ylab='Recombination Rate (cM/Mb)',xlab='',
     main='All autosomes', cex.lab=1.5,
     axes=F)
title(xlab='Distance to Nearest\n Chromosome end (Mb)',line=4,cex.lab=1.5)
axis(1,at=seq(0,6e7,1e7),labels=seq(0,60,10))
axis(2)
ss.telo=smooth.spline(dat.1m$dtelo,dat.1m$m_cj,df=5)
lines(ss.telo,lwd=3,col=2)
abline(v=4e6,lty=3)


## remove regions close to chromosome ends
chr.size=do.call(cbind,list(by(dat.snp$right,dat.snp$chr,max)))
rec.1m=dat.1m
rec.snp=dat.snp
for (ichr in 1:26) {
    subset=(rec.1m$chr==ichr)
    subset= subset & ((rec.1m$left<4e6) | (rec.1m$right > (chr.size[ichr]-4e6)))
    rec.1m=rec.1m[!subset,]
    subset=(rec.snp$chr==ichr)  
    subset=subset & ((rec.snp$left<4e6) | (rec.snp$right > (chr.size[ichr]-4e6)))
    rec.snp=rec.snp[!subset,]
}
rec.1m$chr=as.factor(rec.1m$chr)
rec.snp$chr=as.factor(rec.snp$chr)

library(Hmisc)
library(xtable)
### CHROMOSOME REC RATE 1 Mb WINDOWS
mod=lm(m_cj~chr-1,data=rec.1m,weights=1/s_cj^2)
res=summary(mod)
rec=res$coefficients[,1]
rec.std.err=res$coefficients[,2]
err=confint(mod)
mod.size.log=lm(rec~log(chr.size),weights=1/rec.std.err^2)
mod.size.inv=lm(rec~I(1/chr.size),weights=1/rec.std.err^2)
mod.size.log.b=coefficients(mod.size.log)
mod.size.inv.b=coefficients(mod.size.inv)

size2rec.log=function(s) {
    mod.size.log.b[1]+mod.size.log.b[2]*log(s)
}
size2rec.inv=function(s) {
    mod.size.inv.b[1]+mod.size.inv.b[2]/s
}

print(xtable(summary(mod),caption='Chromosome effect model estimates'
      ,label='tab:chr.rec.1m'),include.rownames=FALSE)

tw=summary(mod)$coefficients[,c(1,2)]
colnames(tw)=c('c','se')
write.table(tw,file='results/family/chromosome_rates.txt',quote=F)

par(mfcol=c(1,2))
par(mar=c(5,4,4,1))
errbar(1:26,rec,yplus=err[,2],yminus=err[,1],
       xlab='Chromosome',ylab='Mean recombination rate (cM/Mb)')
title(main='One Megabase Windows')
par(mar=c(5,4,1,1))
errbar(chr.size,rec,yplus=err[,2],yminus=err[,1],
       xlab='Chromosome Size (Mb)',ylab='Mean recombination rate (cM/Mb)',
       lwd=2,axes=F,col=1,pch='.',errbar.col='gray')
points(chr.size,rec,col=1,pch=16,cex=0.7)
axis(1,at=seq(50,250,25)*1e6,labels=seq(50,250,25))
axis(2)
text(chr.size,rec,1:26,adj=-0.2)
xx=seq(30e6,300e6,1e4)
lines(xx,size2rec.log(xx),col=4,lty=3,lwd=2)
lines(xx,size2rec.inv(xx),col=4,lty=2,lwd=2)
legend(200e6,1.6,legend=c('log(size)','1/size'),lty=c(3,2),col=4)

### CHROMOSOME REC RATE SNP ARRAY
mod=lm(m_cj~chr-1,data=rec.snp,weights=1/s_cj^2)
res=summary(mod)
rec=res$coefficients[,1]
rec.std.err=res$coefficients[,2]
err=confint(mod)
mod.size.log=lm(rec~log(chr.size),weights=1/rec.std.err^2)
mod.size.inv=lm(rec~I(1/chr.size),weights=1/rec.std.err^2)

mod.size.log.b=coefficients(mod.size.log)
mod.size.inv.b=coefficients(mod.size.inv)
size2rec.log=function(s) { mod.size.log.b[1]+mod.size.log.b[2]*log(s) }
size2rec.inv=function(s) { mod.size.inv.b[1]+mod.size.inv.b[2]/s }

print(xtable(summary(mod),caption='Chromosome effect model estimates (SNP array)'
        ,label='tab:chr.rec.snp'),include.rownames=FALSE)

par(mfcol=c(1,2))
par(mar=c(5,4,4,1))
errbar(1:26,rec,yplus=err[,2],yminus=err[,1],
       xlab='Chromosome',ylab='Mean recombination rate (cM/Mb)')
title(main='SNP array Windows')
par(mar=c(5,4,1,1))
errbar(chr.size,rec,yplus=err[,2],yminus=err[,1],
       xlab='Chromosome Size (Mb)',ylab='Mean recombination rate (cM/Mb)',
       lwd=2,axes=F,col=1,pch='.',errbar.col='gray')
points(chr.size,rec,col=1,pch=16,cex=0.7)
axis(1,at=seq(50,250,25)*1e6,labels=seq(50,250,25))
axis(2)
text(chr.size,rec,1:26,adj=-0.2)
xx=seq(30e6,300e6,1e4)
lines(xx,size2rec.log(xx),col=4,lty=3,lwd=2)
lines(xx,size2rec.inv(xx),col=4,lty=2,lwd=2)
legend(200e6,1.5,legend=c('log(size)','1/size'),lty=c(3,2),col=4)

## One Mb windows
gc1=read.table('data/genome/GC_Content_1Mb.txt',head=TRUE)
## SNP array intervals
gcs=read.table('data/genome/GC_Content_SNP_array.txt',head=TRUE)
## check the data are aligned
require(assertthat)
assert_that(mean(gcs$left==dat.snp$left)==1)
assert_that(mean(gc1$left==dat.1m$left)==1)

gc1$cj=dat.1m$m_cj
gcs$cj=dat.snp$m_cj


## remove extreme regions
for (ichr in 1:26) {
    subset=(gc1$chr==ichr)
    subset= subset & ((gc1$left<4e6) | (gc1$right > (chr.size[ichr]-4e6)))
    gc1=gc1[!subset,]
    subset=(gcs$chr==ichr)  
    subset=subset & ((gcs$left<4e6) | (gcs$right > (chr.size[ichr]-4e6)))
    gcs=gcs[!subset,]
}

par(mfrow=c(1,2))
plot(gc1$gc,gc1$cj,main='One Mb Windows',
     xlab='GC content (%)',
     ylab='Rec. rate (cM/Mb)',pch=16,
     col=rgb(0,0,0,0.1),
     axes=F)
axis(1)
axis(2)
plot(gcs$gc,gcs$cj,main='SNP array intervals',
     xlab='GC content (%)',ylab='Rec. rate (cM/Mb)',
     col=rgb(0,0,0,0.1),
     pch=16,axes=F)
axis(1)
axis(2)

dist.gc1=ecdf(gc1$gc)
dist.gcs=ecdf(gcs$gc)
## quantiles of the gc distribution
gc1$qgc=dist.gc1(gc1$gc)
gcs$qgc=dist.gcs(gcs$gc)

par(mar=c(5,4,1,1),mgp=c(3,1.5,0),mfrow=c(1,2))

q1=quantile(gc1$gc,p=c(0.1,0.25,0.5,0.75,0.9))
ll1=paste(round(q1,digits=1),'\n(',names(q1),')',sep='')

qs=quantile(gcs$gc,p=c(0.1,0.25,0.5,0.75,0.9))
lls=paste(round(qs,digits=1),'\n(',names(qs),')',sep='')

par(mfrow=c(1,2)) 
  plot(gc1$qgc,gc1$cj,main='One Mb Windows',
       xlab='GC content (quantile scale)',
       ylab='Rec. rate (cM/Mb)',pch=16,
       col=rgb(0,0,0,0.1),
       axes=F)
  axis(1,at=dist.gc1(q1),labels=ll1,adj=1)
  axis(2)
  plot(gcs$qgc,gcs$cj,main='SNP array intervals',
       xlab='GC content (quantile scale)',ylab='Rec. rate (cM/Mb)',
       col=rgb(0,0,0,0.1),
       pch=16,axes=F)
  axis(1,at=dist.gcs(qs),labels=lls,adj=1)
  axis(2)

## Model on 1 Mb windows, raw GC covariate:
gc1$chr=as.factor(gc1$chr)
lm.1m.chr=lm(cj~chr,data=gc1)
lm.1m.chr.gc=lm(cj~chr+gc,data=gc1)
aov.1m.gc=anova(lm.1m.chr,lm.1m.chr.gc)
pval.1m.gc=aov.1m.gc['Pr(>F)'][2,1]

## Model on 1Mb windows, quantile-transformed gc
lm.1m.chr.qgc=lm(cj~chr+qgc,data=gc1)
aov.1m.qgc=anova(lm.1m.chr,lm.1m.chr.qgc)
pval.1m.qgc=aov.1m.qgc['Pr(>F)'][2,1]

## Model on SNP windows, raw GC covariate:
gcs$chr=as.factor(gcs$chr)
lm.snp.chr=lm(cj~chr,data=gcs)
lm.snp.chr.gc=lm(cj~chr+gc,data=gcs)
aov.snp.gc=anova(lm.snp.chr,lm.snp.chr.gc)
pval.snp.gc=aov.snp.gc['Pr(>F)'][2,1]

## Model on SNP windows, quantile-transformed gc
lm.snp.chr.qgc=lm(cj~chr+qgc,data=gcs)
aov.snp.qgc=anova(lm.snp.chr,lm.snp.chr.qgc)
pval.snp.qgc=aov.snp.qgc['Pr(>F)'][2,1]

res.gc=data.frame("Intervals"=c('1 Mb','1 Mb','SNP','SNP'),
      "GC Model"=c('raw','transformed','raw','transformed'),
      "logp"=-log10(c(pval.1m.gc,pval.1m.qgc,pval.snp.gc,pval.snp.qgc)))
  colnames(res.gc)[2]=c('GC content')
  colnames(res.gc)[3]=c('$-\\log_{10}(p)$')
  print(xtable(res.gc,
caption='Significance of GC content effect on recombination rate'),
type='latex',include.rownames=FALSE,sanitize.text.function=function(x){x})

write.table(gc1,file='results/family/1Mb_map_annotated.txt',
            col.names=T,row.names=F,quote=FALSE)
write.table(gcs,file='results/family/SNP_array_map_annotated.txt',
            col.names=T,row.names=F,quote=FALSE)

hsmap=read.table('results/population/hotspots51.txt',h=T)



## Get number of hotspots in intervals
getnhs=function(v,hsdat=hsmap) {
    v=as.integer(v)
    hs.loc=hsdat[(hsdat[,1]==v[1])&(hsdat[,2]>=v[2])&(hsdat[,3]<=v[3]),]
    return(sum(hs.loc[,5]))
}

gc1$nhs=apply(gc1[,c(1,2,3)],1,getnhs)
gcs$nhs=apply(gcs[,c(1,2,3)],1,getnhs)

## calculate hotspot density (in HS / 10Kb) for SNP array map
gcs$len=gcs$right-gcs$left
gcs$hsdens=1e4*gcs$nhs/gcs$len

csnp.nhs=cor.test(gcs$cj,gcs$nhs)
csnp.hsdens=cor.test(gcs$cj,gcs$hsdens)
c1.nhs=cor.test(gc1$cj,gc1$nhs)
pval.vec=c(c1.nhs$p.value,csnp.nhs$p.value,csnp.hsdens$p.value)
res.cor=data.frame("Intervals"=c('1 Mb','SNP','SNP'),
    "HS model"=c('number/density','number','density'),
    "corr."=c(c1.nhs$estimate,csnp.nhs$estimate,csnp.hsdens$estimate),
    "p-value"=-log10(pval.vec))

colnames(res.cor)[2]=c('HS effect')
colnames(res.cor)[3]=c('Correlation')
colnames(res.cor)[4]=c('$-\\log_{10}(p)$')

print(xtable(res.cor,
  caption='Correlation between historical hotspots and meiotic recombination rate.%
  Number:  number of hotspots.%
  Density:  density of hotspots (in HS/10Kb).'),
  type='latex',include.rownames=FALSE,sanitize.text.function=function(x){x})

lm.1m.chr.qgc.nhs=lm(cj ~ chr + qgc + nhs, data=gc1)
aov.1m.qgc.nhs=anova(lm.1m.chr.qgc,lm.1m.chr.qgc.nhs)
pval.1m.qgc.nhs=aov.1m.qgc.nhs['Pr(>F)'][2,1]

lm.snp.chr.qgc.nhs=lm(cj ~ chr + qgc + nhs, data=gcs)
aov.snp.qgc.nhs=anova(lm.snp.chr.qgc,lm.snp.chr.qgc.nhs)
pval.snp.qgc.nhs=aov.snp.qgc.nhs['Pr(>F)'][2,1]


lm.snp.chr.qgc.hsdens=lm(cj ~ chr + qgc + hsdens, data=gcs)
aov.snp.qgc.hsdens=anova(lm.snp.chr.qgc,lm.snp.chr.qgc.hsdens)
pval.snp.qgc.hsdens=aov.snp.qgc.hsdens['Pr(>F)'][2,1]

res.hs=data.frame("Intervals"=c('1 Mb','SNP','SNP'),
        "HS model"=c('number/density','number','density'),
    "logp"=-log10(c(pval.1m.qgc.nhs,pval.snp.qgc.nhs,pval.snp.qgc.hsdens)))


colnames(res.hs)[2]=c('HS effect')
colnames(res.hs)[3]=c('$-\\log_{10}(p)$')

print(xtable(res.hs,
caption='Significance of hotspot effect on meiotic recombination rate.%
Number: effect of the number of hotspots.%
Density: effect of the density of hotspots (in HS/10Kb).'),
type='latex',include.rownames=FALSE,sanitize.text.function=function(x){x})

csnp.nhs=cor.test(residuals(lm.snp.chr.qgc),gcs$nhs)
csnp.hsdens=cor.test(residuals(lm.snp.chr.qgc),gcs$hsdens)
c1.nhs=cor.test(residuals(lm.1m.chr.qgc),gc1$nhs)
pval.vec=c(c1.nhs$p.value,csnp.nhs$p.value,csnp.hsdens$p.value)
res.cor.res=data.frame("Intervals"=c('1 Mb','SNP','SNP'),
    "HS model"=c('number/density','number','density'),
    "corr."=c(c1.nhs$estimate,csnp.nhs$estimate,csnp.hsdens$estimate),
    "p-value"=-log10(pval.vec))

colnames(res.cor.res)[2]=c('HS effect')
colnames(res.cor.res)[3]=c('Correlation')
colnames(res.cor.res)[4]=c('$-\\log_{10}(p)$')

print(xtable(res.cor.res,
  caption='Correlation between historical hotspots and meiotic recombination rate%
  corrected for chromosome and GC content effects.%
  Number:  number of hotspots.%
  Density:  density of hotspots (in HS/10Kb).'),
  type='latex',include.rownames=FALSE,sanitize.text.function=function(x){x})
