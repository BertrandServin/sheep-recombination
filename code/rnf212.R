rnf212=read.table('results/gwas/rnf212.assoc.txt',h=TRUE)
pred.pos=c(116426000,116448000)
  mut=sapply(rnf212$rs,function(x){substr(x,start=1,stop=3)})=='RNF'
plot(rnf212$ps,-log10(rnf212$p_wald),
     xlim=c(116.2e6,117e6),
     xlab='Position (Mb)',
     ylab='-log(p)',
     axes=FALSE,col=0)
polygon(x=c(pred.pos,rev(pred.pos)),y=c(0,0,16,16),col='gray',border='gray')
points(rnf212$ps,-log10(rnf212$p_wald),pch=16,col=1)
points(rnf212$ps[mut],-log10(rnf212$p_wald[mut]),pch=19,col=2)
axis(1,at=seq(116.2e6,117e6,1e5),labels=seq(116.2,117,0.1))
axis(2)

assoc.soay=read.table('data/Soay/assoc_oar6.txt',h=T,sep=',')
colnames(assoc.soay)

assoc.soay.m=subset(assoc.soay,Model=='Male')
assoc.soay.f=subset(assoc.soay,Model=='Female')
assoc.soay.a=subset(assoc.soay,Model=='All')

head(rnf212)

layout(matrix(c(1,1,2,3),nrow=4))

pred.pos=c(116426000,116448000)
plot(rnf212$ps,-log10(rnf212$p_wald),
       xlim=c(115.8e6,117e6),
       xlab='Position (Mb)',
       ylab='-log(p)',main='Lacaune GWAS (Males)',
     axes=FALSE,col=0)
## RNF 212 gene
polygon(x=c(pred.pos,rev(pred.pos)),y=c(0,0,16,16),col='gray',border='gray')
points(rnf212$ps,-log10(rnf212$p_wald),pch=16,col='gray')
## Common markers with Soay
common=rnf212$ps%in%assoc.soay.a$Position
points(rnf212$ps[common],-log10(rnf212$p_wald)[common],lwd=2)
## typed mutations
points(rnf212$ps[mut],-log10(rnf212$p_wald[mut]),pch=19,col=2)
axis(1,at=seq(115.8e6,117e6,2e5),labels=seq(115.8,117,0.2))
axis(2)

plot(assoc.soay.m$Position,-log10(assoc.soay.m$Pr.Chisq.),
       pch=16,xlim=c(115.8e6,117e6),
       xlab='Position (Mb)',ylim=c(0,16),
       ylab='-log(p)',main='Soay GWAS (Males)',col=0,
     axes=FALSE)
## RNF 212 gene
polygon(x=c(pred.pos,rev(pred.pos)),y=c(0,0,16,16),col='gray',border='gray')
points(assoc.soay.m$Position,-log10(assoc.soay.m$Pr.Chisq.),col='gray')
common=assoc.soay.m$Position%in%rnf212$ps
points(assoc.soay.m$Position[common],-log10(assoc.soay.m$Pr.Chisq.)[common],lwd=2)
##points(assoc.soay.m$Position,-log10(assoc.soay.m$Pr.Chisq.),col=4)
axis(1,at=seq(115.8e6,117e6,2e5),labels=seq(115.8,117,0.2))
axis(2)

plot(assoc.soay.f$Position,-log10(assoc.soay.f$Pr.Chisq.),
       pch=16,xlim=c(115.8e6,117e6),
       xlab='Position (Mb)',
       ylab='-log(p)',main='Soay GWAS (Females)',col=0,
     axes=FALSE)
## RNF 212 gene
polygon(x=c(pred.pos,rev(pred.pos)),y=c(0,0,16,16),col='gray',border='gray')
points(assoc.soay.f$Position,-log10(assoc.soay.f$Pr.Chisq.),col='gray')
common=assoc.soay.f$Position%in%rnf212$ps
points(assoc.soay.f$Position[common],-log10(assoc.soay.f$Pr.Chisq.)[common],lwd=2)
##points(assoc.soay.m$Position,-log10(assoc.soay.m$Pr.Chisq.),col=4)
axis(1,at=seq(115.8e6,117e6,2e5),labels=seq(115.8,117,0.2))
axis(2)
range(assoc.soay.f$Pr.Chisq.[assoc.soay.f$Pr.Chisq.>0])
