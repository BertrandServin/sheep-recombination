library(xtable)
grr=read.table('results/family/parent_recombination.txt',head=T)
mod.nmeio=lm(Ri ~ nbMeio,data=grr)
res=summary(mod.nmeio)
coef(res)[2,c(1,4)]

require(xtable)
print(xtable(summary(mod.nmeio),label='tab:nmeio',
       caption='Estimating and testing the effect of %
the number of observed meioses on individual recombination rates.'),include.rownames=FALSE)

plot(grr$nbMeio,grr$Ri,
     xlab='Number of Meioses',
     ylab='Average number of CO / meiosis',
     axes=F,pch=16)
axis(1)
axis(2)

grr=read.table('results/family/nco_meioses.txt')
colnames(grr)=c('offspring','sire','year','insem','nco')
grr$insem.mo=with(grr,month.abb[insem])
grr$insem.mo=ordered(grr$insem.mo,levels=month.abb)
grr$sire=as.factor(grr$sire)
grr$year=as.factor(grr$year)
grr$offspring=as.factor(grr$offspring)
head(grr)

library(lme4)

grr=grr[complete.cases(grr$year)&complete.cases(grr$insem),]

mod.null=lmer(nco ~ (1|sire),data=grr,REML=FALSE)

mod.year=lmer(nco ~ (1|sire)+year,data=grr,REML=FALSE)
## No Year of birth effect
anova(mod.null,mod.year)
## A small insemination month effect
mod.insem=lmer(nco ~ (1|sire)+insem.mo-1,data=grr,REML=FALSE)
ci=confint(mod.insem)[-c(1,2),]

anova(mod.null,mod.insem)

## fit using REML
mod.insem.reml=update(mod.insem,REML=TRUE)
ci=confint(mod.insem.reml)[-c(1,2),]
layout(matrix(c(1,1,1,1,2,2),ncol=2,byrow=T))
plot(2:8,fixef(mod.insem.reml),axes=F,pch=19,
     ylab='Mean number of Crossovers / meiosis',
     xlab='Insemination Month',
     type='b',lty=2,lwd=2,
     ylim=c(33,38.5),cex=2)
axis(1,at=2:8,labels=month.abb[2:8],las=2)
axis(2) 
head(ci)
segments(x0=2:8,x1=2:8,y0=ci[,1],y1=ci[,2],lwd=2)
barplot(table(grr$insem.mo)[2:8],ylab='Number of Inseminations')

parec=read.table('results/family/parent_recombination.txt',h=T)
sol=read.table('results/family/solutions.aireml',skip=1)
## overall mean
mu=sol[1,4]
## Sire effects
blup=sol[sol$V2==4,]
## Correspondance between levels and sires
corresp.ped=read.table('results/family/renadd04.ped')

parec$id=corresp.ped$V1[match(parec$parent,corresp.ped$V10)]
parec$blup=blup$V4[parec$id]

m=mean(parec$Ri)
plot(parec$Ri-m,parec$blup,
     xlab='Mean number of Crossovers per meiosis',
     ylab='Additive Genetic Value (BLUP)',pch=19,axes=F)
axis(1,at=axTicks(1),labels=round(axTicks(1)+m,digits=1))
abline(0,1,lwd=2,col=2)
axis(2)

write.table(parec,file='results/family/parent_recombination_blup.txt',quote=F,row.names=F)
