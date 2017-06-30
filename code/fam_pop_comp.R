dat.1m=read.table('results/family/1Mb_map.txt',h=T,fill=NA)
dat.snp=read.table('results/family/SNP_array_map.txt',h=T,fill=NA)
dat.comb=read.table('results/combined/compare_family_60K.txt',h=T)
pop.hd=read.table('results/population/HD_SNP_array_map.txt',h=T)

pop.hd$rho=1e3*pop.hd$delta/(pop.hd$right-pop.hd$left) ## per Kb

## Whole chromosome
ichr = 24
thechr=dat.snp[ dat.snp$chr==ichr,]
thechr$mids=apply(thechr[,c(2,3)],1,mean)
thechr.1m=dat.1m[ dat.1m$chr==ichr,]
thechr.1m$mids=apply(thechr.1m[,c(2,3)],1,mean)
mychr=thechr.1m 

png(file='figures/chr24_pop_fam_comp.png',w=14,h=20,unit='cm',res=200.)
M = matrix(c(1,1,2,5,3,6,4,7), 4, 2, byrow = TRUE)
layout(M)

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
     main=paste('Chromosome',ichr),axes=F,
     xlab='Position (Mb)',ylab='Recombination rate c')
axis(2)
axis(1,at=seq(0,max(mychr$right),1e6),label=NA)
axis(1,at=seq(0,max(mychr$right),10e6),
     label=round(1e-6*seq(0,max(mychr$right),10e6),digit=0))
polygon(c(xx,rev(xx)),c(yy.q5,rev(yy.q95)),col='gray',border='gray')
points(mychr$mids,mychr$m_cj,pch=16,type='b',lty=3)

## highlight windows
hi.win=6
xx=c(mychr$left[hi.win],mychr$right[hi.win])
yy.q5=mychr$q5_cj[hi.win]
yy.q5=rep(yy.q5,each=2)
yy.q95=mychr$q95_cj[hi.win]
yy.q95=rep(yy.q95,each=2)
polygon(c(xx,rev(xx)),c(yy.q5,rev(yy.q95)),col=NULL,border='darkred',lwd=2)


low.win=19
xx=c(mychr$left[low.win],mychr$right[low.win])
yy.q5=mychr$q5_cj[low.win]
yy.q5=rep(yy.q5,each=2)
yy.q95=mychr$q95_cj[low.win]
yy.q95=rep(yy.q95,each=2)
polygon(c(xx,rev(xx)),c(yy.q5,rev(yy.q95)),col=NULL,border='darkblue',lwd=2)

####### Hi Window

mysnp=(dat.comb$left>5e6) & (dat.comb$right<6e6) &(dat.comb$chr==24)
tt=seq(5e6,6e6,2e5)
b=seq(5e6,6e6,1e5)
mc=mychr$m_cj[hi.win]

#### meiotic 50K
plot(dat.comb$left[mysnp],dat.comb$c[mysnp],type='n',ylab='c (cM/Mb)',
     xlab = "Position (Mb)", main = "c on 50K",
     ylim=c(0,9),axes=F,xlim=c(5e6,6e6))
abline(v=b, col="grey", lwd = 1, lty = 2)
points(dat.comb$left[mysnp],dat.comb$c[mysnp],type='s',
       pch = 16, lwd = 3, col='darkred')
axis(2, lwd = 2)
axis(1,at=tt,labels=round(tt*1e-6,digits=1),lwd=2)

#### LD based 50K
plot(dat.comb$left[mysnp],dat.comb$rho[mysnp],type='n',
     main = expression(paste(rho," ", "on 50K")),
     xlab = "Position (Mb)", ylab = expression(paste(rho, "(/kb)")),
     ylim=c(0,1.5)*1e-3,axes=F,xlim=c(5e6,6e6))
abline(v=b, col="grey", lwd = 1, lty = 2)
points(dat.comb$left[mysnp],dat.comb$rho[mysnp],type='s',
       pch = 16, lwd = 3, col='darkred')
axis(2, lwd = 2,at=seq(0,1.5,0.5)*1e-3,labels=seq(0,1.5,0.5))
axis(1,at=tt,labels=round(tt*1e-6,digits=1),lwd=2)

##### LD based 600K
mysnp=(pop.hd$chr==24) & (pop.hd$left>5e6) & (pop.hd$right<6e6)
plot(pop.hd$left[mysnp],pop.hd$rho[mysnp],type='n',
     main = expression(paste(rho," ", "on 600K")),
     xlab = "Position (Mb)", ylab = expression(paste(rho, "(/kb)")),
     ylim=c(0,7),axes=F,xlim=c(5e6,6e6))
abline(v=b, col="grey", lwd = 1, lty = 2)
points(pop.hd$left[mysnp],pop.hd$rho[mysnp],type='s',
       pch = 16, lwd = 3, col='darkred')
axis(2, lwd = 2,at=seq(0,7,1))
axis(1,at=tt,labels=round(tt*1e-6,digits=1),lwd=2)

####### Lo Window

mysnp=(dat.comb$left>18e6) & (dat.comb$right<19.1e6) &(dat.comb$chr==24)
tt=seq(18e6,19e6,2e5)
b=seq(18e6,19e6,1e5)
mc=mychr$m_cj[low.win]

#### meiotic 50K
plot(dat.comb$left[mysnp],dat.comb$c[mysnp],type='n',ylab='c (cM/Mb)',
     xlab = "Position (Mb)", main = "c on 50K",
     ylim=c(0,9),axes=F,xlim=c(18e6,19e6))
abline(v=b, col="grey", lwd = 1, lty = 2)
points(dat.comb$left[mysnp],dat.comb$c[mysnp],type='s',
       pch = 16, lwd = 3, col='darkblue')
axis(2, lwd = 2)
axis(1,at=tt,labels=round(tt*1e-6,digits=1),lwd=2)

#### LD based 50K
plot(dat.comb$left[mysnp],dat.comb$rho[mysnp],type='n',
     main = expression(paste(rho," ", "on 50K")),
     xlab = "Position (Mb)", ylab = expression(paste(rho, "(/kb)")),
     ylim=c(0,1.5)*1e-3,axes=F,xlim=c(18e6,19e6))
abline(v=b, col="grey", lwd = 1, lty = 2)
points(dat.comb$left[mysnp],dat.comb$rho[mysnp],type='s',
       pch = 16, lwd = 3, col='darkblue')
axis(2, lwd = 2,at=seq(0,1.5,0.5)*1e-3,labels=seq(0,1.5,0.5))
axis(1,at=tt,labels=round(tt*1e-6,digits=1),lwd=2)

##### LD based 600K
mysnp=(pop.hd$chr==24) & (pop.hd$left>18e6) & (pop.hd$right<19e6)
plot(pop.hd$left[mysnp],pop.hd$rho[mysnp],type='n',
     main = expression(paste(rho," ", "on 600K")),
     xlab = "Position (Mb)", ylab = expression(paste(rho, "(/kb)")),
     ylim=c(0,7),axes=F,xlim=c(18e6,19e6))
abline(v=b, col="grey", lwd = 1, lty = 2)
points(pop.hd$left[mysnp],pop.hd$rho[mysnp],type='s',
       pch = 16, lwd = 3, col='darkblue')
axis(2, lwd = 2,at=seq(0,7,1))
axis(1,at=tt,labels=round(tt*1e-6,digits=1),lwd=2)
dev.off()
