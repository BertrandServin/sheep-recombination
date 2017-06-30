ped=read.table('data/family/Lacaunes/lacaune_genotypes.fam')
ped$V2=as.factor(ped$V2)
ped$V3=as.factor(ped$V3)
ped$V4=as.factor(ped$V4)
## remove individuals with no sex information
ped=ped[ped$V5!=0,]

indivs=as.data.frame(levels(ped$V2))
colnames(indivs)=c('code')
indivs$idx=match(indivs$code,ped$V2)
indivs$isfather=sapply(indivs$code,function(x){ x%in%ped$V3 })
indivs$ismother=sapply(indivs$code,function(x){ x%in%ped$V4 })
indivs$sex=ped$V5[indivs$idx]

## father / mother
peres=ped$V3[indivs$idx]
indivs$hasfather=(peres!='0')&(peres%in%indivs$code)
meres=ped$V4[indivs$idx]
indivs$hasmother=(meres!='0')&(meres%in%indivs$code)

## number of offspring
indivs$noff=sapply(indivs$code,function(x) { sum(ped$V3==as.character(x))+sum(ped$V4==as.character(x))})

## Keep male individuals with known father and at least two offspring
sel1=indivs[(indivs$hasfather)&(indivs$noff>1)&(indivs$sex==1),]

## keep male individuals with unknown father and at least four offspring
sel2=indivs[(!indivs$hasfather)&(indivs$noff>3)&(indivs$sex==1),]

sel=rbind(sel1,sel2)
write.table(sel$code,quote=F,col.names=F,row.names=F,file='data/family/FIDs.txt')

