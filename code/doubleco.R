carte=read.table('map')
colnames(carte)=c('idx','name','pos')
co=read.table('recombinations_hmm')
colnames(co)=c('off','sir','l','r')
co$dl=carte$pos[co$l]
co$dr=carte$pos[co$r]
## should already be ordered correctly but just in case ...
co=co[order(co$sir,co$off,co$l),]
nco=nrow(co)
## is the next observation from the same meiosis
co$meioadj=FALSE
co$meioadj[-nco]=(co$off[-1]==co$off[-nco])
## what is the distance to the next co
co$dadj=100
co$dadj[-nco]=co$dl[-1]-co$dr[-nco]

## identification of CO that are distant by less than 3cM from
## the next, for the same meiosis
doubleidx=which((co$dadj<3)&(co$meioadj))
print(doubleidx)
## Read in genotype data
dat=read.table('typ',row.names=1)
for (gidx in doubleidx) {
    iidx=which(rownames(dat)==co$off[gidx])
    lmk=2*co$l[gidx]-1
    rmk=2*co$r[gidx+1]
    dat[iidx,lmk:rmk]=0
}
write.table(dat,file='typ_cor',col.names=F,row.names=T,quote=F)
