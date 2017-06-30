famdat=read.table('results/Soay/Johnston_2016/2_FamilyPedigree_FullClean_g.txt',h=T)
corresp=read.table('data/Soay/20150129merged1_66nodups.QC2.corresp')
id2keep=unique(c(as.character(famdat$ANIMAL),as.character(famdat$FATHER),as.character(famdat$MOTHER)))
numid2keep=corresp$V2[match(id2keep,corresp$V1)]
fam=read.table('data/Soay/Soays_autosomes.fam')
fam2keep=fam[match(numid2keep,fam$V2),]
write.table(fam2keep,'data/Soay/keep_indivs.txt',row.names=F,col.names=F,quote=F)
