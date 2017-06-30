''' Reads a plink .fam files and renumber invididuals so that 
    older individuals have a smaller integer id than youger ones'''

infile='data/Soay/20150129merged1_66nodups.QC2.fam'
prefix=infile[:-4]

nullid='0'
curid=0

ancestors={nullid:curid}


## get data
rawindivs={}
fathers={}
mothers={}
with open(infile) as f:
   for ligne in f:
       buf=ligne.split()
       ## key = individ, value = fatherid,motherid
       rawindivs[buf[1]]=(buf[2],buf[3])
       try:
           fathers[buf[2]]+=1
       except KeyError:
           fathers[buf[2]]=1
       try:
           mothers[buf[3]]+=1
       except KeyError:
           mothers[buf[3]]=1
           
## remove individuals that only appear in fatherid, motherid
newindiv={}
for k,v in rawindivs.items():
    dad=v[0]
    mom=v[1]
    ## look for dad in individuals
    if v[0]!=nullid:
        try:
            dadhere=rawindivs[v[0]]
            dad=v[0]
        except KeyError:
            #print "Individual Father",v[0],"not in data"
            dad=nullid
    ## look for mom in individuals
    if v[1]!=nullid:
        try:
            momhere=rawindivs[v[1]]
            mom=v[1]
        except KeyError:
            #print "Individual Mother",v[1],"not in data"
            mom=nullid
    newindiv[k]=(dad,mom)

listindiv=dict(newindiv.items())
while len(listindiv)>0:
    nanc={}
    for k,v in listindiv.items():
        n=0
        try:
            dadid=ancestors[v[0]]
            n+=1
        except KeyError:
            pass
        try:
            momid=ancestors[v[1]]
            n+=1
        except KeyError:
            pass
        nanc[k]=n
    newanc=[k for k in nanc.keys() if nanc[k]==2]
    for i,k in enumerate(newanc):
        ancestors[k]=curid+i+1
        ##print k,ancestors[k]
        del listindiv[k]
    curid=curid+len(newanc)

## write corresp
with open(prefix+'.corresp','w') as fout:
    for k,v in ancestors.items():
        print >>fout,k,v

## identify individuals to remove
with open(prefix+'.rm','w') as fout:
    for indiv in ancestors.keys():
        if indiv is nullid:
            continue
        dad,mom=newindiv[indiv]
        if (dad==nullid) and (mom==nullid):
            try:
                tmp=fathers[indiv]
                isdad=True
            except KeyError:
                isdad=False
            try:
                tmp=mothers[indiv]
                ismom=True
            except KeyError:
                ismom=False
            if not (isdad or ismom):
                print >>fout,ancestors[indiv],ancestors[indiv],dad,mom
    
## now everyone should be numbered, we can write a new fam file with numbers
with open(infile) as f:
    with open(prefix+'.num.fam','w') as fout:
        for ligne in f:
            buf=ligne.split()
            indiv=ancestors[buf[0]]
            try:
                dad=ancestors[buf[2]]
            except KeyError:
                dad=ancestors[nullid]
            try:
                mom=ancestors[buf[3]]
            except KeyError:
                mom=ancestors[nullid]
            print>>fout,indiv,indiv,dad,mom,buf[4],buf[5]

