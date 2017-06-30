from family_map import Parent, CrossingOver

def main():
    main_wd='data/family/' 
    ## dictionary of parents
    pere = {}
    ## Get CO from LINKPHASE
    physmap=[]
    for i in range(1,27):
        print i
        ## get marker information
        carte_chr={} 
        with open(main_wd+'LINKPHASE/chr_'+str(i)+'/map','r') as mapfile:
            for ligne in mapfile:
                buf=ligne.split()
                carte_chr[buf[0]]=int(float(buf[2])*1e6) ## converts LINKPHASE input in Mb to bp.
        physmap.append(carte_chr)
        ## get crossing overs
        with open(main_wd+'LINKPHASE/chr_'+str(i)+'/recombinations_hmm','r') as f:
            for ligne in f:
                buf = ligne.split()
                try:
                    lepere = pere[buf[1]]
                except KeyError:
                    lepere = Parent(buf[1]) 
                    pere[buf[1]] = lepere
                lepere.add_offspring_CO(buf[0], i, carte_chr[buf[2]], carte_chr[buf[3]])
    ## get data on yob and insemination month
    covar={}
    with open('data/family/sheep_covariates.txt') as f:
        for i,ligne in enumerate(f):
            if i<1:
                continue
            buf=ligne.split()
            ## yob,insem
            covar[buf[0]]=(buf[1],buf[2])
            
    ## get selected sires
    FIDs={}
    with open(main_wd+"FIDs.txt") as f:
        for i,ligne in enumerate(f):
            buf=ligne.split()
            FIDs[buf[0]]=1

    with open('results/family/nco_meioses.txt','w') as fout:
        ##print >>fout,'sire','offspring','year','insem','nco'
        for nom in FIDs:
            try:
                p=pere[nom]
            except KeyError:
                print nom,'Not Found'
                raise
            for m in p.meioses:
                print>>fout,m,nom,covar[nom][0],covar[m][1],len(p.meioses[m])

            
if __name__ == '__main__':
    main()
