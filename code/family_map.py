#/usr/bin/env python2.7
# -*-coding:utf-8 -*
from __future__ import print_function
import sys
import random
import numpy as np
from scipy.stats import gamma,lognorm,binom
from scipy.special import gamma as gammaf
from scipy.special import gammaln as lgammaf

import matplotlib.pyplot as plt

''' 
Recombination rates estimation for each marker interval on the map from LINKPHASE output
'''

class Parent():
    '''
    Class for storing information for each parent

    Attributes:
    
    name : identifier for the parent
    meioses : dictionary with keys : offspring name, value : list of crossing over objects
    '''
    
    def __init__(self,name): 
        self.name = name
        self.meioses = {}
    def get_offspring_CO(self, name):
        '''
        Get the list of CO for offspring with name *name*
        If offspring is new creates a new empty list in meioses dictionary and returns it.
        '''
        try:
            ## search first, avoiding overwriting values in dict.
            L = self.meioses[name] 
        except KeyError:
            self.meioses[name] = []
            L=self.meioses[name]
        return L
    def add_offspring_CO(self, name, chro, left, right):
        '''
        Add a crossing over in offspring *name* on chromosome chro between *left* and *right*
        if offspring is new, will create new list of CO for this offspring in the meioses dictionary
        '''
        the_meiosis=self.get_offspring_CO(name)
        myco = CrossingOver(chro, left, right)
        the_meiosis.append(myco)
            
    def nb_CO_meioses(self):
        ''' 
        Counts number of crossing over for each meiosis
        '''
        res = [len(v) for v in self.meioses.values()]
        return res

    def nb_CO_tot(self):
        '''
        Counts total number of crossing over in all meioses
        '''
        return np.sum(self.nb_CO_meioses())
    
    def oi_xi(self, chrom, w_left, w_right):
        '''
        Computes probabilities that each crossing over in the parent occurs in
        genomic region on chromosome *chrom* between positions *w_gauche* and *w_droite*.
        
        Returns a tuple with entries:
        -- list of contributions for each CO
        -- number of meioses for the parent
        '''
        contrib = []
        for m in self.meioses.values():
            for co in m:
                if co.chro == chrom and not( (w_right < co.left ) or (w_left > co.right)):
                    contrib += [co.oi_xi(w_left, w_right)]
        return (contrib, len(self.meioses))

        
class CrossingOver():
    '''
    Class to store crossing over information
    
    Attributes:
    -- chro : chromosome
    -- left : position of the marker on the left side
    -- right : position of the marker on the right side
    '''
    def __init__(self, chro, gauche, droite):
        self.chro = chro
        self.left = gauche
        self.right = droite
                
    def oi_xi(self, w_left, w_right):
        '''
        Computes the probability that the crossing over occured between in the window between w_left and w_right
        '''
        if (w_right < self.left ) or (w_left > self.right):
            ## no overlap
            ## wl------wr               wl----------wr
            ##             sl-------sr
            return 0
        elif (w_left <= self.left) and (self.right <= w_right):
            ## co included in window
            ## wl---------------------------wr
            ##      sl---------------sr
            return 1
        elif (self.left <= w_left) and (w_right <= self.right):
            ## window is included in co
            ##          wl------wr
            ## sl---------------------sr
            return float(w_right - w_left)/(self.right-self.left)
        elif (w_left < self.left):
            ## we know w_right<self.right as other case is treated above
            ## wl-----------------wr
            ##       sl------------------sr
            return float(w_right-self.left)/(self.right-self.left)
        else:
            ## only case left
            ##           wl-----------------wr
            ##    sl--------------sr
            try: 
                assert (self.left <= w_left) and (self.right < w_right)
            except AssertionError:
                print(self.right, w_right, self.left, w_left)
                raise
            return float(self.right-w_left)/(self.right-self.left)
                
def main():
    ## set working directory
    main_wd='data/family/' 
    ## dictionary of parents
    pere = {}
    ## chromosome loop
    physmap=[]
    for i in range(1,27):
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

    ## Compute Ri and \bar(Ri)
    Ribar=0
    S_Ri_mi=0
    ## get selected sires
    FIDs={}
    with open(main_wd+"FIDs.txt") as f:
        for i,ligne in enumerate(f):
            buf=ligne.split()
            FIDs[buf[0]]=1
    with open("results/family/parent_recombination.txt",'w') as fout:
        print('parent','nbCO','nbMeio','Ri',file=fout)
        todel=[]
        for nom,p in pere.items():
            try:
                myfid=FIDs[p.name]
                p.Ri=float(p.nb_CO_tot())/len(p.meioses)
                print(p.name, p.nb_CO_tot(),len(p.meioses),p.Ri,file=fout)
                Ribar+=p.Ri
                S_Ri_mi+=p.nb_CO_tot() ## = Rs x ms
            except KeyError:
                todel.append(nom)
        for p in todel:
            del pere[p]
    Ribar/=len(pere)
    S_Ri_mi/=Ribar
    calibrate_prior=False
    if calibrate_prior:
        print('CALIBRATE PRIOR on c')
        ## Calibrate genome-wide prior on w_size windows
        cgw=[]
        w_size=1e5
        g_size=0
        for i in range(26):
            chrom=i+1
            print('CHROMOSOME',chrom)
            xlim=max(physmap[i].values())                ## right-most position
            g_size+=xlim
            for w_left in range(0,xlim,int(w_size)):
                pc=100.0*w_left/xlim
                sys.stdout.write("%3.2f %% \r"%pc)
                sys.stdout.flush()
                w_right=w_left+w_size
                M=0.0 ## total number of meioses
                oi_xi=0.0 ## expected number of CO in interval
                for p in pere.values():
                    res=p.oi_xi(chrom,w_left,w_right)
                    oi_xi+=np.sum(res[0])
                    M+=res[1]
                cgw.append(1e8*oi_xi/(w_size*M)) ## in cM/Mb
        print("Genome Size",g_size)
        cgw=np.array(cgw)
        cgw=cgw[cgw>0]
        prior_params=gamma.fit(cgw,floc=0)
        cj_alpha=prior_params[0]
        cj_beta=1.0/prior_params[2]
        cj_prior=gamma(prior_params[0],0,prior_params[2])
        hh=plt.hist(cgw,normed=1,bins=50, facecolor='g',alpha=0.75)
        ii=np.arange(0.0,max(cgw),0.01)
        ligne,=plt.plot(ii,cj_prior.pdf(ii),'r',linewidth=3)
        plt.legend([ligne],[ 'Gamma( '+str(np.round(cj_alpha,decimals=1))+', '+str(np.round(cj_beta,decimals=1))+')'])
        plt.title("Distribution of Raw c estimates")
        plt.savefig('results/family/Prior_cj.pdf')
        plt.close()
    else:
        cj_alpha=2.7
        cj_beta=1.8
        cj_prior=gamma(cj_alpha,0,1.0/cj_beta)
    ## 1Mb map
    with open("results/family/1Mb_map.txt",'w') as fout,\
         open("results/combined/comb_1Mb_map_family.txt",'w') as fout_comb:
        print("CONSTRUCTING 1Mb RECOMBINATION MAPS")
        print('chr','left','right','m_cj','s_cj','q5_cj','q95_cj',file=fout)
        print('chr','left','right','method','rep','value',file=fout_comb)
        w_size_rec=1e6
        for i in range(26):
            chrom=i+1
            print('CHROMOSOME',chrom)
            xlim=max(physmap[i].values())                ## left-most position
            for w_left in range(0,xlim,int(w_size_rec)):
                pc=100.0*w_left/xlim
                sys.stdout.write("%3.2f %% \r"%pc)
                sys.stdout.flush()
                w_right=w_left+w_size_rec
                p_xi=[] ## expected number of CO in interval
                for p in pere.values():
                    res=p.oi_xi(chrom,w_left,w_right)
                    p_xi+=res[0]
                sxi=np.array([np.sum(np.random.binomial(1,p_xi,len(p_xi))) for i in range(1000)])
                cj_dist=gamma(cj_alpha+sxi,0,1.0/(cj_beta+S_Ri_mi*1e-8*w_size_rec)).rvs(len(sxi))
                pc=np.percentile(cj_dist,[5,95])
                print(chrom,w_left,w_right,np.average(cj_dist),np.std(cj_dist),pc[0],pc[1],file=fout)
                fout.flush()
                ## random samples from posterior distribution for combined maps 
                pseudo_obs=random.sample(cj_dist,20)
                for i,o in enumerate(pseudo_obs):
                    print(chrom,w_left,w_right,'family',i,o,file=fout_comb)
    ## 60K intervals map
    with open("results/family/SNP_array_map.txt",'w') as fout:
        print("CONSTRUCTING SNP ARRAY RECOMBINATION MAPS")
        print('chr','left','right','m_cj','s_cj','q5_cj','q95_cj',file=fout)
        peres=pere.values()
        Ri_mi={}
        for p in peres:
                Ri_mi[p.name]=p.Ri*len(p.meioses)/Ribar
        for i in range(26):
            snp_pos=sorted(physmap[i].values())
            xlim=max(snp_pos) 
            chrom=i+1
            print('CHROMOSOME',chrom)
            for i,w_left in enumerate(snp_pos[:-1]):
                w_right=snp_pos[i+1]
                pc=100.0*w_left/xlim
                w_size_rec=w_right-w_left
                sys.stdout.write("%3.2f %% \r"%pc)
                sys.stdout.flush()
                p_xi=[] ## list of P(CO occurs in interval) for each CO
                nco={}
                for p in peres:
                    res=p.oi_xi(chrom,w_left,w_right)
                    p_xi+=res[0]
                    nco[p.name]=len(res[0])
                ## get samples from the distribution of number of CO across meioses in the interval
                nsamp=1000
                sxi=np.zeros(nsamp,dtype=np.float64)
                for samp in range(nsamp):
                    oxi=np.random.binomial(1,p_xi,len(p_xi))
                    sxi[samp]=np.sum(oxi)
                ## get posterior samples of rec. rate giv. sxi
                cj_dist=gamma(cj_alpha+sxi,0,1.0/(cj_beta+S_Ri_mi*1e-8*w_size_rec)).rvs(len(sxi))
                ## empirical quantiles on rec. rate.
                pc=np.percentile(cj_dist,[5,95])
                ## print out summary stats on c_j
                print(chrom,w_left,w_right,np.average(cj_dist),np.std(cj_dist),pc[0],pc[1],file=fout)
                fout.flush()

if __name__ == '__main__':
    main()
