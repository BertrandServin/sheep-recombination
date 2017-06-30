#!/bin/env python
from __future__ import print_function
import os
import argparse
import numpy as np

class Marker_Interval(object):
    
    def __init__(self,c,l,r):
        '''
        Create a Marker interval on chromosome c, from l to r.
        '''
        self.chrom=c
        if l>r or l<0 or r<0:
            raise ValueError
        self.left=l
        self.right=r
        self.length=self.right-self.left
        self.lambda_hat=None
        self.lbd_dist=None
        self.rho_dist=None
        
    def set_lambda(self,lbd_dist,estimator=np.median):
        self.lbd_dist=lbd_dist
        self.lambda_hat=estimator(lbd_dist)
        
    def get_rho_hat(self,rho):
        '''
        Calculates the scaled recombination rate:
        rho_j = rho * lambda_j * l
        where rho is assumed to be 4Nc
        '''
        return self.lambda_hat*self.length*rho

    def get_delta_i_k(self,rhodist=None):
        if rhodist is not None:
            self.rho_dist=rhodist
        return self.lbd_dist*self.rho_dist*self.length
    
    def get_delta_i(self,rhodist,estimator=np.median):
        delta_i_k=self.get_delta_i_k(rhodist)
        return estimator(delta_i_k),np.std(delta_i_k),np.percentile(delta_i_k,[5,95])
    
class Chromosome_Window(object):
    
    def __init__(self,c,l,r,pad=0):
        self.chrom=c
        if l>0:
            self.left=l+pad
        else:
            self.left=0
        self.right=r-pad
        self.mk_int=[]
        self.rho=None
        self.rhodist=None
        
    def add_mk_interval(self,mi):
        if mi.chrom==self.chrom:
            # if mi.left < self.left and mi.right > self.left:
            #     ## intervalle chevauchant a gauche
            #     #self.left=mi.left
            #     self.mk_int.append(mi)
            if mi.left < self.right and mi.right > self.right:
                ## overlapping interval on the right
                self.mk_int.append(mi)
            elif mi.left >= self.left and mi.right<=self.right:
                ## interval is included in window
                self.mk_int.append(mi)
    def true_left(self):
        return np.min([x.left for x in self.mk_int])
    def true_right(self):
        return np.max([x.right for x in self.mk_int])
    def set_rho(self,rho_dist,estimator=np.median):
        self.rhodist=rho_dist
        self.rho=estimator(rho_dist)

    def length(self):
        return np.sum([x.length for x in self.mk_int])
    
    def get_rho_hat(self):
        '''
        Sum_j rho lambda_j l_j
        '''
        return np.sum([x.get_rho_hat(self.rho) for x in self.mk_int])

    def get_local_rho(self):
        '''
        (1/L) Sum_j rho lambda_j l_j
        '''
        return self.get_rho_hat()/self.length()
    
    def get_mcmc_rho(self):
        myvec=np.array([x.get_delta_i_k(self.rhodist) for x in self.mk_int])
        ## myvec is an array of shape N_MK_INT x N_MCMC_ITER
        rho_w_dot=np.sum(myvec,axis=0)/self.length()
        return np.mean(rho_w_dot)
    
 
def main():
    
    parser=argparse.ArgumentParser()
    parser.add_argument('--dir',dest='wd')
    parser.add_argument('--pad',dest='pad',help='padding length',metavar='L',type=float)
    myopts=parser.parse_args()
    with open('results/population/HD_SNP_array_map51.txt','w') as fout:
        print('chr','left','right','lambda','delta','s_delta','q5_delta','q95_delta',file=fout)
        for chrom in range(1,27):
            thedir=myopts.wd+'/'+str(chrom)
            windows=[]
            for subdir in os.listdir(thedir):
                dat=np.genfromtxt('/'.join([thedir,subdir,'phase_output_recom']))
                ## pos=dat[0,:] SNP pos are truncated in PHASE output
                pos=[int(x) for x in open('/'.join([thedir,subdir,'map.txt'])).readline().split()]
                rho=dat[1:,0]
                lbd_dist=dat[1:,1:]
                beg,end=[int(x) for x in subdir.split('-')]
                thewindow=Chromosome_Window(chrom,beg,end,myopts.pad)
                thewindow.set_rho(rho)
                for ix,rg in enumerate(zip(pos[:-1],pos[1:])):
                    g,d=rg
                    mi=Marker_Interval(c=chrom,l=g,r=d)
                    mi.set_lambda(lbd_dist[:,ix])
                    thewindow.add_mk_interval(mi)
                windows.append(thewindow)
            for win in sorted(windows,key=lambda x: x.left):
                if len(win.mk_int) > 0:
                    for m in win.mk_int:
                        try:
                            mu,sd,pc=m.get_delta_i(win.rhodist)
                            print(m.chrom,int(m.left),int(m.right),m.lambda_hat,mu,sd,pc[0],pc[1],file=fout)
                            fout.flush()
                        except ValueError:
                            print(win.chrom,win.left,win.right,len(win.mk_int))
                            raise


if __name__=='__main__':
    main()


