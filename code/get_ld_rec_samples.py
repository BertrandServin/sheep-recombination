#!/bin/env python
from __future__ import print_function
import os
import numpy as np
from make_pop_map import Marker_Interval

def get_PHASE_results(wd):
    marker_intervals={}
    for chrom in range(1,27):
        marker_intervals[chrom]=[]
        thedir=wd+'/'+str(chrom)
        for subdir in os.listdir(thedir):
            print(chrom,subdir)
            dat=np.genfromtxt('/'.join([thedir,subdir,'phase_output_recom']))
            pos=[int(x) for x in open('/'.join([thedir,subdir,'map.txt'])).readline().split()]
            rho=dat[1:,0]
            lbd_dist=dat[1:,1:]
            beg,end=[int(x) for x in subdir.split('-')]
            for ix,rg in enumerate(zip(pos[:-1],pos[1:])):
                g,d=rg
                mi=Marker_Interval(c=chrom,l=g,r=d)
                mi.set_lambda(lbd_dist[:,ix])
                mi.rho_dist=rho
                marker_intervals[chrom].append(mi)
    return marker_intervals

def main():
    mi=get_PHASE_results('./data/population/PHASE51')
    print(len(mi))
    nsamp=20
    one_mb_rec=np.genfromtxt('results/family/1Mb_map.txt',names=True)
    with open('results/combined/comb_1Mb_map_pop51.txt','w') as fout:
        print('chr','left','right','method','rep','value',file=fout)
        for win in one_mb_rec:
            try:
                win_mi=[x for x in mi[win['chr']] if x.left>=win['left'] and x.right<=win['right']]
            except KeyError:
                break
            delta_tot=np.zeros(nsamp,dtype=np.float64)
            for m in win_mi:
                delta_dist=m.get_delta_i_k()
                delta_tot+=np.random.choice(delta_dist,nsamp)
            delta_tot/=(win['right']-win['left'])
            delta_tot*=1e8
            for i,val in enumerate(delta_tot):
                print(win['chr'],win['left'],win['right'],'pop',i,val,file=fout)
                fout.flush()
            
if __name__=='__main__':
    main()
