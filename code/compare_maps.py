#!/usr/bin/env python2.8
from __future__ import print_function
import numpy as np
import sys

def main():
    hdrec=np.genfromtxt('results/population/HD_SNP_array_map.txt',names=True)

    one_mb_rec=np.genfromtxt('results/family/1Mb_map.txt',names=True)
    with open('results/combined/compare_family_1Mb.txt','w') as fout:
        print('chr','left','right','rho','srho','c','sc',file=fout)
        for win in one_mb_rec:
            print(win['chr'],win['left'],win['right'],sep='-',end='\r')
            sys.stdout.flush()
            sub=(hdrec['chr']==win['chr']) & (hdrec['left']>=win['left']) & (hdrec['right']<=win['right'])
            hdwin=hdrec[sub]
            l=np.sum(hdwin['right']-hdwin['left'])
            d=np.sum(hdwin['delta'])
            if l>0:
                rho=d/l
                sd_rho=np.sqrt(np.sum(hdwin['s_delta']**2))/l
                print(win['chr'],win['left'],win['right'],rho,sd_rho,win['m_cj'],win['s_cj'],file=fout)

    print()
    SNP_rec=np.genfromtxt('results/family/SNP_array_map.txt',names=True)
    with open('results/combined/compare_family_60K.txt','w') as fout:
        print('chr','left','right','rho','srho','c','sc',file=fout)
        for win in SNP_rec:
            print(win['chr'],win['left'],win['right'],sep='-',end='\r')
            sys.stdout.flush()
            sub=(hdrec['chr']==win['chr']) & (hdrec['left']>=win['left']) & (hdrec['right']<=win['right'])
            hdwin=hdrec[sub]
            l=np.sum(hdwin['right']-hdwin['left'])
            d=np.sum(hdwin['delta'])
            if l>0:
                rho=d/l
                sd_rho=np.sqrt(np.sum(hdwin['s_delta']**2))/l
                print(win['chr'],win['left'],win['right'],rho,sd_rho,win['m_cj'],win['s_cj'],file=fout)
            
    print()
if __name__=='__main__':
    main()
