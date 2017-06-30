#!/usr/bin/env python
''' 
Converts fastphase files created by plink to PHASE files
'''
import sys

def usage():
    print 'Usage:',sys.argv[0],'input_file'

    
def main():
    try:
        f=open(sys.argv[1])
        print 'Processing file',sys.argv[1]
    except:
        usage()
        sys.exit(1)
    fout=open(sys.argv[1]+'.phase','w')
    ## ligne nb individus
    ligne=f.readline()
    ni=int(ligne.split()[0])
    print >>fout,ligne,
    ## ligne nb marqueurs
    ligne=f.readline()
    ns=int(ligne.split()[0])
    print >>fout,ligne,
    ## ligne P ...
    print >>fout,f.readline(),
    ## Ajout de la ligne SSS
    print >>fout,'S'*ns
    ## reste
    data=f.readlines()
    for ind in range(ni):
        idl=data[3*ind]
        idh1=data[3*ind+1]
        idh2=data[3*ind+2]
        print >>fout,''.join([x for x in idl.split()])
        print >>fout,' '.join([c for c in idh1]),
        print >>fout,' '.join([c for c in idh2]),
    f.close()
    fout.close()
        
    
if __name__=='__main__':
    main()
