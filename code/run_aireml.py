import os
import tempfile
from subprocess import Popen, PIPE, STDOUT


months=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

## create input file for renumf90
with open('results/family/nco_meioses.txt') as f:
    fpheno=tempfile.NamedTemporaryFile(suffix='.pheno',delete=False)
    print 'Phenotype File:',fpheno.name
    for ligne in f:
        buf=ligne.split()
        try: 
            print >>fpheno,buf[1],buf[1],1,buf[3],buf[4]
        except ValueError:
            continue
    fpheno.close()

## create parameter file for renumf90
fpar=tempfile.NamedTemporaryFile(suffix='.par',delete=False)
print 'Parameter File:',fpar.name

print >>fpar,'DATAFILE'
print >>fpar,fpheno.name
print >>fpar,'TRAITS'
print >>fpar,'5'
print >>fpar,'FIELDS_PASSED TO OUTPUT'
print >>fpar,'1'
print >>fpar,'WEIGHT(S)\n'
print >>fpar,'RESIDUAL_VARIANCE'
print >>fpar,'24.0'
print >>fpar,'EFFECT'
print >>fpar,'3 cross numer'
print >>fpar,'EFFECT'
print >>fpar,'4 cross numer'
print >>fpar,'EFFECT'
print >>fpar,'2 cross alpha'
print >>fpar,'RANDOM'
print >>fpar,'diagonal'
print >>fpar,'EFFECT'
print >>fpar,'1 cross alpha'
print >>fpar,'RANDOM'
print >>fpar,'animal'
print >>fpar,'FILE'
print >>fpar,'data/family/lacaune_pedigree_4G.txt'
print >>fpar,'PED_DEPTH'
print >>fpar,'4'
print >>fpar,'(CO)VARIANCES'
print >>fpar,'4.0000'
print >>fpar,'OPTION sol se\n'
fpar.close()

## Run renumf90
p=Popen(['renumf90'],stdout=PIPE,stdin=PIPE,stderr=STDOUT)
renumf90_out=p.communicate(input=fpar.name)[0]
with open('results/family/renumf90.out','w') as f:
    print >>f,renumf90_out

## Get heritability estimates
with open('renf90.par','a') as frenf90:
    print>>frenf90,"OPTION se_covar_function h2 g_4_4_1_1/(g_4_4_1_1+g_3_3_1_1+r_1_1)"

## run airemlf90
p=Popen(['airemlf90'],stdout=PIPE,stdin=PIPE,stderr=STDOUT)
aireml_out=p.communicate(input=b'renf90.par')[0]
with open('results/family/aireml.out','w') as f:
    print >>f,aireml_out

## store results and clean up files
os.rename('solutions','results/family/solutions.aireml')
os.rename('renadd04.ped','results/family/renadd04.ped')

os.remove('renf90.par')
os.remove('renf90.dat')
os.remove('renf90.tables')
os.remove('fspak90.ord')
os.remove('airemlf90.log')
os.remove('fort.8765')
os.remove(fpar.name)
os.remove(fpheno.name)
