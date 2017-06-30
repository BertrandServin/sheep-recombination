import os
import tempfile
import subprocess

## extract the 345 sires from Lacaune genotypes (store BLUP as phenotype)
cohort={}
fsel=tempfile.NamedTemporaryFile(suffix='.sel',delete=False)
with open('results/family/parent_recombination_blup.txt') as f:
    for i,ligne in enumerate(f):
        if i<1:
            continue
        buf=ligne.split()
        cohort[buf[0]]=buf[-1]
        print >>fsel,10,buf[0]
fsel.close()
print 'Selected animals file:',fsel.name

plink_cmd=['plink','--sheep']
plink_cmd+=['--bfile','data/family/Lacaunes/lacaune_genotypes']
plink_cmd+=['--keep',fsel.name]
plink_cmd+=['--make-bed','--out','data/gwas/cohort']

print '*** Running',' '.join(plink_cmd)
subprocess.check_call(plink_cmd)
os.remove(fsel.name)

## Some individuals are in the HD panel, create a cohort file with HD genotypes
plink_cmd=['plink','--sheep']
plink_cmd+=['--bfile','data/population/HDpanel/frenchsheep_HD']
plink_cmd+=['--keep','data/population/HDpanel/shared_indivs.txt']
plink_cmd+=['--make-bed','--out','data/gwas/HDcohort']
print '*** Running',' '.join(plink_cmd)
subprocess.check_call(plink_cmd)
os.rename('data/gwas/HDcohort.fam','data/gwas/HDcohort.fam.orig')

corresp={}
with open('data/population/HDpanel/shared_indivs.txt') as f:
    for ligne in f:
        buf=ligne.split()
        corresp[buf[1]]=buf[2]

with open('data/gwas/HDcohort.fam.orig') as f:
    with open('data/gwas/HDcohort.fam','w') as fout:
        for ligne in f:
            buf=ligne.split()
            nom=corresp[buf[1]]
            buf[0]='10'
            buf[1]=nom
            print >>fout,' '.join(buf)

## merge datasets
plink_cmd=['plink','--sheep']
plink_cmd+=['--bfile','data/gwas/cohort']
plink_cmd+=['--bmerge','data/gwas/HDcohort']
plink_cmd+=['--out','data/gwas/FinalCohort']
print '*** Running',' '.join(plink_cmd)
subprocess.check_call(plink_cmd)

## run BIMBAM
jobfile=open('bimbamjobs.sh','w')
for chrom in range(1,27):
    ## create cohort file
    plink_cmd=['plink','--sheep']
    plink_cmd+=['--bfile','data/gwas/FinalCohort']
    plink_cmd+=['--chr',str(chrom)]
    plink_cmd+=['--recode','bimbam']
    plink_cmd+=['--out','data/gwas/cohort_chr'+str(chrom)]
    print '*** Running',' '.join(plink_cmd)
    subprocess.check_call(plink_cmd)

    ## create panel file
    plink_cmd=['plink','--sheep']
    plink_cmd+=['--family']
    plink_cmd+=['--bfile','data/population/HDpanel/frenchsheep_HD']
    plink_cmd+=['--keep-cluster-names','LAC','LAM']
    plink_cmd+=['--chr',str(chrom)]
    plink_cmd+=['--recode','bimbam']
    plink_cmd+=['--out','data/gwas/panel_chr'+str(chrom)]
    print '*** Running',' '.join(plink_cmd)
    subprocess.check_call(plink_cmd)

    ## Impute
    bb_cmd=['bimbam']
    bb_cmd+=['-g','data/gwas/panel_chr'+str(chrom)+'.recode.geno.txt','-p','0']
    bb_cmd+=['-g','data/gwas/cohort_chr'+str(chrom)+'.recode.geno.txt','-p','data/gwas/cohort_chr'+str(chrom)+'.recode.pheno.txt']
    bb_cmd+=['-pos','data/gwas/cohort_chr'+str(chrom)+'.recode.pos.txt']
    bb_cmd+=['-e','10','-w','20','-s','1','-c','15','-o','imput_lacaune_chr'+str(chrom),'-wmg','-wgd']
    ## To run the analysis uncomment the following
    # print '*** Running',' '.join(bb_cmd)
    # subprocess.check_call(bb_cmd)
    ## print jobs to file
    print >>jobfile,' '.join(bb_cmd)
