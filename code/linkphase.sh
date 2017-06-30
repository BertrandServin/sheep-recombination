#!/bin/bash
# -*-coding:utf-8 -* 


#Creation of a directory for each chromosome
#Add each file used by LINKPHASE : pedigree, SNPs, genotypes
basedir=`pwd`
datadir='data/family/Lacaunes'
tmpdir='data/family/linkphase.in'
mkdir -p $tmpdir

#Ped file creation

awk '{print $2, $3, $4}' $datadir/lacaune_genotypes.fam > $tmpdir/ped


#Creation of the files for each chromosome
for ((chr=1 ; chr<=26 ; chr++)) 
do
    echo "CHROMOSOME" $chr
    wdir=$tmpdir/chr_$chr
    mkdir -p $wdir
    #Genotype file creation in each directory
    #Genotypes conversion in "AG" to "12"
    plink --sheep --bfile $datadir/lacaune_genotypes --chr $chr --recode 12 --out $wdir/typage
    #Format to linkphase input
    cat $wdir/typage.ped| awk '{$1=$3=$4=$5=$6=""; print $0}' | tr -s ' ' | sed -e 's/^[ \t]*//g' > $wdir/typ

    #SNPs file creation in each directory
    #Elimination of useless informations : cM size, chromosome number
    awk -v r=1000000 '{printf  "%d %s %.6f\n", NR, $2, $4/r}' $wdir/typage.map >  $wdir/map
    ln -s ../ped $wdir/
    ## create parameder file
    echo "#PEDIGREE_FILE" > $wdir/linkin.txt
    echo "ped" >> $wdir/linkin.txt

    echo "#GENOTYPE_FILE" >> $wdir/linkin.txt
    echo "typ" >> $wdir/linkin.txt
    
    echo "#MARKER_FILE" >> $wdir/linkin.txt
    echo "map" >> $wdir/linkin.txt
    
    echo "#HALFSIB_PHASING" >> $wdir/linkin.txt
    echo "yes" >> $wdir/linkin.txt
    
    echo "#HMM_PHASING" >> $wdir/linkin.txt
    echo "yes" >> $wdir/linkin.txt
    
    echo "#N_TEMPLATES" >> $wdir/linkin.txt
    echo "50" >> $wdir/linkin.txt
    
    echo "#CHECK_PREPHASING" >> $wdir/linkin.txt
    echo "yes" >> $wdir/linkin.txt
    ## Get in directory to run linkphase analysis
    cd $wdir
    echo "First Run of LINKPHASE"
    ## run linkphase once
    linkphase | tee run1.log
    echo "Identify double crossovers"
    ## identify double crossovers
    R CMD BATCH $basedir/code/doubleco.R
    mv recombinations_hmm recombinations_hmm.run1
    ## run linkphase again
    echo "Rerun LINKPHASE"
    mv typ typ.1
    mv typ_cor typ
    linkphase | tee run2.log
    sleep 1
    mv recombinations_hmm recombinations_hmm.run2
    ## Test if we improved anything
    nco1=`wc recombinations_hmm.run1 | tr -s ' '| cut -f 2 -d ' '`
    nco2=`wc recombinations_hmm.run2 | tr -s ' '| cut -f 2 -d ' '`
    if (( nco1 < nco2 ))
    then
        echo "Not improved: " $nco1 " " $nco2
        cp recombinations_hmm.run1 recombinations_hmm
    else
        echo "Improved: " $nco1 " " $nco2
        cp recombinations_hmm.run2 recombinations_hmm
    fi
    cd $basedir
    echo "DONE"
done
## Keep only relevant files for the analysis
mkdir -p data/family/LINKPHASE
rsync -av --include="*/" --include="map" --include="recombinations_hmm" --include="ped" --exclude="*" data/family/linkphase.in/ data/family/LINKPHASE/

 
