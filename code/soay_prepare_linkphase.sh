## number animals using integers for linkphase
## this creates several files:
## -- data/Soay/20150129merged1_66nodups.QC2.num.fam: new fam file with integer codes replacing chars
## -- data/Soay/20150129merged1_66nodups.QC2.rm: individuals with no parents nor offspring => to discard
## -- data/Soay/20150129merged1_66nodups.QC2.corresp: correspondance between original IDs and linkphase IDs
python code/soay_renumfam.py

## recode Soay plink files with Linkphase IDs
plink --sheep --chr 1-26 --bed data/Soay/20150129merged1_66nodups.QC2.bed --fam data/Soay/20150129merged1_66nodups.QC2.num.fam --bim data/Soay/20150129merged1_66nodups.QC2.bim --make-bed --out data/Soay/Soays_autosomes

## remove individuals with no offspring and no parents in the data
plink --sheep --bfile data/Soay/Soays_autosomes --remove data/Soay/20150129merged1_66nodups.QC2.rm --make-bed --out data/Soay/Soays_autosomes

## keep the same individuals as those use in the Soay study
plink --sheep --bfile data/Soay/Soays_autosomes --keep data/Soay/keep_indivs.txt --make-bed --out data/Soay/Soays_autosomes_clean
