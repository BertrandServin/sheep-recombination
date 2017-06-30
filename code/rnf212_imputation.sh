## exctract the end of OAR6
plink --sheep --bfile data/gwas/FinalCohort --chr 6 --from-kb 115000 \
      --out data/gwas/cohort-chr6qtl --make-bed
plink --sheep --bfile data/gwas/cohort-chr6qtl --bmerge data/gwas/rnf212mut \
      --make-bed --out data/gwas/cohort-chr6qtl-mut

## get panel
plink --sheep --family --bfile data/population/HDpanel/frenchsheep_HD \
  --keep-cluster-names LAC LAM \
  --chr 6 --from-kb 115 \
  --out data/gwas/panel-chr6qtl-mut --recode-bimbam

## Recode bimbam
plink --sheep --bfile data/gwas/cohort-chr6qtl-mut --chr 6 \
      --recode-bimbam --out data/gwas/cohort-chr6qtl-mut

## Run bimbam
bimbam -g data/gwas/panel-chr6qtl-mut.recode.geno.txt -p 0 \
       -pos data/gwas/panel-chr6qtl-mut.recode.pos.txt \
       -g data/gwas/cohort-chr6qtl-mut.recode.geno.txt \
       -p data/gwas/cohort-chr6qtl-mut.recode.pheno.txt\
       -pos data/gwas/cohort-chr6qtl-mut.recode.pos.txt  \
       -e 10 -w 20 -s 1 -c 15 -o imput_rnf212 -wmg 

mv output/imput_rnf212.mean.genotype.txt data/gwas/
awk '{print $1, $6 ,$5}' output/imput_rnf212.snpinfo.txt  | sed -e s/0$/6/g >data/gwas/imput_rnf212.snpinfo.txt
