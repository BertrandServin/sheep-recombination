## Get the Lacaune individuals from HD panel data
plink --sheep --family --keep-cluster-names LAC LAM --bfile data/population/HDpanel/frenchsheep_HD  --make-bed --out data/population/HDpanel/Lacaune
## remove individuals that are present in the cohort
plink --sheep --bfile data/population/HDpanel/Lacaune --remove data/population/HDpanel/shared_indivs.txt --make-bed --out data/population/HDpanel/Lacaune

## Running PHASE on a window (these needs to be automated to
## analyse the whole genome)
mkdir -p tmpres/PHASE
plink --bfile data/population/HDpanel/Lacaune --chr 1 --from-kb 0 --to-kb 2100 --recode 01 fastphase-1chr --out tmpres/PHASE/window
## convert fastPHASE file to PHASE
python code/fph2ph.py tmpres/PHASE/window.recode.phase.inp
## run PHASE to estimate recombination rates
PHASE -X10 tmpres/PHASE/window.recode.phase.inp.phase tmpres/PHASE/window
