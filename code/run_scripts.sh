mkdir -p data/population/HDpanel
mkdir -p data/family/Lacaunes
## Get Panel data
baseurl=https://www.zenodo.org/record/237116/files
wget -P data/population/HDpanel $baseurl/frenchsheep_HD.bim
wget -P data/population/HDpanel $baseurl/frenchsheep_HD.bed
wget -P data/population/HDpanel $baseurl/frenchsheep_HD.fam
## Get Pedigree data
baseurl=https://www.zenodo.org/record/804264/files
wget -P data/family/Lacaunes $baseurl/lacaune_genotypes.fam
wget -P data/family/Lacaunes $baseurl/lacaune_genotypes.bed
wget -P data/family/Lacaunes $baseurl/lacaune_genotypes.bim

python code/family_map.py

python code/make_pop_map.py --dir data/population/PHASE/ --pad 1e5

python code/get_ld_rec_samples.py

python code/list_meioses.py

python code/run_aireml.py
