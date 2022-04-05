## Running capice using VIP singularity containers
Hardcoded paths in sh files need to be replaced with actual ones. vip/images/bcftools-1.14.sif and 
vip/images/vep-105.0.sif should be present in order to run.

Run vep
```commandline
sbatch -t 60 --mincpus 8 --mem 8GB /groups/solve-rd/tmp10/mslofstra/run_vep_capice_p37.sh
sbatch -t 60 --mincpus 8 --mem 8GB /groups/solve-rd/tmp10/mslofstra/run_vep_capice_p38.sh
```

Run bcftools
```commandline
chmod u+rwx run_bcftools_capice_p37.sh
chmod u+rwx run_bcftools_capice_p38.sh
./run_bcftools_capice_p37.sh
./run_bcftools_capice_p38.sh
```

Get capice v3.0.0 and models:
```commandline
wget https://github.com/molgenis/capice/archive/refs/tags/v3.0.0.zip
wget https://github.com/molgenis/capice/releases/download/v3.0.0/v3.0.0-v1_grch37.pickle.dat
wget https://github.com/molgenis/capice/releases/download/v3.0.0/v3.0.0-v1_grch38.pickle.dat
```
Install capice:
```commandline
unzip v3.0.0.zip
cd capice-3.0.0/
chmod u+rwx setup.py
ml Python3.9
pip install .
```

Run:
```commandline
capice predict -i capice_input_p37.tsv -m v3.0.0-v1_grch37.pickle.dat
capice predict -i capice_input_p38.tsv -m v3.0.0-v1_grch38.pickle.dat
```