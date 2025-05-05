

#POUR CREER UN ENVIRONNEMENT CONDA AVEC LES OUTILS :
#executer la commandes ci dessous dans le meme repertoire que environment.yml
conda env create -n mon_env_test -f environment.yml

#alternative sans utiliser environment.yml :

conda create --yes -n test_env -c bioconda -c conda-forge samtools=1.13 minimap2=2.22 hisat2=2.2.1 subread=2.0.1




