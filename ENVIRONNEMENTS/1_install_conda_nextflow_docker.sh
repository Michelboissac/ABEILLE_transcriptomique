#1)installer conda
#https://www.anaconda.com/docs/getting-started/miniconda/install#linux
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh
source ~/miniconda3/bin/activate
conda init --all

#2)installer nextflow et docker dans un environnement conda : "nextflow_docker"
yes | conda create -n nextflow_docker
conda install bioconda::nextflow
sudo apt install docker-compose
