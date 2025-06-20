#1)installer conda
#https://www.anaconda.com/docs/getting-started/miniconda/install#linux
#mkdir -p ~/miniconda3
#wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
#bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
#rm ~/miniconda3/miniconda.sh
#source ~/miniconda3/bin/activate
#conda init --all

#2)creer un env conda avec les outils
#conda create --yes -n test_env2 -c bioconda -c conda-forge samtools=1.13 minimap2=2.22 hisat2=2.2.1 subread=2.0.1


#bash ANALYSE_ALIGNEMENTS_NANOPORE.sh "/media/mboissac/T7_Shield/Michel_boissac/DATAS/genome/" "/media/mboissac/T7_Shield/Michel_boissac/DATAS/20250429_skm_/" GCF_003254395.2_Amel_HAv3.1_genomic.fna genomic.gtf "/media/mboissac/T7_Shield/Michel_boissac/Pipeline/20250429_skm_/"


DOSSIER_GENOME=$1
DOSSIER_READS=$2
nom_genome=$3
nom_genome_annotation=$4
DOSSIER_TRAVAIL=$5


#########################################################################################
DOSSIER_SAM=${DOSSIER_TRAVAIL}SAM/
DOSSIER_BAM=${DOSSIER_TRAVAIL}BAM/
DOSSIER_SORTED_BAM_and_INDEX=${DOSSIER_TRAVAIL}SORTED_BAM_AND_INDEX/
DOSSIER_GENOME_TRAVAIL=${DOSSIER_TRAVAIL}GENOME/
DOSSIER_COUNTS=${DOSSIER_TRAVAIL}COUNTS/
genome=${DOSSIER_GENOME_TRAVAIL}${nom_genome}
genome_index=${DOSSIER_GENOME_TRAVAIL}${nom_genome}.idx
liste_de_reads_nanopore=$(ls "${DOSSIER_READS[@]}")
tableau_alignements=${DOSSIER_COUNTS}counts.txt
genome_annotation_gtf=${DOSSIER_GENOME_TRAVAIL}genomic.gtf
liste_alignements=()
mkdir "${DOSSIER_TRAVAIL[@]}"
mkdir "${DOSSIER_SAM[@]}"
mkdir "${DOSSIER_BAM[@]}"
mkdir "${DOSSIER_SORTED_BAM_and_INDEX[@]}"
mkdir "${DOSSIER_GENOME_TRAVAIL[@]}"
mkdir "${DOSSIER_COUNTS[@]}"
cp "${DOSSIER_GENOME[@]}${nom_genome}" "${genome[@]}"
cp "${DOSSIER_GENOME[@]}${nom_genome_annotation}" "${genome_annotation_gtf[@]}"


minimap2 -d "${genome_index[@]}" "${genome[@]}"

for nom_reads in $liste_de_reads_nanopore;
do
    reads=${DOSSIER_READS[@]}${nom_reads}
    alignement=${nom_genome}.VS.${nom_reads}
    
    minimap2 -ax splice "${genome[@]}" "${reads[@]}" > "${DOSSIER_SAM[@]}${alignement}.sam"

    samtools view -Sb "${DOSSIER_SAM[@]}${alignement}.sam" > "${DOSSIER_BAM[@]}${alignement}.bam"

    samtools sort "${DOSSIER_BAM[@]}${alignement}.bam" -o "${DOSSIER_SORTED_BAM_and_INDEX[@]}${alignement}.sorted.bam"

    samtools index "${DOSSIER_SORTED_BAM_and_INDEX[@]}${alignement}.sorted.bam"

    liste_alignements+=("${DOSSIER_SORTED_BAM_and_INDEX[@]}${alignement}.sorted.bam")
done

featureCounts -T 4 -t exon -g gene_id -a "${genome_annotation_gtf[@]}" -o "${tableau_alignements[@]}" "${liste_alignements[@]}"
#########################################################################################


