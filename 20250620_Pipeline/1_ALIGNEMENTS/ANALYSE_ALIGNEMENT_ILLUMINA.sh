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

#bash ANALYSE_ALIGNEMENT_ILLUMINA.sh "/media/mboissac/T7_Shield/Michel_boissac/DATAS/genome/" "/media/mboissac/T7_Shield/Michel_boissac/Pipeline/test_analyse_illumina/" GCF_003254395.2_Amel_HAv3.1_genomic.fna genomic.gtf "/media/mboissac/T7_Shield/Michel_boissac/Pipeline/test_analyse_illumina_out/"

DOSSIER_GENOME=$1
DOSSIER_READS=$2
nom_genome=$3
nom_genome_annotation=$4
DOSSIER_TRAVAIL=$5

#########################################################################################


DOSSIER_READS_TRAVAIL=${DOSSIER_TRAVAIL}SRA/
DOSSIER_SAM=${DOSSIER_TRAVAIL}SAM/
DOSSIER_BAM=${DOSSIER_TRAVAIL}BAM/
DOSSIER_SORTED_BAM_and_INDEX=${DOSSIER_TRAVAIL}SORTED_BAM_AND_INDEX/
DOSSIER_GENOME_TRAVAIL=${DOSSIER_TRAVAIL}GENOME/
DOSSIER_COUNTS=${DOSSIER_TRAVAIL}COUNTS/
genome=${DOSSIER_GENOME_TRAVAIL}${nom_genome}
genome_index=${DOSSIER_GENOME_TRAVAIL}${nom_genome}.idx
tableau_alignements=${DOSSIER_COUNTS}counts.txt
genome_annotation_gtf=${DOSSIER_GENOME_TRAVAIL}${nom_genome_annotation}
liste_alignements=()
splicesites=${DOSSIER_GENOME_TRAVAIL}splicesites.txt
exons=${DOSSIER_GENOME_TRAVAIL}exons.txt
mkdir "${DOSSIER_TRAVAIL[@]}"
mkdir "${DOSSIER_SAM[@]}"
mkdir "${DOSSIER_BAM[@]}"
mkdir "${DOSSIER_SORTED_BAM_and_INDEX[@]}"
mkdir "${DOSSIER_GENOME_TRAVAIL[@]}"
mkdir "${DOSSIER_COUNTS[@]}"
mkdir "${DOSSIER_READS_TRAVAIL[@]}"
cp "${DOSSIER_GENOME[@]}${nom_genome}" "${genome[@]}"
cp "${DOSSIER_GENOME[@]}${nom_genome_annotation}" "${genome_annotation_gtf[@]}"
mv "${DOSSIER_READS[@]}/"* "${DOSSIER_READS_TRAVAIL[@]}/"

#"${DOSSIER_READS_TRAVAIL[@]}"


hisat2_extract_splice_sites.py "${genome_annotation_gtf[@]}" > "${splicesites[@]}"
hisat2_extract_exons.py "${genome_annotation_gtf[@]}" > "${exons[@]}"
hisat2-build --ss "${splicesites[@]}" --exon "${exons[@]}" "${genome[@]}" "${genome_index[@]}"


readarray -t liste_de_reads_illumina < <(
  ls "${DOSSIER_READS_TRAVAIL[@]}"*.fastq|sed "s|${DOSSIER_READS_TRAVAIL[@]}||"|sed "s|[1-2]\.fastq||"|sort -u 
)
echo "${liste_de_reads_illumina[@]}"

for nom_reads in "${liste_de_reads_illumina[@]}"; do
    reads1="${nom_reads}1.fastq"
    reads2="${nom_reads}2.fastq"
    
    alignement=${nom_genome}.VS.${nom_reads}

    cd "${DOSSIER_TRAVAIL[@]}"

    hisat2 -x GENOME/${nom_genome}.idx --dta --rna-strandness RF -1 SRA/"${reads1[@]}" -2 SRA/"${reads2[@]}" -S SAM/${alignement}.sam
    
    samtools view -Sb "${DOSSIER_SAM[@]}${alignement}.sam" > "${DOSSIER_BAM[@]}${alignement}.bam"
    
    samtools sort -n "${DOSSIER_BAM[@]}${alignement}.bam" -o "${DOSSIER_BAM[@]}${alignement}.name_sorted.bam"  # tri par nom pour fixmate
    
    samtools fixmate -m "${DOSSIER_BAM[@]}${alignement}.name_sorted.bam" "${DOSSIER_BAM[@]}${alignement}.name_sorted.fixmate.bam"
    
    samtools sort "${DOSSIER_BAM[@]}${alignement}.name_sorted.fixmate.bam" -o "${DOSSIER_SORTED_BAM_and_INDEX[@]}${alignement}.name_sorted.fixmate.sorted.bam"
    
    samtools index "${DOSSIER_SORTED_BAM_and_INDEX[@]}${alignement}.name_sorted.fixmate.sorted.bam"

    liste_alignements+=("${DOSSIER_SORTED_BAM_and_INDEX[@]}${alignement}.name_sorted.fixmate.sorted.bam")
    
done

featureCounts -T 4 -t exon -g gene_id -a "${genome_annotation_gtf[@]}" -o "${tableau_alignements[@]}" "${liste_alignements[@]}"



