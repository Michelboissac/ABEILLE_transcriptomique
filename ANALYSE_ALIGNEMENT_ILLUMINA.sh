DOSSIER_GENOME="/media/mboissac/T7 Shield/Michel_boissac/genome/"
DOSSIER_READS="/media/mboissac/T7 Shield/Michel_boissac/SRA/"
nom_genome=GCF_003254395.2_Amel_HAv3.1_genomic
nom_genome_annotation=genomic.gtf
#########################################################################################
DOSSIER_TRAVAIL="$(pwd)/ANALYSE_ILLUMINA/"

DOSSIER_READS_TRAVAIL=${DOSSIER_TRAVAIL}SRA/
DOSSIER_SAM=${DOSSIER_TRAVAIL}SAM/
DOSSIER_BAM=${DOSSIER_TRAVAIL}BAM/
DOSSIER_SORTED_BAM_and_INDEX=${DOSSIER_TRAVAIL}SORTED_BAM_AND_INDEX/
DOSSIER_GENOME_TRAVAIL=${DOSSIER_TRAVAIL}GENOME/
DOSSIER_COUNTS=${DOSSIER_TRAVAIL}COUNTS/
genome=${DOSSIER_GENOME_TRAVAIL}GCF_003254395.2_Amel_HAv3.1_genomic.fna
genome_index=${DOSSIER_GENOME_TRAVAIL}GCF_003254395.2_Amel_HAv3.1_genomic.idx
tableau_alignements=${DOSSIER_COUNTS}counts.txt
genome_annotation_gtf=${DOSSIER_GENOME_TRAVAIL}genomic.gtf
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
#cp "${DOSSIER_GENOME[@]}GCF_003254395.2_Amel_HAv3.1_genomic.fna" "${genome[@]}"
#cp "${DOSSIER_GENOME[@]}genomic.gtf" "${genome_annotation_gtf[@]}"
#mv "${DOSSIER_READS[@]}" "${DOSSIER_TRAVAIL[@]}"


#hisat2_extract_splice_sites.py "${genome_annotation_gtf[@]}" > "${splicesites[@]}"
#hisat2_extract_exons.py "${genome_annotation_gtf[@]}" > "${exons[@]}"
#hisat2-build --ss "${splicesites[@]}" --exon "${exons[@]}" "${genome[@]}" "${genome_index[@]}"


readarray -t liste_de_reads_illumina < <(
  ls "${DOSSIER_READS_TRAVAIL[@]}"*.fastq|sed "s|${DOSSIER_READS_TRAVAIL[@]}||"|sed "s|[1-2]\.fastq||"|sort -u 
)
echo "${liste_de_reads_illumina[@]}"

for nom_reads in "${liste_de_reads_illumina[@]}"; do
    reads1="${nom_reads}1.fastq"
    reads2="${nom_reads}2.fastq"
    
    alignement=${nom_genome}.VS.${nom_reads}

    cd "${DOSSIER_TRAVAIL[@]}"

    hisat2 -x GENOME/GCF_003254395.2_Amel_HAv3.1_genomic.idx --dta --rna-strandness RF -1 SRA/"${reads1[@]}" -2 SRA/"${reads2[@]}" -S SAM/${alignement}.sam
    
    samtools view -Sb "${DOSSIER_SAM[@]}${alignement}.sam" > "${DOSSIER_BAM[@]}${alignement}.bam"
    
    samtools sort -n "${DOSSIER_BAM[@]}${alignement}.bam" -o "${DOSSIER_BAM[@]}${alignement}.name_sorted.bam"  # tri par nom pour fixmate
    
    samtools fixmate -m "${DOSSIER_BAM[@]}${alignement}.name_sorted.bam" "${DOSSIER_BAM[@]}${alignement}.name_sorted.fixmate.bam"
    
    samtools sort "${DOSSIER_BAM[@]}${alignement}.name_sorted.fixmate.bam" -o "${DOSSIER_SORTED_BAM_and_INDEX[@]}${alignement}.name_sorted.fixmate.sorted.bam"
    
    samtools index "${DOSSIER_SORTED_BAM_and_INDEX[@]}${alignement}.name_sorted.fixmate.sorted.bam"

    liste_alignements+=("${DOSSIER_SORTED_BAM_and_INDEX[@]}${alignement}.name_sorted.fixmate.sorted.bam")
    
done

#illumina



#featureCounts -T 4 -t exon -g gene_id -a "${genome_annotation_gtf[@]}" -o "${tableau_alignements[@]}" "${liste_alignements[@]}"



