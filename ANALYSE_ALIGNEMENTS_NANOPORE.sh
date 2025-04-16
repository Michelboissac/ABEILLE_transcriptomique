DOSSIER_GENOME="/media/mboissac/T7 Shield/Michel_boissac/genome/"
DOSSIER_READS="/media/mboissac/T7 Shield/Michel_boissac/READS/"
nom_genome=GCF_003254395.2_Amel_HAv3.1_genomic
nom_genome_annotation=genomic.gtf
#########################################################################################
DOSSIER_TRAVAIL="$(pwd)/ANALYSE_NANOPORE/"
DOSSIER_SAM=${DOSSIER_TRAVAIL}SAM/
DOSSIER_BAM=${DOSSIER_TRAVAIL}BAM/
DOSSIER_SORTED_BAM_and_INDEX=${DOSSIER_TRAVAIL}SORTED_BAM_AND_INDEX/
DOSSIER_GENOME_TRAVAIL=${DOSSIER_TRAVAIL}GENOME/
DOSSIER_COUNTS=${DOSSIER_TRAVAIL}COUNTS/
genome=${DOSSIER_GENOME_TRAVAIL}GCF_003254395.2_Amel_HAv3.1_genomic.fna
genome_index=${DOSSIER_GENOME_TRAVAIL}GCF_003254395.2_Amel_HAv3.1_genomic.idx
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
cp "${DOSSIER_GENOME[@]}GCF_003254395.2_Amel_HAv3.1_genomic.fna" "${genome[@]}"
cp "${DOSSIER_GENOME[@]}genomic.gtf" "${genome_annotation_gtf[@]}"


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
#OUTILS

#minimap2 2.24-r1122
#samtools 1.13
#featureCounts 2.0.1



