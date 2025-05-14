rm(list=ls())
setwd(getwd())
#################################################################
#FONCTIONS
repertoire_fonctions = paste0(getwd(),"/FONCTIONS_R")
source(paste0(repertoire_fonctions,"/WORKFLOW.R"))
#TXT
repertoire_fichiers_txt = paste0(getwd(),"/FICHIERS_TXT")
#OUTPUT PNG
OUTPUTS = paste0(getwd(),"/OUTPUT")

#################################################################
#WORKFLOW = function(tableau_counts,
#                    normalisation = c("pas_de_normalisation","vst","vst_conditions", "vst_peu_de_genes","log_rpkm", "log_tpm", "rlog"),
#                    noms_genes="",
#                    liste_conditions_a_garder_etoile=NULL,
#                    nom_repertoire_output="")
#############################################################


counts = WORKFLOW(tableau_counts = "2014_jasper_counts_illumina",
                  normalisation = "vst",
                  noms_genes = "OCTOPAMINE_TRP")



###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
#autres :
#TABLEAU POUR IGV : creer un tableau avec nom gene + coordonnee dans une df a partir de counts.
#function_creationtableau_pour_igv_screenshot_alignement_gene_start_end(counts,"start_end_genes_channels_nanopore_nachr.txt") #TABLEAU POUR IGV :
###############################################################
#counts1=function_recupere_tab_counts_txt("counts_paired_zhang_2022_bulk")
#counts2=function_recupere_tab_counts_txt("2014_jasper_counts_illumina")
#counts1 = function_supprimer_x_premieres_colonnes(counts1,5)                                 #supprime les premiers colonnes du tableau (longueur,chromosome, etc ..)
#counts2 = function_supprimer_x_premieres_colonnes(counts2,5)                                 #supprime les premiers colonnes du tableau (longueur,chromosome, etc ..)
#counts1 = function_change_NOMS_ALIGNEMENTS(counts1, "counts_paired_zhang_2022_bulk")                           #change le nom des alignements (colonnes) par cexu du fichiers *names.txt
#counts2 = function_change_NOMS_ALIGNEMENTS(counts2, "2014_jasper_counts_illumina")                           #change le nom des alignements (colonnes) par cexu du fichiers *names.txt
# Vérifie si les rownames sont identiques et dans le même ordre
#if (identical(rownames(counts1), rownames(counts2))) {
#  counts12=cbind(counts1,counts2)
#} else {
#  stop("Les noms de lignes ne correspondent pas ou ne sont pas dans le même ordre.")
#}
#write.table(counts12,"2014_jasper_counts_illumina_AND_zhang_2022_bulk.txt",row.names=TRUE,col.names = TRUE, quote = FALSE, sep = "\t")


