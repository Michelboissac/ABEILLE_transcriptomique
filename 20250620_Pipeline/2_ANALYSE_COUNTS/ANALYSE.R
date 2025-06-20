#Pour utiliser cette fonction, il faut que ANALYSE.R soit dans le meme repertoire que le dossier FONCTIONS_R et le dossier FICHIERS_TXT, 
#ainsi que le dossier OUTPUT (à créer).

#FONCTIONS_R contient l'ensemble des fonctions necessaire (FONCTIONS.R) utilisée dans la fonction principale (WORKFLOW.R).

#FICHIERS_TXT contient le tableau counts 2014_jasper_counts_illumina.txt, ansi que counts_paired_zhang_2022_bulk.txt.
#Il contient egalement respectivment 2014_jasper_counts_illumina_names.txt et counts_paired_zhang_2022_bulk_names.txt,
#permettant de modifier les noms des analyses (colonnes de chaque tableau counts).
#il y a egalement le fichiers OCTOPAMINE_TRP_CAV_NAV.txt contenant une liste de genes que l'on veut etudier, et renommer.
#Genes.txt contient l'ensemble des 12000 genes avec leurs noms officiel dans le tableau counts.

#OUTPUT contient les dossiers de sortie du script.
#############################################################
rm(list=ls())
setwd(getwd())
gc()
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
                  normalisation = "log_tpm",
                  noms_genes = "OCTOPAMINE_TRP_CAV_NAV")
#############################################################




