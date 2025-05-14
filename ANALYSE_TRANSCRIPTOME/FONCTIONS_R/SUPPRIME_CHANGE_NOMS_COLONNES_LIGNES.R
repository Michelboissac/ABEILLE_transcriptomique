function_recupere_tab_counts_txt = function(tableau_counts){
  tableau_counts_chemin=paste0(repertoire_fichiers_txt,"/",tableau_counts,".txt")
  counts <- read.table(tableau_counts_chemin, header=TRUE, row.names=1, sep="\t", comment.char="#")
}

function_change_NOMS_ALIGNEMENTS = function(counts,tableau_counts){
  noms_alignements=paste0(repertoire_fichiers_txt,"/",tableau_counts,"_names.txt")
  # Vérifier que le répertoire existe
  if (file.exists(noms_alignements)) {
    illumina_nom <- read.table(noms_alignements, header = FALSE, sep = "\t", comment.char = "#")
    for (echantillon in colnames(counts)) {
      if (echantillon %in% illumina_nom[,1]) {
        nom_a_remplacer <- illumina_nom[illumina_nom[,1] == echantillon, 2]
        print(paste("Remplacement :", echantillon, "->", nom_a_remplacer))
        colnames(counts)[colnames(counts) == echantillon] <- nom_a_remplacer
      } else {
        print(paste("Pas de correspondance pour :", echantillon))
      }
    }
    # action(s) à effectuer si le répertoire existe
  } else {
    message("Le répertoire n'existe pas : ", noms_alignements)
    # éventuellement arrêter ou proposer une alternative
  }


  return(counts)
}

function_selection_gene_dans_counts = function(liste_de_gene_a_selectionner,counts){
  true_false_channels_recep <- rownames(counts) %in% liste_de_gene_a_selectionner
  counts <- counts[true_false_channels_recep, ]
  return(counts)
}

function_CONVERSION_NOM_DES_GENES = function(df_nom_actuel_conversion,counts){
  noms_genes_actuel=df_nom_actuel_conversion[[1]]
  noms_genes_conversion=df_nom_actuel_conversion[[2]]
  mapping <- setNames(noms_genes_conversion, noms_genes_actuel)
  current_names <- rownames(counts)
  new_names <- ifelse(current_names %in% names(mapping),
                      mapping[current_names],
                      current_names)
  new_names <- make.unique(new_names)
  rownames(counts) <- new_names
  return(counts)
}

function_selection_gene_dans_counts_AND_CONVERSION_NOM_DES_GENES = function(noms_genes="",counts){
  if(noms_genes != ""){
    chemin_tab_noms_genes = paste0(repertoire_fichiers_txt,"/",noms_genes,".txt")                            #repertoire_fichiers_txt
    
    df_genes_names = read.table(chemin_tab_noms_genes, header=FALSE,  sep="\t", comment.char="#")
    counts = function_selection_gene_dans_counts(df_genes_names[[1]],counts)
    counts = function_CONVERSION_NOM_DES_GENES(df_genes_names,counts)
    
  }
  return(counts)
}


function_supprimer_x_premieres_colonnes = function(counts,x){
  for(i in 1:x){
    counts <- counts[, -c(1)]
  }
  return(counts)
}