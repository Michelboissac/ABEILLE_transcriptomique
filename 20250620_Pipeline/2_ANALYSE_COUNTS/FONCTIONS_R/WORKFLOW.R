WORKFLOW = function(tableau_counts,normalisation = c("pas_de_normalisation","vst","vst_conditions", "vst_peu_de_genes","log_rpkm", "log_tpm", "rlog","binaire"),noms_genes="",liste_conditions_a_garder_etoile=NULL,nom_repertoire_output=""){
  #sources
  source(paste0(repertoire_fonctions,"/FONCTIONS.R"))

  #LIBRARY
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  
  #Creer le repertoire output
  if(nom_repertoire_output == ""){
    chemin_repertoire_output = paste0(OUTPUTS,"/",tableau_counts,"_",normalisation,"_",noms_genes)
    nom_repertoire_output=tableau_counts}
  else{chemin_repertoire_output = paste0(OUTPUTS,"/",nom_repertoire_output)}
  dir.create(chemin_repertoire_output)

  #PREPARATION TABLEAU COUNTS
  #recupere le tableau "counts"  (alignement vs 12000 genes)
  counts = function_recupere_tab_counts_txt(tableau_counts)                                  #Tableau alignements VS tout les genes
  #recupere longueurs des genes
  longueurs_genes = recupere_longueurs_genes(counts)                                         #faire avant "function_selection_gene_dans_counts_AND_CONVERSION_NOM_DES_GENES"
  #supprime les 1eres colonnes
  counts = function_supprimer_x_premieres_colonnes(counts,5)                                 #supprime les premiers colonnes du tableau (longueur,chromosome, etc ..)
  #change le nom des alignements
  counts = function_change_NOMS_ALIGNEMENTS(counts,tableau_counts)                           #change le nom des alignements (colonnes) par cexu du fichiers *names.txt
  #alignements que l'on veut garder
  counts = function_conditions_a_garder_dans_counts(counts,liste_conditions_a_garder_etoile) #pour recuperer uniquement certains alignements (colonnes)

  #NORMALISATION 
  #NORMALISATION de counts (a faire sur le plus grands jeux de donnée, tableau entiers avec 12000 genes)
  print("normalisation :")
  counts_norm = NORMALISATION(counts = counts,tableau_counts = tableau_counts,normalisation = normalisation,filename=paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,"_dispEsts.png"),liste_conditions_a_garder_etoile=liste_conditions_a_garder_etoile,longueurs_genes=longueurs_genes)

  #SELECTION GENES
  #selectionne une sous partie du tableau normalise avec les genes d'interets 
  print("genes :")
  counts_norm_genes = function_selection_gene_dans_counts_AND_CONVERSION_NOM_DES_GENES(noms_genes,counts_norm)

  #QUALITE PROFONDEUR DE SEQUENCAGE
  #creer des barplots non triés de l'expression du genome
  dir.create(paste0(chemin_repertoire_output,"/barplot_genome_expression_",normalisation))
  function_barplot_genome_expression(counts = counts_norm,filename = paste0(chemin_repertoire_output,"/barplot_genome_expression_",normalisation,"/",nom_repertoire_output,".",normalisation,"_barplot_genome_expression"))  
  #nbr de reads moyen/gene exprimés (>0reads) 
  function_barplot_nbr_reads_moyen_par_gene_exprime(counts = counts_norm,filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,"_all_genes.barplot_nbr_reads_moyen_par_gene_exprime"))
  #nbr de genes exprimés par alignements (>0 reads/gene)
  function_barplot_nbr_genes_expr(counts = counts,filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,"_all_genes.barplot_nbr_genes_expr"))
  
  function_boxplot_taille_reads(counts,
                                plot_titre = nom_repertoire_output,
                                filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".boxplot_nbr_reads.png")
  )
  function_boxplot_taille_reads(counts_norm,
                                plot_titre = nom_repertoire_output,
                                filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,".boxplot_nbr_reads.png")
  )
  function_boxplot_taille_reads(counts_norm_genes,
                                plot_titre = nom_repertoire_output,
                                filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,".",noms_genes,".boxplot_nbr_reads.png")
  )
  
  
  #QUALITE SIGNATURE TRANSCRIPTIONNELLE DISCRIMINANTE
  #dendrogramme bootstrap
  function_cluster_hierarchique_bootsrap(counts = counts_norm,filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,"_all_genes.dendrogramme_bootsrap"))
  #dendrogramme
  function_cluster_hierarchique(counts = counts_norm,filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,"_all_genes.dendrogramme"))
  #dendrogramme
  function_cluster_hierarchique(counts = counts_norm_genes,filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,".",noms_genes,".dendrogramme"))
  #ACP sur le tableau entiers normalise 
  function_ACP(counts = counts_norm,plot_titre = paste0(nom_repertoire_output,".",normalisation,"_all_genes"), filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,"_all_genes.acp.png"))
  #ACP sur le tableau genes normalise 
  function_ACP(counts = counts_norm_genes,plot_titre = paste0(nom_repertoire_output,".",normalisation,"_all_genes"), filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,".",noms_genes,".acp.png"))
  #ACP 3D
  function_3D_PCA(counts = counts_norm,filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,"_all_genes.acp3D.html"))
  
      #CORRELATIONS
  chemin_tab_noms_genes = paste0(repertoire_fichiers_txt,"/",noms_genes,".txt")                            #repertoire_fichiers_txt
  df_genes_names = read.table(chemin_tab_noms_genes, header=FALSE,  sep="\t", comment.char="#")
  R2=0.99
  function_CNN_genes(counts = counts_norm,df_genes_names = df_genes_names[[1]],R2 = R2,filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,".",noms_genes,".",R2))
  
  #
  function_matrice_correlation(counts = counts_norm_genes,plot_titre = paste0(nom_repertoire_output,".",normalisation,".",noms_genes),filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,".",noms_genes,".correlation.png"))
  #
  function_matrice_correlation_courbes(counts = counts_norm_genes,plot_titre = paste0(nom_repertoire_output,".",normalisation,".",noms_genes),filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,".",noms_genes,".correlation_courbes.png"))
  #graphe de correlation des genes d'interets
  function_CNN(counts = counts_norm_genes,R2= 0.5,filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,"_CNN"))
  function_CNN(counts = counts_norm_genes,R2= 0.6,filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,"_CNN"))
  function_CNN(counts = counts_norm_genes,R2= 0.7,filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,"_CNN"))
  function_CNN(counts = counts_norm_genes,R2= 0.8,filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,"_CNN"))
  function_CNN(counts = counts_norm_genes,R2= 0.9,filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,"_CNN"))
  print("coorelation :")
  function_CNN(counts = counts_norm_genes,R2= 0,filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,"_CNN"))
  #trie la matrice par ordre alphabetique (colonnes et lignes) pour la visualisation heatmap
  counts_norm_genes = counts_norm_genes[order(rownames(counts_norm_genes)), order(colnames(counts_norm_genes))] 
  
    #SAUVEGARDE TABLEAU COUNTS
  #sauvegarde le sous tableau normalise avec nos genes
  write.table(counts_norm_genes,paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,".",noms_genes,".counts.txt"),row.names=TRUE,col.names = TRUE, quote = FALSE, sep = "\t")
  #sauvegarde le tableau counts normalise
  write.table(counts_norm,paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,".","all_genes.counts.txt"),row.names=TRUE,col.names = TRUE, quote = FALSE, sep = "\t")
  
    #HEATMAP

  #heatmap sur le tableau entiers normalise
  function_heatmap(counts = counts_norm,plot_titre = paste0(nom_repertoire_output,".",normalisation,"_all_genes"),filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,"_all_genes.png"))
  #heatmap sur le sous tableau normalise avec nos genes
  function_heatmap(counts = counts_norm_genes,plot_titre = paste0(nom_repertoire_output,".",normalisation,".",noms_genes),filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,".",noms_genes,".heatmap.png"))
  #pheatmap sur le sous tableau normalise avec nos genes
  function_pheatmap(counts = counts_norm_genes,plot_titre = paste0(nom_repertoire_output,".",normalisation,".",noms_genes),filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,".",noms_genes,".pheatmap.png"))
}
