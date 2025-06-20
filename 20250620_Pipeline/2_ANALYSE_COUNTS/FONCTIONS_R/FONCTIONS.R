function_pheatmap <- function(counts,plot_titre = "pheatmap", filename = "pheatmap.png", width = 1200, height = 1000, res = 150) {
  library(pheatmap)
  
  # Transposition et nettoyage
  counts <- t(counts)
  counts <- counts[rowSums(is.na(counts)) < ncol(counts), ]
  counts <- as.matrix(counts)
  counts[is.na(counts)] <- 0
  #counts <- counts[, colSums(counts) != 0]
   
  # Créer la matrice de labels
  number_matrix <- matrix(sprintf("%.1f", counts), 
                          nrow = nrow(counts), ncol = ncol(counts))
  
  # Mettre les "0.0" en texte blanc
  number_color <- matrix("black", nrow = nrow(counts), ncol = ncol(counts))
  number_color[counts == 0] <- "white"
  print(filename)
  png(filename, width = width, height = height, res = res)
  
  # Affichage du heatmap
  pheatmap(counts,
           main =plot_titre,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           display_numbers = number_matrix,
           number_color = number_color,
           fontsize_number = 10,
           col = colorRampPalette(c("white", "yellow", "red"))(100)
           )
  dev.off()
}

function_heatmap = function(counts,plot_titre = "heatmap", filename = "heatmap.png", width = 1200, height = 1000, res = 150){
  #HEATMPAP:
  counts <- t(counts)
  counts <- counts[rowSums(is.na(counts)) < ncol(counts), ]
  counts=as.matrix(counts)
  counts[is.na(counts)] <- 0
  #counts <- counts[, colSums(counts) != 0]
  counts <- counts[order(rownames(counts), decreasing = TRUE), order(colnames(counts))]
  
  
  png(filename, width = width, height = height, res = res)
  heatmap(counts,Rowv = NA, Colv = NA, margins = c(5, 5),scale = "none",keep.dendro = FALSE,main = plot_titre,
         col = colorRampPalette(c("white", "yellow", "red"))(100)
  )
  dev.off()
  
  
}

#sorte de diagramme de venn
function_genes_communs = function(counts){
  genes = rownames(counts)
  for(experience in colnames(counts)){
    for(gene in genes){
      nbr_de_reads = counts[gene,experience]
      if(nbr_de_reads>0){
        counts[gene,experience]=1 #gene
      }
    }
  }
  # Install if needed
  #install.packages("UpSetR")
  library(UpSetR)
  # Affichage du diagramme avec des options pour compacité
  upset(counts,
        sets = colnames(counts),
        keep.order = TRUE,
        sets.bar.color = "#56B4E9",
        order.by = "freq",   # trie les combinaisons par fréquence
        mb.ratio = c(0.6, 0.4),  # réduit la taille des barres du bas (main bar)
        text.scale = 1.2)     # ajuste la taille du texte pour lisibilité
}


function_ACP <- function(counts, plot_titre = "ACP", filename = "ACP_plot.png", width = 1200, height = 1000, res = 150){  #ACP :
  data_t <- t(counts)
  pca <- prcomp(data_t, scale. = F)
  summary(pca)
  
  png(filename, width = width, height = height, res = res)
  plot(pca$x[,1:2], col=1:nrow(data_t), pch=19, main = "plot_titre")
  text(pca$x[,1:2], labels=rownames(pca$x), pos=3)
  dev.off()
  
  
  
  
  
  
  
  
  
  # Voir les contributions aux composantes principales
  contributions <- pca$rotation
  
  # Quelle variable contribue le plus à PC1 ?
  #abs(contributions[,1])  # valeurs absolues pour l'importance
  # Classement décroissant
  liste_axe_gene_contribution = contributions[,1]
  liste_axe_gene_contribution = liste_axe_gene_contribution[liste_axe_gene_contribution>0.01]
  
  importance_PC1 <- sort(abs(liste_axe_gene_contribution), decreasing = TRUE)
  #print(importance_PC1)
  #boxplot(importance_PC1)
  write.csv(importance_PC1, paste0(filename,"contributions1.txt"))
  
  
  #importance_PC2 <- sort(abs(contributions[,2]), decreasing = TRUE)
  #print(importance_PC2)
  #write.csv(importance_PC2, paste0(filename,"contributions2.txt"))
  
  #importance_PC3 <- sort(abs(contributions[,3]), decreasing = TRUE)
  #print(importance_PC3)
  #write.csv(importance_PC3, paste0(filename,"contributions3.txt"))
  
  
  
}

function_acp_sur_illumina_enplusnano =function(){
  # Spécifie les noms à exclure
  colonnes_a_exclure <- c("b.6", "b.10", "b.14","h.6","h.10","h.14","s.6","s.10","s.14","p.6")  # adapte selon ton cas
  
  # 1. Séparer les données
  data_incluse <- counts[, !(colnames(counts) %in% colonnes_a_exclure)]
  data_projeter <- counts[, colnames(counts) %in% colonnes_a_exclure]
  
  # 2. Transposer pour l'ACP (échantillons = lignes)
  data_t <- t(data_incluse)
  pca <- prcomp(data_t, scale. = FALSE)
  
  # 3. Projeter les colonnes exclues
  # Transposer + s'assurer que les lignes (gènes) sont dans le même ordre
  proj_t <- t(data_projeter)
  proj_t <- proj_t[, match(names(pca$center), colnames(proj_t))]
  
  stopifnot(identical(rownames(data_projeter), rownames(data_incluse)))
  
  # Centrer avec le même centre que l'ACP
  proj_centered <- as.matrix(proj_t) - pca$center
  proj_coords <- proj_centered %*% pca$rotation  # projection dans l'espace PCA
  
  # 4. Plot
  x_all <- c(pca$x[,1], proj_coords[,1])
  y_all <- c(pca$x[,2], proj_coords[,2])
  png("pca_with_projection.png", width = 1200, height = 1000, res = 150)
  
  plot(pca$x[,1:2], col="blue", pch=19, main="PCA avec projections",
       xlim = range(x_all), ylim = range(y_all))
  text(pca$x[,1:2], labels=rownames(pca$x), pos=3, col="blue")
  
  points(proj_coords[,1:2], col="red", pch=17)
  text(proj_coords[,1:2], labels=rownames(proj_coords), pos=3, col="red")
  dev.off()
  
}


function_UMAP = function(counts){
  
  #install.packages("uwot")
  
  # Charger le package
  library(uwot)
  counts <- t(counts)
  
  
  # UMAP (par défaut en 2D)
  umap_result <- umap(counts)
  
  # Résultat = une matrice avec les coordonnées projetées
  head(umap_result)
  noms=colnames(counts)
  noms=as.factor(noms)
  
  # Optionnel : visualisation
  plot(umap_result, col = as.numeric(noms), pch = 19,
       main = "Projection UMAP de iris")
  legend("topright", legend = levels(noms),
         col = 1:3, pch = 19) 
  
}


function_3D_PCA = function(counts, filename = "ACP_3D.html"){
  #PCA 3d : 
  library(plotly)
  counts <- t(counts)  
  pca <- prcomp(counts, scale. = FALSE)
  pca_df <- as.data.frame(pca$x)
  
  p <- plot_ly(data = pca_df, 
               
               x = ~PC1, y = ~PC2, z = ~PC3, 
               type = 'scatter3d', 
               mode = 'markers+text',
               text = rownames(pca_df),
               textposition = 'top center',
               marker = list(size = 5,
                             color = as.numeric(as.factor(rownames(pca_df))),
                             colorscale = 'Viridis'))
  
  htmlwidgets::saveWidget(p, filename)
  #browseURL(filename)
}

function_cluster_hierarchique = function(counts,filename){
  counts <- t(counts)  # maintenant [40 samples x 12000 gènes]
  d <- dist(counts, method = "euclidean")
  hc <- hclust(d, method = "ward.D2")  # méthode ward.D2 souvent bonne pour clustering transcriptomique
  
  png(paste0(filename,".png"), width = 800, height = 600)  # tu peux ajuster la taille
  plot(hc, main = "Hierarchical Clustering des échantillons", xlab = "", sub = "", cex = 0.9)
  dev.off() 
}


function_cluster_hierarchique_bootsrap = function(counts, filename){
  #DENDROGRAMME AVEC BOOTSTRAP POUR VERIFIER QUALTIE DATA ?
  library("pvclust")
  result <- pvclust(counts, method.hclust = "ward.D2", method.dist = "euclidean", nboot = 100)
  png(paste0(filename,".png"), width = 800, height = 600)
  
  plot(result)
  pvrect(result, alpha=0.95)  # entoure les clusters avec p-value > 95%
  
  dev.off() 
}

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
  counts <- counts[true_false_channels_recep, , drop = FALSE]
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

function_supprimer_x_premieres_colonnes <- function(counts, x) {
  for (i in 1:x) {
    counts <- counts[, -1, drop = FALSE]
  }
  return(counts)
}

function_collage_tableau = function(counts1,counts2,counts12_name){
  # Vérifie si les rownames sont identiques et dans le même ordre
  if (identical(rownames(counts1), rownames(counts2))) {
    counts12=cbind(counts1,counts2)
  } else {
    stop("Les noms de lignes ne correspondent pas ou ne sont pas dans le même ordre.")
  }
  write.table(counts12,counts12_name,row.names=TRUE,col.names = TRUE, quote = FALSE, sep = "\t")
}


function_conditions_a_garder_dans_counts = function(counts,liste_conditions_a_garder_etoile){
  if(!is.null(liste_conditions_a_garder_etoile)){
    liste_conditions_a_garder_etoile=paste0("^",liste_conditions_a_garder_etoile)     #enlever le ^si on veut grep n'importe ou et pas que au debut
    liste_conditions_a_garder_entiere=c()
    for(condition in liste_conditions_a_garder_etoile){
      liste_conditions_a_garder_entiere=c(liste_conditions_a_garder_entiere,grep(condition, colnames(counts)))
    }  
    counts <- counts[, liste_conditions_a_garder_entiere]
  }
  return(counts)
}

function_ANALYSE_EXPRESSION_DIFFERENTIELLE = function(counts,listes_conditions_a_garder){
  counts = function_conditions_a_garder_dans_counts(counts,listes_conditions_a_garder)
  colonnes_conditions=c()
  
  for(condition in listes_conditions_a_garder){
    colonne = grep(paste0("^",condition), names(counts))
    for(i in colonne){colonnes_conditions=c(colonnes_conditions,condition)}
  }
  
  metadata <- data.frame(
    row.names = colnames(counts),
    condition = colonnes_conditions  # Adapte selon ton design expérimental
  )
  metadata$condition <- factor(metadata$condition, levels = listes_conditions_a_garder)
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ condition)
  
  
  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- DESeq(dds)    #ANALYSE EXPRESSION DIFF
  res <- results(dds)
  res <- res[order(res$padj), ] 
  dev.new()
  plotMA(res, ylim=c(-2,2))
  
  
  
  condition_VS=resultsNames(dds)[2]
  resLFC <- lfcShrink(dds, coef=condition_VS, type="apeglm") #Log fold change shrinkage for visualization and ranking
  resLFC
  resOrdered <- res[order(res$pvalue),] #p-values and adjusted p-values
  summary(resLFC)
  
  #sum(res$padj < 0.1, na.rm=TRUE) #How many adjusted p-values were less than 0.1?
  
  #results function automatically performs independent filtering based on the mean of normalized counts for each gene, optimizing the number of genes which will have an adjusted p value below a given FDR cutoff, alpha. Independent filtering is further discussed below. By default the argument alpha is set to 0.1. If the adjusted p value cutoff will be a value other than 0.1, alpha should be set to that value:
  #res05 <- results(dds, alpha=0.05)  
  #summary(res05)    
  #sum(res05$padj < 0.05, na.rm=TRUE)
  
  
  #Independent hypothesis weighting
  #A generalization of the idea of p value filtering is to weight hypotheses to optimize power
  library("IHW")
  resIHW <- results(dds, filterFun=ihw)
  summary(resIHW)
  sum(resIHW$padj < 0.1, na.rm=TRUE)
  metadata(resIHW)$ihwResult
  
  #MA-plot
  #In DESeq2, the function plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet. Points will be colored blue if the adjusted p value is less than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.
  dev.new()
  plotMA(resLFC, ylim=c(-2,2))
  
  #It is more useful to visualize the MA-plot for the shrunken log2 fold changes, which remove the noise associated with log2 fold changes from low count genes without requiring arbitrary filtering thresholds.
  dev.new()
  plotMA(resLFC, ylim=c(-2,2))
  
  #After calling plotMA, one can use the function identify to interactively detect the row number of individual genes by clicking on the plot. One can then recover the gene identifiers by saving the resulting indices:
  #idx <- identify(res$baseMean, res$log2FoldChange)
  #rownames(res)[idx]
  
  #Alternative shrinkage estimators
  resNorm <- lfcShrink(dds, coef=2, type="normal")
  resAsh <- lfcShrink(dds, coef=2, type="ashr")
  par(mfrow=c(1,3), mar=c(4,4,2,1))
  xlim <- c(1,1e5); ylim <- c(-3,3)
  plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
  plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
  plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
  
  #plot counts
  dev.new()
  plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
  
  #more information on resultas columns
  #mcols(res)$description
  
  #Exporting results to CSV files
  #write.csv(as.data.frame(resOrdered), file="condition_treated_results.csv")
  
  #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
  #resSig <- subset(resOrdered, padj < 0.1)
  #resSig
  
  #Multi-factor designs : pour effet de lot etc .. ou diff entre pairend end single end etc ..
  
  #Effects of transformations on the variance :
  # this gives log2(n + 1)
  ntd <- normTransform(dds)
  library("vsn")
  dev.new()
  meanSdPlot(assay(ntd))
  #meanSdPlot(assay(vsd))
  #meanSdPlot(assay(rld))
  
  
  #Heatmap of the count matrix
  library("pheatmap")
  
  select <- order(rowMeans(counts(dds,normalized=TRUE)),
                  decreasing=TRUE)[1:20]
  print(colData(dds))
  df <- as.data.frame(colData(dds)[,c("condition","sizeFactor")])  #sizeFactor = type ?
  
  dev.new()
  pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
           cluster_cols=FALSE, annotation_col=df)
  
  
  #
  #pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
  #         cluster_cols=FALSE, annotation_col=df)
  
  #
  #pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
  #         cluster_cols=FALSE, annotation_col=df)
  
  #Heatmap of the sample-to-sample distances
  #sampleDists <- dist(t(assay(vsd)))
  sampleDists <- dist(t(assay(ntd)))
  library("RColorBrewer")
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(ntd$condition, ntd$type, sep="-")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  
  dev.new()
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
  
  #Principal component plot of the samples
  #plotPCA(vsd, intgroup=c("condition", "type"))
  dev.new()
  plotPCA(ntd, intgroup=c("condition", "sizeFactor")) #sizeFactor = type ?
  
  #dispersion
  dev.new()
  plotDispEsts(dds)
  
  #independent filtering of results
  dev.new()
  metadata(res)$alpha
  metadata(res)$filterThreshold
  plot(metadata(res)$filterNumRej, 
       type="b", ylab="number of rejections",
       xlab="quantiles of filter")
  lines(metadata(res)$lo.fit, col="red")
  abline(v=metadata(res)$filterTheta)
}

#######################################################################
recupere_longueurs_genes = function(counts){
  longueurs_genes=counts$Length #recupere les longueurs des genes pour la noramlisation
  return(longueurs_genes)
}
function_normalise_longueur_gene = function(counts,longueurs_genes){
  
  counts=counts/longueurs_genes
  return(counts)
}
function_normalise_log_plus_1 = function(counts){
  counts=log(counts+1)
  return(counts)
}
function_normalise_Reads_Per_Million = function(counts){
  for (col in colnames(counts)) {
    sample = counts[[col]]
    nbr_reads_total_samples = sum(sample)
    facteur_mise_a_echelle_per_million= nbr_reads_total_samples/1000000
    RPM = sample/facteur_mise_a_echelle_per_million #reads per million
    counts[[col]]=RPM
    
  }  
  return(counts)
}
#######################################################################
function_normalise_log_TPM = function(counts,tableau_counts,longueurs_genes){
  counts = function_normalise_longueur_gene(counts,longueurs_genes)
  counts = function_normalise_Reads_Per_Million(counts)
  counts = function_normalise_log_plus_1(counts)
  return(counts)
}

function_normalise_log_RPKM = function(counts,tableau_counts,longueurs_genes){
  counts = function_normalise_Reads_Per_Million(counts)
  counts = function_normalise_longueur_gene(counts,longueurs_genes)
  counts = function_normalise_log_plus_1(counts)
  return(counts)
}

function_normalisation_VST = function(counts,tableau_counts,filename = "plotDispEsts.png", width = 1200, height = 1000, res = 150){
  #Utilise une version plus rapide, mais nécessite suffisamment de gènes exprimés (sinon erreur)
  colnames(counts)=make.unique(colnames(counts),sep = ".") 
  metadata <- data.frame(
    row.names = colnames(counts),
    condition = colnames(counts) 
  )
  
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ 1)
  vst_data_condition <- vst(dds, blind = TRUE) #laisser sur TRUE, calcul la VST en ne prenant pas compte des conditions. sinon normalise difference ?
  vst_matrix_condition <- assay(vst_data_condition)
  #dispersion
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  
  png(filename, width = width, height = height, res = res)
  plotDispEsts(dds)
  dev.off()
  
  return(vst_matrix_condition)
}

function_normalisation_VST_conditions = function(counts,tableau_counts,filename = "plotDispEsts.png", width = 1200, height = 1000, res = 150){
  #Utilise une version plus rapide que function_normalisation_VST_peu_de_genes , mais nécessite suffisamment de gènes exprimés (sinon erreur)
  print(colnames(counts))
  metadata <- data.frame(
    row.names = colnames(counts),
    condition = sapply(strsplit(colnames(counts), split = "\\."), `[`, 1)    # passe de "dt_fo.trucaenelever"  à "dt_fo"
  )
  print(metadata)
  
  colnames(counts)= make.unique(colnames(counts),sep = ".")
  print(colnames(counts))
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~condition)
  vst_data_condition <- vst(dds, blind = FALSE) #laisser sur TRUE, calcul la VST en ne prenant pas compte des conditions. sinon normalise difference ?
  vst_matrix_condition <- assay(vst_data_condition)
  #dispersion
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  
  png(filename, width = width, height = height, res = res)
  plotDispEsts(dds)
  dev.off()
  
  return(vst_matrix_condition)
}

function_normalisation_VST_peu_de_genes<- function(counts,tableau_counts,filename = "plotDispEsts.png", width = 1200, height = 1000, res = 150){
  
  #Plus robuste, fonctionne même avec peu de gènes exprimés ou des données peu denses
  colnames(counts)=make.unique(colnames(counts),sep = ".") 
  metadata <- data.frame(
    row.names = colnames(counts),
    condition = colnames(counts) 
  )
  
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ 1)
  
  # Utilisation de varianceStabilizingTransformation à la place de vst()
  vst_data_condition <- varianceStabilizingTransformation(dds, blind = TRUE)
  vst_matrix_condition <- assay(vst_data_condition)
  
  # Affichage de la dispersion
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  
  png(filename, width = width, height = height, res = res)
  plotDispEsts(dds)
  dev.off()
  
  return(vst_matrix_condition)
}

function_normalisation_rlog = function(counts,tableau_counts,filename = "plotDispEsts.png", width = 1200, height = 1000, res = 150){
  
  colnames(counts)=make.unique(colnames(counts),sep = ".") 
  metadata <- data.frame(
    row.names = colnames(counts),
    condition =  colnames(counts) 
  )
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata,design = ~ 1)
  rlog_data <- rlog(dds, blind = TRUE)
  rlog_matrix <- assay(rlog_data)
  #dispersion
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  plotDispEsts(dds)
  
  return(rlog_matrix)
}

function_conversion_expression_binaire = function(counts){
  counts[counts > 0] = 1
  return(counts)
}

#######################################################################
NORMALISATION = function(counts, tableau_counts, normalisation = c("pas_de_normalisation","vst","vst_conditions", "vst_peu_de_genes","log_rpkm", "log_tpm", "rlog","binaire"),filename,liste_conditions_a_garder_etoile=NULL,longueurs_genes){
  normalisation <- match.arg(normalisation)
  if(normalisation=="pas_de_normalisation"){print("pas_de_normalisation")}
  if(normalisation=="log_rpkm"){counts = function_normalise_log_RPKM(counts,tableau_counts,longueurs_genes)}
  if(normalisation=="log_tpm"){counts = function_normalise_log_TPM(counts,tableau_counts,longueurs_genes)}
  if(normalisation=="rlog"){counts = function_normalisation_rlog(counts,tableau_counts)}
  if(normalisation=="vst"){counts = function_normalisation_VST(counts,tableau_counts,filename)}
  if(normalisation=="vst_peu_de_genes"){counts = function_normalisation_VST_peu_de_genes(counts,tableau_counts,filename)}
  if(normalisation=="vst_conditions"){counts = function_normalisation_VST_conditions(counts,tableau_counts,filename)}
  if(normalisation=="binaire"){counts = function_conversion_expression_binaire(counts)}
  
  return(counts)
}
#######################################################################

function_creationtableau_pour_igv_screenshot_alignement_gene_start_end = function(counts,output_txt){
  chromosomes = counts$Chr
  debuts = counts$Start
  fins = counts$End
  liste_de_genes_start_stop=c()
  for(gene in 1:length(chromosomes)){
    nom_gene=rownames(counts[gene,])
    isoformes=counts[gene,1]
    isoformes_start=counts[gene,2]
    isoformes_end=counts[gene,3]
    isoformes=strsplit(isoformes,";"  )[[1]]
    isoformes_start=strsplit(isoformes_start,";"  )[[1]]
    isoformes_end=strsplit(isoformes_end,";"  )[[1]]
    isoformes_start = as.numeric(isoformes_start)
    isoformes_end = as.numeric(isoformes_end)
    if(length(unique(isoformes)) == 1){
      iso=unique(isoformes)
      iso_start=min(isoformes_start)
      iso_end=max(isoformes_end)
      gene_start_end = paste0(iso,":",iso_start,"-",iso_end)
      names(gene_start_end)=nom_gene
      liste_de_genes_start_stop=c(liste_de_genes_start_stop,gene_start_end)
    }
    
  }
  
  write.table(liste_de_genes_start_stop, output_txt, sep = "\t", quote = FALSE, row.names = TRUE)
  
}

function_barplot_nbr_genes_expr = function(counts, filename){
  noms_colonnes = colnames(counts)
  noms_genes = rownames(counts)
  liste_nbr_genes_exprimes=c()
  for(colonne in noms_colonnes){
    col = counts[,colonne]
    #
    nbr_genes_exprimes = sum(col != 0)  #nbr genes avec + de 1 reads
    names(nbr_genes_exprimes)=colonne
    liste_nbr_genes_exprimes = c(liste_nbr_genes_exprimes,nbr_genes_exprimes)
  }
  png(paste0(filename,".png"), width = 800, height = 600) 
  barplot(liste_nbr_genes_exprimes,main ="Nombre de genes exprimés",las = 2 )
  dev.off()  # ferme le fichier
}

function_barplot_nbr_reads_moyen_par_gene_exprime = function(counts, filename){
  noms_colonnes = colnames(counts)
  noms_genes = rownames(counts)
  liste_nbr_reads_moyen_par_gene_exprime=c()
  liste_nbr_genes_exprimes=c()
  
  for(colonne in noms_colonnes){
    col = counts[,colonne]
    nbr_reads_total = sum(col)
    nbr_genes = length(col)
    nbr_reads_moyen_par_gene_all_genome = nbr_reads_total/nbr_genes
    
    #
    nbr_genes_exprimes = sum(col != 0)  #nbr genes avec + de 1 reads
    names(nbr_genes_exprimes)=colonne
    liste_nbr_genes_exprimes = c(liste_nbr_genes_exprimes,nbr_genes_exprimes)
    
    #
    nbr_reads_moyen_par_gene_exprime = nbr_reads_total/nbr_genes_exprimes
    names(nbr_reads_moyen_par_gene_exprime)=colonne
    liste_nbr_reads_moyen_par_gene_exprime = c(liste_nbr_reads_moyen_par_gene_exprime,nbr_reads_moyen_par_gene_exprime)
    
  }
  png(paste0(filename,".png"), width = 800, height = 600) 
  barplot(liste_nbr_reads_moyen_par_gene_exprime,main = "Nombre de reads moyens /genes exprimés",las = 2)
  dev.off()  # ferme le fichier
}

function_barplot_genome_expression = function(counts, filename){
  noms_colonnes = colnames(counts)
  for(colonne in noms_colonnes){
    col = counts[,colonne]
    #
    png(paste0(filename,colonne,".png"), width = 800, height = 600)  # tu peux ajuster la taille
    barplot(col,main = colonne,las = 2)   #permet de garder l'ordre des genes et de voir des profils
    dev.off()  # ferme le fichier
  }
}

function_garde_500_genes__les_plus_variables = function(counts){
  gene_sd <- apply(counts, 1, sd)
  top_genes <- names(sort(gene_sd, decreasing = TRUE)[1:500])
  counts = counts[top_genes,]
  
  return(counts)
}

function_CNN_genes = function(counts,df_genes_names=df_genes_names,R2 = 0.99,filename){
  library("igraph")
  library("corrr")
  # Supprimer les gènes constants avant la corrélation
  counts_filtered <- counts[apply(counts, 1, function(x) sd(x) != 0), ]
  
  # Puis recalculer la corrélation
  cor_matrix <- correlate(t(counts_filtered), method = "pearson")
  
  # 2. Transformer la matrice en format utilisable
  cor_df <- stretch(cor_matrix) # passer en format "long"
  
  # 3. Filtrer les fortes corrélations
  cor_df <- subset(cor_df, abs(r) >= R2)  
  
  # 4. Construire le graph
  g <- graph_from_data_frame(cor_df, directed = FALSE)
  
  
  # 5. Dessiner le graph
  # Définir les couleurs selon le signe de la corrélation
  edge_colors <- ifelse(E(g)$r < 0, "red", "green")
  
  # Générer le PNG
  png(paste0(filename,".",R2,".png"), width = 1200, height = 1000, res = 150)
  
  plot(g, 
       vertex.label = V(g)$name,
       vertex.label.cex = 0.7,
       vertex.size = 5,
       edge.width = abs(E(g)$r) * 5,  # épaisseur proportionnelle à la force
       edge.color = edge_colors,      # couleur selon le signe
       main = "CNN",
       layout = layout_with_fr)
  
  dev.off()
  
  
  
  
  
  # 1. Filtrer les gènes qui existent dans le graphe
  genes_cibles_valides <- intersect(df_genes_names, V(g)$name)
  
  # 2. Vérifier qu'au moins un gène est présent
  if (length(genes_cibles_valides) > 0) {
    # 3. Trouver les voisins directs (ordre 1) + inclure les gènes eux-mêmes
    voisins <- unlist(neighborhood(g, order = 1, nodes = genes_cibles_valides, mode = "all"))
    
    # 4. Extraire les noms des sommets à inclure dans le sous-graphe
    sommets_sousgraphe <- unique(V(g)[voisins]$name)
    
    # 5. Créer le sous-graphe
    g_sous <- induced_subgraph(g, vids = sommets_sousgraphe)
    
    # 6. Afficher ou enregistrer le sous-graphe
    png(filename = paste0(filename,".png"), width = 8000, height = 8000, res = 400)
    # Créer un vecteur de couleurs : bleu pour les gènes cibles, orange pour les autres
    couleurs_sommets <- ifelse(V(g_sous)$name %in% genes_cibles_valides, "blue", "orange")
    layout_graphopt <- layout_with_graphopt(g_sous, charge = 0.30, niter = 2000)
    plot(g_sous,
         vertex.label = V(g_sous)$name,
         vertex.label.cex = 0.8,
         vertex.size = 6,
         vertex.color = couleurs_sommets,  # <<<< couleurs ici
         edge.width = abs(E(g_sous)$r) * 5,
         main = "Gènes cibles (bleu) + voisins (orange)",
         layout = layout_graphopt)
    
    dev.off()
  } else {
    cat("Aucun des gènes cibles n'est présent dans le graphe.\n")
  }
  
  
  
  
  
}

function_matrice_correlation = function(counts,plot_titre = "pheatmap", filename = "pheatmap.png", width = 1200, height = 1000, res = 150){
  library("igraph")
  library("corrr")
  
  cor_matrix <- correlate(t(counts), method = "pearson")
  cor_matrix = as.data.frame(cor_matrix,row.names = TRUE)
  rownames(cor_matrix) = cor_matrix$term
  cor_matrix = cor_matrix[,-1]
  cor_matrix[is.na(cor_matrix)] <- 0
  cor_matrix =as.matrix(cor_matrix)
  
  
  png(filename, width = width, height = height, res = res)
  
  # Créer la matrice de labels
  number_matrix <- matrix(sprintf("%.1f", cor_matrix), 
                          nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))
  
  # Mettre les "0.0" en texte blanc
  number_color <- matrix("black", nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))
  number_color[cor_matrix == 0] <- "white"
  # Affichage du heatmap
  pheatmap(cor_matrix,
           main =plot_titre,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           display_numbers = number_matrix,
           number_color = number_color,
           fontsize_number = 10,
           col = colorRampPalette(c("white", "yellow", "red"))(100)
  )
  dev.off()  
}

function_matrice_correlation_courbes = function(counts,plot_titre = "pheatmap", filename = "pheatmap.png", width = 1200, height = 1000, res = 150){
  library("igraph")
  library("corrr")
  
  #install.packages("GGally")     # à faire une seule fois
  library(GGally)
  library(ggplot2)
  # Affiche les nuages de points pour toutes les paires de variables
  
  
  p = ggpairs(
    t(counts),
    
    upper = list(continuous = wrap("cor", size = 3)),
    
    lower = list(
      continuous = wrap("smooth", 
                        method = "lm", 
                        se = FALSE,
                        color = "red", 
                        fullrange = TRUE,
                        alpha = 0.8)
    ),
    
    diag = list(continuous = wrap("densityDiag", alpha = 0.5)),
    title = plot_titre
  ) +
    theme_minimal()
  ggsave(filename = filename, plot = p,
         width = 12, height = 12, dpi = 300, units = "in")
  
}

function_CNN_seurat = function(counts,R2 = 0.99,filename){
  
  library("igraph")
  library("corrr")
  library("dplyr")
  counts_filtered <- counts[apply(counts, 1, function(x) sd(x) != 0), ]
  cor_matrix <- correlate(t(counts_filtered), method = "pearson")
  cor_df <- stretch(cor_matrix) # passer en format "long"
  cor_df <- cor_df %>% filter(cor_df$x != cor_df$y)
  
  # 5. Réordonner les paires pour ne garder qu’une direction (évite doublons)
  cor_df <- cor_df %>%
    mutate(pair_id = paste0(pmin(cor_df$x, cor_df$y), "_", pmax(cor_df$x, cor_df$y))) %>%
    distinct(pair_id, .keep_all = TRUE)
  
  #cor_df <- subset(cor_df, r*r >= R2)         #R2 varieance
  cor_df <- subset(cor_df, abs(r) >= R2)     #coef de correlation
  
  g <- graph_from_data_frame(cor_df, directed = FALSE)
  edge_colors <- ifelse(E(g)$r < 0, "red", "green")
  png(paste0(filename,".",R2,".png"), width = 1200, height = 1000, res = 150)
  plot(g, 
       vertex.label = V(g)$name,
       vertex.label.cex = 0.7,
       vertex.size = 5,
       edge.width = abs(E(g)$r) * 5,  # épaisseur proportionnelle à la force
       edge.color = edge_colors,      # couleur selon le signe
       main = "CNN",
       layout = layout_with_fr)
  dev.off()
  
}

function_CNN = function(counts,R2 = 0.99,filename){
  library("igraph")
  library("corrr")
  library("dplyr")
  # Supprimer les gènes constants avant la corrélation
  counts_filtered <- counts[apply(counts, 1, function(x) sd(x) != 0), ]
  
  # Puis recalculer la corrélation
  cor_matrix <- correlate(t(counts_filtered), method = "pearson")
  
  # 2. Transformer la matrice en format utilisable
  cor_df <- stretch(cor_matrix) # passer en format "long"
  
  # 4. Supprimer les auto-corrélations (x == y)
  cor_df <- cor_df %>% filter(cor_df$x != cor_df$y)
  
  # 5. Réordonner les paires pour ne garder qu’une direction (évite doublons)
  cor_df <- cor_df %>%
    mutate(pair_id = paste0(pmin(cor_df$x, cor_df$y), "_", pmax(cor_df$x, cor_df$y))) %>%
    distinct(pair_id, .keep_all = TRUE)
  
  # 3. Filtrer les fortes corrélations
  cor_df <- subset(cor_df, abs(r) >= R2)  
  
  # 4. Construire le graph
  g <- graph_from_data_frame(cor_df, directed = FALSE)
  
  # 5. Dessiner le graph
  # Définir les couleurs selon le signe de la corrélation
  edge_colors <- ifelse(E(g)$r < 0, "red", "green")
  
  # Générer le PNG
  png(paste0(filename,".",R2,".png"), width = 1200, height = 1000, res = 150)
  
  plot(g, 
       vertex.label = V(g)$name,
       vertex.label.cex = 0.7,
       vertex.size = 5,
       edge.width = abs(E(g)$r) * 5,  # épaisseur proportionnelle à la force
       edge.color = edge_colors,      # couleur selon le signe
       main = "CNN",
       layout = layout_with_fr)
  
  dev.off()
  
}

function_CNN_seurat_multiple = function(counts,liste_conditions_a_garder_etoile,chemin_repertoire_output,nom_repertoire_output,normalisation,noms_genes){
  
  chemin_repertoire_output = paste0(chemin_repertoire_output,"/",liste_conditions_a_garder_etoile[1])
  dir.create(chemin_repertoire_output)
  
  function_CNN_seurat(counts = counts,R2= 0,filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,"_",noms_genes,"_",liste_conditions_a_garder_etoile[1],"_CNN"))
  function_CNN_seurat(counts = counts,R2= 0.001,filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,"_",noms_genes,"_",liste_conditions_a_garder_etoile[1],"_CNN"))
  function_CNN_seurat(counts = counts,R2= 0.1,filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,"_",noms_genes,"_",liste_conditions_a_garder_etoile[1],"_CNN"))
  function_CNN_seurat(counts = counts,R2= 0.2,filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,"_",noms_genes,"_",liste_conditions_a_garder_etoile[1],"_CNN"))
  function_CNN_seurat(counts = counts,R2= 0.3,filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,"_",noms_genes,"_",liste_conditions_a_garder_etoile[1],"_CNN"))
  function_CNN_seurat(counts = counts,R2= 0.4,filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,"_",noms_genes,"_",liste_conditions_a_garder_etoile[1],"_CNN"))
  function_CNN_seurat(counts = counts,R2= 0.5,filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,"_",noms_genes,"_",liste_conditions_a_garder_etoile[1],"_CNN"))
  function_CNN_seurat(counts = counts,R2= 0.6,filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,"_",noms_genes,"_",liste_conditions_a_garder_etoile[1],"_CNN"))
  function_CNN_seurat(counts = counts,R2= 0.7,filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,"_",noms_genes,"_",liste_conditions_a_garder_etoile[1],"_CNN"))
  function_CNN_seurat(counts = counts,R2= 0.8,filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,"_",noms_genes,"_",liste_conditions_a_garder_etoile[1],"_CNN"))
  function_CNN_seurat(counts = counts,R2= 0.9,filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,"_",noms_genes,"_",liste_conditions_a_garder_etoile[1],"_CNN"))
  function_CNN_seurat(counts = counts,R2= 0.99,filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,"_",noms_genes,"_",liste_conditions_a_garder_etoile[1],"_CNN"))
  
}

function_creation_graph_coexpression_cellule_expression_binaire_condition = function(counts,condition,noms_genes){
  counts_condition = function_selection_gene_dans_counts_AND_CONVERSION_NOM_DES_GENES(noms_genes,counts)
  liste_genes = rownames(counts_condition)
  df <- data.frame(matrix(0, nrow = length(liste_genes), ncol = length(liste_genes)))
  rownames(df) = liste_genes
  colnames(df) = liste_genes
  counts_condition = function_conditions_a_garder_dans_counts(counts_condition,c(condition)) #pour recuperer uniquement certains alignements (colonnes)
  counts_condition = function_conversion_expression_binaire(counts_condition)
  counts_condition = counts_condition[order(rownames(counts_condition)), order(colnames(counts_condition))] 
  liste_nom_cellules = colnames(counts_condition)
  for(nom_cellule in liste_nom_cellules){
    genes = rownames(counts_condition)
    cellule = counts_condition[[nom_cellule]]
    names(cellule) = genes
    cellule = cellule[cellule != 0]     
    genes_coexprime_dans_cellule = names(cellule)
    for(gene1 in genes_coexprime_dans_cellule){
      for(gene2 in genes_coexprime_dans_cellule){
        if(gene1 != gene2){
          df[gene1,gene2] = df[gene1,gene2] + 1
        }
        
      }
    }
  }
  library(igraph)
  
  # Conversion du data frame en matrice d'adjacence
  adj_matrix <- as.matrix(df)
  
  # Création du graphe avec poids (valeurs de la matrice)
  g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  # Récupérer les poids attribués automatiquement
  weights <- E(g)$weight
  
  # Normaliser les poids pour les utiliser comme largeur des arêtes
  E(g)$width <- weights/max(weights)*10 + 1  # Ajuste *5 selon ce que tu veux visuellement
  
  
  
  # graph % cellule coexpr/cellule totale
  nbr_cellules=length(colnames(counts_condition))
  label_ =  round(weights/nbr_cellules, 2)
  
  png(paste0(chemin_repertoire_output,"/proportion_coexpression_cellule_expression_binaire_",condition,"_.png"), width = 800, height = 600)
  
  plot(g,
       main = paste0("Proportion de coexpression des ",condition),
       layout = layout_in_circle,
       vertex.label.cex = 0.8,
       vertex.size = 15,
       edge.color = "grey",
       vertex.label.color = "black",
       vertex.color = "lightblue",
       edge.width = E(g)$width,
       edge.label = label_,
       edge.label.cex = 0.8,
       edge.label.color = "black")
  # Ajouter du texte en bas à droite
  usr <- par("usr")  # Limites du plot: c(xmin, xmax, ymin, ymax)
  x <- usr[2] - 0.05 * (usr[2] - usr[1])  # Légèrement à gauche du bord droit
  y <- usr[3] + 0.05 * (usr[4] - usr[3])  # Légèrement au-dessus du bas
  texte_a_rajouter = paste0("Nombre de cellules = ",nbr_cellules)
  text(x, y, labels = texte_a_rajouter, adj = c(1, 0), cex = 0.9, col = "black")
  dev.off()
  
  
  png(paste0(chemin_repertoire_output,"/nombre_coexpression_cellule_expression_binaire_",condition,"_.png"), width = 800, height = 600)
  # Affichage final
  plot(g,
       main = paste0("Nombre de coexpression des ",condition),
       layout = layout_in_circle,
       vertex.label.cex = 0.8,
       vertex.size = 15,
       edge.color = "grey",
       vertex.label.color = "black",
       vertex.color = "lightblue",
       edge.width = E(g)$width,
       edge.label = round(weights, 2),
       edge.label.cex = 0.8,
       edge.label.color = "black")
  # Ajouter du texte en bas à droite
  usr <- par("usr")  # Limites du plot: c(xmin, xmax, ymin, ymax)
  x <- usr[2] - 0.05 * (usr[2] - usr[1])  # Légèrement à gauche du bord droit
  y <- usr[3] + 0.05 * (usr[4] - usr[3])  # Légèrement au-dessus du bas
  texte_a_rajouter = paste0("Nombre de cellules = ",nbr_cellules)
  text(x, y, labels = texte_a_rajouter, adj = c(1, 0), cex = 0.9, col = "black")
  dev.off()
  #heatmap(as.matrix(df))
  
  
  
}

function_creation_graph_coexpression_specifique_cellule_expression_binaire_condition = function(counts,condition,noms_genes){
  
  counts_condition = function_selection_gene_dans_counts_AND_CONVERSION_NOM_DES_GENES(noms_genes,counts)
  liste_genes = rownames(counts_condition)
  df <- data.frame(matrix(0, nrow = length(liste_genes), ncol = length(liste_genes)))
  rownames(df) = liste_genes
  colnames(df) = liste_genes
  counts_condition = function_conditions_a_garder_dans_counts(counts_condition,c(condition)) #pour recuperer uniquement certains alignements (colonnes)
  counts_condition = function_conversion_expression_binaire(counts_condition)
  counts_condition = counts_condition[order(rownames(counts_condition)), order(colnames(counts_condition))] 
  
  liste_nom_cellules = colnames(counts_condition)
  liste_arete_possible_genes = c()
  for(gene1 in liste_genes){
    for(gene2 in liste_genes){
      if(gene1 != gene2){
        arete_gene1_gene2=c(gene1,gene2)
        arete_gene1_gene2 = sort(unlist(arete_gene1_gene2))
        liste_arete_possible_genes=append(liste_arete_possible_genes,list(arete_gene1_gene2))    
      }
    }
  }
  liste_arete_possible_genes = unique(liste_arete_possible_genes)
  
  nbr_cellules=length(colnames(counts_condition))
  for(arete in liste_arete_possible_genes){
    nbr11 = 0
    nbr00 = 0
    nbr01_10_11 = 0
    for(nom_cellule in liste_nom_cellules){
      nom_gene1 = arete[1]
      nom_gene2 = arete[2]
      
      gene1 = counts_condition[nom_gene1,nom_cellule]
      gene2 = counts_condition[nom_gene2,nom_cellule]
      
      if((gene1 == 1) && (gene2 == 1)){nbr11 = nbr11 +1}
      
      if((gene1 == 0) && (gene2 == 0)){nbr00 = nbr00 +1}
    }
    nbr01_10_11 = nbr_cellules - nbr00
    prop11sur011011 = nbr11/nbr01_10_11
    print(prop11sur011011)
    df[nom_gene1,nom_gene2] = prop11sur011011
  }
  df[] <- lapply(df, function(col) {
    if (is.numeric(col)) {
      col[is.na(col) | is.nan(col)] <- 0
    }
    return(col)
  })
  library(igraph)
  
  # Conversion du data frame en matrice d'adjacence
  adj_matrix <- as.matrix(df)
  heatmap(adj_matrix)
  # Création du graphe avec poids (valeurs de la matrice)
  g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  # Récupérer les poids attribués automatiquement
  weights <- E(g)$weight
  
  # Normaliser les poids pour les utiliser comme largeur des arêtes
  E(g)$width <- weights * 30 # Ajuste *5 selon ce que tu veux visuellement
  
  # graph % cellule coexpr/cellule totale
  nbr_cellules=length(colnames(counts_condition))
  label_ =  round(weights/nbr_cellules, 2)
  png(paste0(chemin_repertoire_output,"/proportion_coexpression_specifique_cellule_expression_binaire_",condition,"_.png"), width = 800, height = 600)
  # Affichage final
  plot(g,
       main = paste0("Nombre de coexpression spécifique des ",condition),
       layout = layout_in_circle,
       vertex.label.cex = 0.8,
       vertex.size = 15,
       edge.color = "grey",
       vertex.label.color = "black",
       vertex.color = "lightblue",
       edge.width = E(g)$width,
       edge.label = round(weights, 2),
       edge.label.cex = 0.8,
       edge.label.color = "black")
  # Ajouter du texte en bas à droite
  usr <- par("usr")  # Limites du plot: c(xmin, xmax, ymin, ymax)
  x <- usr[2] - 0.05 * (usr[2] - usr[1])  # Légèrement à gauche du bord droit
  y <- usr[3] + 0.05 * (usr[4] - usr[3])  # Légèrement au-dessus du bas
  texte_a_rajouter = paste0("Nombre de cellules = ",nbr_cellules)
  text(x, y, labels = texte_a_rajouter, adj = c(1, 0), cex = 0.9, col = "black")
  dev.off()
  #heatmap(as.matrix(df))
  
}



function_boxplot_taille_reads = function(counts,plot_titre,filename){
  
  png(filename, width = 800, height = 600)
  
  # Fonction pour calculer le N50 d'un vecteur numérique
  calc_N50 <- function(x) {
    x <- sort(x, decreasing = TRUE)
    cumsum_x <- cumsum(x)
    total <- sum(x)
    N50_val <- x[min(which(cumsum_x >= total / 2))]
    return(N50_val)
  }
  
  # Créons une figure avec boxplot et annotations
  boxplot(counts, 
          main = plot_titre,
          xlab = "Colonnes", ylab = "Valeurs",
          las = 2, col = "lightblue", border = "darkblue",
          outline = FALSE)  # Optionnel pour cacher les outliers
  
  # Ajout des statistiques pour chaque colonne
  for (i in seq_along(counts)) {
    col_data <- counts[[i]]
    
    # Calculs
    moy <- mean(col_data, na.rm = TRUE)
    med <- median(col_data, na.rm = TRUE)
    n50 <- calc_N50(col_data)
    quartiles <- quantile(col_data, probs = c(0.25, 0.75), na.rm = TRUE)
    
    # Position verticale pour texte
    ymax <- max(col_data, na.rm = TRUE)
    
    # Ajouter la moyenne (en rouge, triangle)
    points(i, moy, col = "red", pch = 17, cex = 1.5)
    
    # Ajouter la médiane (en bleu, cercle)
    points(i, med, col = "blue", pch = 19, cex = 1.5)
    
    # Ajouter le N50 (en vert, carré)
    points(i, n50, col = "darkgreen", pch = 15, cex = 1.5)
    
    # Ajouter les quartiles (en violet, croix)
    points(rep(i, 2), quartiles, col = "purple", pch = 4, cex = 1.5)
    
    # Ajouter texte descriptif au-dessus
    text(i, ymax, 
         labels = paste0("M:", round(moy, 1), 
                         "\nMd:", round(med, 1), 
                         "\nN50:", round(n50, 1)), 
         pos = 3, cex = 0.7)
  }
  
  # Légende
  legend("topright", legend = c("Moyenne", "Médiane", "N50", "Quartiles"),
         col = c("red", "blue", "darkgreen", "purple"),
         pch = c(17, 19, 15, 4), cex = 0.8)
  
  dev.off()
  
  
  
}


