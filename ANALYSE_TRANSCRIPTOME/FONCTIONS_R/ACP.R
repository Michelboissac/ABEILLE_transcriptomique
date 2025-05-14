function_ACP <- function(counts, plot_titre = "ACP", filename = "ACP_plot.png", width = 1200, height = 1000, res = 150){  #ACP :
  data_t <- t(counts)
  pca <- prcomp(data_t, scale. = F)
  summary(pca)
  
  png(filename, width = width, height = height, res = res)
  plot(pca$x[,1:2], col=1:nrow(data_t), pch=19, main = "plot_titre")
  text(pca$x[,1:2], labels=rownames(pca$x), pos=3)
  dev.off()
  
  
  
  
  

  
  
  
  # Voir les contributions aux composantes principales
 # contributions <- pca$rotation
  
  # Quelle variable contribue le plus à PC1 ?
  #abs(contributions[,1])  # valeurs absolues pour l'importance
  # Classement décroissant
 # importance_PC1 <- sort(abs(contributions[,1]), decreasing = TRUE)
 # print(importance_PC1)
  #boxplot(importance_PC1)
  
  #importance_PC2 <- sort(abs(contributions[,2]), decreasing = TRUE)
  #print(importance_PC2)
  
 # importance_PC3 <- sort(abs(contributions[,3]), decreasing = TRUE)
  #print(importance_PC3)

 # head(importance_PC1)
 # head(importance_PC2)
  
  #head(importance_PC3)
  #biplot(pca)
  
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
