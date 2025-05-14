


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


function_CNN = function(counts,R2 = 0.99,filename){
  library("igraph")
  library("corrr")
  # Supprimer les gènes constants avant la corrélation
  counts_filtered <- counts[apply(counts, 1, function(x) sd(x) != 0), ]
  
  # Puis recalculer la corrélation
  cor_matrix <- correlate(t(counts_filtered), method = "pearson")

  # 2. Transformer la matrice en format utilisable
  cor_df <- stretch(cor_matrix) # passer en format "long"
  
  # 4. Supprimer les auto-corrélations (x == y)
  cor_df <- cor_df %>% filter(x != y)
  
  # 5. Réordonner les paires pour ne garder qu’une direction (évite doublons)
  cor_df <- cor_df %>%
    mutate(pair_id = paste0(pmin(x, y), "_", pmax(x, y))) %>%
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
