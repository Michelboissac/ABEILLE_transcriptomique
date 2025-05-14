
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
