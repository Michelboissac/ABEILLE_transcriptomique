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