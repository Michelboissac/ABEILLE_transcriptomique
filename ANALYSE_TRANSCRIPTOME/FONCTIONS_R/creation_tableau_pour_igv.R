
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