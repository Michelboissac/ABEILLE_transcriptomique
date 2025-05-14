
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
