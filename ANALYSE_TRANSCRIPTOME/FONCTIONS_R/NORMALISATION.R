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
#######################################################################
NORMALISATION = function(counts, tableau_counts, normalisation = c("pas_de_normalisation","vst","vst_conditions", "vst_peu_de_genes","log_rpkm", "log_tpm", "rlog"),filename,liste_conditions_a_garder_etoile=NULL,longueurs_genes){
  normalisation <- match.arg(normalisation)
  if(normalisation=="pas_de_normalisation"){print("pas_de_normalisation")}
  if(normalisation=="log_rpkm"){counts = function_normalise_log_RPKM(counts,tableau_counts,longueurs_genes)}
  if(normalisation=="log_tpm"){counts = function_normalise_log_TPM(counts,tableau_counts,longueurs_genes)}
  if(normalisation=="rlog"){counts = function_normalisation_rlog(counts,tableau_counts)}
  if(normalisation=="vst"){counts = function_normalisation_VST(counts,tableau_counts,filename)}
  if(normalisation=="vst_peu_de_genes"){counts = function_normalisation_VST_peu_de_genes(counts,tableau_counts,filename)}
  if(normalisation=="vst_conditions"){counts = function_normalisation_VST_conditions(counts,tableau_counts,filename)}
  return(counts)
}
#######################################################################