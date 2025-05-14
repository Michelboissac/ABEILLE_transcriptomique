
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