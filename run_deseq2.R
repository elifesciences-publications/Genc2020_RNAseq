runDeSeq2 <- function(samps2, outtag, gt, crossWT=F){
  require(DESeq2)
  require(tximport)
  require(ggplot2)
  require(cowplot)
  
  # Update outtag
  outtag = paste(outtag, gt, sep='_')
  
  
  # Load gene information for merging
  genes = read.delim('TxMatrix.BDGP6.V1.txt')
  tx2gene <- genes[,c('transcript_id', 'gene_id')]
  geneMat = genes[, c('gene_id', 'gene_type', 'gene_biotype', 'gene_symbol')]
  geneMat = geneMat[!duplicated(geneMat$gene_id),]
  
  # tximport
  txObj <- tximport(samps2$file_path,
                    type = "salmon",
                    tx2gene = tx2gene,
                    ignoreTxVersion=FALSE, txIn = TRUE, txOut = FALSE)
  
  # Construct DESeq model
  sampMat = sampMat[,c('batch', 'condition')]
  if (crossWT){
    print ('design =  ~ condition')
    dds <- DESeqDataSetFromTximport(txObj, sampMat, design =  ~ condition)
    } else {
      print ('design =  ~ batch + condition')
      dds <- DESeqDataSetFromTximport(txObj, sampMat, design =  ~ batch + condition)
    }
  
  
  # pre-filtering the dataset
  print (paste('[Progress] The number of genes used from the result table : ', nrow(dds), sep=''))
  dds <- dds[ rowSums(counts(dds)) != 0, ]
  # dds <- dds[rowSums(counts(dds) == 0) == 0,]
  print (paste('[Progress] The number of genes after removing low count genes: ', nrow(dds), sep=''))
  
  # add size factors
  print ('Perform estimateSizeFactors on dds.')
  dds <- estimateSizeFactors(dds)
  print ('Perform estimateDispersions on dds.')
  dds <- estimateDispersions(dds)
  
  # Distribution of gene expression by samples and gene types
  print ('Distribution of gene expression by samples and gene types.')
  require(tidyr)
  df = as.data.frame(counts(dds, normalize=T))
  df$gene_id <- as.character(do.call(rbind.data.frame, strsplit(rownames(df), '.', fixed = T))[[1]])
  df = df %>% gather(sample, count, -gene_id)
  genes1 = geneMat[,c('gene_id','gene_biotype')]
  genes1 = genes1[!duplicated(genes1),]
  df = merge(df, genes1, by='gene_id', all.x=T)
    
  df_unnorm = as.data.frame(counts(dds, normalize=F))
  df_unnorm$gene_id <- as.character(do.call(rbind.data.frame, strsplit(rownames(df_unnorm), '.', fixed = T))[[1]])
  df_unnorm = df_unnorm %>% gather(sample, count, -gene_id)
  genes1 = geneMat[,c('gene_id','gene_biotype')]
  genes1 = genes1[!duplicated(genes1),]
  df_unnorm = merge(df_unnorm, genes1, by='gene_id', all.x=T)
  
  df$type <- ifelse(df$gene_biotype == 'protein_coding', 'protein_coding', 'pseudogene')
  df$type <- factor(df$type, levels = c('protein_coding','pseudogene'))
  df_unnorm$type <- ifelse(df_unnorm$gene_biotype == 'protein_coding', 'protein_coding', 'pseudogene')
  df_unnorm$type <- factor(df_unnorm$type, levels = c('protein_coding','pseudogene'))
  
  p1 <- ggplot(df, aes(x = sample, y = log2(count), fill=type)) + 
    geom_boxplot(outlier.size = 0.25) + facet_wrap(~type, ncol = 1) + 
    theme_bw(base_size = 7) + 
    labs(x='', title = 'After estimateSizeFactors') +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    guides(color=guide_legend(override.aes=list(fill=NA)))
  
  p2 <- ggplot(df_unnorm, aes(x = sample, y = log2(count), fill=type)) + 
    geom_boxplot(outlier.size = 0.25) + facet_wrap(~type, ncol = 1) + 
    theme_bw(base_size = 7) + 
    labs(x='', title = 'Before estimateSizeFactors') +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    guides(color=guide_legend(override.aes=list(fill=NA)))
  
  p = plot_grid(p2, p1, labels = c("A", "B"), ncol=2)
  
  save_plot(paste("Figures/plot.counts_dist", outtag,"pdf", sep='.'), 
            p,
            base_height = 5*1, 
            base_width = 5*1)
  
  # Similarity between samples
  require(RColorBrewer)
  require(pheatmap)
  rld <- rlog(dds, blind = FALSE)
  sampleDists <- dist(t(assay(rld)))
  sampleDistMatrix <- as.matrix( sampleDists )
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors,
           cellwidth = 10, cellheight = 10, fontsize = 4, 
           fontsize_row = 4, fontsize_col = 4, height = 2.5,
           filename = paste("Figures/plot.sample_similarity", outtag, "png", sep='.')) 
  print ('[Progress] Similarity matrix is created and plotted.')
  
  # Remove batch effect using limma
  require(limma)
  vsd <- vst(dds)
  p1 <- plotPCA(vsd, "batch")
  p2 <- plotPCA(vsd, "condition")
  assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$batch)
  p3 <- plotPCA(vsd, "batch")
  p4 <- plotPCA(vsd, "condition")
  p = plot_grid(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), ncol=2)
  ggsave(paste("Figures/plot.pca_removedBatch", outtag,"png", sep='.'), 
         plot = p,  scale = 1, width = 4*2, height = 4*2, dpi = 300)
  print ('[Progress] PCA is created and plotted')
  
  # DEX
  dds <- DESeq(dds)
  resultsNames(dds)
  res <- results(dds, name='condition_1_vs_0')
  res <- lfcShrink(dds, coef = 'condition_1_vs_0')
  
  mcols(res, use.names = TRUE)
  summary(res)
  
  print ('Create a gene-based matrix.')
  data <- data.frame(gene_id = row.names(res),
                     baseMean = res$baseMean,
                     lfc = res$log2FoldChange,
                     lfcSE = res$lfcSE,
                     stat = res$stat,
                     pvalue = res$pvalue,
                     padj = res$padj)
  
  
  data <- na.omit(data)
  print ('[Progress] DEX model is created')  
  
  # fix the gene names
  genes1 = geneMat[,c('gene_id', 'gene_symbol', 'gene_biotype')]
  genes1 = genes1[!duplicated(genes1),]
  data$gene_id <- as.character(do.call(rbind.data.frame, strsplit(as.character(data$gene_id), '.', fixed = T))[[1]])
  data = merge(data, genes1, by='gene_id', all.x=TRUE)
  print ('[Progress] Fix the gene name for Ensembl Gene IDs')
  
  # Add changes (Fold Change = 1.5 and adjusted p-value < 0.05)
  print ('Add changes (Any Fold Change and adjusted p-value < 0.05)')
  data$change <- ifelse(data$padj < 0.05, 
                        ifelse(data$lfc > 0, 'Up', 'Down'), 'None')
  
  # MA plot
  pdf(paste("Figures/plot.MA", outtag,"pdf", sep='.'))
  print ('MA plot for protein coding genes.')
  par(mfrow = c(1,2))
  nGenes = sum(rownames(res) %in% geneMat[geneMat$gene_biotype=='protein_coding',]$gene_id & !is.na(res$padj))
  DESeq2::plotMA(res[rownames(res) %in% geneMat[geneMat$gene_biotype=='protein_coding',]$gene_id & !is.na(res$padj),], 
                 ylim = c(-5,5), alpha=0.05, main=paste('Protein Coding Genes', paste('(n=', nGenes, ')', sep=''), sep='\n'))
  abline(h=c(log2(1.5), log2(1/1.5) ),col="dodgerblue",lwd=2)
  
  print ('MA plot for non-coding genes.')
  nGenes = sum(rownames(res) %in% geneMat[geneMat$gene_biotype!='protein_coding',]$gene_id & !is.na(res$padj))
  DESeq2::plotMA(res[rownames(res) %in% geneMat[geneMat$gene_biotype!='protein_coding',]$gene_id,], 
                 ylim = c(-5,5), alpha=0.05, main=paste('Protein Coding Genes', paste('(n=', nGenes, ')', sep=''), sep='\n'))
  abline(h=c(log2(1.5), log2(1/1.5) ),col="dodgerblue",lwd=2)
  dev.off()
  
  # Volcano plot
  min_baseMean = 0
  plotMat = data[data$gene_biotype=='protein_coding' & data$baseMean >= min_baseMean,]
  n_up_dex = sum(plotMat$change=='Up')
  n_down_dex = sum(plotMat$change=='Down')
  
  p = ggplot(plotMat, 
             aes(x = lfc, y = -log10(padj) )) + 
    geom_point(aes(color = change), size = 0.75, alpha = 0.8, na.rm = T) + # Make dots bigger
    theme_bw(base_size = 7) + # change theme
    ggtitle(label = paste( paste(strsplit(gt, '_vs_', fixed = T)[[1]][1], 'vs.', strsplit(gt, '_vs_', fixed = T)[[1]][2], sep=' '), 
                          '\n(Up=', n_up_dex, ';',
                          ' Down=',n_down_dex, ')', sep=''
    )
    ) + # Add a title
    xlab(expression(log[2]("Het" / "WT"))) + # x-axis label
    ylab(expression(-log[10]("adjusted p-value"))) + # y-axis label
    xlim(-4.5, 4.5) + ylim(0,62) + 
    # geom_vline(xintercept = c(-1.5,1.5), colour = "darkgrey") + # Add cutoffs
    geom_hline(yintercept = 1.3, colour = "darkgrey", lwd=0.20) + # Add cutoffs
    geom_vline(xintercept = 0, colour = "black", lwd=0.20) + # Add 0 lines
    scale_color_manual(values = c("Up" = "#E64B35", 
                                  "Down" = "#3182bd", 
                                  "None" = "grey")) # change colors
  
  save_plot(paste("Figures/plot.volcano_ProteinCoding", outtag,"pdf", sep='.'), 
            p,
            base_height = 2.7*1, 
            base_width = 3.5*1)
  
  print ('return the volcano plot')
  return(p)
}
 