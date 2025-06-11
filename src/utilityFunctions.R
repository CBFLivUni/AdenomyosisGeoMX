

# Utility function to summary the QC flags for a NanoStringGeoMxSet
QCSummary <- function(object) {
  # Collate QC Results
  QCResults <- protocolData(object)[["QCFlags"]]
  flag_columns <- colnames(QCResults)
  QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                           Warning = colSums(QCResults[, flag_columns]))
  
  QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
    ifelse(sum(x) == 0L, "PASS", "WARNING")
  })
  
  QC_Summary["TOTAL FLAGS", ] <-
    c(sum(QCResults[, "QCStatus"] == "PASS"),
      sum(QCResults[, "QCStatus"] == "WARNING"))
  
  
  warn_formatter <- formatter("span", 
                              style = x ~ style( "background-color" = ifelse(x > 0 , "yellow", "white")))
  
  pass_formatter <- formatter("span", 
                              style = x ~ style( "background-color" = ifelse(x > 0 , "lightgreen", "white")))
  
  QCTable <- formattable(QC_Summary, list(Pass=pass_formatter,Warning=warn_formatter),caption = "Summary of QC flags")
  
  return(list(QCResults=QCResults,QCTable=QCTable))
}


# Graphical summaries of QC statistics plot function- taken from the GeoMxWorkflows vignette
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = "segment",
                         thr = NULL,
                         scale_trans = NULL,
                         fill_cols=NULL) {
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = fill_by)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
    labs(x = annotation, y = "Segments, #", title = annotation)
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  if(!is.null(fill_cols)) {
    plt <- plt +
      scale_color_manual(values=fill_cols) +
      scale_fill_manual(values=fill_cols)
  }
  plt
}

#  Function to filter a geomxdataset to keep samples suitable for running a random effects model
# i.e contains multiple observations per slide
filterDataWithinSlide <- function(dataset){
  
  pheno <- as.data.frame(pData(dataset))%>% rownames_to_column("ID") 
  
  #keep samples where the staining is present on multiple slides per Disease.
  phenoFiltered <- pheno %>% group_by(Disease,slide) %>% filter(n_distinct(segment)>1) %>%
    group_by(Disease,segment) %>% filter(n() > 1)  %>% group_by(Disease) %>%
    mutate(segmentTypes=n_distinct(segment)) %>% group_by(Disease,slide) %>%
    filter(!segmentTypes > n_distinct(segment))
  
  ind <- which(pheno$ID %in% phenoFiltered$ID)
  
  return(dataset[,ind])
  
}

# Function to allow group split but keeping the names of the groups in the resulting list
# modified from https://github.com/tidyverse/dplyr/issues/4223#issuecomment-469269857
named_group_split <- function(.tbl, ...) {
  grouped <- group_by(.tbl, ...)
  names <- rlang::inject(paste(!!!group_keys(grouped), sep = "_"))
  
  grouped %>% 
    group_split() %>% 
    rlang::set_names(names)
}


# Function to plot a MAplot given the differential expression results and the negative probes
MAPlot <- function(diffExp,contrasts=1,FCColumn="Estimate",FDRColumn="FDR",negProbes="NegProbe-WTX") {
  
  diffExp$Color <- ifelse(diffExp[,FDRColumn]<=0.05,"yes","no")
    
  cutoff <- quantile(diffExp$meanExp, 0.9)
  diffExp$absFC <- abs(diffExp[,FCColumn])
  
  FCColumn <- sym(FCColumn)
  FDRColumn <- sym(FDRColumn)
  
  if(contrasts == 1){
    labelData <- diffExp %>% group_by(Subset) %>% arrange(!!FDRColumn) %>%
      filter(!!FDRColumn<0.01 & absFC >1 & meanExp > cutoff ) %>% slice_head(n=8)
  } else {
    labelData <- diffExp %>% group_by(Subset,Contrast) %>% arrange(!!FDRColumn) %>%
      filter(!!FDRColumn<0.01 & absFC >1 & meanExp > cutoff ) %>% slice_head(n=8)
  }
  
  
  
  g <- ggplot(subset(diffExp, !Gene %in% negProbes),
         aes(x = meanExp, y = !!FCColumn,
             color = Color, label = Gene)) +
    geom_hline(yintercept = c(0.5, -0.5), lty = "dashed") +
    scale_x_continuous(trans = "log2") +
    geom_point(alpha = 0.5) + 
    labs(y = "log2 Fold Change",
         x = "Mean Expression",
         color = "FDR ≤ 0.05") +
    scale_color_manual(values = c(yes = "dodgerblue",
                                  no = "grey")) +
    geom_text_repel(data = labelData,
                    size = 4, point.padding = 0.15, color = "black",
                    min.segment.length = .1, box.padding = .2, lwd = 2,max.overlaps = Inf,seed=42) +
    theme_bw(base_size = 16)
    
  
  if(contrasts == 1){
    g <- g +  facet_wrap(~Subset,ncol = 1)
  } else {
      g <- g + facet_wrap(~Subset + Contrast)
  }
  
  return(g)
}



data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}



# Function to get the results table for a limma contrast
getResultsDataFrame <- function(fit2, contrast, numerator, 
                                denominator) {
  data <- topTable(fit2, coef = contrast, number = Inf, 
                   sort.by = "P")
  #data <- independentFiltering(data, filter = data$AveExpr, objectType = "limma")
  #colnames(data) <- paste(paste(numerator, denominator, sep = "vs"), colnames(data), sep = "_")
  return(data)
}




plotBoxPlotSegment <- function(GOI, spatialData, elt="normmat", normalisation="GeoDiff", segment="") {
  
  pal  <- c("#2dc4cc", "#e0690d", "#7570B3", "#E7298A", "#66A61E", "#E6AB02")
  
  DF_ <- data.frame(t(assayDataElement(spatialData[GOI, ], elt = elt)), pData(spatialData))
  DF_ <- DF_ %>% rename("Gene" = all_of(GOI[1]), "Region" = "region")
  
  DF_ <- DF_ %>%
    mutate(segment=replace(segment, grepl(pattern = "Myometrium", x = Region), "Myo")) 
  DF_$Region <- gsub("Myometrium_near_adeno", "Myo 2", gsub("Myometrium_near_endo", "Myo 1", DF_$Region))
  DF_ <-DF_[str_detect(DF_$segment, all_of(segment)), ]
  
  plot <- ggplot(DF_, aes(y=Gene, x=Region))+
    geom_boxplot(width = 0.4, alpha = 0.3, na.rm = TRUE)+
    geom_point(aes(color = Region), position = ggplot2::position_jitterdodge(dodge.width = 0.60),
               alpha = 0.5,
               size = 3,
               stroke = 0,
               na.rm = TRUE)+
    theme_cowplot()+
    scale_color_manual(values=pal) +
    theme(axis.title.x=element_blank(),
          axis.text.x = element_text(angle = 40,  hjust=1))+
    labs(y = paste(GOI,"Expression"), x = "Region")+
    ggtitle(segment)
  
  return(plot)
  
}


# Function to produce boxplot for GOI #
# https://github.com/IndrajeetPatil/ggstatsplot/blob/main/R/ggbetweenstats.R

plotBoxPlot <- function(GOI, spatialData, elt="normmat", normalisation="GeoDiff") {
  
  pal  <- c("#E7298A", "#7570B3", "#e0690d", "#2dc4cc", "#66A61E", "#E6AB02")
  
  if (normalisation=="GeoDiff") {
    
    DF_ <- data.frame(t(assayDataElement(spatialData[GOI, ], elt = elt)), pData(spatialData))
    DF_ <- DF_ %>% rename("Gene" = all_of(GOI[1]), "Region" = "region")
    
  } else if (normalisation=="Q3") {
    
    DF_ <- data.frame(t(log2(assayDataElement(spatialData[GOI, ], elt = elt))), pData(spatialData))
    DF_ <- DF_ %>% rename("Gene" = all_of(GOI[1]), "Region" = "region")
    
  } else if (normalisation=="TMM") {
    
    d <- DGEList(na.omit(assayDataElement(spatialData, elt = elt)))
    d <- calcNormFactors(d)
    cpm <- cpm(d, log=TRUE)
    DF_ <- data.frame(cpm[GOI,], pData(spatialData))
    DF_ <- DF_ %>% rename("Gene" = "cpm.GOI...", "Region" = "region")
    
  } else {
    
    print("Not one of the available normalisation methods, please choose GeoDiff, Q3 or TMM...")
    
    
  }  
  
  DF_ <- DF_ %>%
    mutate(segment=replace(segment, grepl(pattern = "Myometrium", x = Region), "Myo")) 
  
  DF_$segment <- factor(DF_$segment, ordered = TRUE, levels = c("PanCK", "Stroma", "CD45", "CD31", "Myo"))
  
  DF_$Region <- gsub("Myometrium_near_adeno", "Myo 2", gsub("Myometrium_near_endo", "Myo 1", DF_$Region))
  DF_$Region <- factor(DF_$Region, ordered = TRUE, levels = c("LE", "Functionalis", "Basalis", "Adenomyosis", "Myo 1", "Myo 2"))
  
  plot <- ggplot(DF_, aes(y=Gene, x=Region))+
    geom_boxplot(width = 0.4, alpha = 0.3, na.rm = TRUE)+
    geom_point(aes(color = Region), position = ggplot2::position_jitterdodge(dodge.width = 0.60),
               alpha = 0.5,
               size = 3,
               stroke = 0,
               na.rm = TRUE)+
    facet_wrap(~segment, ncol = 5, scales = "free_x")+
    theme_cowplot()+
    scale_color_manual(values=pal) +
    theme(axis.title.x=element_blank(),
          axis.text.x = element_text(angle = 40,  hjust=1))+
    labs(y = paste(GOI,"Expression"), x = "Region") 
  
  return(plot)
}



# Function for PCA correlation plot 

library(PCAtools)

eigencorplotPCA <- function(x, metavars="") {
  
  # Wrapper around the PCAtools eigencorplot to use preferred aesthetics
  # Takes as input PCA object from PCAtools and a vector of metavars to measure correlation with
  # The metavars have to be present in the original metadata provided to the PCA object
  
  use_cex <- 16/16
  
  corplot <- PCAtools::eigencorplot(x,
                                    metavars = metavars,
                                    col = c( "blue2", "blue1", "black", "red1", "red2"),
                                    colCorval = "white",
                                    scale = TRUE,
                                    main = "",
                                    plotRsquared = FALSE,
                                    cexTitleX= use_cex,
                                    cexTitleY= use_cex,
                                    cexLabX = use_cex,
                                    cexLabY = use_cex,
                                    cexMain = use_cex,
                                    cexLabColKey = use_cex,
                                    cexCorval = use_cex)
  
  return(corplot)
  
}



#Function to create PCA plot #

library(PCAtools)
library(viridis)
library(ggrepel)

plotPCA <- function(x, PCs="", colours=NULL, colour.data=NULL, shape.data=NULL, sample.lab=F, colour.lab="", shape.lab="", borderWidth = 0.8) {
  
  df_out <- data.frame(PC1=x$rotated[,PCs[1]], PC2=x$rotated[,PCs[2]])
  row.names(df_out) <- row.names(x$rotated)
  
  
  if(!is.null(colour.data) & !is.null(shape.data)){
    
    if(!is.null(colours)){
      
      if(is.factor(colour.data)){
        
        pca_plot <- ggplot(df_out, aes(x=PC1, y=PC2, colour = colour.data, shape=shape.data)) +
          geom_point(size = 3.5) +
          scale_color_manual(values=colours) +
          scale_shape_manual(values=c(19,15,17,18)) +
          labs(colour = colour.lab, shape = shape.lab)
        
        
      } else if(is.numeric(colour.data)){
        
        pca_plot <- ggplot(df_out, aes(x=PC1, y=PC2, colour = colour.data, shape=shape.data)) +
          geom_point(size = 3.5) +
          scale_color_gradientn(colours = colours) +
          scale_shape_manual(values=c(19,15,17,18)) +
          labs(colour = colour.lab, shape = shape.lab)
        
      } else {
        
        stop("Error: Colour data must be type 'numeric' or 'factor'.")
        
      }
      
    } else {
      
      if(is.factor(colour.data)){
        
        pca_plot <- ggplot(df_out, aes(x=PC1, y=PC2, colour = colour.data, shape=shape.data)) +
          geom_point(size = 3.5) +
          scale_colour_brewer(palette = "Set1") +
          scale_shape_manual(values=c(19,15,17,18)) +
          labs(colour = colour.lab, shape = shape.lab)
        
      } else if(is.numeric(colour.data)){
        
        pca_plot <- ggplot(df_out, aes(x=PC1, y=PC2, colour = colour.data, shape=shape.data)) +
          geom_point(size = 3.5) +
          scale_color_viridis(option = "plasma") +
          scale_shape_manual(values=c(19,15,17,18)) +
          labs(colour = colour.lab, shape = shape.lab)
        
      } else {
        
        stop("Error: Colour data must be type 'numeric' or 'factor'.")
        
      }
      
    }
    
    
    
  } else if(!is.null(colour.data) & is.null(shape.data)) {
    
    if(!is.null(colours)){
      
      if(is.factor(colour.data)){
        
        pca_plot <- ggplot(df_out, aes(x=PC1, y=PC2, colour = colour.data)) +
          geom_point(size = 3.5) +
          scale_color_manual(values=colours) +
          labs(colour = colour.lab)
        
      } else if(is.numeric(colour.data)){
        
        pca_plot <- ggplot(df_out, aes(x=PC1, y=PC2, colour = colour.data)) +
          geom_point(size = 3.5) +
          scale_color_gradientn(colours = colours) +
          labs(colour = colour.lab)
        
      } else {
        
        stop("Error: Colour data must be type 'numeric' or 'factor'.")
        
      }
      
    } else {
      
      if(is.factor(colour.data)){
        
        pca_plot <- ggplot(df_out, aes(x=PC1, y=PC2, colour = colour.data)) +
          geom_point(size = 3.5) +
          scale_colour_brewer(palette = "Set1") +
          labs(colour = colour.lab)
        
      } else if(is.numeric(colour.data)){
        
        pca_plot <- ggplot(df_out, aes(x=PC1, y=PC2, colour = colour.data)) +
          geom_point(size = 3.5) +
          scale_color_viridis(option = "plasma") +
          labs(colour = colour.lab)
        
      } else {
        
        stop("Error: Colour data must be type 'numeric' or 'factor'.")
        
      }
      
    }
    
    
    
  } else {
    
    pca_plot <- ggplot(df_out, aes(x=PC1, y=PC2)) +
      geom_point(size = 3.5)
    
  }
  
  
  pca_plot <- pca_plot +
    xlab(paste0('PC', PCs[1], ': ', round(as.numeric(x$variance[PCs[1]])), '% expl.var')) +
    ylab(paste0('PC', PCs[2], ': ', round(as.numeric(x$variance[PCs[2]])), '% expl.var')) +
    #lims(x= c(min(df_out[,"PC1"])-20, max(df_out[,"PC1"])+20), y = c(min(df_out[,"PC2"])-20, max(df_out[,"PC2"])+20)) +
    theme_cowplot()
  
  
  if(sample.lab){
    
    pca_plot <- pca_plot + geom_text_repel(aes(label = row.names(df_out), size = NULL, color = NULL),
                                           check_overlap = T,
                                           size = 3.0
    )
    
  }
  
  return(pca_plot)
  
}






#Function to run RUV-seq with DEseq2

RUV.total <- function(raw,pData,fData,cIdx,k,exclude = NULL){
  
  library(RUVSeq)
  library(DESeq2)
  library(limma)
  library(matrixStats)
  
  
  fData = fData[rownames(raw),]
  int = intersect(rownames(raw),rownames(fData))
  fData = fData[int,]
  raw = raw[int,]
  
  ## USE DESEQ2 FORMULATION TO INTEGRATE RAW EXPRESSION
  ## PDATA AND FDATA
  set <- newSeqExpressionSet(as.matrix(round(raw)),
                             phenoData=pData,
                             featureData=fData)
  
  ## UPPER QUANTILE NORMALIZATION (BULLARD 2010)
  set <- betweenLaneNormalization(set, which="upper")
  
  ## RUVg USING HOUSKEEPING GENES
  set <- RUVg(set, cIdx, k=k)
  dds <- DESeqDataSetFromMatrix(counts(set),colData=pData(set),design=~1)
  rowData(dds) <- fData
  
  ## SIZE FACTOR ESTIMATIONS
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersionsGeneEst(dds)
  cts <- counts(dds, normalized=TRUE)
  disp <- pmax((rowVars(cts) - rowMeans(cts)),0)/rowMeans(cts)^2
  mcols(dds)$dispGeneEst <- disp
  dds <- estimateDispersionsFit(dds, fitType="mean")
  
  ## TRANSFORMATION TO THE LOG SPACE WITH A VARIANCE STABILIZING TRANSFORMATION
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  mat <- assay(vsd)
  
  ## REMOVE THE UNWANTED VARIATION ESTIMATED BY RUVg
  covars <- as.matrix(colData(dds)[,grep("W",colnames(colData(dds))),drop=FALSE])
  mat <- removeBatchEffect(mat, covariates=covars)
  assay(vsd) <- mat
  
  return(vsd = vsd)

}


# Function to create volcano plot

plotVolcano <- function(x, title="", top.genes = NULL, point.size = 0.9, lab.size=3, quadrants=F, foldChangeColumn="logFC",FDRColumn="adj.P.Val") {
  
  x <- x %>%
    rownames_to_column(var = "Genes") %>%
    # Rename columns to be consistent
    rename("logFC" = all_of(foldChangeColumn), "adj.P.Val" = all_of(FDRColumn)) %>%
    mutate(
      Expression = case_when(logFC >= 0.5 & adj.P.Val <= 0.05 ~ "Up-regulated",
                             logFC <= -0.5 & adj.P.Val <= 0.05 ~ "Down-regulated",
                             TRUE ~ "Unchanged")
    )
  
  x$Expression <- factor(x$Expression, levels=c("Down-regulated", "Unchanged", "Up-regulated"))
  myColors <- c("dodgerblue3", "black", "firebrick3")
  names(myColors) <- levels(x$Expression)

  
  
  volcano_plot <- ggplot(x, aes(logFC, -log(adj.P.Val, 10))) + # -log10 conversion
    geom_point(aes(color = Expression), size = point.size) +
    xlab(expression("log"[2]*"FC")) +
    ylab(expression("-log"[10]*"FDR")) + 
    scale_color_manual(name = "Expression",values = myColors) +
    theme_cowplot() + 
    ggtitle(title) + 
    theme(legend.position="none")
  
  
  if(quadrants){
    
    volcano_plot <- volcano_plot +
      geom_hline(yintercept = -log10(0.05),
                 linetype = "dashed") + 
      geom_vline(xintercept = c(-0.4, 0.4),
                 linetype = "dashed")  
    
  }
  
  
  if(!is.null(top.genes)){
    
    top_genes <- bind_rows(
      x %>% 
        filter(Expression == 'Up-regulated') %>% 
        arrange(adj.P.Val, desc(abs(logFC))) %>% 
        head(top.genes),
      x %>% 
        filter(Expression == 'Down-regulated') %>% 
        arrange(adj.P.Val, desc(abs(logFC))) %>% 
        head(top.genes)
    )
    
    
    volcano_plot <- volcano_plot +
      #geom_label_repel(data = top_genes, mapping = aes(logFC, -log(adj.P.Val,10), label = Genes), size = lab.size, force = 1,
      #                 nudge_y = 0.2)
      geom_text_repel(data = top_genes, aes(label=Genes),max.overlaps = Inf,min.segment.length = 0,seed = 42) 
    
  }
  
  
  return(volcano_plot)
  
  
}




# Function to create barplot and table of # of DE genes per contrast 
# Partially adapted from code provided by Lauren Mee (CBF)
# res_list is a list of results generated from a series of contrasts of interest
# 'foldChangeColumn' and 'FDRColumn' are the relevant columns containing the log2FC and adjusted p.values for each object in the results list
# stains is a character vector relevant for spatial transcriptomics experiments
#   It contains the stains you want to facet the plot by
#   Each stain must be seperated by "|" e.g. "CD45|CD31|PanCK|Stroma"
# Example:
#   DEResPlot(limma_GeoDiff_results_table, foldChangeColumn="logFC", FDRColumn="adj.P.Val", stains=c("CD45|CD31|PanCK|Stroma"))

DEResPlot <- function(res_list, foldChangeColumn="logFC", FDRColumn="adj.P.Val", stains=NULL, title=""){
  
  # Create dataframe from results list
  for (i in 1:length(res_list)){
    contrast <- names(res_list)[[i]]
    
    res_list[[i]] <- res_list[[i]] %>%
      # Rename columns to be consistent
      rename("logFC" = all_of(foldChangeColumn), "adj.P.Val" = all_of(FDRColumn)) %>%
      # Record whether genes are up or down regulated
      mutate(Direction = ifelse(logFC > 0, "Up", "Down"), Contrast = contrast) %>%
      rownames_to_column(var = "Genes")
  }
  
  # Combine all results
  res <- bind_rows(res_list) %>%
    # Convert contrast to factor so any cases of there being no sig genes
    # Are not lost once genes are filtered by significance
    mutate(Contrast = as.factor(Contrast)) %>%
    # Filter by signficance
    filter(adj.P.Val < 0.05) %>%
    # Make note of numbers of DEGs per contrast per direction of expression
    select(Contrast, Direction) %>%
    table(.)
  
  # Make a dataframe from the above information
  results <- data.frame("Up" = res[,"Up"],
                        "Down" = res[,"Down"]) %>%
    rownames_to_column(var = "Contrast")
  
  
  # Create barplots from the dataframe containing number of DE genes
  # If no stain information is provided all contrasts will be plotted together
  if(is.null(stains)){
    
    plotdf <- suppressMessages(reshape2::melt(results)) %>%
      filter(variable != "Not_Significant")
    for (i in 1:nrow(plotdf)){
      if(plotdf$variable[i] == "Down"){
        if(plotdf$value[i] > 0){
          plotdf$value[i] <- plotdf$value[i] * -1
        }
      }
    }
    
    p <- ggplot(plotdf, aes(x = Contrast, y = value)) +
      geom_bar(stat = "identity", aes(fill = variable)) +
      theme_bw(base_size = 13) +
      labs(y = "Number Significantly DE Genes",
           fill = "Expression Direction") +
      scale_fill_manual(values = c("Up" = "firebrick3",
                                   "Down" = "dodgerblue3")) +
      theme(legend.position = "bottom") +
      ggtitle(title)
  }
  
  # If stain information is provided the plots will be faceted depending on the stain
  else{
    
    # Add new column based off stain
    plotdf <- results %>% 
      mutate("Stain" = str_extract(.$Contrast, stains)) 
    
    # Remove rows that don't have any genes up/down regulated or don't have a matching stain
    plotdf <- suppressMessages(reshape2::melt(plotdf)) %>%
      filter(variable != "Not_Significant", Stain !="NA")
    
    # Make number of down-regulated genes negative in value
    for (i in 1:nrow(plotdf)){
      if(plotdf$variable[i] == "Down"){
        if(plotdf$value[i] > 0){
          plotdf$value[i] <- plotdf$value[i] * -1
        }
      }
    }
    
    # Tidy contrasts
    plotdf$Contrast <- str_remove_all(plotdf$Contrast, paste0(stains, "|_"))
    
    # Create plot
    p <- ggplot(plotdf, aes(x = Contrast, y = value)) +
      geom_bar(stat = "identity", aes(fill = variable)) +
      theme_bw(base_size = 13) +
      labs(y = "Number Significantly DE Genes",
           fill = "Expression Direction",
           x = "") +
      scale_fill_manual(values = c("Up" = "firebrick3",
                                   "Down" = "dodgerblue3")) +
      #theme(legend.position = "bottom",
      #      axis.text.x = element_text(angle = 15, vjust = 0.7)) +
      theme(legend.position = "bottom") +
      facet_wrap(~ Stain, ncol = 1, scales="free")
  }
  
  return(list(DEplot = p, DEtable = results))
  
}



# Function to check correlation between fold changes between models 

modelCorrelation <- function(res_list, filter.sig=T, by.pval=F, title="", subtitle=""){
  res_list <- list()
  
  for (i in 1:length(all_results)){
    res_list[[i]] <- do.call(rbind, all_results[[i]]) %>% mutate("gene_name" = gsub('.*\\.', '', row.names(.)))
    
    # Rename columns to be consistent with limma column names
    colnames(res_list[[i]])[which(names(res_list[[i]]) == "log2FC")] <- "logFC"
    colnames(res_list[[i]])[which(names(res_list[[i]]) == "adjp")] <- "adj.P.Val"
    
    
    # Choose whether to filter to only include significant genes
    if(filter.sig){
      res_list[[i]] <- res_list[[i]][res_list[[i]]$adj.P.Val <0.05,]
    }
    
    
    # Either look at the correlation between fold-change or p-value
    if(by.pval){
      res_list[[i]] <- res_list[[i]] %>% rownames_to_column(var = "Names") %>% pull(adj.P.Val, Names)
    
    } else {
      res_list[[i]] <- res_list[[i]] %>% rownames_to_column(var = "Names") %>% pull(logFC, Names)
      
    }
    
  }
  
  res_list <- t(tibble(V = res_list) %>% 
                  unnest_wider(V, names_sep = ""))
  
  
  colnames(res_list) <- names(all_results)
  
  corMat <- cor(res_list,use = "pairwise.complete.obs")
  p <- ggcorrplot::ggcorrplot(corMat,type = "upper",ggtheme = theme_cowplot(),lab=TRUE, colors = c("#6D9EC1", "white", "firebrick3"))
  p <- p + ggtitle(title, subtitle = subtitle)
  return(p)

}


# Function to create an enrichment dataframe from a list of contrasts containing DE results 
# The genes are ranked by signed -log10(pvalue)

compareClustersDF <- function(res_list, pvalueCutoff = 1, gmt=gmt){
  
  # Create empty list to store ranked genes
  res_list_ranked_pval <- list()
  
  # Iterate through results to generate gene rank metric per contrast
  for(i in 1:length(res_list)){
    
    # Rename columns to be consistent with limma column names
    colnames(res_list[[i]])[which(names(res_list[[i]]) == "log2FC")] <- "logFC"
    colnames(res_list[[i]])[which(names(res_list[[i]]) == "adjp")] <- "adj.P.Val"
    colnames(res_list[[i]])[which(names(res_list[[i]]) == "pvalue")] <- "P.Value"
    
    ranked_list <- res_list[[i]] %>% mutate(rank_metric = -log10(P.Value) * sign(logFC)) %>%
      arrange(desc(rank_metric)) %>%
      rownames_to_column() %>%
      pull(rank_metric, rowname)
    
    ranked_list[!is.finite(ranked_list)] <- NA
    ranked_list <- na.omit(ranked_list)
    res_list_ranked_pval[[i]] <- ranked_list
    
  }
  
  # Name list of ranked genes to match contrasts
  names(res_list_ranked_pval) <- names(res_list)
  
  
  # Pass list of ranked FC values to clusterProfiler's compareCluster function
  # This carries out GSEA for each contrast and adjusts p-values accordingly
  enrichment_results <- compareCluster(geneCluster = res_list_ranked_pval,
                                       fun = GSEA,
                                       exponent = 1,
                                       minGSSize = 8,
                                       maxGSSize = 500,
                                       pvalueCutoff = pvalueCutoff,
                                       pAdjustMethod = "BH",
                                       TERM2GENE = gmt, 
                                       by = "fgsea",
                                       eps = 0)
  
  
  # Return dataframe of enrichment results
  enrichment_results_df <- enrichment_results@compareClusterResult
  return(enrichment_results_df)
  
}


# Function to create an enrichment object from a list of contrasts containing DE results 
# The genes are ranked by signed -log10(pvalue)

compareClustersGSEA <- function(res_list, pvalueCutoff = 0.05, gmt=gmt){
  
  # Create empty list to store ranked genes
  res_list_ranked_pval <- list()
  
  # Iterate through results to generate gene rank metric per contrast
  for(i in 1:length(res_list)){
    
    # Rename columns to be consistent with limma column names
    colnames(res_list[[i]])[which(names(res_list[[i]]) == "log2FC")] <- "logFC"
    colnames(res_list[[i]])[which(names(res_list[[i]]) == "adjp")] <- "adj.P.Val"
    colnames(res_list[[i]])[which(names(res_list[[i]]) == "pvalue")] <- "P.Value"
    
    ranked_list <- res_list[[i]] %>% mutate(rank_metric = -log10(P.Value) * sign(logFC)) %>%
      arrange(desc(rank_metric)) %>%
      rownames_to_column() %>%
      pull(rank_metric, rowname)
    
    ranked_list[!is.finite(ranked_list)] <- NA
    ranked_list <- na.omit(ranked_list)
    res_list_ranked_pval[[i]] <- ranked_list
    
  }
  
  # Name list of ranked genes to match contrasts
  names(res_list_ranked_pval) <- names(res_list)
  
  
  # Pass list of ranked FC values to clusterProfiler's compareCluster function
  # This carries out GSEA for each contrast and adjusts p-values accordingly
  enrichment_results <- compareCluster(geneCluster = res_list_ranked_pval,
                                       fun = GSEA,
                                       exponent = 1,
                                       minGSSize = 8,
                                       maxGSSize = 500,
                                       pvalueCutoff = pvalueCutoff,
                                       pAdjustMethod = "BH",
                                       TERM2GENE = gmt, 
                                       by = "fgsea",
                                       eps = 0)
  
  
  # Return enrichment results
  return(enrichment_results)
  
}



# Function to create an enrichment object from a list of contrasts containing DE results 
# The genes are ranked by signed -log10(pvalue)

compareClustersGO <- function(res_list, pvalueCutoff = 0.05, ont="ALL", OrgDb = org.Hs.eg.db){
  
  # Create empty list to store ranked genes
  res_list_ranked_pval <- list()
  
  # Iterate through results to generate gene rank metric per contrast
  for(i in 1:length(res_list)){
    
    # Rename columns to be consistent with limma column names
    colnames(res_list[[i]])[which(names(res_list[[i]]) == "log2FC")] <- "logFC"
    colnames(res_list[[i]])[which(names(res_list[[i]]) == "adjp")] <- "adj.P.Val"
    colnames(res_list[[i]])[which(names(res_list[[i]]) == "pvalue")] <- "P.Value"
    
    ranked_list <- res_list[[i]] %>% mutate(rank_metric = -log10(P.Value) * sign(logFC)) %>%
      arrange(desc(rank_metric)) %>%
      rownames_to_column() %>%
      pull(rank_metric, rowname)
    
    ranked_list[!is.finite(ranked_list)] <- NA
    ranked_list <- na.omit(ranked_list)
    res_list_ranked_pval[[i]] <- ranked_list
    
  }
  
  # Name list of ranked genes to match contrasts
  names(res_list_ranked_pval) <- names(res_list)
  
  
  # Pass list of ranked FC values to clusterProfiler's compareCluster function
  # This carries out GSEA for each contrast and adjusts p-values accordingly
  enrichment_results <- compareCluster(geneCluster = res_list_ranked_pval,
                                       fun = gseGO,
                                       OrgDb = org.Hs.eg.db,
                                       keyType = 'SYMBOL',
                                       minGSSize = 8,
                                       maxGSSize = 500,
                                       pvalueCutoff = pvalueCutoff,
                                       pAdjustMethod = "BH",
                                       ont=ont)
  
  
  # Return enrichment results
  return(enrichment_results)
  
}

#  Create dotplot 

createDotplot <- function(res_list, title){
clusterProfiler::dotplot(res_list, color="enrichmentScore") +
  scale_color_gradientn(colours = c("#6D9EC1", "white", "firebrick3")) +
  theme_hc() +
  ggtitle(title)
}



# Function to run nichenet analysis 
runNicheNetAnalysis <- function(expressed_genes_sender, expressed_genes_receiver, geneset_oi, senderName, receiverName, organism="human"){
  
  if(organism == "human"){
    lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
    lr_network = lr_network %>% mutate(from = make.names(from), to = make.names(to))
    ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
    colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% make.names()
    rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% make.names()
  } else if(organism == "mouse"){
    lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
    lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor) %>% mutate(ligand = make.names(ligand), receptor = make.names(receptor))
    ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
    colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% make.names()
    rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% make.names()
  }
  
  background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
  
  # If wanted, users can remove ligand-receptor interactions that were predicted based on protein-protein interactions and only keep ligand-receptor interactions that are described in curated databases. To do this: uncomment following line of code:
  # lr_network = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
  
  ligands = lr_network %>% pull(from) %>% unique()
  expressed_ligands = intersect(ligands,expressed_genes_sender)
  
  receptors = lr_network %>% pull(to) %>% unique()
  expressed_receptors = intersect(receptors,expressed_genes_receiver)
  lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
  head(lr_network_expressed)
  
  potential_ligands = lr_network_expressed %>% pull(from) %>% unique()
  head(potential_ligands)
  
  ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
  ligand_activities %>% arrange(-pearson) 
  
  best_upstream_ligands = ligand_activities %>% top_n(10, aupr) %>% arrange(-aupr) %>% pull(test_ligand)
  
  # Changed histogram to be aupr score instead
  p_hist_lig_activity = ggplot(ligand_activities, aes(x=aupr)) + 
    geom_histogram(color="black", fill="darkorange")  + 
    # geom_density(alpha=.1, fill="orange") +
    geom_vline(aes(xintercept=min(ligand_activities %>% top_n(10, aupr) %>% pull(aupr))), color="red", linetype="dashed", size=1) + 
    labs(x="ligand activity (AUPR)", y = "# ligands") +
    theme_classic()
  p_hist_lig_activity
  
  active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 100) %>% bind_rows()
  
  active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = na.omit(active_ligand_target_links_df), ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)
  
  order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
  order_targets = active_ligand_target_links_df$target %>% unique()
  vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
  
  p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot(senderName,paste0("Differentially expressed genes in ",receiverName), color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))
  
  p_ligand_target_network
  
  
  # get the ligand-receptor network of the top-ranked ligands
  lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
  best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
  
  # get the weights of the ligand-receptor interactions as used in the NicheNet model
  # weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
  weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds?download=1"))
  lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
  
  # convert to a matrix
  lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
  lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
  
  # perform hierarchical clustering to order the ligands and receptors
  dist_receptors = dist(lr_network_top_matrix, method = "binary")
  hclust_receptors = hclust(dist_receptors, method = "ward.D2")
  order_receptors = hclust_receptors$labels[hclust_receptors$order]
  
  dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
  hclust_ligands = hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
  
  vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
  p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot(senderName,receiverName, color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
  p_ligand_receptor_network
  
  ligand_pearson_matrix = ligand_activities %>% dplyr::select(aupr) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
  
  vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("aupr")
  p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot(senderName,"Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "aupr \ntarget gene prediction ability)")
  p_ligand_pearson
  
  return(list(ligandActivityPlot=p_hist_lig_activity,ligandTargetLinksTable=active_ligand_target_links_df,ligandActivityTable=ligand_activities,ligandTarget=p_ligand_target_network,ligandReceptor=p_ligand_receptor_network,ligandPearson=p_ligand_pearson))
  
  
}


# Function to produce circos plot from the outputs of a nichenet analysis
# Code adapted from multinichenetR circos plot
# https://github.com/saeyslab/multinichenetr/blob/main/R/plotting.R
# circos_links_oi should be a dataframe containing the columns, ligand, target, weight, id, sender & receiver
# id should be a column that serves as a unique identifier for each interaction 
#  
produceCircosPlot = function(circos_links_oi, colors_sender, colors_receiver, title=""){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  requireNamespace("circlize")
  
  
  # Link each cell type to a color
  grid_col_tbl_ligand = tibble::tibble(sender = colors_sender %>% names(), color_ligand_type = colors_sender)
  grid_col_tbl_receptor = tibble::tibble(receiver = colors_receiver %>% names(), color_receptor_type = colors_receiver)
  
  
  # Make the plot for condition of interest - title of the plot
  title = title
  
  
  # Rename 'target' column to receptor
  circos_links = circos_links_oi %>% dplyr::rename(receptor = "target")
  df = circos_links
  
  
  # Each pair of ligand-receptors needs to be unique for circos plot to work
  # Code to make each pair unique by adding spaces after name
  ligand.uni = unique(df$ligand)
  for (i in 1:length(ligand.uni)) {
    df.i = df[df$ligand == ligand.uni[i], ]
    sender.uni = unique(df.i$sender)
    for (j in 1:length(sender.uni)) {
      df.i.j = df.i[df.i$sender == sender.uni[j], ]
      df.i.j$ligand = paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
      df$ligand[df$id %in% df.i.j$id] = df.i.j$ligand
    }
  }
  receptor.uni = unique(df$receptor)
  for (i in 1:length(receptor.uni)) {
    df.i = df[df$receptor == receptor.uni[i], ]
    receiver.uni = unique(df.i$receiver)
    for (j in 1:length(receiver.uni)) {
      df.i.j = df.i[df.i$receiver == receiver.uni[j], ]
      df.i.j$receptor = paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
      df$receptor[df$id %in% df.i.j$id] = df.i.j$receptor
    }
  }
  
  intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))
  
  while(length(intersecting_ligands_receptors) > 0){
    df_unique = df %>% dplyr::filter(!receptor %in% intersecting_ligands_receptors)
    df_duplicated = df %>% dplyr::filter(receptor %in% intersecting_ligands_receptors)
    df_duplicated = df_duplicated %>% dplyr::mutate(receptor = paste(" ",receptor, sep = ""))
    df = dplyr::bind_rows(df_unique, df_duplicated)
    intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))
  }
  
  circos_links = df
  
  # Link ligands/Receptors to the colors of senders/receivers
  circos_links = circos_links %>% dplyr::inner_join(grid_col_tbl_ligand) %>% dplyr::inner_join(grid_col_tbl_receptor)
  links_circle = circos_links %>% dplyr::distinct(ligand,receptor, weight)
  ligand_color = circos_links %>% dplyr::distinct(ligand,color_ligand_type)
  grid_ligand_color = ligand_color$color_ligand_type %>% magrittr::set_names(ligand_color$ligand)
  receptor_color = circos_links %>% dplyr::distinct(receptor,color_receptor_type)
  grid_receptor_color = receptor_color$color_receptor_type %>% magrittr::set_names(receptor_color$receptor)
  grid_col =c(grid_ligand_color,grid_receptor_color)
  
  # give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
  transparency = circos_links %>% dplyr::mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% dplyr::mutate(transparency = 1-weight) %>% .$transparency
  
  # Define order of the ligands and receptors and the gaps
  ligand_order = circos_links_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
    ligands = circos_links %>% dplyr::filter(sender == sender_oi) %>%  dplyr::arrange(ligand) %>% dplyr::distinct(ligand)
  }) %>% unlist()
  
  receptor_order = circos_links_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
    receptors = circos_links %>% dplyr::filter(receiver == receiver_oi) %>%  dplyr::arrange(receptor) %>% dplyr::distinct(receptor)
  }) %>% unlist()
  
  order = c(ligand_order,receptor_order)
  
  width_same_cell_same_ligand_type = 0.275
  width_different_cell = 3
  width_ligand_receptor = 9
  width_same_cell_same_receptor_type = 0.275
  
  sender_gaps = circos_links_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
    sector = rep(width_same_cell_same_ligand_type, times = (circos_links %>% dplyr::filter(sender == sender_oi) %>% dplyr::distinct(ligand) %>% nrow() -1))
    gap = width_different_cell
    return(c(sector,gap))
  }) %>% unlist()
  sender_gaps = sender_gaps[-length(sender_gaps)]
  
  receiver_gaps = circos_links_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
    sector = rep(width_same_cell_same_receptor_type, times = (circos_links %>% dplyr::filter(receiver == receiver_oi) %>% dplyr::distinct(receptor) %>% nrow() -1))
    gap = width_different_cell
    return(c(sector,gap))
  }) %>% unlist()
  receiver_gaps = receiver_gaps[-length(receiver_gaps)]
  
  gaps = c(sender_gaps, width_ligand_receptor, receiver_gaps, width_ligand_receptor)
  
  if(length(gaps) != length(union(circos_links$ligand, circos_links$receptor) %>% unique())){
    warning("Specified gaps have different length than combined total of ligands and receptors - This is probably due to duplicates in ligand-receptor names")
  }
  
  links_circle$weight[links_circle$weight == 0] = 0.01
  circos.clear()
  circos.par(gap.degree = gaps)
  chordDiagram(links_circle,
               directional = 1,
               order=order,
               link.sort = TRUE,
               link.decreasing = TRUE,
               grid.col = grid_col,
               transparency = transparency,
               diffHeight = 0.0075,
               direction.type = c("diffHeight", "arrows"),
               link.visible = links_circle$weight > 0.01,
               annotationTrack = "grid",
               preAllocateTracks = list(track.height = 0.175),
               link.arr.length = 0.05, link.arr.type = "big.arrow",  link.lwd = 1.25, link.lty = 1,
               reduce = 0,
               scale = TRUE)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1)
  }, bg.border = NA) #
  
  title(title)
  p_circos = recordPlot()
  
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  # grid_col_all = c(colors_receiver, colors_sender)
  legend = ComplexHeatmap::Legend(at = circos_links_oi$receiver %>% unique() %>% sort(),
                                  type = "grid",
                                  legend_gp = grid::gpar(fill = colors_receiver[circos_links_oi$receiver %>% unique() %>% sort()]),
                                  title_position = "topleft",
                                  title = "Receiver")
  ComplexHeatmap::draw(legend, just = c("left", "bottom"))
  
  legend = ComplexHeatmap::Legend(at = circos_links_oi$sender %>% unique() %>% sort(),
                                  type = "grid",
                                  legend_gp = grid::gpar(fill = colors_sender[circos_links_oi$sender %>% unique() %>% sort()]),
                                  title_position = "topleft",
                                  title = "Sender")
  ComplexHeatmap::draw(legend, just = c("left", "top"))
  
  p_legend = grDevices::recordPlot()
  
  circos_plot <- list(circos = p_circos,
                      legend = p_legend)
  
  return(circos_plot)
}


# Plot nichenet res
nicheNetPlot <- function(res, orientation="long"){
  
  if (orientation=="long") {
    
    p1 <- res[["ligandActivityPlot"]]
    p2 <- res[["ligandPearson"]]
    p4 <- res[["ligandTarget"]]
    p3 <- res[["ligandReceptor"]]
    
    
    layout <- "#AAAAAA
    #BCCCCC
    #DDDDDD"
    
    niche_net_plot <- p1 + p2 + p3 + p4 + 
      plot_layout(design = layout)
    
    return(niche_net_plot)

  } else if (orientation=="wide") {  
    
    p1 <- res[["ligandActivityPlot"]]
    p2 <- res[["ligandPearson"]]
    p4 <- res[["ligandTarget"]]
    p3 <- res[["ligandReceptor"]]
    
    
    niche_net_plot <- p1 + p2 + p3 + p4
    
    return(niche_net_plot)
    
  } else {
    
    "Choose wide or long"
    
    
  }  
  
  
}


# Function to plot a dotplot of enrichment results 
plotPathwayDotPlot <- function(enrichmentTable,title){
  enrichmentTable <- result(enrichmentTable) %>% slice_max(abs(NES),n=10) %>%
    mutate(Description = fct_reorder(Description,NES))
  
  ggplot(enrichmentTable, aes(x=NES, y=Description, size=setSize,color=p.adjust)) +
    geom_point() + theme_cowplot(18) +
    scale_size(range=c(3, 8)) +
    guides(size  = guide_legend(order = 1),
           color = guide_colorbar(order = 2)) +
    ggtitle(title)
  
}



# Function to carry out dsea 
calculateDSEA <- function(signatureSearch){
  
  drugs <- result(signatureSearch)%>% group_by(pert) %>%
    slice_max(abs(cor_score))
  drugsList <- drugs$cor_score
  names(drugsList) <- drugs$pert
  drugsList <- sort(drugsList,decreasing = TRUE)
  
  #dsea_GSEA_KEGG <- dsea_GSEA(drugList = drugsList, type="KEGG",exponent=1, 
                              #nPerm=10000, pvalueCutoff=1, minGSSize=2)
  dsea_GSEA_GOMF <- dsea_GSEA(drugList = drugsList, type="GO",ont="MF",exponent=1, 
                              nPerm=10000, pvalueCutoff=1, minGSSize=2)
  #dsea_GSEA_GOBP <- dsea_GSEA(drugList = drugsList, type="KEGG",ont="BP",exponent=1, 
                              #nPerm=10000, pvalueCutoff=1, minGSSize=2)
  
  GOMFPlot <- plotPathwayDotPlot(dsea_GSEA_GOMF,"GO Molecular function")
  #return(list(dsea_GSEA_KEGG, dsea_GSEA_GOMF, dsea_GSEA_GOBP))
  return(list(dsea_GSEA_GOMF, GOMFPlot))
  
}


# Function to create plot of t-score correlation between two comparisons
library(smplot2)

plotCorrelation <- function(x, y, title="", xlab ="", ylab ="") {
  
  DF_ <- merge(x, y, by = 'row.names')
  DF <- DF_ %>% dplyr::select("t.x", "t.y")
  
  cor_plot <- ggplot(DF, aes(x=t.x, y=t.y)) + 
    geom_point() + 
    geom_smooth(method = "lm", alpha = 0.15) +
    sm_statCorr(corr_method = "pearson") + 
    labs(x=xlab, y=ylab) +
    ggtitle(title) +
    theme_cowplot()
  
  return(cor_plot)
  
  
}


# Function to carry out differential cell proportion analysis for a set of contrasts
proportionDEanalysis <- function(contrasts_vector, immune_propslist, design) {
  
  contrast.matrix <- makeContrasts(contrasts=contrasts_vector,
                                   levels=design)
  
  
  outs <- speckle::propeller.ttest(immune_propslist, design, contrast.matrix, robust=TRUE,trend=FALSE, sort=TRUE)
  
  return(outs)
  
  
}
