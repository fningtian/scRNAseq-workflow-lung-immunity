rm(list=ls())

library(DropletUtils) # read 10x
library(hdf5r) # read Read10X_h5
library(Seurat)
library(dplyr)
library(cowplot)
library(patchwork)
library(gridExtra)
library(DoubletDecon)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(tibble) # rownmaes to colnames 


lib = "NS-DD-1s-DEC-1" 
samples <- c("699-unstim-0h", "699-Ig-4h", "699-CD3-4h")

lib = "NS-DD-1s-DEC-2" 
samples <- c("640-unstim-0h", "640-Ig-4h", "640-CD3-4h",
             "678-unstim-0h", "678-Ig-4h", "678-CD3-4h")

lib = "NS-DD-1s-DEC-3" 
samples <- c("548-unstim-0h", "548-Ig-4h", "548-CD3-4h",
             "569-unstim-0h", "569-Ig-4h", "569-CD3-4h")

lib = "NS-DD-1s-DEC-4" 
samples <- c("589-unstim-0h", "589-Ig-4h", "589-CD3-4h",
             "633-unstim-0h", "633-Ig-4h", "633-CD3-4h")

margin <- 1
# margin <- 2

# c("#f7ede2", "#e4c1f9", "#f5cac3", "#a9def9", "#b8f2e6")

check.joint.bcs <- FALSE


sample_colors <- RColorBrewer::brewer.pal(n = length(samples), name = "Pastel2")
color_vector <- c(setNames(sample_colors, samples),
                  "Negative" = "#3a0ca3", "Doublet" = "#005f73")
resDir <- "qc/cellranger"

### ----------------------------------------------------------------------------------------
## I HTODEMUX ---------------------------------------------------------------------
## 0.0 5' demultiplexing setup | 1.0 is dependent on 0.0 ----------------------
  
data_path <- file.path(paste0(resDir, "/", lib, "/count/sample_filtered_feature_bc_matrix.h5"))


data <- Read10X_h5(data_path)
umis_data <- data$`Gene Expression`
hto_data <- data$`Antibody Capture`

if (check.joint.bcs){
  joint.bcs <- intersect(colnames(umis_data), colnames(hto_data))
  umis_data <- umis_data[, joint.bcs]
  hto_data <- as.matrix(hto_data[, joint.bcs])
  print(length(colnames(umis_data)))
  print(length(colnames(hto_data)))
  print(length(joint.bcs))
}



seuratObj_multiplex <- CreateSeuratObject(counts = umis_data)
seuratObj_multiplex <- NormalizeData(seuratObj_multiplex)
seuratObj_multiplex <- FindVariableFeatures(seuratObj_multiplex, selection.method = "mean.var.plot")
seuratObj_multiplex <- ScaleData(seuratObj_multiplex, features = VariableFeatures(seuratObj_multiplex))
seuratObj_multiplex[["HTO"]] <- CreateAssayObject(counts = hto_data)

positive.quantile <- c(0.90, 0.95, 0.99)


# change to defualt clara to better deal with noise and outliers
plot_HTO_tag_counts <- function(seuratObj, positive.quantile, i, margin = margin, 
                                check.joint.bcs, plots = TRUE){
  # seuratObj[["HTO"]] <- CreateAssayObject(counts = data$`Antibody Capture`)
  # rownames(seuratObj[["HTO"]])  #"699-unstim-0h" "699-Ig-4h"     "699-CD3-4h"  
  seuratObj <- NormalizeData(seuratObj, assay = "HTO", normalization.method = "CLR", margin = margin)
  seuratObj <- HTODemux(seuratObj, assay = "HTO", kfunc = "clara", 
                        positive.quantile = positive.quantile[i]) 
  
  
  HTOHeatmap(seuratObj)
  
  # Idents(seuratObj) <- "HTO_maxID"
  # 
  # seuratObj <- SetIdent(seuratObj,
  #                       value = factor(Idents(seuratObj), levels = c("699-unstim-0h",
  #                                                                    "699-Ig-4h",
  #                                                                    "699-CD3-4h")))
  
  
  Idents(seuratObj) <- "hash.ID"
  
  
  seuratObj <- SetIdent(seuratObj,
                        value = factor(Idents(seuratObj), levels = c(samples,
                                                                     "Negative",
                                                                     "Doublet")))
  
  # RidgePlot(seuratObj, assay = "HTO", features = rownames(seuratObj[["HTO"]])[1:3], ncol = 3) + 
  #   scale_fill_manual(values = c("#f7ede2", "#e4c1f9", "#f5cac3", "#a9def9", "#b8f2e6")) +
  #   theme_bw()
  
  
  ridge_plots <- lapply(1:length(samples), function(i) {
    RidgePlot(seuratObj, assay = "HTO", features = rownames(seuratObj[["HTO"]])[i], ncol = 1) +
      scale_fill_manual(values = color_vector) +
      theme_bw() +
      theme(axis.text.y = element_text(size = 8)) +
      guides(fill = "none")
  })
  
  
  # Combine the Ridge plots using patchwork
  RidgePlot <- wrap_plots(ridge_plots, ncol = 3)
  
  VlnPlot1 <- VlnPlot(seuratObj, features = "nFeature_RNA", pt.size = 0, log = TRUE) +
    geom_jitter(aes(y = nFeature_RNA), color = "grey", size = 0.05, alpha = 1, position = position_jitter(width = 0.3)) +  
    scale_fill_manual(values = color_vector) +
    theme_bw() +
    guides(fill = "none")
  
  VlnPlot2 <- VlnPlot(seuratObj, features = "nCount_RNA", pt.size = 0, log = TRUE) +
    geom_jitter(aes(y = nCount_RNA), color = "grey", size = 0.05, alpha = 1, position = position_jitter(width = 0.3)) +  
    scale_fill_manual(values = color_vector) +
    theme_bw() +
    guides(fill = "none")
  
  # VlnPlot2 <- VlnPlot(seuratObj, features = "nFeature_RNA", pt.size = 0.1, log = TRUE) +
  #   scale_fill_manual(values = c("#f7ede2", "#e4c1f9", "#f5cac3", "#a9def9", "#b8f2e6")) +
  #   theme_bw()
  
  # counts of classification
  table(seuratObj@meta.data$hash.ID) 
  
  DefaultAssay(seuratObj) <- "HTO"
  seuratObj <- ScaleData(seuratObj, features = rownames(seuratObj),
                         verbose = FALSE)
  seuratObj <- RunPCA(seuratObj, features = rownames(seuratObj), approx = FALSE)
  seuratObj <- RunTSNE(seuratObj, dims = 1:3, check_duplicates = FALSE, perplexity = 100)
  
  
  seuratObj@meta.data$hash.ID <-  factor(seuratObj@meta.data$hash.ID, levels = c(samples,
                                                                                 "Negative",
                                                                                 "Doublet"))
  
  if (plots==FALSE) {
    return(seuratObj)
  }
  
  
  
  DimPlot1 <- DimPlot(seuratObj, group.by = "hash.ID") +
    scale_color_manual(values = color_vector) +
    theme_bw()
  
  # DimPlot <- (DimPlot1 + DimPlot1) + plot_layout(ncol = 2)
  
  
  hash.ID.df <- as.data.frame(table(seuratObj@meta.data$hash.ID))
  colnames(hash.ID.df) <- c("hash.ID", "count")
  
  RidgePlot_rel_height <- 1 + 1.3 * ceiling((length(samples) - 3) / 3) 
  
  
  output_dir <- paste0("qc/qc_", lib, "/")
  
  if (check.joint.bcs) {
    output_dir <- paste0(output_dir, "htodemux_joint_bcs/")
  }else{
    output_dir <- paste0(output_dir, "htodemux/")
  }
  
  # Create the directory structure if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  
  pdf(paste0(output_dir, positive.quantile[i], "_per_tag_margin_", margin, ".pdf"), width = 7.6, height = 9.76)
  
  p1 <- plot_grid(RidgePlot, VlnPlot1, VlnPlot2,
                  labels = c("A", "B", "C", "D"), 
                  nrow = 3,
                  rel_heights = c(RidgePlot_rel_height,1,1,1))
  
  
  p2 <- plot_grid(DimPlot1, tableGrob(hash.ID.df),
                  labels = c("E", "F"),
                  nrow = 1,
                  rel_widths = c(1.4, 1))
  
  print(plot_grid(p1, p2,
                  nrow = 2,
                  rel_heights = c(2,1)))
  
  dev.off()
  
  # normalize and scale the RNA layer 
  
  DefaultAssay(seuratObj) <- "RNA"
  seuratObj <- RunPCA(object = seuratObj, features = VariableFeatures(object = seuratObj))
  seuratObj <- RunUMAP(seuratObj, reduction = "pca", dims = 1:20)
  seuratObj <- FindNeighbors(seuratObj, reduction = "pca", dims = 1:20)
  seuratObj <- FindClusters(seuratObj, resolution = 0.5)
  seuratObj <- RunTSNE(object = seuratObj, dims = 1:20)
  
  return(seuratObj)
  
}

## 1.0 5' demultiplexing  -----------------------

# plot_HTO_tag_counts(positive.quantile, 1)
# plot_HTO_tag_counts(positive.quantile, 2)
# plot_HTO_tag_counts(positive.quantile, 3)

# subset seurat object based on hash IDs
# the object has everything, not just singlets

# based on the gene expression profile, where does each hashtag go? 


seuratObj_0.99 <- plot_HTO_tag_counts(seuratObj_multiplex, positive.quantile, 3, margin, check.joint.bcs)


DimPlot1_0.99 <- DimPlot(seuratObj_0.99, reduction = "umap", label = T, repel = T) + 
  labs(title = 'UMAP clustering 0.99', x = "UMAP 1", y = 'UMAP 2') +
  theme_bw() +
  theme(legend.position = "none")

DimPlot2_0.99 <- DimPlot(seuratObj_0.99, group.by = "hash.ID") +
  scale_color_manual(values = color_vector) +
  labs(title = 'hash.ID 0.99', x = "UMAP 1", y = 'UMAP 2') +
  theme_bw() 


seuratObj_0.95 <- plot_HTO_tag_counts(seuratObj_multiplex, positive.quantile, 2, margin, check.joint.bcs)

DimPlot1_0.95 <- DimPlot(seuratObj_0.95, reduction = "umap", label = T, repel = T) + 
  labs(title = 'UMAP clustering 0.95', x = "UMAP 1", y = 'UMAP 2') +
  theme_bw() +
  theme(legend.position = "none")

DimPlot2_0.95 <- DimPlot(seuratObj_0.95, group.by = "hash.ID") +
  scale_color_manual(values = color_vector) +
  labs(title = 'hash.ID 0.95', x = "UMAP 1", y = 'UMAP 2') +
  theme_bw()


seuratObj_0.9 <- plot_HTO_tag_counts(seuratObj_multiplex, positive.quantile, 1, margin, check.joint.bcs)


DimPlot1_0.9 <- DimPlot(seuratObj_0.9, reduction = "umap", label = T, repel = T) + 
  labs(title = 'UMAP clustering 0.9', x = "UMAP 1", y = 'UMAP 2') +
  theme_bw() +
  theme(legend.position = "none")

DimPlot2_0.9 <- DimPlot(seuratObj_0.9, group.by = "hash.ID") +
  scale_color_manual(values = color_vector) +
  labs(title = 'hash.ID 0.9', x = "UMAP 1", y = 'UMAP 2') +
  theme_bw()

output_dir <- paste0("qc/qc_", lib, "/")

if (check.joint.bcs) {
  output_dir <- paste0(output_dir, "htodemux_joint_bcs/")
}else{
  output_dir <- paste0(output_dir, "htodemux/")
}

# Create the directory structure if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}


pdf(paste0(output_dir, "gex_0.9_0.95_0.99_per_tag_margin_", margin,  ".pdf"), width = 7.35, height = 8.88)

DimPlot <- plot_grid(DimPlot1_0.9, DimPlot2_0.9, 
                     DimPlot1_0.95, DimPlot2_0.95, 
                     DimPlot1_0.99, DimPlot2_0.99, 
                     ncol = 2, 
                     labels = c("A", "", "B", "", "C", ""),
                     rel_widths = c(0.7, 1, 0.7, 1, 0.7, 1))
print(DimPlot)

dev.off()


# further exploration on demutiplexing


# run till the plot_HTO_tag_counts function
# we are not generating any plots 

# demultiplex_seuratObj <- function(seuratObj, positive.quantile, quantile_index, margin) {
#   # Apply plot_HTO_tag_counts with specified quantile index
#   seuratObj <- plot_HTO_tag_counts(seuratObj_multiplex, positive.quantile, quantile_index, margin = margin, plots = FALSE)
#   
#   # Set default assay to RNA for dimensional reduction and clustering
#   DefaultAssay(seuratObj) <- "RNA"
#   
#   # Run PCA
#   seuratObj <- RunPCA(object = seuratObj, features = VariableFeatures(object = seuratObj))
#   
#   # Run UMAP, neighbors, clustering, and t-SNE
#   seuratObj <- RunUMAP(seuratObj, reduction = "pca", dims = 1:20)
#   seuratObj <- FindNeighbors(seuratObj, reduction = "pca", dims = 1:20)
#   seuratObj <- FindClusters(seuratObj, resolution = 0.5)
#   seuratObj <- RunTSNE(object = seuratObj, dims = 1:20)
#   
#   # Optional: set default assay back to HTO if needed
#   DefaultAssay(seuratObj) <- "HTO"
#   
#   return(seuratObj)
# }

# Process objects for different quantile levels
# seuratObj_0.9 <- demultiplex_seuratObj(seuratObj_multiplex, positive.quantile, 1, margin = margin)
# seuratObj_0.95 <- demultiplex_seuratObj(seuratObj_multiplex, positive.quantile, 2, margin = margin)
# seuratObj_0.99 <- demultiplex_seuratObj(seuratObj_multiplex, positive.quantile, 3, margin = margin)

DefaultAssay(seuratObj_0.9) <- "HTO"
DefaultAssay(seuratObj_0.95) <- "HTO"
DefaultAssay(seuratObj_0.99) <- "HTO"



# Generate FeaturePlots for each quantile
plot_0.9 <- FeaturePlot(seuratObj_0.9, features = rownames(seuratObj_0.9[["HTO"]]), 
                        cols = c("lightgrey", "blue"), ncol = 3,  order = T) +
  plot_annotation(title = paste0(lib, " 0.9"))

plot_0.95 <- FeaturePlot(seuratObj_0.95, features = rownames(seuratObj_0.95[["HTO"]]), 
                         cols = c("lightgrey", "blue"), ncol = 3,  order = T) +
  plot_annotation(title = paste0(lib, " 0.95"))

plot_0.99 <- FeaturePlot(seuratObj_0.99, features = rownames(seuratObj_0.99[["HTO"]]), 
                         cols = c("lightgrey", "blue"), ncol = 3,  order = T) +
  plot_annotation(title = paste0(lib, " 0.99"))


if(length(samples) == 3){
  plot_height = 8
}else{
  plot_height = 16
}


if (check.joint.bcs) {
  output_dir <- paste0(output_dir, "htodemux_joint_bcs/")
}else{
  output_dir <- paste0(output_dir, "htodemux/")
}

# Create the directory structure if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}


pdf(paste0(output_dir, "hto_exp_0.9_0.95_0.99_per_tag_margin", margin, ".pdf"), width = 8.94, height = plot_height)

plot_grid(plot_0.9, plot_0.95, plot_0.99,
          ncol = 1)

dev.off()



plot_HTO_expression_pairs <- function(seuratObj, positive.quantile, quantile_index, margin,
                                      check.joint.bcs) {
  hashtags <- rownames(seuratObj[["HTO"]])
  plots_list <- list()
  
  # Loop over all pairs of hashtags
  for (i in 1:(length(hashtags) - 1)) {
    for (j in (i + 1):length(hashtags)) {
      
      # Extract the expression data for the two hashtags
      hashtag1 <- hashtags[i]
      hashtag2 <- hashtags[j]
      
      exp1 <- seuratObj[["HTO"]]@data[hashtag1, ]
      exp2 <- seuratObj[["HTO"]]@data[hashtag2, ]
      
      # Create a data frame for ggplot
      plot_data <- data.frame(
        Exp1 = exp1,
        Exp2 = exp2,
        HashID = seuratObj$hash.ID  # Metadata for coloring
      )
      
      # Create the scatter plot for this pair of hashtags
      p <- ggplot(plot_data, aes(x = Exp1, y = Exp2, color = HashID)) +
        geom_point(alpha = 1) +
        labs(
          x = paste(hashtag1, "exp"),
          y = paste(hashtag2, "exp"),
          title = paste(hashtag1, "vs", hashtag2)
        ) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # Add 1:1 line
        theme_bw() +
        theme(plot.title = element_text(size = 10)) +
        scale_color_manual(values = color_vector, name = paste("p.q.", positive.quantile[quantile_index]))  # Updated to use positive.quantile[3]
      
      # Add the plot to the list
      plots_list[[paste(hashtag1, hashtag2, sep = "_vs_")]] <- p
    }
  }
  
  plot_height <- ifelse(length(hashtags) == 3, 2.6, 13)
  
  # Combine all the plots into one figure and only display the legend once
  final_plot <- wrap_plots(plots_list, ncol = 3) +
    plot_layout(guides = "collect")
  
  if (check.joint.bcs) {
    output_dir <- paste0(output_dir, "htodemux_joint_bcs/")
  }else{
    output_dir <- paste0(output_dir, "htodemux/")
  }
  
  # Create the directory structure if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  
  pdf(paste0(output_dir, "hto_exp_pairwise_", positive.quantile[quantile_index], "_margin_",  margin, ".pdf"), width = 9, height = plot_height)
  print(final_plot)
  dev.off()
}


plot_HTO_expression_pairs(
  seuratObj = seuratObj_0.99,
  positive.quantile = positive.quantile,
  quantile_index = 3,
  margin, check.joint.bcs)


plot_HTO_expression_pairs(
  seuratObj = seuratObj_0.95,
  positive.quantile = positive.quantile,
  quantile_index = 2,
  margin, check.joint.bcs)

plot_HTO_expression_pairs(
  seuratObj = seuratObj_0.9,
  positive.quantile = positive.quantile,
  quantile_index = 1,
  margin, check.joint.bcs)














### ----------------------------------------------------------------------------------------
### ----------------------------------------------------------------------------------------
### ----------------------------------------------------------------------------------------
### ----------------------------------------------------------------------------------------
### ----------------------------------------------------------------------------------------
## II MULTIseqDemux ----------------------------------------------------------------
## 0.0 5' demultiplexing setup ----------------------

resDir <- "qc/cellranger"

data_path <- file.path(paste0(resDir, "/", lib, "/count/sample_filtered_feature_bc_matrix.h5"))


data <- Read10X_h5(data_path)
umis_data <- data$`Gene Expression`
hto_data <- data$`Antibody Capture`

if (check.joint.bcs){
  joint.bcs <- intersect(colnames(umis_data), colnames(hto_data))
  umis_data <- umis_data[, joint.bcs]
  hto_data <- as.matrix(hto_data[, joint.bcs])
  print(length(colnames(umis_data)))
  print(length(colnames(hto_data)))
  print(length(joint.bcs))
}



seuratObj_multiplex <- CreateSeuratObject(counts = umis_data)
seuratObj_multiplex <- NormalizeData(seuratObj_multiplex)
seuratObj_multiplex <- FindVariableFeatures(seuratObj_multiplex, selection.method = "mean.var.plot")
seuratObj_multiplex <- ScaleData(seuratObj_multiplex, features = VariableFeatures(seuratObj_multiplex))
seuratObj_multiplex[["HTO"]] <- CreateAssayObject(counts = hto_data)


# change to defualt clara to better deal with noise and outliers
plot_HTO_tag_counts <- function(seuratObj, margin = margin, 
                                check.joint.bcs, plots = TRUE){
  # seuratObj[["HTO"]] <- CreateAssayObject(counts = data$`Antibody Capture`)
  # rownames(seuratObj[["HTO"]])  #"699-unstim-0h" "699-Ig-4h"     "699-CD3-4h"  
  seuratObj <- NormalizeData(seuratObj, assay = "HTO", normalization.method = "CLR", margin = margin)
  seuratObj <- MULTIseqDemux(seuratObj,
                      assay = "HTO",
                      autoThresh = TRUE,
                      maxiter = 10,
                      qrange = seq(from = 0.1, to = 0.9, by = 0.05),
                      verbose = TRUE)
  
  
  # HTOHeatmap(seuratObj)
  
  # Idents(seuratObj) <- "HTO_maxID"
  # 
  # seuratObj <- SetIdent(seuratObj,
  #                       value = factor(Idents(seuratObj), levels = c("699-unstim-0h",
  #                                                                    "699-Ig-4h",
  #                                                                    "699-CD3-4h")))
  
  
  Idents(seuratObj) <- "MULTI_ID"
  
  
  seuratObj <- SetIdent(seuratObj,
                        value = factor(Idents(seuratObj), levels = c(samples,
                                                                     "Negative",
                                                                     "Doublet")))
  
  # RidgePlot(seuratObj, assay = "HTO", features = rownames(seuratObj[["HTO"]])[1:3], ncol = 3) + 
  #   scale_fill_manual(values = c("#f7ede2", "#e4c1f9", "#f5cac3", "#a9def9", "#b8f2e6")) +
  #   theme_bw()
  
  
  ridge_plots <- lapply(1:length(samples), function(i) {
    RidgePlot(seuratObj, assay = "HTO", features = rownames(seuratObj[["HTO"]])[i], ncol = 1) +
      scale_fill_manual(values = color_vector) +
      theme_bw() +
      theme(axis.text.y = element_text(size = 8)) +
      guides(fill = "none")
  })
  
  
  # Combine the Ridge plots using patchwork
  RidgePlot <- wrap_plots(ridge_plots, ncol = 3)
  
  VlnPlot1 <- VlnPlot(seuratObj, features = "nFeature_RNA", pt.size = 0, log = TRUE) +
    geom_jitter(aes(y = nFeature_RNA), color = "grey", size = 0.05, alpha = 1, position = position_jitter(width = 0.3)) +  
    scale_fill_manual(values = color_vector) +
    theme_bw() +
    guides(fill = "none")
  
  VlnPlot2 <- VlnPlot(seuratObj, features = "nCount_RNA", pt.size = 0, log = TRUE) +
    geom_jitter(aes(y = nCount_RNA), color = "grey", size = 0.05, alpha = 1, position = position_jitter(width = 0.3)) +  
    scale_fill_manual(values = color_vector) +
    theme_bw() +
    guides(fill = "none")
  
  # VlnPlot2 <- VlnPlot(seuratObj, features = "nFeature_RNA", pt.size = 0.1, log = TRUE) +
  #   scale_fill_manual(values = c("#f7ede2", "#e4c1f9", "#f5cac3", "#a9def9", "#b8f2e6")) +
  #   theme_bw()
  
  # counts of classification
  table(seuratObj@meta.data$MULTI_ID) 
  
  DefaultAssay(seuratObj) <- "HTO"
  seuratObj <- ScaleData(seuratObj, features = rownames(seuratObj),
                         verbose = FALSE)
  seuratObj <- RunPCA(seuratObj, features = rownames(seuratObj), approx = FALSE)
  seuratObj <- RunTSNE(seuratObj, dims = 1:3, check_duplicates = FALSE, perplexity = 100)
  
  
  seuratObj@meta.data$MULTI_ID <-  factor(seuratObj@meta.data$MULTI_ID, levels = c(samples,
                                                                                 "Negative",
                                                                                 "Doublet"))
  
  if (plots==FALSE) {
    return(seuratObj)
  }
  
  
  
  DimPlot1 <- DimPlot(seuratObj, group.by = "MULTI_ID") +
    scale_color_manual(values = color_vector) +
    theme_bw()
  
  # DimPlot <- (DimPlot1 + DimPlot1) + plot_layout(ncol = 2)
  
  
  MULTI_ID.df <- as.data.frame(table(seuratObj@meta.data$MULTI_ID))
  colnames(MULTI_ID.df) <- c("MULTI_ID", "count")
  
  RidgePlot_rel_height <- 1 + 1.3 * ceiling((length(samples) - 3) / 3) 
  
  
  output_dir <- paste0("qc/qc_", lib, "/")
  
  if (check.joint.bcs) {
    output_dir <- paste0(output_dir, "multiseqdemux_joint_bcs/")
  }else{
    output_dir <- paste0(output_dir, "multiseqdemux/")
  }
  
  # Create the directory structure if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  
  pdf(paste0(output_dir, "optimal_pq_per_tag_margin_", margin, ".pdf"), width = 7.6, height = 9.76)
  
  p1 <- plot_grid(RidgePlot, VlnPlot1, VlnPlot2,
                  labels = c("A", "B", "C", "D"), 
                  nrow = 3,
                  rel_heights = c(RidgePlot_rel_height,1,1,1))
  
  
  p2 <- plot_grid(DimPlot1, tableGrob(MULTI_ID.df),
                  labels = c("E", "F"),
                  nrow = 1,
                  rel_widths = c(1.4, 1))
  
  print(plot_grid(p1, p2,
                  nrow = 2,
                  rel_heights = c(2,1)))
  
  dev.off()
  
  # normalize and scale the RNA layer 
  
  DefaultAssay(seuratObj) <- "RNA"
  seuratObj <- RunPCA(object = seuratObj, features = VariableFeatures(object = seuratObj))
  seuratObj <- RunUMAP(seuratObj, reduction = "pca", dims = 1:20)
  seuratObj <- FindNeighbors(seuratObj, reduction = "pca", dims = 1:20)
  seuratObj <- FindClusters(seuratObj, resolution = 0.5)
  seuratObj <- RunTSNE(object = seuratObj, dims = 1:20)
  
  return(seuratObj)
  
}

## run 0.0 frist -> 1.0 5' demultiplexing  -----------------------

# plot_HTO_tag_counts(positive.quantile, 1)
# plot_HTO_tag_counts(positive.quantile, 2)
# plot_HTO_tag_counts(positive.quantile, 3)

# subset seurat object based on hash IDs
# the object has everything, not just singlets

# based on the gene expression profile, where does each hashtag go? 


seuratObj_demultiplex <- plot_HTO_tag_counts(seuratObj_multiplex, margin, check.joint.bcs)

DimPlot1_demultiplex <- DimPlot(seuratObj_demultiplex, reduction = "umap", label = T, repel = T) + 
  labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2') +
  theme_bw() +
  theme(legend.position = "none")

DimPlot2_demultiplex <- DimPlot(seuratObj_demultiplex, group.by = "MULTI_ID") +
  scale_color_manual(values = color_vector) +
  labs(title = 'MULTI_ID', x = "UMAP 1", y = 'UMAP 2') +
  theme_bw()

output_dir <- paste0("qc/qc_", lib, "/")

if (check.joint.bcs) {
  output_dir <- paste0(output_dir, "multiseqdemux_joint_bcs/")
}else{
  output_dir <- paste0(output_dir, "multiseqdemux/")
}

# # Create the directory structure if it doesn't exist
# if (!dir.exists(output_dir)) {
#   dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
# }
# 

pdf(paste0(output_dir, "gex_per_tag_margin_", margin,  ".pdf"), width = 7.35, height = 8.88/3)

DimPlot <- plot_grid(DimPlot1_demultiplex, DimPlot2_demultiplex, 
                     ncol = 2, 
                     labels = c("A", ""),
                     rel_widths = c(0.7, 1))
print(DimPlot)

dev.off()


# further exploration on demutiplexing
DefaultAssay(seuratObj_demultiplex) <- "HTO"



# Generate FeaturePlots for each quantile
plot_demultiplex <- FeaturePlot(seuratObj_demultiplex, features = rownames(seuratObj_demultiplex[["HTO"]]), 
                        cols = c("lightgrey", "blue"), ncol = 3,  order = T) +
  plot_annotation(title = paste0(lib, " 0.9"))

if(length(samples) == 3){
  plot_height = 8/3
}else{
  plot_height = 16/3
}

pdf(paste0(output_dir, "hto_exp_per_tag_margin", margin, ".pdf"), width = 8.94, height = plot_height)

plot_grid(plot_demultiplex,
          ncol = 1)

dev.off()



plot_HTO_expression_pairs <- function(seuratObj, margin,
                                      check.joint.bcs) {
  hashtags <- rownames(seuratObj[["HTO"]])
  plots_list <- list()
  
  # Loop over all pairs of hashtags
  for (i in 1:(length(hashtags) - 1)) {
    for (j in (i + 1):length(hashtags)) {
      
      # Extract the expression data for the two hashtags
      hashtag1 <- hashtags[i]
      hashtag2 <- hashtags[j]
      
      exp1 <- seuratObj[["HTO"]]@data[hashtag1, ]
      exp2 <- seuratObj[["HTO"]]@data[hashtag2, ]
      
      # Create a data frame for ggplot
      plot_data <- data.frame(
        Exp1 = exp1,
        Exp2 = exp2,
        MULTIID = seuratObj$MULTI_ID  # Metadata for coloring
      )
      
      # Create the scatter plot for this pair of hashtags
      p <- ggplot(plot_data, aes(x = Exp1, y = Exp2, color = MULTIID)) +
        geom_point(alpha = 1) +
        labs(
          x = paste(hashtag1, "exp"),
          y = paste(hashtag2, "exp"),
          title = paste(hashtag1, "vs", hashtag2)
        ) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # Add 1:1 line
        theme_bw() +
        theme(plot.title = element_text(size = 10)) +
        scale_color_manual(values = color_vector, name = paste("optimal p.q."))  # Updated to use positive.quantile[3]
      
      # Add the plot to the list
      plots_list[[paste(hashtag1, hashtag2, sep = "_vs_")]] <- p
    }
  }
  
  plot_height <- ifelse(length(hashtags) == 3, 2.6, 13)
  
  # Combine all the plots into one figure and only display the legend once
  final_plot <- wrap_plots(plots_list, ncol = 3) +
    plot_layout(guides = "collect")
  
  output_dir <- paste0("qc/qc_", lib, "/")
  
  if (check.joint.bcs) {
    output_dir <- paste0(output_dir, "multiseqdemux_joint_bcs/")
  }else{
    output_dir <- paste0(output_dir, "multiseqdemux/")
  }
  
  # Create the directory structure if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  
  pdf(paste0(output_dir, "hto_exp_pairwise_optimal_pq_margin_",  margin, ".pdf"), width = 9, height = plot_height)
  print(final_plot)
  dev.off()
}


plot_HTO_expression_pairs(
  seuratObj = seuratObj_demultiplex,
  margin, check.joint.bcs)







## run 0.0 frist -> 1.0.1 look into NS-DD-1s-DEC-2 ----

doublets <- seuratObj_demultiplex@meta.data %>%
  filter(MULTI_ID == "Doublet")

unique(doublets$MULTI_classification)

# remove leading and trailing spaces 
doublets$MULTI_classification <- trimws(doublets$MULTI_classification)

doublets$bar_color <- ifelse(grepl("678-CD3-4h", doublets$MULTI_classification), "red", "skyblue")

unique(doublets$bar_color)

p <- ggplot(doublets, aes(x = MULTI_classification, fill = bar_color)) +
  geom_bar(color = "black") +  # Bar plot with border color
  labs(
    title = "Distribution of MULTI_classification for Doublets",
    x = "MULTI_classification",
    y = "Count"
  ) +
  scale_fill_identity() +  # Use the colors defined in the 'bar_color' column
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x labels
    strip.text = element_text(size = 12, face = "bold") # Adjust strip text size
  )

pdf(paste0("qc/qc_NS-DD-1s-DEC-2/multiseqdemux/doublets.pdf"), width = 7.25, height = 5.18)
print(p)
dev.off()



# doublets exp of 678-CD3-4h - all else
# seuratObj_demultiplex@assays$HTO@data
exp_diff <- t(seuratObj_demultiplex@assays$HTO@data) %>%
  as.data.frame() %>%
  mutate(
    diff = `678-CD3-4h` - (rowSums(.[, -which(names(.) == "678-CD3-4h")])) # Subtract sum of other antibodies
  ) %>%
  rownames_to_column("barcode") 

doublets_exp <- doublets %>%
  rownames_to_column("barcode") %>% 
  left_join(exp_diff, by = "barcode")

any(is.na(doublets_exp$diff))

doublets_exp_long <- doublets_exp %>%
  filter(grepl("678-CD3-4h", MULTI_classification)) %>%
  select(colnames(exp_diff)) %>%
  pivot_longer(cols = colnames(exp_diff)[-1], names_to = "antibody", values_to = "norm.exp")

p <- ggplot(doublets_exp_long, aes(x = norm.exp, color = antibody)) +
  geom_density(linewidth = 1.2) +
  scale_color_manual(values = c(color_vector[1:6], "diff" = "darkred")) +
  theme_bw() +
  labs(
    x = "Antibody",
    y = "Normalized expression"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(paste0("qc/qc_NS-DD-1s-DEC-2/multiseqdemux/doublets_subtraction.pdf"), width = 4.5, height = 3)
print(p)
dev.off()

seuratObj_demultiplex@meta.data$MULTI_ID_PREV <- seuratObj_demultiplex@meta.data$MULTI_ID


# adjust parameters below
antibody_barcodes <- doublets_exp %>%
  filter(grepl("678-CD3-4h", MULTI_classification))

# Apply quantile threshold only to filtered barcodes
quantiles <- quantile(antibody_barcodes$diff, c(0.025, 0.975))

rescued_cells <- antibody_barcodes %>%
  filter(diff >= quantiles[1] & diff <= quantiles[2]) %>%
  pull(barcode)

seuratObj_demultiplex@meta.data$MULTI_ID <- seuratObj_demultiplex@meta.data$MULTI_ID_PREV

seuratObj_demultiplex@meta.data[rescued_cells, "MULTI_ID"] <- "678-CD3-4h"


Idents(seuratObj_demultiplex) <- "MULTI_ID"


seuratObj_demultiplex <- SetIdent(seuratObj_demultiplex,
                      value = factor(Idents(seuratObj_demultiplex), levels = c(samples,
                                                                   "Negative",
                                                                   "Doublet")))


ridge_plots <- lapply(1:length(samples), function(i) {
  RidgePlot(seuratObj_demultiplex, assay = "HTO", features = rownames(seuratObj_demultiplex[["HTO"]])[i], ncol = 1) +
    scale_fill_manual(values = color_vector) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 8)) +
    guides(fill = "none")
})


# Combine the Ridge plots using patchwork
RidgePlot <- wrap_plots(ridge_plots, ncol = 3)

VlnPlot1 <- VlnPlot(seuratObj_demultiplex, features = "nFeature_RNA", pt.size = 0, log = TRUE) +
  geom_jitter(aes(y = nFeature_RNA), color = "grey", size = 0.05, alpha = 1, position = position_jitter(width = 0.3)) +  
  scale_fill_manual(values = color_vector) +
  theme_bw() +
  guides(fill = "none")

VlnPlot2 <- VlnPlot(seuratObj_demultiplex, features = "nCount_RNA", pt.size = 0, log = TRUE) +
  geom_jitter(aes(y = nCount_RNA), color = "grey", size = 0.05, alpha = 1, position = position_jitter(width = 0.3)) +  
  scale_fill_manual(values = color_vector) +
  theme_bw() +
  guides(fill = "none")


# counts of classification
table(seuratObj_demultiplex@meta.data$MULTI_ID) 

DefaultAssay(seuratObj_demultiplex) <- "HTO"
seuratObj_demultiplex <- ScaleData(seuratObj_demultiplex, features = rownames(seuratObj_demultiplex),
                       verbose = FALSE)
seuratObj_demultiplex <- RunPCA(seuratObj_demultiplex, features = rownames(seuratObj_demultiplex), approx = FALSE)
seuratObj_demultiplex <- RunTSNE(seuratObj_demultiplex, dims = 1:3, check_duplicates = FALSE, perplexity = 100)


seuratObj_demultiplex@meta.data$MULTI_ID <-  factor(seuratObj_demultiplex@meta.data$MULTI_ID, levels = c(samples,
                                                                                 "Negative",
                                                                                 "Doublet"))
DimPlot1 <- DimPlot(seuratObj_demultiplex, group.by = "MULTI_ID") +
  scale_color_manual(values = color_vector) +
  theme_bw()

# DimPlot <- (DimPlot1 + DimPlot1) + plot_layout(ncol = 2)


MULTI_ID.df <- as.data.frame(table(seuratObj_demultiplex@meta.data$MULTI_ID))
colnames(MULTI_ID.df) <- c("MULTI_ID", "count")

RidgePlot_rel_height <- 1 + 1.3 * ceiling((length(samples) - 3) / 3) 


output_dir <- paste0("qc/qc_", lib, "/")

pdf(paste0(output_dir, "multiseqdemux/adjusted_optimal_pq_per_tag_margin_", margin, ".pdf"), width = 7.6, height = 9.76)

p1 <- plot_grid(RidgePlot, VlnPlot1, VlnPlot2,
                labels = c("A", "B", "C", "D"), 
                nrow = 3,
                rel_heights = c(RidgePlot_rel_height,1,1,1))


p2 <- plot_grid(DimPlot1, tableGrob(MULTI_ID.df),
                labels = c("E", "F"),
                nrow = 1,
                rel_widths = c(1.4, 1))

print(plot_grid(p1, p2,
                nrow = 2,
                rel_heights = c(2,1)))

dev.off()


## run 0.0 frist -> 1.0.2 look into NS-DD-1s-DEC-2 after souporcell ----
seuratObj_demultiplex <- plot_HTO_tag_counts(seuratObj_multiplex, margin, check.joint.bcs, plots=FALSE)

saveRDS(seuratObj_demultiplex, "qc/qc_NS-DD-1s-DEC-2/multiseqdemux/seuratObj_demultiplex.rds")
seuratObj_demultiplex <- readRDS("qc/qc_NS-DD-1s-DEC-2/multiseqdemux/seuratObj_demultiplex.rds")


souporcell_path <- "qc/qc_NS-DD-1s-DEC-2/souporcell_count/NS-DD-1s-DEC-2"
souporcell_clts <- read.table(paste0(souporcell_path, "/clusters.tsv"), header = TRUE, sep = "\t")

metadata <- seuratObj_demultiplex@meta.data
metadata$barcode <- rownames(metadata) 
metadata <- merge(metadata, souporcell_clts, by = "barcode", all = TRUE)  
all(rownames(seuratObj_demultiplex@meta.data) == metadata$barcode)

rownames(metadata) <- metadata$barcode
seuratObj_demultiplex@meta.data <- metadata %>% select(-barcode)


# 1.0.2.1 next, ask the question, when (MULTI_ID is not Negative or Doublet) & (status is singlet),
# stacked barchart of assignemnt per MULTI_ID


singlets <- seuratObj_demultiplex@meta.data %>%
  filter(MULTI_ID != "Negative", MULTI_ID != "Doublet", status == "singlet")

p <- ggplot(singlets, aes(x = MULTI_ID, fill = assignment)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("#E3F2FD", "#DB5461")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  labs(title = "% Souporcell assignment \nper MULTI_ID (singlets only)", 
       x = "MULTI_ID", 
       y = "Percentage",
       fill = "Assignment")


pdf(paste0("qc/qc_NS-DD-1s-DEC-2/multiseqdemux/souporcell_assignment_per_MULTI_ID.pdf"), width = 4, height = 4)
print(p)
dev.off()



# 1.0.2.2 focus on doublets in MULTI_ID 
# for Doublet that contains "678-CD3-4h" in MULTI_ID: if they are singlet in status: demultiplex by assignment 
# store the orginal demultiplex results 

seuratObj_demultiplex@meta.data$MULTI_ID_original <- seuratObj_demultiplex@meta.data$MULTI_ID

seuratObj_demultiplex@meta.data <- seuratObj_demultiplex@meta.data %>%
  mutate(barcode = rownames(seuratObj_demultiplex@meta.data)) %>%
  rowwise() %>%
  mutate(
    MULTI_classification_split = strsplit(as.character(MULTI_classification), "_"),
    MULTI_ID = if_else(
      MULTI_ID == "Doublet" & 
        grepl("678-CD3-4h", MULTI_classification) & 
        grepl("640", MULTI_classification) & # change 640 to 640-unstim-0h
        status == "singlet",
      if_else(
        assignment == "1",
        first(grep("^640", MULTI_classification_split, value = TRUE)),
        first(grep("^678", MULTI_classification_split, value = TRUE))
      ),
      MULTI_ID
    )
  ) %>%
  ungroup() %>% # Removes row-wise grouping
  select(-MULTI_classification_split) %>% # Remove the temporary column
  column_to_rownames("barcode") 



underrep_sample1 <-seuratObj_demultiplex@meta.data %>%
  filter(MULTI_ID_original == "Doublet", 
         grepl("678-CD3-4h", MULTI_classification)) %>%
  select(MULTI_ID_original, MULTI_classification, status, assignment) %>%
  mutate(assignment = case_when(
    assignment == "1" ~ "640",
    assignment == "0" ~ "678",
    assignment == "0/1" ~ "Doublet",
    assignment == "1/0" ~ "Doublet",
  ))

max_cnt <- underrep_sample1 %>%
  group_by(MULTI_classification) %>%
  summarize(count = n(), .groups = "drop") %>%
  summarize(max_count = max(count)) %>%
  pull(max_count)

p1 <- underrep_sample1 %>%
  group_by(MULTI_classification, assignment) %>%
  summarize(count = n(), .groups = "drop") %>%
  ggplot(aes(x = MULTI_classification, y = count, fill = assignment)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#ef476f", "#ffd166", "#118ab2")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "678-CD3-4h doublets\nMulti-seq",
       x = "Classification",
       y = "Count",
       fill = "Souporcell\nassignment") 

underrep_sample2 <-seuratObj_demultiplex@meta.data %>%
  filter(MULTI_ID == "Doublet", 
         grepl("678-CD3-4h", MULTI_classification)) %>%
  select(MULTI_ID, MULTI_classification, status, assignment) %>%
  mutate(assignment = case_when(
    assignment == "1" ~ "640",
    assignment == "0" ~ "678",
    assignment == "0/1" ~ "Doublet",
    assignment == "1/0" ~ "Doublet",
  ))

p2 <- underrep_sample2 %>%
  group_by(MULTI_classification, assignment) %>%
  summarize(count = n(), .groups = "drop") %>%
  ggplot(aes(x = MULTI_classification, y = count, fill = assignment)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#ef476f", "#ffd166", "#118ab2")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "678-CD3-4h doublets\nMulti-seq & Souporcell",
       x = "Classification",
       y = "Count",
       fill = "Souporcell\nassignment") +
    ylim(0, max_cnt)
  



pdf(paste0("qc/qc_NS-DD-1s-DEC-2/multiseqdemux/678-CD3-4h_doublets_souporcell_multi-seq.pdf"), width = 6.21, height = 4)
p1 + p2 +
  plot_layout(guides = "collect") & # Shared legend
  theme(legend.position = "bottom", axis.title.x = element_blank())
dev.off()


Idents(seuratObj_demultiplex) <- "MULTI_ID"


seuratObj_demultiplex <- SetIdent(seuratObj_demultiplex,
                                  value = factor(Idents(seuratObj_demultiplex), levels = c(samples,
                                                                                           "Negative",
                                                                                           "Doublet")))

ridge_plots <- lapply(1:length(samples), function(i) {
  RidgePlot(seuratObj_demultiplex, assay = "HTO", features = rownames(seuratObj_demultiplex[["HTO"]])[i], ncol = 1) +
    scale_fill_manual(values = color_vector) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 8)) +
    guides(fill = "none")
})


# Combine the Ridge plots using patchwork
RidgePlot <- wrap_plots(ridge_plots, ncol = 3)

VlnPlot1 <- VlnPlot(seuratObj_demultiplex, features = "nFeature_RNA", pt.size = 0, log = TRUE) +
  geom_jitter(aes(y = nFeature_RNA), color = "grey", size = 0.05, alpha = 1, position = position_jitter(width = 0.3)) +  
  scale_fill_manual(values = color_vector) +
  theme_bw() +
  guides(fill = "none")

VlnPlot2 <- VlnPlot(seuratObj_demultiplex, features = "nCount_RNA", pt.size = 0, log = TRUE) +
  geom_jitter(aes(y = nCount_RNA), color = "grey", size = 0.05, alpha = 1, position = position_jitter(width = 0.3)) +  
  scale_fill_manual(values = color_vector) +
  theme_bw() +
  guides(fill = "none")


# counts of classification
table(seuratObj_demultiplex@meta.data$MULTI_ID) 

DefaultAssay(seuratObj_demultiplex) <- "HTO"
seuratObj_demultiplex <- ScaleData(seuratObj_demultiplex, features = rownames(seuratObj_demultiplex),
                                   verbose = FALSE)
seuratObj_demultiplex <- RunPCA(seuratObj_demultiplex, features = rownames(seuratObj_demultiplex), approx = FALSE)
seuratObj_demultiplex <- RunTSNE(seuratObj_demultiplex, dims = 1:3, check_duplicates = FALSE, perplexity = 100)


seuratObj_demultiplex@meta.data$MULTI_ID <-  factor(seuratObj_demultiplex@meta.data$MULTI_ID, levels = c(samples,
                                                                                                         "Negative",
                                                                                                         "Doublet"))
DimPlot1 <- DimPlot(seuratObj_demultiplex, group.by = "MULTI_ID") +
  scale_color_manual(values = color_vector) +
  theme_bw()

# DimPlot <- (DimPlot1 + DimPlot1) + plot_layout(ncol = 2)


MULTI_ID.df <- as.data.frame(table(seuratObj_demultiplex@meta.data$MULTI_ID))
colnames(MULTI_ID.df) <- c("MULTI_ID", "count")

RidgePlot_rel_height <- 1 + 1.3 * ceiling((length(samples) - 3) / 3) 


output_dir <- paste0("qc/qc_", lib, "/")


pdf(paste0(output_dir, "multiseqdemux/souporcell_adjusted_optimal_pq_per_tag_margin_", margin, ".pdf"), width = 7.6, height = 9.76)
# pdf(paste0(output_dir, "multiseqdemux/640-unstim-0h_678-CD3-4h_souporcell_adjusted_optimal_pq_per_tag_margin_", margin, ".pdf"), width = 7.6, height = 9.76)


p1 <- plot_grid(RidgePlot, VlnPlot1, VlnPlot2,
                labels = c("A", "B", "C", "D"), 
                nrow = 3,
                rel_heights = c(RidgePlot_rel_height,1,1,1))


p2 <- plot_grid(DimPlot1, tableGrob(MULTI_ID.df),
                labels = c("E", "F"),
                nrow = 1,
                rel_widths = c(1.4, 1))

print(plot_grid(p1, p2,
                nrow = 2,
                rel_heights = c(2,1)))

dev.off()




## 1.0.3 look into NS-DD-1s-DEC-2 after souporcell + manual curation ----
# seuratObj_demultiplex <- plot_HTO_tag_counts(seuratObj_multiplex, margin, check.joint.bcs, plots=FALSE)
# 
# saveRDS(seuratObj_demultiplex, "qc/qc_NS-DD-1s-DEC-2/multiseqdemux/seuratObj_demultiplex.rds")
seuratObj_demultiplex <- readRDS("qc/qc_NS-DD-1s-DEC-2/multiseqdemux/seuratObj_demultiplex.rds")


souporcell_path <- "qc/qc_NS-DD-1s-DEC-2/souporcell_count/NS-DD-1s-DEC-2"
souporcell_clts <- read.table(paste0(souporcell_path, "/clusters.tsv"), header = TRUE, sep = "\t")

metadata <- seuratObj_demultiplex@meta.data
metadata$barcode <- rownames(metadata) 
metadata <- merge(metadata, souporcell_clts, by = "barcode", all = TRUE)  
all(rownames(seuratObj_demultiplex@meta.data) == metadata$barcode)

rownames(metadata) <- metadata$barcode
seuratObj_demultiplex@meta.data <- metadata %>% select(-barcode)

# step 1: min exp of doublets in 678-CD3-4h that can possibly be consider as 678-CD3-4h
selected_cells <- rownames(seuratObj_demultiplex@meta.data[
  seuratObj_demultiplex@meta.data$MULTI_ID == "678-CD3-4h", 
])

min_exp <- min(seuratObj_demultiplex[["HTO"]]@data["678-CD3-4h", selected_cells], na.rm = TRUE)

# step 2, apply filters on current doublets (n = 3274)
doublets <- rownames(seuratObj_demultiplex@meta.data[
  seuratObj_demultiplex@meta.data$MULTI_ID == "Doublet", 
])

# length(rownames(seuratObj_demultiplex@meta.data[
#   grep("678-CD3-4h", seuratObj_demultiplex@meta.data$MULTI_classification),
# ]))
# 
# df <- as.data.frame(table(seuratObj_demultiplex@meta.data$MULTI_classification))
# 
# length(rownames(seuratObj_demultiplex@meta.data[
#   is.na(seuratObj_demultiplex@meta.data$MULTI_classification), 
# ]))


hto_data <- seuratObj_demultiplex[["HTO"]]@data

# Step 2: meet the conditions (n = 2279)
filtered_doublets <- doublets[
  hto_data["678-CD3-4h", doublets] > min_exp & # Above min_exp
    hto_data["640-Ig-4h", doublets] < 2.5 &     # <2.5 in 640-Ig-4h
    hto_data["640-CD3-4h", doublets] < 2 &     # <2 in 640-CD3-4h
    hto_data["678-Ig-4h", doublets] < 2 &      # <2 in 678-Ig-4h
    hto_data["678-unstim-0h", doublets] < 2    # <2 in 678-unstim-0h
]

excluded_doublets <- seuratObj_demultiplex@meta.data %>%
  filter(MULTI_ID == "Doublet", 
         grepl("678-CD3-4h", MULTI_classification),
         grepl("640-unstim-0h", MULTI_classification)) %>%
  mutate(assignment = case_when(
    assignment == "1" ~ "640",
    assignment == "0" ~ "678",
    assignment == "0/1" ~ "Doublet",
    assignment == "1/0" ~ "Doublet",
  )) %>%
  filter(assignment == "640")



doublets1 <- doublets[
  hto_data["678-CD3-4h", doublets] >= min_exp &
    hto_data["640-Ig-4h", doublets] >= 2.5
]

doublets2 <- doublets[
  hto_data["678-CD3-4h", doublets] >= min_exp &
    hto_data["640-CD3-4h", doublets] >= 2
]

doublets3 <- doublets[
  hto_data["678-CD3-4h", doublets] >= min_exp &
    hto_data["678-Ig-4h", doublets] >= 2
]

doublets4 <- doublets[
  hto_data["678-CD3-4h", doublets] >= min_exp &
    hto_data["678-unstim-0h", doublets] >= 2
]

doublets5 <- doublets[
  hto_data["678-CD3-4h", doublets] <= min_exp
]

#12, 74, 57, 72, 784
length(unique(c(doublets1, doublets2, doublets3, doublets4, doublets5)))


p <- data.frame(value = hto_data["640-unstim-0h", doublets5]) %>%
  ggplot(aes(x = value)) +
  geom_histogram(binwidth = 1, fill = color_vector["640-unstim-0h"], color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(x = "Exp of doublets <=min(678-CD3-4h) \nin 640-unstim-0h",
       y = "Frequency") +
  theme_bw()



pdf(paste0("qc/qc_NS-DD-1s-DEC-2/multiseqdemux/expression_doublets_high_678-CD3-4h_in_640-unstim-0h.pdf"), width = 2.97, height = 3.57)
print(p)
dev.off()

# Step 3: meet the conditions (n = 2279)
filtered_doublets <- setdiff(filtered_doublets, rownames(excluded_doublets))
seuratObj_demultiplex@meta.data[filtered_doublets, "MULTI_ID"] <- "678-CD3-4h"


Idents(seuratObj_demultiplex) <- "MULTI_ID"


seuratObj_demultiplex <- SetIdent(seuratObj_demultiplex,
                                  value = factor(Idents(seuratObj_demultiplex), levels = c(samples,
                                                                                           "Negative",
                                                                                           "Doublet")))

ridge_plots <- lapply(1:length(samples), function(i) {
  RidgePlot(seuratObj_demultiplex, assay = "HTO", features = rownames(seuratObj_demultiplex[["HTO"]])[i], ncol = 1) +
    scale_fill_manual(values = color_vector) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 8)) +
    guides(fill = "none")
})


# Combine the Ridge plots using patchwork
RidgePlot <- wrap_plots(ridge_plots, ncol = 3)

VlnPlot1 <- VlnPlot(seuratObj_demultiplex, features = "nFeature_RNA", pt.size = 0, log = TRUE) +
  geom_jitter(aes(y = nFeature_RNA), color = "grey", size = 0.05, alpha = 1, position = position_jitter(width = 0.3)) +  
  scale_fill_manual(values = color_vector) +
  theme_bw() +
  guides(fill = "none")

VlnPlot2 <- VlnPlot(seuratObj_demultiplex, features = "nCount_RNA", pt.size = 0, log = TRUE) +
  geom_jitter(aes(y = nCount_RNA), color = "grey", size = 0.05, alpha = 1, position = position_jitter(width = 0.3)) +  
  scale_fill_manual(values = color_vector) +
  theme_bw() +
  guides(fill = "none")


# counts of classification
table(seuratObj_demultiplex@meta.data$MULTI_ID) 

DefaultAssay(seuratObj_demultiplex) <- "HTO"
seuratObj_demultiplex <- ScaleData(seuratObj_demultiplex, features = rownames(seuratObj_demultiplex),
                                   verbose = FALSE)
seuratObj_demultiplex <- RunPCA(seuratObj_demultiplex, features = rownames(seuratObj_demultiplex), approx = FALSE)
seuratObj_demultiplex <- RunTSNE(seuratObj_demultiplex, dims = 1:3, check_duplicates = FALSE, perplexity = 100)


seuratObj_demultiplex@meta.data$MULTI_ID <-  factor(seuratObj_demultiplex@meta.data$MULTI_ID, levels = c(samples,
                                                                                                         "Negative",
                                                                                                         "Doublet"))
DimPlot1 <- DimPlot(seuratObj_demultiplex, group.by = "MULTI_ID") +
  scale_color_manual(values = color_vector) +
  theme_bw()

# DimPlot <- (DimPlot1 + DimPlot1) + plot_layout(ncol = 2)


MULTI_ID.df <- as.data.frame(table(seuratObj_demultiplex@meta.data$MULTI_ID))
colnames(MULTI_ID.df) <- c("MULTI_ID", "count")

RidgePlot_rel_height <- 1 + 1.3 * ceiling((length(samples) - 3) / 3) 


output_dir <- paste0("qc/qc_", lib, "/")


pdf(paste0(output_dir, "multiseqdemux/souporcell_manual_adjusted_optimal_pq_per_tag_margin_", margin, ".pdf"), width = 7.6, height = 9.76)


p1 <- plot_grid(RidgePlot, VlnPlot1, VlnPlot2,
                labels = c("A", "B", "C", "D"), 
                nrow = 3,
                rel_heights = c(RidgePlot_rel_height,1,1,1))


p2 <- plot_grid(DimPlot1, tableGrob(MULTI_ID.df),
                labels = c("E", "F"),
                nrow = 1,
                rel_widths = c(1.4, 1))

print(plot_grid(p1, p2,
                nrow = 2,
                rel_heights = c(2,1)))

dev.off()


# Stacked bar plot of souporcell calling for negatives and doublets 
negatives <- seuratObj_demultiplex@meta.data %>%
  filter(MULTI_ID == "Negative") %>%
  select(MULTI_ID, status, assignment) %>%
  mutate(assignment = case_when(
    assignment == "1" ~ "640",
    assignment == "0" ~ "678",
    assignment == "0/1" ~ "Doublet",
    assignment == "1/0" ~ "Doublet",
  ))

p_negatives <- negatives %>%
  group_by(MULTI_ID, assignment) %>%
  summarize(count = n(), .groups = "drop") %>%
  ggplot(aes(x = MULTI_ID, y = count, fill = assignment)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#ef476f", "#ffd166", "#118ab2")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(title = "Negatives Multi-seq",
       x = "Classification",
       y = "Count",
       fill = "Souporcell\nassignment") 

doublets_before <- seuratObj_demultiplex@meta.data %>%
  filter(grepl("_", MULTI_classification) | is.na(MULTI_classification)) %>%
  select(MULTI_ID, MULTI_classification, status, assignment) %>%
  mutate(assignment = case_when(
    assignment == "1" ~ "640",
    assignment == "0" ~ "678",
    assignment == "0/1" ~ "Doublet",
    assignment == "1/0" ~ "Doublet",
  ))

p_doublets_before <-  doublets_before %>%
  group_by(MULTI_classification, assignment) %>%
  summarize(count = n(), .groups = "drop") %>%
  ggplot(aes(x = MULTI_classification, y = count, fill = assignment)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#ef476f", "#ffd166", "#118ab2")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(title = "Doublets Multi-seq",
       x = "Classification",
       y = "Count",
       fill = "Souporcell assignment") 

doublets_after <- seuratObj_demultiplex@meta.data %>%
  filter(MULTI_ID == "Doublet") %>%
  select(MULTI_ID, MULTI_classification, status, assignment) %>%
  mutate(assignment = case_when(
    assignment == "1" ~ "640",
    assignment == "0" ~ "678",
    assignment == "0/1" ~ "Doublet",
    assignment == "1/0" ~ "Doublet",
  ))

p_doublets_after <- doublets_after %>%
  group_by(MULTI_classification, assignment) %>%
  summarize(count = n(), .groups = "drop") %>%
  ggplot(aes(x = MULTI_classification, y = count, fill = assignment)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#ef476f", "#ffd166", "#118ab2")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(title = "Doublets Multi-seq & curation",
       x = "Classification",
       y = "Count",
       fill = "Souporcell assignment")
# +
#   ylim(0, 400)


p_doublets_after_same_scale <- doublets_after %>%
  group_by(MULTI_classification, assignment) %>%
  summarize(count = n(), .groups = "drop") %>%
  ggplot(aes(x = MULTI_classification, y = count, fill = assignment)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#ef476f", "#ffd166", "#118ab2")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(title = "Doublets Multi-seq & curation",
       x = "Classification",
       y = "Count",
       fill = "Souporcell assignment") +
  ylim(0, 400)


pdf(paste0("qc/qc_NS-DD-1s-DEC-2/multiseqdemux/doublets_curation_multi-seq.pdf"), width = 10.94, height = 5.77)
p_doublets_before + p_doublets_after +
  plot_layout(guides = "collect") & # Shared legend
  theme(legend.position = "bottom", axis.title.x = element_blank())
dev.off()

pdf(paste0("qc/qc_NS-DD-1s-DEC-2/multiseqdemux/doublets_curation_multi-seq_same_y-axis.pdf"), width = 10.94, height = 5.77)
p_doublets_before + p_doublets_after_same_scale +
  plot_layout(guides = "collect") & # Shared legend
  theme(legend.position = "bottom", axis.title.x = element_blank())
dev.off()

pdf(paste0("qc/qc_NS-DD-1s-DEC-2/multiseqdemux/negatives_curation_multi-seq.pdf"), width = 2.97, height = 3.57)
p_negatives
dev.off()





# 
# negative <- subset(seuratObj_demultiplex, subset = MULTI_ID == "Negative")
# hto_data_negative <- as.data.frame(t(negative[["HTO"]]@data))
# 
# # Melt the data into long format
# hto_long <- melt(hto_data, variable.name = "Antibody", value.name = "Expression")
# 
# # Create a ridge plot
# ggplot(hto_long, aes(x = Expression, y = Antibody, fill = Antibody)) +
#   geom_density_ridges(scale = 2) +
#   theme_minimal() +
#   labs(title = "Expression of Antibodies in Negative Cells",
#        x = "Expression Level",
#        y = "Antibody") +
#   theme(legend.position = "none")



## 1.0.4 look into NS-DD-1s-DEC-2 after souporcell + manual curation + negative and small 640-unstim-0h clusters ----
# seuratObj_demultiplex <- plot_HTO_tag_counts(seuratObj_multiplex, margin, check.joint.bcs, plots=FALSE)
# 
# saveRDS(seuratObj_demultiplex, "qc/qc_NS-DD-1s-DEC-2/multiseqdemux/seuratObj_demultiplex.rds")
seuratObj_demultiplex <- readRDS("qc/qc_NS-DD-1s-DEC-2/multiseqdemux/seuratObj_demultiplex.rds")


souporcell_path <- "qc/qc_NS-DD-1s-DEC-2/souporcell_count/NS-DD-1s-DEC-2"
souporcell_clts <- read.table(paste0(souporcell_path, "/clusters.tsv"), header = TRUE, sep = "\t")

metadata <- seuratObj_demultiplex@meta.data
metadata$barcode <- rownames(metadata) 
metadata <- merge(metadata, souporcell_clts, by = "barcode", all = TRUE)  
all(rownames(seuratObj_demultiplex@meta.data) == metadata$barcode)

rownames(metadata) <- metadata$barcode
seuratObj_demultiplex@meta.data <- metadata %>% select(-barcode)

# step 1: min exp of doublets in 678-CD3-4h that can possibly be consider as 678-CD3-4h
selected_cells <- rownames(seuratObj_demultiplex@meta.data[
  seuratObj_demultiplex@meta.data$MULTI_ID == "678-CD3-4h", 
])

min_exp <- min(seuratObj_demultiplex[["HTO"]]@data["678-CD3-4h", selected_cells], na.rm = TRUE)

# step 2, apply filters on current doublets (n = 3274)
doublets <- rownames(seuratObj_demultiplex@meta.data[
  seuratObj_demultiplex@meta.data$MULTI_ID == "Doublet", 
])

# length(rownames(seuratObj_demultiplex@meta.data[
#   grep("678-CD3-4h", seuratObj_demultiplex@meta.data$MULTI_classification),
# ]))
# 
# df <- as.data.frame(table(seuratObj_demultiplex@meta.data$MULTI_classification))
# 
# length(rownames(seuratObj_demultiplex@meta.data[
#   is.na(seuratObj_demultiplex@meta.data$MULTI_classification), 
# ]))


hto_data <- seuratObj_demultiplex[["HTO"]]@data

# Step 2: meet the conditions (n = 2279)
filtered_doublets <- doublets[
  hto_data["678-CD3-4h", doublets] > min_exp & # Above min_exp
    hto_data["640-Ig-4h", doublets] < 2.5 &     # <2.5 in 640-Ig-4h
    hto_data["640-CD3-4h", doublets] < 2 &     # <2 in 640-CD3-4h
    hto_data["678-Ig-4h", doublets] < 2 &      # <2 in 678-Ig-4h
    hto_data["678-unstim-0h", doublets] < 2    # <2 in 678-unstim-0h
]

excluded_doublets <- seuratObj_demultiplex@meta.data %>%
  filter(MULTI_ID == "Doublet", 
         grepl("678-CD3-4h", MULTI_classification),
         grepl("640-unstim-0h", MULTI_classification)) %>%
  mutate(assignment = case_when(
    assignment == "1" ~ "640",
    assignment == "0" ~ "678",
    assignment == "0/1" ~ "Doublet",
    assignment == "1/0" ~ "Doublet",
  )) %>%
  filter(assignment == "640")



doublets1 <- doublets[
  hto_data["678-CD3-4h", doublets] >= min_exp &
    hto_data["640-Ig-4h", doublets] >= 2.5
]

doublets2 <- doublets[
  hto_data["678-CD3-4h", doublets] >= min_exp &
    hto_data["640-CD3-4h", doublets] >= 2
]

doublets3 <- doublets[
  hto_data["678-CD3-4h", doublets] >= min_exp &
    hto_data["678-Ig-4h", doublets] >= 2
]

doublets4 <- doublets[
  hto_data["678-CD3-4h", doublets] >= min_exp &
    hto_data["678-unstim-0h", doublets] >= 2
]

doublets5 <- doublets[
  hto_data["678-CD3-4h", doublets] <= min_exp
]

#12, 74, 57, 72, 784
length(unique(c(doublets1, doublets2, doublets3, doublets4, doublets5)))


p <- data.frame(value = hto_data["640-unstim-0h", doublets5]) %>%
  ggplot(aes(x = value)) +
  geom_histogram(binwidth = 1, fill = color_vector["640-unstim-0h"], color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(x = "Exp of doublets <=min(678-CD3-4h) \nin 640-unstim-0h",
       y = "Frequency") +
  theme_bw()



pdf(paste0("qc/qc_NS-DD-1s-DEC-2/multiseqdemux/expression_doublets_high_678-CD3-4h_in_640-unstim-0h.pdf"), width = 2.97, height = 3.57)
print(p)
dev.off()

# Step 3: meet the conditions (n = 2279)
filtered_doublets <- setdiff(filtered_doublets, rownames(excluded_doublets))
seuratObj_demultiplex@meta.data[filtered_doublets, "MULTI_ID"] <- "678-CD3-4h"


# step 4 run celltypist
rna_data <- GetAssayData(seuratObj_demultiplex, assay = "RNA", layer = "data")
write.csv(as.matrix(rna_data), file = "qc/qc_NS-DD-1s-DEC-2/multiseqdemux/rna_data.csv", quote = F)


# step 5 figure out the coordinates of negative and small 640-unstim-0h clusters -> cluster1

DefaultAssay(seuratObj_demultiplex) <- "HTO"
seuratObj_demultiplex <- ScaleData(seuratObj_demultiplex, features = rownames(seuratObj_demultiplex),
                                   verbose = FALSE)
seuratObj_demultiplex <- RunPCA(seuratObj_demultiplex, features = rownames(seuratObj_demultiplex), approx = FALSE)
seuratObj_demultiplex <- RunTSNE(seuratObj_demultiplex, dims = 1:3, check_duplicates = FALSE, perplexity = 100)


seuratObj_demultiplex@meta.data$MULTI_ID <-  factor(seuratObj_demultiplex@meta.data$MULTI_ID, levels = c(samples,
                                                                                                         "Negative",
                                                                                                         "Doublet"))
DimPlot1 <- DimPlot(seuratObj_demultiplex, group.by = "MULTI_ID") +
  scale_color_manual(values = color_vector) +
  theme_bw()


p <- DimPlot1 + coord_cartesian(xlim = c(-40.5, -30),
                           ylim = c(-20, -5)) +
  scale_x_continuous(breaks = seq(floor(min(seuratObj_demultiplex@reductions$tsne@cell.embeddings[, 1])),
                                  ceiling(max(seuratObj_demultiplex@reductions$tsne@cell.embeddings[, 1])), 2)) +
  scale_y_continuous(breaks = seq(floor(min(seuratObj_demultiplex@reductions$tsne@cell.embeddings[, 2])),
                                  ceiling(max(seuratObj_demultiplex@reductions$tsne@cell.embeddings[, 2])),2))
output_dir <- paste0("qc/qc_", lib, "/")
pdf(paste0(output_dir, "multiseqdemux/souporcell_manual_adjusted_optimal_pq_per_tag_margin_", margin, "_tsne_cluster2.pdf"), width = 3.80, height = 2.54)
print(p)
dev.off()

# step 6 look into the negative -> cluster 1 and small 640-unstim-0h clusters -> cluster2
# definition of cluster1 cells: "Negative" in MULTI_ID & TSNE2 > 0
# definition of cluster2 cells: -41 <= TSNE1 <= -30 & -20 <= TSNE2 <= -5 & not doublets
cell_annot <- read.table("qc/qc_NS-DD-1s-DEC-2/multiseqdemux/celltypist_predicted_labels.csv", header = TRUE, sep = ",")
cell_annot <- cell_annot %>% 
  rename("barcode" = "X")

metadata <- seuratObj_demultiplex@meta.data
metadata$barcode <- rownames(metadata) 
metadata <- merge(metadata, cell_annot, by = "barcode", all = TRUE)  
all(rownames(seuratObj_demultiplex@meta.data) == metadata$barcode)

rownames(metadata) <- metadata$barcode
seuratObj_demultiplex@meta.data <- metadata %>% select(-barcode)

seuratObj_demultiplex@meta.data$MULTI_ID_2 <- seuratObj_demultiplex@meta.data$MULTI_ID
seuratObj_demultiplex@meta.data$MULTI_ID_2 <- as.character(seuratObj_demultiplex@meta.data$MULTI_ID_2)

seuratObj_demultiplex@meta.data$MULTI_ID_2[seuratObj_demultiplex@meta.data$MULTI_ID == "Negative" & 
                                             seuratObj_demultiplex@reductions$tsne@cell.embeddings[, 2] > 0] <- "cluster1"

seuratObj_demultiplex@meta.data$MULTI_ID_2[seuratObj_demultiplex@reductions$tsne@cell.embeddings[, 1] >= -41 & 
                                             seuratObj_demultiplex@reductions$tsne@cell.embeddings[, 1] <= -30 & 
                                             seuratObj_demultiplex@reductions$tsne@cell.embeddings[, 2] >= -20 & 
                                             seuratObj_demultiplex@reductions$tsne@cell.embeddings[, 2] <= -5 & 
                                             seuratObj_demultiplex@meta.data$MULTI_ID != "Doublet"] <- "cluster2"



library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


p <- seuratObj_demultiplex@meta.data %>%
  filter(MULTI_ID_2 %in% c("cluster1", "cluster2")) %>%
  group_by(MULTI_ID_2, predicted_labels) %>%
  count(name = "count")  %>%
  # mutate(predicted_labels = factor(predicted_labels, 
  #                                  levels = unique(predicted_labels[order(-count)]))) %>%  # Reorder by count
  
  ggplot(aes(x = MULTI_ID_2, y = count, fill = factor(predicted_labels, 
                                                      levels = unique(predicted_labels[order(-count)])))) +
  geom_bar(stat = "identity") +
  labs(x = "Cluster", y = "Count", 
       title = "Predicted Labels for cluster1 (negative, possibly Ig) & cluster2 (unstim)",
       fill = "Predicted labels") +
  theme_bw() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  guides(fill = guide_legend(ncol = 2))  


pdf(paste0(output_dir, "multiseqdemux/celltypist_cluster1_cluster2.pdf"), width = 6.59, height = 4.14)
print(p)
dev.off()

seuratObj_demultiplex@meta.data$MULTI_ID_3 <- seuratObj_demultiplex@meta.data$MULTI_ID
seuratObj_demultiplex@meta.data$MULTI_ID_3 <- as.character(seuratObj_demultiplex@meta.data$MULTI_ID_3)




seuratObj_demultiplex@meta.data$MULTI_ID_3[seuratObj_demultiplex@meta.data$MULTI_ID_2 == "cluster1" &
                                             seuratObj_demultiplex@meta.data$predicted_labels == "B cells" &
                                             seuratObj_demultiplex@meta.data$assignment == "1"] <- "640-Ig-4h"
  

seuratObj_demultiplex@meta.data$MULTI_ID <-  factor(seuratObj_demultiplex@meta.data$MULTI_ID, levels = c(samples,
                                                                                                         "Negative",
                                                                                                         "Doublet"))

seuratObj_demultiplex@meta.data$MULTI_ID_3 <-  factor(seuratObj_demultiplex@meta.data$MULTI_ID_3, levels = c(samples,
                                                                                                         "Negative",
                                                                                                         "Doublet"))
DimPlot1 <- DimPlot(seuratObj_demultiplex, group.by = "MULTI_ID") +
  coord_cartesian(
    xlim = c(min(seuratObj_demultiplex@reductions$tsne@cell.embeddings[, 1]), 25),
                  ylim = c(0, 40)) +
  scale_color_manual(values = color_vector) +
  theme_bw()


DimPlot2 <- DimPlot(seuratObj_demultiplex, group.by = "MULTI_ID_3") +
  coord_cartesian(
    xlim = c(min(seuratObj_demultiplex@reductions$tsne@cell.embeddings[, 1]), 25),
             ylim = c(0, 40)) +
  scale_color_manual(values = color_vector) +
  theme_bw()




pdf(paste0(output_dir, "multiseqdemux/tsne_cluster1_MULTI_ID_3.pdf"), width = 10.82, height = 4.39)
plot_grid(DimPlot1, DimPlot2,
          labels = c("A", "B"),
          nrow = 1,
          rel_widths = c(1, 1))
dev.off()


# plot the expression level of a specific antibody
FeaturePlot1 <- FeaturePlot(
  seuratObj_demultiplex, 
  features = samples[1], 
  reduction = "tsne"
)

library(patchwork)



plots <- lapply(samples, function(feature) {
  FeaturePlot(
    seuratObj_demultiplex, 
    features = feature, 
    reduction = "tsne"
  ) +
    # scale_color_gradient(low = "#F5F5F5", high = darkened_color_vector[feature]) +
    scale_color_gradient(low = "#D3D3D3", high = "blue") +
    
    ggtitle(paste("Expression of", feature)) +
    theme_bw()
})

pdf(paste0(output_dir, "multiseqdemux/tsne_per_antibody_exp.pdf"), width = 10.72, height = 6.13)

wrap_plots(plots)

dev.off()


seuratObj_demultiplex@meta.data$range_nCount_RNA[seuratObj_demultiplex@meta.data$nCount_RNA >= 25000] <- ">=25K"
seuratObj_demultiplex@meta.data$range_nCount_RNA[seuratObj_demultiplex@meta.data$nCount_RNA >= 20000 & 
                                                   seuratObj_demultiplex@meta.data$nCount_RNA < 25000] <- "20K"
seuratObj_demultiplex@meta.data$range_nCount_RNA[seuratObj_demultiplex@meta.data$nCount_RNA >= 15000 & 
                                                   seuratObj_demultiplex@meta.data$nCount_RNA < 20000] <- "15K"
seuratObj_demultiplex@meta.data$range_nCount_RNA[seuratObj_demultiplex@meta.data$nCount_RNA >= 10000 & 
                                                   seuratObj_demultiplex@meta.data$nCount_RNA < 15000] <- "10K"
seuratObj_demultiplex@meta.data$range_nCount_RNA[seuratObj_demultiplex@meta.data$nCount_RNA >= 7500 & 
                                                   seuratObj_demultiplex@meta.data$nCount_RNA < 10000] <- "7.5K"
seuratObj_demultiplex@meta.data$range_nCount_RNA[seuratObj_demultiplex@meta.data$nCount_RNA >= 5000 & 
                                                   seuratObj_demultiplex@meta.data$nCount_RNA < 7500] <- "5K"
seuratObj_demultiplex@meta.data$range_nCount_RNA[seuratObj_demultiplex@meta.data$nCount_RNA >= 2500 & 
                                                   seuratObj_demultiplex@meta.data$nCount_RNA < 5000] <- "2.5K"
seuratObj_demultiplex@meta.data$range_nCount_RNA[seuratObj_demultiplex@meta.data$nCount_RNA < 2500] <- "<2.5K"

seuratObj_demultiplex@meta.data$range_nCount_RNA <- factor(
  seuratObj_demultiplex@meta.data$range_nCount_RNA,
  levels = c("<2.5K", "2.5K", "5K", "7.5K", "10K", "15K", "20K", ">=25K"))


p2 <- DimPlot(
  seuratObj_demultiplex, 
  group.by = "range_nCount_RNA",  # Specify the categorical variable
  cols = c("<2.5K" = "lightgrey", "2.5K" = "blue", "5K" = "green", 
           "7.5K" = "yellow", "10K" = "orange", "15K" = "red", 
           "20K" = "purple", ">=25K" = "darkred")  # Assign colors to each range
) +
  ggtitle("nCount RNA") +
  theme_bw()


p1 <- FeaturePlot(
  seuratObj_demultiplex,  
  features = "nCount_RNA",  
  reduction = "tsne"  
) +
  # scale_color_viridis() +
  ggtitle("nCount RNA") + 
  theme_bw()  #


pdf(paste0(output_dir, "multiseqdemux/tsne_nCount_RNA.pdf"), width = 10.82, height = 4.39)
plot_grid(p1, p2,
          labels = c("A", "B"),
          nrow = 1,
          rel_widths = c(1, 1))
dev.off()



# step 7 souporcell of negatives MULTI_ID_4

negatives <-seuratObj_demultiplex@meta.data %>%
  filter(MULTI_ID == "Negative") %>%
  select(MULTI_ID, status, assignment, predicted_labels) %>%
  mutate(assignment = case_when(
    assignment == "1" ~ "640",
    assignment == "0" ~ "678",
    assignment == "0/1" ~ "Doublet",
    assignment == "1/0" ~ "Doublet",
  ))

negatives %>%
  group_by(predicted_labels, assignment) %>%
  summarize(count = n(), .groups = "drop")  %>%
  ggplot(aes(x = assignment, y = count, fill = factor(predicted_labels, 
                                                      levels = unique(predicted_labels[order(-count)])))) +
  geom_bar(stat = "identity") +
  # scale_fill_manual(values = c("#ef476f", "#ffd166", "#118ab2")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Predicted Labels for cluster1 (negative, possibly Ig)",
       x = "Souporcell classification",
       y = "Count",
       fill = "Predicted labels") +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  guides(fill = guide_legend(ncol = 2))  



# assign negatives that are 640 by souporcell as 640_Ig_4h
seuratObj_demultiplex@meta.data$MULTI_ID_4 <- seuratObj_demultiplex@meta.data$MULTI_ID
seuratObj_demultiplex@meta.data$MULTI_ID_4 <- as.character(seuratObj_demultiplex@meta.data$MULTI_ID_4)
seuratObj_demultiplex@meta.data$MULTI_ID_4[seuratObj_demultiplex@meta.data$MULTI_ID == "Negative" &
                                             seuratObj_demultiplex@meta.data$assignment == "1"] <- "640-Ig-4h"
seuratObj_demultiplex@meta.data$MULTI_ID <-  factor(seuratObj_demultiplex@meta.data$MULTI_ID, levels = c(samples,
                                                                                                         "Negative",
                                                                                                         "Doublet"))

seuratObj_demultiplex@meta.data$MULTI_ID_4 <-  factor(seuratObj_demultiplex@meta.data$MULTI_ID_4, levels = c(samples,
                                                                                                             "Negative",
                                                                                                             "Doublet"))

seuratObj_demultiplex@meta.data %>%
  filter(MULTI_ID_4 == "Negative") %>%
  select(MULTI_ID_4, status, assignment) %>%
  mutate(assignment = case_when(
    assignment == "1" ~ "640",
    assignment == "0" ~ "678",
    assignment == "0/1" ~ "Doublet",
    assignment == "1/0" ~ "Doublet",
  ))%>%
  group_by(MULTI_ID_4, assignment) %>%
  summarize(count = n(), .groups = "drop") %>%
  ggplot(aes(x = MULTI_ID_4, y = count, fill = assignment)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#ef476f", "#ffd166", "#118ab2")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "678-CD3-4h doublets\nMulti-seq",
       x = "Classification",
       y = "Count",
       fill = "Souporcell\nassignment") 



DimPlot1 <- DimPlot(seuratObj_demultiplex, group.by = "MULTI_ID") +
  coord_cartesian(
    xlim = c(min(seuratObj_demultiplex@reductions$tsne@cell.embeddings[, 1]), 25),
    ylim = c(0, 40)) +
  scale_color_manual(values = color_vector) +
  theme_bw()


DimPlot2 <- DimPlot(seuratObj_demultiplex, group.by = "MULTI_ID_4") +
  coord_cartesian(
    xlim = c(min(seuratObj_demultiplex@reductions$tsne@cell.embeddings[, 1]), 25),
    ylim = c(0, 40)) +
  scale_color_manual(values = color_vector) +
  theme_bw()



pdf(paste0(output_dir, "multiseqdemux/tsne_cluster1_MULTI_ID_4.pdf"), width = 10.82, height = 4.39)
plot_grid(DimPlot1, DimPlot2,
          labels = c("A", "B"),
          nrow = 1,
          rel_widths = c(1, 1))
dev.off()



# step 8 souporcell of negatives MULTI_ID_5, exclduing 640-CD3-4h and UMI counts < 2.5k -> final -> check out 1.0.4 

## FINAL 1.0.5 look into NS-DD-1s-DEC-2 final -> check out 1.0.4 ####
# many of the codes are the same with the previous step
# put them all together here 

seuratObj_demultiplex <- readRDS("qc/qc_NS-DD-1s-DEC-2/multiseqdemux/seuratObj_demultiplex.rds")


souporcell_path <- "qc/qc_NS-DD-1s-DEC-2/souporcell_count/NS-DD-1s-DEC-2"
souporcell_clts <- read.table(paste0(souporcell_path, "/clusters.tsv"), header = TRUE, sep = "\t")

metadata <- seuratObj_demultiplex@meta.data
metadata$barcode <- rownames(metadata) 
metadata <- merge(metadata, souporcell_clts, by = "barcode", all = TRUE)  
all(rownames(seuratObj_demultiplex@meta.data) == metadata$barcode)

rownames(metadata) <- metadata$barcode
seuratObj_demultiplex@meta.data <- metadata %>% select(-barcode)

# step 1: min exp of doublets in 678-CD3-4h that can possibly be consider as 678-CD3-4h
selected_cells <- rownames(seuratObj_demultiplex@meta.data[
  seuratObj_demultiplex@meta.data$MULTI_ID == "678-CD3-4h", 
])

min_exp <- min(seuratObj_demultiplex[["HTO"]]@data["678-CD3-4h", selected_cells], na.rm = TRUE)

# step 2, apply filters on current doublets (n = 3274)
doublets <- rownames(seuratObj_demultiplex@meta.data[
  seuratObj_demultiplex@meta.data$MULTI_ID == "Doublet", 
])

# length(rownames(seuratObj_demultiplex@meta.data[
#   grep("678-CD3-4h", seuratObj_demultiplex@meta.data$MULTI_classification),
# ]))
# 
# df <- as.data.frame(table(seuratObj_demultiplex@meta.data$MULTI_classification))
# 
# length(rownames(seuratObj_demultiplex@meta.data[
#   is.na(seuratObj_demultiplex@meta.data$MULTI_classification), 
# ]))


hto_data <- seuratObj_demultiplex[["HTO"]]@data

# Step 2: meet the conditions (n = 2279)
filtered_doublets <- doublets[
  hto_data["678-CD3-4h", doublets] > min_exp & # Above min_exp
    hto_data["640-Ig-4h", doublets] < 2.5 &     # <2.5 in 640-Ig-4h
    hto_data["640-CD3-4h", doublets] < 2 &     # <2 in 640-CD3-4h
    hto_data["678-Ig-4h", doublets] < 2 &      # <2 in 678-Ig-4h
    hto_data["678-unstim-0h", doublets] < 2    # <2 in 678-unstim-0h
]

excluded_doublets <- seuratObj_demultiplex@meta.data %>%
  filter(MULTI_ID == "Doublet", 
         grepl("678-CD3-4h", MULTI_classification),
         grepl("640-unstim-0h", MULTI_classification)) %>%
  mutate(assignment = case_when(
    assignment == "1" ~ "640",
    assignment == "0" ~ "678",
    assignment == "0/1" ~ "Doublet",
    assignment == "1/0" ~ "Doublet",
  )) %>%
  filter(assignment == "640")



doublets1 <- doublets[
  hto_data["678-CD3-4h", doublets] >= min_exp &
    hto_data["640-Ig-4h", doublets] >= 2.5
]

doublets2 <- doublets[
  hto_data["678-CD3-4h", doublets] >= min_exp &
    hto_data["640-CD3-4h", doublets] >= 2
]

doublets3 <- doublets[
  hto_data["678-CD3-4h", doublets] >= min_exp &
    hto_data["678-Ig-4h", doublets] >= 2
]

doublets4 <- doublets[
  hto_data["678-CD3-4h", doublets] >= min_exp &
    hto_data["678-unstim-0h", doublets] >= 2
]

doublets5 <- doublets[
  hto_data["678-CD3-4h", doublets] <= min_exp
]

#12, 74, 57, 72, 784
length(unique(c(doublets1, doublets2, doublets3, doublets4, doublets5)))


p <- data.frame(value = hto_data["640-unstim-0h", doublets5]) %>%
  ggplot(aes(x = value)) +
  geom_histogram(binwidth = 1, fill = color_vector["640-unstim-0h"], color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(x = "Exp of doublets <=min(678-CD3-4h) \nin 640-unstim-0h",
       y = "Frequency") +
  theme_bw()



pdf(paste0("qc/qc_NS-DD-1s-DEC-2/multiseqdemux/FINAL/expression_doublets_high_678-CD3-4h_in_640-unstim-0h.pdf"), width = 2.97, height = 3.57)
print(p)
dev.off()

# Step 3: meet the conditions (n = 2279)
filtered_doublets <- setdiff(filtered_doublets, rownames(excluded_doublets))
seuratObj_demultiplex@meta.data[filtered_doublets, "MULTI_ID"] <- "678-CD3-4h"

# step 4 run celltypist
write.csv(as.matrix(seuratObj_demultiplex@assays$RNA$counts), file = "qc/qc_NS-DD-1s-DEC-2/multiseqdemux/rna_counts.csv", quote = F)


# step 5 figure out the coordinates of negative and small 640-unstim-0h clusters -> cluster1

DefaultAssay(seuratObj_demultiplex) <- "HTO"
seuratObj_demultiplex <- ScaleData(seuratObj_demultiplex, features = rownames(seuratObj_demultiplex),
                                   verbose = FALSE)
seuratObj_demultiplex <- RunPCA(seuratObj_demultiplex, features = rownames(seuratObj_demultiplex), approx = FALSE)
seuratObj_demultiplex <- RunTSNE(seuratObj_demultiplex, dims = 1:3, check_duplicates = FALSE, perplexity = 100)


seuratObj_demultiplex@meta.data$MULTI_ID <-  factor(seuratObj_demultiplex@meta.data$MULTI_ID, levels = c(samples,
                                                                                                         "Negative",
                                                                                                         "Doublet"))
DimPlot1 <- DimPlot(seuratObj_demultiplex, group.by = "MULTI_ID") +
  scale_color_manual(values = color_vector) +
  theme_bw()


p <- DimPlot1 + coord_cartesian(xlim = c(-40.5, -30),
                                ylim = c(-20, -5)) +
  scale_x_continuous(breaks = seq(floor(min(seuratObj_demultiplex@reductions$tsne@cell.embeddings[, 1])),
                                  ceiling(max(seuratObj_demultiplex@reductions$tsne@cell.embeddings[, 1])), 2)) +
  scale_y_continuous(breaks = seq(floor(min(seuratObj_demultiplex@reductions$tsne@cell.embeddings[, 2])),
                                  ceiling(max(seuratObj_demultiplex@reductions$tsne@cell.embeddings[, 2])),2))
output_dir <- paste0("qc/qc_", lib, "/")
pdf(paste0(output_dir, "multiseqdemux/FINAL/souporcell_manual_adjusted_optimal_pq_per_tag_margin_", margin, "_tsne_cluster2.pdf"), width = 3.80, height = 2.54)
print(p)
dev.off()

# step 6 look into the negative -> cluster 1 and small 640-unstim-0h clusters -> cluster2
# definition of cluster1 cells: "Negative" in MULTI_ID & TSNE2 > 0
# definition of cluster2 cells: -41 <= TSNE1 <= -30 & -20 <= TSNE2 <= -5 & not doublets
cell_annot <- read.table("qc/qc_NS-DD-1s-DEC-2/multiseqdemux/celltypist_predicted_labels.csv", header = TRUE, sep = ",")
cell_annot <- cell_annot %>% 
  rename("barcode" = "X")

metadata <- seuratObj_demultiplex@meta.data
metadata$barcode <- rownames(metadata) 
metadata <- merge(metadata, cell_annot, by = "barcode", all = TRUE)  
all(rownames(seuratObj_demultiplex@meta.data) == metadata$barcode)

rownames(metadata) <- metadata$barcode
seuratObj_demultiplex@meta.data <- metadata %>% select(-barcode)

# seuratObj_demultiplex@meta.data$MULTI_ID_2 <- seuratObj_demultiplex@meta.data$MULTI_ID
# seuratObj_demultiplex@meta.data$MULTI_ID_2 <- as.character(seuratObj_demultiplex@meta.data$MULTI_ID_2)
# 
# seuratObj_demultiplex@meta.data$MULTI_ID_2[seuratObj_demultiplex@meta.data$MULTI_ID == "Negative" & 
#                                              seuratObj_demultiplex@reductions$tsne@cell.embeddings[, 2] > 0] <- "cluster1"
# 
# seuratObj_demultiplex@meta.data$MULTI_ID_2[seuratObj_demultiplex@reductions$tsne@cell.embeddings[, 1] >= -41 & 
#                                              seuratObj_demultiplex@reductions$tsne@cell.embeddings[, 1] <= -30 & 
#                                              seuratObj_demultiplex@reductions$tsne@cell.embeddings[, 2] >= -20 & 
#                                              seuratObj_demultiplex@reductions$tsne@cell.embeddings[, 2] <= -5 & 
#                                              seuratObj_demultiplex@meta.data$MULTI_ID != "Doublet"] <- "cluster2"
# 


# library(RColorBrewer)
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# 
# 
# p <- seuratObj_demultiplex@meta.data %>%
#   filter(MULTI_ID_2 %in% c("cluster1", "cluster2")) %>%
#   group_by(MULTI_ID_2, predicted_labels) %>%
#   count(name = "count")  %>%
#   # mutate(predicted_labels = factor(predicted_labels, 
#   #                                  levels = unique(predicted_labels[order(-count)]))) %>%  # Reorder by count
#   
#   ggplot(aes(x = MULTI_ID_2, y = count, fill = factor(predicted_labels, 
#                                                       levels = unique(predicted_labels[order(-count)])))) +
#   geom_bar(stat = "identity") +
#   labs(x = "Cluster", y = "Count", 
#        title = "Predicted Labels for cluster1 (negative, possibly Ig) & cluster2 (unstim)",
#        fill = "Predicted labels") +
#   theme_bw() +
#   scale_fill_manual(values = cols) +
#   theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
#   guides(fill = guide_legend(ncol = 2))  
# 
# 
# pdf(paste0(output_dir, "multiseqdemux/celltypist_cluster1_cluster2.pdf"), width = 6.59, height = 4.14)
# print(p)
# dev.off()
# 
# seuratObj_demultiplex@meta.data$MULTI_ID_3 <- seuratObj_demultiplex@meta.data$MULTI_ID
# seuratObj_demultiplex@meta.data$MULTI_ID_3 <- as.character(seuratObj_demultiplex@meta.data$MULTI_ID_3)
# 
# 
# 
# 
# seuratObj_demultiplex@meta.data$MULTI_ID_3[seuratObj_demultiplex@meta.data$MULTI_ID_2 == "cluster1" &
#                                              seuratObj_demultiplex@meta.data$predicted_labels == "B cells" &
#                                              seuratObj_demultiplex@meta.data$assignment == "1"] <- "640-Ig-4h"
# 
# 
# seuratObj_demultiplex@meta.data$MULTI_ID <-  factor(seuratObj_demultiplex@meta.data$MULTI_ID, levels = c(samples,
#                                                                                                          "Negative",
#                                                                                                          "Doublet"))
# 
# seuratObj_demultiplex@meta.data$MULTI_ID_3 <-  factor(seuratObj_demultiplex@meta.data$MULTI_ID_3, levels = c(samples,
#                                                                                                              "Negative",
#                                                                                                              "Doublet"))
# DimPlot1 <- DimPlot(seuratObj_demultiplex, group.by = "MULTI_ID") +
#   coord_cartesian(
#     xlim = c(min(seuratObj_demultiplex@reductions$tsne@cell.embeddings[, 1]), 25),
#     ylim = c(0, 40)) +
#   scale_color_manual(values = color_vector) +
#   theme_bw()
# 
# 
# DimPlot2 <- DimPlot(seuratObj_demultiplex, group.by = "MULTI_ID_3") +
#   coord_cartesian(
#     xlim = c(min(seuratObj_demultiplex@reductions$tsne@cell.embeddings[, 1]), 25),
#     ylim = c(0, 40)) +
#   scale_color_manual(values = color_vector) +
#   theme_bw()
# 
# 
# 
# 
# pdf(paste0(output_dir, "multiseqdemux/tsne_cluster1_MULTI_ID_3.pdf"), width = 10.82, height = 4.39)
# plot_grid(DimPlot1, DimPlot2,
#           labels = c("A", "B"),
#           nrow = 1,
#           rel_widths = c(1, 1))
# dev.off()


# plot the expression level of a specific antibody
FeaturePlot1 <- FeaturePlot(
  seuratObj_demultiplex, 
  features = samples[1], 
  reduction = "tsne"
)

library(patchwork)



plots <- lapply(samples, function(feature) {
  FeaturePlot(
    seuratObj_demultiplex, 
    features = feature, 
    reduction = "tsne"
  ) +
    # scale_color_gradient(low = "#F5F5F5", high = darkened_color_vector[feature]) +
    scale_color_gradient(low = "#D3D3D3", high = "blue") +
    
    ggtitle(paste("Expression of", feature)) +
    theme_bw()
})

pdf(paste0(output_dir, "multiseqdemux/FINAL/tsne_per_antibody_exp.pdf"), width = 10.72, height = 6.13)

wrap_plots(plots)

dev.off()


seuratObj_demultiplex@meta.data$range_nCount_RNA[seuratObj_demultiplex@meta.data$nCount_RNA >= 25000] <- ">=25K"
seuratObj_demultiplex@meta.data$range_nCount_RNA[seuratObj_demultiplex@meta.data$nCount_RNA >= 20000 & 
                                                   seuratObj_demultiplex@meta.data$nCount_RNA < 25000] <- "20K"
seuratObj_demultiplex@meta.data$range_nCount_RNA[seuratObj_demultiplex@meta.data$nCount_RNA >= 15000 & 
                                                   seuratObj_demultiplex@meta.data$nCount_RNA < 20000] <- "15K"
seuratObj_demultiplex@meta.data$range_nCount_RNA[seuratObj_demultiplex@meta.data$nCount_RNA >= 10000 & 
                                                   seuratObj_demultiplex@meta.data$nCount_RNA < 15000] <- "10K"
seuratObj_demultiplex@meta.data$range_nCount_RNA[seuratObj_demultiplex@meta.data$nCount_RNA >= 7500 & 
                                                   seuratObj_demultiplex@meta.data$nCount_RNA < 10000] <- "7.5K"
seuratObj_demultiplex@meta.data$range_nCount_RNA[seuratObj_demultiplex@meta.data$nCount_RNA >= 5000 & 
                                                   seuratObj_demultiplex@meta.data$nCount_RNA < 7500] <- "5K"
seuratObj_demultiplex@meta.data$range_nCount_RNA[seuratObj_demultiplex@meta.data$nCount_RNA >= 2500 & 
                                                   seuratObj_demultiplex@meta.data$nCount_RNA < 5000] <- "2.5K"
seuratObj_demultiplex@meta.data$range_nCount_RNA[seuratObj_demultiplex@meta.data$nCount_RNA < 2500] <- "<2.5K"

seuratObj_demultiplex@meta.data$range_nCount_RNA <- factor(
  seuratObj_demultiplex@meta.data$range_nCount_RNA,
  levels = c("<2.5K", "2.5K", "5K", "7.5K", "10K", "15K", "20K", ">=25K"))


p2 <- DimPlot(
  seuratObj_demultiplex, 
  group.by = "range_nCount_RNA",  # Specify the categorical variable
  cols = c("<2.5K" = "lightgrey", "2.5K" = "blue", "5K" = "green", 
           "7.5K" = "yellow", "10K" = "orange", "15K" = "red", 
           "20K" = "purple", ">=25K" = "darkred")  # Assign colors to each range
) +
  ggtitle("nCount RNA") +
  theme_bw()


p1 <- FeaturePlot(
  seuratObj_demultiplex,  
  features = "nCount_RNA",  
  reduction = "tsne"  
) +
  # scale_color_viridis() +
  ggtitle("nCount RNA") + 
  theme_bw()  #


pdf(paste0(output_dir, "multiseqdemux/tsne_nCount_RNA.pdf"), width = 10.82, height = 4.39)
plot_grid(p1, p2,
          labels = c("A", "B"),
          nrow = 1,
          rel_widths = c(1, 1))
dev.off()

# step 7 deal with negatives




# final step to make all relevant plots 


qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


negatives <- seuratObj_demultiplex@meta.data %>%
  filter(MULTI_ID == "Negative") %>%
  select(MULTI_ID, status, assignment, predicted_labels) %>%
  mutate(assignment = case_when(
    assignment == "1" ~ "640",
    assignment == "0" ~ "678",
    assignment == "0/1" ~ "Doublet",
    assignment == "1/0" ~ "Doublet",
  ))

negatives %>%
  group_by(predicted_labels, assignment) %>%
  summarize(count = n(), .groups = "drop")  %>%
  ggplot(aes(x = assignment, y = count, fill = factor(predicted_labels, 
                                                      levels = unique(predicted_labels[order(-count)])))) +
  geom_bar(stat = "identity") +
  # scale_fill_manual(values = c("#ef476f", "#ffd166", "#118ab2")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Predicted Labels for cluster1 (negative, possibly Ig)",
       x = "Souporcell classification",
       y = "Count",
       fill = "Predicted labels") +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  guides(fill = guide_legend(ncol = 2))  



# assign some of the negatives to 640_Ig_4h
seuratObj_demultiplex@meta.data$MULTI_ID <- as.character(seuratObj_demultiplex@meta.data$MULTI_ID)
seuratObj_demultiplex@meta.data$MULTI_ID[seuratObj_demultiplex@meta.data$MULTI_ID == "Negative" &
                                             seuratObj_demultiplex@meta.data$assignment == "1" &
                                             seuratObj_demultiplex@meta.data$nCount_RNA >= 2500 &
                                             !(seuratObj_demultiplex@reductions$tsne@cell.embeddings[, 1] < -25 & 
                                                 seuratObj_demultiplex@reductions$tsne@cell.embeddings[, 2] < 7.5)] <- "640-Ig-4h"

seuratObj_demultiplex@meta.data$MULTI_ID <-  factor(seuratObj_demultiplex@meta.data$MULTI_ID, levels = c(samples,
                                                                                                         "Negative",
                                                                                                         "Doublet"))



negative_p1 <- seuratObj_demultiplex@meta.data %>%
  filter(MULTI_classification == "Negative") %>%
  select(MULTI_classification, status, assignment) %>%
  mutate(assignment = case_when(
    assignment == "1" ~ "640",
    assignment == "0" ~ "678",
    assignment == "0/1" ~ "Doublet",
    assignment == "1/0" ~ "Doublet",
  ))%>%
  group_by(MULTI_classification, assignment) %>%
  summarize(count = n(), .groups = "drop") %>%
  ggplot(aes(x = MULTI_classification, y = count, fill = assignment)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#ef476f", "#ffd166", "#118ab2")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(title = "Negatives \nbefore curation",
       x = "Classification",
       y = "Count",
       fill = "Souporcell\nassignment") +
  ylim(0, 1000)


negative_p2 <- seuratObj_demultiplex@meta.data %>%
  filter(MULTI_ID == "Negative") %>%
  select(MULTI_ID, status, assignment) %>%
  mutate(assignment = case_when(
    assignment == "1" ~ "640",
    assignment == "0" ~ "678",
    assignment == "0/1" ~ "Doublet",
    assignment == "1/0" ~ "Doublet",
  ))%>%
  group_by(MULTI_ID, assignment) %>%
  summarize(count = n(), .groups = "drop") %>%
  ggplot(aes(x = MULTI_ID, y = count, fill = assignment)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#ef476f", "#ffd166", "#118ab2")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(title = "Negatives \nafter curation",
       x = "Classification",
       y = "Count",
       fill = "Souporcell\nassignment") +
  ylim(0, 1000)



doublets_before <- seuratObj_demultiplex@meta.data %>%
  filter(grepl("_", MULTI_classification) | is.na(MULTI_classification)) %>%
  select(MULTI_ID, MULTI_classification, status, assignment) %>%
  mutate(assignment = case_when(
    assignment == "1" ~ "640",
    assignment == "0" ~ "678",
    assignment == "0/1" ~ "Doublet",
    assignment == "1/0" ~ "Doublet",
  ))

doublet_p1 <-  doublets_before %>%
  group_by(MULTI_classification, assignment) %>%
  summarize(count = n(), .groups = "drop") %>%
  ggplot(aes(x = MULTI_classification, y = count, fill = assignment)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#ef476f", "#ffd166", "#118ab2")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(title = "Doublets \nbefore curation",
       x = "Classification",
       y = "Count",
       fill = "Souporcell assignment")+
  ylim(0, 400)

doublets_after <- seuratObj_demultiplex@meta.data %>%
  filter(MULTI_ID == "Doublet") %>%
  select(MULTI_ID, MULTI_classification, status, assignment) %>%
  mutate(assignment = case_when(
    assignment == "1" ~ "640",
    assignment == "0" ~ "678",
    assignment == "0/1" ~ "Doublet",
    assignment == "1/0" ~ "Doublet",
  ))

doublet_p2 <- doublets_after %>%
  group_by(MULTI_classification, assignment) %>%
  summarize(count = n(), .groups = "drop") %>%
  ggplot(aes(x = MULTI_classification, y = count, fill = assignment)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#ef476f", "#ffd166", "#118ab2")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(title = "Doublets \nafter curation",
       x = "Classification",
       y = "Count",
       fill = "Souporcell assignment")+
  ylim(0, 400)


pdf(paste0("qc/qc_NS-DD-1s-DEC-2/multiseqdemux/FINAL/negatives_doublets_curation_multi-seq.pdf"), width = 6.07, height = 8.56)
(negative_p1 + doublet_p1 + plot_layout(widths = c(0.2, 1))) /
  (negative_p2 + doublet_p2 + plot_layout(widths = c(0.2, 1))) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", axis.title.x = element_blank())

dev.off()


# 
# 
# 
# DimPlot1 <- DimPlot(seuratObj_demultiplex, group.by = "MULTI_ID") +
#   coord_cartesian(
#     xlim = c(min(seuratObj_demultiplex@reductions$tsne@cell.embeddings[, 1]), 25),
#     ylim = c(0, 40)) +
#   scale_color_manual(values = color_vector) +
#   theme_bw()
# 
# 
# DimPlot2 <- DimPlot(seuratObj_demultiplex, group.by = "MULTI_ID_5") +
#   coord_cartesian(
#     xlim = c(min(seuratObj_demultiplex@reductions$tsne@cell.embeddings[, 1]), 25),
#     ylim = c(0, 40)) +
#   scale_color_manual(values = color_vector) +
#   theme_bw()
# 
# 
# 
# pdf(paste0(output_dir, "multiseqdemux/tsne_cluster1_MULTI_ID_5.pdf"), width = 10.82, height = 4.39)
# plot_grid(DimPlot1, DimPlot2,
#           labels = c("A", "B"),
#           nrow = 1,
#           rel_widths = c(1, 1))
# dev.off()


Idents(seuratObj_demultiplex) <- "MULTI_ID"


seuratObj_demultiplex <- SetIdent(seuratObj_demultiplex,
                                  value = factor(Idents(seuratObj_demultiplex), levels = c(samples,
                                                                                           "Negative",
                                                                                           "Doublet")))

ridge_plots <- lapply(1:length(samples), function(i) {
  RidgePlot(seuratObj_demultiplex, assay = "HTO", features = rownames(seuratObj_demultiplex[["HTO"]])[i], ncol = 1) +
    scale_fill_manual(values = color_vector) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 8)) +
    guides(fill = "none")
})


# Combine the Ridge plots using patchwork
RidgePlot <- wrap_plots(ridge_plots, ncol = 3)

VlnPlot1 <- VlnPlot(seuratObj_demultiplex, features = "nFeature_RNA", pt.size = 0, log = TRUE) +
  geom_jitter(aes(y = nFeature_RNA), color = "grey", size = 0.05, alpha = 1, position = position_jitter(width = 0.3)) +  
  scale_fill_manual(values = color_vector) +
  theme_bw() +
  guides(fill = "none")

VlnPlot2 <- VlnPlot(seuratObj_demultiplex, features = "nCount_RNA", pt.size = 0, log = TRUE) +
  geom_jitter(aes(y = nCount_RNA), color = "grey", size = 0.05, alpha = 1, position = position_jitter(width = 0.3)) +  
  scale_fill_manual(values = color_vector) +
  theme_bw() +
  guides(fill = "none")


# counts of classification
table(seuratObj_demultiplex@meta.data$MULTI_ID) 

DefaultAssay(seuratObj_demultiplex) <- "HTO"
seuratObj_demultiplex <- ScaleData(seuratObj_demultiplex, features = rownames(seuratObj_demultiplex),
                                   verbose = FALSE)
seuratObj_demultiplex <- RunPCA(seuratObj_demultiplex, features = rownames(seuratObj_demultiplex), approx = FALSE)
seuratObj_demultiplex <- RunTSNE(seuratObj_demultiplex, dims = 1:3, check_duplicates = FALSE, perplexity = 100)


seuratObj_demultiplex@meta.data$MULTI_ID <-  factor(seuratObj_demultiplex@meta.data$MULTI_ID, levels = c(samples,
                                                                                                         "Negative",
                                                                                                         "Doublet"))
DimPlot1 <- DimPlot(seuratObj_demultiplex, group.by = "MULTI_ID") +
  scale_color_manual(values = color_vector) +
  theme_bw()



# DimPlot <- (DimPlot1 + DimPlot1) + plot_layout(ncol = 2)


MULTI_ID.df <- as.data.frame(table(seuratObj_demultiplex@meta.data$MULTI_ID))
colnames(MULTI_ID.df) <- c("MULTI_ID", "count")

RidgePlot_rel_height <- 1 + 1.3 * ceiling((length(samples) - 3) / 3) 


output_dir <- paste0("qc/qc_", lib, "/")


pdf(paste0(output_dir, "multiseqdemux/FINAL/final_optimal_pq_per_tag_margin_", margin, ".pdf"), width = 7.6, height = 9.76)


p1 <- plot_grid(RidgePlot, VlnPlot1, VlnPlot2,
                labels = c("A", "B", "C", "D"), 
                nrow = 3,
                rel_heights = c(RidgePlot_rel_height,1,1,1))


p2 <- plot_grid(DimPlot1, tableGrob(MULTI_ID.df),
                labels = c("E", "F"),
                nrow = 1,
                rel_widths = c(1.4, 1))

print(plot_grid(p1, p2,
                nrow = 2,
                rel_heights = c(2,1)))

dev.off()



# Generate FeaturePlots for each quantile
plot_demultiplex <- FeaturePlot(seuratObj_demultiplex, features = rownames(seuratObj_demultiplex[["HTO"]]), 
                                cols = c("lightgrey", "blue"), ncol = 3,  order = T) +
  plot_annotation(title = paste0(lib, " multi-seq"))

if(length(samples) == 3){
  plot_height = 8/3
}else{
  plot_height = 16/3
}

pdf(paste0(output_dir, "multiseqdemux/FINAL/final_hto_exp_per_tag_margin", margin, ".pdf"), width = 8.94, height = plot_height)

plot_grid(plot_demultiplex,
          ncol = 1)

dev.off()



plot_HTO_expression_pairs <- function(seuratObj, margin,
                                      check.joint.bcs) {
  hashtags <- rownames(seuratObj[["HTO"]])
  plots_list <- list()
  
  # Loop over all pairs of hashtags
  for (i in 1:(length(hashtags) - 1)) {
    for (j in (i + 1):length(hashtags)) {
      
      # Extract the expression data for the two hashtags
      hashtag1 <- hashtags[i]
      hashtag2 <- hashtags[j]
      
      exp1 <- seuratObj[["HTO"]]@data[hashtag1, ]
      exp2 <- seuratObj[["HTO"]]@data[hashtag2, ]
      
      # Create a data frame for ggplot
      plot_data <- data.frame(
        Exp1 = exp1,
        Exp2 = exp2,
        MULTIID = seuratObj$MULTI_ID  # Metadata for coloring
      )
      
      # Create the scatter plot for this pair of hashtags
      p <- ggplot(plot_data, aes(x = Exp1, y = Exp2, color = MULTIID)) +
        geom_point(alpha = 1) +
        labs(
          x = paste(hashtag1, "exp"),
          y = paste(hashtag2, "exp"),
          title = paste(hashtag1, "vs", hashtag2)
        ) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # Add 1:1 line
        theme_bw() +
        theme(plot.title = element_text(size = 10)) +
        scale_color_manual(values = color_vector, name = paste("optimal p.q."))  # Updated to use positive.quantile[3]
      
      # Add the plot to the list
      plots_list[[paste(hashtag1, hashtag2, sep = "_vs_")]] <- p
    }
  }
  
  plot_height <- ifelse(length(hashtags) == 3, 2.6, 13)
  
  # Combine all the plots into one figure and only display the legend once
  final_plot <- wrap_plots(plots_list, ncol = 3) +
    plot_layout(guides = "collect")
  
  output_dir <- paste0("qc/qc_", lib, "/")
  
  if (check.joint.bcs) {
    output_dir <- paste0(output_dir, "multiseqdemux_joint_bcs/")
  }else{
    output_dir <- paste0(output_dir, "multiseqdemux/")
  }
  
  # Create the directory structure if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  
  pdf(paste0(output_dir, "FINAL/hto_exp_pairwise_optimal_pq_margin_",  margin, ".pdf"), width = 9, height = plot_height)
  print(final_plot)
  dev.off()
}


plot_HTO_expression_pairs(
  seuratObj = seuratObj_demultiplex,
  margin, check.joint.bcs)















## 1.2 table summarzing counts per treatment + experiment + positive.quantile threshold ----
## htodemux_norm_by_features ----
calculate_ratios <- function(df) {
  # Calculate column sums (ignoring NA values)
  col_sum <- colSums(df[1:(nrow(df)-2), c("cnt.0.9", "cnt.0.95", "cnt.0.99", "exp.cnt")], na.rm = TRUE)

  # Initialize the ratio columns
  df$ratio.0.9 <- NA
  df$ratio.0.95 <- NA
  df$ratio.0.99 <- NA
  df$ratio.exp <- NA

  # Calculate ratios for cnt.0.9, cnt.0.95, cnt.0.99, and exp.cnt
  df$ratio.0.9[1:(nrow(df)-2)] <- round(df$cnt.0.9[1:(nrow(df)-2)] / col_sum[1], 2)
  df$ratio.0.95[1:(nrow(df)-2)] <- round(df$cnt.0.95[1:(nrow(df)-2)] / col_sum[2], 2)
  df$ratio.0.99[1:(nrow(df)-2)] <- round(df$cnt.0.99[1:(nrow(df)-2)] / col_sum[3], 2)
  df$ratio.exp[1:(nrow(df)-2)] <- round(df$exp.cnt[1:(nrow(df)-2)] / col_sum[4], 2)

  return(df)
}

exp_obs_cell_cnts1 <- data.frame(
  hash.label = c("699−unstim−0h", "699−Ig−4h", "699−CD3−4h", "Negative", "Doublet"),
  hash.id = c("0251", "0252", "0253", NA, NA),
  cnt.0.9 = c(8090, 4630, 4350, 664, 1578),
  cnt.0.95 = c(7186, 4671, 4452, 1710, 1293),
  cnt.0.99 = c(4941, 4864, 4590, 4205, 712),
  exp.cnt = c(42000, 22000, 20050, NA, NA)
)

exp_obs_cell_cnts2 <- data.frame(
  hash.label = c("640-unstim-0h", "640-Ig-4h", "640-CD3-4h",
                 "678-unstim-0h", "678-Ig-4h", "678-CD3-4h", "Negative", "Doublet"),
  hash.id = c("0251", "0252", "0253", "0254", "0255", "0256", NA, NA),
  cnt.0.9 = c(1255, 617, 1954, 3675, 2365, 193, 757, 5146),
  cnt.0.95 = c(1421, 698, 2120, 4055, 2573, 218, 828, 4049),
  cnt.0.99 = c(785, 753, 2316, 4367, 2742, 724, 1756, 2519),
  exp.cnt = c(26500, 6700, 22000, 44000, 24000, 22000, NA, NA)
)


exp_obs_cell_cnts3 <- data.frame(
  hash.label = c("548-unstim-0h", "548-Ig-4h", "548-CD3-4h",
                 "569-unstim-0h", "569-Ig-4h", "569-CD3-4h", "Negative", "Doublet"),
  hash.id = c("0251", "0252", "0253", "0254", "0255", "0256", NA, NA),
  cnt.0.9 = c(1403, 1112, 2052, 2034, 2278, 574, 46, 5701),
  cnt.0.95 = c(1880, 1260, 2296, 2241, 2599, 728, 136, 4060),
  cnt.0.99 = c(2707, 1370, 2471, 2167, 2953, 991, 538, 2003),
  exp.cnt = c(22000, 8500, 18400, 22000, 21000, 9550, NA, NA)
)

exp_obs_cell_cnts4 <- data.frame(
  hash.label = c("589-unstim-0h", "589-Ig-4h", "589-CD3-4h",
                 "633-unstim-0h", "633-Ig-4h", "633-CD3-4h", "Negative", "Doublet"),
  hash.id = c("0251", "0252", "0253", "0254", "0255", "0256", NA, NA),
  cnt.0.9 = c(2127, 3290, 2675, 2311, 2526, 2555, 289, 2701),
  cnt.0.95 = c(2494, 3330, 2704, 2315, 2533, 2622, 388, 2088),
  cnt.0.99 = c(2624, 3340, 2686, 2285, 2479, 2654, 648, 1758),
  exp.cnt = c(22000, 24000, 24000, 22000, 24000, 24000, NA, NA)
)


exp_obs_cell_cnts1 <- calculate_ratios(exp_obs_cell_cnts1)
exp_obs_cell_cnts2 <- calculate_ratios(exp_obs_cell_cnts2)
exp_obs_cell_cnts3 <- calculate_ratios(exp_obs_cell_cnts3)
exp_obs_cell_cnts4 <- calculate_ratios(exp_obs_cell_cnts4)

# add additional information (experiment and demultiplexing)
exp_obs_cell_cnts1$lib <- "NS-DD-1s-DEC-1"
exp_obs_cell_cnts2$lib <- "NS-DD-1s-DEC-2"
exp_obs_cell_cnts3$lib <- "NS-DD-1s-DEC-3"
exp_obs_cell_cnts4$lib <- "NS-DD-1s-DEC-4"

exp_obs_cell_cnts1$demux <- "htodemux_norm_by_features"
exp_obs_cell_cnts2$demux <- "htodemux_norm_by_features"
exp_obs_cell_cnts3$demux <- "htodemux_norm_by_features"
exp_obs_cell_cnts4$demux <- "htodemux_norm_by_features"

df1 <- rbind(
  exp_obs_cell_cnts1,
  exp_obs_cell_cnts2,
  exp_obs_cell_cnts3,
  exp_obs_cell_cnts4
)


custom_theme <- ttheme_minimal(
  core = list(bg_params = list(fill = "white", col = "black", lwd = 1.5)), # white cells with black borders
  colhead = list(bg_params = list(fill = "white", col = "black", lwd = 1.5)) # white header with black borders
)

p <- plot_grid(tableGrob(exp_obs_cell_cnts1, theme = custom_theme),
               tableGrob(exp_obs_cell_cnts2, theme = custom_theme),
               tableGrob(exp_obs_cell_cnts3, theme = custom_theme),
               tableGrob(exp_obs_cell_cnts4, theme = custom_theme),
               labels = c("NS-DD-1s-DEC-1",
                          "NS-DD-1s-DEC-2",
                          "NS-DD-1s-DEC-3",
                          "NS-DD-1s-DEC-4"),
               nrow = 4,
               rel_heights = c(0.7,1,1,1))

pdf(paste0("qc/qc_NS-DD-1s-DEC_summary/exp_obs_cell_cnts_htodemux_margin_1_per_tag.pdf"), width = 9, height = 12)
print(p)
dev.off()

## htodemux_norm_by_cells ----
exp_obs_cell_cnts1 <- data.frame(
  hash.label = c("699−unstim−0h", "699−Ig−4h", "699−CD3−4h", "Negative", "Doublet"),
  hash.id = c("0251", "0252", "0253", NA, NA),
  cnt.0.9 = c(7846, 4331, 4175, 70, 2890),
  cnt.0.95 = c(8634, 4599, 4323, 108, 1648),
  cnt.0.99 = c(8714, 4707, 4368, 282, 1241),
  exp.cnt = c(42000, 22000, 20050, NA, NA)
)

exp_obs_cell_cnts2 <- data.frame(
  hash.label = c("640-unstim-0h", "640-Ig-4h", "640-CD3-4h",
                 "678-unstim-0h", "678-Ig-4h", "678-CD3-4h", "Negative", "Doublet"),
  hash.id = c("0251", "0252", "0253", "0254", "0255", "0256", NA, NA),
  cnt.0.9 = c(620, 677, 2741, 4116, 2565, 343, 875, 4025),
  cnt.0.95 = c(648, 727, 2831, 4268, 2673, 360, 939, 3516),
  cnt.0.99 = c(486, 748, 2408, 4414, 2752, 790, 1965, 2399),
  exp.cnt = c(26500, 6700, 22000, 44000, 24000, 22000, NA, NA)
)


exp_obs_cell_cnts3 <- data.frame(
  hash.label = c("548-unstim-0h", "548-Ig-4h", "548-CD3-4h",
                 "569-unstim-0h", "569-Ig-4h", "569-CD3-4h", "Negative", "Doublet"),
  hash.id = c("0251", "0252", "0253", "0254", "0255", "0256", NA, NA),
  cnt.0.9 = c(2135, 1282, 2311, 2351, 2651, 817, 63, 3590),
  cnt.0.95 = c(2718, 1359, 2443, 2468, 2869, 974, 183, 2186),
  cnt.0.99 = c(3048, 1402, 2508, 2155, 2837, 1121, 837, 1292),
  exp.cnt = c(22000, 8500, 18400, 22000, 21000, 9550, NA, NA)
)

exp_obs_cell_cnts4 <- data.frame(
  hash.label = c("589-unstim-0h", "589-Ig-4h", "589-CD3-4h",
                 "633-unstim-0h", "633-Ig-4h", "633-CD3-4h", "Negative", "Doublet"),
  hash.id = c("0251", "0252", "0253", "0254", "0255", "0256", NA, NA),
  cnt.0.9 = c(2004, 3275, 2572, 2286, 2506, 2379, 362, 3090),
  cnt.0.95 = c(2398, 3322, 2557, 2298, 2531, 2535, 534, 2299),
  cnt.0.99 = c(2653, 3354, 2376, 2225, 2454, 2630, 1017, 1765),
  exp.cnt = c(22000, 24000, 24000, 22000, 24000, 24000, NA, NA)
)


exp_obs_cell_cnts1 <- calculate_ratios(exp_obs_cell_cnts1)
exp_obs_cell_cnts2 <- calculate_ratios(exp_obs_cell_cnts2)
exp_obs_cell_cnts3 <- calculate_ratios(exp_obs_cell_cnts3)
exp_obs_cell_cnts4 <- calculate_ratios(exp_obs_cell_cnts4)


# add additional information (experiment and demultiplexing)
exp_obs_cell_cnts1$lib <- "NS-DD-1s-DEC-1"
exp_obs_cell_cnts2$lib <- "NS-DD-1s-DEC-2"
exp_obs_cell_cnts3$lib <- "NS-DD-1s-DEC-3"
exp_obs_cell_cnts4$lib <- "NS-DD-1s-DEC-4"

exp_obs_cell_cnts1$demux <- "htodemux_norm_by_cells"
exp_obs_cell_cnts2$demux <- "htodemux_norm_by_cells"
exp_obs_cell_cnts3$demux <- "htodemux_norm_by_cells"
exp_obs_cell_cnts4$demux <- "htodemux_norm_by_cells"

df2 <- rbind(
  exp_obs_cell_cnts1,
  exp_obs_cell_cnts2,
  exp_obs_cell_cnts3,
  exp_obs_cell_cnts4
)

custom_theme <- ttheme_minimal(
  core = list(bg_params = list(fill = "white", col = "black", lwd = 1.5)), # white cells with black borders
  colhead = list(bg_params = list(fill = "white", col = "black", lwd = 1.5)) # white header with black borders
)

p <- plot_grid(tableGrob(exp_obs_cell_cnts1, theme = custom_theme),
               tableGrob(exp_obs_cell_cnts2, theme = custom_theme),
               tableGrob(exp_obs_cell_cnts3, theme = custom_theme),
               tableGrob(exp_obs_cell_cnts4, theme = custom_theme),
               labels = c("NS-DD-1s-DEC-1",
                          "NS-DD-1s-DEC-2",
                          "NS-DD-1s-DEC-3",
                          "NS-DD-1s-DEC-4"),
               nrow = 4,
               rel_heights = c(0.7,1,1,1))

pdf(paste0("qc/qc_NS-DD-1s-DEC_summary/exp_obs_cell_cnts_htodemux_margin_2_per_tag.pdf"), width = 9, height = 12)
print(p)
dev.off()

## multiseqdemux ----
calculate_ratios_multiseqdemux <- function(df) {
  # Calculate column sums (ignoring NA values)
  col_sum <- colSums(df[1:(nrow(df)-2), c("cnt.optimal", "exp.cnt")], na.rm = TRUE)

  # Initialize the ratio columns
  df$ratio.optimal <- NA  
  df$ratio.exp <- NA
  
  # Calculate ratios for cnt.optimal, and exp.cnt
  df$ratio.optimal[1:(nrow(df)-2)] <- round(df$cnt.optimal[1:(nrow(df)-2)] / col_sum[1], 2)
  df$ratio.exp[1:(nrow(df)-2)] <- round(df$exp.cnt[1:(nrow(df)-2)] / col_sum[2], 2)
  
  return(df)
}


exp_obs_cell_cnts1 <- data.frame(
  hash.label = c("699−unstim−0h", "699−Ig−4h", "699−CD3−4h", "Negative", "Doublet"),
  hash.id = c("0251", "0252", "0253", NA, NA),
  cnt.optimal = c(9013, 4691, 4433, 142, 1033),
  exp.cnt = c(42000, 22000, 20050, NA, NA)
)

exp_obs_cell_cnts2 <- data.frame(
  hash.label = c("640-unstim-0h", "640-Ig-4h", "640-CD3-4h",
                 "678-unstim-0h", "678-Ig-4h", "678-CD3-4h", "Negative", "Doublet"),
  hash.id = c("0251", "0252", "0253", "0254", "0255", "0256", NA, NA),
  cnt.optimal = c(1698, 739, 2114, 4278, 2713, 176, 970, 3274),
  exp.cnt = c(26500, 6700, 22000, 44000, 24000, 22000, NA, NA)
)


exp_obs_cell_cnts3 <- data.frame(
  hash.label = c("548-unstim-0h", "548-Ig-4h", "548-CD3-4h",
                 "569-unstim-0h", "569-Ig-4h", "569-CD3-4h", "Negative", "Doublet"),
  hash.id = c("0251", "0252", "0253", "0254", "0255", "0256", NA, NA),
  cnt.optimal = c(3072, 1379, 2511, 2669, 3078, 1246, 300, 945),
  exp.cnt = c(22000, 8500, 18400, 22000, 21000, 9550, NA, NA)
)

exp_obs_cell_cnts4 <- data.frame(
  hash.label = c("589-unstim-0h", "589-Ig-4h", "589-CD3-4h",
                 "633-unstim-0h", "633-Ig-4h", "633-CD3-4h", "Negative", "Doublet"),
  hash.id = c("0251", "0252", "0253", "0254", "0255", "0256", NA, NA),
  cnt.optimal = c(2575, 3320, 2710, 2336, 2556, 2642, 451, 1884),
  exp.cnt = c(22000, 24000, 24000, 22000, 24000, 24000, NA, NA)
)


exp_obs_cell_cnts1 <- calculate_ratios_multiseqdemux(exp_obs_cell_cnts1)
exp_obs_cell_cnts2 <- calculate_ratios_multiseqdemux(exp_obs_cell_cnts2)
exp_obs_cell_cnts3 <- calculate_ratios_multiseqdemux(exp_obs_cell_cnts3)
exp_obs_cell_cnts4 <- calculate_ratios_multiseqdemux(exp_obs_cell_cnts4)


# add additional information (experiment and demultiplexing)
exp_obs_cell_cnts1$lib <- "NS-DD-1s-DEC-1" 
exp_obs_cell_cnts2$lib <- "NS-DD-1s-DEC-2" 
exp_obs_cell_cnts3$lib <- "NS-DD-1s-DEC-3" 
exp_obs_cell_cnts4$lib <- "NS-DD-1s-DEC-4" 

exp_obs_cell_cnts1$demux <- "multiseqdemux" 
exp_obs_cell_cnts2$demux <- "multiseqdemux" 
exp_obs_cell_cnts3$demux <- "multiseqdemux" 
exp_obs_cell_cnts4$demux <- "multiseqdemux" 

df3 <- rbind(
  exp_obs_cell_cnts1,
  exp_obs_cell_cnts2,
  exp_obs_cell_cnts3,
  exp_obs_cell_cnts4
)

custom_theme <- ttheme_minimal(
  core = list(bg_params = list(fill = "white", col = "black", lwd = 1.5)), # white cells with black borders
  colhead = list(bg_params = list(fill = "white", col = "black", lwd = 1.5)) # white header with black borders
)

p <- plot_grid(tableGrob(exp_obs_cell_cnts1, theme = custom_theme),
               tableGrob(exp_obs_cell_cnts2, theme = custom_theme),
               tableGrob(exp_obs_cell_cnts3, theme = custom_theme),
               tableGrob(exp_obs_cell_cnts4, theme = custom_theme),
               labels = c("NS-DD-1s-DEC-1", 
                          "NS-DD-1s-DEC-2",
                          "NS-DD-1s-DEC-3",
                          "NS-DD-1s-DEC-4"),
               nrow = 4,
               rel_heights = c(0.7,1,1,1))

pdf(paste0("qc/qc_NS-DD-1s-DEC_summary/exp_obs_cell_cnts_multiseqdemux_per_tag.pdf"), width = 9, height = 12)
print(p)
dev.off()

##summary----
df <- rbind(df1, df2)

df_long <- df %>%
  gather(key = "cnt_type", value = "cnt", cnt.0.9, cnt.0.95, cnt.0.99) %>%
  gather(key = "ratio_type", value = "ratio", ratio.0.9, ratio.0.95, ratio.0.99) %>%
  filter(sub("cnt.", "", cnt_type) == sub("ratio.", "", ratio_type)) %>%  # Match cnt with ratio type
  mutate(
    demux = paste0(demux, "_", sub("cnt.", "", cnt_type)),  # Create demux value
    cnt_type = NULL,  # Drop cnt_type column
    ratio_type = NULL  # Drop ratio_type column
  )


df3_long <- df3 %>%
  gather(key = "cnt_type", value = "cnt", cnt.optimal) %>%
  gather(key = "ratio_type", value = "ratio", ratio.optimal) %>%
  mutate(
    cnt_type = NULL,  # Drop cnt_type column
    ratio_type = NULL  # Drop ratio_type column
  )

exp_obs_cell_cnts <- rbind(df_long, df3_long)


p <- ggplot(exp_obs_cell_cnts, aes(x = hash.label, y = cnt, color = demux, group = interaction(demux, lib))) +
  geom_point(size = 3) +
  geom_line(size = 0.8) +
  facet_wrap(~lib, scales = "free_x") +
  labs(
    title = "Counts by hashtags for each experiment",
    x = "Hashtags",
    y = "Counts",
    color = "Demux"
  ) +
  # scale_color_manual(
  #   values = brewer.pal(n = length(unique(exp_obs_cell_cnts$demux)), "Set2")
  # ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "right",
    strip.background = element_rect(fill = "white", color = "black"),
  )

pdf(paste0("qc/qc_NS-DD-1s-DEC_summary/exp_obs_cell_cnts_comparison.pdf"), width = 10.22, height = 7.42)
print(p)
dev.off()






## FINAL multiseqdemux ----
calculate_ratios_multiseqdemux <- function(df) {
  # Calculate column sums (ignoring NA values)
  col_sum <- colSums(df[1:(nrow(df)-2), c("cnt.optimal", "exp.cnt")], na.rm = TRUE)
  
  # Initialize the ratio columns
  df$ratio.optimal <- NA  
  df$ratio.exp <- NA
  
  # Calculate ratios for cnt.optimal, and exp.cnt
  df$ratio.optimal[1:(nrow(df)-2)] <- round(df$cnt.optimal[1:(nrow(df)-2)] / col_sum[1], 2)
  df$ratio.exp[1:(nrow(df)-2)] <- round(df$exp.cnt[1:(nrow(df)-2)] / col_sum[2], 2)
  
  return(df)
}


exp_obs_cell_cnts1 <- data.frame(
  hash.label = c("699−unstim−0h", "699−Ig−4h", "699−CD3−4h", "Negative", "Doublet"),
  hash.id = c("0251", "0252", "0253", NA, NA),
  cnt.optimal = c(9013, 4691, 4433, 142, 1033),
  exp.cnt = c(42000, 22000, 20050, NA, NA)
)

exp_obs_cell_cnts2 <- data.frame(
  hash.label = c("640-unstim-0h", "640-Ig-4h", "640-CD3-4h",
                 "678-unstim-0h", "678-Ig-4h", "678-CD3-4h", "Negative", "Doublet"),
  hash.id = c("0251", "0252", "0253", "0254", "0255", "0256", NA, NA),
  cnt.optimal = c(1698, 1497, 2114, 4278, 2713, 2455, 212, 995),
  exp.cnt = c(26500, 6700, 22000, 44000, 24000, 22000, NA, NA)
)


exp_obs_cell_cnts3 <- data.frame(
  hash.label = c("548-unstim-0h", "548-Ig-4h", "548-CD3-4h",
                 "569-unstim-0h", "569-Ig-4h", "569-CD3-4h", "Negative", "Doublet"),
  hash.id = c("0251", "0252", "0253", "0254", "0255", "0256", NA, NA),
  cnt.optimal = c(3072, 1379, 2511, 2669, 3078, 1246, 300, 945),
  exp.cnt = c(22000, 8500, 18400, 22000, 21000, 9550, NA, NA)
)

exp_obs_cell_cnts4 <- data.frame(
  hash.label = c("589-unstim-0h", "589-Ig-4h", "589-CD3-4h",
                 "633-unstim-0h", "633-Ig-4h", "633-CD3-4h", "Negative", "Doublet"),
  hash.id = c("0251", "0252", "0253", "0254", "0255", "0256", NA, NA),
  cnt.optimal = c(2575, 3320, 2710, 2336, 2556, 2642, 451, 1884),
  exp.cnt = c(22000, 24000, 24000, 22000, 24000, 24000, NA, NA)
)


exp_obs_cell_cnts1 <- calculate_ratios_multiseqdemux(exp_obs_cell_cnts1)
exp_obs_cell_cnts2 <- calculate_ratios_multiseqdemux(exp_obs_cell_cnts2)
exp_obs_cell_cnts3 <- calculate_ratios_multiseqdemux(exp_obs_cell_cnts3)
exp_obs_cell_cnts4 <- calculate_ratios_multiseqdemux(exp_obs_cell_cnts4)


# add additional information (experiment and demultiplexing)
exp_obs_cell_cnts1$lib <- "NS-DD-1s-DEC-1" 
exp_obs_cell_cnts2$lib <- "NS-DD-1s-DEC-2" 
exp_obs_cell_cnts3$lib <- "NS-DD-1s-DEC-3" 
exp_obs_cell_cnts4$lib <- "NS-DD-1s-DEC-4" 

exp_obs_cell_cnts1$demux <- "multiseqdemux" 
exp_obs_cell_cnts2$demux <- "multiseqdemux" 
exp_obs_cell_cnts3$demux <- "multiseqdemux" 
exp_obs_cell_cnts4$demux <- "multiseqdemux" 

df3 <- rbind(
  exp_obs_cell_cnts1,
  exp_obs_cell_cnts2,
  exp_obs_cell_cnts3,
  exp_obs_cell_cnts4
)

custom_theme <- ttheme_minimal(
  core = list(bg_params = list(fill = "white", col = "black", lwd = 1.5)), # white cells with black borders
  colhead = list(bg_params = list(fill = "white", col = "black", lwd = 1.5)) # white header with black borders
)

p <- plot_grid(tableGrob(exp_obs_cell_cnts1, theme = custom_theme),
               tableGrob(exp_obs_cell_cnts2, theme = custom_theme),
               tableGrob(exp_obs_cell_cnts3, theme = custom_theme),
               tableGrob(exp_obs_cell_cnts4, theme = custom_theme),
               labels = c("NS-DD-1s-DEC-1", 
                          "NS-DD-1s-DEC-2",
                          "NS-DD-1s-DEC-3",
                          "NS-DD-1s-DEC-4"),
               nrow = 4,
               rel_heights = c(0.7,1,1,1))

pdf(paste0("qc/qc_NS-DD-1s-DEC_summary_FINAL/exp_obs_cell_cnts_multiseqdemux_per_tag.pdf"), width = 9, height = 12)
print(p)
dev.off()


######
######
## 1.3 which threshold to use  -----------------------

# if (lib == "NS-DD-1s-DEC-3"){
#   seuratObj_singlets <- plot_HTO_tag_counts(seuratObj_multiplex, positive.quantile, 3)
# }else if(lib == "NS-DD-1s-DEC-1"){
#   seuratObj_singlets <- plot_HTO_tag_counts(seuratObj_multiplex, positive.quantile, 1)
# }

# seuratObj_tags <- list()
# 
# # Loop through each unique hash.ID and subset the Seurat object
# for (id in samples) {
#   seuratObj_tag <- subset(seuratObj_singlets, subset = hash.ID == id)
#   seuratObj_tags[[id]] <- seuratObj_tag
#   saveRDS(seuratObj_tag, file = paste0(resDir, "/", lib, "/count/sample_filtered_feature_bc_matrix_", id, ".rds"))
#   
# }


#--
# go with multiseqdeux
# run till seuratObj_demultiplex <- plot_HTO_tag_counts(seuratObj_multiplex, margin, check.joint.bcs)


seuratObj_tags <- list()

# Loop through each unique MULTI_ID and subset the Seurat object
for (id in samples) {
  seuratObj_tag <- subset(seuratObj_demultiplex, subset = MULTI_ID == id)
  seuratObj_tags[[id]] <- seuratObj_tag
  saveRDS(seuratObj_tag, file = paste0(resDir, "/", lib, "/count/multiseqdemux_sample_filtered_feature_bc_matrix_", id, ".rds"))
  
}


# spot check
seuratObj_tags[["678-CD3-4h"]]@meta.data %>%
  dim()





