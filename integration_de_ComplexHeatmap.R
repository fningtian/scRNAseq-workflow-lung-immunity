rm(list = ls())


library(Seurat)
library(readxl) 
library(writexl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(patchwork)
library(data.table)
library(tibble) # rownames_to_column
library(tidyr) # pivot_longer
library(stringr) # str_replace_all

# some packages
require(RColorBrewer)
require(ComplexHeatmap)
require(circlize)
require(digest)
require(cluster)

excel.dir <- "integration_2/integration_2_leiden/results_wOrgClusterAnnotation_DEGs/"

rds <- readRDS('integration_2/integration_2_leiden/RDS_Dir/integration_2_leiden_annot.rds')
names(table(Idents(rds)))
# de.dir <- "stimuli_time_MAST/"
# the complete path is this, results_wOrgClusterAnnotation_DEGs/stimuli_time_MAST/
# that contains de outputs and the figures below 


# Read GTF file as tab-delimited text
gtf_path <- "integration_2/integration_2_leiden/refdata-gex-GRCh38-2024-A/genes/protein_coding_genes.gtf"
gtf <- fread(gtf_path, sep = "\t", header = FALSE)

colnames(gtf) <- c("seqname", "source", "feature", "start", "end", 
                        "score", "strand", "frame", "attribute")

gtf[, gene_id := sub('.*gene_id "([^"]+)".*', '\\1', attribute)]
gtf[, gene_name := sub('.*gene_name "([^"]+)".*', '\\1', attribute)]
gtf[, gene_type := sub('.*gene_type "([^"]+)".*', '\\1', attribute)]

# gtf <- gtf[feature == "gene" & gtf == "protein_coding"]

pcgs <- unique(gtf$gene_name)


############################### COMPLEX HEATMAP ###############################
##-- skip cellType CD4 T, CD8 T, NK, Mono/Mph | unstim_0h CD3_4h CD3_18h | 5 donors | 80 top de genes ----
de.dir <- "stimuli_time_MAST/"
for (cellType in c("NK cells", "CD4 T cells", "CD8 T cells")){
  print(cellType)
  rds.sub <- subset(rds, idents = cellType)
  rds.sub <- subset(rds.sub, subset = expCond.stimuli %in% c("unstim", "CD3"))
  
  file.all <- list.files(paste0(excel.dir, de.dir))
  file <- file.all[grepl("CD3.*top", file.all, ignore.case = TRUE)]
  # excel_data <- lapply(paste0(excel.dir, de.dir, file), read_excel)
  
  expression <- function(filename) {
    # Remove the .xlsx extension
    filename <- sub("\\.xlsx$", "", filename)
    filename <- sub("expCondCompDeMarkers_", "", filename)
    filename <- sub("-unstim_0h", "", filename)
    return(filename)
  }
  
  # Read each Excel file and add a new column with the last two underscore-separated parts of the filename
  excel_data <- lapply(file, function(f) {
    # Read the Excel file
    df <- read_excel(paste0(excel.dir, de.dir, f))
    
    # Extract the last two parts of the filename
    expression <- expression(f)
    
    # Add a new column with the extracted filename part
    df$expression <- expression
    
    return(df)
  })
  
  
  combined_data <- do.call(rbind, excel_data) %>% 
    filter(geneType %in% cellType)
  
  geneNames <- c(combined_data$Gene)
  features <- geneNames
  resDir <- paste0(excel.dir, de.dir)
  seuratObjFinal <- rds.sub
  expCondReorderLevels <- c("602_unstim_0h", "602_CD3_4h", "602_CD3_18h",
                            "616_unstim_0h", "616_CD3_4h", "616_CD3_18h",
                            "662_unstim_0h", "662_CD3_4h", "662_CD3_18h",
                            "673_unstim_0h", "673_CD3_4h", "673_CD3_18h",
                            "681_unstim_0h", "681_CD3_4h", "681_CD3_18h")
  slot <- "scale.data"
  Seurat::DefaultAssay(seuratObjFinal) <- "RNA"
  
  # by default, both do.scale/center are on, 
  # will scale the expression level for each gene feature by dividing the centered feature expression levels by their standard deviations
  seuratObjFinal <- Seurat::ScaleData(object = seuratObjFinal, do.scale = T, do.center = T, features = features) 
  
  # 1. input for heatmap
  data <- as.data.frame(as.matrix(t(SeuratObject::GetAssayData(object = seuratObjFinal, 
                                                               layer = slot)[features, colnames(x = seuratObjFinal), drop = FALSE]))) 
  
  metadata <- seuratObjFinal@meta.data[rownames(data), ]
  
  # focus on genes of interest
  data <- t(data[,geneNames])
  data[data < -2.5] <- -2.5
  data[data > 2.5] <- 2.5
  
  
  ## set the color scheme
  # myCol <- colorRampPalette(c("#31d3ff", "black", "#fff200"))(100)
  # myCol <- colorRampPalette(c(rev(brewer.pal(3, "Blues")), "white", brewer.pal(3, "Reds")))(100)
  # myCol <- colorRampPalette(c("blue", "white", "red"))(100)
  myCol <- colorRampPalette(c("#3B4CC0", "#FFFFFF", "#B40426"))(100)
  # myCol <- colorRampPalette(c("#3B4CC0", "#89A0D0", "#FFFFFF", "#D08080", "#B40426"))(100)
  
  # myCol <- colorRampPalette(brewer.pal(5, "RdYlBu"))(100)
  
  myBreaks <- seq(-2.5, 2.5, length.out = 100)
  
  col_ann <- data.frame(
    Donor = factor(metadata$expCond.donor,levels=c('602', '616', '662', '673', '681')),
    Stimuli.Time = factor(metadata$expCond.stimuli.time,
                          levels=c('unstim_0h',
                                   'CD3_4h',
                                   'CD3_18h')),
    stringsAsFactors = FALSE) #do not have to set factor 
  
  row_ann <- data.frame(
    # Genes = rownames(data),
    Stimuli.Time = rev(c(rep('CD3_4h',20), rep('CD3_4h',20),
                         rep('CD3_18h',20), rep('CD3_18h',20))),
    Expression = rev(c(rep('Top20_upDe',20), rep('Top20_downDe',20),
                       rep('Top20_upDe',20), rep('Top20_downDe',20))),
    stringsAsFactors = FALSE)
  
  col_colours <- list(
    Donor = c("602" = "#f2a65a",
              "616" = "#b5c2b7", 
              "662" = "#62466b", 
              "673" = "#9ed8db", 
              "681" = "#394787"),
    Stimuli.Time = c('unstim_0h' = '#eff9f0', 
                     'CD3_4h' = '#FF88DC', 
                     'CD3_18h' = '#91A6FF'))
  
  row_colours <- list(
    Expression = c('Top20_downDe' = '#a8dadc',
                   'Top20_upDe' = '#e63946',
                   'Top20_downDe' = '#a8dadc',
                   'Top20_upDe' = '#e63946'),
    Stimuli.Time = c('CD3_4h' = '#FF88DC', 
                     'CD3_18h' = '#91A6FF'))
  
  
  colAnn <- HeatmapAnnotation(
    df = col_ann,
    which = 'col', # 'col' (samples) or 'row' (gene) annotation?
    na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
    col = col_colours,
    annotation_height = 0.6,
    annotation_width = unit(1, 'cm'),
    simple_anno_size = unit(3, "mm"), # the height of the bar 
    gap = unit(1, 'mm'),
    annotation_legend_param = list(
      Donor = list(
        nrow = 2, # number of rows across which the legend will be arranged
        title = 'Donor',
        title_position = 'topcenter',
        legend_direction = 'vertical',
        title_gp = gpar(fontsize = 11),  # legend on the side
        labels_gp = gpar(fontsize = 11)),  # legend on the side
      Stimuli.Time = list(
        nrow = 3,
        title = 'Stimuli\nTime',
        title_position = 'topcenter',
        legend_direction = 'vertical',
        title_gp = gpar(fontsize = 11), # legend on the side
        labels_gp = gpar(fontsize = 11)))) # legend on the side
  
  rowAnn <- HeatmapAnnotation(
    df = row_ann,
    which = 'row', # 'col' (samples) or 'row' (gene) annotation?
    na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
    col = row_colours,
    annotation_height = 0.6,
    annotation_width = unit(1, 'cm'),
    simple_anno_size = unit(3, "mm"),
    gap = unit(1, 'mm'),
    annotation_legend_param = list(
      Stimuli.Time = list(
        ncol = 4, 
        title = 'Stimuli.Time',
        title_position = 'topcenter',
        legend_direction = 'vertical',
        title_gp = gpar(fontsize = 11),   # row annot legend title (side) 
        labels_gp = gpar(fontsize = 11)),  # row annot legend label (side)
      Expression = list(
        ncol = 4, 
        title = 'Expression',
        title_position = 'topcenter',
        legend_direction = 'vertical',
        title_gp = gpar(fontsize = 11),  # row annot legend title (side)
        labels_gp = gpar(fontsize = 11))  # row annot legend label (side)
    ))
  
  
  
  hmap <- Heatmap(data,
                  # split the genes / rows according to the PAM clusters
                  column_split = factor(metadata$expCond.donor.stimuli.time,
                                        levels = expCondReorderLevels),
                  row_split = factor(c(rep('CD3_18h_top20_downDe',20), rep('CD3_18h_top20_upDe',20),
                                       rep('CD3_4h_top20_downDe',20), rep('CD3_4h_top20_upDe',20)),
                                     levels=c('CD3_4h_top20_upDe','CD3_4h_top20_downDe',
                                              'CD3_18h_top20_upDe','CD3_18h_top20_downDe')),
                  
                  # cluster_row_slices = FALSE,
                  name = 'Scaled\ngene\nexpression',
                  
                  col = colorRamp2(myBreaks, myCol),
                  
                  # parameters for the colour-bar that represents gradient of expression
                  heatmap_legend_param = list(
                    at = c(-2.5, -2, -1, 0, 1, 2, 2.5),
                    color_bar = 'continuous',
                    legend_direction = 'vertical',
                    legend_width = unit(6, 'cm'),
                    legend_height = unit(4, 'cm'),
                    title_position = 'topcenter',
                    title_gp=gpar(fontsize = 11), # heatmap y axis title 
                    labels_gp=gpar(fontsize = 11)), # heatmap gradient legend
                  
                  # row (gene) parameters
                  cluster_rows = FALSE,
                  show_row_dend = FALSE,
                  row_title = 'Differentially expressed genes',
                  row_title_side = 'left',
                  row_title_gp = gpar(fontsize = 11),
                  row_title_rot = 90,
                  show_row_names = TRUE,
                  row_names_gp = gpar(fontsize = 8), # y-axis ticks
                  row_names_side = 'left',
                  # row_dend_width = unit(25,'mm'),
                  
                  # column (sample) parameters
                  cluster_columns = FALSE,
                  show_column_dend = FALSE,
                  # column_title = '',
                  column_title_side = 'top',
                  column_title_gp = gpar(fontsize = 0), # column ticks
                  column_title_rot = 90,
                  show_column_names = FALSE,
                  column_names_gp = gpar(fontsize = 9),
                  column_names_max_height = unit(10, 'cm'),
                  # column_dend_height = unit(25,'mm')
                  
                  # specify top and bottom annotations
                  top_annotation = colAnn,
                  left_annotation = rowAnn)
  
  pdf(paste0(excel.dir, de.dir, cellType, "_complexheatmap.pdf"), height=10, width=12)
  draw(hmap, 
       # + genelabels,
       heatmap_legend_side = 'right',
       annotation_legend_side = 'right')
  dev.off()
}








##-- skip cellType Mono/Mph | unstim_0h Ig_4h Ig_18h | 5 donors | 80 top de genes ----
de.dir <- "stimuli_time_MAST/"
cellTypes <- c("Mono/Mph")

for (cellType in cellTypes){
  print(cellType)
  rds.sub <- subset(rds, idents = cellType)
  rds.sub <- subset(rds.sub, subset = expCond.stimuli %in% c("unstim", "Ig"))
  
  file.all <- list.files(paste0(excel.dir, de.dir))
  file <- file.all[grepl("Ig.*top", file.all, ignore.case = TRUE)]
  # excel_data <- lapply(paste0(excel.dir, de.dir, file), read_excel)
  
  expression <- function(filename) {
    # Remove the .xlsx extension
    filename <- sub("\\.xlsx$", "", filename)
    filename <- sub("expCondCompDeMarkers_", "", filename)
    filename <- sub("-unstim_0h", "", filename)
    return(filename)
  }
  
  # Read each Excel file and add a new column with the last two underscore-separated parts of the filename
  excel_data <- lapply(file, function(f) {
    # Read the Excel file
    df <- read_excel(paste0(excel.dir, de.dir, f))
    
    # Extract the last two parts of the filename
    expression <- expression(f)
    
    # Add a new column with the extracted filename part
    df$expression <- expression
    
    return(df)
  })
  
  
  combined_data <- do.call(rbind, excel_data) %>% 
    filter(geneType %in% cellType)
  
  geneNames <- c(combined_data$Gene)
  features <- geneNames
  resDir <- paste0(excel.dir, de.dir)
  seuratObjFinal <- rds.sub
  
  # Use paste() and expand.grid() to create all combinations
  expCondReorderLevels <- apply(expand.grid(c("unstim_0h", "Ig_4h", "Ig_18h"), 
                                            sort(unique(seuratObjFinal@meta.data$expCond.donor))), 1, function(x) paste(x[2], x[1], sep = "_"))
  # 
  # expCondReorderLevels <- c("602_unstim_0h", "602_Ig_4h", "602_Ig_18h",
  #                           "616_unstim_0h", "616_Ig_4h", "616_Ig_18h",
  #                           "662_unstim_0h", "662_Ig_4h", "662_Ig_18h",
  #                           "673_unstim_0h", "673_Ig_4h", "673_Ig_18h",
  #                           "681_unstim_0h", "681_Ig_4h", "681_Ig_18h")
  slot <- "scale.data"
  Seurat::DefaultAssay(seuratObjFinal) <- "RNA"
  
  # by default, both do.scale/center are on, 
  # will scale the expression level for each gene feature by dividing the centered feature expression levels by their standard deviations
  seuratObjFinal <- Seurat::ScaleData(object = seuratObjFinal, do.scale = T, do.center = T, features = features) 
  
  # 1. input for heatmap
  data <- as.data.frame(as.matrix(t(SeuratObject::GetAssayData(object = seuratObjFinal, 
                                                               layer = slot)[features, colnames(x = seuratObjFinal), drop = FALSE]))) 
  
  metadata <- seuratObjFinal@meta.data[rownames(data), ]
  
  # focus on genes of interest
  data <- t(data[,geneNames])
  data[data < -2.5] <- -2.5
  data[data > 2.5] <- 2.5
  
  
  ## set the color scheme
  # myCol <- colorRampPalette(c("#31d3ff", "black", "#fff200"))(100)
  # myCol <- colorRampPalette(c(rev(brewer.pal(3, "Blues")), "white", brewer.pal(3, "Reds")))(100)
  # myCol <- colorRampPalette(c("blue", "white", "red"))(100)
  myCol <- colorRampPalette(c("#3B4CC0", "#FFFFFF", "#B40426"))(100)
  # myCol <- colorRampPalette(c("#3B4CC0", "#89A0D0", "#FFFFFF", "#D08080", "#B40426"))(100)
  
  # myCol <- colorRampPalette(brewer.pal(5, "RdYlBu"))(100)
  
  myBreaks <- seq(-2.5, 2.5, length.out = 100)
  
  col_ann <- data.frame(
    Donor = factor(metadata$expCond.donor,levels=sort(unique(metadata$expCond.donor))),
    Stimuli.Time = factor(metadata$expCond.stimuli.time,
                          levels=c('unstim_0h',
                                   'Ig_4h',
                                   'Ig_18h')),
    stringsAsFactors = FALSE) #do not have to set factor 
  
  row_ann <- data.frame(
    # Genes = rownames(data),
    Stimuli.Time = rev(c(rep('Ig_4h',20), rep('Ig_4h',20),
                         rep('Ig_18h',20), rep('Ig_18h',20))),
    Expression = rev(c(rep('Top20_upDe',20), rep('Top20_downDe',20),
                       rep('Top20_upDe',20), rep('Top20_downDe',20))),
    stringsAsFactors = FALSE)
  
  col_colours <- list(
    Donor = c("602" = "#f2a65a",
              "616" = "#b5c2b7", 
              "662" = "#62466b", 
              "673" = "#9ed8db", 
              "681" = "#394787"),
    Stimuli.Time = c('unstim_0h' = '#eff9f0', 
                     'Ig_4h' = '#FFBCB5', 
                     'Ig_18h' = '#03978E'))
  
  row_colours <- list(
    Expression = c('Top20_downDe' = '#a8dadc',
                   'Top20_upDe' = '#e63946',
                   'Top20_downDe' = '#a8dadc',
                   'Top20_upDe' = '#e63946'),
    Stimuli.Time = c('Ig_4h' = '#FFBCB5', 
                     'Ig_18h' = '#03978E'))
  
  
  colAnn <- HeatmapAnnotation(
    df = col_ann,
    which = 'col', # 'col' (samples) or 'row' (gene) annotation?
    na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
    col = col_colours,
    annotation_height = 0.6,
    annotation_width = unit(1, 'cm'),
    simple_anno_size = unit(3, "mm"), # the height of the bar 
    gap = unit(1, 'mm'),
    annotation_legend_param = list(
      Donor = list(
        nrow = 4, # number of rows across which the legend will be arranged
        title = 'Donor',
        title_position = 'topcenter',
        legend_direction = 'vertical',
        title_gp = gpar(fontsize = 11),  # legend on the side
        labels_gp = gpar(fontsize = 11)),  # legend on the side
      Stimuli.Time = list(
        nrow = 3,
        title = 'Stimuli\nTime',
        title_position = 'topcenter',
        legend_direction = 'vertical',
        title_gp = gpar(fontsize = 11), # legend on the side
        labels_gp = gpar(fontsize = 11)))) # legend on the side
  
  rowAnn <- HeatmapAnnotation(
    df = row_ann,
    which = 'row', # 'col' (samples) or 'row' (gene) annotation?
    na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
    col = row_colours,
    annotation_height = 0.6,
    annotation_width = unit(1, 'cm'),
    simple_anno_size = unit(3, "mm"),
    gap = unit(1, 'mm'),
    annotation_legend_param = list(
      Stimuli.Time = list(
        ncol = 4, 
        title = 'Stimuli.Time',
        title_position = 'topcenter',
        legend_direction = 'vertical',
        title_gp = gpar(fontsize = 11),   # row annot legend title (side) 
        labels_gp = gpar(fontsize = 11)),  # row annot legend label (side)
      Expression = list(
        ncol = 4, 
        title = 'Expression',
        title_position = 'topcenter',
        legend_direction = 'vertical',
        title_gp = gpar(fontsize = 11),  # row annot legend title (side)
        labels_gp = gpar(fontsize = 11))  # row annot legend label (side)
    ))
  
  
  
  hmap <- Heatmap(data,
                  # split the genes / rows according to the PAM clusters
                  column_split = factor(metadata$expCond.donor.stimuli.time,
                                        levels = expCondReorderLevels),
                  row_split = factor(c(rep('Ig_18h_top20_downDe',20), rep('Ig_18h_top20_upDe',20),
                                       rep('Ig_4h_top20_downDe',20), rep('Ig_4h_top20_upDe',20)),
                                     levels=c('Ig_4h_top20_upDe','Ig_4h_top20_downDe',
                                              'Ig_18h_top20_upDe','Ig_18h_top20_downDe')),
                  
                  # cluster_row_slices = FALSE,
                  name = 'Scaled\ngene\nexpression',
                  
                  col = colorRamp2(myBreaks, myCol),
                  
                  # parameters for the colour-bar that represents gradient of expression
                  heatmap_legend_param = list(
                    at = c(-2.5, -2, -1, 0, 1, 2, 2.5),
                    color_bar = 'continuous',
                    legend_direction = 'vertical',
                    legend_width = unit(6, 'cm'),
                    legend_height = unit(4, 'cm'),
                    title_position = 'topcenter',
                    title_gp=gpar(fontsize = 11), # heatmap y axis title 
                    labels_gp=gpar(fontsize = 11)), # heatmap gradient legend
                  
                  # row (gene) parameters
                  cluster_rows = FALSE,
                  show_row_dend = FALSE,
                  row_title = 'Differentially expressed genes',
                  row_title_side = 'left',
                  row_title_gp = gpar(fontsize = 11),
                  row_title_rot = 90,
                  show_row_names = TRUE,
                  row_names_gp = gpar(fontsize = 8), # y-axis ticks
                  row_names_side = 'left',
                  # row_dend_width = unit(25,'mm'),
                  
                  # column (sample) parameters
                  cluster_columns = FALSE,
                  show_column_dend = FALSE,
                  # column_title = '',
                  column_title_side = 'top',
                  column_title_gp = gpar(fontsize = 0), # column ticks
                  column_title_rot = 90,
                  show_column_names = FALSE,
                  column_names_gp = gpar(fontsize = 9),
                  column_names_max_height = unit(10, 'cm'),
                  # column_dend_height = unit(25,'mm')
                  
                  # specify top and bottom annotations
                  top_annotation = colAnn,
                  left_annotation = rowAnn)
  
  pdf(paste0(excel.dir, de.dir, cellType, "_complexheatmap.pdf"), height=10, width=12)
  draw(hmap, 
       # + genelabels,
       heatmap_legend_side = 'right',
       annotation_legend_side = 'right')
  dev.off()
}









##-- skip cellType B | unstim_0h Ig_4h Ig_18h | 5 donors | 80 top de genes ####
de.dir <- "stimuli_time_MAST/"


cellTypes <- c("B")

for (cellType in cellTypes){
  print(cellType)
  rds.sub <- subset(rds, idents = cellType)
  rds.sub <- subset(rds.sub, subset = expCond.stimuli %in% c("unstim", "Ig"))
  
  rds.sub@meta.data$expCond.stimuli.time.asthma <- paste(rds.sub@meta.data$expCond.stimuli.time, 
                                                      rds.sub@meta.data$asthma, sep = "_")
  
  file.all <- list.files(paste0(excel.dir, de.dir))
  file <- file.all[grepl("Ig.*top", file.all, ignore.case = TRUE)]
  # excel_data <- lapply(paste0(excel.dir, de.dir, file), read_excel)
  
  expression <- function(filename) {
    # Remove the .xlsx extension
    filename <- sub("\\.xlsx$", "", filename)
    filename <- sub("expCondCompDeMarkers_", "", filename)
    filename <- sub("-unstim_0h", "", filename)
    return(filename)
  }
  
  # Read each Excel file and add a new column with the last two underscore-separated parts of the filename
  excel_data <- lapply(file, function(f) {
    # Read the Excel file
    df <- read_excel(paste0(excel.dir, de.dir, f))
    
    # Extract the last two parts of the filename
    expression <- expression(f)
    
    # Add a new column with the extracted filename part
    df$expression <- expression
    
    return(df)
  })
  
  
  combined_data <- do.call(rbind, excel_data) %>% 
    filter(geneType %in% cellType)
  
  geneNames <- c(combined_data$Gene)
  features <- geneNames
  resDir <- paste0(excel.dir, de.dir)
  seuratObjFinal <- rds.sub
  
  # Use paste() and expand.grid() to create all combinations
  expCondReorderLevels <- apply(expand.grid(c("unstim_0h", "Ig_4h", "Ig_18h"), 
                                            sort(unique(seuratObjFinal@meta.data$expCond.donor))), 1, function(x) paste(x[2], x[1], sep = "_"))
  # expCondReorderLevels <- c("602_unstim_0h", "602_Ig_4h", "602_Ig_18h",
  #                           "616_unstim_0h", "616_Ig_4h", "616_Ig_18h",
  #                           "662_unstim_0h", "662_Ig_4h", "662_Ig_18h",
  #                           "673_unstim_0h", "673_Ig_4h", "673_Ig_18h",
  #                           "681_unstim_0h", "681_Ig_4h", "681_Ig_18h")
  slot <- "scale.data"
  Seurat::DefaultAssay(seuratObjFinal) <- "RNA"
  
  # by default, both do.scale/center are on, 
  # will scale the expression level for each gene feature by dividing the centered feature expression levels by their standard deviations
  seuratObjFinal <- Seurat::ScaleData(object = seuratObjFinal, do.scale = T, do.center = T, features = features) 
  
  # 1. input for heatmap
  data <- as.data.frame(as.matrix(t(SeuratObject::GetAssayData(object = seuratObjFinal, 
                                                               layer = slot)[features, colnames(x = seuratObjFinal), drop = FALSE]))) 
  
  metadata <- seuratObjFinal@meta.data[rownames(data), ]
  
  # focus on genes of interest
  data <- t(data[,geneNames])
  data[data < -2.5] <- -2.5
  data[data > 2.5] <- 2.5
  
  
  ## set the color scheme
  # myCol <- colorRampPalette(c("#31d3ff", "black", "#fff200"))(100)
  # myCol <- colorRampPalette(c(rev(brewer.pal(3, "Blues")), "white", brewer.pal(3, "Reds")))(100)
  # myCol <- colorRampPalette(c("blue", "white", "red"))(100)
  myCol <- colorRampPalette(c("#3B4CC0", "#FFFFFF", "#B40426"))(100)
  # myCol <- colorRampPalette(c("#3B4CC0", "#89A0D0", "#FFFFFF", "#D08080", "#B40426"))(100)
  
  # myCol <- colorRampPalette(brewer.pal(5, "RdYlBu"))(100)
  
  myBreaks <- seq(-2.5, 2.5, length.out = 100)
  
  col_ann <- data.frame(
    # Donor = factor(metadata$expCond.donor,levels=sort(unique(metadata$expCond.donor))),
    Asthma = factor(metadata$expCond.asthma,levels=sort(unique(metadata$expCond.asthma))),
    Stimuli.Time = factor(metadata$expCond.stimuli.time,
                          levels=c('unstim_0h',
                                   'Ig_4h',
                                   'Ig_18h')),
    stringsAsFactors = FALSE) #do not have to set factor 
  
  row_ann <- data.frame(
    # Genes = rownames(data),
    Stimuli.Time = rev(c(rep('Ig_4h',20), rep('Ig_4h',20),
                         rep('Ig_18h',20), rep('Ig_18h',20))),
    Expression = rev(c(rep('Top20_upDe',20), rep('Top20_downDe',20),
                       rep('Top20_upDe',20), rep('Top20_downDe',20))),
    stringsAsFactors = FALSE)
  
  
  donors <- sort(unique(metadata$expCond.donor))
  
  # Generate a color gradient
  gradient_colors <- colorRampPalette(c("#f2a65a", "#313715"))(length(donors))
  
  # Create the named vector
  donor_colors <- setNames(gradient_colors, donors)
  
  col_colours <- list(
    # Donor = donor_colors,
    Asthma = c("Asthmatic" = "#CEB992",
               "Non-asthmatic" = "#471323"),
    Stimuli.Time = c('unstim_0h' = '#eff9f0', 
                     'Ig_4h' = '#FFBCB5', 
                     'Ig_18h' = '#03978E'))
  
  row_colours <- list(
    Expression = c('Top20_downDe' = '#a8dadc',
                   'Top20_upDe' = '#e63946',
                   'Top20_downDe' = '#a8dadc',
                   'Top20_upDe' = '#e63946'),
    Stimuli.Time = c('Ig_4h' = '#FFBCB5', 
                     'Ig_18h' = '#03978E'))
  
  
  colAnn <- HeatmapAnnotation(
    df = col_ann,
    which = 'col', # 'col' (samples) or 'row' (gene) annotation?
    na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
    col = col_colours,
    annotation_height = 0.6,
    annotation_width = unit(1, 'cm'),
    simple_anno_size = unit(3, "mm"), # the height of the bar 
    gap = unit(1, 'mm'),
    annotation_legend_param = list(
      # Donor = list(
      #   nrow = 4, # number of rows across which the legend will be arranged
      #   title = 'Donor',
      #   title_position = 'topcenter',
      #   legend_direction = 'vertical',
      #   title_gp = gpar(fontsize = 11),  # legend on the side
      #   labels_gp = gpar(fontsize = 11)),  # legend on the side
      Asthma = list(
        nrow = 2, # number of rows across which the legend will be arranged
        title = 'Asthma',
        title_position = 'topcenter',
        legend_direction = 'vertical',
        title_gp = gpar(fontsize = 11),  # legend on the side
        labels_gp = gpar(fontsize = 11)),  # legend on the side
      Stimuli.Time = list(
        nrow = 3,
        title = 'Stimuli\nTime',
        title_position = 'topcenter',
        legend_direction = 'vertical',
        title_gp = gpar(fontsize = 11), # legend on the side
        labels_gp = gpar(fontsize = 11)))) # legend on the side
  
  rowAnn <- HeatmapAnnotation(
    df = row_ann,
    which = 'row', # 'col' (samples) or 'row' (gene) annotation?
    na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
    col = row_colours,
    annotation_height = 0.6,
    annotation_width = unit(1, 'cm'),
    simple_anno_size = unit(3, "mm"),
    gap = unit(1, 'mm'),
    annotation_legend_param = list(
      Stimuli.Time = list(
        ncol = 4, 
        title = 'Stimuli.Time',
        title_position = 'topcenter',
        legend_direction = 'vertical',
        title_gp = gpar(fontsize = 11),   # row annot legend title (side) 
        labels_gp = gpar(fontsize = 11)),  # row annot legend label (side)
      Expression = list(
        ncol = 4, 
        title = 'Expression',
        title_position = 'topcenter',
        legend_direction = 'vertical',
        title_gp = gpar(fontsize = 11),  # row annot legend title (side)
        labels_gp = gpar(fontsize = 11))  # row annot legend label (side)
    ))
  
  
  
  hmap <- Heatmap(data,
                  # split the genes / rows according to the PAM clusters
                  # column_split = factor(metadata$expCond.donor.stimuli.time,
                  #                       levels = expCondReorderLevels), # donor-specific
                  column_split = factor(metadata$expCond.stimuli.time.asthma,
                                        levels = apply(expand.grid(c("unstim_0h", "Ig_4h", "Ig_18h"),
                                                                   sort(unique(seuratObjFinal@meta.data$expCond.asthma))),
                                                       1, function(x) paste(x[1], x[2], sep = "_"))),
                  row_split = factor(c(rep('Ig_18h_top20_downDe',20), rep('Ig_18h_top20_upDe',20),
                                       rep('Ig_4h_top20_downDe',20), rep('Ig_4h_top20_upDe',20)),
                                     levels=c('Ig_4h_top20_upDe','Ig_4h_top20_downDe',
                                              'Ig_18h_top20_upDe','Ig_18h_top20_downDe')),
                  
                  # cluster_row_slices = FALSE,
                  name = 'Scaled\ngene\nexpression',
                  
                  col = colorRamp2(myBreaks, myCol),
                  
                  # parameters for the colour-bar that represents gradient of expression
                  heatmap_legend_param = list(
                    at = c(-2.5, -2, -1, 0, 1, 2, 2.5),
                    color_bar = 'continuous',
                    legend_direction = 'vertical',
                    legend_width = unit(6, 'cm'),
                    legend_height = unit(4, 'cm'),
                    title_position = 'topcenter',
                    title_gp=gpar(fontsize = 11), # heatmap y axis title 
                    labels_gp=gpar(fontsize = 11)), # heatmap gradient legend
                  
                  # row (gene) parameters
                  cluster_rows = FALSE,
                  show_row_dend = FALSE,
                  row_title = 'Differentially expressed genes',
                  row_title_side = 'left',
                  row_title_gp = gpar(fontsize = 11),
                  row_title_rot = 90,
                  show_row_names = TRUE,
                  row_names_gp = gpar(fontsize = 8), # y-axis ticks
                  row_names_side = 'left',
                  # row_dend_width = unit(25,'mm'),
                  
                  # column (sample) parameters
                  cluster_columns = FALSE,
                  show_column_dend = FALSE,
                  # column_title = '',
                  column_title_side = 'top',
                  column_title_gp = gpar(fontsize = 0), # column ticks
                  column_title_rot = 90,
                  show_column_names = FALSE,
                  column_names_gp = gpar(fontsize = 9),
                  column_names_max_height = unit(10, 'cm'),
                  # column_dend_height = unit(25,'mm')
                  
                  # specify top and bottom annotations
                  top_annotation = colAnn,
                  left_annotation = rowAnn,
                  
                  use_raster =FALSE)
  
  pdf(paste0(excel.dir, de.dir, cellType, "_complexheatmap.pdf"), height=10, width=20)
  draw(hmap, 
       # + genelabels,
       heatmap_legend_side = 'left',
       annotation_legend_side = 'top')
  dev.off()
}








##-- SETTING UP: subset rds for cellType B | Under Ig and unstim, & HLA ----
cellType <- "B"


rds.sub <- subset(rds, idents = cellType)

rds.sub <- subset(rds.sub, subset = (expCond.stimuli == "Ig" | expCond.stimuli == "unstim"))

rds.sub@meta.data$HLA <- 'None'

rds.sub@meta.data$HLA[which(rds.sub$RNA$data[rownames(rds.sub$RNA$data)[6458], ] != 0)] <- 'HLA-DQA2'
rds.sub@meta.data$HLA[which(rds.sub$RNA$data[rownames(rds.sub$RNA$data)[6459], ] != 0)] <- 'HLA-DQB2'

rds.sub@meta.data$HLA[
  which(
    rds.sub$RNA$data[rownames(rds.sub$RNA$data)[6458], ] != 0 & 
      rds.sub$RNA$data[rownames(rds.sub$RNA$data)[6459], ] != 0
  )
] <- 'HLA-DQA2&DQB2'

rds.sub@meta.data$expCond.stimuli.time.HLA <- paste(rds.sub@meta.data$expCond.stimuli.time, 
                                                    rds.sub@meta.data$HLA, sep = "_")
# rds.sub$expCond.stimuli.time.HLA <- gsub("unstim_.*", "unstim", rds.sub$expCond.stimuli.time.HLA)

rds.sub@meta.data$expCond.stimuli.time.HLA.donor <- paste(rds.sub@meta.data$expCond.stimuli.time.HLA, 
                                                          rds.sub@meta.data$expCond.donor, sep = "_")



saveRDS(rds.sub, file = paste0(excel.dir, "B_ig_unstim_hla.rds"))

##-- cellType B | Under Ig and unstim, what are the de genes for B cells expressing HLA-DQA|B vs None | group by asthma, hla & treatment ----
de.dir <- "HLA-DQA2|DQB2_ig_unstim_MAST/"
cellType <- "B"

rds.sub <- readRDS(paste0(excel.dir, "B_ig_unstim_hla.rds"))
rds.sub@meta.data$HLA_tmp <- rds.sub@meta.data$HLA

rds.sub@meta.data$HLA <- ifelse(grepl("None", rds.sub@meta.data$HLA), "None", 'HLA-DQA2|DQB2')

# rds.sub@meta.data %>% count(HLA)

file.all <- list.files(paste0(excel.dir, de.dir))
# file <- file.all[grepl("DQA2.*_top", file.all, ignore.case = TRUE)]

# [1] "expCondCompDeMarkers_DQA2|DQB2-None_top20_downDe.xlsx"
# [2] "expCondCompDeMarkers_DQA2|DQB2-None_top20_upDe.xlsx"   

# excel_data <- lapply(paste0(excel.dir, de.dir, file), read_excel)
# combined_data <- do.call(rbind, excel_data)

file <- file.all[grepl("*_adjSig_1SelClusters.xlsx", file.all, ignore.case = TRUE)]

# expCondCompDeMarkers_Asthmatic-Non-asthmatic_adjSig_1SelClusters.xlsx

# excel_data <- lapply(paste0(excel.dir, de.dir, file), read_excel)
# combined_data <- do.call(rbind, excel_data)

excel_data <- read_excel(paste0(excel.dir, de.dir, file))


top_down <- excel_data %>%
  filter(avg_log2FC < 0 & ...1 %in% pcgs) %>%
  arrange(-abs(avg_log2FC)) %>%  
  head(20) %>%
  rename(Gene = ...1) %>%
  mutate(geneType = "B") %>%
  select(c(Gene, geneType))


top_up <- excel_data %>%
  filter(avg_log2FC > 0 & ...1 %in% pcgs) %>%
  arrange(-abs(avg_log2FC)) %>%  
  head(20) %>%
  rename(Gene = ...1) %>%
  mutate(geneType = "B") %>%
  select(c(Gene, geneType))

combined_data <- bind_rows(top_down, top_up)


# library(biomaRt)
# 
# # Connect to Ensembl
# ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")  # Adjust species as needed
# 
# 
# # Query biotypes
# gene_info <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
#                    filters = "hgnc_symbol",
#                    values = combined_data$Gene,
#                    mart = ensembl)

# skip the line below for now 
# rds.sub <- subset(rds.sub, 
#                   subset = (expCond.donor == "602" | 
#                               expCond.donor == "681"))

geneNames <- c(combined_data$Gene, "HLA-DOA", "HLA-DOB", "IGHD")
features <- geneNames
resDir <- paste0(excel.dir, de.dir)
seuratObjFinal <- rds.sub
# expCondReorderLevels <- c("unstim_602", "Ig_4h_None_602", "Ig_18h_None_602", 
#                           "Ig_4h_DQA2|DQB2_602", "Ig_18h_DQA2|DQB2_602",
#                           "unstim_681", "Ig_4h_None_681","Ig_18h_None_681",
#                           "Ig_4h_DQA2|DQB2_681","Ig_18h_DQA2|DQB2_681")


expCondReorderLevels <- apply(expand.grid(c("unstim_0h", "Ig_4h_None", "Ig_18h_None", "Ig_4h_HLA-DQA2|DQB2", "Ig_18h_HLA-DQA2|DQB2"), 
                                          sort(unique(seuratObjFinal@meta.data$expCond.donor))), 1, function(x) paste(x[2], x[1], sep = "_"))

slot <- "scale.data"

# clusterLevels <- levels(Seurat::Idents(seuratObjFinal)) # cell types
# seuratObjFinal <- subset(seuratObjFinal, idents = "B cells")
Seurat::DefaultAssay(seuratObjFinal) <- "RNA"

# by default, both do.scale/center are on, 
# will scale the expression level for each gene feature by dividing the centered feature expression levels by their standard deviations
seuratObjFinal <- Seurat::ScaleData(object = seuratObjFinal, do.scale = T, do.center = T, features = features) 

# 1. input for heatmap
data <- as.data.frame(as.matrix(t(SeuratObject::GetAssayData(object = seuratObjFinal, 
                                                             layer = slot)[features, colnames(x = seuratObjFinal), drop = FALSE]))) 


seuratObjFinal@meta.data$expCond.stimuli.time.HLA <- paste(seuratObjFinal@meta.data$expCond.stimuli.time, 
                                                           seuratObjFinal@meta.data$HLA, sep = "_")

seuratObjFinal@meta.data$expCond.stimuli.time.HLA <- ifelse(grepl("unstim_0h", seuratObjFinal@meta.data$expCond.stimuli.time.HLA), 
                                                            "unstim_0h", 
                                                            seuratObjFinal@meta.data$expCond.stimuli.time.HLA)

seuratObjFinal@meta.data$expCond.stimuli.time.HLA.asthma <- paste(seuratObjFinal@meta.data$expCond.stimuli.time.HLA, 
                                                                  seuratObjFinal@meta.data$expCond.asthma, sep = "_")

metadata <- seuratObjFinal@meta.data[rownames(data), ]

# focus on genes of interest
geneNames = c(combined_data$Gene, "HLA-DOA", "HLA-DOB", "IGHD")
data <- t(data[,geneNames])

data[data < -2.5] <- -2.5
data[data > 2.5] <- 2.5


# 2. set the color scheme
# myCol <- colorRampPalette(c("#31d3ff", "black", "#fff200"))(100)
# myCol <- colorRampPalette(c( "#003049", "#669bbc", "#fdf0d5", "#c1121f", "#780000"))(100)
myCol <- colorRampPalette(c("#3B4CC0", "#FFFFFF", "#B40426"))(100)

# myCol <- colorRampPalette(brewer.pal(5, "RdYlBu"))(100)

myBreaks <- seq(-2.5, 2.5, length.out = 100)

col_ann <- data.frame(
  # Donor = factor(metadata$expCond.donor,levels=sort(unique(metadata$expCond.donor))),
  Asthma = factor(metadata$expCond.asthma,levels=sort(unique(metadata$expCond.asthma))),
  Stimuli.Time.HLA = factor(metadata$expCond.stimuli.time.HLA,
                            levels=c('unstim_0h',
                                     'Ig_4h_None',
                                     'Ig_18h_None',
                                     'Ig_4h_HLA-DQA2|DQB2',
                                     'Ig_18h_HLA-DQA2|DQB2')),
  stringsAsFactors = FALSE) #do not have to set factor 

row_ann <- data.frame(
  # Genes = rownames(data),
  Expression = c(rep('Down-regulated',20), rep('Up-regulated',20),
                 rep('Genes of interest',3)),
  stringsAsFactors = FALSE)

donors <- sort(unique(metadata$expCond.donor))

# Generate a color gradient
gradient_colors <- colorRampPalette(c("#f2a65a", "#313715"))(length(donors))

# Create the named vector
donor_colors <- setNames(gradient_colors, donors)

col_colours <- list(
  # Donor = donor_colors,
  Asthma = c("Asthmatic" = "#471323",
             "Non-asthmatic" = "#CEB992"),
  Stimuli.Time.HLA = c('unstim_0h' = '#eff9f0', 
                       'Ig_4h_None' = '#FFBCB5', 
                       'Ig_18h_None' = '#03978E', 
                       'Ig_4h_HLA-DQA2|DQB2' = '#e29578', 
                       'Ig_18h_HLA-DQA2|DQB2' = '#006d77'))


row_colours <- list(
  Expression = c('Down-regulated' = '#6C8FC8',
                 'Up-regulated' = '#EE4687',
                 'Genes of interest' = '#f1faee'))


colAnn <- HeatmapAnnotation(
  df = col_ann,
  which = 'col', # 'col' (samples) or 'row' (gene) annotation?
  na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
  col = col_colours,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  simple_anno_size = unit(3, "mm"),
  gap = unit(1, 'mm'),
  annotation_legend_param = list(
    # Donor = list(
    #   nrow = 4, # number of rows across which the legend will be arranged
    #   title = 'Donor',
    #   title_position = 'topcenter',
    #   legend_direction = 'vertical',
    #   title_gp = gpar(fontsize = 11),  # legend on the side
    #   labels_gp = gpar(fontsize = 11)),  # legend on the side
    Asthma = list(
      nrow = 2, # number of rows across which the legend will be arranged
      title = 'Asthma',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 11),  # legend on the side
      labels_gp = gpar(fontsize = 11)),  # legend on the side
    Stimuli.Time.HLA = list(
      nrow = 3,
      title = 'Stimuli.Time.HLA',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 11), # legend on the side
      labels_gp = gpar(fontsize = 11)))) # legend on the side


rowAnn <- HeatmapAnnotation(
  df = row_ann,
  which = 'row', # 'col' (samples) or 'row' (gene) annotation?
  na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
  col = row_colours,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  simple_anno_size = unit(3, "mm"),
  gap = unit(1, 'mm'),
  annotation_legend_param = list(
    Expression = list(
      nrow = 3, 
      title = 'Expression w/\n HLA-DQA2/DQB2',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 11),  # legend on the side
      labels_gp = gpar(fontsize = 11))  # legend on the side
  ))



hmap <- Heatmap(data,
                # split the genes / rows according to the PAM clusters
                # column_split = factor(metadata$expCond.stimuli.time.HLA.donor,
                #                       levels = expCondReorderLevels),
                # column_split = factor(metadata$expCond.asthma,
                #                       levels = c("Asthmatic", "Non-asthmatic")),
                column_split = factor(metadata$expCond.stimuli.time.HLA.asthma,
                                      levels = apply(expand.grid(c("unstim_0h", "Ig_4h_None", "Ig_18h_None",
                                                                   "Ig_4h_HLA-DQA2|DQB2", "Ig_18h_HLA-DQA2|DQB2"),
                                                                 sort(unique(seuratObjFinal@meta.data$expCond.asthma))),
                                                     1, function(x) paste(x[1], x[2], sep = "_"))),
                # column_split = factor(metadata$expCond.stimuli.time.HLA,
                #                       levels = c("unstim", "Ig_4h_None", "Ig_18h_None",
                #                                  "Ig_4h_HLA-DQA2|DQB2", "Ig_18h_HLA-DQA2|DQB2")),
                row_split = factor(c(rep('Down-regulated',20), rep('Up-regulated',20),
                                     rep('Genes of interest',3)),
                                   levels=c('Up-regulated','Down-regulated','Genes of interest')),
                
                # cluster_row_slices = FALSE,
                name = 'Scaled\ngene\nexpression',
                
                col = colorRamp2(myBreaks, myCol),
                use_raster = FALSE,
                # parameters for the colour-bar that represents gradient of expression
                heatmap_legend_param = list(
                  at = c(-2.5, -2, -1, 0, 1, 2, 2.5),
                  color_bar = 'continuous',
                  legend_direction = 'vertical',
                  legend_width = unit(6, 'cm'),
                  legend_height = unit(4, 'cm'),
                  title_position = 'topcenter',
                  title_gp=gpar(fontsize = 11),
                  labels_gp=gpar(fontsize = 11)),
                
                # row (gene) parameters
                cluster_rows = FALSE,
                show_row_dend = FALSE,
                row_title = 'Differentially expressed genes',
                row_title_side = 'left',
                row_title_gp = gpar(fontsize = 11),
                row_title_rot = 90,
                show_row_names = TRUE,
                row_names_gp = gpar(fontsize = 9), # y-axis ticks
                row_names_side = 'left',
                # row_dend_width = unit(25,'mm'),
                
                # column (sample) parameters
                cluster_columns = FALSE,
                show_column_dend = FALSE,
                # column_title = '',
                column_title_side = 'top',
                column_title_gp = gpar(fontsize = 0), # column ticks
                column_title_rot = 90,
                show_column_names = FALSE,
                column_names_gp = gpar(fontsize = 7),
                column_names_max_height = unit(10, 'cm'),
                # column_dend_height = unit(25,'mm')
                # specify top and bottom annotations
                top_annotation = colAnn,
                left_annotation = rowAnn)


pdf(paste0(excel.dir, de.dir, cellType, "_complexheatmap.pdf"), height=6.67, width=7.28)
draw(hmap, 
     # + genelabels,
     heatmap_legend_side = 'right',
     annotation_legend_side = 'top')
dev.off()

# tetst if there are more B cells expressing HLA-DQ2 genes in asthma 
contingency_table <- table(metadata$HLA, metadata$expCond.asthma)

#                 Asthmatic Non-asthmatic
# HLA-DQA2|DQB2     11184          4376
# None               7048          9677



# Perform Fisher's exact test
fisher_test_result <- fisher.test(contingency_table)
# Fisher's Exact Test for Count Data
# 
# data:  contingency_table
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  3.348830 3.677406
# sample estimates:
# odds ratio 
#   3.508898 

# strong association between the gene being expressed and asthmatic
# B cells from asthmatic donors are 3.51 times more likely to express HLA-DQA2/DQB2


metadata_unstim<- metadata %>%
  filter(expCond.stimuli.time %in% "unstim_0h")
contingency_table <- table(metadata_unstim$HLA, metadata_unstim$expCond.asthma)
res <- fisher.test(contingency_table) #p-value < 2.2e-16, odds ratio 3.810912

metadata_ig_4h <- metadata %>%
  filter(expCond.stimuli.time %in% "Ig_4h")
contingency_table <- table(metadata_ig_4h$HLA, metadata_ig_4h$expCond.asthma)
res <- fisher.test(contingency_table) #p-value < 2.2e-16, odds ratio 2.062551 




##-- cellType B | Under Ig and unstim, what are the de genes for B cells expressing HLA-DQA|B vs None | group by hla & treatment ----
de.dir <- "HLA-DQA2|DQB2_ig_unstim_MAST/"
cellType <- "B"

rds.sub <- readRDS(paste0(excel.dir, "B_ig_unstim_hla.rds"))
rds.sub@meta.data$HLA_tmp <- rds.sub@meta.data$HLA

rds.sub@meta.data$HLA <- ifelse(grepl("None", rds.sub@meta.data$HLA), "None", 'HLA-DQA2|DQB2')

# rds.sub@meta.data %>% count(HLA)

file.all <- list.files(paste0(excel.dir, de.dir))
# file <- file.all[grepl("DQA2.*_top", file.all, ignore.case = TRUE)]

# [1] "expCondCompDeMarkers_DQA2|DQB2-None_top20_downDe.xlsx"
# [2] "expCondCompDeMarkers_DQA2|DQB2-None_top20_upDe.xlsx"   

# excel_data <- lapply(paste0(excel.dir, de.dir, file), read_excel)
# combined_data <- do.call(rbind, excel_data)

file <- file.all[grepl("*_adjSig_1SelClusters.xlsx", file.all, ignore.case = TRUE)]

# expCondCompDeMarkers_Asthmatic-Non-asthmatic_adjSig_1SelClusters.xlsx

# excel_data <- lapply(paste0(excel.dir, de.dir, file), read_excel)
# combined_data <- do.call(rbind, excel_data)

excel_data <- read_excel(paste0(excel.dir, de.dir, file))


top_down <- excel_data %>%
  filter(avg_log2FC < 0 & ...1 %in% pcgs) %>%
  arrange(-abs(avg_log2FC)) %>%  
  head(20) %>%
  rename(Gene = ...1) %>%
  mutate(geneType = "B") %>%
  select(c(Gene, geneType))


top_up <- excel_data %>%
  filter(avg_log2FC > 0 & ...1 %in% pcgs) %>%
  arrange(-abs(avg_log2FC)) %>%  
  head(20) %>%
  rename(Gene = ...1) %>%
  mutate(geneType = "B") %>%
  select(c(Gene, geneType))

combined_data <- bind_rows(top_down, top_up)


# library(biomaRt)
# 
# # Connect to Ensembl
# ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")  # Adjust species as needed
# 
# 
# # Query biotypes
# gene_info <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
#                    filters = "hgnc_symbol",
#                    values = combined_data$Gene,
#                    mart = ensembl)

# skip the line below for now 
# rds.sub <- subset(rds.sub, 
#                   subset = (expCond.donor == "602" | 
#                               expCond.donor == "681"))

geneNames <- c(combined_data$Gene, "HLA-DOA", "HLA-DOB", "IGHD")
features <- geneNames
resDir <- paste0(excel.dir, de.dir)
seuratObjFinal <- rds.sub
# expCondReorderLevels <- c("unstim_602", "Ig_4h_None_602", "Ig_18h_None_602", 
#                           "Ig_4h_DQA2|DQB2_602", "Ig_18h_DQA2|DQB2_602",
#                           "unstim_681", "Ig_4h_None_681","Ig_18h_None_681",
#                           "Ig_4h_DQA2|DQB2_681","Ig_18h_DQA2|DQB2_681")


expCondReorderLevels <- apply(expand.grid(c("unstim_0h", "Ig_4h_None", "Ig_18h_None", "Ig_4h_HLA-DQA2|DQB2", "Ig_18h_HLA-DQA2|DQB2"), 
                                          sort(unique(seuratObjFinal@meta.data$expCond.donor))), 1, function(x) paste(x[2], x[1], sep = "_"))

slot <- "scale.data"

# clusterLevels <- levels(Seurat::Idents(seuratObjFinal)) # cell types
# seuratObjFinal <- subset(seuratObjFinal, idents = "B cells")
Seurat::DefaultAssay(seuratObjFinal) <- "RNA"

# by default, both do.scale/center are on, 
# will scale the expression level for each gene feature by dividing the centered feature expression levels by their standard deviations
seuratObjFinal <- Seurat::ScaleData(object = seuratObjFinal, do.scale = T, do.center = T, features = features) 

# 1. input for heatmap
data <- as.data.frame(as.matrix(t(SeuratObject::GetAssayData(object = seuratObjFinal, 
                                                             layer = slot)[features, colnames(x = seuratObjFinal), drop = FALSE]))) 


seuratObjFinal@meta.data$expCond.stimuli.time.HLA <- paste(seuratObjFinal@meta.data$expCond.stimuli.time, 
                                                           seuratObjFinal@meta.data$HLA, sep = "_")

# seuratObjFinal@meta.data$expCond.stimuli.time.HLA <- ifelse(grepl("unstim_0h", seuratObjFinal@meta.data$expCond.stimuli.time.HLA), 
#                                                             "unstim_0h", 
#                                                             seuratObjFinal@meta.data$expCond.stimuli.time.HLA)

# seuratObjFinal@meta.data$expCond.stimuli.time.HLA.asthma <- paste(seuratObjFinal@meta.data$expCond.stimuli.time.HLA, 
#                                                                   seuratObjFinal@meta.data$expCond.asthma, sep = "_")
# 
metadata <- seuratObjFinal@meta.data[rownames(data), ]

# focus on genes of interest
geneNames = c(combined_data$Gene, "HLA-DOA", "HLA-DOB", "IGHD")
data <- t(data[,geneNames])

data[data < -2.5] <- -2.5
data[data > 2.5] <- 2.5


# 2. set the color scheme
# myCol <- colorRampPalette(c("#31d3ff", "black", "#fff200"))(100)
# myCol <- colorRampPalette(c( "#003049", "#669bbc", "#fdf0d5", "#c1121f", "#780000"))(100)
myCol <- colorRampPalette(c("#3B4CC0", "#FFFFFF", "#B40426"))(100)

# myCol <- colorRampPalette(brewer.pal(5, "RdYlBu"))(100)

myBreaks <- seq(-2.5, 2.5, length.out = 100)

col_ann <- data.frame(
  # Donor = factor(metadata$expCond.donor,levels=sort(unique(metadata$expCond.donor))),
  # Asthma = factor(metadata$expCond.asthma,levels=sort(unique(metadata$expCond.asthma))),
  Stimuli.Time.HLA = factor(metadata$expCond.stimuli.time.HLA,
                            levels=c('unstim_0h_None',
                                     'Ig_4h_None',
                                     'Ig_18h_None',
                                     'unstim_0h_HLA-DQA2|DQB2',
                                     'Ig_4h_HLA-DQA2|DQB2',
                                     'Ig_18h_HLA-DQA2|DQB2')),
  stringsAsFactors = FALSE) #do not have to set factor 

row_ann <- data.frame(
  # Genes = rownames(data),
  Expression = c(rep('Down-regulated',20), rep('Up-regulated',20),
                 rep('Genes of interest',3)),
  stringsAsFactors = FALSE)

donors <- sort(unique(metadata$expCond.donor))

# Generate a color gradient
gradient_colors <- colorRampPalette(c("#f2a65a", "#313715"))(length(donors))

# Create the named vector
donor_colors <- setNames(gradient_colors, donors)

col_colours <- list(
  # Donor = donor_colors,
  # Asthma = c("Asthmatic" = "#471323",
  #            "Non-asthmatic" = "#CEB992"),
  Stimuli.Time.HLA = c('unstim_0h_None' = '#eff9f0', 
                       'Ig_4h_None' = '#FFBCB5', 
                       'Ig_18h_None' = '#03978E', 
                       'unstim_0h_HLA-DQA2|DQB2' = '#b2d4ca',
                       'Ig_4h_HLA-DQA2|DQB2' = '#e29578', 
                       'Ig_18h_HLA-DQA2|DQB2' = '#006d77'))


row_colours <- list(
  Expression = c('Down-regulated' = '#6C8FC8',
                 'Up-regulated' = '#EE4687',
                 'Genes of interest' = '#f1faee'))


colAnn <- HeatmapAnnotation(
  df = col_ann,
  which = 'col', # 'col' (samples) or 'row' (gene) annotation?
  na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
  col = col_colours,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  simple_anno_size = unit(3, "mm"),
  gap = unit(1, 'mm'),
  annotation_legend_param = list(
    # Donor = list(
    #   nrow = 4, # number of rows across which the legend will be arranged
    #   title = 'Donor',
    #   title_position = 'topcenter',
    #   legend_direction = 'vertical',
    #   title_gp = gpar(fontsize = 11),  # legend on the side
    #   labels_gp = gpar(fontsize = 11)),  # legend on the side
    # Asthma = list(
    #   nrow = 2, # number of rows across which the legend will be arranged
    #   title = 'Asthma',
    #   title_position = 'topcenter',
    #   legend_direction = 'vertical',
    #   title_gp = gpar(fontsize = 11),  # legend on the side
    #   labels_gp = gpar(fontsize = 11)),  # legend on the side
    Stimuli.Time.HLA = list(
      nrow = 3,
      title = 'Stimuli.Time.HLA',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 11), # legend on the side
      labels_gp = gpar(fontsize = 11)))) # legend on the side


rowAnn <- HeatmapAnnotation(
  df = row_ann,
  which = 'row', # 'col' (samples) or 'row' (gene) annotation?
  na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
  col = row_colours,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  simple_anno_size = unit(3, "mm"),
  gap = unit(1, 'mm'),
  annotation_legend_param = list(
    Expression = list(
      nrow = 3, 
      title = 'Expression w/\n HLA-DQA2/DQB2',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 11),  # legend on the side
      labels_gp = gpar(fontsize = 11))  # legend on the side
  ))



hmap <- Heatmap(data,
                # split the genes / rows according to the PAM clusters
                # column_split = factor(metadata$expCond.stimuli.time.HLA.donor,
                #                       levels = expCondReorderLevels),
                # column_split = factor(metadata$expCond.asthma,
                #                       levels = c("Asthmatic", "Non-asthmatic")),
                # column_split = factor(metadata$expCond.stimuli.time.HLA.asthma,
                #                       levels = apply(expand.grid(c("unstim_0h", "Ig_4h_None", "Ig_18h_None",
                #                                                    "Ig_4h_HLA-DQA2|DQB2", "Ig_18h_HLA-DQA2|DQB2"),
                #                                                  sort(unique(seuratObjFinal@meta.data$expCond.asthma))),
                #                                      1, function(x) paste(x[1], x[2], sep = "_"))),
                column_split = factor(metadata$expCond.stimuli.time.HLA,
                                      levels = c("unstim_0h_None", "Ig_4h_None", "Ig_18h_None",
                                                 'unstim_0h_HLA-DQA2|DQB2', "Ig_4h_HLA-DQA2|DQB2", "Ig_18h_HLA-DQA2|DQB2")),
                row_split = factor(c(rep('Down-regulated',20), rep('Up-regulated',20),
                                     rep('Genes of interest',3)),
                                   levels=c('Up-regulated','Down-regulated','Genes of interest')),
                
                # cluster_row_slices = FALSE,
                name = 'Scaled\ngene\nexpression',
                
                col = colorRamp2(myBreaks, myCol),
                use_raster = FALSE,
                # parameters for the colour-bar that represents gradient of expression
                heatmap_legend_param = list(
                  at = c(-2.5, -2, -1, 0, 1, 2, 2.5),
                  color_bar = 'continuous',
                  legend_direction = 'vertical',
                  legend_width = unit(6, 'cm'),
                  legend_height = unit(4, 'cm'),
                  title_position = 'topcenter',
                  title_gp=gpar(fontsize = 11),
                  labels_gp=gpar(fontsize = 11)),
                
                # row (gene) parameters
                cluster_rows = FALSE,
                show_row_dend = FALSE,
                row_title = 'Differentially expressed genes',
                row_title_side = 'left',
                row_title_gp = gpar(fontsize = 11),
                row_title_rot = 90,
                show_row_names = TRUE,
                row_names_gp = gpar(fontsize = 9), # y-axis ticks
                row_names_side = 'left',
                # row_dend_width = unit(25,'mm'),
                
                # column (sample) parameters
                cluster_columns = FALSE,
                show_column_dend = FALSE,
                # column_title = '',
                column_title_side = 'top',
                column_title_gp = gpar(fontsize = 0), # column ticks
                column_title_rot = 90,
                show_column_names = FALSE,
                column_names_gp = gpar(fontsize = 7),
                column_names_max_height = unit(10, 'cm'),
                # column_dend_height = unit(25,'mm')
                # specify top and bottom annotations
                top_annotation = colAnn,
                left_annotation = rowAnn)


pdf(paste0(excel.dir, de.dir, cellType, "_hla_cond_complexheatmap.pdf"), height=6.67, width=7.28)
draw(hmap, 
     # + genelabels,
     heatmap_legend_side = 'right',
     annotation_legend_side = 'top')
dev.off()

# tetst if there are more B cells expressing HLA-DQ2 genes in asthma 
contingency_table <- table(metadata$HLA, metadata$expCond.asthma)

#                 Asthmatic Non-asthmatic
# HLA-DQA2|DQB2     11184          4376
# None               7048          9677



# Perform Fisher's exact test
fisher_test_result <- fisher.test(contingency_table)
# Fisher's Exact Test for Count Data
# 
# data:  contingency_table
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  3.348830 3.677406
# sample estimates:
# odds ratio 
#   3.508898 

# strong association between the gene being expressed and asthmatic
# B cells from asthmatic donors are 3.51 times more likely to express HLA-DQA2/DQB2


metadata_unstim<- metadata %>%
  filter(expCond.stimuli.time %in% "unstim_0h")
contingency_table <- table(metadata_unstim$HLA, metadata_unstim$expCond.asthma)
res <- fisher.test(contingency_table) #p-value < 2.2e-16, odds ratio 3.810912

metadata_ig_4h <- metadata %>%
  filter(expCond.stimuli.time %in% "Ig_4h")
contingency_table <- table(metadata_ig_4h$HLA, metadata_ig_4h$expCond.asthma)
res <- fisher.test(contingency_table) #p-value < 2.2e-16, odds ratio 2.062551 





##-- cellType B | Under Ig and unstim, what are the de genes for B cells expressing HLA-DQA only ----
de.dir <- "HLA-DQA2_ig_unstim_MAST/"
cellType <- "B"

rds.sub <- readRDS(paste0(excel.dir, "B_ig_unstim_hla.rds"))

rds.sub@meta.data$HLA <- ifelse(grepl("DQA2", rds.sub@meta.data$HLA), "HLA-DQA2", "None")


# rds.sub@meta.data %>% count(HLA)

file.all <- list.files(paste0(excel.dir, de.dir))
# file <- file.all[grepl("DQA2.*_top", file.all, ignore.case = TRUE)]

# [1] "expCondCompDeMarkers_DQA2|DQB2-None_top20_downDe.xlsx"
# [2] "expCondCompDeMarkers_DQA2|DQB2-None_top20_upDe.xlsx"   

# excel_data <- lapply(paste0(excel.dir, de.dir, file), read_excel)
# combined_data <- do.call(rbind, excel_data)
file <- file.all[grepl("*_adjSig_1SelClusters.xlsx", file.all, ignore.case = TRUE)]
excel_data <- read_excel(paste0(excel.dir, de.dir, file))


top_down <- excel_data %>%
  filter(avg_log2FC < 0 & ...1 %in% pcgs) %>%
  arrange(-abs(avg_log2FC)) %>%  
  head(20) %>%
  rename(Gene = ...1) %>%
  mutate(geneType = "B") %>%
  select(c(Gene, geneType))


top_up <- excel_data %>%
  filter(avg_log2FC > 0 & ...1 %in% pcgs) %>%
  arrange(-abs(avg_log2FC)) %>%  
  head(20) %>%
  rename(Gene = ...1) %>%
  mutate(geneType = "B") %>%
  select(c(Gene, geneType))

combined_data <- bind_rows(top_down, top_up)


# library(biomaRt)
# 
# # Connect to Ensembl
# ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")  # Adjust species as needed
# 
# 
# # Query biotypes
# gene_info <- getBM(attributes = c("hgnc_symbol", "gene_biotype"),
#                    filters = "hgnc_symbol",
#                    values = combined_data$Gene,
#                    mart = ensembl)

# skip the line below for now 
# rds.sub <- subset(rds.sub, 
#                   subset = (expCond.donor == "602" | 
#                               expCond.donor == "681"))

geneNames <- c(combined_data$Gene, "HLA-DOA", "HLA-DOB", "IGHD")
features <- geneNames
resDir <- paste0(excel.dir, de.dir)
seuratObjFinal <- rds.sub
# expCondReorderLevels <- c("unstim_602", "Ig_4h_None_602", "Ig_18h_None_602", 
#                           "Ig_4h_DQA2|DQB2_602", "Ig_18h_DQA2|DQB2_602",
#                           "unstim_681", "Ig_4h_None_681","Ig_18h_None_681",
#                           "Ig_4h_DQA2|DQB2_681","Ig_18h_DQA2|DQB2_681")


expCondReorderLevels <- apply(expand.grid(c("unstim_0h", "Ig_4h_None", "Ig_18h_None", "Ig_4h_HLA-DQA2", "Ig_18h_HLA-DQA2"), 
                                          sort(unique(seuratObjFinal@meta.data$expCond.donor))), 
                              1, function(x) paste(x[2], x[1], sep = "_"))

slot <- "scale.data"

# clusterLevels <- levels(Seurat::Idents(seuratObjFinal)) # cell types
# seuratObjFinal <- subset(seuratObjFinal, idents = "B cells")
Seurat::DefaultAssay(seuratObjFinal) <- "RNA"

# by default, both do.scale/center are on, 
# will scale the expression level for each gene feature by dividing the centered feature expression levels by their standard deviations
seuratObjFinal <- Seurat::ScaleData(object = seuratObjFinal, do.scale = T, do.center = T, features = features) 

# 1. input for heatmap
data <- as.data.frame(as.matrix(t(SeuratObject::GetAssayData(object = seuratObjFinal, 
                                                             layer = slot)[features, colnames(x = seuratObjFinal), drop = FALSE]))) 


seuratObjFinal@meta.data$expCond.stimuli.time.HLA <- paste(seuratObjFinal@meta.data$expCond.stimuli.time, 
                                                           seuratObjFinal@meta.data$HLA, sep = "_")

seuratObjFinal@meta.data$expCond.stimuli.time.HLA <- ifelse(grepl("unstim_0h", seuratObjFinal@meta.data$expCond.stimuli.time.HLA), 
                                                            "unstim_0h", 
                                                            seuratObjFinal@meta.data$expCond.stimuli.time.HLA)

seuratObjFinal@meta.data$expCond.stimuli.time.HLA.asthma <- paste(seuratObjFinal@meta.data$expCond.stimuli.time.HLA, 
                                                                  seuratObjFinal@meta.data$expCond.asthma, sep = "_")

metadata <- seuratObjFinal@meta.data[rownames(data), ]



# focus on genes of interest
geneNames = c(combined_data$Gene, "HLA-DOA", "HLA-DOB", "IGHD")
data <- t(data[,geneNames])

data[data < -2.5] <- -2.5
data[data > 2.5] <- 2.5


# 2. set the color scheme
# myCol <- colorRampPalette(c("#31d3ff", "black", "#fff200"))(100)
# myCol <- colorRampPalette(c( "#003049", "#669bbc", "#fdf0d5", "#c1121f", "#780000"))(100)
myCol <- colorRampPalette(c("#3B4CC0", "#FFFFFF", "#B40426"))(100)

# myCol <- colorRampPalette(brewer.pal(5, "RdYlBu"))(100)

myBreaks <- seq(-2.5, 2.5, length.out = 100)

col_ann <- data.frame(
  # Donor = factor(metadata$expCond.donor,levels=sort(unique(metadata$expCond.donor))),
  Asthma = factor(metadata$expCond.asthma,levels=sort(unique(metadata$expCond.asthma))),
  Stimuli.Time.HLA = factor(metadata$expCond.stimuli.time.HLA,
                            levels=c('unstim_0h',
                                     'Ig_4h_None',
                                     'Ig_18h_None',
                                     'Ig_4h_HLA-DQA2',
                                     'Ig_18h_HLA-DQA2')),
  stringsAsFactors = FALSE) #do not have to set factor 

row_ann <- data.frame(
  # Genes = rownames(data),
  Expression = c(rep('Down-regulated',20), rep('Up-regulated',20),
                 rep('Genes of interest',3)),
  stringsAsFactors = FALSE)

donors <- sort(unique(metadata$expCond.donor))

# Generate a color gradient
gradient_colors <- colorRampPalette(c("#f2a65a", "#313715"))(length(donors))

# Create the named vector
donor_colors <- setNames(gradient_colors, donors)

col_colours <- list(
  # Donor = donor_colors,
  Asthma = c("Asthmatic" = "#471323",
             "Non-asthmatic" = "#CEB992"),
  Stimuli.Time.HLA = c('unstim_0h' = '#eff9f0', 
                       'Ig_4h_None' = '#FFBCB5', 
                       'Ig_18h_None' = '#03978E', 
                       'Ig_4h_HLA-DQA2' = '#e29578', 
                       'Ig_18h_HLA-DQA2' = '#006d77'))


row_colours <- list(
  Expression = c('Down-regulated' = '#6C8FC8',
                 'Up-regulated' = '#EE4687',
                 'Genes of interest' = '#f1faee'))


colAnn <- HeatmapAnnotation(
  df = col_ann,
  which = 'col', # 'col' (samples) or 'row' (gene) annotation?
  na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
  col = col_colours,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  simple_anno_size = unit(3, "mm"),
  gap = unit(1, 'mm'),
  annotation_legend_param = list(
    # Donor = list(
    #   nrow = 4, # number of rows across which the legend will be arranged
    #   title = 'Donor',
    #   title_position = 'topcenter',
    #   legend_direction = 'vertical',
    #   title_gp = gpar(fontsize = 11),  # legend on the side
    #   labels_gp = gpar(fontsize = 11)),  # legend on the side
    Asthma = list(
      nrow = 2, # number of rows across which the legend will be arranged
      title = 'Asthma',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 11),  # legend on the side
      labels_gp = gpar(fontsize = 11)),  # legend on the side
    Stimuli.Time.HLA = list(
      nrow = 3,
      title = 'Stimuli.Time.HLA',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 11), # legend on the side
      labels_gp = gpar(fontsize = 11)))) # legend on the side


rowAnn <- HeatmapAnnotation(
  df = row_ann,
  which = 'row', # 'col' (samples) or 'row' (gene) annotation?
  na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
  col = row_colours,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  simple_anno_size = unit(3, "mm"),
  gap = unit(1, 'mm'),
  annotation_legend_param = list(
    Expression = list(
      nrow = 3, 
      title = 'Expression w/\n HLA-DQA2',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 11),  # legend on the side
      labels_gp = gpar(fontsize = 11))  # legend on the side
  ))



hmap <- Heatmap(data,
                # split the genes / rows according to the PAM clusters
                # column_split = factor(metadata$expCond.stimuli.time.HLA.donor,
                #                       levels = expCondReorderLevels),
                # column_split = factor(metadata$expCond.asthma,
                #                       levels = c("Asthmatic", "Non-asthmatic")),
                column_split = factor(metadata$expCond.stimuli.time.HLA.asthma,
                                      levels = apply(expand.grid(c("unstim_0h", "Ig_4h_None", "Ig_18h_None",
                                                                   "Ig_4h_HLA-DQA2", "Ig_18h_HLA-DQA2"),
                                                                 sort(unique(seuratObjFinal@meta.data$expCond.asthma))),
                                                     1, function(x) paste(x[1], x[2], sep = "_"))),
                # column_split = factor(metadata$expCond.stimuli.time.HLA,
                #                       levels = c("unstim", "Ig_4h_None", "Ig_18h_None",
                #                                  "Ig_4h_HLA-DQA2", "Ig_18h_HLA-DQA2")),
                row_split = factor(c(rep('Down-regulated',20), rep('Up-regulated',20),
                                     rep('Genes of interest',3)),
                                   levels=c('Up-regulated','Down-regulated','Genes of interest')),
                
                # cluster_row_slices = FALSE,
                name = 'Scaled\ngene\nexpression',
                
                col = colorRamp2(myBreaks, myCol),
                use_raster = FALSE,
                # parameters for the colour-bar that represents gradient of expression
                heatmap_legend_param = list(
                  at = c(-2.5, -2, -1, 0, 1, 2, 2.5),
                  color_bar = 'continuous',
                  legend_direction = 'vertical',
                  legend_width = unit(6, 'cm'),
                  legend_height = unit(4, 'cm'),
                  title_position = 'topcenter',
                  title_gp=gpar(fontsize = 11),
                  labels_gp=gpar(fontsize = 11)),
                
                # row (gene) parameters
                cluster_rows = FALSE,
                show_row_dend = FALSE,
                row_title = 'Differentially expressed genes',
                row_title_side = 'left',
                row_title_gp = gpar(fontsize = 11),
                row_title_rot = 90,
                show_row_names = TRUE,
                row_names_gp = gpar(fontsize = 9), # y-axis ticks
                row_names_side = 'left',
                # row_dend_width = unit(25,'mm'),
                
                # column (sample) parameters
                cluster_columns = FALSE,
                show_column_dend = FALSE,
                # column_title = '',
                column_title_side = 'top',
                column_title_gp = gpar(fontsize = 0), # column ticks
                column_title_rot = 90,
                show_column_names = FALSE,
                column_names_gp = gpar(fontsize = 7),
                column_names_max_height = unit(10, 'cm'),
                # column_dend_height = unit(25,'mm')
                # specify top and bottom annotations
                top_annotation = colAnn,
                left_annotation = rowAnn)


pdf(paste0(excel.dir, de.dir, cellType, "_complexheatmap.pdf"), height=6.67, width=7.28)
draw(hmap, 
     # + genelabels,
     heatmap_legend_side = 'right',
     annotation_legend_side = 'top')
dev.off()

# tetst if there are more B cells expressing HLA-DQA2 genes in asthma 
contingency_table <- table(metadata$HLA, metadata$expCond.asthma)

# Asthmatic Non-asthmatic
# HLA-DQA2     10179          3257
# None          8053         10796



# Perform Fisher's exact test
fisher_test_result <- fisher.test(contingency_table)
# Fisher's Exact Test for Count Data
# 
# data:  contingency_table
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  3.988483 4.400806
# sample estimates:
# odds ratio 
#   4.189537 


# strong association between the gene being expressed and asthmatic
# B cells from asthmatic donors are 4.189537 times more likely to express HLA-DQA2\





##-- cellType B | Under Ig and unstim, what are the de genes for B cells expressing HLA-DQB only ----
de.dir <- "HLA-DQB2_ig_unstim_MAST/"
cellType <- "B"

rds.sub <- readRDS(paste0(excel.dir, "B_ig_unstim_hla.rds"))

rds.sub@meta.data$HLA <- ifelse(grepl("DQB2", rds.sub@meta.data$HLA), "HLA-DQB2", "None")


# rds.sub@meta.data %>% count(HLA)

file.all <- list.files(paste0(excel.dir, de.dir))
# file <- file.all[grepl("DQB2.*_top", file.all, ignore.case = TRUE)]

# [1] "expCondCompDeMarkers_DQA2|DQB2-None_top20_downDe.xlsx"
# [2] "expCondCompDeMarkers_DQA2|DQB2-None_top20_upDe.xlsx"   

# excel_data <- lapply(paste0(excel.dir, de.dir, file), read_excel)
# combined_data <- do.call(rbind, excel_data)

file <- file.all[grepl("*_adjSig_1SelClusters.xlsx", file.all, ignore.case = TRUE)]
excel_data <- read_excel(paste0(excel.dir, de.dir, file))


top_down <- excel_data %>%
  filter(avg_log2FC < 0 & ...1 %in% pcgs) %>%
  arrange(-abs(avg_log2FC)) %>%  
  head(20) %>%
  rename(Gene = ...1) %>%
  mutate(geneType = "B") %>%
  select(c(Gene, geneType))


top_up <- excel_data %>%
  filter(avg_log2FC > 0 & ...1 %in% pcgs) %>%
  arrange(-abs(avg_log2FC)) %>%  
  head(20) %>%
  rename(Gene = ...1) %>%
  mutate(geneType = "B") %>%
  select(c(Gene, geneType))

combined_data <- bind_rows(top_down, top_up)



# library(biomaRt)
# 
# # Connect to Ensembl
# ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")  # Adjust species as needed
# 
# 
# # Query biotypes
# gene_info <- getBM(attributes = c("hgnc_symbol", "gene_biotype"),
#                    filters = "hgnc_symbol",
#                    values = combined_data$Gene,
#                    mart = ensembl)

# skip the line below for now 
# rds.sub <- subset(rds.sub, 
#                   subset = (expCond.donor == "602" | 
#                               expCond.donor == "681"))

geneNames <- c(combined_data$Gene, "HLA-DOA", "HLA-DOB", "IGHD")
features <- geneNames
resDir <- paste0(excel.dir, de.dir)
seuratObjFinal <- rds.sub
# expCondReorderLevels <- c("unstim_602", "Ig_4h_None_602", "Ig_18h_None_602", 
#                           "Ig_4h_DQA2|DQB2_602", "Ig_18h_DQA2|DQB2_602",
#                           "unstim_681", "Ig_4h_None_681","Ig_18h_None_681",
#                           "Ig_4h_DQA2|DQB2_681","Ig_18h_DQA2|DQB2_681")


expCondReorderLevels <- apply(expand.grid(c("unstim_0h", "Ig_4h_None", "Ig_18h_None", "Ig_4h_HLA-DQB2", "Ig_18h_HLA-DQB2"), 
                                          sort(unique(seuratObjFinal@meta.data$expCond.donor))), 
                              1, function(x) paste(x[2], x[1], sep = "_"))

slot <- "scale.data"

# clusterLevels <- levels(Seurat::Idents(seuratObjFinal)) # cell types
# seuratObjFinal <- subset(seuratObjFinal, idents = "B cells")
Seurat::DefaultAssay(seuratObjFinal) <- "RNA"

# by default, both do.scale/center are on, 
# will scale the expression level for each gene feature by dividing the centered feature expression levels by their standard deviations
seuratObjFinal <- Seurat::ScaleData(object = seuratObjFinal, do.scale = T, do.center = T, features = features) 

# 1. input for heatmap
data <- as.data.frame(as.matrix(t(SeuratObject::GetAssayData(object = seuratObjFinal, 
                                                             layer = slot)[features, colnames(x = seuratObjFinal), drop = FALSE]))) 


seuratObjFinal@meta.data$expCond.stimuli.time.HLA <- paste(seuratObjFinal@meta.data$expCond.stimuli.time, 
                                                           seuratObjFinal@meta.data$HLA, sep = "_")

seuratObjFinal@meta.data$expCond.stimuli.time.HLA <- ifelse(grepl("unstim_0h", seuratObjFinal@meta.data$expCond.stimuli.time.HLA), 
                                                            "unstim_0h", 
                                                            seuratObjFinal@meta.data$expCond.stimuli.time.HLA)

seuratObjFinal@meta.data$expCond.stimuli.time.HLA.asthma <- paste(seuratObjFinal@meta.data$expCond.stimuli.time.HLA, 
                                                                  seuratObjFinal@meta.data$expCond.asthma, sep = "_")

metadata <- seuratObjFinal@meta.data[rownames(data), ]



# focus on genes of interest
geneNames = c(combined_data$Gene, "HLA-DOA", "HLA-DOB", "IGHD")
data <- t(data[,geneNames])

data[data < -2.5] <- -2.5
data[data > 2.5] <- 2.5


# 2. set the color scheme
# myCol <- colorRampPalette(c("#31d3ff", "black", "#fff200"))(100)
# myCol <- colorRampPalette(c( "#003049", "#669bbc", "#fdf0d5", "#c1121f", "#780000"))(100)
myCol <- colorRampPalette(c("#3B4CC0", "#FFFFFF", "#B40426"))(100)

# myCol <- colorRampPalette(brewer.pal(5, "RdYlBu"))(100)

myBreaks <- seq(-2.5, 2.5, length.out = 100)

col_ann <- data.frame(
  # Donor = factor(metadata$expCond.donor,levels=sort(unique(metadata$expCond.donor))),
  Asthma = factor(metadata$expCond.asthma,levels=sort(unique(metadata$expCond.asthma))),
  Stimuli.Time.HLA = factor(metadata$expCond.stimuli.time.HLA,
                            levels=c('unstim_0h',
                                     'Ig_4h_None',
                                     'Ig_18h_None',
                                     'Ig_4h_HLA-DQB2',
                                     'Ig_18h_HLA-DQB2')),
  stringsAsFactors = FALSE) #do not have to set factor 

row_ann <- data.frame(
  # Genes = rownames(data),
  Expression = c(rep('Down-regulated',20), rep('Up-regulated',20),
                 rep('Genes of interest',3)),
  stringsAsFactors = FALSE)

donors <- sort(unique(metadata$expCond.donor))

# Generate a color gradient
gradient_colors <- colorRampPalette(c("#f2a65a", "#313715"))(length(donors))

# Create the named vector
donor_colors <- setNames(gradient_colors, donors)

col_colours <- list(
  # Donor = donor_colors,
  Asthma = c("Asthmatic" = "#471323",
             "Non-asthmatic" = "#CEB992"),
  Stimuli.Time.HLA = c('unstim_0h' = '#eff9f0', 
                       'Ig_4h_None' = '#FFBCB5', 
                       'Ig_18h_None' = '#03978E', 
                       'Ig_4h_HLA-DQB2' = '#e29578', 
                       'Ig_18h_HLA-DQB2' = '#006d77'))


row_colours <- list(
  Expression = c('Down-regulated' = '#6C8FC8',
                 'Up-regulated' = '#EE4687',
                 'Genes of interest' = '#f1faee'))


colAnn <- HeatmapAnnotation(
  df = col_ann,
  which = 'col', # 'col' (samples) or 'row' (gene) annotation?
  na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
  col = col_colours,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  simple_anno_size = unit(3, "mm"),
  gap = unit(1, 'mm'),
  annotation_legend_param = list(
    # Donor = list(
    #   nrow = 4, # number of rows across which the legend will be arranged
    #   title = 'Donor',
    #   title_position = 'topcenter',
    #   legend_direction = 'vertical',
    #   title_gp = gpar(fontsize = 11),  # legend on the side
    #   labels_gp = gpar(fontsize = 11)),  # legend on the side
    Asthma = list(
      nrow = 2, # number of rows across which the legend will be arranged
      title = 'Asthma',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 11),  # legend on the side
      labels_gp = gpar(fontsize = 11)),  # legend on the side
    Stimuli.Time.HLA = list(
      nrow = 3,
      title = 'Stimuli.Time.HLA',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 11), # legend on the side
      labels_gp = gpar(fontsize = 11)))) # legend on the side


rowAnn <- HeatmapAnnotation(
  df = row_ann,
  which = 'row', # 'col' (samples) or 'row' (gene) annotation?
  na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
  col = row_colours,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  simple_anno_size = unit(3, "mm"),
  gap = unit(1, 'mm'),
  annotation_legend_param = list(
    Expression = list(
      nrow = 3, 
      title = 'Expression w/\n HLA-DQB2',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 11),  # legend on the side
      labels_gp = gpar(fontsize = 11))  # legend on the side
  ))



hmap <- Heatmap(data,
                # split the genes / rows according to the PAM clusters
                # column_split = factor(metadata$expCond.stimuli.time.HLA.donor,
                #                       levels = expCondReorderLevels),
                # column_split = factor(metadata$expCond.asthma,
                #                       levels = c("Asthmatic", "Non-asthmatic")),
                column_split = factor(metadata$expCond.stimuli.time.HLA.asthma,
                                      levels = apply(expand.grid(c("unstim_0h", "Ig_4h_None", "Ig_18h_None",
                                                                   "Ig_4h_HLA-DQB2", "Ig_18h_HLA-DQB2"),
                                                                 sort(unique(seuratObjFinal@meta.data$expCond.asthma))),
                                                     1, function(x) paste(x[1], x[2], sep = "_"))),
                # column_split = factor(metadata$expCond.stimuli.time.HLA,
                #                       levels = c("unstim", "Ig_4h_None", "Ig_18h_None",
                #                                  "Ig_4h_HLA-DQB2", "Ig_18h_HLA-DQB2")),
                row_split = factor(c(rep('Down-regulated',20), rep('Up-regulated',20),
                                     rep('Genes of interest',3)),
                                   levels=c('Up-regulated','Down-regulated','Genes of interest')),
                
                # cluster_row_slices = FALSE,
                name = 'Scaled\ngene\nexpression',
                
                col = colorRamp2(myBreaks, myCol),
                use_raster = FALSE,
                # parameters for the colour-bar that represents gradient of expression
                heatmap_legend_param = list(
                  at = c(-2.5, -2, -1, 0, 1, 2, 2.5),
                  color_bar = 'continuous',
                  legend_direction = 'vertical',
                  legend_width = unit(6, 'cm'),
                  legend_height = unit(4, 'cm'),
                  title_position = 'topcenter',
                  title_gp=gpar(fontsize = 11),
                  labels_gp=gpar(fontsize = 11)),
                
                # row (gene) parameters
                cluster_rows = FALSE,
                show_row_dend = FALSE,
                row_title = 'Differentially expressed genes',
                row_title_side = 'left',
                row_title_gp = gpar(fontsize = 11),
                row_title_rot = 90,
                show_row_names = TRUE,
                row_names_gp = gpar(fontsize = 9), # y-axis ticks
                row_names_side = 'left',
                # row_dend_width = unit(25,'mm'),
                
                # column (sample) parameters
                cluster_columns = FALSE,
                show_column_dend = FALSE,
                # column_title = '',
                column_title_side = 'top',
                column_title_gp = gpar(fontsize = 0), # column ticks
                column_title_rot = 90,
                show_column_names = FALSE,
                column_names_gp = gpar(fontsize = 7),
                column_names_max_height = unit(10, 'cm'),
                # column_dend_height = unit(25,'mm')
                # specify top and bottom annotations
                top_annotation = colAnn,
                left_annotation = rowAnn)


pdf(paste0(excel.dir, de.dir, cellType, "_complexheatmap.pdf"), height=6.67, width=7.28)
draw(hmap, 
     # + genelabels,
     heatmap_legend_side = 'right',
     annotation_legend_side = 'top')
dev.off()

# tetst if there are more B cells expressing HLA-DQB2 genes in asthma 
contingency_table <- table(metadata$HLA, metadata$expCond.asthma)

# Asthmatic Non-asthmatic
# HLA-DQB2      2992          1706
# None         15240         12347



# Perform Fisher's exact test
fisher_test_result <- fisher.test(contingency_table)
# Fisher's Exact Test for Count Data
# 
# data:  contingency_table
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.332181 1.515794
# sample estimates:
# odds ratio 
#   1.420865 

# strong association between the gene being expressed and asthmatic
# B cells from asthmatic donors are 1.420865 times more likely to express HLA-DQB2\






##-- cellType B | Under Ig and unstim, what are the de genes for B cells expressing HLA-DQA&B only ----
de.dir <- "HLA-DQA2&DQB2_ig_unstim_MAST/"
cellType <- "B"

rds.sub <- readRDS(paste0(excel.dir, "B_ig_unstim_hla.rds"))

rds.sub@meta.data$HLA <- ifelse(grepl("HLA-DQA2&DQB2", rds.sub@meta.data$HLA), "HLA-DQA2&DQB2", "None")


# rds.sub@meta.data %>% count(HLA)

file.all <- list.files(paste0(excel.dir, de.dir))
# file <- file.all[grepl("DQB2.*_top", file.all, ignore.case = TRUE)]

# [1] "expCondCompDeMarkers_DQA2|DQB2-None_top20_downDe.xlsx"
# [2] "expCondCompDeMarkers_DQA2|DQB2-None_top20_upDe.xlsx"   

# excel_data <- lapply(paste0(excel.dir, de.dir, file), read_excel)
# combined_data <- do.call(rbind, excel_data)

file <- file.all[grepl("*_adjSig_1SelClusters.xlsx", file.all, ignore.case = TRUE)]
excel_data <- read_excel(paste0(excel.dir, de.dir, file))


top_down <- excel_data %>%
  filter(avg_log2FC < 0 & ...1 %in% pcgs) %>%
  arrange(-abs(avg_log2FC)) %>%  
  head(20) %>%
  rename(Gene = ...1) %>%
  mutate(geneType = "B") %>%
  select(c(Gene, geneType))


top_up <- excel_data %>%
  filter(avg_log2FC > 0 & ...1 %in% pcgs) %>%
  arrange(-abs(avg_log2FC)) %>%  
  head(20) %>%
  rename(Gene = ...1) %>%
  mutate(geneType = "B") %>%
  select(c(Gene, geneType))

combined_data <- bind_rows(top_down, top_up)





# library(biomaRt)
# 
# # Connect to Ensembl
# ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")  # Adjust species as needed
# 
# 
# # Query biotypes
# gene_info <- getBM(attributes = c("hgnc_symbol", "gene_biotype"),
#                    filters = "hgnc_symbol",
#                    values = combined_data$Gene,
#                    mart = ensembl)

# skip the line below for now 
# rds.sub <- subset(rds.sub, 
#                   subset = (expCond.donor == "602" | 
#                               expCond.donor == "681"))

geneNames <- c(combined_data$Gene, "HLA-DOA", "HLA-DOB", "IGHD")
features <- geneNames
resDir <- paste0(excel.dir, de.dir)
seuratObjFinal <- rds.sub
# expCondReorderLevels <- c("unstim_602", "Ig_4h_None_602", "Ig_18h_None_602", 
#                           "Ig_4h_DQA2|DQB2_602", "Ig_18h_DQA2|DQB2_602",
#                           "unstim_681", "Ig_4h_None_681","Ig_18h_None_681",
#                           "Ig_4h_DQA2|DQB2_681","Ig_18h_DQA2|DQB2_681")


expCondReorderLevels <- apply(expand.grid(c("unstim_0h", "Ig_4h_None", "Ig_18h_None", 
                                            "Ig_4h_HLA-DQA2&DQB2", "Ig_18h_HLA-DQA2&DQB2"), 
                                          sort(unique(seuratObjFinal@meta.data$expCond.donor))), 
                              1, function(x) paste(x[2], x[1], sep = "_"))

slot <- "scale.data"

# clusterLevels <- levels(Seurat::Idents(seuratObjFinal)) # cell types
# seuratObjFinal <- subset(seuratObjFinal, idents = "B cells")
Seurat::DefaultAssay(seuratObjFinal) <- "RNA"

# by default, both do.scale/center are on, 
# will scale the expression level for each gene feature by dividing the centered feature expression levels by their standard deviations
seuratObjFinal <- Seurat::ScaleData(object = seuratObjFinal, do.scale = T, do.center = T, features = features) 

# 1. input for heatmap
data <- as.data.frame(as.matrix(t(SeuratObject::GetAssayData(object = seuratObjFinal, 
                                                             layer = slot)[features, colnames(x = seuratObjFinal), drop = FALSE]))) 


seuratObjFinal@meta.data$expCond.stimuli.time.HLA <- paste(seuratObjFinal@meta.data$expCond.stimuli.time, 
                                                           seuratObjFinal@meta.data$HLA, sep = "_")

seuratObjFinal@meta.data$expCond.stimuli.time.HLA <- ifelse(grepl("unstim_0h", seuratObjFinal@meta.data$expCond.stimuli.time.HLA), 
                                                            "unstim_0h", 
                                                            seuratObjFinal@meta.data$expCond.stimuli.time.HLA)

seuratObjFinal@meta.data$expCond.stimuli.time.HLA.asthma <- paste(seuratObjFinal@meta.data$expCond.stimuli.time.HLA, 
                                                                  seuratObjFinal@meta.data$expCond.asthma, sep = "_")

metadata <- seuratObjFinal@meta.data[rownames(data), ]



# focus on genes of interest
geneNames = c(combined_data$Gene, "HLA-DOA", "HLA-DOB", "IGHD")
data <- t(data[,geneNames])

data[data < -2.5] <- -2.5
data[data > 2.5] <- 2.5


# 2. set the color scheme
# myCol <- colorRampPalette(c("#31d3ff", "black", "#fff200"))(100)
# myCol <- colorRampPalette(c( "#003049", "#669bbc", "#fdf0d5", "#c1121f", "#780000"))(100)
myCol <- colorRampPalette(c("#3B4CC0", "#FFFFFF", "#B40426"))(100)

# myCol <- colorRampPalette(brewer.pal(5, "RdYlBu"))(100)

myBreaks <- seq(-2.5, 2.5, length.out = 100)

col_ann <- data.frame(
  # Donor = factor(metadata$expCond.donor,levels=sort(unique(metadata$expCond.donor))),
  Asthma = factor(metadata$expCond.asthma,levels=sort(unique(metadata$expCond.asthma))),
  Stimuli.Time.HLA = factor(metadata$expCond.stimuli.time.HLA,
                            levels=c('unstim_0h',
                                     'Ig_4h_None',
                                     'Ig_18h_None',
                                     'Ig_4h_HLA-DQA2&DQB2',
                                     'Ig_18h_HLA-DQA2&DQB2')),
  stringsAsFactors = FALSE) #do not have to set factor 

row_ann <- data.frame(
  # Genes = rownames(data),
  Expression = c(rep('Down-regulated',20), rep('Up-regulated',20),
                 rep('Genes of interest',3)),
  stringsAsFactors = FALSE)

donors <- sort(unique(metadata$expCond.donor))

# Generate a color gradient
gradient_colors <- colorRampPalette(c("#f2a65a", "#313715"))(length(donors))

# Create the named vector
donor_colors <- setNames(gradient_colors, donors)

col_colours <- list(
  # Donor = donor_colors,
  Asthma = c("Asthmatic" = "#471323",
             "Non-asthmatic" = "#CEB992"),
  Stimuli.Time.HLA = c('unstim_0h' = '#eff9f0', 
                       'Ig_4h_None' = '#FFBCB5', 
                       'Ig_18h_None' = '#03978E', 
                       'Ig_4h_HLA-DQA2&DQB2' = '#e29578', 
                       'Ig_18h_HLA-DQA2&DQB2' = '#006d77'))


row_colours <- list(
  Expression = c('Down-regulated' = '#6C8FC8',
                 'Up-regulated' = '#EE4687',
                 'Genes of interest' = '#f1faee'))


colAnn <- HeatmapAnnotation(
  df = col_ann,
  which = 'col', # 'col' (samples) or 'row' (gene) annotation?
  na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
  col = col_colours,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  simple_anno_size = unit(3, "mm"),
  gap = unit(1, 'mm'),
  annotation_legend_param = list(
    # Donor = list(
    #   nrow = 4, # number of rows across which the legend will be arranged
    #   title = 'Donor',
    #   title_position = 'topcenter',
    #   legend_direction = 'vertical',
    #   title_gp = gpar(fontsize = 11),  # legend on the side
    #   labels_gp = gpar(fontsize = 11)),  # legend on the side
    Asthma = list(
      nrow = 2, # number of rows across which the legend will be arranged
      title = 'Asthma',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 11),  # legend on the side
      labels_gp = gpar(fontsize = 11)),  # legend on the side
    Stimuli.Time.HLA = list(
      nrow = 3,
      title = 'Stimuli.Time.HLA',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 11), # legend on the side
      labels_gp = gpar(fontsize = 11)))) # legend on the side


rowAnn <- HeatmapAnnotation(
  df = row_ann,
  which = 'row', # 'col' (samples) or 'row' (gene) annotation?
  na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
  col = row_colours,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  simple_anno_size = unit(3, "mm"),
  gap = unit(1, 'mm'),
  annotation_legend_param = list(
    Expression = list(
      nrow = 3, 
      title = 'Expression w/\n HLA-DQA2&DQB2',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 11),  # legend on the side
      labels_gp = gpar(fontsize = 11))  # legend on the side
  ))



hmap <- Heatmap(data,
                # split the genes / rows according to the PAM clusters
                # column_split = factor(metadata$expCond.stimuli.time.HLA.donor,
                #                       levels = expCondReorderLevels),
                # column_split = factor(metadata$expCond.asthma,
                #                       levels = c("Asthmatic", "Non-asthmatic")),
                column_split = factor(metadata$expCond.stimuli.time.HLA.asthma,
                                      levels = apply(expand.grid(c("unstim_0h", "Ig_4h_None", "Ig_18h_None",
                                                                   "Ig_4h_HLA-DQA2&DQB2", "Ig_18h_HLA-DQA2&DQB2"),
                                                                 sort(unique(seuratObjFinal@meta.data$expCond.asthma))),
                                                     1, function(x) paste(x[1], x[2], sep = "_"))),
                # column_split = factor(metadata$expCond.stimuli.time.HLA,
                #                       levels = c("unstim", "Ig_4h_None", "Ig_18h_None",
                #                                  "Ig_4h_HLA-DQA2&DQB2", "Ig_18h_HLA-DQA2&DQB2")),
                row_split = factor(c(rep('Down-regulated',20), rep('Up-regulated',20),
                                     rep('Genes of interest',3)),
                                   levels=c('Up-regulated','Down-regulated','Genes of interest')),
                
                # cluster_row_slices = FALSE,
                name = 'Scaled\ngene\nexpression',
                
                col = colorRamp2(myBreaks, myCol),
                use_raster = FALSE,
                # parameters for the colour-bar that represents gradient of expression
                heatmap_legend_param = list(
                  at = c(-2.5, -2, -1, 0, 1, 2, 2.5),
                  color_bar = 'continuous',
                  legend_direction = 'vertical',
                  legend_width = unit(6, 'cm'),
                  legend_height = unit(4, 'cm'),
                  title_position = 'topcenter',
                  title_gp=gpar(fontsize = 11),
                  labels_gp=gpar(fontsize = 11)),
                
                # row (gene) parameters
                cluster_rows = FALSE,
                show_row_dend = FALSE,
                row_title = 'Differentially expressed genes',
                row_title_side = 'left',
                row_title_gp = gpar(fontsize = 11),
                row_title_rot = 90,
                show_row_names = TRUE,
                row_names_gp = gpar(fontsize = 9), # y-axis ticks
                row_names_side = 'left',
                # row_dend_width = unit(25,'mm'),
                
                # column (sample) parameters
                cluster_columns = FALSE,
                show_column_dend = FALSE,
                # column_title = '',
                column_title_side = 'top',
                column_title_gp = gpar(fontsize = 0), # column ticks
                column_title_rot = 90,
                show_column_names = FALSE,
                column_names_gp = gpar(fontsize = 7),
                column_names_max_height = unit(10, 'cm'),
                # column_dend_height = unit(25,'mm')
                # specify top and bottom annotations
                top_annotation = colAnn,
                left_annotation = rowAnn)


pdf(paste0(excel.dir, de.dir, cellType, "_complexheatmap.pdf"), height=6.67, width=7.28)
draw(hmap, 
     # + genelabels,
     heatmap_legend_side = 'right',
     annotation_legend_side = 'top')
dev.off()

# tetst if there are more B cells expressing HLA-DQB2 genes in asthma 
contingency_table <- table(metadata$HLA, metadata$expCond.asthma)

# Asthmatic Non-asthmatic
# HLA-DQA2&DQB2      1987           587
# None              16245         13466


# Perform Fisher's exact test
fisher_test_result <- fisher.test(contingency_table)

# Fisher's Exact Test for Count Data
# 
# data:  contingency_table
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  2.550525 3.090214
# sample estimates:
# odds ratio 
#   2.805662 


# strong association between the gene being expressed and asthmatic
# B cells from asthmatic donors are 2.805662 times more likely to express HLA-DQB2\







##-- cellType B | Under Ig and unstim, what are the de genes for asthmatic vs non-asthmatic----
de.dir <- "asthma_ig_unstim_MAST/"
cellType <- "B"

rds.sub <- readRDS(paste0(excel.dir, "B_ig_unstim_hla.rds"))
rds.sub@meta.data$HLA_tmp <- rds.sub@meta.data$HLA

rds.sub@meta.data$HLA <- ifelse(grepl("None", rds.sub@meta.data$HLA), "None", 'HLA-DQA2|DQB2')

file.all <- list.files(paste0(excel.dir, de.dir))
file <- file.all[grepl("*_adjSig_1SelClusters.xlsx", file.all, ignore.case = TRUE)]

# expCondCompDeMarkers_Asthmatic-Non-asthmatic_adjSig_1SelClusters.xlsx

# excel_data <- lapply(paste0(excel.dir, de.dir, file), read_excel)
# combined_data <- do.call(rbind, excel_data)

excel_data <- read_excel(paste0(excel.dir, de.dir, file))


top_down <- excel_data %>%
  filter(avg_log2FC < 0 & ...1 %in% pcgs) %>%
  arrange(-abs(avg_log2FC)) %>%  
  head(20) %>%
  rename(Gene = ...1) %>%
  mutate(geneType = "B") %>%
  select(c(Gene, geneType))


top_up <- excel_data %>%
  filter(avg_log2FC > 0 & ...1 %in% pcgs) %>%
  arrange(-abs(avg_log2FC)) %>%  
  head(20) %>%
  rename(Gene = ...1) %>%
  mutate(geneType = "B") %>%
  select(c(Gene, geneType))

combined_data <- bind_rows(top_down, top_up)


geneNames <- c(combined_data$Gene, "HLA-DOA", "HLA-DOB", "IGHD")
# geneNames <- combined_data$Gene

features <- geneNames
resDir <- paste0(excel.dir, de.dir)
seuratObjFinal <- rds.sub
# expCondReorderLevels <- c("unstim_602", "Ig_4h_None_602", "Ig_18h_None_602", 
#                           "Ig_4h_DQA2|DQB2_602", "Ig_18h_DQA2|DQB2_602",
#                           "unstim_681", "Ig_4h_None_681","Ig_18h_None_681",
#                           "Ig_4h_DQA2|DQB2_681","Ig_18h_DQA2|DQB2_681")


expCondReorderLevels <- apply(expand.grid(c("unstim_0h", "Ig_4h", "Ig_18h"), 
                                          sort(unique(seuratObjFinal@meta.data$expCond.donor))), 1, function(x) paste(x[2], x[1], sep = "_"))

slot <- "scale.data"

# clusterLevels <- levels(Seurat::Idents(seuratObjFinal)) # cell types
# seuratObjFinal <- subset(seuratObjFinal, idents = "B cells")
Seurat::DefaultAssay(seuratObjFinal) <- "RNA"

# by default, both do.scale/center are on, 
# will scale the expression level for each gene feature by dividing the centered feature expression levels by their standard deviations
seuratObjFinal <- Seurat::ScaleData(object = seuratObjFinal, do.scale = T, do.center = T, features = features) 

# 1. input for heatmap
data <- as.data.frame(as.matrix(t(SeuratObject::GetAssayData(object = seuratObjFinal, 
                                                             layer = slot)[features, colnames(x = seuratObjFinal), drop = FALSE]))) 


seuratObjFinal@meta.data$expCond.stimuli.time.asthma <- paste(seuratObjFinal@meta.data$expCond.stimuli.time, 
                                                              seuratObjFinal@meta.data$expCond.asthma, sep = "_")

seuratObjFinal@meta.data$expCond.stimuli.time.asthma.hla <- paste(seuratObjFinal@meta.data$expCond.stimuli.time.asthma, 
                                                                  seuratObjFinal@meta.data$HLA, sep = "_")

metadata <- seuratObjFinal@meta.data[rownames(data), ]

# focus on genes of interest
data <- t(data[,geneNames])

data[data < -2.5] <- -2.5
data[data > 2.5] <- 2.5


# 2. set the color scheme
# myCol <- colorRampPalette(c("#31d3ff", "black", "#fff200"))(100)
# myCol <- colorRampPalette(c( "#003049", "#669bbc", "#fdf0d5", "#c1121f", "#780000"))(100)
myCol <- colorRampPalette(c("#3B4CC0", "#FFFFFF", "#B40426"))(100)

# myCol <- colorRampPalette(brewer.pal(5, "RdYlBu"))(100)

myBreaks <- seq(-2.5, 2.5, length.out = 100)

col_ann <- data.frame(
  # Donor = factor(metadata$expCond.donor,levels=sort(unique(metadata$expCond.donor))),
  Asthma = factor(metadata$expCond.asthma,levels=sort(unique(metadata$expCond.asthma))),
  Stimuli.Time = factor(metadata$expCond.stimuli.time,
                        levels=c('unstim_0h',
                                 'Ig_4h',
                                 'Ig_18h')),
  HLA = factor(metadata$HLA,
                        levels=c("None",
                                 "HLA-DQA2|DQB2")),
  stringsAsFactors = FALSE) #do not have to set factor 

row_ann <- data.frame(
  # Genes = rownames(data),
  Expression = c(rep('Down-regulated',20), rep('Up-regulated',20),
                 rep('Genes of interest',3)),
  stringsAsFactors = FALSE)

donors <- sort(unique(metadata$expCond.donor))

# Generate a color gradient
gradient_colors <- colorRampPalette(c("#f2a65a", "#313715"))(length(donors))

# Create the named vector
donor_colors <- setNames(gradient_colors, donors)

col_colours <- list(
  # Donor = donor_colors,
  Asthma = c("Asthmatic" = "#471323",
             "Non-asthmatic" = "#CEB992"),
  Stimuli.Time = c('unstim_0h' = '#eff9f0', 
                   'Ig_4h' = '#FFBCB5', 
                   'Ig_18h' = '#03978E'),
  HLA = c("None" = "#bcb8b1",
          "HLA-DQA2|DQB2" = "#ffc300"))


row_colours <- list(
  Expression = c('Down-regulated' = '#6C8FC8',
                 'Up-regulated' = '#EE4687',
                 'Genes of interest' = '#f1faee'))


colAnn <- HeatmapAnnotation(
  df = col_ann,
  which = 'col', # 'col' (samples) or 'row' (gene) annotation?
  na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
  col = col_colours,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  simple_anno_size = unit(3, "mm"),
  gap = unit(1, 'mm'),
  annotation_legend_param = list(
    # Donor = list(
    #   nrow = 4, # number of rows across which the legend will be arranged
    #   title = 'Donor',
    #   title_position = 'topcenter',
    #   legend_direction = 'vertical',
    #   title_gp = gpar(fontsize = 11),  # legend on the side
    #   labels_gp = gpar(fontsize = 11)),  # legend on the side
    Asthma = list(
      nrow = 2, # number of rows across which the legend will be arranged
      title = 'Asthma',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 11),  # legend on the side
      labels_gp = gpar(fontsize = 11)),  # legend on the side
    Stimuli.Time = list(
      nrow = 3,
      title = 'Stimuli.Time',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 11), # legend on the side
      labels_gp = gpar(fontsize = 11)), # legend on the side
    HLA = list(
      nrow = 2,
      title = 'HLA',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 11), # legend on the side
      labels_gp = gpar(fontsize = 11))) # legend on the side
  ) 


rowAnn <- HeatmapAnnotation(
  df = row_ann,
  which = 'row', # 'col' (samples) or 'row' (gene) annotation?
  na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
  col = row_colours,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  simple_anno_size = unit(3, "mm"),
  gap = unit(1, 'mm'),
  annotation_legend_param = list(
    Expression = list(
      nrow = 3, 
      title = 'Expression',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 11),  # legend on the side
      labels_gp = gpar(fontsize = 11))  # legend on the side
  ))



hmap <- Heatmap(data,
                # split the genes / rows according to the PAM clusters
                # column_split = factor(metadata$expCond.stimuli.time.HLA.donor,
                #                       levels = expCondReorderLevels),
                # column_split = factor(metadata$expCond.asthma,
                #                       levels = c("Asthmatic", "Non-asthmatic")),
                # column_split = factor(metadata$expCond.stimuli.time.asthma,
                #                       levels = apply(expand.grid(c("unstim_0h", "Ig_4h", "Ig_18h"),
                #                                                  sort(unique(seuratObjFinal@meta.data$expCond.asthma))),
                #                                      1, function(x) paste(x[1], x[2], sep = "_"))),
                # column_split = factor(metadata$expCond.stimuli.time.HLA,
                #                       levels = c("unstim", "Ig_4h_None", "Ig_18h_None",
                #                                  "Ig_4h_HLA-DQA2|DQB2", "Ig_18h_HLA-DQA2|DQB2")),
                column_split = factor(metadata$expCond.stimuli.time.asthma.hla,
                                      levels = c("unstim_0h_Asthmatic_None", 
                                                 "unstim_0h_Asthmatic_HLA-DQA2|DQB2",
                                                 "Ig_4h_Asthmatic_None", 
                                                 "Ig_4h_Asthmatic_HLA-DQA2|DQB2",
                                                 "Ig_18h_Asthmatic_None",
                                                 "Ig_18h_Asthmatic_HLA-DQA2|DQB2",
                                                 "unstim_0h_Non-asthmatic_None", 
                                                 "unstim_0h_Non-asthmatic_HLA-DQA2|DQB2",
                                                 "Ig_4h_Non-asthmatic_None", 
                                                 "Ig_4h_Non-asthmatic_HLA-DQA2|DQB2",
                                                 "Ig_18h_Non-asthmatic_None",
                                                 "Ig_18h_Non-asthmatic_HLA-DQA2|DQB2"
                                                 )),
                row_split = factor(c(rep('Down-regulated',20), rep('Up-regulated',20),
                                     rep('Genes of interest',3)),
                                   levels=c('Up-regulated','Down-regulated','Genes of interest')),
                
                # cluster_row_slices = FALSE,
                name = 'Scaled\ngene\nexpression',
                
                col = colorRamp2(myBreaks, myCol),
                use_raster = FALSE,
                # parameters for the colour-bar that represents gradient of expression
                heatmap_legend_param = list(
                  at = c(-2.5, -2, -1, 0, 1, 2, 2.5),
                  color_bar = 'continuous',
                  legend_direction = 'vertical',
                  legend_width = unit(6, 'cm'),
                  legend_height = unit(4, 'cm'),
                  title_position = 'topcenter',
                  title_gp=gpar(fontsize = 11),
                  labels_gp=gpar(fontsize = 11)),
                
                # row (gene) parameters
                cluster_rows = FALSE,
                show_row_dend = FALSE,
                row_title = 'Differentially expressed genes',
                row_title_side = 'left',
                row_title_gp = gpar(fontsize = 11),
                row_title_rot = 90,
                show_row_names = TRUE,
                row_names_gp = gpar(fontsize = 9), # y-axis ticks
                row_names_side = 'left',
                # row_dend_width = unit(25,'mm'),
                
                # column (sample) parameters
                cluster_columns = FALSE,
                show_column_dend = FALSE,
                # column_title = '',
                column_title_side = 'top',
                column_title_gp = gpar(fontsize = 0), # column ticks
                column_title_rot = 90,
                show_column_names = FALSE,
                column_names_gp = gpar(fontsize = 7),
                column_names_max_height = unit(10, 'cm'),
                # column_dend_height = unit(25,'mm')
                # specify top and bottom annotations
                top_annotation = colAnn,
                left_annotation = rowAnn)


pdf(paste0(excel.dir, de.dir, cellType, "_complexheatmap.pdf"), height=6.67, width=7.28)
draw(hmap, 
     # + genelabels,
     heatmap_legend_side = 'right',
     annotation_legend_side = 'top')
dev.off()

# # tetst if there are more B cells expressing HLA-DQ2 genes in asthma 
# contingency_table <- table(metadata$HLA, metadata$expCond.asthma)
# 
# #                 Asthmatic Non-asthmatic
# # HLA-DQA2|DQB2     11184          4376
# # None               7048          9677
# 
# 
# 
# # Perform Fisher's exact test
# fisher_test_result <- fisher.test(contingency_table)
# # Fisher's Exact Test for Count Data
# # 
# # data:  contingency_table
# # p-value < 2.2e-16
# # alternative hypothesis: true odds ratio is not equal to 1
# # 95 percent confidence interval:
# #  3.348830 3.677406
# # sample estimates:
# # odds ratio 
# #   3.508898 
# 
# # strong association between the gene being expressed and asthmatic
# # B cells from asthmatic donors are 3.51 times more likely to express HLA-DQA2/DQB2
# 
# 
# metadata_unstim<- metadata %>%
#   filter(expCond.stimuli.time %in% "unstim_0h")
# contingency_table <- table(metadata_unstim$HLA, metadata_unstim$expCond.asthma)
# res <- fisher.test(contingency_table) #p-value < 2.2e-16, odds ratio 3.810912
# 
# metadata_ig_4h <- metadata %>%
#   filter(expCond.stimuli.time %in% "Ig_4h")
# contingency_table <- table(metadata_ig_4h$HLA, metadata_ig_4h$expCond.asthma)
# res <- fisher.test(contingency_table) #p-value < 2.2e-16, odds ratio 2.062551 






##-- cellType B | Under unstim, what are the de genes for asthmatic vs non-asthmatic----
de.dir <- "asthma_unstim_MAST/"
cellType <- "B"

rds.sub <- readRDS(paste0(excel.dir, "B_ig_unstim_hla.rds"))

file.all <- list.files(paste0(excel.dir, de.dir))
file <- file.all[grepl("*_adjSig_1SelClusters.xlsx", file.all, ignore.case = TRUE)]

# expCondCompDeMarkers_Asthmatic-Non-asthmatic_adjSig_1SelClusters.xlsx

# excel_data <- lapply(paste0(excel.dir, de.dir, file), read_excel)
# combined_data <- do.call(rbind, excel_data)

excel_data <- read_excel(paste0(excel.dir, de.dir, file))


top_down <- excel_data %>%
  filter(avg_log2FC < 0 & ...1 %in% pcgs) %>%
  arrange(-abs(avg_log2FC)) %>%  
  head(20) %>%
  rename(Gene = ...1) %>%
  mutate(geneType = "B") %>%
  select(c(Gene, geneType))


top_up <- excel_data %>%
  filter(avg_log2FC > 0 & ...1 %in% pcgs) %>%
  arrange(-abs(avg_log2FC)) %>%  
  head(20) %>%
  rename(Gene = ...1) %>%
  mutate(geneType = "B") %>%
  select(c(Gene, geneType))

combined_data <- bind_rows(top_down, top_up)


geneNames <- c(combined_data$Gene, "HLA-DOA", "HLA-DOB", "IGHD")
features <- geneNames
resDir <- paste0(excel.dir, de.dir)
seuratObjFinal <- rds.sub
# expCondReorderLevels <- c("unstim_602", "Ig_4h_None_602", "Ig_18h_None_602", 
#                           "Ig_4h_DQA2|DQB2_602", "Ig_18h_DQA2|DQB2_602",
#                           "unstim_681", "Ig_4h_None_681","Ig_18h_None_681",
#                           "Ig_4h_DQA2|DQB2_681","Ig_18h_DQA2|DQB2_681")


expCondReorderLevels <- apply(expand.grid(c("unstim_0h", "Ig_4h", "Ig_18h"), 
                                          sort(unique(seuratObjFinal@meta.data$expCond.donor))), 1, function(x) paste(x[2], x[1], sep = "_"))

slot <- "scale.data"

# clusterLevels <- levels(Seurat::Idents(seuratObjFinal)) # cell types
# seuratObjFinal <- subset(seuratObjFinal, idents = "B cells")
Seurat::DefaultAssay(seuratObjFinal) <- "RNA"

# by default, both do.scale/center are on, 
# will scale the expression level for each gene feature by dividing the centered feature expression levels by their standard deviations
seuratObjFinal <- Seurat::ScaleData(object = seuratObjFinal, do.scale = T, do.center = T, features = features) 

# 1. input for heatmap
data <- as.data.frame(as.matrix(t(SeuratObject::GetAssayData(object = seuratObjFinal, 
                                                             layer = slot)[features, colnames(x = seuratObjFinal), drop = FALSE]))) 


seuratObjFinal@meta.data$expCond.stimuli.time.asthma <- paste(rds.sub@meta.data$expCond.stimuli.time, 
                                                              rds.sub@meta.data$expCond.asthma, sep = "_")

metadata <- seuratObjFinal@meta.data[rownames(data), ]

# focus on genes of interest
geneNames = c(combined_data$Gene, "HLA-DOA", "HLA-DOB", "IGHD")
data <- t(data[,geneNames])

data[data < -2.5] <- -2.5
data[data > 2.5] <- 2.5


# 2. set the color scheme
# myCol <- colorRampPalette(c("#31d3ff", "black", "#fff200"))(100)
# myCol <- colorRampPalette(c( "#003049", "#669bbc", "#fdf0d5", "#c1121f", "#780000"))(100)
myCol <- colorRampPalette(c("#3B4CC0", "#FFFFFF", "#B40426"))(100)

# myCol <- colorRampPalette(brewer.pal(5, "RdYlBu"))(100)

myBreaks <- seq(-2.5, 2.5, length.out = 100)

col_ann <- data.frame(
  # Donor = factor(metadata$expCond.donor,levels=sort(unique(metadata$expCond.donor))),
  Asthma = factor(metadata$expCond.asthma,levels=sort(unique(metadata$expCond.asthma))),
  Stimuli.Time = factor(metadata$expCond.stimuli.time,
                        levels=c('unstim_0h',
                                 'Ig_4h',
                                 'Ig_18h')),
  stringsAsFactors = FALSE) #do not have to set factor 

row_ann <- data.frame(
  # Genes = rownames(data),
  Expression = c(rep('Down-regulated',20), rep('Up-regulated',20),
                 rep('Genes of interest',3)),
  stringsAsFactors = FALSE)

donors <- sort(unique(metadata$expCond.donor))

# Generate a color gradient
gradient_colors <- colorRampPalette(c("#f2a65a", "#313715"))(length(donors))

# Create the named vector
donor_colors <- setNames(gradient_colors, donors)

col_colours <- list(
  # Donor = donor_colors,
  Asthma = c("Asthmatic" = "#471323",
             "Non-asthmatic" = "#CEB992"),
  Stimuli.Time = c('unstim_0h' = '#eff9f0', 
                   'Ig_4h' = '#FFBCB5', 
                   'Ig_18h' = '#03978E'))


row_colours <- list(
  Expression = c('Down-regulated' = '#6C8FC8',
                 'Up-regulated' = '#EE4687',
                 'Genes of interest' = '#f1faee'))


colAnn <- HeatmapAnnotation(
  df = col_ann,
  which = 'col', # 'col' (samples) or 'row' (gene) annotation?
  na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
  col = col_colours,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  simple_anno_size = unit(3, "mm"),
  gap = unit(1, 'mm'),
  annotation_legend_param = list(
    # Donor = list(
    #   nrow = 4, # number of rows across which the legend will be arranged
    #   title = 'Donor',
    #   title_position = 'topcenter',
    #   legend_direction = 'vertical',
    #   title_gp = gpar(fontsize = 11),  # legend on the side
    #   labels_gp = gpar(fontsize = 11)),  # legend on the side
    Asthma = list(
      nrow = 2, # number of rows across which the legend will be arranged
      title = 'Asthma',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 11),  # legend on the side
      labels_gp = gpar(fontsize = 11)),  # legend on the side
    Stimuli.Time = list(
      nrow = 3,
      title = 'Stimuli.Time',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 11), # legend on the side
      labels_gp = gpar(fontsize = 11)))) # legend on the side


rowAnn <- HeatmapAnnotation(
  df = row_ann,
  which = 'row', # 'col' (samples) or 'row' (gene) annotation?
  na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
  col = row_colours,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  simple_anno_size = unit(3, "mm"),
  gap = unit(1, 'mm'),
  annotation_legend_param = list(
    Expression = list(
      nrow = 3, 
      title = 'Expression',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 11),  # legend on the side
      labels_gp = gpar(fontsize = 11))  # legend on the side
  ))



hmap <- Heatmap(data,
                # split the genes / rows according to the PAM clusters
                # column_split = factor(metadata$expCond.stimuli.time.HLA.donor,
                #                       levels = expCondReorderLevels),
                # column_split = factor(metadata$expCond.asthma,
                #                       levels = c("Asthmatic", "Non-asthmatic")),
                column_split = factor(metadata$expCond.stimuli.time.asthma,
                                      levels = apply(expand.grid(c("unstim_0h", "Ig_4h", "Ig_18h"),
                                                                 sort(unique(seuratObjFinal@meta.data$expCond.asthma))),
                                                     1, function(x) paste(x[1], x[2], sep = "_"))),
                # column_split = factor(metadata$expCond.stimuli.time.HLA,
                #                       levels = c("unstim", "Ig_4h_None", "Ig_18h_None",
                #                                  "Ig_4h_HLA-DQA2|DQB2", "Ig_18h_HLA-DQA2|DQB2")),
                row_split = factor(c(rep('Down-regulated',20), rep('Up-regulated',20),
                                     rep('Genes of interest',3)),
                                   levels=c('Up-regulated','Down-regulated','Genes of interest')),
                
                # cluster_row_slices = FALSE,
                name = 'Scaled\ngene\nexpression',
                
                col = colorRamp2(myBreaks, myCol),
                use_raster = FALSE,
                # parameters for the colour-bar that represents gradient of expression
                heatmap_legend_param = list(
                  at = c(-2.5, -2, -1, 0, 1, 2, 2.5),
                  color_bar = 'continuous',
                  legend_direction = 'vertical',
                  legend_width = unit(6, 'cm'),
                  legend_height = unit(4, 'cm'),
                  title_position = 'topcenter',
                  title_gp=gpar(fontsize = 11),
                  labels_gp=gpar(fontsize = 11)),
                
                # row (gene) parameters
                cluster_rows = FALSE,
                show_row_dend = FALSE,
                row_title = 'Differentially expressed genes',
                row_title_side = 'left',
                row_title_gp = gpar(fontsize = 11),
                row_title_rot = 90,
                show_row_names = TRUE,
                row_names_gp = gpar(fontsize = 9), # y-axis ticks
                row_names_side = 'left',
                # row_dend_width = unit(25,'mm'),
                
                # column (sample) parameters
                cluster_columns = FALSE,
                show_column_dend = FALSE,
                # column_title = '',
                column_title_side = 'top',
                column_title_gp = gpar(fontsize = 0), # column ticks
                column_title_rot = 90,
                show_column_names = FALSE,
                column_names_gp = gpar(fontsize = 7),
                column_names_max_height = unit(10, 'cm'),
                # column_dend_height = unit(25,'mm')
                # specify top and bottom annotations
                top_annotation = colAnn,
                left_annotation = rowAnn)


pdf(paste0(excel.dir, de.dir, cellType, "_complexheatmap.pdf"), height=6.67, width=7.28)
draw(hmap, 
     # + genelabels,
     heatmap_legend_side = 'right',
     annotation_legend_side = 'top')
dev.off()

# # tetst if there are more B cells expressing HLA-DQ2 genes in asthma 
# contingency_table <- table(metadata$HLA, metadata$expCond.asthma)
# 
# #                 Asthmatic Non-asthmatic
# # HLA-DQA2|DQB2     11184          4376
# # None               7048          9677
# 
# 
# 
# # Perform Fisher's exact test
# fisher_test_result <- fisher.test(contingency_table)
# # Fisher's Exact Test for Count Data
# # 
# # data:  contingency_table
# # p-value < 2.2e-16
# # alternative hypothesis: true odds ratio is not equal to 1
# # 95 percent confidence interval:
# #  3.348830 3.677406
# # sample estimates:
# # odds ratio 
# #   3.508898 
# 
# # strong association between the gene being expressed and asthmatic
# # B cells from asthmatic donors are 3.51 times more likely to express HLA-DQA2/DQB2
# 
# 
# metadata_unstim<- metadata %>%
#   filter(expCond.stimuli.time %in% "unstim_0h")
# contingency_table <- table(metadata_unstim$HLA, metadata_unstim$expCond.asthma)
# res <- fisher.test(contingency_table) #p-value < 2.2e-16, odds ratio 3.810912
# 
# metadata_ig_4h <- metadata %>%
#   filter(expCond.stimuli.time %in% "Ig_4h")
# contingency_table <- table(metadata_ig_4h$HLA, metadata_ig_4h$expCond.asthma)
# res <- fisher.test(contingency_table) #p-value < 2.2e-16, odds ratio 2.062551 







##-- cellType B | Under Ig_4h, what are the de genes for asthmatic vs non-asthmatic----
de.dir <- "asthma_ig_4h_MAST/"
cellType <- "B"

rds.sub <- readRDS(paste0(excel.dir, "B_ig_unstim_hla.rds"))

file.all <- list.files(paste0(excel.dir, de.dir))
file <- file.all[grepl("*_adjSig_1SelClusters.xlsx", file.all, ignore.case = TRUE)]

# expCondCompDeMarkers_Asthmatic-Non-asthmatic_adjSig_1SelClusters.xlsx

# excel_data <- lapply(paste0(excel.dir, de.dir, file), read_excel)
# combined_data <- do.call(rbind, excel_data)

excel_data <- read_excel(paste0(excel.dir, de.dir, file))


top_down <- excel_data %>%
  filter(avg_log2FC < 0 & ...1 %in% pcgs) %>%
  arrange(-abs(avg_log2FC)) %>%  
  head(20) %>%
  rename(Gene = ...1) %>%
  mutate(geneType = "B") %>%
  select(c(Gene, geneType))


top_up <- excel_data %>%
  filter(avg_log2FC > 0 & ...1 %in% pcgs) %>%
  arrange(-abs(avg_log2FC)) %>%  
  head(20) %>%
  rename(Gene = ...1) %>%
  mutate(geneType = "B") %>%
  select(c(Gene, geneType))

combined_data <- bind_rows(top_down, top_up)


geneNames <- c(combined_data$Gene, "HLA-DOA", "HLA-DOB", "IGHD")
features <- geneNames
resDir <- paste0(excel.dir, de.dir)
seuratObjFinal <- rds.sub
# expCondReorderLevels <- c("unstim_602", "Ig_4h_None_602", "Ig_18h_None_602", 
#                           "Ig_4h_DQA2|DQB2_602", "Ig_18h_DQA2|DQB2_602",
#                           "unstim_681", "Ig_4h_None_681","Ig_18h_None_681",
#                           "Ig_4h_DQA2|DQB2_681","Ig_18h_DQA2|DQB2_681")


expCondReorderLevels <- apply(expand.grid(c("unstim_0h", "Ig_4h", "Ig_18h"), 
                                          sort(unique(seuratObjFinal@meta.data$expCond.donor))), 1, function(x) paste(x[2], x[1], sep = "_"))

slot <- "scale.data"

# clusterLevels <- levels(Seurat::Idents(seuratObjFinal)) # cell types
# seuratObjFinal <- subset(seuratObjFinal, idents = "B cells")
Seurat::DefaultAssay(seuratObjFinal) <- "RNA"

# by default, both do.scale/center are on, 
# will scale the expression level for each gene feature by dividing the centered feature expression levels by their standard deviations
seuratObjFinal <- Seurat::ScaleData(object = seuratObjFinal, do.scale = T, do.center = T, features = features) 

# 1. input for heatmap
data <- as.data.frame(as.matrix(t(SeuratObject::GetAssayData(object = seuratObjFinal, 
                                                             layer = slot)[features, colnames(x = seuratObjFinal), drop = FALSE]))) 


seuratObjFinal@meta.data$expCond.stimuli.time.asthma <- paste(rds.sub@meta.data$expCond.stimuli.time, 
                                                              rds.sub@meta.data$expCond.asthma, sep = "_")

metadata <- seuratObjFinal@meta.data[rownames(data), ]

# focus on genes of interest
geneNames = c(combined_data$Gene, "HLA-DOA", "HLA-DOB", "IGHD")
data <- t(data[,geneNames])

data[data < -2.5] <- -2.5
data[data > 2.5] <- 2.5


# 2. set the color scheme
# myCol <- colorRampPalette(c("#31d3ff", "black", "#fff200"))(100)
# myCol <- colorRampPalette(c( "#003049", "#669bbc", "#fdf0d5", "#c1121f", "#780000"))(100)
myCol <- colorRampPalette(c("#3B4CC0", "#FFFFFF", "#B40426"))(100)

# myCol <- colorRampPalette(brewer.pal(5, "RdYlBu"))(100)

myBreaks <- seq(-2.5, 2.5, length.out = 100)

col_ann <- data.frame(
  # Donor = factor(metadata$expCond.donor,levels=sort(unique(metadata$expCond.donor))),
  Asthma = factor(metadata$expCond.asthma,levels=sort(unique(metadata$expCond.asthma))),
  Stimuli.Time = factor(metadata$expCond.stimuli.time,
                        levels=c('unstim_0h',
                                 'Ig_4h',
                                 'Ig_18h')),
  stringsAsFactors = FALSE) #do not have to set factor 

row_ann <- data.frame(
  # Genes = rownames(data),
  Expression = c(rep('Down-regulated',20), rep('Up-regulated',20),
                 rep('Genes of interest',3)),
  stringsAsFactors = FALSE)

donors <- sort(unique(metadata$expCond.donor))

# Generate a color gradient
gradient_colors <- colorRampPalette(c("#f2a65a", "#313715"))(length(donors))

# Create the named vector
donor_colors <- setNames(gradient_colors, donors)

col_colours <- list(
  # Donor = donor_colors,
  Asthma = c("Asthmatic" = "#471323",
             "Non-asthmatic" = "#CEB992"),
  Stimuli.Time = c('unstim_0h' = '#eff9f0', 
                   'Ig_4h' = '#FFBCB5', 
                   'Ig_18h' = '#03978E'))


row_colours <- list(
  Expression = c('Down-regulated' = '#6C8FC8',
                 'Up-regulated' = '#EE4687',
                 'Genes of interest' = '#f1faee'))


colAnn <- HeatmapAnnotation(
  df = col_ann,
  which = 'col', # 'col' (samples) or 'row' (gene) annotation?
  na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
  col = col_colours,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  simple_anno_size = unit(3, "mm"),
  gap = unit(1, 'mm'),
  annotation_legend_param = list(
    # Donor = list(
    #   nrow = 4, # number of rows across which the legend will be arranged
    #   title = 'Donor',
    #   title_position = 'topcenter',
    #   legend_direction = 'vertical',
    #   title_gp = gpar(fontsize = 11),  # legend on the side
    #   labels_gp = gpar(fontsize = 11)),  # legend on the side
    Asthma = list(
      nrow = 2, # number of rows across which the legend will be arranged
      title = 'Asthma',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 11),  # legend on the side
      labels_gp = gpar(fontsize = 11)),  # legend on the side
    Stimuli.Time = list(
      nrow = 3,
      title = 'Stimuli.Time',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 11), # legend on the side
      labels_gp = gpar(fontsize = 11)))) # legend on the side


rowAnn <- HeatmapAnnotation(
  df = row_ann,
  which = 'row', # 'col' (samples) or 'row' (gene) annotation?
  na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
  col = row_colours,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  simple_anno_size = unit(3, "mm"),
  gap = unit(1, 'mm'),
  annotation_legend_param = list(
    Expression = list(
      nrow = 3, 
      title = 'Expression',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 11),  # legend on the side
      labels_gp = gpar(fontsize = 11))  # legend on the side
  ))



hmap <- Heatmap(data,
                # split the genes / rows according to the PAM clusters
                # column_split = factor(metadata$expCond.stimuli.time.HLA.donor,
                #                       levels = expCondReorderLevels),
                # column_split = factor(metadata$expCond.asthma,
                #                       levels = c("Asthmatic", "Non-asthmatic")),
                column_split = factor(metadata$expCond.stimuli.time.asthma,
                                      levels = apply(expand.grid(c("unstim_0h", "Ig_4h", "Ig_18h"),
                                                                 sort(unique(seuratObjFinal@meta.data$expCond.asthma))),
                                                     1, function(x) paste(x[1], x[2], sep = "_"))),
                # column_split = factor(metadata$expCond.stimuli.time.HLA,
                #                       levels = c("unstim", "Ig_4h_None", "Ig_18h_None",
                #                                  "Ig_4h_HLA-DQA2|DQB2", "Ig_18h_HLA-DQA2|DQB2")),
                row_split = factor(c(rep('Down-regulated',20), rep('Up-regulated',20),
                                     rep('Genes of interest',3)),
                                   levels=c('Up-regulated','Down-regulated','Genes of interest')),
                
                # cluster_row_slices = FALSE,
                name = 'Scaled\ngene\nexpression',
                
                col = colorRamp2(myBreaks, myCol),
                use_raster = FALSE,
                # parameters for the colour-bar that represents gradient of expression
                heatmap_legend_param = list(
                  at = c(-2.5, -2, -1, 0, 1, 2, 2.5),
                  color_bar = 'continuous',
                  legend_direction = 'vertical',
                  legend_width = unit(6, 'cm'),
                  legend_height = unit(4, 'cm'),
                  title_position = 'topcenter',
                  title_gp=gpar(fontsize = 11),
                  labels_gp=gpar(fontsize = 11)),
                
                # row (gene) parameters
                cluster_rows = FALSE,
                show_row_dend = FALSE,
                row_title = 'Differentially expressed genes',
                row_title_side = 'left',
                row_title_gp = gpar(fontsize = 11),
                row_title_rot = 90,
                show_row_names = TRUE,
                row_names_gp = gpar(fontsize = 9), # y-axis ticks
                row_names_side = 'left',
                # row_dend_width = unit(25,'mm'),
                
                # column (sample) parameters
                cluster_columns = FALSE,
                show_column_dend = FALSE,
                # column_title = '',
                column_title_side = 'top',
                column_title_gp = gpar(fontsize = 0), # column ticks
                column_title_rot = 90,
                show_column_names = FALSE,
                column_names_gp = gpar(fontsize = 7),
                column_names_max_height = unit(10, 'cm'),
                # column_dend_height = unit(25,'mm')
                # specify top and bottom annotations
                top_annotation = colAnn,
                left_annotation = rowAnn)


pdf(paste0(excel.dir, de.dir, cellType, "_complexheatmap.pdf"), height=6.67, width=7.28)
draw(hmap, 
     # + genelabels,
     heatmap_legend_side = 'right',
     annotation_legend_side = 'top')
dev.off()

# # tetst if there are more B cells expressing HLA-DQ2 genes in asthma 
# contingency_table <- table(metadata$HLA, metadata$expCond.asthma)
# 
# #                 Asthmatic Non-asthmatic
# # HLA-DQA2|DQB2     11184          4376
# # None               7048          9677
# 
# 
# 
# # Perform Fisher's exact test
# fisher_test_result <- fisher.test(contingency_table)
# # Fisher's Exact Test for Count Data
# # 
# # data:  contingency_table
# # p-value < 2.2e-16
# # alternative hypothesis: true odds ratio is not equal to 1
# # 95 percent confidence interval:
# #  3.348830 3.677406
# # sample estimates:
# # odds ratio 
# #   3.508898 
# 
# # strong association between the gene being expressed and asthmatic
# # B cells from asthmatic donors are 3.51 times more likely to express HLA-DQA2/DQB2
# 
# 
# metadata_unstim<- metadata %>%
#   filter(expCond.stimuli.time %in% "unstim_0h")
# contingency_table <- table(metadata_unstim$HLA, metadata_unstim$expCond.asthma)
# res <- fisher.test(contingency_table) #p-value < 2.2e-16, odds ratio 3.810912
# 
# metadata_ig_4h <- metadata %>%
#   filter(expCond.stimuli.time %in% "Ig_4h")
# contingency_table <- table(metadata_ig_4h$HLA, metadata_ig_4h$expCond.asthma)
# res <- fisher.test(contingency_table) #p-value < 2.2e-16, odds ratio 2.062551 






##-- cellType B | Under Ig_18h, what are the de genes for asthmatic vs non-asthmatic----
de.dir <- "asthma_ig_18h_MAST/"
cellType <- "B"

rds.sub <- readRDS(paste0(excel.dir, "B_ig_unstim_hla.rds"))

file.all <- list.files(paste0(excel.dir, de.dir))
file <- file.all[grepl("*_adjSig_1SelClusters.xlsx", file.all, ignore.case = TRUE)]

# expCondCompDeMarkers_Asthmatic-Non-asthmatic_adjSig_1SelClusters.xlsx

# excel_data <- lapply(paste0(excel.dir, de.dir, file), read_excel)
# combined_data <- do.call(rbind, excel_data)

excel_data <- read_excel(paste0(excel.dir, de.dir, file))


top_down <- excel_data %>%
  filter(avg_log2FC < 0 & ...1 %in% pcgs) %>%
  arrange(-abs(avg_log2FC)) %>%  
  head(20) %>%
  rename(Gene = ...1) %>%
  mutate(geneType = "B") %>%
  select(c(Gene, geneType))


top_up <- excel_data %>%
  filter(avg_log2FC > 0 & ...1 %in% pcgs) %>%
  arrange(-abs(avg_log2FC)) %>%  
  head(20) %>%
  rename(Gene = ...1) %>%
  mutate(geneType = "B") %>%
  select(c(Gene, geneType))

combined_data <- bind_rows(top_down, top_up)


geneNames <- c(combined_data$Gene, "HLA-DOA", "HLA-DOB", "IGHD")
features <- geneNames
resDir <- paste0(excel.dir, de.dir)
seuratObjFinal <- rds.sub
# expCondReorderLevels <- c("unstim_602", "Ig_4h_None_602", "Ig_18h_None_602", 
#                           "Ig_4h_DQA2|DQB2_602", "Ig_18h_DQA2|DQB2_602",
#                           "unstim_681", "Ig_4h_None_681","Ig_18h_None_681",
#                           "Ig_4h_DQA2|DQB2_681","Ig_18h_DQA2|DQB2_681")


expCondReorderLevels <- apply(expand.grid(c("unstim_0h", "Ig_4h", "Ig_18h"), 
                                          sort(unique(seuratObjFinal@meta.data$expCond.donor))), 1, function(x) paste(x[2], x[1], sep = "_"))

slot <- "scale.data"

# clusterLevels <- levels(Seurat::Idents(seuratObjFinal)) # cell types
# seuratObjFinal <- subset(seuratObjFinal, idents = "B cells")
Seurat::DefaultAssay(seuratObjFinal) <- "RNA"

# by default, both do.scale/center are on, 
# will scale the expression level for each gene feature by dividing the centered feature expression levels by their standard deviations
seuratObjFinal <- Seurat::ScaleData(object = seuratObjFinal, do.scale = T, do.center = T, features = features) 

# 1. input for heatmap
data <- as.data.frame(as.matrix(t(SeuratObject::GetAssayData(object = seuratObjFinal, 
                                                             layer = slot)[features, colnames(x = seuratObjFinal), drop = FALSE]))) 


seuratObjFinal@meta.data$expCond.stimuli.time.asthma <- paste(rds.sub@meta.data$expCond.stimuli.time, 
                                                              rds.sub@meta.data$expCond.asthma, sep = "_")

metadata <- seuratObjFinal@meta.data[rownames(data), ]

# focus on genes of interest
geneNames = c(combined_data$Gene, "HLA-DOA", "HLA-DOB", "IGHD")
data <- t(data[,geneNames])

data[data < -2.5] <- -2.5
data[data > 2.5] <- 2.5


# 2. set the color scheme
# myCol <- colorRampPalette(c("#31d3ff", "black", "#fff200"))(100)
# myCol <- colorRampPalette(c( "#003049", "#669bbc", "#fdf0d5", "#c1121f", "#780000"))(100)
myCol <- colorRampPalette(c("#3B4CC0", "#FFFFFF", "#B40426"))(100)

# myCol <- colorRampPalette(brewer.pal(5, "RdYlBu"))(100)

myBreaks <- seq(-2.5, 2.5, length.out = 100)

col_ann <- data.frame(
  # Donor = factor(metadata$expCond.donor,levels=sort(unique(metadata$expCond.donor))),
  Asthma = factor(metadata$expCond.asthma,levels=sort(unique(metadata$expCond.asthma))),
  Stimuli.Time = factor(metadata$expCond.stimuli.time,
                        levels=c('unstim_0h',
                                 'Ig_4h',
                                 'Ig_18h')),
  stringsAsFactors = FALSE) #do not have to set factor 

row_ann <- data.frame(
  # Genes = rownames(data),
  Expression = c(rep('Down-regulated',20), rep('Up-regulated',20),
                 rep('Genes of interest',3)),
  stringsAsFactors = FALSE)

donors <- sort(unique(metadata$expCond.donor))

# Generate a color gradient
gradient_colors <- colorRampPalette(c("#f2a65a", "#313715"))(length(donors))

# Create the named vector
donor_colors <- setNames(gradient_colors, donors)

col_colours <- list(
  # Donor = donor_colors,
  Asthma = c("Asthmatic" = "#471323",
             "Non-asthmatic" = "#CEB992"),
  Stimuli.Time = c('unstim_0h' = '#eff9f0', 
                   'Ig_4h' = '#FFBCB5', 
                   'Ig_18h' = '#03978E'))


row_colours <- list(
  Expression = c('Down-regulated' = '#6C8FC8',
                 'Up-regulated' = '#EE4687',
                 'Genes of interest' = '#f1faee'))


colAnn <- HeatmapAnnotation(
  df = col_ann,
  which = 'col', # 'col' (samples) or 'row' (gene) annotation?
  na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
  col = col_colours,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  simple_anno_size = unit(3, "mm"),
  gap = unit(1, 'mm'),
  annotation_legend_param = list(
    # Donor = list(
    #   nrow = 4, # number of rows across which the legend will be arranged
    #   title = 'Donor',
    #   title_position = 'topcenter',
    #   legend_direction = 'vertical',
    #   title_gp = gpar(fontsize = 11),  # legend on the side
    #   labels_gp = gpar(fontsize = 11)),  # legend on the side
    Asthma = list(
      nrow = 2, # number of rows across which the legend will be arranged
      title = 'Asthma',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 11),  # legend on the side
      labels_gp = gpar(fontsize = 11)),  # legend on the side
    Stimuli.Time = list(
      nrow = 3,
      title = 'Stimuli.Time',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 11), # legend on the side
      labels_gp = gpar(fontsize = 11)))) # legend on the side


rowAnn <- HeatmapAnnotation(
  df = row_ann,
  which = 'row', # 'col' (samples) or 'row' (gene) annotation?
  na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
  col = row_colours,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  simple_anno_size = unit(3, "mm"),
  gap = unit(1, 'mm'),
  annotation_legend_param = list(
    Expression = list(
      nrow = 3, 
      title = 'Expression',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 11),  # legend on the side
      labels_gp = gpar(fontsize = 11))  # legend on the side
  ))



hmap <- Heatmap(data,
                # split the genes / rows according to the PAM clusters
                # column_split = factor(metadata$expCond.stimuli.time.HLA.donor,
                #                       levels = expCondReorderLevels),
                # column_split = factor(metadata$expCond.asthma,
                #                       levels = c("Asthmatic", "Non-asthmatic")),
                column_split = factor(metadata$expCond.stimuli.time.asthma,
                                      levels = apply(expand.grid(c("unstim_0h", "Ig_4h", "Ig_18h"),
                                                                 sort(unique(seuratObjFinal@meta.data$expCond.asthma))),
                                                     1, function(x) paste(x[1], x[2], sep = "_"))),
                # column_split = factor(metadata$expCond.stimuli.time.HLA,
                #                       levels = c("unstim", "Ig_4h_None", "Ig_18h_None",
                #                                  "Ig_4h_HLA-DQA2|DQB2", "Ig_18h_HLA-DQA2|DQB2")),
                row_split = factor(c(rep('Down-regulated',20), rep('Up-regulated',20),
                                     rep('Genes of interest',3)),
                                   levels=c('Up-regulated','Down-regulated','Genes of interest')),
                
                # cluster_row_slices = FALSE,
                name = 'Scaled\ngene\nexpression',
                
                col = colorRamp2(myBreaks, myCol),
                use_raster = FALSE,
                # parameters for the colour-bar that represents gradient of expression
                heatmap_legend_param = list(
                  at = c(-2.5, -2, -1, 0, 1, 2, 2.5),
                  color_bar = 'continuous',
                  legend_direction = 'vertical',
                  legend_width = unit(6, 'cm'),
                  legend_height = unit(4, 'cm'),
                  title_position = 'topcenter',
                  title_gp=gpar(fontsize = 11),
                  labels_gp=gpar(fontsize = 11)),
                
                # row (gene) parameters
                cluster_rows = FALSE,
                show_row_dend = FALSE,
                row_title = 'Differentially expressed genes',
                row_title_side = 'left',
                row_title_gp = gpar(fontsize = 11),
                row_title_rot = 90,
                show_row_names = TRUE,
                row_names_gp = gpar(fontsize = 9), # y-axis ticks
                row_names_side = 'left',
                # row_dend_width = unit(25,'mm'),
                
                # column (sample) parameters
                cluster_columns = FALSE,
                show_column_dend = FALSE,
                # column_title = '',
                column_title_side = 'top',
                column_title_gp = gpar(fontsize = 0), # column ticks
                column_title_rot = 90,
                show_column_names = FALSE,
                column_names_gp = gpar(fontsize = 7),
                column_names_max_height = unit(10, 'cm'),
                # column_dend_height = unit(25,'mm')
                # specify top and bottom annotations
                top_annotation = colAnn,
                left_annotation = rowAnn)


pdf(paste0(excel.dir, de.dir, cellType, "_complexheatmap.pdf"), height=6.67, width=7.28)
draw(hmap, 
     # + genelabels,
     heatmap_legend_side = 'right',
     annotation_legend_side = 'top')
dev.off()

# # tetst if there are more B cells expressing HLA-DQ2 genes in asthma 
# contingency_table <- table(metadata$HLA, metadata$expCond.asthma)
# 
# #                 Asthmatic Non-asthmatic
# # HLA-DQA2|DQB2     11184          4376
# # None               7048          9677
# 
# 
# 
# # Perform Fisher's exact test
# fisher_test_result <- fisher.test(contingency_table)
# # Fisher's Exact Test for Count Data
# # 
# # data:  contingency_table
# # p-value < 2.2e-16
# # alternative hypothesis: true odds ratio is not equal to 1
# # 95 percent confidence interval:
# #  3.348830 3.677406
# # sample estimates:
# # odds ratio 
# #   3.508898 
# 
# # strong association between the gene being expressed and asthmatic
# # B cells from asthmatic donors are 3.51 times more likely to express HLA-DQA2/DQB2
# 
# 
# metadata_unstim<- metadata %>%
#   filter(expCond.stimuli.time %in% "unstim_0h")
# contingency_table <- table(metadata_unstim$HLA, metadata_unstim$expCond.asthma)
# res <- fisher.test(contingency_table) #p-value < 2.2e-16, odds ratio 3.810912
# 
# metadata_ig_4h <- metadata %>%
#   filter(expCond.stimuli.time %in% "Ig_4h")
# contingency_table <- table(metadata_ig_4h$HLA, metadata_ig_4h$expCond.asthma)
# res <- fisher.test(contingency_table) #p-value < 2.2e-16, odds ratio 2.062551 








############################### % HLA-DQA2/DQB2 per donor and treatment ###############################
# figures_dir <- "results_wOrgClusterAnnotation_DEGs/MAST_HLA_ig_unstim/results_wOrgClusterAnnotation/heatmap_expCond.stimuli.time.HLA.donor/"
deResDir <- "HLA-DQA2|DQB2_ig_unstim_MAST/"

rds.sub <- readRDS(paste0(excel.dir, "B_ig_unstim_hla.rds"))


rds.sub@meta.data$HLA_tmp <- rds.sub@meta.data$HLA

rds.sub@meta.data$HLA <- ifelse(
  rds.sub@meta.data$HLA_tmp != "None", 
  'HLA-DQA2|DQB2', 
  "None"
)
# 616, 673, 662 have 18, 21, 23 B cells in total for each, exclude 

tmp <- rds.sub@meta.data %>% 
  count(expCond.donor) %>%
  group_by(expCond.donor) %>%
  mutate(total_n = sum(n)) 

tmp %>% 
  filter(!expCond.donor %in% c(616, 673, 662)) %>% 
  pull(expCond.donor)

p1 <- rds.sub@meta.data %>%
  count(expCond.donor, expCond.asthma, expCond.stimuli.time, HLA)  %>%
  filter(!expCond.donor %in% c(616, 673, 662)) %>%
  ggplot(aes(x = factor(expCond.stimuli.time,
                        levels = c('unstim_0h', 'Ig_4h', 'Ig_18h'),
                        labels = c('Unstim_0h', 'Ig_4h', 'Ig_18h')), y = n, fill = HLA)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(
    factor(expCond.asthma, levels = c("Asthmatic", "Non-asthmatic")) +
      factor(expCond.donor, levels = unique(expCond.donor[order(expCond.asthma)])) ~., 
    space ="free") +  # Reorder donors by asthma
  scale_fill_manual(values = c("HLA-DQA2|DQB2" = "#B40426",
                               "None" = "#3B4CC0")) +
  theme_bw() +
  labs(
    x = "",
    y = "Number of B cells", 
    fill = "HLA-DQA2/DQB2"
  ) +
  theme(
    strip.background = element_rect(fill ="white"),
    axis.text.x = element_text(angle = 0, hjust = 0.5), # Adjust x-axis labels
  )

p2 <- rds.sub@meta.data %>%
  count(expCond.asthma, expCond.stimuli.time, HLA)  %>%
  # filter(!expCond.donor %in% c(616, 673, 662)) %>%
  ggplot(aes(x = factor(expCond.stimuli.time,
                        levels = c('unstim_0h', 'Ig_4h', 'Ig_18h'),
                        labels = c('Unstim_0h', 'Ig_4h', 'Ig_18h')), y = n, fill = HLA)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~ expCond.asthma,
             space = "free") +  # Reorder donors by asthma
  scale_fill_manual(values = c("HLA-DQA2|DQB2" = "#B40426",
                               "None" = "#3B4CC0")) +
  theme_bw() +
  labs(
    x = "",
    y = "Number of B cells", 
    fill = "HLA-DQA2/DQB2"
  ) +
  theme(
    strip.background = element_rect(fill ="white"),
    axis.text.x = element_text(angle = 0, hjust = 0.5), # Adjust x-axis labels
    
  )


# unstim ig B cell counts
p3 <- rds.sub@meta.data %>%
  count(expCond.donor, expCond.asthma, HLA)  %>%
  ggplot(aes(x = factor(expCond.donor), y = n, fill = HLA)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("HLA-DQA2|DQB2" = "#B40426",
                               "None" = "#3B4CC0")) +
  theme_bw() +
  labs(
    x = "",
    y = "Number of B cells", 
    fill = "Unstim & Ig"
  ) +
  theme(
    strip.background = element_rect(fill ="white"),
    axis.text.x = element_text(angle = 90, hjust = 0.5), # Adjust x-axis labels
  ) +
  facet_grid(~expCond.asthma, scales = "free_x", space = "free_x")


# unstim ig B cell freq
p4 <- rds.sub@meta.data %>%
  count(expCond.donor, expCond.asthma, HLA)  %>%
  ggplot(aes(x = factor(expCond.donor), y = n, fill = HLA)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("HLA-DQA2|DQB2" = "#B40426",
                               "None" = "#3B4CC0")) +
  theme_bw() +
  labs(
    x = "",
    y = "Proportion of B cells", 
    fill = "Unstim & Ig"
  ) +
  theme(
    strip.background = element_rect(fill ="white"),
    axis.text.x = element_text(angle = 90, hjust = 0.5), # Adjust x-axis labels
  ) +
  facet_grid(~expCond.asthma, scales = "free_x", space = "free_x")


# unstim B cell counts
p5 <- rds.sub@meta.data %>%
  filter(expCond.stimuli.time %in% "unstim_0h") %>%
  count(expCond.donor, expCond.asthma, HLA)  %>%
  ggplot(aes(x = factor(expCond.donor), y = n, fill = HLA)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("HLA-DQA2|DQB2" = "#B40426",
                               "None" = "#3B4CC0")) +
  theme_bw() +
  labs(
    x = "",
    y = "Number of B cells", 
    fill = "Unstim"
  ) +
  theme(
    strip.background = element_rect(fill ="white"),
    axis.text.x = element_text(angle = 90, hjust = 0.5), # Adjust x-axis labels
  ) +
  facet_grid(~expCond.asthma, scales = "free_x", space = "free_x")


# unstim B cell freq
p6 <- rds.sub@meta.data %>%
  filter(expCond.stimuli.time %in% "unstim_0h") %>%
  count(expCond.donor, expCond.asthma, HLA)  %>%
  ggplot(aes(x = factor(expCond.donor), y = n, fill = HLA)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("HLA-DQA2|DQB2" = "#B40426",
                               "None" = "#3B4CC0")) +
  theme_bw() +
  labs(
    x = "",
    y = "Proportion of B cells", 
    fill = "Unstim"
  ) +
  theme(
    strip.background = element_rect(fill ="white"),
    axis.text.x = element_text(angle = 90, hjust = 0.5), # Adjust x-axis labels
  ) +
  facet_grid(~expCond.asthma, scales = "free_x", space = "free_x")

# Ig 4h B cell counts
p7 <- rds.sub@meta.data %>%
  filter(expCond.stimuli.time %in% "Ig_4h") %>%
  count(expCond.donor, expCond.asthma, HLA)  %>%
  ggplot(aes(x = factor(expCond.donor), y = n, fill = HLA)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("HLA-DQA2|DQB2" = "#B40426",
                               "None" = "#3B4CC0")) +
  theme_bw() +
  labs(
    x = "",
    y = "Number of B cells", 
    fill = "Ig_4h"
  ) +
  theme(
    strip.background = element_rect(fill ="white"),
    axis.text.x = element_text(angle = 90, hjust = 0.5), # Adjust x-axis labels
  ) +
  facet_grid(~expCond.asthma, scales = "free_x", space = "free_x")


# Ig 4h B cell freq
p8 <- rds.sub@meta.data %>%
  filter(expCond.stimuli.time %in% "Ig_4h") %>%
  count(expCond.donor, expCond.asthma, HLA)  %>%
  ggplot(aes(x = factor(expCond.donor), y = n, fill = HLA)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("HLA-DQA2|DQB2" = "#B40426",
                               "None" = "#3B4CC0")) +
  theme_bw() +
  labs(
    x = "",
    y = "Proportion of B cells", 
    fill = "Ig_4h"
  ) +
  theme(
    strip.background = element_rect(fill ="white"),
    axis.text.x = element_text(angle = 90, hjust = 0.5), # Adjust x-axis labels
  ) +
  facet_grid(~expCond.asthma, scales = "free_x", space = "free_x")



p3 / p4
p5 / p6
p7 / p8

pdf(paste0(excel.dir, "HLA-DQ2_B_unstim.ig_asthma_donor_barchart.pdf"),
    width = 4.18,
    height = 7.98)
print(p1)
dev.off()

pdf(paste0(excel.dir, "HLA-DQ2_B_unstim.ig_asthma_barchart.pdf"),
    width = 5.63,
    height = 2.41)
print(p2)
dev.off()

pdf(paste0(excel.dir, "HLA-DQ2_B_unstim.ig_by_asthma_donor_barcharts.pdf"),
    width = 8.31,
    height = 6.41)
(p3 + p4) / (p5 + p6) / (p7 + p8)
dev.off()















HLA_DQ2_percent_increase <- rds.sub@meta.data %>% 
  count(expCond.donor, expCond.asthma, expCond.stimuli.time, HLA)%>%
  group_by(expCond.donor, expCond.stimuli.time) %>%
  mutate(total_n = sum(n)) %>%
  filter(all(total_n > 20)) %>% # same with excluding donors, one donor & time combination less than 20, all removed
  mutate(percentage = n/total_n*100) %>%
  ungroup() %>%
  filter(HLA != "None") %>%
  group_by(expCond.donor) %>%
  mutate(baseline = percentage[expCond.stimuli.time == "unstim_0h"],
         percent_increase = (percentage - baseline) / baseline * 100) %>%
  ungroup()

rds.sub@meta.data %>% 
  select(expCond.donor, expCond.asthma) %>%
  unique()


HLA_DQ2_expressed <- rds.sub@meta.data %>% 
  count(expCond.donor, expCond.asthma, HLA) %>%
  filter(!expCond.donor %in% c(616, 673, 662)) %>%
  filter(HLA != "None")

HLA_DQ2_expressed_time <- rds.sub@meta.data %>% 
  count(expCond.donor, expCond.asthma, expCond.stimuli.time, HLA) %>%
  filter(!expCond.donor %in% c(616, 673, 662)) %>%
  filter(HLA != "None")



HLA_DQ2_expressed_percent <- rds.sub@meta.data %>%
  count(expCond.donor, expCond.asthma, HLA) %>%
  filter(!expCond.donor %in% c(616, 673, 662)) %>%
  group_by(expCond.donor) %>%
  mutate(total_n = sum(n)) %>%
  # filter(all(n > 10)) %>%
  mutate(percentage = n/total_n*100) %>%
  ungroup() %>%
  filter(HLA != "None")


HLA_DQ2_expressed_time_percent <- rds.sub@meta.data %>%
  count(expCond.donor, expCond.asthma, expCond.stimuli.time, HLA) %>%
  filter(!expCond.donor %in% c(616, 673, 662)) %>%
  group_by(expCond.donor, expCond.stimuli.time) %>%
  mutate(total_n = sum(n)) %>%
  # filter(all(n > 10)) %>%
  mutate(percentage = n/total_n*100) %>%
  ungroup() %>%
  filter(!HLA %in% "None") %>%
  filter(expCond.stimuli.time != "Ig_18h") 

  
p1 <- ggplot(HLA_DQ2_expressed, aes(x = expCond.asthma, y = n, 
             fill = as.character(expCond.asthma))) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  ylim(0, 12000) +
  geom_jitter(width = 0.2, size = 1, shape = 21, alpha = 0.8) +  # Add jitter
  scale_fill_manual(values = c("Asthmatic" = "#471323",
                                "Non-asthmatic" = "#CEB992")) +
  labs(
    x = "", 
    y = "Number of \nHLA-DQA2|DQB2-expressing\nB cells ", 
    fill = "Asthma diagnosis"
  ) +
  theme_bw() +
  stat_compare_means(aes(group = expCond.asthma),
                     method = "t.test", # Use "t.test" for pairwise comparison or "anova" for multiple groups
                     label = "p.signif")

t.test(n ~ expCond.asthma, 
       data = HLA_DQ2_expressed,
       alternative = 'greater')$p.value



t.test(n ~ expCond.asthma, 
       data = HLA_DQA2_expressed,
       alternative = 'greater')$p.value


p2 <- ggplot(HLA_DQ2_expressed_percent, aes(x = expCond.asthma, y = percentage,
                              fill = as.character(expCond.asthma))) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  ylim(0, 100) +
  geom_jitter(width = 0.2, size = 1, shape = 21, alpha = 0.8) +  # Add jitter
  scale_fill_manual(values = c("Asthmatic" = "#471323",
                               "Non-asthmatic" = "#CEB992")) +
  labs(
    x = "",
    y = "% of \nHLA-DQA2|DQB2-expressing\nB cells ",
    fill = "Asthma diagnosis"
  ) +
  theme_bw() +
  stat_compare_means(aes(group = expCond.asthma),
                     method = "t.test", # Use "t.test" for pairwise comparison or "anova" for multiple groups
                     label = "p.signif")

t.test(percentage ~ expCond.asthma, 
       data = HLA_DQ2_expressed_percent,
       alternative = 'greater')$p.value







p3 <- ggplot(HLA_DQ2_expressed_time_percent, aes(x = expCond.asthma, y = percentage,
                                                 fill = as.character(expCond.asthma))) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  ylim(0, 100) +
  geom_jitter(width = 0.2, size = 1, shape = 21, alpha = 0.8) +  # Add jitter
  scale_fill_manual(values = c("Asthmatic" = "#471323",
                               "Non-asthmatic" = "#CEB992")) +
  labs(
    x = "",
    y = "% of \nHLA-DQA2|DQB2-expressing\nB cells ",
    fill = "Asthma diagnosis"
  ) +
  theme_bw() +
  facet_grid( ~ factor(expCond.stimuli.time,
                       levels = c('unstim_0h', 'Ig_4h'),
                       labels = c('Unstim_0h', 'Ig_4h')), 
              scale = "free", space = "free") +
  theme(
    strip.background = element_rect(fill = "white", color = "black"))+
  stat_compare_means(aes(group = expCond.asthma),
                     method = "t.test", # Use "t.test" for pairwise comparison or "anova" for multiple groups
                     label = "p.signif")



HLA_DQ2_expressed_time_percent %>%
  group_by(expCond.stimuli.time) %>%  #
  summarise(
    p_value = t.test(percentage ~ expCond.asthma, 
                     data = pick(everything()), 
                     alternative = 'greater')$p.value
  )


p4 <- ggplot(HLA_DQ2_percent_increase, aes(x = factor(expCond.stimuli.time,
                                        levels = c("unstim_0h", "Ig_4h", "Ig_18h")), 
             y = percent_increase, 
             color = as.character(expCond.donor),
             group = as.character(expCond.donor))) +
  geom_line(linewidth=0.3) +
  geom_point(size=2) +
  geom_label_repel(data = HLA_DQ2_percent_increase %>% filter(expCond.stimuli.time == "Ig_4h"),
                   aes(label = as.character(expCond.donor),
                       fill = expCond.asthma),
                   size = 3,
                   color = "black",
                   # label.padding = 0.5,  # Padding inside label box
                   label.r = 0.15,  # Round label corners
                   label.size = 0.25,  # Border thickness
                   alpha = 0.9,
                   box.padding = 0.3,  # Adjust box padding
                   point.padding = 0.5,
                   max.overlaps = 20) +
  theme_bw() +
  # scale_color_manual(values = c("#f2a65a","#394787")) +
  scale_color_manual(values = c('#ce0000', '#008080', '#ccb100', '#4363d8', 
                                '#bcf60c', '#f58231', '#46f0f0', '#f032e6', 
                                '#911eb4')) +
  scale_fill_manual(values = c("Asthmatic" = "#be929a", "Non-asthmatic" = "#f6f2e2")) +
  scale_x_discrete(labels = c("unstim_0h" = "Unstim_0h")) +
  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, by = 50)) +
  labs(x = "", 
       y = "% increase of \nHLA-DQA2|DQB2-expressing\nB cells", 
       color = "Donors (n=9)",
       fill = "Asthma diagnosis") +
  guides(color = guide_legend(ncol = 2), fill = guide_legend(ncol = 1))
# theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  # facet_wrap(~ expCond.donor)


p5 <- ggplot(HLA_DQ2_percent_increase, aes(x = factor(expCond.stimuli.time,
                                                      levels = c("unstim_0h", 
                                                                 "Ig_4h", 
                                                                 "Ig_18h")), 
                                           y = percent_increase, 
                                           color = as.character(expCond.asthma),
                                           group = as.character(expCond.donor))) +
  geom_line(linewidth=0.3) +
  geom_point(size=2) +
  theme_bw() +
  scale_color_manual(values = c("Asthmatic" = "#471323", "Non-asthmatic" = "#CEB992")) +
  scale_x_discrete(labels = c("unstim_0h" = "Unstim_0h")) +
  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, by = 50)) +
  labs(x = "", 
       y = "% increase of \nHLA-DQA2|DQB2-\nexpressing\nB cells", 
       color = "Donors (n=9)",
       fill = "Asthma diagnosis")



pdf(paste0(excel.dir, "HLA-DQ2_B_unstim.ig_asthma.pdf"),
    width = 6.80,
    height = 2.77)

(p1 | p2)  +
  plot_layout(guides = "collect")

dev.off()




pdf(paste0(excel.dir, "HLA-DQ2_B_by_unstim.ig_asthma.pdf"),
    width = 5.64,
    height = 2.77)

p3

dev.off()


pdf(paste0(excel.dir, "HLA-DQ2_B_unstim.ig_asthma_percent_increase_donors.pdf"),
    width = 4.28,
    height = 2.77)

p4

dev.off()




pdf(paste0(excel.dir, "HLA-DQ2_B_unstim.ig_asthma_percent_increase.pdf"),
    width = 4.32,
    height = 2.01)

p5

dev.off()












############################### avg exp HLA-DQA2/DQB2 per donor and treatment ####
deResDir <- "HLA-DQA2|DQB2_ig_unstim_MAST/"

rds.sub <- readRDS(paste0(excel.dir, "B_ig_unstim_hla.rds"))

DefaultAssay(rds.sub) <- "RNA" # do not use the integrated layer 

rds.sub@meta.data$expCond.donor.stimuli.time.asthma <- paste0(
  rds.sub@meta.data$expCond.donor.stimuli.time, 
  "|", 
  rds.sub@meta.data$expCond.asthma
)


sel_genes <- c("HLA-DQA2", "HLA-DQB2")
# existing_genes <- intersect(unique(sel_genes), rownames(GetAssayData(object = rds.sub, layer = "data")))

Idents(rds.sub) <- "expCond.donor.stimuli.time"

dotplot <- DotPlot(rds.sub,
                   features = c(sel_genes),
                   cols = c('#D3D3D3', '#CC0000'),
                   scale = T, scale.by = 'size',
                   dot.min = 0) + RotatedAxis()


sel_exp <- as.data.frame(t(GetAssayData(object = rds.sub, layer = "data")[c(sel_genes), ]))




# Calculate the average expression per cluster
sel_avg_exp <- sel_exp %>%
  rownames_to_column(var = "cell") %>%
  inner_join(rds.sub@meta.data %>% rownames_to_column(var = "cell"), by = "cell") %>%
  group_by(expCond.donor.stimuli.time.asthma) %>%
  summarise(across(2:(ncol(sel_exp) + 1), \(x) mean(x, na.rm = TRUE)))

# Calculate the percentage of cells expressing each gene per cluster
sel_pct_exp <- sel_exp %>%
  mutate(both = rowSums(across(all_of(sel_genes)))) %>% 
  rownames_to_column(var = "cell") %>%
  inner_join(rds.sub@meta.data %>% rownames_to_column(var = "cell"), by = "cell") %>%
  group_by(expCond.donor.stimuli.time.asthma) %>%
  # summarise(across(2:(ncol(sel_exp)+1), ~mean(. > 0) * 100))
  summarise(
    across(2:(ncol(sel_exp) + 1 + 1), list(
      percent = ~mean(. > 0) * 100,
      count = ~sum(. > 0)
    ), .names = "{.col}_{.fn}")
  ) # HLA-DQA2, HLA_DQB2, both -> ncol(sel_exp) + 1

# Combine average expression and percentage expression data
sel_avg_pct_exp <- sel_avg_exp %>%
  pivot_longer(-expCond.donor.stimuli.time.asthma, names_to = "gene", 
               values_to = "expression") %>%
  # inner_join(sel_pct_exp %>% pivot_longer(-expCond.donor.stimuli.time.asthma, 
  #                                        names_to = "gene", values_to = "pct.exp"), 
  #            by = c("expCond.donor.stimuli.time.asthma", "gene"))
  full_join(
    sel_pct_exp %>%
      pivot_longer(-expCond.donor.stimuli.time.asthma, 
                   names_to = c("gene", ".value"), 
                   names_sep = "_"),  # Separates 'percent' and 'count' into different columns
    by = c("expCond.donor.stimuli.time.asthma", "gene")
  )



sel_avg_pct_exp <- sel_avg_pct_exp %>%
  separate(
    expCond.donor.stimuli.time.asthma,
    into = c("expCond.donor", "expCond.stimuli.time.asthma"),
    sep = "_",
    extra = "merge" 
  ) %>%
  separate(
    expCond.stimuli.time.asthma,  
    into = c("expCond.stimuli.time", "expCond.asthma"), 
    sep = "\\|"
  )

# sel_avg_pct_exp <- sel_avg_pct_exp %>%
#   mutate(cluster_annotation_plot = case_when(
#     # cluster_annotation == "Alveolar Mph MT-positive" ~ "Alveolar Mph\nMT-positive",
#     # cluster_annotation == "EC general capillary" ~ "EC general\ncapillary",
#     # cluster_annotation == "Monocyte-derived Mph" ~ "Monocyte-\nderived Mph",
#     TRUE ~ cluster_annotation  # Default case for other annotations
#   ))

# sel_avg_pct_exp <- sel_avg_pct_exp %>%
#   mutate(pct.exp.shape = case_when(
#     pct.exp < 25 ~ "<25%",
#     pct.exp >= 25 & pct.exp < 50 ~ "25-50%",
#     pct.exp >= 50 & pct.exp < 75 ~ "50-75%",
#     pct.exp >= 75 ~ "75-100%"
#   ))


# sel_avg_pct_exp <- sel_avg_pct_exp %>%
#   mutate(gene_cat = case_when(
#     # gene == "HEY2" ~ "Up-regulated",
#     # gene == "HES7" ~ "Up-regulated",
#     # gene == "SP1" ~ "Dowm-regulated",
#     # gene == "E2F4" ~ "Dowm-regulated",
#     TRUE ~ gene  # Default case for other annotations
#   ))


sel_avg_pct_exp <- sel_avg_pct_exp %>%
  mutate(expCond.stimuli.time = str_replace_all(expCond.stimuli.time, 
                                                c("unstim_0h" = "Unstim_0h")))

B_donor_cnt <- sel_avg_pct_exp %>%
  group_by(expCond.donor, gene) %>%
  summarise(count_sum = sum(count, na.rm = TRUE)) # good



t.test(expression ~ expCond.asthma, 
       data = sel_avg_pct_exp %>%
         filter(expCond.stimuli.time %in% "Unstim_0h",
                gene %in% "HLA-DQA2",
                !expCond.donor %in% c(662, 673, 616)),
       alternative = "greater")

t.test(expression ~ expCond.asthma, 
       data = sel_avg_pct_exp %>%
         filter(expCond.stimuli.time %in% "Unstim_0h",
                gene %in% "HLA-DQB2",
                !expCond.donor %in% c(662, 673, 616)),
       alternative = "greater")


t.test(expression ~ expCond.asthma, 
       data = sel_avg_pct_exp %>%
         filter(expCond.stimuli.time %in% "Ig_4h",
                gene %in% "HLA-DQA2",
                !expCond.donor %in% c(662, 673, 616)),
       alternative = "greater")

t.test(expression ~ expCond.asthma, 
       data = sel_avg_pct_exp %>%
         filter(expCond.stimuli.time %in% "Ig_4h",
                gene %in% "HLA-DQB2",
                !expCond.donor %in% c(662, 673, 616)),
       alternative = "greater")



is_grouped_df(sel_avg_pct_exp)



sel_avg_pct_exp_fc <- sel_avg_pct_exp %>%
  filter(expCond.stimuli.time %in% c("Unstim_0h", "Ig_4h"),
         !expCond.donor %in% c(662, 673, 616)) %>%
  pivot_wider(names_from = expCond.stimuli.time, values_from = c(expression, count, percent)) %>%
  mutate(
    expression_fc = `expression_Ig_4h` / `expression_Unstim_0h`,
    percent_fc = `percent_Ig_4h` / `percent_Unstim_0h`
  ) 


sel_avg_pct_exp_fc <- sel_avg_pct_exp_fc %>%
  mutate(expCond.genotype = case_when(
    expCond.donor %in% c(569, 673, 548, 633, 640, 681) ~ "Heterozygous",
    expCond.donor %in% c(616) ~ "Homozygous non-risk/non-risk",
    expCond.donor %in% c(589, 602, 662, 678, 699) ~ "Homozygous risk/risk",
    TRUE ~ NA_character_ 
  ))


sel_avg_pct_exp_fc %>%
  select(expCond.asthma, expCond.donor, expCond.genotype) %>%
  group_by(expCond.asthma) %>%
  distinct()



sel_avg_pct_exp_fc_indiv <- sel_avg_pct_exp_fc %>% filter(!gene %in% "both") 
sel_avg_pct_exp_fc_both <- sel_avg_pct_exp_fc %>% filter(gene %in% "both") 


p1 <-ggplot(sel_avg_pct_exp_fc_indiv, aes(x = gene, y = percent_fc)) +
  # stat_boxplot(geom = "errorbar", position = position_dodge(width = 0.75), width = 0.4) +  #  error bars
  geom_boxplot(outlier.shape = NA, 
               position = position_dodge(width = 0.75),
               alpha = 0.5) +
  geom_jitter(aes(fill = expCond.genotype,
                  shape = expCond.asthma), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              alpha = 1,
              size = 2,
              color = "black") + 
  scale_y_continuous(breaks = seq(floor(min(sel_avg_pct_exp_fc_indiv$percent_fc)), 
                                  ceiling(max(sel_avg_pct_exp_fc_indiv$percent_fc)), 
                                  by = 0.25)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +  
  scale_shape_manual(values = c("Asthmatic" = 24, "Non-asthmatic" = 21)) +  
  scale_fill_manual(values = c("Heterozygous"="#bbd7f1",
                               "Homozygous non-risk/non-risk"= "#90d3d0",
                               "Homozygous risk/risk" = "#8b7ab8",
                               "Asthmatic" = "#471323",
                               "Non-asthmatic" = "#CEB992")) +
  # scale_fill_manual(values = c("Non-asthmatic.Heterozygous"="#bbd7f1",
  #                              "Asthmatic.Heterozygous"= "#bbd7f1",
  #                              "Non-asthmatic.Homozygous risk/risk" = "#8b7ab8",
  #                              "Asthmatic.Homozygous risk/risk" = "#8b7ab8",
  #                              "Asthmatic" = "#471323",
  #                              "Non-asthmatic" = "#CEB992")) +
  labs(x = "Genes", y = "Fold change of \n%HLA-DQA2/DQB2-expressing\n B cells", 
       shape = "expCond.asthma",
       fill = "expCond.asthma")

p2 <- ggplot(sel_avg_pct_exp_fc_indiv, aes(x = gene, y = expression_fc)) +
  # stat_boxplot(geom = "errorbar", position = position_dodge(width = 0.75), width = 0.4) +  #  error bars
  geom_boxplot(outlier.shape = NA, 
               position = position_dodge(width = 0.75),
               alpha = 0.5) +
  geom_jitter(aes(fill = expCond.genotype,
                  shape = expCond.asthma), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              alpha = 1,
              size = 2,
              color = "black") + 
  scale_y_continuous(breaks = seq(floor(min(sel_avg_pct_exp_fc_indiv$expression_fc)), 
                                  ceiling(max(sel_avg_pct_exp_fc_indiv$expression_fc)), 
                                  by = 0.25)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +  
  scale_fill_manual(values = c("Heterozygous"="#bbd7f1",
                               "Homozygous non-risk/non-risk"= "#90d3d0",
                               "Homozygous risk/risk" = "#8b7ab8",
                               "Asthmatic" = "#471323",
                               "Non-asthmatic" = "#CEB992")) +
  # scale_fill_manual(values = c("Non-asthmatic.Heterozygous"="#bbd7f1",
  #                              "Asthmatic.Heterozygous"= "#bbd7f1",
  #                              "Non-asthmatic.Homozygous risk/risk" = "#8b7ab8",
  #                              "Asthmatic.Homozygous risk/risk" = "#8b7ab8",
  #                              "Asthmatic" = "#471323",
  #                              "Non-asthmatic" = "#CEB992")) +
  scale_shape_manual(values = c("Asthmatic" = 24, "Non-asthmatic" = 21)) +  
  labs(x = "Genes", y = "Fold change of expression\n \n", 
       shape = "expCond.asthma",
       fill = "expCond.asthma")


p3 <-ggplot(sel_avg_pct_exp_fc_both, aes(x = gene, y = percent_fc)) +
  # stat_boxplot(geom = "errorbar", position = position_dodge(width = 0.75), width = 0.4) +  #  error bars
  geom_boxplot(outlier.shape = NA, 
               position = position_dodge(width = 0.75),
               alpha = 0.5) +
  geom_jitter(aes(fill = expCond.genotype,
                  shape = expCond.asthma), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              alpha = 1,
              size = 2,
              color = "black") + 
  scale_y_continuous(breaks = seq(floor(min(sel_avg_pct_exp_fc_both$percent_fc)), 
                                  ceiling(max(sel_avg_pct_exp_fc_both$percent_fc)), 
                                  by = 0.25)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +  
  scale_shape_manual(values = c("Asthmatic" = 24, "Non-asthmatic" = 21)) +  
  scale_fill_manual(values = c("Heterozygous"="#bbd7f1",
                               "Homozygous non-risk/non-risk"= "#90d3d0",
                               "Homozygous risk/risk" = "#8b7ab8",
                               "Asthmatic" = "#471323",
                               "Non-asthmatic" = "#CEB992")) +
  scale_x_discrete(labels = c("both" = "HLA-DQA2/DQB2")) + 
  labs(x = "Genes", y = "Fold change of \n%HLA-DQA2/DQB2-expressing\n B cells", 
       shape = "expCond.asthma",
       fill = "expCond.asthma")


pdf(paste0(excel.dir, deResDir, "GOEnrichment/HLA-DQ2_boxplot_indiv.pdf"), height=2.46, width=8.46)

print(p1 + p2)

dev.off()

pdf(paste0(excel.dir, deResDir, "GOEnrichment/HLA-DQ2_boxplot_both.pdf"), height=2.46, width=3.68)

print(p3)

dev.off()



# Scale the expression data if needed
sel_avg_pct_exp <- sel_avg_pct_exp %>%
  group_by(gene) %>%
  mutate(expression = scale(expression)) %>%
  ungroup()

# make sure to use scaled expression
p <- ggplot(sel_avg_pct_exp, aes(x = factor(expCond.stimuli.time, 
                                            levels=c("Unstim_0h", "Ig_4h", "Ig_18h")), 
                                 y = factor(gene, levels=unique(sel_avg_pct_exp$gene)),
                                 size = pct.exp, 
                                 fill = expression)) +
  # scale_shape_manual(values =  c("<25%" = 21, "25-50%" = 22, "50-75%" = 23, "75-100%" = 24)) +
  geom_point(shape=21) +
  # scale_color_manual(values =  c("<0.01" = "#0000ff", "Not significant" = "white", "NA" = "white")) +
  # scale_discrete_manual(aesthetics = "stroke", 
  #                       values =  c("<0.01" = 1, "Not significant" = 0, "NA" = 0), 
  #                       guide = "none") +
  scale_fill_gradientn(colors = colorRampPalette(c('#f9dbbd', '#ffa5ab', "#da627d", "#a53860", "#450920"))(100)) +
  facet_grid(. ~ factor(expCond.donor,
                        levels=c(569,589,602,662,673,678,548,616,633,640,681,699)), scales="free", space="free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "right",
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
        strip.text.y = element_blank(),
        strip.background.y = element_blank()) +
  labs(x = "Conditions",
       y = "Genes",
       size = "Percent expressed",
       fill = "Scaled expression")


pdf(paste0(excel.dir, deResDir, "GOEnrichment/HLA-DQ2_dotplot.pdf"), height=1.90, width=9.66)

print(p)

dev.off()








