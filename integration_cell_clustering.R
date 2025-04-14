library(Seurat)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(ggpubr) #cor
library(reshape2)
library(tibble) # for 5.2
library(tidyr) # for 5.2
library(corrplot) # for 5.2

seuratObj_annot <- readRDS('integration_2/integration_2_leiden/RDS_Dir/integration_2_leiden_annot.rds')

seuratObj_annot@meta.data %>%
  select(expCond.stimuli.time, expCond.media) %>%
  distinct() %>%
  arrange(expCond.media)

seuratObj_annot@meta.data %>%
  select(expCond.lib, expCond.media) %>%
  count(expCond.lib, expCond.media)

# if not running hippo, jump to 4.0!!!!
# for intergration_2_leiden, go from 4.3
#### 1.0 customize tsne and umap visualizations ####
# setup custom theme for plotting
theme1noLegend <- theme(plot.title = element_text(size = 16, hjust = 0.5),
                        # legend.key.size = unit(0.7, "cm"),
                        axis.title = element_text(size = 20),
                        axis.text = element_text(size = 25),
                        legend.position="bottom",
                        legend.text = element_text(size = 14) ) + Seurat::NoLegend()
theme1wLegend <- theme(plot.title = element_text(size = 16, hjust = 0.5),
                       # legend.key.size = unit(0.7, "cm"),
                       axis.title = element_text(size = 15),
                       axis.text = element_text(size = 25),
                       legend.position="bottom",
                       legend.text = element_text(size = 15) )
theme1wLegendRight <- theme(plot.title = element_text(size = 16, hjust = 0.5),
                            # legend.key.size = unit(0.7, "cm"),
                            axis.title = element_text(size = 20),
                            axis.text = element_text(size = 25),
                            legend.position="right",
                            legend.text = element_text(size = 14) )

# getClusterSummaryReplot(resDir = 'integration/integration_1/integration_1_louvain', 
#                         rds = seuratObj_annot, 
#                         newAnnotation = F, 
#                         expCondCheck = 'cluster_annotation', 
#                         expCondCheckFname = 'cluster_annotation')



selectedCol <- DiscretePalette(n = length(levels(Seurat::Idents(seuratObj_annot))), palette = 'alphabet')
workdir <- "integration_2/integration_2_Leiden"

# umap clusters with and without labeling and split by donors 
umapCluster <- DimPlot(seuratObj_annot, reduction = "umap", cols = selectedCol, label = T, label.size = 6, repel = T) + 
  labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')
umapClusterNolabel <- DimPlot(seuratObj_annot, reduction = "umap", cols = selectedCol, label = F, repel = T) + 
  labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')
umapSplit <- DimPlot(seuratObj_annot, reduction = "umap", cols = selectedCol, label = T, 
                     label.size = 4, repel = T, split.by = 'expCond.donor') +
  labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')

umapCluster.stimuli.time <- DimPlot(seuratObj_annot, reduction = "umap", group.by = "expCond.stimuli.time",
                       cols = selectedCol,
                       label = T, label.size = 6, repel = T) + 
  labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')

umapCluster.media <- DimPlot(seuratObj_annot, reduction = "umap", group.by = "expCond.media",
                                    cols = selectedCol,
                                    label = T, label.size = 6, repel = T) + 
  labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')

pdf(file = paste0(workdir, "/umap_plot_noLabel_integrate_orgAnnotation.pdf"), 
    width = 5.7, height = 7)
print(umapClusterNolabel + theme1wLegend + guides(colour = guide_legend(nrow=4, byrow=TRUE, override.aes = list(size=6))))
dev.off()

pdf(file = paste0(workdir, "/umap_plot_wLabel_integrate_orgAnnotation.pdf"), 
    width = 5.7, height = 7)
print(umapCluster + theme1wLegend + guides(colour = guide_legend(nrow=4, byrow=TRUE, override.aes = list(size=6))) )
dev.off()

pdf(file = paste0(workdir, "/umap_plot_wLabel_integrate_orgAnnotation_stimuli.time.pdf"), 
    width = 5.7, height = 7)
print(umapCluster.stimuli.time + theme1wLegend + guides(colour = guide_legend(nrow=4, byrow=TRUE, override.aes = list(size=6))) )
dev.off()

pdf(file = paste0(workdir, "/umap_plot_wLabel_integrate_orgAnnotation_media.pdf"), 
    width = 5.7, height = 7)
print(umapCluster.media + theme1wLegend + guides(colour = guide_legend(nrow=4, byrow=TRUE, override.aes = list(size=6))) )
dev.off()

plotName = paste(workdir, sprintf('/umap_plot_wLabel_orgAnnotation_%s.pdf', 'donor'), sep = '/')
if ( length(levels(as.factor(seuratObj_annot$expCond.donor))) > 1 ) {
  if ( length(levels(as.factor(seuratObj_annot$expCond.donor))) == 2) {
    pdf(file = plotName, width = 11, height = 7)
  } else if ( length(levels(as.factor(seuratObj_annot$expCond.donor))) == 3) {
    pdf(file = plotName, width = 13, height = 7)
  } else if ( length(levels(as.factor(seuratObj_annot$expCond.donor))) == 4) {
    pdf(file = plotName, width = 21, height = 7)
  } else if ( length(levels(as.factor(seuratObj_annot$expCond.donor))) > 4) {
    pdf(file = plotName, width = 5.5*length(levels(as.factor(seuratObj_annot$expCond.donor))), height = 7)
  }
  print(umapSplit + theme1wLegend + guides(colour = guide_legend(nrow=4, byrow=TRUE, override.aes = list(size=6))) )
  dev.off()
}



# tsne clusters with and without labeling and split by donors 
tsneCluster <- DimPlot(seuratObj_annot, reduction = "tsne", cols = selectedCol, label = T, label.size = 6, repel = T) + labs(title = 'tSNE clustering', x = "tSNE 1", y = 'tSNE 2')
tsneClusterNolabel <- DimPlot(seuratObj_annot, reduction = "tsne", cols = selectedCol, label = F, repel = T) + labs(title = 'tSNE clustering', x = "tSNE 1", y = 'tSNE 2')
tsneSplit <- DimPlot(seuratObj_annot, reduction = "tsne", cols = selectedCol, label = T, 
                     label.size = 4, repel = T, split.by = 'expCond.donor') +
  labs(title = 'tSNE clustering', x = "tSNE 1", y = 'tSNE 2')


pdf(file = paste0(workdir,"/tsne_plot_noLabel_integrate_orgAnnotation.pdf"), width = 5.7, height = 6.7)
print(tsneClusterNolabel + theme1wLegend + guides(colour = guide_legend(nrow=4, byrow=TRUE, override.aes = list(size=6))))
dev.off()

pdf(file = paste0(workdir,"/tsne_plot_wLabel_integrate_orgAnnotation.pdf"), width = 5.7, height = 6.7)
print(tsneCluster + theme1wLegend + guides(colour = guide_legend(nrow=4, byrow=TRUE, override.aes = list(size=6))) )
dev.off()

plotName = paste(workdir,sprintf('tsne_plot_wLabel_orgAnnotation_%s.pdf', 'donor'), sep = '/')
if ( length(levels(as.factor(seuratObj_annot$expCond.donor))) > 1 ) {
  if ( length(levels(as.factor(seuratObj_annot$expCond.donor))) == 2) {
    pdf(file = plotName, width = 11, height = 7)
  } else if ( length(levels(as.factor(seuratObj_annot$expCond.donor))) == 3) {
    pdf(file = plotName, width = 13, height = 7)
  } else if ( length(levels(as.factor(seuratObj_annot$expCond.donor))) == 4) {
    pdf(file = plotName, width = 21, height = 7)
  } else if ( length(levels(as.factor(seuratObj_annot$expCond.donor))) > 4) {
    pdf(file = plotName, width = 5.5*length(levels(as.factor(seuratObj_annot$expCond.donor))), height = 7)
  }
  # print(tsneSplit + theme1noLegend)
  print(tsneSplit + theme1wLegend + guides(colour = guide_legend(nrow=4, byrow=TRUE, override.aes = list(size=6))) )
  dev.off()
}

#### 2.0 T cells and NK cells | hippo  ####
getHippoRes <- function(resDir=NULL, rds=NULL, newAnnotation=F, newAnnotationRscriptName=NULL, expCondCheck='all',
                        cellcluster = NULL , expCond = NULL, noClusters = 3, sparseMatrix = F, initial.label.on = F,
                        hippoResNamePrefix = 'hippo_cluster_test', topN = 100, debug = as.logical(F)) {
  ##--------------------------------------------------------------------------------------##
  if (is.null(cellcluster)) stop("Please provide 'cellcluster' for hippo analysis")
  sel.expConds <- expCond ## to avoid 'expCond' otpion with metadata 'expCond' rename this paratmer into sel.expConds
  ## ----
  newAnnotation           <- as.logical(newAnnotation)
  if (newAnnotation & is.null(newAnnotationRscriptName)) stop("Option 'newAnnotation' is on, please provide corresponding option 'newAnnotationRscriptName'.")
  ##--------------------------------------------------------------------------------------##
  if (is.null(resDir) & !is.null(rds)) {
    if (class(rds)=='Seurat') {
      seuratObjFinal      <<- rds
      print('RDS is provided with rds option')
    } else {
      rdsFname            <- rds
      ## ---
      if (!file.exists(rdsFname)) stop("Please execute getClusterMarker() to conduct integration analysis before running getClusterSummaryReplot().")
      seuratObjFinal      <<- readRDS(file = as.character(rdsFname))
      print('Done for RDS read in')
    }
    resDir                <- getwd()
  } else if (!is.null(resDir) & is.null(rds)) {
    rdsFname              <- sprintf('%s/RDS_Dir/%s.rds', resDir, basename(resDir))
    resDir                <- resDir
    ## ---
    if (!file.exists(rdsFname)) stop("Please execute getClusterMarker() to conduct integration analysis before running getClusterSummaryReplot().")
    seuratObjFinal          <<- readRDS(file = as.character(rdsFname))
    print('Done for RDS read in')
  } else if (is.null(resDir) & is.null(rds)){
    stop("Error: please provide either option 'resDir' or 'rds', or both. ")
  } else if (!is.null(resDir) & !is.null(rds)){
    if (class(rds)=='Seurat') {
      seuratObjFinal      <<- rds
      print('RDS is provided with rds option')
    } else {
      rdsFname            <- rds
      ## ---
      if (!file.exists(rdsFname)) stop("Please execute getClusterMarker() to conduct integration analysis before running getClusterSummaryReplot().")
      seuratObjFinal      <<- readRDS(file = as.character(rdsFname))
      print('Done for RDS read in')
    }
    resDir                <- resDir
  }
  ##--------------------------------------------------------------------------------------##
  ## update results directory if new annotation is used
  resDir                <- paste(resDir, 'hippo_results', sep = '/')
  if (!dir.exists(resDir)) dir.create(resDir)
  resDir                <- paste(resDir, hippoResNamePrefix, sep = '/')
  if (!dir.exists(resDir)) dir.create(resDir)
  ##--------------------------------------------------------------------------------------##
  # if (expCondCheck == 'sample') {
  #   if (is.null(expCondCheckFname)) {
  #     expCondCheckFname        <- 'expCond_sample'
  #   } else {
  #     expCondCheckFname        <- expCondCheckFname
  #   }
  # } else {
  #   if (is.null(expCondCheckFname)) {
  #     expCondCheckFname        <- as.character(expCondCheck)
  #   } else {
  #     expCondCheckFname        <- expCondCheckFname
  #   }
  # }
  ##--------------------------------------------------------------------------------------##
  if (newAnnotation) {
    ## Assign cell type identity to clusters
    ## redefine the level of Idents on the y-axis can be adjusted here by inputting order for cell annotation
    source(as.character(newAnnotationRscriptName))
    print('-=-=-=-')
    print('updated Idents are as below:')
    print(table(Idents(seuratObjFinal)))
    print('-=-=-=-')
  }
  ##--------------------------------------------------------------------------------------##
  ## update 'seuratObjFinal@meta.data$expCond'
  if (expCondCheck == 'all') {
    seuratObjFinal                     <- seuratObjFinal
  } else {
    if (!expCondCheck%in%colnames(seuratObjFinal@meta.data)) {
      stop("ERROR: 'expCondCheck' does not exist in your 'rds' metadata.")
    } else {
      seuratObjFinal@meta.data$expCond <- seuratObjFinal@meta.data[, grep(sprintf('^%s$', as.character(expCondCheck)), colnames(seuratObjFinal@meta.data))]
    }
  }
  ##--------------------------------------------------------------------------------------##
  ## if provided, subset on 'cellcluster'
  orgClusterLevels <- levels(Seurat::Idents(seuratObjFinal))
  if (debug) print(sprintf("%s idents levels are %s", length(orgClusterLevels), paste(orgClusterLevels, collapse = ', ')))
  if (!is.null(cellcluster)) {
    if (any(!cellcluster %in% orgClusterLevels ) ) stop('Please provide the corresponding cell clusters in identfied idents for hippo analysis.')
    print(sprintf('Subsetting %s specific cell clusters: %s', length(cellcluster), paste(cellcluster, collapse = ',')))
    seuratObjFinal        <- subset(seuratObjFinal, idents = as.character(cellcluster) )
  }
  ## in provided, subset on 'expCond'
  expCondLevels           <- levels(factor(seuratObjFinal@meta.data$expCond))
  if (debug) print(sprintf("%s expCond levels are %s", length(expCondLevels), paste(expCondLevels, collapse = ', ')))
  if (!is.null(sel.expConds)) {
    if (any(!sel.expConds %in% expCondLevels ) ) stop("Please provide the corresponding experimental condition levels specified in 'expCondCheck' option.")
    print(sprintf('Subsetting %s specific expCond: %s', length(sel.expConds), paste(sel.expConds, collapse = ',')))
    
    if (length(sel.expConds)==1) {
      seuratObjFinal          <- subset(seuratObjFinal, expCond == sel.expConds)
    } else {
      for (i in 1:length(sel.expConds)) {
        seuratObjFinalPrep    <- subset(seuratObjFinal, subset = expCond == as.character(sel.expConds[i]))
        if (i ==1) {
          seuratObjFinalPrep2 = seuratObjFinalPrep
        } else {
          seuratObjFinalPrep2 <- merge(seuratObjFinalPrep2, seuratObjFinalPrep)
        }
      }
      seuratObjFinal          <- seuratObjFinalPrep2
    }
  }
  # print(sprintf("hippo input is as below"))
  # seuratObjFinal
  ##--------------------------------------------------------------------------------------##
  inputDataPrep <- seuratObjFinal
  if (sparseMatrix) {
    print("lighthippo input is a sparse matrix")
    inputData               <- inputDataPrep@assays$RNA$counts ## sparse matrix
  } else {
    print("lighthippo input is a dense matrix")
    inputData               <- as.matrix(inputDataPrep@assays$RNA$counts) ##dense matrix
  }
  # print(head(inputData[,1:4]))
  # print(dim(inputData))
  # print('7574655647578282-=-=-=-=-=-=-=-=-=-')
  print(sprintf('lightHippo: %s cells in combined cell cluster (%s) with a total of %s genes expressed', dim(inputData)[2], paste(cellcluster, collapse = '; '), dim(inputData)[1] ))
  print('-=-=-=-=-=')
  ##--------------------------------------------------------------------------------------##
  ## 1. run light hippo, returned
  systime1         <- Sys.time()
  if (initial.label.on) {
    inputData.cluster  <- as.numeric(Idents(inputDataPrep))
    print(sprintf("Processing analysis with 'initial.labels' on for 'K.round' = %s", noClusters))
    set.seed(1234)
    K.round            <- noClusters - length(unique(inputData.cluster)) + 1
    lightHippoRes      <- lightHippo::lightHIPPO(dat = inputData, K.round = K.round, initial.labels = inputData.cluster)
    final_feature_list <- lightHippo::organizing_hippo_features(lightHippoRes)
    save(inputData, lightHippoRes, inputData.cluster, noClusters, cellcluster, hippoResNamePrefix, file = file.path(resDir, sprintf('lightHIPPOres_%s_k%s.Rdata', hippoResNamePrefix, noClusters) ) )
  } else {
    print(sprintf("Processing analysis with 'initial.labels' off for 'K.round' = %s", noClusters))
    set.seed(1234)
    lightHippoRes      <- lightHippo::lightHIPPO(dat = inputData, K.round = noClusters, initial.round = 0)
    final_feature_list <- lightHippo::organizing_hippo_features(lightHippoRes)
    save(inputData, lightHippoRes, noClusters, cellcluster, hippoResNamePrefix, file = file.path(resDir, sprintf('lightHIPPOres_%s_k%s.Rdata', hippoResNamePrefix, noClusters) ) )
  }
  systime2         <- Sys.time()
  print(sprintf("lighthippo complete used %s %s.", round(difftime(systime2, systime1), digits = 2), attr(difftime(systime2, systime1), "units") ))
  print("lighthippo next round ID is shown as below:")
  print(table(lightHippoRes$next_round_IDs))
  print('-=-=-=-')
  ##--------------------------------------------------------------------------------------##
  ## 2. make diagnostic plot
  total.num.gene  <- nrow(inputData)
  set.seed(20200610)
  randomIDs       <- sample(1:total.num.gene, 5000)
  summarizing_dat <- lightHippo::summarize_current_zero_proportions(inputData[randomIDs, ], lightHippoRes$next_round_IDs)
  plot_dat_per_cluster_inflation <- lightHippo::visualize_current_zero_proportions(summarizing_dat)
  if (noClusters <= 3) {
    pdf(file = file.path(resDir, sprintf('diagnosticPlot_cellCluster_%s_k%s.pdf', hippoResNamePrefix, noClusters) ), width = 2*(noClusters+1), height = ceiling((noClusters+1)/4)*3)
  } else {
    pdf(file = file.path(resDir, sprintf('diagnosticPlot_cellCluster_%s_k%s.pdf', hippoResNamePrefix, noClusters) ), width = 8, height = ceiling(noClusters/3.9)*3)
  }
  # pdf("lightHIPPO_counts_inflation_check.png")
  print(plot_dat_per_cluster_inflation)
  dev.off()
  ##--------------------------------------------------------------------------------------##
  ttID       <- lightHippo::cut_hierarchy(lightHippoRes, K = noClusters)
  print(sprintf("%s clusters table is:", noClusters))
  print(table(ttID))
  print('-=-=-=-')
  check_these <- final_feature_list$ID[[noClusters-1]]
  tt.selected <- lightHippo::summarize_for_feature_dot(inputData[check_these , ], ttID)
  p <-  lightHippo::makeDotplot(input = tt.selected, topN = topN)
  if (noClusters <= 4) {
    ggsave(filename = file.path(resDir, sprintf('topFeatures_cellCluster_%s_k%s_dotplot.pdf', hippoResNamePrefix, noClusters) ), plot = p, width = ceiling(topN/4), height = 5)
  } else if (noClusters>4 & noClusters <=6) {
    ggsave(filename = file.path(resDir, sprintf('topFeatures_cellCluster_%s_k%s_dotplot.pdf', hippoResNamePrefix, noClusters) ), plot = p, width = ceiling(topN/4), height = noClusters)
  } else {
    ggsave(filename = file.path(resDir, sprintf('topFeatures_cellCluster_%s_k%s_dotplot.pdf', hippoResNamePrefix, noClusters) ), plot = p, width = ceiling(topN/4), height = ceiling(noClusters/1.5))
  }
  ##--------------------------------------------------------------------------------------##
  viewRes  <- lightHippo::cut_hierarchy(lightHippoRes, K = noClusters, cut_sequence = T)
  pdf(file = file.path(resDir, sprintf('%s_k%s_clustersTree.pdf', hippoResNamePrefix, noClusters) ), width = 5, height = 4 )
  lightHippo::visualize_hippo_hierarchy(viewRes)
  dev.off()
  ##--------------------------------------------------------------------------------------##
  if (initial.label.on) {
    noStart = length(unique(inputData.cluster))
  } else {
    noStart = 2
  }
  for (j in noStart:noClusters) {
    res          <- lightHippo::cut_hierarchy(lightHippoRes, K = j)
    cellHippoPer <- as.data.frame( table(res))
    colnames(cellHippoPer) <- c('hippo', 'No cells')
    fname <- file.path(resDir, sprintf('Per_cellCluster_%s_k%s.xlsx', hippoResNamePrefix, noClusters) )
    # print(cellHippoPer)
    if (j ==1) {
      xlsx::write.xlsx(x = cellHippoPer, file = fname, sheetName = paste('k', j, sep = '_'), row.names = F, append = F)
    } else {
      xlsx::write.xlsx(x = cellHippoPer, file = fname, sheetName = paste('k', j, sep = '_'), row.names = F, append = T)
    }
  }
  print(sprintf("light Hippo analysis complete for k=%s", noClusters))
  print('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-')
  ##--------------------------------------------------------------------------------------##
}

workdir <- "integration/integration_1/integration_1_louvain/"

getHippoRes(resDir = workdir,
            rds = seuratObj_annot,
            cellcluster = "NK cells",
            noClusters = 3,
            initial.label.on = F,
            sparseMatrix = F,
            hippoResNamePrefix = 'hippo_NK')

getHippoRes(resDir = workdir,
            rds = seuratObj_annot,
            cellcluster = "CD4 T cells",
            noClusters = 3,
            initial.label.on = F,
            sparseMatrix = F,
            hippoResNamePrefix = 'hippo_CD4T')

getHippoRes(resDir = workdir,
            rds = seuratObj_annot,
            cellcluster = "CD8 T cells",
            noClusters = 3,
            initial.label.on = F,
            sparseMatrix = F,
            hippoResNamePrefix = 'hippo_CD8T')


# update idents by hippo clusters
hippoRds  <- updateHippoIdents(
  resDir = workdir,
  rds = seuratObj_annot,
  newAnnotation = F, 
  hippoResList = list(
    'NK cells' = paste0(workdir, 'hippo_results/hippo_NK/lightHIPPOres_hippo_NK_k3.Rdata'),
    'CD4 T cells' =  paste0(workdir, 'hippo_results/hippo_CD4T/lightHIPPOres_hippo_CD4T_k3.Rdata'),
    'CD8 T cells' =  paste0(workdir, 'hippo_results/hippo_CD8T/lightHIPPOres_hippo_CD8T_k3.Rdata')
  ),
  hippoResK = c(3, 3, 3)
)

#### 3.0 after 2.0, customize tsne and umap visualizations ####


Tcells <- unique(Idents(hippoRds))[grepl("T cells", unique(Idents(hippoRds)))]
NKcells <- unique(Idents(hippoRds))[grepl("NK cells", unique(Idents(hippoRds)))]


hippoRds_sub <- subset(hippoRds, idents = c(Tcells, NKcells))

# setup custom theme for plotting
theme1noLegend <- theme(plot.title = element_text(size = 16, hjust = 0.5),
                        # legend.key.size = unit(0.7, "cm"),
                        axis.title = element_text(size = 20),
                        axis.text = element_text(size = 25),
                        legend.position="bottom",
                        legend.text = element_text(size = 14) ) + Seurat::NoLegend()
theme1wLegend <- theme(plot.title = element_text(size = 16, hjust = 0.5),
                       # legend.key.size = unit(0.7, "cm"),
                       axis.title = element_text(size = 15),
                       axis.text = element_text(size = 25),
                       legend.position="bottom",
                       legend.text = element_text(size = 15) )
theme1wLegendRight <- theme(plot.title = element_text(size = 16, hjust = 0.5),
                            # legend.key.size = unit(0.7, "cm"),
                            axis.title = element_text(size = 20),
                            axis.text = element_text(size = 25),
                            legend.position="right",
                            legend.text = element_text(size = 14) )

# getClusterSummaryReplot(resDir = 'integration/integration_1/integration_1_louvain', 
#                         rds = seuratObj_annot, 
#                         newAnnotation = F, 
#                         expCondCheck = 'cluster_annotation', 
#                         expCondCheckFname = 'cluster_annotation')



selectedCol <- DiscretePalette(n = length(levels(Seurat::Idents(hippoRds_sub))), palette = 'glasbey')
workdir <- "integration/integration_1/integration_1_louvain"
hippodir <- "hippo_results"

# umap clusters with and without labeling and split by donors 
umapCluster <- DimPlot(hippoRds_sub, reduction = "umap", cols = selectedCol, label = T, label.size = 6, repel = T) + 
  labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')
umapClusterNolabel <- DimPlot(hippoRds_sub, reduction = "umap", cols = selectedCol, label = F, repel = T) + 
  labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')
umapSplit <- DimPlot(hippoRds_sub, reduction = "umap", cols = selectedCol, label = T, 
                     label.size = 4, repel = T, split.by = 'expCond.donor') +
  labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')

umapCluster.stimuli.time <- DimPlot(hippoRds_sub, reduction = "umap", group.by = "expCond.stimuli.time",
                                    cols = selectedCol,
                                    label = T, label.size = 6, repel = T) +
  labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')

pdf(file = paste0(workdir, "/", hippodir, "/umap_plot_noLabel_integrate_orgAnnotation.pdf"), 
    width = 5.7, height = 7)
print(umapClusterNolabel + theme1wLegend + guides(colour = guide_legend(nrow=4, byrow=TRUE, override.aes = list(size=6))))
dev.off()

pdf(file = paste0(workdir, "/", hippodir, "/umap_plot_wLabel_integrate_orgAnnotation.pdf"), 
    width = 5.7, height = 7)
print(umapCluster + theme1wLegend + guides(colour = guide_legend(nrow=4, byrow=TRUE, override.aes = list(size=6))) )
dev.off()

pdf(file = paste0(workdir, "/", hippodir, "/umap_plot_wLabel_integrate_orgAnnotation_stimuli.time.pdf"), 
    width = 5.7, height = 7)
print(umapCluster.stimuli.time + theme1wLegend + guides(colour = guide_legend(nrow=4, byrow=TRUE, override.aes = list(size=6))) )
dev.off()

plotName = paste(workdir, hippodir, sprintf('/umap_plot_wLabel_orgAnnotation_%s.pdf', 'donor'), sep = '/')
if ( length(levels(as.factor(hippoRds_sub$expCond.donor))) > 1 ) {
  if ( length(levels(as.factor(hippoRds_sub$expCond.donor))) == 2) {
    pdf(file = plotName, width = 11, height = 7)
  } else if ( length(levels(as.factor(hippoRds_sub$expCond.donor))) == 3) {
    pdf(file = plotName, width = 13, height = 7)
  } else if ( length(levels(as.factor(hippoRds_sub$expCond.donor))) == 4) {
    pdf(file = plotName, width = 21, height = 7)
  } else if ( length(levels(as.factor(hippoRds_sub$expCond.donor))) > 4) {
    pdf(file = plotName, width = 5.5*length(levels(as.factor(hippoRds_sub$expCond.donor))), height = 7)
  }
  print(umapSplit + theme1wLegend + guides(colour = guide_legend(nrow=4, byrow=TRUE, override.aes = list(size=6))) )
  dev.off()
}



# tsne clusters with and without labeling and split by donors 
tsneCluster <- DimPlot(hippoRds_sub, reduction = "tsne", cols = selectedCol, label = T, label.size = 6, repel = T) + labs(title = 'tSNE clustering', x = "tSNE 1", y = 'tSNE 2')
tsneClusterNolabel <- DimPlot(hippoRds_sub, reduction = "tsne", cols = selectedCol, label = F, repel = T) + labs(title = 'tSNE clustering', x = "tSNE 1", y = 'tSNE 2')
tsneSplit <- DimPlot(hippoRds_sub, reduction = "tsne", cols = selectedCol, label = T, 
                     label.size = 4, repel = T, split.by = 'expCond.donor') +
  labs(title = 'tSNE clustering', x = "tSNE 1", y = 'tSNE 2')


pdf(file = paste0(workdir, "/", hippodir, "/tsne_plot_noLabel_integrate_orgAnnotation.pdf"), width = 5.7, height = 6.7)
print(tsneClusterNolabel + theme1wLegend + guides(colour = guide_legend(nrow=4, byrow=TRUE, override.aes = list(size=6))))
dev.off()

pdf(file = paste0(workdir, "/", hippodir, "/tsne_plot_wLabel_integrate_orgAnnotation.pdf"), width = 5.7, height = 6.7)
print(tsneCluster + theme1wLegend + guides(colour = guide_legend(nrow=4, byrow=TRUE, override.aes = list(size=6))) )
dev.off()

plotName = paste(workdir, hippodir, sprintf('tsne_plot_wLabel_orgAnnotation_%s.pdf', 'donor'), sep = '/')
if ( length(levels(as.factor(hippoRds_sub$expCond.donor))) > 1 ) {
  if ( length(levels(as.factor(hippoRds_sub$expCond.donor))) == 2) {
    pdf(file = plotName, width = 11, height = 7)
  } else if ( length(levels(as.factor(hippoRds_sub$expCond.donor))) == 3) {
    pdf(file = plotName, width = 13, height = 7)
  } else if ( length(levels(as.factor(hippoRds_sub$expCond.donor))) == 4) {
    pdf(file = plotName, width = 21, height = 7)
  } else if ( length(levels(as.factor(hippoRds_sub$expCond.donor))) > 4) {
    pdf(file = plotName, width = 5.5*length(levels(as.factor(hippoRds_sub$expCond.donor))), height = 7)
  }
  # print(tsneSplit + theme1noLegend)
  print(tsneSplit + theme1wLegend + guides(colour = guide_legend(nrow=4, byrow=TRUE, override.aes = list(size=6))) )
  dev.off()
}


#### if cell type per cluster, go with 4.1 customize tsne and umap visualizations and subplots for each treatment | dir majority_voting ####
workdir <- "integration_2/integration_2_leiden/cell_clustering"

cell_types <- c("Mono Mph",
                "B", 
                "CD4 T", 
                "CD8 T", 
                "NK",
                "DC")


# selectedCol <- c("#f77f00", "#075a84", "#9163cb", "#ff4d6d", "#679436", "#f07167", "#603808", "#07beb8")
# 
# selectedCol <- c("#ff9b54", "#6d98ba", "#d4c2fc", "#ff758f", "#98c9a3", "#ffb4a2", "#ddb892", "#b2f7ef")
# 


selectedCol <- c("#075a84", "#9163cb", "#ff4d6d", "#679436", "#f77f00","#07beb8")

selectedCol <- c("#6d98ba", "#d4c2fc", "#ff758f", "#98c9a3", "#ff9b54", "#b2f7ef")


names(selectedCol) <- cell_types

# umap clustering
pdf(file = paste0(workdir, "/umap_all.pdf"), width = 5.38, height = 4.14)

DimPlot(seuratObj_annot, reduction = "umap", cols = selectedCol, label = T, label.size = 4, repel = T) + 
  labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')

dev.off()


pdf(file = paste0(workdir, "/umap_expCond.stimuli.time.pdf"), width = 21, height = 3.5)

DimPlot(seuratObj_annot, reduction = "umap", cols = selectedCol, label = F, label.size = 4, repel = T, split.by = 'expCond.stimuli.time') +
  labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')

dev.off()

pdf(file = paste0(workdir, "/umap_expCond.media.pdf"), width = 7, height = 3.5)

DimPlot(seuratObj_annot, reduction = "umap", cols = selectedCol, label = F, label.size = 4, repel = T, split.by = 'expCond.media') +
  labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')

dev.off()

# umap clustering
pdf(file = paste0(workdir, "/umap_all.darker.pdf"), width = 5.38, height = 4.14)

DimPlot(seuratObj_annot, reduction = "umap", cols = selectedCol, label = T, label.size = 4, repel = T) + 
  labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')

dev.off()


pdf(file = paste0(workdir, "/umap_expCond.stimuli.time.darker.pdf"), width = 21, height = 3.5)

DimPlot(seuratObj_annot, reduction = "umap", cols = selectedCol, label = F, label.size = 4, repel = T, split.by = 'expCond.stimuli.time') +
  labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')

dev.off()


#### if cell type per cell, go with 4.2 customize tsne and umap visualizations and subplots for each treatment | dir per_cell_and_majority_voting ####
workdir <- "integration_2/integration_2_leiden/cell_clustering/per_cell_and_majority_voting"

# unique(seuratObj_annot@meta.data$predicted_labels)
seuratObj_annot@meta.data$cluster_annotation2 <- seuratObj_annot@meta.data$predicted_labels


cell_labels <- c(
  "EC general capillary" = "EC",
  
  "B cells" = "B",
  
  "CD4 T cells" = "CD4 T",
  "CD8 T cells" = "CD8 T",
  
  "NK cells" = "NK",
  
  "Monocyte-derived Mph" = "Mono/Mph",
  
  "Alveolar Mph CCL3+" = "Mono/Mph",
  "Alveolar Mph MT-positive" = "Mono/Mph",
  
  "Classical monocytes" = "Mono/Mph",
  "Non-classical monocytes" = "Mono/Mph",
  
  "Migratory DCs" = "DC"
)


cluster_labels <- c(
  `2` = "B", `3` = "B", `13` = "B",
  `9` = "NK",
  `6` = "CD4 T", `12` = "CD4 T", `15` = "CD4 T",
  `5` = "CD8 T",
  `16` = "DC",
  `1` = "Mono/Mph", `4` = "Mono/Mph", `7` = "Mono/Mph",
  `8` = "Mono/Mph", `10` = "Mono/Mph", `11` = "Mono/Mph",
  `14` = "Mono/Mph", `17` = "Mono/Mph"
)


seuratObj_annot@meta.data$cluster_annotation2 <- mapply(
  function(annotation, cluster) {
    cluster <- as.character(cluster)
    if (annotation %in% names(cell_labels)) {
      # Replace using cell_labels
      cell_labels[annotation]
    } else if (cluster %in% names(cluster_labels)) {
      # Replace using cluster_labels
      cluster_labels[[cluster]]
    } else {
      # If no match, keep original annotation
      annotation
    }
  },
  seuratObj_annot@meta.data$cluster_annotation2,
  seuratObj_annot@meta.data$seurat_clusters
)

Idents(seuratObj_annot) <- seuratObj_annot@meta.data$cluster_annotation2

cell_types <- c(
  "Mono/Mph", # From "Monocyte-derived Mph" "Classical monocytes" and "Non-classical monocytes" "Alveolar Mph CCL3+" and "Alveolar Mph MT-positive"
  "B",        # From "B cells"
  "CD4 T",    # From "CD4 T cells"
  "CD8 T",    # From "CD8 T cells"
  "NK",       # From "NK cells"
  "DC",       # From "Migratory DCs"
  "EC"       # From "EC general capillary"
)


saveRDS(seuratObj_annot, 'integration_2/integration_2_leiden/RDS_Dir/integration_2_leiden_annot.rds')
names(table(Idents(seuratObj_annot)))


# selectedCol <- DiscretePalette(n = length(levels(Seurat::Idents(seuratObj_annot))), palette = 'parade')

selectedCol<-c( 'cornflowerblue', 
                'mediumvioletred',
                'wheat3', 
                'lightsalmon',  
                "#7CE3D8", 
                'yellowgreen',
                'seagreen')
# ,
#                 'indianred1',
#                 "tan1",  '#DDAD4B',
#                 'seagreen1',  'seagreen',
#                 "darkolivegreen1" ,"tan2" , "cadetblue1", "darkolivegreen4")

# selectedCol <- c("#f77f00", "#075a84", "#9163cb", "#ff4d6d", "#679436", "#f07167", "#603808", "#07beb8")
# 
# selectedCol <- c("#ff9b54", "#6d98ba", "#d4c2fc", "#ff758f", "#98c9a3", "#ffb4a2", "#ddb892", "#b2f7ef")
# 
# 
# 
# selectedCol <- c("#075a84", "#9163cb", "#ff4d6d", "#679436", "#f77f00","#07beb8")
# 
# selectedCol <- c("#6d98ba", "#d4c2fc", "#ff758f", "#98c9a3", "#ff9b54", "#b2f7ef")


names(selectedCol) <- cell_types


# umap clustering
pdf(file = paste0(workdir, "/umap_all.pdf"), width = 5.38, height = 4.14)

DimPlot(seuratObj_annot, reduction = "umap", cols = selectedCol, label = T, label.size = 4, repel = T) + 
  labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')

dev.off()


pdf(file = paste0(workdir, "/umap_expCond.stimuli.time.pdf"), width = 21, height = 3.5)

DimPlot(seuratObj_annot, reduction = "umap", cols = selectedCol, label = F, label.size = 4, repel = T, split.by = 'expCond.stimuli.time') +
  labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')

dev.off()

pdf(file = paste0(workdir, "/umap_expCond.media.pdf"), width = 7, height = 3.5)

DimPlot(seuratObj_annot, reduction = "umap", cols = selectedCol, label = F, label.size = 4, repel = T, split.by = 'expCond.media') +
  labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')

dev.off()

# umap clustering
pdf(file = paste0(workdir, "/umap_all.darker.pdf"), width = 5.38, height = 4.14)

DimPlot(seuratObj_annot, reduction = "umap", cols = selectedCol, label = T, label.size = 4, repel = T) + 
  labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')

dev.off()


pdf(file = paste0(workdir, "/umap_expCond.stimuli.time.darker.pdf"), width = 21, height = 3.5)

DimPlot(seuratObj_annot, reduction = "umap", cols = selectedCol, label = F, label.size = 4, repel = T, split.by = 'expCond.stimuli.time') +
  labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')

dev.off()



#### if cell type per cell, go with 4.3 customize tsne and umap visualizations and subplots for each treatment | dir per_cell_and_majority_voting_curation ####
workdir <- "integration_2/integration_2_leiden/cell_clustering/per_cell_and_majority_voting"

# unique(seuratObj_annot@meta.data$predicted_labels)
seuratObj_annot@meta.data$cluster_annotation2 <- seuratObj_annot@meta.data$predicted_labels


cell_labels <- c(
  "EC general capillary" = "EC",
  
  "B cells" = "B",
  
  "CD4 T cells" = "CD4 T",
  "CD8 T cells" = "CD8 T",
  
  "NK cells" = "NK",
  
  "Monocyte-derived Mph" = "Mono/Mph",
  
  "Alveolar Mph CCL3+" = "Mono/Mph",
  "Alveolar Mph MT-positive" = "Mono/Mph",
  
  "Classical monocytes" = "Mono/Mph",
  "Non-classical monocytes" = "Mono/Mph",
  
  "Migratory DCs" = "DC"
)


cluster_labels <- c(
  `2` = "B", `3` = "B", `13` = "B",
  `9` = "NK",
  `6` = "CD4 T", `12` = "CD4 T", `15` = "CD4 T",
  `5` = "CD8 T",
  `16` = "DC",
  `1` = "Mono/Mph", `4` = "Mono/Mph", `7` = "Mono/Mph",
  `8` = "Mono/Mph", `10` = "Mono/Mph", `11` = "Mono/Mph",
  `14` = "Mono/Mph", `17` = "Mono/Mph"
)


seuratObj_annot@meta.data$cluster_annotation2 <- mapply(
  function(annotation, cluster) {
    cluster <- as.character(cluster)
    if (annotation %in% names(cell_labels)) {
      # Replace using cell_labels
      cell_labels[annotation]
    } else if (cluster %in% names(cluster_labels)) {
      # Replace using cluster_labels
      cluster_labels[[cluster]]
    } else {
      # If no match, keep original annotation
      annotation
    }
  },
  seuratObj_annot@meta.data$cluster_annotation2,
  seuratObj_annot@meta.data$seurat_clusters
)

Idents(seuratObj_annot) <- seuratObj_annot@meta.data$cluster_annotation2

cell_types <- c(
  "Mono/Mph", # From "Monocyte-derived Mph" "Classical monocytes" and "Non-classical monocytes" "Alveolar Mph CCL3+" and "Alveolar Mph MT-positive"
  "B",        # From "B cells"
  "CD4 T",    # From "CD4 T cells"
  "CD8 T",    # From "CD8 T cells"
  "NK",       # From "NK cells"
  "DC",       # From "Migratory DCs"
  "EC"       # From "EC general capillary"
)


names(table(Idents(seuratObj_annot)))


# selectedCol <- DiscretePalette(n = length(levels(Seurat::Idents(seuratObj_annot))), palette = 'parade')

selectedCol<-c( 'cornflowerblue', 
                'mediumvioletred',
                'wheat3', 
                'lightsalmon',  
                "#7CE3D8", 
                'yellowgreen',
                'seagreen')


names(selectedCol) <- cell_types


seuratObj_annot$expCond.stimuli.time <- factor(
  seuratObj_annot$expCond.stimuli.time,
  levels = c("unstim_0h", "Ig_4h", "Ig_18h", "CD3_4h", "CD3_18h", "LPS_4h", "LPS_18h")
)


# umap clustering
pdf(file = paste0(workdir, "/umap_all.pdf"), width = 5.38, height = 4.14)

DimPlot(seuratObj_annot, reduction = "umap", cols = selectedCol, label = T, label.size = 4, repel = T) + 
  labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')

dev.off()


pdf(file = paste0(workdir, "/umap_expCond.stimuli.time.pdf"), width = 18, height = 3.5)

DimPlot(seuratObj_annot, reduction = "umap", cols = selectedCol, label = F, label.size = 4, repel = T, split.by = 'expCond.stimuli.time') +
  labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')

dev.off()


# pdf(file = paste0(workdir, "/umap_expCond.stimuli.time.unstim_CD3.pdf"), width = 9, height = 3.5)
# 
# DimPlot(subset(seuratObj_annot, subset = expCond.stimuli.time %in%  c("unstim_0h", "CD3_4h", "CD3_18h")), reduction = "umap", 
#         cols = selectedCol, label = F, label.size = 4, repel = T, split.by = 'expCond.stimuli.time') +
#   labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')
# 
# dev.off()

pdf(file = paste0(workdir, "/umap_expCond.media.pdf"), width = 7, height = 3.5)

DimPlot(seuratObj_annot, reduction = "umap", cols = selectedCol, label = F, label.size = 4, repel = T, split.by = 'expCond.media') +
  labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')

dev.off()

# umap clustering
pdf(file = paste0(workdir, "/umap_all.darker.pdf"), width = 5.38, height = 4.14)

DimPlot(seuratObj_annot, reduction = "umap", cols = selectedCol, label = T, label.size = 4, repel = T) + 
  labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')

dev.off()


pdf(file = paste0(workdir, "/umap_expCond.stimuli.time.darker.pdf"), width = 16, height = 3.5)

DimPlot(seuratObj_annot, reduction = "umap", cols = selectedCol, label = F, label.size = 4, repel = T, split.by = 'expCond.stimuli.time') +
  labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')

dev.off()



# CD4T: UMAP 1 >0 mark them as CD4T fp
# CD8 T: UMAP 1 >0  mark them as CD8T fp
# CD8 T: UMAP 1 >0  mark them as CD8T fp
# B: seurat cluster only 2 3 13 are considered as B, all else B mark as B fp 
workdir2 <- "integration_2/integration_2_leiden/cell_clustering/per_cell_and_majority_voting_curation"

# table(seuratObj_annot@meta.data$cluster_annotation3)
# seuratObj_annot@meta.data %>%
#   filter(cluster_annotation3 == "B") %>%
#   filter(rownames(.) %in% rownames(umap_coords[umap_coords$umap_2< -5, ])) %>%
#   filter(rownames(.) %in% rownames(umap_coords[umap_coords$umap_1< 0, ])) %>%
#   dim()

seuratObj_annot@meta.data$cluster_annotation3 <-seuratObj_annot@meta.data$cluster_annotation2

umap_coords <- as.data.frame(Embeddings(seuratObj_annot, "umap"))

seuratObj_annot@meta.data <- seuratObj_annot@meta.data %>%
  mutate(cluster_annotation3 = ifelse(cluster_annotation2 == "B" & !(seurat_clusters %in% c(2, 3, 13)), 
                                      "B fp", cluster_annotation3),
         cluster_annotation3 = ifelse(cluster_annotation2 == "CD4 T" &
                                               !(rownames(.) %in% rownames(umap_coords[umap_coords$umap_1 < 0, ])), 
                                                 "CD4 T fp", cluster_annotation3),
         cluster_annotation3 = ifelse(cluster_annotation2 == "CD8 T" &
                                        !(rownames(.) %in% rownames(umap_coords[umap_coords$umap_1 < 0, ])), 
                                      "CD8 T fp", cluster_annotation3),
         cluster_annotation3 = ifelse(cluster_annotation2 == "NK" &
                                        !(rownames(.) %in% rownames(umap_coords[umap_coords$umap_1 < 0, ])), 
                                      "NK fp", cluster_annotation3),
         # cluster 15 being B CD4T misc
         cluster_annotation3 = ifelse(seurat_clusters %in% 15, 
                                      "B CD4 T misc", cluster_annotation3))

Idents(seuratObj_annot) <- seuratObj_annot@meta.data$cluster_annotation3


selectedCol[c("B fp", "CD4 T fp",  "CD8 T fp", "NK fp", "B CD4 T misc")] <- "#4C4C4C"

pdf(file = paste0(workdir2, "/umap_all.pdf"), width = 5.38, height = 4.14)

DimPlot(seuratObj_annot, reduction = "umap", cols = selectedCol, label = T, label.size = 4, repel = T) + 
  labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')

dev.off()



seuratObj_annot$expCond.stimuli.time <- factor(
  seuratObj_annot$expCond.stimuli.time,
  levels = c("unstim_0h", "Ig_4h", "Ig_18h", "CD3_4h", "CD3_18h", "LPS_4h", "LPS_18h")
)

pdf(file = paste0(workdir2, "/umap_expCond.stimuli.time.pdf"), width = 18, height = 3.5)

DimPlot(seuratObj_annot, reduction = "umap", cols = selectedCol, label = F, label.size = 4, repel = T, split.by = 'expCond.stimuli.time') +
  labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')

dev.off()


pdf(file = paste0(workdir2, "/umap_expCond.stimuli.time.unstim_CD3.pdf"), width = 7, height = 3.5)

DimPlot(subset(seuratObj_annot, subset = expCond.stimuli.time %in%  c("unstim_0h", "CD3_4h", "CD3_18h")), reduction = "umap", 
        cols = selectedCol, label = F, label.size = 4, repel = T, split.by = 'expCond.stimuli.time') +
  labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')

dev.off()


table(seuratObj_annot@meta.data$cluster_annotation3)
table(Idents(seuratObj_annot)) 

saveRDS(seuratObj_annot, 'integration_2/integration_2_leiden/RDS_Dir/integration_2_leiden_annot.rds')




















#### 5.0 correlation heatmap #####
# 5.1 randi prepare input-----
setwd("/gpfs/data/schoettler-lab/ftian1/robi")
.libPaths(c("/ess/home/home1/ftian1/R/x86_64-pc-linux-gnu-library/4.2", 
            "/gpfs/data/icelake-apps/software/gcc-12.1.0/R/4.2.1/lib64/R/library"))

# add libraries here

workdir <- "integration_2/integration_2_leiden/cell_clustering/per_cell_and_majority_voting_curation"

seuratObj_annot <- readRDS('integration_2/integration_2_leiden/RDS_Dir/integration_2_leiden_annot.rds')

cell_types <- c(
  "Mono/Mph", # From "Monocyte-derived Mph" "Classical monocytes" and "Non-classical monocytes" "Alveolar Mph CCL3+" and "Alveolar Mph MT-positive"
  "B",        # From "B cells"
  "CD4 T",    # From "CD4 T cells"
  "CD8 T",    # From "CD8 T cells"
  "NK",       # From "NK cells"
  "DC",       # From "Migratory DCs"
  "EC"       # From "EC general capillary"
)

seuratObj_annot <- subset(seuratObj_annot, idents = cell_types)

seuratObj_annot$expCond.celltype.stimuli.time <- paste0(Idents(seuratObj_annot), " ", seuratObj_annot$expCond.stimuli.time)
combinations <- combn(names(table(seuratObj_annot$expCond.celltype.stimuli.time)), 2, simplify = FALSE)
comp_group_list <- lapply(combinations, as.list)

expCondCheck <- 'expCond.celltype.stimuli.time'
correlation_values_df <- data.frame()

getScatterPlot <- function (resDir = NULL, rds = NULL, newAnnotation = F, newAnnotationRscriptName = NULL, 
                            subsetOn = NULL, subset = NULL, subsetGenes = NULL, expCondCheck = "idents", 
                            selectedGroups, interactive = T) 
{
  if (length(selectedGroups) != 2) 
    stop("Please provide corresponding x/y groups used for scatter plot")
  if (is.null(resDir) & !is.null(rds)) {
    if (class(rds) == "Seurat") {
      seuratObjFinal <<- rds
      print("RDS is provided with rds option")
    }
    else {
      rdsFname <- rds
      if (!file.exists(rdsFname)) 
        stop("Please execute getClusterMarker() to conduct integration analysis before running getClusterSummaryReplot().")
      seuratObjFinal <<- readRDS(file = as.character(rdsFname))
      print("Done for RDS read in")
    }
    resDir <- getwd()
  }
  else if (!is.null(resDir) & is.null(rds)) {
    rdsFname <- sprintf("%s/RDS_Dir/%s.rds", resDir, basename(resDir))
    resDir <- resDir
    if (!file.exists(rdsFname)) 
      stop("Please execute getClusterMarker() to conduct integration analysis before running getClusterSummaryReplot().")
    seuratObjFinal <<- readRDS(file = as.character(rdsFname))
    print("Done for RDS read in")
  }
  else if (is.null(resDir) & is.null(rds)) {
    stop("Error: please provide either option 'resDir' or 'rds', or both. ")
  }
  else if (!is.null(resDir) & !is.null(rds)) {
    if (class(rds) == "Seurat") {
      seuratObjFinal <<- rds
      print("RDS is provided with rds option")
    }
    else {
      rdsFname <- rds
      if (!file.exists(rdsFname)) 
        stop("Please execute getClusterMarker() to conduct integration analysis before running getClusterSummaryReplot().")
      seuratObjFinal <<- readRDS(file = as.character(rdsFname))
      print("Done for RDS read in")
    }
    resDir <- resDir
  }
  if (newAnnotation) {
    resDir <- paste(resDir, "scatterPlot_wNewAnnotation", 
                    sep = "/")
  }
  else {
    resDir <- paste(resDir, "scatterPlot_wOrgClusterAnnotation", 
                    sep = "/")
  }
  if (!dir.exists(resDir)) 
    dir.create(resDir)
  if (newAnnotation) {
    source(as.character(newAnnotationRscriptName))
  }
  if (!is.null(subsetOn)) {
    if (subsetOn == "idents") {
      cellcluster = subset
      clusterLevels <- levels(Seurat::Idents(seuratObjFinal))
      if (!is.null(cellcluster)) {
        if (any(!cellcluster %in% clusterLevels)) 
          stop("Please provide the corresponding cell clusters in identfied idents.")
        seuratObjFinal <- subset(seuratObjFinal, idents = as.character(cellcluster))
        print(sprintf("Note: Only these cell clusters (%s) can be explored", 
                      paste(cellcluster, collapse = ", ")))
      }
      else {
        print("Note: All cell clusters definided in idents can be explored.")
      }
    }
    else {
      sel.expConds = subset
      if (!subsetOn %in% colnames(seuratObjFinal@meta.data)) {
        stop("ERROR: 'subsetOn' does not exist in your 'rds' metadata.")
      }
      else {
        seuratObjFinal@meta.data$expCond <- seuratObjFinal@meta.data[, 
                                                                     grep(as.character(subsetOn), colnames(seuratObjFinal@meta.data))]
      }
      expCondLevels <- levels(factor(seuratObjFinal@meta.data$expCond))
      if (!is.null(sel.expConds)) {
        if (any(!sel.expConds %in% expCondLevels)) 
          stop("Please provide the corresponding experimental condition levesl specified in 'subset' option.")
        print(sprintf("Subsetting %s specific expCond: %s", 
                      length(sel.expConds), paste(sel.expConds, 
                                                  collapse = ",")))
        if (length(sel.expConds) == 1) {
          seuratObjFinal <- subset(seuratObjFinal, expCond == 
                                     sel.expConds)
        }
        else {
          for (i in 1:length(sel.expConds)) {
            seuratObjFinalPrep <- subset(seuratObjFinal, 
                                         subset = expCond == as.character(sel.expConds[i]))
            if (i == 1) {
              seuratObjFinalPrep2 = seuratObjFinalPrep
            }
            else {
              seuratObjFinalPrep2 <- merge(seuratObjFinalPrep2, 
                                           seuratObjFinalPrep)
            }
          }
          seuratObjFinal <- seuratObjFinalPrep2
        }
      }
    }
  }
  if (expCondCheck == "idents") {
    seuratObjFinal@meta.data$expCond <- Seurat::Idents(seuratObjFinal)
  }
  else {
    if (!expCondCheck %in% colnames(seuratObjFinal@meta.data)) {
      stop("ERROR: 'expCondCheck' does not exist in your 'rds' metadata.")
    }
    else {
      seuratObjFinal@meta.data$expCond <- seuratObjFinal@meta.data[, 
                                                                   grep(sprintf("^%s$", as.character(expCondCheck)), 
                                                                        colnames(seuratObjFinal@meta.data))]
    }
  }
  DefaultAssay(seuratObjFinal) <- "RNA"
  selected.groups <- as.character(unlist(selectedGroups))
  print("-=-=-=-=-=-=-=-=-")
  print("Start to extract expressions for scatter plot")
  res4plotPrep <- sapply(selected.groups, FUN = function(x) {
    subset.res <- subset(seuratObjFinal, expCond == x)
    norm.res <- subset.res@assays$RNA$data
    print(sprintf("%s cells for %s mean calculations", dim(norm.res)[2], 
                  x))
    norm.res.mean <- Matrix::rowMeans(x = norm.res)
    return(norm.res.mean)
  })
  if (!is.null(subsetGenes)) {
    res4plotPrep2 <- res4plotPrep[match(subsetGenes, rownames(res4plotPrep)), 
    ]
    res4plotPrep <- res4plotPrep2[!is.na(rownames(res4plotPrep2)), 
    ]
    print(sprintf("subsetGenes is on, %s gene expression are extracted for scatter plot.", 
                  dim(res4plotPrep)[1]))
  }
  else {
    print(sprintf("subsetGenes is off, %s gene expression are extracted for scatter plot.", 
                  dim(res4plotPrep)[1]))
  }
  print("END to extract expressions for scatter plot")
  print("-=-=-=-=-=-=-=-=-")
  res4plot <- sapply(selectedGroups, function(x) {
    if (length(x) == 1) {
      res <- res4plotPrep[, match(x, colnames(res4plotPrep))]
    }
    else {
      res <- rowMeans(res4plotPrep[, match(x, colnames(res4plotPrep))])
    }
    return(res)
  })
  colnames(res4plot) <- gsub(pattern = "-|[ ]|/", replacement = "_", 
                             colnames(res4plot))
  if (is.null(colnames(res4plot))) {
    colnames(res4plot) <- c("Group1", "Group2")
  }
  if (interactive) {
    easylabel::easylabel(res4plot, x = colnames(res4plot)[1], 
                         y = colnames(res4plot)[2], output_shiny = T)
  }
  else {
  }
  return(res4plot)
}


for (n in 1:length(comp_group_list)) {
  tryCatch({
    
    # get input
    mydata <- getScatterPlot(
      rds = seuratObj_annot, 
      expCondCheck = expCondCheck, 
      selectedGroups = comp_group_list[[n]],
      interactive = F)
    
    # plot scatters
    myplot <- ggscatter(
      as.data.frame(mydata), x='Group1', y='Group2',
      add = "reg.line",
      add.params = list(color = "blue", fill = "lightgray"),
      conf.int = TRUE, conf.int.level = 0.9) +
      # coord_equal(ratio = 1) +
      scale_x_continuous(breaks=seq(0,20,1)) +
      scale_y_continuous(breaks=seq(0,20,1)) +
      labs(title = paste0("scatter ", n), x = comp_group_list[[n]][1], y = comp_group_list[[n]][2]) +
      stat_regline_equation(label.y = 5.4) +
      stat_cor(method = "pearson", label.x.npc = "left", label.y = 5, r.accuracy = 0.01)
    # ggsave(paste0(resDir, "/scatterPlot_wOrgClusterAnnotation/", n, ".jpg"), plot = myplot, width = 4, height = 4, units = "in")
    
    # Extract R and p values
    r_value <- ggplot2::ggplot_build(myplot)$data[[4]]$r
    p_value <- ggplot2::ggplot_build(myplot)$data[[4]]$p
    
    # Append new row to the dataframe
    new_row <- data.frame(R = r_value, p = p_value)
    correlation_values_df <- rbind(correlation_values_df, new_row)
    
    print(paste0(comp_group_list[[n]], " r: ", r_value, " p: ", p_value))
  }, error = function(e) {
    message("Error occurred in comp_group ", n, ": ", e$message)
  })
  next
}

correlation_values_df <- cbind(as.data.frame(do.call(rbind, combinations)) , correlation_values_df$R, correlation_values_df$p)
colnames(correlation_values_df) <- c("x", "y", "r", "p")
write.table(correlation_values_df, paste0(workdir, "/correlation_values.txt"), sep = "\t", quote = FALSE, row.names = FALSE)


# 5.2 local visulization -----
workdir <- "integration_2/integration_2_leiden/cell_clustering/per_cell_and_majority_voting_curation"

corr_df <- read.table(paste0(workdir, "/correlation_values.txt"), sep = "\t", header = 1)


# corr_df$x <- sub(" r$", "", corr_df$x)  
# corr_df$y <- sub(" r$", "", corr_df$y) 
# write.table(corr_df, paste0(workdir, "/correlation_values.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

corr_df$r <- as.numeric(corr_df$r)
corr_df$p_value <- as.numeric(corr_df$p_value)



# Create the reverse version of corr_df
reversed_corr_df <- corr_df %>%
  rename(x = y, y = x) %>%
  mutate(r = r)  # No change to V3

# Combine the original and reversed data frames
combined_corr_df <- bind_rows(corr_df, reversed_corr_df)


# Convert the long-format data frame to a wide-format matrix
corr_matrix <- combined_corr_df[c(1,2,3)] %>%
  pivot_wider(names_from = y, values_from = r, values_fill = list(r = NA)) %>%
  column_to_rownames(var = "x") %>%
  as.matrix()


na_count <- sum(is.na(corr_matrix))
print(paste("Number of NA values:", na_count))

# Fill NA values with 1
corr_matrix[is.na(corr_matrix)] <- 1

cell_types <- c(
  "B",      
  "CD4 T",   
  "CD8 T",    
  "NK", 
  "Mono/Mph", 
  "EC",
  "DC"  )

stimuli.time <- c("unstim_0h", 
                          "LPS_4h", "LPS_18h",
                          "CD3_4h", "CD3_18h",
                          "Ig_4h", "Ig_18h")      
# Initialize an empty vector to store the ordered combinations
order_vector <- character()

# Use nested loops to create the order_vector
for (cell_type in cell_types) {
  for (stimuli_time in stimuli.time) {
    order_vector <- c(order_vector, paste0(cell_type, " ", stimuli_time))
  }
}

corr_matrix_ordered <- corr_matrix[order_vector, order_vector]

min_value <- min(corr_matrix_ordered, na.rm = TRUE)
max_value <- max(corr_matrix_ordered, na.rm = TRUE)

pdf(paste0(workdir, "/correlation.pdf"), width = 10.73, height = 6.59)
# Plotting the matrix
corrplot(corr_matrix_ordered, type = "lower", order = "original", 
         tl.col = "black", 
         tl.srt = 45, 
         tl.cex = 0.6,
         # tl.pos="n",
         # col = colorRampPalette(c(rep("#62a1db", 10),
         #                          "white", "#b01111"))(400),
         col = colorRampPalette(c(rep("#62a1db", 8),
                                  "white", "#b01111"))(100),
         method = "circle",
         col.lim = c(0.58, 1)
         )

dev.off()
# 
# all_levels <- unique(c(as.character(corr_df$V1), as.character(corr_df$V2)))
# 
# # Convert V1 and V2 to factors with the same levels
# corr_df$V1 <- factor(corr_df$V1, levels = all_levels)
# corr_df$V2 <- factor(corr_df$V2, levels = all_levels)
# 
# ggplot(corr_df, aes(y = V1, 
#                     x = V2,
#                     size = V3, 
#                     fill = V3)) +
#   geom_tile(color = "white") +
#   scale_fill_gradientn(colors = colorRampPalette(c('#f9dbbd', '#ffa5ab', "#da627d", "#a53860", "#450920"))(100)) +
#   theme_bw() +
#   geom_point(shape=21) +
#   labs(x = "Cell Type", y = "Cell Type", fill = "Correlation")+
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1),
#         legend.position = "left") 
# 
# # Plot using ggplot2
# ggplot(corr_df, aes(x = V1, y = V2, fill = V3)) +
#   geom_tile(color = "white") +
#   geom_point(aes(size = V3), shape = 21, color = "black") +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   coord_fixed() +
#   labs(x = "Cell Type", y = "Cell Type", fill = "Correlation")


#### 6.0 counts of cells ####


cell_types <- c(
  "Mono/Mph", # From "Monocyte-derived Mph" "Classical monocytes" and "Non-classical monocytes" "Alveolar Mph CCL3+" and "Alveolar Mph MT-positive"
  "B",        # From "B cells"
  "CD4 T",    # From "CD4 T cells"
  "CD8 T",    # From "CD8 T cells"
  "NK",       # From "NK cells"
  "DC",       # From "Migratory DCs"
  "EC"       # From "EC general capillary"
)

counts <- seuratObj_annot@meta.data %>% 
  group_by(cluster_annotation3) %>% 
  count() %>%
  filter(cluster_annotation3 %in% cell_types)
# selectedCol <- c("#f77f00", "#075a84", "#9163cb", "#ff4d6d", "#679436", "#f07167", "#603808", "#07beb8")
# 
# selectedCol <- c("#ff9b54", "#6d98ba", "#d4c2fc", "#ff758f", "#98c9a3", "#ffb4a2", "#ddb892", "#b2f7ef")
# 



names(table(Idents(seuratObj_annot)))


# selectedCol <- DiscretePalette(n = length(levels(Seurat::Idents(seuratObj_annot))), palette = 'parade')



selectedCol <- c("#075a84", "#9163cb", "#ff4d6d", "#679436", "#f77f00","#07beb8")

selectedCol <- c("#6d98ba", "#d4c2fc", "#ff758f", "#98c9a3", "#ff9b54", "#b2f7ef")

selectedCol<-c( 'cornflowerblue', 
                'mediumvioletred',
                'wheat3', 
                'lightsalmon',  
                "#7CE3D8", 
                'yellowgreen',
                'seagreen')

names(selectedCol) <- cell_types



pdf("integration_2/integration_2_leiden/cell_clustering/per_cell_and_majority_voting_curation/celltype_cnts.pdf", width = 5.9, height = 3)

ggplot(counts, aes(x = 2, y = n, fill = factor(cluster_annotation3, 
                                               levels=c( "B", "CD4 T", 
                                                         "CD8 T", "NK",
                                                         "Mono/Mph", "EC", "DC")))) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +    # Make it circular (polar coordinates)
  xlim(0.5, 2.5) +              # Create space for the hole in the middle
  theme_void() +                # Remove all background, gridlines, and axes
  # theme(legend.position = "none") +  # Position the legend
  labs(fill = "Cluster Annotation") +  # Label the legend
  scale_fill_manual(values = selectedCol) +
  guides(fill = guide_legend(title = "Cluster Annotation")) +
  geom_text(aes(label = ifelse(n > 3000, n, n)), position = position_stack(vjust = 0.5)) 

dev.off()

# #### 7.0 counts of cells by donor and treatment -----
# celltype_donor_stimuli.time_cnts <- rds@meta.data %>%
#   group_by(cluster_annotation, expCond.donor, expCond.stimuli.time, expCond.donor.stimuli.time) %>%
#   summarize(count = n(), .groups = 'drop') %>%
#   filter(cluster_annotation %in% c("B cells", "CD4 T cells", "CD8 T cells", "NK cells"))
# 
# colnames(celltype_donor_stimuli.time_cnts) <- gsub("^expCond\\.", "", colnames(celltype_donor_stimuli.time_cnts))
# 
# write.csv(celltype_donor_stimuli.time_cnts, file = "integration/integration_1/integration_1_louvain/cell_clustering/celltype_cnts.csv", row.names = FALSE)
# 





