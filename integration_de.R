rm(list=ls())
.libPaths(c("/ess/home/home1/ftian1/R/x86_64-pc-linux-gnu-library/4.2",
            "/gpfs/data/icelake-apps/software/gcc-12.1.0/R/4.2.1/lib64/R/library"))

library(Seurat)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(openxlsx) # read.xlsx
library(readxl)
library(writexl)
library(MAST)
library(xlsx)
library(tidyr)
library(data.table)



rds <- readRDS('integration_2/integration_2_leiden/RDS_Dir/integration_2_leiden_annot.rds')
names(table(Idents(rds)))

deResDir <- "integration_2/integration_2_leiden/"


#### 1.0 run on randi ####
## ---------------------------------------------------------------------------------------
getClusterExpCondDe <- function(resDir=NULL, rds=NULL, newAnnotation=F, newAnnotationRscriptName=NULL, expCondCheck='sample',
                                expCondCheckFname = NULL, cellcluster = NULL, compGroup, deMethod = 'wilcox',
                                covariateVarName = NULL, deseq2bulk.metaCovariateInput = NULL, covariateVarLevels = NULL, norm.method = 'UQ', run.dispersion = as.logical(T),
                                min.cells.group = 10, min.pct = 0.1, logfc.threshold = 0.25, pAdjValCutoff = 0.05, topNo = 10, debug = F, outputExcel = as.logical(T)) {
  options(java.parameters = "-Xmx128g")
  ## ---
  if (missing(compGroup)) stop("Please provide option 'compGroup' to specify which 2 groups in your updated experimental condition levels with 'expCondCheck' for comparision.")
  pAdjValCutoff           <- as.numeric(pAdjValCutoff)
  topNo                   <- as.numeric(topNo)
  deMethod                <- as.character(deMethod)
  newAnnotation           <- as.logical(newAnnotation)
  if (newAnnotation & is.null(newAnnotationRscriptName)) print("Option 'newAnnotation' is on, please provide corresponding option 'newAnnotationRscriptName'.")
  if (!is.null(covariateVarName)) {
    if (!deMethod %in% c('LR', 'negbinom', 'poisson','MAST', 'DESeq2')) stop("'covariateVarName' cannot be applied for the specified 'deMethod'.")
  }
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
  if (newAnnotation) {
    resDir <- paste(resDir, 'results_wNewAnnotation_DEGs', sep = '/')
  } else {
    resDir <- paste(resDir, 'results_wOrgClusterAnnotation_DEGs', sep = '/')
  }
  if (!dir.exists(resDir)) dir.create(resDir)
  print(sprintf('DEG identification on cell clusters will be saved at %s', resDir))
  ## -------------------------------------------------------------------------------------
  if (newAnnotation) {
    ## Assign cell type identity to clusters
    ## redefine the level of Idents on the y-axis can be adjusted here by inputting order for cell annotation
    source(as.character(newAnnotationRscriptName))
  }
  ##--------------------------------------------------------------------------------------##
  if (expCondCheck == 'sample') {
    if (is.null(expCondCheckFname)) {
      expCondCheckFname        <- 'expCond_sample'
    } else {
      expCondCheckFname        <- expCondCheckFname
    }
  } else {
    if (is.null(expCondCheckFname)) {
      expCondCheckFname        <- as.character(expCondCheck)
    } else {
      expCondCheckFname        <- expCondCheckFname
    }
  }
  ## update resDir
  resDir                  <- paste(resDir, expCondCheckFname, sep = '/')
  if (!dir.exists(resDir)) dir.create(resDir)
  ##--------------------------------------------------------------------------------------##
  ## update 'seuratObjFinal@meta.data$expCond' and add 'covariateVarName'
  if (expCondCheck == 'sample') {
    seuratObjFinal                     <- seuratObjFinal
    if (covariateVarName == 'sample') stop("DE tests are conducted on 'sample', no need to regress on this variable again with 'covariateVarName'.")
    if (deMethod == 'DESeq2.bulk') stop("DEseq2 bulk RNA tests are combining results on the sample level, no DE analysis can be conducted on the same level.")
  } else {
    if (deMethod == 'MAST') {
      ## covariateVar used for MAST latent.varible option
      if (!is.null(covariateVarName)) {
        if (covariateVarName == 'sample') {
          seuratObjFinal@meta.data$covariateVar = seuratObjFinal@meta.data$orig.ident
        } else {
          seuratObjFinal@meta.data$covariateVar = seuratObjFinal@meta.data[, grep(as.character(covariateVarName), colnames(seuratObjFinal@meta.data))]
        }
      }
      ## ---
    }
    if (deMethod == 'DESeq2.bulk') {
      seuratObjFinal@meta.data$deseq2bulk = seuratObjFinal@meta.data$orig.ident
    }
    if (!expCondCheck%in%colnames(seuratObjFinal@meta.data)) {
      stop("ERROR: 'expCondCheck' does not exist in your 'rds' metadata.")
    } else {
      seuratObjFinal@meta.data$expCond <- seuratObjFinal@meta.data[, grep(sprintf('^%s$', as.character(expCondCheck)), colnames(seuratObjFinal@meta.data))]
    }
  }
  ##--------------------------------------------------------------------------------------##
  Seurat::DefaultAssay(seuratObjFinal)   <- "RNA"
  ##--------------------------------------------------------------------------------------##
  ## if 'cellcluster' provided, only certain cell clusters to conduct DE analysis.
  clusterLevels <- levels(Seurat::Idents(seuratObjFinal))
  if (!is.null(cellcluster)) {
    if (any(!cellcluster %in% clusterLevels ) ) stop('Please provide the corresponding cell clusters in identfied idents.')
    print(sprintf("Start step 1: identifing DEGs for selected clusters (%s) based on experimental conditions in rds metadata table column '%s'.", paste(unlist(cellcluster), collapse = ', '), expCondCheck ))
    noClusters                             <- unique(cellcluster)
  } else {
    print(sprintf("Start step 1: identifing DEGs for all identified clusters based on experimental conditions in rds metadata table column '%s'.", expCondCheck))
    noClusters                             <- levels(Seurat::Idents(seuratObjFinal))
  }
  ## -----
  print(sprintf("A total of %s clusters to process", length(noClusters)))
  ## -
  clusterDeMarkers                       <- list()
  normCounts <- list()
  orgCounts  <- list()
  metaTabs   <- list()
  clusterDeResSummary                    <- matrix(data = NA, nrow = length(noClusters), ncol = 4)
  clusterTopDeMarkers                    <- list()
  clusterTopDeMarkersUp                  <- NA
  clusterTopDeMarkersUpCluster           <- NA
  clusterTopDeMarkersDown                <- NA
  clusterTopDeMarkersDownCluster         <- NA
  clusterDeMarkersFname                  <- paste(sprintf('%s/expCondCompDeMarkers_%s', resDir, gsub('/', '-', compGroup) ))
  ## -
  if (is.null(cellcluster)) {
    if (newAnnotation) {
      resFname1                  <- paste(clusterDeMarkersFname, '_full_wNewAnnotation_allClusters', sep = '')
      resFname2                  <- paste(clusterDeMarkersFname, '_adjSig_wNewAnnotation_allClusters', sep = '')
      resFname3                  <- paste(clusterDeMarkersFname, '_adjSig_up_wNewAnnotation_allClusters', sep = '')
      resFname4                  <- paste(clusterDeMarkersFname, '_adjSig_down_wNewAnnotation_allClusters', sep = '')
    } else {
      resFname1                  <- paste(clusterDeMarkersFname, '_full_allClusters', sep = '')
      resFname2                  <- paste(clusterDeMarkersFname, '_adjSig_allClusters', sep = '')
      resFname3                  <- paste(clusterDeMarkersFname, '_adjSig_up_allClusters', sep = '')
      resFname4                  <- paste(clusterDeMarkersFname, '_adjSig_down_allClusters', sep = '')
    }
  } else {
    if (newAnnotation) {
      resFname1                  <- paste(clusterDeMarkersFname, sprintf("_full_wNewAnnotation_%sSelClusters", length(noClusters)), sep = '')
      resFname2                  <- paste(clusterDeMarkersFname, sprintf("_adjSig_wNewAnnotation_%sSelClusters", length(noClusters)), sep = '')
      resFname3                  <- paste(clusterDeMarkersFname, sprintf("_adjSig_up_wNewAnnotation_%sSelClusters", length(noClusters)), sep = '')
      resFname4                  <- paste(clusterDeMarkersFname, sprintf("_adjSig_down_wNewAnnotation_%sSelClusters", length(noClusters)), sep = '')
    } else {
      resFname1                  <- paste(clusterDeMarkersFname, sprintf("_full_%sSelClusters", length(noClusters)), sep = '')
      resFname2                  <- paste(clusterDeMarkersFname, sprintf("_adjSig_%sSelClusters", length(noClusters)), sep = '')
      resFname3                  <- paste(clusterDeMarkersFname, sprintf("_adjSig_up_%sSelClusters", length(noClusters)), sep = '')
      resFname4                  <- paste(clusterDeMarkersFname, sprintf("_adjSig_down_%sSelClusters", length(noClusters)), sep = '')
    }
  }
  ## -
  if (debug) print(sprintf("clusters to be processed are %s", paste(noClusters, collapse = ', ')))
  for (i in 1:length(noClusters)) {
    # for (i in c(1,3,4)) {
    print(sprintf('Start 1.%s processing cluster %s for DE markers identification', i, noClusters[i]))
    seuratObjFinalSubet                  <- subset(seuratObjFinal, idents = as.character(noClusters[i]) )
    print("-=-=-=")
    print("experimental conditions:")
    print(table(seuratObjFinalSubet$expCond))
    print("-=-=-=")
    ## -------------------- ##
    if (deMethod == 'MAST' & !is.null(covariateVarName)) {
      if (!is.null(covariateVarLevels)) {
        if (debug) {
          print("*********")
          print(table(seuratObjFinalSubet@meta.data$covariateVar))
          print(sprintf("input 'covariateVarLevels': %s", paste(covariateVarLevels, collapse = ', ')))
          print("*********")
        }
        if (length(unique(seuratObjFinalSubet@meta.data$covariateVar))<length(covariateVarLevels)) {
          covariateVarLevels <- covariateVarLevels[which(!is.na(match(covariateVarLevels, unique(seuratObjFinalSubet@meta.data$covariateVar))))]
          if (debug) print(sprintf("covariateVarLevels is updated into %s. ", paste(covariateVarLevels, collapse = ', ')))
        }
        seuratObjFinalSubet@meta.data$covariateVar = factor(seuratObjFinalSubet@meta.data$covariateVar, levels = covariateVarLevels)
        # seuratObjFinalSubet@meta.data$covariateVar = relevel(x =  seuratObjFinalSubet@meta.data$covariateVar, ref = covariateVarLevels[length(covariateVarLevels)] )
      } else {
        covariateVarLevels <- names(sort(table(seuratObjFinalSubet@meta.data$covariateVar), decreasing = F))
        if (debug) {
          print("*********")
          print(table(seuratObjFinalSubet@meta.data$covariateVar))
          print(sprintf("sorted levels: %s", paste(covariateVarLevels, collapse = ', ')))
          print("*********")
        }
        seuratObjFinalSubet@meta.data$covariateVar = factor(seuratObjFinalSubet@meta.data$covariateVar, levels = covariateVarLevels)
        # seuratObjFinalSubet@meta.data$covariateVar = relevel(x =  seuratObjFinalSubet@meta.data$covariateVar, ref = covariateVarLevels[length(covariateVarLevels)] )
      }
    }
    ## -------------------- ##
    compGroupName                      <- strsplit(compGroup, split = '/')[[1]]
    if ( all(compGroupName %in% names(table(seuratObjFinalSubet$expCond))) ) {
      ## ---
      if ( length(compGroupName)>2 ) {
        stop("Too many comparision groups provided in 'compGroup', only 2 names shoulb be provided ")
      } else {
        if (length(compGroupName) == 1) {
          if (all(table(seuratObjFinalSubet$expCond) >= min.cells.group) ) {
            if (deMethod=='DESeq2.bulk') {
              deMarkersRes                 <- runBulkDEseq2(object = seuratObjFinalSubet, ident.1 =  compGroupName, min.pct = min.pct, metaCovariateInput = deseq2bulk.metaCovariateInput, debug = debug, min.cells.group = min.cells.group, norm.method = norm.method, run.dispersion = run.dispersion)
              deMarkers                    <- deMarkersRes$deres
              orgCount                     <- deMarkersRes$orgCount
              normCount                    <- deMarkersRes$normCount
              metaTab                      <- deMarkersRes$metaTab
            } else {
              if (is.null(covariateVarName)) {
                deMarkers                  <- FindMarkers(seuratObjFinalSubet, ident.1 = compGroupName, group.by = 'expCond', test.use = deMethod, min.cells.group = min.cells.group, logfc.threshold = logfc.threshold, min.pct =  min.pct)
              } else {
                if (debug) {
                  print('--------------')
                  print(table(seuratObjFinalSubet@meta.data$covariateVar))
                  print('--------------')
                }
                deMarkers                  <- FindMarkers(seuratObjFinalSubet, ident.1 = compGroupName, group.by = 'expCond', test.use = deMethod, min.cells.group = min.cells.group, logfc.threshold = logfc.threshold, min.pct =  min.pct, latent.vars = 'covariateVar')
              }
            }
          } else {
            print(sprintf('cluster %s has less than 3 cells in one experimental condiction, no DE markers can be identfied for this cluster', noClusters[i]))
          }
        } else if (length(compGroupName) == 2) {
          if (all(table(seuratObjFinalSubet$expCond)[match(compGroupName, names(table(seuratObjFinalSubet$expCond)))] >= min.cells.group) ) {
            if (deMethod=='DESeq2.bulk') {
              deMarkersRes                 <- runBulkDEseq2(object = seuratObjFinalSubet, ident.1 = compGroupName[1], ident.2 = compGroupName[2], min.pct = min.pct, metaCovariateInput = deseq2bulk.metaCovariateInput, debug = debug, min.cells.group = min.cells.group, norm.method = norm.method, run.dispersion = run.dispersion)
              deMarkers                    <- deMarkersRes$deres
              orgCount                     <- deMarkersRes$orgCount
              normCount                    <- deMarkersRes$normCount
              metaTab                      <- deMarkersRes$metaTab
            } else {
              if (is.null(covariateVarName)) {
                deMarkers                  <- FindMarkers(seuratObjFinalSubet, ident.1 = compGroupName[1], ident.2 = compGroupName[2], group.by = 'expCond', test.use = deMethod, min.cells.group = min.cells.group, logfc.threshold = logfc.threshold, min.pct =  min.pct)
              } else {
                if (debug) {
                  print('--------------')
                  print(table(seuratObjFinalSubet@meta.data$covariateVar))
                  print('--------------')
                }
                deMarkers                  <- FindMarkers(seuratObjFinalSubet, ident.1 = compGroupName[1], ident.2 = compGroupName[2], group.by = 'expCond', test.use = deMethod, min.cells.group = min.cells.group, logfc.threshold = logfc.threshold, min.pct =  min.pct, latent.vars = 'covariateVar')
              }
            }
          } else {
            print(sprintf('cluster %s has less than 3 cells in one experimental condiction, no DE markers can be identfied for this cluster', noClusters[i]))
          }
        }
        if(debug) print('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-')
        if(debug) print('111111111111111')
        clusterDeMarkers[[i]]          <- deMarkers
        if (deMethod!="DESeq2.bulk") {
          # print(sprintf('Maximum p_value is %s, Maximum adjusted p_value is %s', round(max(deMarkers$pvalue), digits = 4), round(max(deMarkers$padj, na.rm = T), digits = 4)))
          deMarkersAdjSig                <- deMarkers %>% dplyr::filter(p_val_adj <= pAdjValCutoff) %>% dplyr::filter(abs(avg_log2FC) >= logfc.threshold) %>% dplyr::mutate(perDiff = pct.1-pct.2) %>% dplyr::mutate(FC = ifelse(avg_log2FC>0, 2^avg_log2FC, -2^(-avg_log2FC)) )
          deMarkersAdjSigUp              <- deMarkersAdjSig %>% dplyr::filter(avg_log2FC > 0) %>% dplyr::arrange(desc(FC))
          deMarkersAdjSigDown            <- deMarkersAdjSig %>% dplyr::filter(avg_log2FC < 0) %>% dplyr::arrange(FC)
        } else {
          if(debug) print('22222222222222')
          normCounts[[i]] <- normCount
          orgCounts[[i]]  <- orgCount
          metaTabs[[i]]   <- metaTab
          # print(sprintf('Maximum p_value is %s, Maximum adjusted p_value is %s', round(max(deMarkers$p_val), digits = 4), round(max(deMarkers$p_val_adj), digits = 4)))
          deMarkersAdjSig                <- deMarkers %>% dplyr::filter(padj <= pAdjValCutoff) %>% dplyr::filter(abs(log2FoldChange) > logfc.threshold) %>% dplyr::mutate(FC = ifelse(log2FoldChange>0, 2^log2FoldChange, -2^(-log2FoldChange)) )
          deMarkersAdjSigUp              <- deMarkersAdjSig %>% dplyr::filter(log2FoldChange > 0) %>% dplyr::arrange(desc(FC))
          deMarkersAdjSigDown            <- deMarkersAdjSig %>% dplyr::filter(log2FoldChange < 0) %>% dplyr::arrange(FC)
        }
        if(debug) print('333333333333')
        if (dim(deMarkersAdjSigUp)[1] > topNo) {
          topNo1 = topNo
        } else {
          if (dim(deMarkersAdjSigUp)[1]!=0) {
            topNo1 = dim(deMarkersAdjSigUp)[1]
          } else {
            topNo1 = 0
          }
        }
        if(debug) print('44444444444444')
        if (dim(deMarkersAdjSigDown)[1] > topNo) {
          topNo2 = topNo
        } else {
          if (dim(deMarkersAdjSigDown)[1]!=0) {
            topNo2 = dim(deMarkersAdjSigDown)[1]
          } else {
            topNo2 = 0
          }
        }
        if(debug) print('55555555555555')
        if (topNo1!=0 & topNo2!=0) {
          clusterTopDeMarkers[[i]]       <- list('up' = rownames(deMarkersAdjSigUp)[1:topNo1],
                                                 'down' = rownames(deMarkersAdjSigDown)[1:topNo2])
        } else if (topNo1!=0 & topNo2==0) {
          clusterTopDeMarkers[[i]]       <- list('up' = rownames(deMarkersAdjSigUp)[1:topNo1])
        } else if (topNo1==0 & topNo2!=0) {
          clusterTopDeMarkers[[i]]       <- list('down' = rownames(deMarkersAdjSigDown)[1:topNo2])
        } else {
          clusterTopDeMarkers[[i]]       <- NA
        }
        
        if(debug) print('666666666666666')
        if (dim(deMarkersAdjSigUp)[1]!=0) {
          clusterTopDeMarkersUp          <- c(clusterTopDeMarkersUp, rownames(deMarkersAdjSigUp)[1:topNo1])
          clusterTopDeMarkersUpCluster   <- c(clusterTopDeMarkersUpCluster, rep(noClusters[i], length(rownames(deMarkersAdjSigUp)[1:topNo1])) )
        }
        if (dim(deMarkersAdjSigDown)[1]!=0){
          clusterTopDeMarkersDown        <- c(clusterTopDeMarkersDown, rownames(deMarkersAdjSigDown)[1:topNo2])
          clusterTopDeMarkersDownCluster <- c(clusterTopDeMarkersDownCluster, rep(noClusters[i], length(rownames(deMarkersAdjSigDown)[1:topNo2])))
        }
        clusterDeResSummary[i,]        <- c(dim(deMarkers)[1], dim(deMarkersAdjSig)[1], dim(deMarkersAdjSigUp)[1], dim(deMarkersAdjSigDown)[1] )
        print(sprintf('out of %s DE markers, %s are significantly DE at adjusted p-value of %s', dim(deMarkers)[1], dim(deMarkersAdjSig)[1], pAdjValCutoff))
        print(sprintf('out of %s significantly DE markers, %s are positively expressed, %s are negative expressed.', dim(deMarkersAdjSig)[1], dim(deMarkersAdjSigUp)[1], dim(deMarkersAdjSigDown)[1] ))
        ## -
        if(debug) print('777777777777777')
        if (outputExcel) {
          if (i ==1) {
            if (dim(deMarkers)[1] > 0) write.xlsx(x = deMarkers, file = sprintf('%s.xlsx', resFname1), sheetName = paste('cluster', gsub('/|:|[ ]|-', '', noClusters[i]), sep = '_'), row.names = T, append = F )
            if (dim(deMarkersAdjSig)[1] > 0) write.xlsx(x = deMarkersAdjSig, file = sprintf('%s.xlsx', resFname2), sheetName = paste('cluster', gsub('/|:|[ ]|-', '', noClusters[i]), sep = '_'), row.names = T, append = F )
            if (dim(deMarkersAdjSigUp)[1] > 0) write.xlsx(x = deMarkersAdjSigUp, file = sprintf('%s.xlsx', resFname3), sheetName = paste('cluster', gsub('/|:|[ ]|-', '', noClusters[i]), sep = '_'), row.names = T, append = F )
            if (dim(deMarkersAdjSigDown)[1] > 0) write.xlsx(x = deMarkersAdjSigDown, file = sprintf('%s.xlsx', resFname4), sheetName = paste('cluster', gsub('/|:|[ ]|-', '', noClusters[i]), sep = '_'), row.names = T, append = F )
          } else {
            if (dim(deMarkers)[1] > 0) write.xlsx(x = deMarkers, file = sprintf('%s.xlsx', resFname1), sheetName = paste('cluster', gsub('/|:|[ ]|-', '', noClusters[i]), sep = '_'), row.names = T, append = T )
            if (dim(deMarkersAdjSig)[1] > 0)write.xlsx(x = deMarkersAdjSig, file = sprintf('%s.xlsx', resFname2), sheetName = paste('cluster', gsub('/|:|[ ]|-', '', noClusters[i]), sep = '_'), row.names = T, append = T )
            if (dim(deMarkersAdjSigUp)[1] > 0) write.xlsx(x = deMarkersAdjSigUp, file = sprintf('%s.xlsx', resFname3), sheetName = paste('cluster', gsub('/|:|[ ]|-', '', noClusters[i]), sep = '_'), row.names = T, append = T )
            if (dim(deMarkersAdjSigDown)[1] > 0) write.xlsx(x = deMarkersAdjSigDown, file = sprintf('%s.xlsx', resFname4), sheetName = paste('cluster', gsub('/|:|[ ]|-', '', noClusters[i]), sep = '_'), row.names = T, append = T )
          }
        } else {
          if (dim(deMarkers)[1] > 0) write.table(x = deMarkers, file = sprintf('%s_cluster%s.txt', resFname1, gsub('/|:|[ ]|-', '', noClusters[i]) ), row.names = T, quote = F, col.names = NA, sep = '\t')
          if (dim(deMarkersAdjSig)[1] > 0) write.table(x = deMarkersAdjSig, file = sprintf('%s_cluster%s.txt', resFname2, gsub('/|:|[ ]|-', '', noClusters[i])), row.names = T, quote = F, col.names = NA, sep = '\t')
          if (dim(deMarkersAdjSigUp)[1] > 0) write.table(x = deMarkersAdjSigUp, file = sprintf('%s_cluster%s.txt', resFname3, gsub('/|:|[ ]|-', '', noClusters[i])), row.names = T, quote = F, col.names = NA, sep = '\t')
          if (dim(deMarkersAdjSigDown)[1] > 0) write.table(x = deMarkersAdjSigDown, file = sprintf('%s_cluster%s.txt', resFname4, gsub('/|:|[ ]|-', '', noClusters[i])), row.names = T, quote = F, col.names = NA, sep = '\t')
        }
        if(debug) print('8888888888888888')
      }
      ## ---
    } else {
      stop("Provided comparision group names in 'compGroup' does not match experimental conditions ")
    }
    print(sprintf('End 1.%s processing cluster %s for DE markers identification', i, noClusters[i]))
    print('=========')
  }
  ## -
  clusterTopDeMarkersUpComb          <- data.frame('Gene' = clusterTopDeMarkersUp[-1], 'geneType'= clusterTopDeMarkersUpCluster[-1])
  clusterTopDeMarkersDownComb        <- data.frame('Gene' = clusterTopDeMarkersDown[-1], 'geneType' = clusterTopDeMarkersDownCluster[-1])
  if(debug) print(head(clusterTopDeMarkersUpComb))
  if (dim(clusterTopDeMarkersUpComb[!is.na(clusterTopDeMarkersUpComb$Gene),])[1]!=0) {
    write.xlsx( x = clusterTopDeMarkersUpComb[!is.na(clusterTopDeMarkersUpComb$Gene),], file = sprintf('%s_top%s_upDe.xlsx', clusterDeMarkersFname, topNo1), sheetName = 'up', row.names = F, col.names = T, append = F)
  }
  if (dim(clusterTopDeMarkersDownComb[!is.na(clusterTopDeMarkersDownComb$Gene),])[1] !=0) {
    write.xlsx( x = clusterTopDeMarkersDownComb[!is.na(clusterTopDeMarkersDownComb$Gene),], file = sprintf('%s_top%s_downDe.xlsx', clusterDeMarkersFname, topNo2), sheetName = 'down', row.names = F, col.names = T, append = F)
  }
  ## -
  rownames(clusterDeResSummary)      <- noClusters
  colnames(clusterDeResSummary)      <- c('Total', 'DE', 'Up', 'Down')
  clusterDeResSummaryFname           <- paste(clusterDeMarkersFname, '_NoDEmarkers_summary.xlsx', sep = '')
  write.xlsx( x = clusterDeResSummary, file = clusterDeResSummaryFname, sheetName = 'No. summary', row.names = T, col.names = T, append = F)
  ## ---
  names(clusterDeMarkers)            <- noClusters
  names(clusterTopDeMarkers)         <- noClusters
  if (deMethod=='DESeq2.bulk') {
    names(normCounts)                  <- noClusters
    names(orgCounts)                   <- noClusters
    names(metaTabs)                    <- noClusters
  }
  print(sprintf('End step 1: identifing DEGs for each identified clusters based on sample name experimental conditions'))
  print('-=-=-=-=-=-=-=-=-=-=-=-=-=-=')
  if (deMethod=='DESeq2.bulk') {
    return(list('orgCounts' = orgCounts, 'normCounts' = normCounts, 'metaTabs' = metaTabs, 'clusterDeMarkers' = clusterDeMarkers, 'clusterDeResSummary' = clusterDeResSummary, 'clusterTopDeMarkers' = clusterTopDeMarkers ))
  }else {
    return(list('clusterDeMarkers' = clusterDeMarkers, 'clusterDeResSummary' = clusterDeResSummary, 'clusterTopDeMarkers' = clusterTopDeMarkers ))
  }
  ## -------------------------------------------------------------------------------------
}

## ----------------- ##
IdentsToCells <- function( object, group.by, ident.1, ident.2 = NULL, cellnames.use) {
  if (is.null(x = ident.1)) {
    stop("Please provide ident.1")
  }
  Seurat::Idents(object = object) <- group.by
  ident.1 <- Seurat::WhichCells(object = object, idents = ident.1)
  if (is.null(ident.2)) {
    ident.2 <- setdiff(x = cellnames.use, y = ident.1)
  } else {
    ident.2 <- Seurat::WhichCells(object = object, idents = ident.2)
  }
  return(list(group1 = ident.1, group2 = ident.2))
}
## ----------------- ##
runBulkDEseq2 <- function(object, ident.1, ident.2=NULL, min.pct, metaCovariateInput = NULL, debug = F, min.cells.group, norm.method = 'TMM', run.dispersion = as.logical(T)) {
  cellGroup = 'expCond' ##updated with 'expCondCheck', compGroupName should be levels shown in 'expCond'
  min.cells.group = as.numeric(min.cells.group)
  if (debug) print(table(object@meta.data[,match(cellGroup, colnames(object@meta.data))]))
  if (is.null(ident.2)) {
    print(sprintf("group1 = %s, group2 = NULL", ident.1))
  } else {
    print(sprintf("group1 = %s, group2 = %s", ident.1, ident.2))
  }
  ## select cells corresponding ident.1/2
  cells      <- IdentsToCells(object = object, group.by = cellGroup, ident.1 = ident.1, ident.2 = ident.2,  cellnames.use = colnames(object) )
  if (debug) print(str(cells))
  # ## calculate normalization factors with respect to all selected cells, either full sets or subsets.
  # cells.meta <- object@meta.data[match(unlist(cells), rownames(object@meta.data)),]
  # cells.no   <- count(cells.meta, deseq2bulk) %>% as.data.frame()
  # cells.no$norm.factor<- cells.no$n/median(cells.no$n)
  ## calculate FC results
  fc.results <- Seurat::FoldChange(object = object, ident.1 = ident.1, ident.2 = ident.2,  group.by = cellGroup )
  if (debug) {
    print(sprintf("Initially, %s genes are used for FC calculation.", dim(fc.results)[1]))
  }
  ## filter based on min.pct option
  fc.results$pmax <- pmax(fc.results$pct.1, fc.results$pct.2)
  fc.results      <- fc.results[fc.results$pmax>min.pct,]
  genes4de        <- rownames(fc.results)
  print(sprintf("After low expression removal with min.pct = %s, %s genes will be used for DE analysis.", min.pct, length(genes4de)))
  ## aggregating counts based on deseq2bulk inherited from 'orig.ident', and meanwhile prepare corresponding metatab.
  countInputs     <- lapply(1:length(cells), function(x) {
    counts             <- Seurat::FetchData(object = object, vars = genes4de, cells = cells[[x]], slot = "count")
    counts$deseq2bulk  <- Seurat::FetchData(object = object, vars = 'deseq2bulk', cells = cells[[x]])$deseq2bulk
    counts.group       <- counts %>% group_by(deseq2bulk) %>% summarise(across(everything(), list(sum))) %>% as.data.frame()
    colnames(counts.group) <- gsub('_1', '', colnames(counts.group))
    ## normalize the sum by median of number of cells weight
    deseq2bulk.no      <- as.data.frame(table(counts$deseq2bulk))
    ## remove samples with low expression cell number at cut-off of 'min.cells.group'
    if (debug) print("before remval")
    if (debug)  print(deseq2bulk.no)
    deseq2bulk.no        <- deseq2bulk.no[deseq2bulk.no$Freq >= min.cells.group,]
    if (debug) print("after remval")
    if (debug) print(deseq2bulk.no)
    if (debug) print('-=-=-=-=-=-=-=-=-')
    deseq2bulk.no$weight <- deseq2bulk.no$Freq/median(deseq2bulk.no$Freq)
    ## match count table with weight table
    counts.group.match   <-counts.group %>% dplyr::filter(deseq2bulk %in% as.character(deseq2bulk.no$Var1) )
    ##match deseq2bulk.no has the same sample order as shown in counts.group.match
    deseq2bulk.no        <- deseq2bulk.no[match(counts.group.match$deseq2bulk, deseq2bulk.no$Var1),]
    if (debug) print(sprintf("min.cells.group = %s", min.cells.group))
    if (debug) print(sprintf("Before: [row=%s %s], After [row=%s %s]", dim(counts.group)[1], dim(counts.group)[2], dim(counts.group.match)[1], dim(counts.group.match)[2] ) )
    counts.group2            <- round(counts.group.match[,-1]/deseq2bulk.no$weight, digits = 0)
    counts.group2$deseq2bulk <- counts.group.match$deseq2bulk
    counts.group             <- counts.group2
    if (debug) print(dim(counts.group))
    if (debug) print(deseq2bulk.no)
    ## -
    counts.group$group <- paste('group', x, sep = '')
    ## match on metaTabPrep
    metaTabPrep               <- data.frame('bulksamp' = deseq2bulk.no %>% dplyr::pull(Var1))
    if (!is.null(metaCovariateInput)) {
      colnames(metaCovariateInput) <- tolower(colnames(metaCovariateInput))
      if (!any(grepl('^bulksamp$', colnames(metaCovariateInput)))) stop("Please provide corresponding 'metaCovariateInput' with column 'bulksamp' to represent bulk pseudo samples. ")
      metaTabPrep <- dplyr::left_join(x = metaTabPrep, y = metaCovariateInput, by = 'bulksamp')
    }
    metaTabPrep$group <- paste('group', x, sep = '')
    ## ----
    return(list('countTab' = counts.group, 'metaTab' = metaTabPrep))
  } )
  ## combining 2 list of cells counts into a combined DF as 'countTab'
  countInputsComb <- rbind(countInputs[[1]]$countTab, countInputs[[2]]$countTab)
  ## ----------------- ##
  if (debug) print("Complete bulk aggregation.")
  if (debug) print(countInputsComb[,1:4])
  countTab <- countInputsComb %>% dplyr::select(-group) %>% dplyr::select(-deseq2bulk) %>% t()
  colnames(countTab) <- countInputsComb %>% dplyr::pull(deseq2bulk)
  if (debug) print("header of countTab")
  if (debug) print(head(countTab))
  if (debug) print(dim(countTab))
  # if (debug) print(sprintf("countTab.txt in '%s'.", paste(getwd(), 'countTab.txt', sep = '/')))
  # if (debug) write.table(x = countTab, file = paste(getwd(), 'countTab.txt', sep = '/'), col.names = NA, row.names = T, quote = F, sep = '\t')
  if (debug) print('-=-=-=-=-=-=-=-=-=-=-')
  ## ----------------- ##
  ## prepare metaTab, if 'metaCovariateInput' != NULL, add corresponding covariates for model DESeq2 model fitting.
  metaTab         <- rbind(countInputs[[1]]$metaTab, countInputs[[2]]$metaTab)
  metaTab2        <- metaTab %>% dplyr::select(-bulksamp) ## remove sample name for model formular establishment
  ## metaTab2 and metaTab are identical, except 'bulksamp' column is removed in metaTab2 for DESeq2 model formular and DESeq2 readin.
  # metaTab2      <- as.data.frame(unclass(metaTab2),stringsAsFactors=TRUE) ##change all character column into factor
  # rownames(metaTab2) <- as.character(metaTab$bulksamp)
  ## ----------------- ##
  if (debug) print("Full metaTab")
  if (debug) print(metaTab)
  print('-=-=-=-=-=-=-=-=-=-=-')
  print(sprintf("A total of %s samples", dim(metaTab)[1]))
  print(table(metaTab$group))
  print('-=-=-=-=-=-=-=-=-=-=-')
  ## ----------------- ##
  design.formula <- formula(paste0(' ~ ', paste(colnames(metaTab2), collapse = '+')))
  print("Design formular is:")
  print(design.formula)
  print('-=-=-=-=-=-=-=-=-=-=-')
  #
  if (norm.method == 'TMM') {
    dge   <- edgeR::DGEList(counts = countTab)
    if (debug) print("using TMM normalization")
    dge           <- edgeR::calcNormFactors(object = dge, method = 'TMM')
    normCount     <- edgeR::cpm(dge)
    dds <- DESeq2::DESeqDataSetFromMatrix(countData=round(edgeR::cpm(dge), digits = 0), colData=S4Vectors::DataFrame(metaTab2), design = design.formula )
  } else if (norm.method == 'UQ'){
    dge   <- edgeR::DGEList(counts = countTab)
    if (debug) print("using upperquartile normalization")
    dge   <- edgeR::calcNormFactors(object = dge, method = 'upperquartile')
    normCount     <- edgeR::cpm(dge)
    dds <- DESeq2::DESeqDataSetFromMatrix(countData=round(edgeR::cpm(dge), digits = 0), colData=S4Vectors::DataFrame(metaTab2), design = design.formula )
  } else {
    normCount     <- edgeR::cpm(countTab)
    dds   <- DESeq2::DESeqDataSetFromMatrix(countData=countTab, colData=S4Vectors::DataFrame(metaTab2), design = design.formula )
  }
  if (norm.method == 'TMM' | norm.method == 'UQ') {
    BiocGenerics::sizeFactors(dds) <- rep(1, length(colnames(edgeR::cpm(dge))))
    if (run.dispersion) {
      dds <- DESeq2::estimateDispersions(dds)
    } else {
      DESeq2::dispersions(dds) <- rep(1, length(rownames(edgeR::cpm(dge))))
    }
    if (debug) print("size factor is")
    if (debug) print(BiocGenerics::sizeFactors(dds))
    if (debug) print("dispersion head is")
    if (debug) print(head(DESeq2::dispersions(dds)))
    dds <- DESeq2::nbinomWaldTest(dds)
  } else {
    dds <- DESeq2::estimateSizeFactors(dds)
    dds <- DESeq2::estimateDispersions(dds)
    dds <- DESeq2::nbinomWaldTest(dds)
    # dds           <- DESeq2::DESeq(dds)
  }
  ## ---
  if (debug) print("resultsNames is")
  if (debug) print(DESeq2::resultsNames(dds))
  ## ---
  res           <- DESeq2::results(dds, contrast=c("group", 'group1', 'group2' ), format="DataFrame")
  res           <- as.data.frame(res)
  tpdeseq2.save <- res[order(res$padj) , ]
  if (debug) print(head(tpdeseq2.save))
  # print('090909090909')
  # print(head(countTab))
  # print('-=-=-=-')
  # print(head(normCount))
  # print('-=-=-=-')
  # print('090909090909')
  return(list(deres = tpdeseq2.save, orgCount = countTab, normCount =  normCount, metaTab = metaTab))
  print("DESeq2 bulk analysis is completed.")
}
################# DEGs (MAST) across major cell types, do "CD3_18h", "CD3_4h", "Ig_18h", "Ig_4h", "LPS_18h", "LPS_4h" vs "unstim_0h"################# 
compGroup1 <- names(table(rds$expCond.stimuli.time))[1:6] 
compGroup2 <- names(table(rds$expCond.stimuli.time))[7]
compGroups <- paste0(compGroup1, "/", compGroup2) # treatment_tp/unstim_0h

for (compGroup in compGroups) {
  # sink("getClusterExpCondDe_MAST.log") #record running log for the codes below
  getClusterExpCondDe(resDir = deResDir,
                      rds = rds,
                      newAnnotation = F,
                      expCondCheck = 'expCond.stimuli.time',
                      expCondCheckFname = paste0('stimuli_time', '_MAST'),
                      compGroup = compGroup,
                      cellcluster = c("B", "CD4 T", "CD8 T", "NK",
                                      "Mono/Mph", "DC", "EC"),
                      deMethod = 'MAST',
                      min.cells.group = 6,
                      topNo = 20
  )
  # sink() #stop recording log
}


################# DEGs (MAST) across major cell types, do "CD3", "Ig", "LPS" vs "unstim_0h"################# 
compGroup1 <- names(table(rds$expCond.stimuli))[1:3] 
compGroup2 <- names(table(rds$expCond.stimuli))[4]
compGroups <- paste0(compGroup1, "/", compGroup2) # treatment_tp/unstim_0h

for (compGroup in compGroups) {
  # sink("getClusterExpCondDe_MAST.log") #record running log for the codes below
  getClusterExpCondDe(resDir = deResDir,
                      rds = rds,
                      newAnnotation = F,
                      expCondCheck = 'expCond.stimuli',
                      expCondCheckFname = paste0('stimuli', '_MAST'),
                      compGroup = compGroup,
                      cellcluster = c("B", "CD4 T", "CD8 T", "NK",
                                      "Mono/Mph","DC", "EC"),
                      deMethod = 'MAST',
                      min.cells.group = 6,
                      topNo = 20
  )
  # sink() #stop recording log
}

################# DEGs (MAST) for B HLA-DQA2|DQB2 against B ################# 
rds.ig.unstim <- subset(rds, subset = (expCond.stimuli == "Ig" | expCond.stimuli == "unstim"))

rds.ig.unstim@meta.data$HLA <- 'None'

grep('DQA2', rownames(rds$RNA$data), value = FALSE)

grep('DQB2', rownames(rds$RNA$data), value = FALSE)

rds.ig.unstim@meta.data$HLA[which(rds.ig.unstim$RNA$data[rownames(rds.ig.unstim$RNA$data)[6458], ] != 0)] <- 'HLA-DQA2|DQB2'
rds.ig.unstim@meta.data$HLA[which(rds.ig.unstim$RNA$data[rownames(rds.ig.unstim$RNA$data)[6459], ] != 0)] <- 'HLA-DQA2|DQB2'

rds.ig.unstim@meta.data$expCond.stimuli.time.HLA <- paste(rds.ig.unstim@meta.data$expCond.stimuli.time, 
                                                          rds.ig.unstim@meta.data$HLA, sep = "_")

rds.ig.unstim$expCond.stimuli.time.HLA <- gsub("unstim_.*", "unstim", rds.ig.unstim$expCond.stimuli.time.HLA)

rds.ig.unstim@meta.data$expCond.stimuli.time.HLA.donor <- paste(rds.ig.unstim@meta.data$expCond.stimuli.time.HLA, 
                                                                rds.ig.unstim@meta.data$expCond.donor, sep = "_")

# rds.ig.unstim@meta.data %>% count(expCond.stimuli.time.HLA)

# names(table(Idents(rds))) # cell types

# only consider ig and unstim -> compare DQA2|DQB2 vs None
names(table(rds.ig.unstim$HLA))


getClusterExpCondDe(resDir = deResDir,
                    rds = rds.ig.unstim,
                    newAnnotation = F,
                    expCondCheck = 'HLA',
                    expCondCheckFname = 'HLA-DQA2|DQB2_ig_unstim_MAST',
                    cellcluster = "B",
                    compGroup = "HLA-DQA2|DQB2/None",
                    deMethod = 'MAST',
                    min.cells.group = 6,
                    topNo = 20
)

################# DEGs (MAST) for B HLA-DQA2&DQB2 against B ################# 
rds.ig.unstim <- subset(rds, subset = (expCond.stimuli == "Ig" | expCond.stimuli == "unstim"))

rds.ig.unstim@meta.data$HLA <- 'None'

grep('DQA2', rownames(rds$RNA$data), value = FALSE)

grep('DQB2', rownames(rds$RNA$data), value = FALSE)

rds.ig.unstim@meta.data$HLA[
  which(
    rds.ig.unstim$RNA$data[rownames(rds.ig.unstim$RNA$data)[6458], ] != 0 & 
      rds.ig.unstim$RNA$data[rownames(rds.ig.unstim$RNA$data)[6459], ] != 0
  )
] <- 'HLA-DQA2&DQB2'

rds.ig.unstim@meta.data$expCond.stimuli.time.HLA <- paste(rds.ig.unstim@meta.data$expCond.stimuli.time, 
                                                          rds.ig.unstim@meta.data$HLA, sep = "_")

rds.ig.unstim$expCond.stimuli.time.HLA <- gsub("unstim_.*", "unstim", rds.ig.unstim$expCond.stimuli.time.HLA)

rds.ig.unstim@meta.data$expCond.stimuli.time.HLA.donor <- paste(rds.ig.unstim@meta.data$expCond.stimuli.time.HLA, 
                                                                rds.ig.unstim@meta.data$expCond.donor, sep = "_")

# rds.ig.unstim@meta.data %>% count(expCond.stimuli.time.HLA)

# names(table(Idents(rds))) # cell types

# only consider ig and unstim -> compare DQA2|DQB2 vs None
names(table(rds.ig.unstim$HLA))


getClusterExpCondDe(resDir = deResDir,
                    rds = rds.ig.unstim,
                    newAnnotation = F,
                    expCondCheck = 'HLA',
                    expCondCheckFname = 'HLA-DQA2&DQB2_ig_unstim_MAST',
                    cellcluster = "B",
                    compGroup = "HLA-DQA2&DQB2/None",
                    deMethod = 'MAST',
                    min.cells.group = 6,
                    topNo = 20
)



################# DEGs (MAST) for B HLA-DQA2 against B ################# 
rds.ig.unstim <- subset(rds, subset = (expCond.stimuli == "Ig" | expCond.stimuli == "unstim"))

rds.ig.unstim@meta.data$HLA <- 'None'

grep('DQA2', rownames(rds$RNA$data), value = FALSE)

grep('DQB2', rownames(rds$RNA$data), value = FALSE)

rds.ig.unstim@meta.data$HLA[which(rds.ig.unstim$RNA$data[rownames(rds.ig.unstim$RNA$data)[6458], ] != 0)] <- 'HLA-DQA2'

rds.ig.unstim@meta.data$expCond.stimuli.time.HLA <- paste(rds.ig.unstim@meta.data$expCond.stimuli.time, 
                                                          rds.ig.unstim@meta.data$HLA, sep = "_")

rds.ig.unstim$expCond.stimuli.time.HLA <- gsub("unstim_.*", "unstim", rds.ig.unstim$expCond.stimuli.time.HLA)

rds.ig.unstim@meta.data$expCond.stimuli.time.HLA.donor <- paste(rds.ig.unstim@meta.data$expCond.stimuli.time.HLA, 
                                                                rds.ig.unstim@meta.data$expCond.donor, sep = "_")

# rds.ig.unstim@meta.data %>% count(expCond.stimuli.time.HLA)

# names(table(Idents(rds))) # cell types

# only consider ig and unstim -> compare DQA2|DQB2 vs None
names(table(rds.ig.unstim$HLA))


getClusterExpCondDe(resDir = deResDir,
                    rds = rds.ig.unstim,
                    newAnnotation = F,
                    expCondCheck = 'HLA',
                    expCondCheckFname = 'HLA-DQA2_ig_unstim_MAST',
                    cellcluster = "B",
                    compGroup = "HLA-DQA2/None",
                    deMethod = 'MAST',
                    min.cells.group = 6,
                    topNo = 20
)



################# DEGs (MAST) for B HLA-DQB2 against B ################# 
rds.ig.unstim <- subset(rds, subset = (expCond.stimuli == "Ig" | expCond.stimuli == "unstim"))

rds.ig.unstim@meta.data$HLA <- 'None'

grep('DQA2', rownames(rds$RNA$data), value = FALSE)

grep('DQB2', rownames(rds$RNA$data), value = FALSE)


rds.ig.unstim@meta.data$HLA[which(rds.ig.unstim$RNA$data[rownames(rds.ig.unstim$RNA$data)[6459], ] != 0)] <- 'HLA-DQB2'

rds.ig.unstim@meta.data$expCond.stimuli.time.HLA <- paste(rds.ig.unstim@meta.data$expCond.stimuli.time, 
                                                          rds.ig.unstim@meta.data$HLA, sep = "_")

rds.ig.unstim$expCond.stimuli.time.HLA <- gsub("unstim_.*", "unstim", rds.ig.unstim$expCond.stimuli.time.HLA)

rds.ig.unstim@meta.data$expCond.stimuli.time.HLA.donor <- paste(rds.ig.unstim@meta.data$expCond.stimuli.time.HLA, 
                                                                rds.ig.unstim@meta.data$expCond.donor, sep = "_")

# rds.ig.unstim@meta.data %>% count(expCond.stimuli.time.HLA)

# names(table(Idents(rds))) # cell types

# only consider ig and unstim -> compare DQA2|DQB2 vs None
names(table(rds.ig.unstim$HLA))


getClusterExpCondDe(resDir = deResDir,
                    rds = rds.ig.unstim,
                    newAnnotation = F,
                    expCondCheck = 'HLA',
                    expCondCheckFname = 'HLA-DQB2_ig_unstim_MAST',
                    cellcluster = "B",
                    compGroup = "HLA-DQB2/None",
                    deMethod = 'MAST',
                    min.cells.group = 6,
                    topNo = 20
)





################# DEGs (MAST) for B asthma against B non-asthma ################# 
rds.ig.unstim <- subset(rds, subset = (expCond.stimuli == "Ig" | expCond.stimuli == "unstim"))

# rds.ig.unstim@meta.data %>% count(expCond.stimuli.time.HLA)

# names(table(Idents(rds))) # cell types

# only consider ig and unstim -> compare Asthmatic vs Non-asthmatic
names(table(rds.ig.unstim$expCond.asthma))


getClusterExpCondDe(resDir = deResDir,
                    rds = rds.ig.unstim,
                    newAnnotation = F,
                    expCondCheck = 'expCond.asthma',
                    expCondCheckFname = 'asthma_ig_unstim_MAST',
                    cellcluster = "B",
                    compGroup = "Asthmatic/Non-asthmatic",
                    deMethod = 'MAST',
                    min.cells.group = 6,
                    topNo = 20
)


################# DEGs (MAST) for B asthma against B non-asthma unstim_0h ################# 
rds.unstim <- subset(rds, subset = (expCond.stimuli == "unstim"))


# only consider ig and unstim -> compare Asthmatic vs Non-asthmatic
names(table(rds.unstim$expCond.asthma))


getClusterExpCondDe(resDir = deResDir,
                    rds = rds.unstim,
                    newAnnotation = F,
                    expCondCheck = 'expCond.asthma',
                    expCondCheckFname = 'asthma_unstim_MAST',
                    cellcluster = "B",
                    compGroup = "Asthmatic/Non-asthmatic",
                    deMethod = 'MAST',
                    min.cells.group = 6,
                    topNo = 20
)



################# DEGs (MAST) for B asthma against B non-asthma Ig_4h ################# 
rds.ig <- subset(rds, subset = (expCond.stimuli.time == "Ig_4h"))


# names(table(Idents(rds))) # cell types

# only consider ig and unstim -> compare Asthmatic vs Non-asthmatic
names(table(rds.ig$expCond.asthma))


getClusterExpCondDe(resDir = deResDir,
                    rds = rds.ig,
                    newAnnotation = F,
                    expCondCheck = 'expCond.asthma',
                    expCondCheckFname = 'asthma_ig_4h_MAST',
                    cellcluster = "B",
                    compGroup = "Asthmatic/Non-asthmatic",
                    deMethod = 'MAST',
                    min.cells.group = 6,
                    topNo = 20
)




################# DEGs (MAST) for B asthma against B non-asthma Ig_18h ################# 
rds.ig <- subset(rds, subset = (expCond.stimuli.time == "Ig_18h"))


# names(table(Idents(rds))) # cell types

# only consider ig and unstim -> compare Asthmatic vs Non-asthmatic
names(table(rds.ig$expCond.asthma))


getClusterExpCondDe(resDir = deResDir,
                    rds = rds.ig,
                    newAnnotation = F,
                    expCondCheck = 'expCond.asthma',
                    expCondCheckFname = 'asthma_ig_18h_MAST',
                    cellcluster = "B",
                    compGroup = "Asthmatic/Non-asthmatic",
                    deMethod = 'MAST',
                    min.cells.group = 6,
                    topNo = 20
)


integration_2_leiden_de_asthma_unstim.R


################# DEGs (MAST) for CD4 expressing IL4 or IL5 or IL13 regradless of conditions ################# 
grep('^IL4$', rownames(rds$RNA$data), value = FALSE) # exaxt match 
grep('^IL5$', rownames(rds$RNA$data), value = FALSE) # exaxt match 
grep('^IL13$', rownames(rds$RNA$data), value = FALSE) # exaxt match 

rds@meta.data$IL <- 'None'

rds@meta.data$IL[which(rds$RNA$data[rownames(rds$RNA$data)[23360], ] != 0)] <- 'IL4|IL5|IL13'
rds@meta.data$IL[which(rds$RNA$data[rownames(rds$RNA$data)[22774], ] != 0)] <- 'IL4|IL5|IL13'
rds@meta.data$IL[which(rds$RNA$data[rownames(rds$RNA$data)[21228], ] != 0)] <- 'IL4|IL5|IL13'


getClusterExpCondDe(
    resDir = deResDir,
    rds = rds,
    newAnnotation = F,
    expCondCheck = "IL",
    expCondCheckFname = "IL_MAST",
    compGroup = "IL4|IL5|IL13/None",
    cellcluster = "CD4 T",
    deMethod = "MAST",
    min.cells.group = 6,
    topNo = 20
  )










## look into seurat clusters ----
################# DEGs (MAST) across major seurat clusters regradless of conditions ################# 
Idents(rds) <- "all_seurat_clusters"

for (cluster in unique(rds@meta.data$seurat_clusters)) {
  
  # Create a new column for the current cluster
  column_name <- paste0("seurat_cluster_", cluster)
  rds@meta.data[[column_name]] <- ifelse(rds@meta.data$seurat_clusters == cluster, as.character(cluster), "other")
  
  
  # Run differential expression analysis
  getClusterExpCondDe(
    resDir = deResDir,
    rds = rds,
    newAnnotation = F,
    expCondCheck = column_name,
    expCondCheckFname = paste0(column_name, "_MAST"),
    compGroup = paste0(as.character(cluster), "/", "other"),
    # cellcluster = as.character(cluster),
    deMethod = "MAST",
    min.cells.group = 6,
    topNo = 20
  )
}










#### 2.0 run locally ####
################# Gene Counts for DEGs (MAST) across major cell types, do "CD3_18h", "CD3_4h", "Ig_18h", "Ig_4h", "LPS_18h", "LPS_4h" vs "unstim_0h" #############
#### prepare de gene cnts df that includes up, down, unique_up, unique_down at 2 time points ####

# Load DEGs for each contrast 
excel.dir <- "results_wOrgClusterAnnotation_DEGs/stimuli_time_MAST/"
file.all <- list.files(paste0(deResDir,excel.dir))
# could be "full_allClusters.xlsx"
files <- file.all[grepl("full_.*Clusters\\.xlsx", file.all, ignore.case = TRUE)]
contrast <- gsub("expCondCompDeMarkers_|_full_allClusters.xlsx", "", files)

# save each cell type to a tibble for each comparison 
for (i in 1:length(files)) {
  excel_file <- paste0(deResDir, excel.dir, files[i])
  sheet_names <- excel_sheets(excel_file)
  for (sheet_name in sheet_names) {
    assign(paste0(contrast[i], "_", sheet_name), read_excel(excel_file, sheet = sheet_name))
  }
}

# get pcgs
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

# parse DEGs 
getUpSetInput <- function(degtbl = "AllCellTypes$", logFC.cutoff = 0.25, fdr.cutoff = 0.05) {
  grp <- ls(pattern = degtbl, envir = parent.frame())
  
  # filter by threshold and get a list of deduplicated gene names for each 
  sig.degs.list <- lapply(grp, function(x) {
    degs <- get(x)
    # sig.degs      <- degs %>% dplyr::filter(p_val_adj<=fdr.cutoff & abs(avg_log2FC)>=logFC.cutoff) %>% dplyr::arrange(desc(avg_log2FC))
    sig.degs      <- degs %>% dplyr::filter(p_val_adj<=fdr.cutoff & abs(avg_log2FC)>=logFC.cutoff & ...1 %in% pcgs)  %>% dplyr::arrange(desc(avg_log2FC))
    sig.degs.up   <- sig.degs %>% dplyr::filter(avg_log2FC>0)
    sig.degs.dn   <- sig.degs %>% dplyr::filter(avg_log2FC<0)
    sig.degs.dst     <- dplyr::rename(distinct(sig.degs[,1]), sig=1)$sig #Note that some up/dn genes are concurrent so the total number of distinct sig genes are smaller than the sum of sig up and sig dn
    sig.degs.up.dst  <- dplyr::rename(distinct(sig.degs.up[,1]), sigUp=1)$sigUp
    sig.degs.dn.dst  <- dplyr::rename(distinct(sig.degs.dn[,1]), sigDn=1)$sigDn
    sig.list <- list(sig.degs.dst, sig.degs.up.dst, sig.degs.dn.dst)
    
    print(x)
    print(paste0("No.DEGs(both): ", dim(sig.degs)[1], " -> ", length(sig.degs.dst)))
    print(paste0("No.DEGs(up): ", dim(sig.degs.up)[1], " -> ", length(sig.degs.up.dst)))
    print(paste0("No.DEGs(dn): ", dim(sig.degs.dn)[1], " -> ", length(sig.degs.dn.dst)))
    
    return(sig.list)
  })
  
  # transform to upset input
  upsetInput.up <- list()
  upsetInput.dn <- list()
  upsetInput.both <- list()
  
  compName <- gsub("-", "/", grp)
  compName <- gsub("_cluster_.*$","",compName)
  
  for (i in 1:length(sig.degs.list)) {
    upsetInput.both[[compName[i]]] <- sig.degs.list[[i]][[1]]
    upsetInput.up[[compName[i]]] <- sig.degs.list[[i]][[2]]
    upsetInput.dn[[compName[i]]] <- sig.degs.list[[i]][[3]]
  }
  
  return(list(both = upsetInput.both, up = upsetInput.up, dn = upsetInput.dn))
  return(sig.degs.list)
}


combined_df <- data.frame()
cell_types <- c("Mono/Mph",
                "B", 
                "CD4 T", 
                "CD8 T", 
                "NK",
                "DC", 
                "EC")

groups <- c("LPS", "Ig", "CD3")



for (celltype in cell_types){
  # print(celltype)
  ct <- gsub("/| |-","",celltype) # replace any / space -
  upsetInput <- getUpSetInput(ct)
  
  # "CD3_18h/unstim_0h_full_5SelClusters.xlsx" will be marked as "CD3_18h/unstim_0h" for selected groups
  names(upsetInput$up) <- sub("_full.*$", "", names(upsetInput$up))
  names(upsetInput$dn) <- sub("_full.*$", "", names(upsetInput$dn))
  names(upsetInput$both) <- sub("_full.*$", "", names(upsetInput$both))
  
  names(upsetInput$up) <- lapply(names(upsetInput$up), function(x){ return(paste0("up_",x)) })
  names(upsetInput$dn) <- lapply(names(upsetInput$dn), function(x){ return(paste0("down_",x)) })
  upsetInput.merge <- c(upsetInput$up, upsetInput$dn)
  for (group in groups){
    # print(group)
    upsetInput.group <- upsetInput.merge[grep(group, names(upsetInput.merge), value = TRUE)]
    
    down_4h <- upsetInput.group[[paste0("down_", group, "_4h/unstim_0h")]]
    down_18h <- upsetInput.group[[paste0("down_", group, "_18h/unstim_0h")]]
    up_4h <- upsetInput.group[[paste0("up_", group, "_4h/unstim_0h")]]
    up_18h <- upsetInput.group[[paste0("up_", group, "_18h/unstim_0h")]]
    
    unique_down_4h <- setdiff(down_4h, down_18h)
    unique_down_18h <- setdiff(down_18h, down_4h)
    unique_up_4h <- setdiff(up_4h, up_18h)
    unique_up_18h <- setdiff(up_18h, up_4h)
    
    # Add these unique sets back to upsetInput.group
    upsetInput.group[[paste0("unique_down_", group, "_4h/unstim_0h")]] <- unique_down_4h
    upsetInput.group[[paste0("unique_down_", group, "_18h/unstim_0h")]] <- unique_down_18h
    upsetInput.group[[paste0("unique_up_", group, "_4h/unstim_0h")]] <- unique_up_4h
    upsetInput.group[[paste0("unique_up_", group, "_18h/unstim_0h")]] <- unique_up_18h
    
    upsetInput.group.cnts <- data.frame(
      up_down = c("up", "up", "down", "down",
                  "unique_up", "unique_up", "unique_down", "unique_down"),
      gene_cnts = c(length(up_4h), length(up_18h), length(down_4h), length(down_18h),
                    length(unique_up_4h), length(unique_up_18h), length(unique_down_4h), length(unique_down_18h)),
      treatment = c(paste0(group, "_4h"), paste0(group, "_18h"), paste0(group, "_4h"), paste0(group, "_18h"),
                    paste0(group, "_4h"), paste0(group, "_18h"), paste0(group, "_4h"), paste0(group, "_18h")),
      celltype = rep(celltype, 8))
    
    combined_df <- rbind(combined_df, upsetInput.group.cnts)
  }
}


write_xlsx(combined_df, path = paste0(deResDir, excel.dir, 'de_gene_cnts.xlsx'), col_names = T)




# stop("DE run on randi stop here.")

#### make a plot for gene cnts ####

combined_df <- read_excel(paste0(deResDir,excel.dir, 'de_gene_cnts.xlsx'), , sheet = "Sheet1")


combined_df <- combined_df %>%
  pivot_wider(names_from = up_down, 
              values_from = gene_cnts) %>%
  mutate(up_share = up - unique_up) %>%
  mutate(down_share = down - unique_down)

"#E63A47"
"#A9DADC"

combined_df$treatment <- factor(combined_df$treatment, 
                                levels =  c("LPS_4h", "Ig_4h", "CD3_4h", 
                                            "LPS_18h", "Ig_18h", "CD3_18h"))


pdf(file = paste0(deResDir, excel.dir, "de_gene_cnts.pdf"), width = 8.35, height = 4.36)

ggplot(combined_df, aes(x = factor(celltype, levels=c("B", 
                                                       "CD4 T", 
                                                       "CD8 T", 
                                                       "NK",
                                                       "Mono/Mph",
                                                      "EC",
                                                      "DC")))) +
  geom_bar(aes(y = up, fill = "up"), stat = "identity", position = "identity") +
  geom_bar(aes(y = -down, fill = "down"), stat = "identity", position = "identity") +
  geom_text(aes(label = up, y = up), vjust = -0.5, color = "#FE4088", size = 2.7) +
  geom_text(aes(label = down, y = -down), vjust = 1.5, color = "#5F96FA", size = 2.7) +
  geom_bar(aes(y = up_share, fill = "up_share"), stat = "identity", position = "identity") +
  geom_bar(aes(y = -down_share, fill = "down_share"), stat = "identity", position = "identity") +
  scale_fill_manual(values = c("up" = "#FE4088", "down" = "#5F96FA", 
                               "up_share" = "#CC336E", "down_share" = "#4B6EAD"),
                    breaks = c("up", "up_share", "down_share", "down"),
                    labels = c("Up-regulated", "Shared up-regulated",
                               "Shared down-regulated", "Down-regulated"),
                    guide = guide_legend(title = "DEGs")) +
  labs(x = "", y = "Number of DEGs") +
  theme_light() +
  theme(panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 11, color = "black"), 
        strip.background = element_rect(fill = "white")) +
  ylim(-6500, 4000) +
  facet_wrap(~ treatment,
             labeller = labeller(treatment = c(
               "LPS_4h" = "LPS_4h vs Unstim_0h",
               "Ig_4h" = "Ig_4h vs Unstim_0h",
               "CD3_4h" = "CD3_4h vs Unstim_0h",
               "LPS_18h" = "LPS_18h vs Unstim_0h",
               "Ig_18h" = "Ig_18h vs Unstim_0h",
               "CD3_18h" = "CD3_18h vs Unstim_0h"
             )))
dev.off()




# change color code
ggplot(combined_df, aes(x = celltype)) +
  geom_bar(aes(y = up, fill = "up"), stat = "identity", position = "identity") +
  geom_bar(aes(y = -down, fill = "down"), stat = "identity", position = "identity") +
  geom_text(aes(label = up, y = up), vjust = -0.5, color = "#c1121f", size = 2.5) +
  geom_text(aes(label = down, y = -down), vjust = 1.5, color = "#669bbc", size = 2.5) +
  geom_bar(aes(y = up_share, fill = "up_share"), stat = "identity", position = "identity") +
  geom_bar(aes(y = -down_share, fill = "down_share"), stat = "identity", position = "identity") +
  scale_fill_manual(values = c("up" = "#c1121f", "down" = "#669bbc", 
                               "up_share" = "#780000", "down_share" = "#003049"),
                    breaks = c("up", "up_share", "down_share", "down"),
                    labels = c("up", "Shared up", "Shared down", "down"),
                    guide = guide_legend(title = "Category")) +
  labs(x = "Cell type", y = "Number of \ndifferentially expressed genes") +
  theme_light() +
  theme(panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 11),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text = element_text(size = 9),
        strip.text = element_text(size = 11, color = "black"), 
        strip.background = element_rect(fill = "white")) +
  ylim(-6500, 3500) +
  facet_wrap(~ factor(treatment, levels = c("LPS_4h", "Ig_4h", "CD3_4h", 
                                            "LPS_18h", "Ig_18h", "CD3_18h")))

