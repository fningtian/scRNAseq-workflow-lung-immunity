rm(list=ls())

library(Seurat)
library(dplyr)
library(cowplot)
library(patchwork)
library(gridExtra)
library(DoubletDecon)
library(ggplot2)
library(tidyr)
library(RColorBrewer)

## qc on samples -----------------------

# par(mar=c(1,1,1,1))
# dev.off()
x11()
##----------------------------------------------------------------------------------------
Improved_Seurat_Pre_Process <- function (seuratObject, num_genes = 50, write_files = FALSE, data_type = "counts") {
  
  version = packageVersion("Seurat")
  # seuratObject = UpdateSeuratObject(object = seuratObject)
  if (data_type == "counts") {
    expression = as.data.frame(seuratObject@assays[["RNA"]]$counts)
  }
  else if (data_type == "data") {
    expression = as.data.frame(seuratObject@assays[["RNA"]]$data)
  }
  else if (data_type == "scaled.data") {
    expression = as.data.frame(seuratObject@assays[["RNA"]]$scale.data)
  }
  seuratObject.markers = FindAllMarkers(object = seuratObject,
                                        only.pos = TRUE, min.pct = 0.25)
  if (version >= package_version(x = "3.9.9")) {
    genes = seuratObject.markers %>% group_by(cluster) %>%
      top_n(n = num_genes, wt = avg_log2FC)
  }
  else if (version >= package_version(x = "3.0.0") && version <
           package_version(x = "3.9.9")) {
    genes = seuratObject.markers %>% group_by(cluster) %>%
      top_n(n = num_genes, wt = avg_logFC)
  }
  else {
    print("This function only works with Seurat 3 or 4. Please update Seurat.")
  }
  clusters = as.data.frame(Idents(object = seuratObject))
  colnames(expression) = gsub("-", ".", colnames(expression))
  if (class(clusters[, 1]) == "character") {
    clusters[, 1] = gsub("-", ".", clusters[, 1])
  }
  else {
    row.names(clusters) = gsub("-", ".", row.names(clusters))
  }
  if (class(genes$cluster) == "factor") {
    if (min(as.numeric(as.character(genes$cluster))) == 0) {
      genes$cluster = as.numeric(as.character(genes$cluster)) +
        1
      clusters[, 1] = as.numeric(as.character(clusters[,
                                                       1])) + 1
    }
    else if (min(as.numeric(as.character(genes$cluster))) ==
             1) {
      genes$cluster = as.numeric(as.character(genes$cluster))
      clusters[, 1] = as.numeric(as.character(clusters[,
                                                       1]))
    }
    else {
      print("Unexpected cluster numbering scheme. Cluster numbers are expected to be continuous numbers starting from either 0 or 1. Please check conversion for correctness following this function.")
    }
  }
  clusters2 = clusters[order(clusters[, 1]), , drop = FALSE]
  expression = expression[row.names(clusters2)]
  genes2 = genes[order(genes$cluster), ]
  allgenes = expression
  expression = expression[row.names(expression) %in% as.character(genes$gene),
  ]
  geneOrder = intersect(genes2$gene, as.character(row.names(expression)))
  expression = expression[match(geneOrder, row.names(expression)),
  ]
  allgenes = rbind(clusters2[, 1], allgenes)
  expression = rbind(clusters2[, 1], expression)
  row.names(allgenes)[1] = "column_clusters-flat"
  row.names(expression)[1] = "column_clusters-flat"
  genes3 = genes2[match(geneOrder, genes2$gene), ]
  rowToAdd = c(NA, genes3$cluster)
  rowToAdd2 = rep(NA, nrow(allgenes))
  expression = cbind(rowToAdd, expression)
  allgenes = cbind(rowToAdd2, allgenes)
  colnames(expression)[1] = "row_clusters-flat"
  colnames(allgenes)[1] = "row_clusters-flat"
  groups = cbind(as.numeric(expression[1, 2:ncol(expression)]),
                 as.numeric(expression[1, 2:ncol(expression)]))
  row.names(groups) = as.character(colnames(expression)[2:ncol(expression)])
  if (write_files == TRUE) {
    write.table(expression, "ICGS_expression.txt", sep = "\t")
    write.table(allgenes, "ICGS_fullExpression.txt", sep = "\t")
    write.table(groups, "ICGS_groups.txt", sep = "\t", col.names = F)
  }
  return(list(newExpressionFile = expression, newFullExpressionFile = allgenes,
              newGroupsFile = groups))
  
}

##----------------------------------------------------------------------------------------
findDoublets <- function(metadata, multiomics = F, extraFilter=F, genomeSpecies=NULL, doubletDeconRhop=0.5, doubletDeconPMF=F, doubletDeconNoCore=-1, resFilename=NULL) {
  ## ---
  if (is.null(genomeSpecies)) genomeSpecies <- 'human'
  if (is.null(resFilename)) resFilename <- 'doublets_results'
  doubletDeconRhop               <- as.numeric(doubletDeconRhop)
  doubletDeconPMF                <- as.logical(doubletDeconPMF)
  multiomics                     <- as.logical(multiomics)
  extraFilter                    <- as.logical(extraFilter)
  ## ---
  cellrangerResList              <- meata2list(metadata = metadata)
  ## ---
  ## prepare results saving directory, if not exist, create one
  resDir  <- sprintf('%s/%s', getwd(), resFilename)
  if(!dir.exists(resDir)) dir.create(resDir)
  print(sprintf('Doublets identification results will be saved in %s', resDir))
  ## intermediate 'resProcessDir' under/inside provided 'resDir'
  resProcessDir                  <- paste(resDir, 'doubletDecon_preProcessed_results', sep = '/')
  if(!dir.exists(resProcessDir)) dir.create(resProcessDir)
  ## ---
  ## loop over provided 'cellrangerResList' to identify/estimate doublets in each list item of provided 'cellrangerResList'
  for (x in 1:length(cellrangerResList)) {
    ## ---
    ## 1. creat seurat object as DoubletDecon suggested
    print('---')
    print(sprintf("Processing sample '%s'.", as.character(cellrangerResList[[x]])))
    
    # funing: change the line below to include the rds format
    # cellrangerCountsOrg          <- Seurat::Read10X(data.dir = cellrangerResList[[x]])
    
    if (tools::file_ext(cellrangerResList[[x]])=='rds') {
      print("read in count data in rds format")
      rdsOrg                     <- readRDS(file = as.character(cellrangerResList[[x]]))
      cellrangerCountsOrg        <- rdsOrg@assays$RNA$counts
    }else {
      cellrangerCountsOrg          <- Seurat::Read10X(data.dir = cellrangerResList[[1]])
    }
    
    # funing: also here
    # if (multiomics) {
    #   cellrangerCounts           <- cellrangerCountsOrg$`Gene Expression`
    # } else {
    #   cellrangerCounts           <- cellrangerCountsOrg
    # }
    
    
    if (multiomics) {
      if(tools::file_ext(cellrangerResList[[x]])=='rds'){
        cellrangerCounts        <- cellrangerCountsOrg
      }else{
        cellrangerCounts           <- cellrangerCountsOrg$`Gene Expression`
      }
    } else {
      cellrangerCounts           <- cellrangerCountsOrg
    }
    
    
    
    
    print(sprintf('Step1: Orignially it has %s cells and %s features originated from cellranger to import into seurat object', length(cellrangerCounts@Dimnames[[2]]), length(cellrangerCounts@Dimnames[[1]]) ))
    ## include feature detected in at least 'min.cells = 3', and include cells where at least 'min.features = 200' detected
    ## ---
    if(extraFilter) {
      if (!'filterFname' %in% colnames(metadata)) stop("Option extraFilter is on, but no filter files is provided in the metadata table column 'filterFname'.")
      if (file_ext(metadata$filterFname[x])=='csv') {
        filterRes   <- read.csv(file = metadata$filterFname[x], header = T)
      } else  if (file_ext(metadata$filterFname[x])=='txt') {
        filterRes   <- read.delim(file = metadata$filterFname[x], header = T, sep = '\t')
      }
      if (sum(grepl('barcode|filter',colnames(filterRes)))!=2) stop("Please make sure columns 'barcode' & 'filter' are inside provided ''.")
      if (!all(filterRes$barcode %in% colnames(cellrangerCounts) )) stop("provided filter files does not include all corresponding cells information.")
      filter.index  <- grep('filter', colnames(filterRes))
      print(sprintf("%s cells will be filtered based on input 'filterFname' in metadata table; %s cells will be used for next steps analysis.",
                    sum(filterRes[,filter.index]==as.logical(T)), sum(filterRes[,filter.index]==as.logical(F)) ))
      cellrangerCountsFilter <- cellrangerCounts[, match(filterRes$barcode[filterRes[,filter.index]==as.logical(F)], colnames(cellrangerCounts))]
      cellrangerCounts       <- cellrangerCountsFilter
    }
    ## ---
    seuratObject                 <- Seurat::CreateSeuratObject(counts = cellrangerCounts,  project = names(cellrangerResList)[x], min.cells = 3, min.features = 200)
    ## add metadata feature into object, here is 'expCond'
    seuratObject                 <- Seurat::AddMetaData(object = seuratObject,  col.name = 'expCond', metadata = as.factor(names(cellrangerResList)[x]))
    ## -
    seuratObject[['percent.mt']] <- Seurat::PercentageFeatureSet(object = seuratObject, pattern = as.character(mtPatten(genomeSpecies)) )
    seuratObject                 <- Seurat::NormalizeData(seuratObject, normalization.method = "LogNormalize", scale.factor = 10000)
    print(sprintf('Complete Log Normlization'))
    seuratObject                 <- Seurat::FindVariableFeatures(seuratObject, selection.method = 'vst', nfeatures = 2000)
    print(sprintf('Complete vst top variable feature identification'))
    if (dim(seuratObject@meta.data)[1] > 10000) {
      seuratObject               <- Seurat::ScaleData( object = seuratObject, features = VariableFeatures(seuratObject) )
    } else {
      seuratObject               <- Seurat::ScaleData(object = seuratObject, features = rownames(seuratObject))
    }
    print(sprintf('Complete data centering/scaling'))
    seuratObject                 <- Seurat::RunPCA(object = seuratObject)
    print(sprintf('Complete PCA clustering'))
    seuratObject                 <- Seurat::RunUMAP(seuratObject, dims = 1:10)
    print(sprintf('Complete UMAP clustering'))
    seuratObject                 <- Seurat::RunTSNE(seuratObject, dims = 1:10)
    print(sprintf('Complete tSNE clustering'))
    seuratObject                 <- Seurat::FindNeighbors(seuratObject, dims = 1:10)
    print(sprintf('Complete finding neighbors'))
    seuratObject                 <- Seurat::FindClusters(seuratObject, resolution = 0.8)
    print(sprintf('Complete finding clustering'))
    print('END step1 for orignal seurat processing')
    ## 2. Improved_Seurat_Pre_Process() from DoubletDecon on established seurat object, and output results to 'resDir' for next step usage
    print(sprintf("Start Step2: 'improve seurat pre process'."))
    filename                     <- names(cellrangerResList)[x]
    newFiles                     <- Improved_Seurat_Pre_Process(seuratObject = seuratObject, num_genes=50, write_files=FALSE)
    write.table(newFiles$newExpressionFile, paste0(resProcessDir, '/', filename, "_expression"), sep="\t")
    write.table(newFiles$newFullExpressionFile, paste0(resProcessDir, '/', filename, "_fullExpression"), sep="\t")
    write.table(newFiles$newGroupsFile, paste0(resProcessDir, '/', filename , "_groups"), sep="\t", col.names = F)
    print(sprintf("END Step2: 'improve seurat pre process'."))
    ## 3.1 DoubletDecon doublets detection with Main_Doublet_Decon() based on centroid method
    print('=========')
    print(sprintf("Start Step3: DoubletDecon (centroids & medoids) use rhop = %s based on %s genome with PMF = %s", doubletDeconRhop, as.character(doubletDeconSpecies(genomeSpecies)), doubletDeconPMF ))
    print('Start Step 3.1: DoubletDecon centroid detection')
    doubletDeconResCentroids     <- DoubletDecon::Main_Doublet_Decon(rawDataFile = paste0(resProcessDir, '/', filename, "_expression"),
                                                                     groupsFile = paste0(resProcessDir, '/', filename , "_groups"),
                                                                     filename = filename,
                                                                     location = paste(resProcessDir, 'Centroids_', sep = '/'),
                                                                     removeCC = F, ## default is FALSE
                                                                     species = as.character(doubletDeconSpecies(genomeSpecies)), ## default is 'mmu'
                                                                     rhop = doubletDeconRhop, ## Default is 1, x in mean+x*SD to determine upper cutoff for correlation in the blacklist.
                                                                     PMF = doubletDeconPMF, ## default = T, Use step 2 (unique gene expression) in doublet determination criteria
                                                                     useFull =  F, ## default = F, Use full gene list for PMF analysis
                                                                     heatmap = F,
                                                                     centroids = TRUE,
                                                                     nCores = doubletDeconNoCore)
    print('Doublet centroids detection results table:')
    print( table(doubletDeconResCentroids$DRS_doublet_table$isADoublet) )
    print('END Step 3.1: DoubletDecon centroid detection')
    ## 3.2 DoubletDecon doublets detection with Main_Doublet_Decon() based on medoids method
    print('Start Step 3.2: DoubletDecon medoids detection')
    print('=========')
    doubletDeconResMedoids       <- DoubletDecon::Main_Doublet_Decon(rawDataFile = paste0(resProcessDir, '/', filename, "_expression"),
                                                                     groupsFile = paste0(resProcessDir, '/', filename , "_groups"),
                                                                     filename = filename,
                                                                     location = paste(resProcessDir, 'Medoids_', sep = '/'),
                                                                     removeCC = F, ## default is FALSE
                                                                     species = as.character(doubletDeconSpecies(genomeSpecies)), ## default is 'mmu'
                                                                     rhop = doubletDeconRhop, ## x in mean+x*SD to determine upper cutoff for correlation in the blacklist. Default is 1
                                                                     PMF = doubletDeconPMF, ## default = T, Use step 2 (unique gene expression) in doublet determination criteria
                                                                     useFull =  F, ## default = F, Use full gene list for PMF analysis
                                                                     heatmap = F,
                                                                     centroids = FALSE,
                                                                     nCores = doubletDeconNoCore)
    print('Doublet medoids detection results table:')
    print( table(doubletDeconResMedoids$DRS_doublet_table$isADoublet) )
    print('=========')
    print('END Step 3.2: DoubletDecon medoids detection')
    print(sprintf("Doublets identification results ('doubletDeconResCentroids' & 'doubletDeconResMedoids') saved in '%s'.", as.character(file.path(resDir, sprintf('%s_doubletDeconRes.Rdata', filename)))))
    save(doubletDeconResCentroids, doubletDeconResMedoids, file = file.path(resDir, sprintf('%s_doubletDeconRes.Rdata', filename)) )
    ## -
    print('Start Step 4: Doublets identification results summary')
    ## 4. summarize doubletDecon detection results.
    ## 4.1 Medoids detection summary
    doubletNoSummaryMedoids             <- as.data.frame(table(doubletDeconResMedoids$DRS_doublet_table$isADoublet))
    doubletNoSummaryMedoids$Per         <- round(x = doubletNoSummaryMedoids[,2] / sum(doubletNoSummaryMedoids[,2]) * 100, digits = 1)
    colnames(doubletNoSummaryMedoids)   <- c('Doublet', 'Medioids detection No', 'Medioids detection Per')
    print('Doublet detection results table with medoids:')
    print( table(doubletDeconResMedoids$DRS_doublet_table$isADoublet) )
    ## Centroids detection
    doubletNoSummaryCentroids           <- as.data.frame(table(doubletDeconResCentroids$DRS_doublet_table$isADoublet))
    doubletNoSummaryCentroids$Per       <- round(x = doubletNoSummaryCentroids[,2] / sum(doubletNoSummaryCentroids[,2]) * 100, digits = 1)
    colnames(doubletNoSummaryCentroids) <- c('Doublet', 'Centroids detection No', 'Centroids detection Per')
    print('Doublet detection results table with centroids:')
    print( table(doubletDeconResCentroids$DRS_doublet_table$isADoublet) )
    if (sum(doubletNoSummaryMedoids[,2]) != sum(doubletNoSummaryCentroids[,2])) stop('Error: total number of cells different from medoids and centroids detection methods.')
    medoidsDoubletCells = rownames(doubletDeconResMedoids$DRS_doublet_table %>% dplyr:: filter(isADoublet == TRUE))
    save(medoidsDoubletCells, file = file.path(resDir, sprintf("%s_medoids_doublet_cells_name.Rdata", filename) ) )
    ## -
    centroidsDoubletCells = rownames(doubletDeconResCentroids$DRS_doublet_table %>% dplyr:: filter(isADoublet == TRUE))
    save(centroidsDoubletCells, file = file.path(resDir, sprintf("%s_centroids_doublet_cells_name.Rdata", filename) ) )
    ## -
    print('Start overlapping medoids and centroids doublet detection results')
    olRes <- gplots::venn(list('medoids'  = rownames(doubletDeconResMedoids$DRS_doublet_table %>% dplyr::filter(isADoublet == TRUE)),
                               'centroid' = rownames(doubletDeconResCentroids$DRS_doublet_table %>% dplyr::filter(isADoublet == TRUE)) ))
    plot(olRes)
    ## -
    pdf(file = file.path(resDir, sprintf("%s_doubletDecon_mediod_centroid_vennComp.pdf", filename) ), width = 3, height = 4 )
    plot(olRes)
    dev.off()
    ## -
    olDoubletCells <- attr(olRes, 'intersection')[[grep('medoids:centroid', names(attr(olRes, 'intersection')) )]]
    print('Finish overlapping medoids and centroids doublet detection results')
    ## -
    doubletNoSummaryComb                <- data.frame(Doublet = c('TRUE', 'FALSE'),
                                                      No = c(length(olDoubletCells), sum(doubletNoSummaryCentroids[,2])-length(olDoubletCells) ))
    doubletNoSummaryComb$Per            <- round(x = doubletNoSummaryComb[,2] / sum(doubletNoSummaryComb[,2]) * 100, digits = 1)
    colnames(doubletNoSummaryComb)      <- c('Doublet', 'Combined detection No', 'Combined detection Per')
    doubletNoSummary                    <- dplyr::left_join(doubletNoSummaryMedoids, doubletNoSummaryCentroids, by = 'Doublet') %>% dplyr::left_join(doubletNoSummaryComb, by = 'Doublet')
    write.table(x = doubletNoSummary, file = file.path(resDir, sprintf("%s_doublets_no_summary.txt", filename)), quote = F, sep = '\t', row.names = F, col.names = T)
    ## -
    save(olDoubletCells, file = file.path(resDir, sprintf("%s_OL_doublet_cells_name.Rdata", filename) ) )
    print('END Step 4: Doublets identification results summary')
    print('=========')
    print(sprintf("COMPLETE doublets identification for sample '%s'.", as.character(cellrangerResList[[x]])))
    print('---END---END---END---')
    ## ---
  }
  return(resDir)
  ## ---
}

## Two minor fns to return mitochodrial content search pattern
## and Main_Doublet_Decon() species based on genome input used in findDoublets()
## ---
mtPatten            <- function(genomeSpecies) {
  if (genomeSpecies == 'human') return('^MT-')
  if (genomeSpecies == 'mouse') return('^mt-')
  if (genomeSpecies == 'human_mouse') return('^MT-|mt-')
}
## ---
doubletDeconSpecies <- function(genomeSpecies) {
  if (genomeSpecies == 'human') return('hsa')
  if (genomeSpecies == 'mouse') return('mmu')
  if (genomeSpecies == 'rat') return('rno')
}

##----------------------------------------------------------------------------------------
processQC <- function(metadata, multiomics = F, extraFilter=F, resDirName=NULL, genomeSpecies=NULL, minCells=3, minFeatures=200, mtFiltering=F, mtPerCutoff=NULL, nfeatures = 5000) {
  ##--------------------------------------------------------------------------------------##
  if (!all(c("sample", "path", 'expCond1',"doubletsRmMethod" ) %in% colnames(metadata))) stop('Please provide metadata table with at least 4 columns: sample, path, expCond1, and doubletsRmMethod')
  ## ---
  if (is.null(resDirName)) resDirName <- as.character('scRICA_results')
  # if (!is.list(cellrangerResList)) stop("Please provide list item of 'cellrangerResList', which is required.")
  minCells                       <- as.numeric(minCells)
  minFeatures                    <- as.numeric(minFeatures)
  mtFiltering                    <- as.logical(mtFiltering)
  multiomics                     <- as.logical(multiomics)
  extraFilter                    <- as.logical(extraFilter)
  if (is.null(genomeSpecies)) genomeSpecies <- as.character('human')
  if (mtFiltering & is.null(mtPerCutoff)) stop("Mitochondrial content filtering option ('mtFiltering') is on, please provide corresponding Mitochondrial content percentage filtering option in 'mtPerCutoff'.  ")
  mtPerCutoff                    <- as.numeric(mtPerCutoff)
  ##--------------------------------------------------------------------------------------##
  # if ( 'expCond1' %in% colnames(metadata) & !'expCond2' %in% colnames(metadata)) {
  #   print(sprintf('Only 1 experimental condition are provided with %s experimental factor levels for comparisons.', length(levels(factor(metadata$expCond1))) ))
  # } else if ( 'expCond1' %in% colnames(metadata) & 'expCond2' %in% colnames(metadata) ) {
  #   print(sprintf('2 experimental conditions are provided, experimental condition 1 has %s experimental factor levels, and experimental condition 2 has %s experimental factor levels for comparisons', length(levels(factor(metadata$expCond1))), length(levels(factor(metadata$expCond2))) ))
  # } else if ( !'expCond1' %in% colnames(metadata) & !'expCond2' %in% colnames(metadata) ) {
  #   print(sprintf('No experimental condition factors are provided, the analysis will be conducted only based samples integration'))
  # } else if ( !'expCond1' %in% colnames(metadata) & 'expCond2' %in% colnames(metadata)) {
  #   print(sprintf('Only 1 experimental condition are provided with %s experimental factor levels for comparisons.', length(levels(factor(metadata$expCond2))) ))
  # }
  print("-=-=-=-=-=-=-")
  expCond.no    <- length(grep('expCond', colnames(metadata)))
  expCond.index <- grep('expCond', colnames(metadata))
  print(sprintf("%s experimental condition are provided in the metadata columns: %s.", expCond.no, paste(colnames(metadata)[expCond.index], collapse = ', ')) )
  for (i in 1:expCond.no) {
    print(sprintf("%s: %s has %s experimental factor levels", i, colnames(metadata)[expCond.index], length(levels(factor(metadata[,expCond.index[i]]))) ))
  }
  print("-=-=-=-=-=-=-")
  ##--------------------------------------------------------------------------------------##
  ## 0. create main directory for results to be save in
  resDir                         <- paste(getwd(), resDirName, sep = '/')
  if (!dir.exists(resDir)) dir.create(resDir)
  ##--------------------------------------------------------------------------------------##
  ## run findDoublet and update metadata to include doublets detection results
  if ( !'doubletsResDir' %in% colnames(metadata) ) {
    metadataOrg                  <- metadata
    metadataOrg$doubletsRmMethod <- gsub('[]|[ ]', 'none', metadataOrg$doubletsRmMethod)
    metadataOrg$doubletsRmMethod[is.na(metadataOrg$doubletsRmMethod)] <- 'none'
    if(all(metadataOrg$doubletsRmMethod=='None'|metadataOrg$doubletsRmMethod=='none'|metadataOrg$doubletsRmMethod=='NONE')) {
      metadata                   <- metadataOrg
      metadata$doubletsResDir    <- as.character('NA')
    } else {
      md4doublet                 <- metadataOrg %>% dplyr::filter(doubletsRmMethod != 'None' ) %>% dplyr::filter(doubletsRmMethod != 'none' ) %>% dplyr::filter(doubletsRmMethod != 'NONE' )
      print('===================================================================')
      print('START: Doublets identification analysis before processing to the next step.')
      doubletsRes                <- findDoublets(metadata = md4doublet, multiomics = multiomics, genomeSpecies = genomeSpecies, resFilename = paste(resDirName, 'doublet_results', sep = '/') )
      print('END: Doublets identification, process to next step QC.')
      print('===================================================================')
      md4doublet$doubletsResDir  <- rep(as.character(doubletsRes), length(md4doublet$sample))
      doubletsMd                 <- md4doublet %>% dplyr::select(c('sample', 'doubletsResDir'))
      metadata                   <- dplyr::left_join(x = metadataOrg, y = doubletsMd, by = 'sample', )
    }
  } else {
    print('===================================================================')
    print('NO doublets identification analysis is needed, because doublets identification either conducted or not needed.')
    metadataOrg                  <- metadata
    metadataOrg$doubletsRmMethod <- gsub('[]|[ ]', 'none', metadataOrg$doubletsRmMethod)
    metadataOrg$doubletsRmMethod[is.na(metadataOrg$doubletsRmMethod)] <- 'none'
    metadataOrg$doubletsRmMethod[metadataOrg$doubletsRmMethod=='NONE'] <- 'none'
    metadataOrg$doubletsRmMethod[metadataOrg$doubletsRmMethod=='None'] <- 'none'
    metadataOrgDbNotNa           <- metadataOrg %>% dplyr::filter(doubletsRmMethod != 'None') %>% dplyr::filter(doubletsRmMethod != 'NONE') %>% dplyr::filter(doubletsRmMethod != 'none')
    if(dim(metadataOrgDbNotNa)[1]>0) {
      if (!all(dir.exists(metadataOrgDbNotNa$doubletsResDir))){
        print(sprintf("Stop at sample %s:", metadata$sample[which(!dir.exists(metadataOrgDbNotNa$doubletsResDir))]))
        stop("please provide correct corresponding 'doubletsResDir' in 'metadata' for column 'doubletsRmMethod' specified NOT as 'none'.")
      }
    }
    metadata                     <- metadataOrg
    print('Process to QC.')
    print('===================================================================')
  }
  # print(metadata)
  # print('*****************************')
  ##--------------------------------------------------------------------------------------##
  cellrangerResList              <- meata2list(metadata = metadata) ## by default names(cellrangerResList) = metadata$sample
  ## ---
  doubletsRmMethods              <- as.character(metadata$doubletsRmMethod)
  doubletsResDirs                <- as.character(metadata$doubletsResDir)
  ## 1. setup the original Seurat object into a list for each item in the 'cellrangerResList'
  # Sys.time()
  print(sprintf('Step 1: read in 10X data into Seurat object at %s', Sys.time()))
  seuratObjList                  <- list()
  for (x in 1:length(cellrangerResList)) {
    print('---===------------')
    print(sprintf("Processing sample %s: '%s'.", x, as.character(cellrangerResList[[x]])))
    if (dir.exists(cellrangerResList[[x]])) {
      print("read in count data in MEX format")
      cellrangerCountsOrg          <- Seurat::Read10X(data.dir = cellrangerResList[[x]])
    } else {
      if (tools::file_ext(cellrangerResList[[x]])=='h5' | tools::file_ext(cellrangerResList[[x]]) == 'hdf5') {
        print("read in count data in h5/hdf5 format")
        cellrangerCountsOrg        <- Seurat::Read10X_h5(filename = as.character(cellrangerResList[[x]]), use.names = TRUE, unique.features = TRUE)
      } else if (tools::file_ext(cellrangerResList[[x]])=='txt' | tools::file_ext(gsub('.gz', '', cellrangerResList[[x]])) == 'txt') {
        print("read in count data in txt format")
        cellrangerCountsOrg        <- read.delim2(file = as.character(cellrangerResList[[x]]))
      } else if (tools::file_ext(cellrangerResList[[x]])=='rds') {
        print("read in count data in rds format") # funing
        rdsOrg                     <- readRDS(file = as.character(cellrangerResList[[x]]))
        cellrangerCountsOrg        <- rdsOrg@assays$RNA$counts
      } else {
        stop("input file in metadata table cannot be read in, it should be any of these 3 formats: txt, hdf5, or MEX in a directory.")
      }
    }
    ##--------------------------------------------------------------------------------------##
    
    # funing: consider the rds format
    # if (multiomics) {
    #   cellrangerCounts           <- cellrangerCountsOrg$`Gene Expression`
    # } else {
    #   cellrangerCounts           <- cellrangerCountsOrg
    # }
    
    if (multiomics) {
      if(tools::file_ext(cellrangerResList[[x]])=='rds'){
        cellrangerCounts        <- cellrangerCountsOrg
      }else{
        cellrangerCounts           <- cellrangerCountsOrg$`Gene Expression`
      }
    } else {
      cellrangerCounts           <- cellrangerCountsOrg
    }
    
    
    
    # print(colnames(cellrangerCounts)[1:5])
    ## -
    if (dir.exists(cellrangerResList[[x]]) | tools::file_ext(cellrangerResList[[x]])=='h5' | tools::file_ext(cellrangerResList[[x]]) == 'hdf5') {
      print(sprintf('Originally it has %s cells and %s features originated from cellranger (MEX or hdf5 format) to import into Seurat object', length(cellrangerCounts@Dimnames[[2]]), length(cellrangerCounts@Dimnames[[1]]) ))
    } else if (tools::file_ext(cellrangerResList[[x]])=='txt' | tools::file_ext(gsub('.gz', '', cellrangerResList[[x]])) == 'txt') {
      print(sprintf('Originally it has %s cells and %s features originated from a txt format file to import into Seurat object', dim(cellrangerCounts)[2], dim(cellrangerCounts)[1] ))
    }
    ##--------------------------------------------------------------------------------------##
    ##--------------------------------------------------------------------------------------##
    seuratObjOrg               <- Seurat::CreateSeuratObject(counts = cellrangerCounts,  project = names(cellrangerResList)[x], min.cells = as.numeric(minCells), min.features = as.numeric(minFeatures))
    ## add metadata feature2 into object, here is 'expCond' for metadata$sample
    seuratObjOrg               <- Seurat::AddMetaData(object = seuratObjOrg,  col.name = 'expCond', metadata = as.factor(metadata$sample[x]))
    # print(head(seuratObjOrg@meta.data)) ##for debug
    ## extra metadata columns in metadata columns 'expCond*' will be added respectively
    for (i in 1:expCond.no ) {
      seuratObjOrg             <- Seurat::AddMetaData(object = seuratObjOrg,  col.name = sprintf("expCond%s",i), metadata = as.factor(metadata[x, expCond.index[i]]) )
    }
    # print(head(seuratObjOrg@meta.data)) ##for debug
    # print(table(seuratObjOrg@meta.data$expCond3)) ##for debug
    ## -------------------------
    # if ( 'expCond1' %in% colnames(metadata) & 'expCond2' %in% colnames(metadata) ) {
    #   seuratObjOrg             <- Seurat::AddMetaData(object = seuratObjOrg,  col.name = 'expCond1', metadata = as.factor(metadata$expCond1[x]))
    #   seuratObjOrg             <- Seurat::AddMetaData(object = seuratObjOrg,  col.name = 'expCond2', metadata = as.factor(metadata$expCond2[x]))
    # } else if ( 'expCond1' %in% colnames(metadata) & !'expCond2' %in% colnames(metadata) ) {
    #   seuratObjOrg             <- Seurat::AddMetaData(object = seuratObjOrg,  col.name = 'expCond1', metadata = as.factor(metadata$expCond1[x]))
    # } else if ( !'expCond1' %in% colnames(metadata) & 'expCond2' %in% colnames(metadata) ) {
    #   seuratObjOrg             <- Seurat::AddMetaData(object = seuratObjOrg,  col.name = 'expCond2', metadata = as.factor(metadata$expCond2[x]))
    # }
    ## -------------------------
    ##--------------------------------------------------------------------------------------##
    orgCellNoSummary           <- data.frame('cellrangeRcellNo' = dim(seuratObjOrg)[2], 'cellrangeRfeatureNo' = dim(seuratObjOrg)[1] )
    if (x == 1) {
      orgCellNoSummarySampComb <- orgCellNoSummary
    } else {
      orgCellNoSummarySampComb <- rbind(orgCellNoSummarySampComb, orgCellNoSummary)
    }
    ##--------------------------------------------------------------------------------------##
    doubletsMethod             <- tolower(doubletsRmMethods[x])
    if (doubletsMethod == 'none') {
      seuratObj                <- seuratObjOrg
      print('No doublet removal is executed here')
      print(seuratObjOrg)
      print('---')
    } else {
      if (is.na(doubletsResDirs[x])) stop('doubletsMethod is on, please provide corresponding full path to the saved doublets removal results')
      print(sprintf('%s doublet removal methods were implemented.', doubletsMethod))
      if (doubletsMethod == 'centroids') {
        doubletFname           <- paste(doubletsResDirs[x], '/', names(cellrangerResList), '_centroids_doublet_cells_name.Rdata', sep = '')
        load(doubletFname[x])
        print(sprintf('%s doublets were estimated from %s sample', length(centroidsDoubletCells), names(cellrangerResList)[x] ))
        doubletCellsUpdate <- gsub(pattern = '[.]', replacement = '-', centroidsDoubletCells)
      } else if (doubletsMethod == 'medoids') {
        doubletFname           <- paste(doubletsResDirs[x], '/', names(cellrangerResList), '_medoids_doublet_cells_name.Rdata', sep = '')
        load(doubletFname[x])
        print(sprintf('%s doublets were estimated from %s sample', length(medoidsDoubletCells), names(cellrangerResList)[x] ))
        doubletCellsUpdate <- gsub(pattern = '[.]', replacement = '-', medoidsDoubletCells)
      } else if (doubletsMethod == 'ol') {
        doubletFname           <- paste(doubletsResDirs[x], '/', names(cellrangerResList), '_OL_doublet_cells_name.Rdata', sep = '')
        load(doubletFname[x])
        print(sprintf('%s doublets were estimated from %s sample', length(olDoubletCells), names(cellrangerResList)[x] ))
        doubletCellsUpdate <- gsub(pattern = '[.]', replacement = '-', olDoubletCells)
      }
      doubletDf                <- data.frame(cells = rownames(seuratObjOrg@meta.data))
      doubletDf                <- doubletDf %>% dplyr::mutate(doublet = ifelse(rownames((seuratObjOrg@meta.data)) %in% doubletCellsUpdate, 'TRUE', 'FASLE') )
      seuratObjOrg             <- Seurat::AddMetaData(object = seuratObjOrg,  col.name = 'doublet', metadata = as.factor(doubletDf$doublet) )
      seuratObj                <- subset(seuratObjOrg, doublet == 'FASLE')
      print('Original Seurat object without doublet removal:')
      print(seuratObjOrg)
      print('---')
      print(sprintf('%s Doublets were removal, %s cells were left', length(doubletCellsUpdate), (dim(seuratObjOrg@meta.data)[1] - length(doubletCellsUpdate)) ))
      print('doublet removed object:')
      print(seuratObj)
      print('---')
    }
    seuratObjList[[x]]         <- seuratObj
  }
  print(sprintf('Step 1: END read in 10X data into Seurat object at %s', Sys.time()))
  print('---===---')
  ##--------------------------------------------------------------------------------------##
  names(seuratObjList)         <- names(cellrangerResList)
  cellNoSummary                <- data.frame('afterDoubletsRmCellNo' = unlist( lapply(seuratObjList, function(x) dim(x)[2])), 'afterDoubletsRmfeatureNo' = unlist(lapply(seuratObjList, function(x) dim(x)[1])) )
  print('---===---')
  print('Original cellrangR processed cell No. and features')
  print(orgCellNoSummarySampComb)
  print('---')
  print('To be processed cell No. and features')
  print(cellNoSummary)
  print('---===---')
  cellNoSummaryComb            <- cbind(orgCellNoSummarySampComb, cellNoSummary)
  rownames(cellNoSummaryComb)  <- names(cellrangerResList)
  write.table(x = cellNoSummaryComb, file = file.path(resDir, 'org_doubletsRemoval_cellNoSummary.txt'), quote = F, row.names = T, col.names = NA, sep = '\t')
  ##--------------------------------------------------------------------------------------##
  ## 2. pre-processing: 1) check/filter the mitochondrial content, 2) normalization, 3) find variable features, and scaling
  print('Step 2: check mitochondrial content & normalization')
  ## 2.1 create 'seuratObjListPlotsDir' for each item in processed 'seuratObjList' based on input 'cellrangerResList'
  qcPlotsDir                  <- paste(resDir, 'QC_plots', sep = '/')
  if (!dir.exists(qcPlotsDir)) dir.create(qcPlotsDir)
  ## ---
  seuratQcProcessObjList      <- list()
  for (x in 1:length(seuratObjList)) {
    ## ---------
    seuratObj                 <- seuratObjList[[x]]
    ## 1).1 calculating mitochondrial content
    # print(sprintf('%s mitochondrial genes are processed', length(grep('^MT-', seuratObj@assays$RNA@data@Dimnames[[1]])) ))
    # seuratObj[['percent.mt']] <- PercentageFeatureSet(object = seuratObj, pattern = '^MT-')
    print(sprintf('%s mitochondrial genes are processed', length(grep(mtPatten(as.character(genomeSpecies)), seuratObj@assays$RNA$data@Dimnames[[1]])) ))
    seuratObj[['percent.mt']] <- Seurat::PercentageFeatureSet(object = seuratObj, pattern = mtPatten(as.character(genomeSpecies)) )
    print(head(seuratObj@meta.data, 5))
    ## calculating no. of cells with certain mitochondrial percentage
    mtPer <- list()
    rRNAper <- list()
    for(k in seq_along(1:4)){
      mtPer[[k]] <- sum(seuratObj@meta.data$percent.mt < 5*k) / length(seuratObj@meta.data$percent.mt)
    }
    mtPer                     <- do.call(rbind, mtPer)
    mtPerTab                  <- data.frame(`MT cutoff` = c("<5%", "<10%", "<15%", "<20%"), `Percentage of cells` =  scales::label_percent()(mtPer[,1]))
    print(mtPerTab)
    if (x == 1) {
      mtPerTabSampComb        <- mtPerTab
    } else {
      mtPerTabSampComb        <- dplyr::full_join(x = mtPerTabSampComb, y = mtPerTab, by = 'MT.cutoff')
    }
    ## ---
    ## calculate rRNA content
    print(sprintf('%s rRNA genes are processed', length(grep(rRNAcontent(as.character(genomeSpecies)), seuratObj@assays$RNA$data@Dimnames[[1]])) ))
    seuratObj[['rRNA.content']] <- Seurat::PercentageFeatureSet(object = seuratObj, pattern = rRNAcontent(as.character(genomeSpecies)) )
    for(k in seq_along(1:20)){
      rRNAper[[k]] <- sum(seuratObj@meta.data$rRNA.content < 5*k) / length(seuratObj@meta.data$rRNA.content)
    }
    rRNAper                     <- do.call(rbind, rRNAper)
    rRNAperTab                  <- data.frame(`rRNA content` = c("<5%", "<10%", "<15%", "<20%",
                                                                 "<25%", "<30%", "<35%", "<40%",
                                                                 "<45%", "<50%", "<55%", "<60%",
                                                                 "<65%", "<70%", "<75%", "<80%",
                                                                 "<85%", "<90%", "<95%", "<100%"),
                                              `Percentage of cells` =  scales::label_percent()(rRNAper[,1]))
    if (x == 1) {
      rRNAperSampComb           <- rRNAperTab
    } else {
      rRNAperSampComb           <- dplyr::full_join(x = rRNAperSampComb, y = rRNAperTab, by = 'rRNA.content')
    }
    ## ------
    # print(rRNAperSampComb)
    # write.table(x = mtPerTab, file = paste(resDir, '/ToDelete_MT_percentage_summary_', names(seuratObjList)[x], '.txt', sep = ''), quote = F, sep = '\t', row.names = F, col.names = T)
    ## 1).2 plot mitochondrial content
    pdf(file = paste(qcPlotsDir, '/featureScatter_', names(seuratObjList)[x], '.pdf', sep = ''), width = 10, height = 6)
    # VlnPlot(object = seuratObj, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
    plot1                     <- Seurat::FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2                     <- Seurat::FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    print(plot1 + plot2)
    dev.off()
    ## 1). 3 visualize QC metrics as a violin plot
    pdf(file = paste(qcPlotsDir, '/featureViolin_', names(seuratObjList)[x], '.pdf', sep = ''), width = 10, height = 6)
    vlnFeaturePlot            <- Seurat::VlnPlot(seuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "rRNA.content"), ncol = 4)
    print(vlnFeaturePlot)
    dev.off()
    ## --
    if (mtFiltering) {
      ## Filtering is on, update 'seuratObj' with subsetted obj
      ## 1). 3 remove unwanted features and cells
      seuratObjBeforeFilter   <- seuratObj
      seuratObj               <- subset(seuratObjBeforeFilter, subset = nFeature_RNA > 200 & percent.mt < mtPerCutoff)
      ## 1).2 plot mitochondrial content
      pdf(file = paste(qcPlotsDir, '/featureScatter_', names(seuratObjList)[x], '_afterFiltering.pdf', sep = ''), width = 10, height = 6)
      # VlnPlot(object = seuratObj, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
      plot1 <- Seurat::FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "percent.mt")
      plot2 <- Seurat::FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
      print(plot1 + plot2)
      dev.off()
      pdf(file = paste(qcPlotsDir, '/featureViolin_', names(seuratObjList)[x], '_afterFiltering.pdf', sep = ''), width = 10, height = 6)
      vlnFeaturePlot <- Seurat::VlnPlot(seuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "rRNA.content"), ncol = 4)
      print(vlnFeaturePlot)
      dev.off()
      print(sprintf('Processing sample %s (%s), After filtering out high percentage contained mitochondrial genes, low and high gene feature expressed cells, cell number reduced from %s to %s.', x, names(seuratObjList)[x], dim(seuratObjBeforeFilter@meta.data)[1], dim(seuratObj@meta.data)[1] ))
      noFilteredCells           <- data.frame(`filtering` = c('before', 'after'), `No cells` = c(dim(seuratObjBeforeFilter@meta.data)[1], dim(seuratObj@meta.data)[1]) )
      
      # write.table(x = noFilteredCells, file = paste(resDir, '/ToDelete_No_mtFiltered_cells_summary_', names(seuratObjList)[x], '.txt', sep = ''), quote = F, sep = '\t', row.names = F, col.names = T)
      ## Complete filtering
    } else {
      noFilteredCells           <- data.frame(`filtering` = c('before'), `No cells` = c(dim(seuratObj@meta.data)[1]) )
    }
    if (x == 1) {
      noFilteredCellsSampComb <- noFilteredCells
    } else {
      noFilteredCellsSampComb <- dplyr::full_join(x = noFilteredCellsSampComb, y = noFilteredCells, by = 'filtering')
    }
    ## -
    ## 2). Normalization: 'normalization.method' & 'scale.factor' are default options
    seuratObj                   <- Seurat::NormalizeData(seuratObj, normalization.method = "LogNormalize", scale.factor = 10000)
    ## 3). Find variable features, by default top 2000
    ## by default, select top 2000 features, vst is also default method, other options are 'mean.var.plot(mvp)' and 'dispersion (disp)'
    seuratObj                   <- Seurat::FindVariableFeatures(seuratObj, selection.method = 'vst', nfeatures = nfeatures)
    topFeatures                 <- head(Seurat::VariableFeatures(seuratObj), n = 10)
    ## 3).2 plot top 100 variable features
    pdf(file = paste(qcPlotsDir, '/topVariableFeature_', names(seuratObjList)[x], '.pdf', sep = ''), width = 8, height = 6)
    plot1 <- Seurat::VariableFeaturePlot(seuratObj)
    plot2 <- Seurat::LabelPoints(plot = plot1, points = topFeatures, repel = TRUE)
    print(plot2)
    dev.off()
    seuratQcProcessObjList[[x]] <- seuratObj
  }
  names(seuratQcProcessObjList) <- names(seuratObjList)
  colnames(mtPerTabSampComb)    <- c(colnames(mtPerTabSampComb)[1], names(seuratObjList))
  write.table(x = mtPerTabSampComb, file = file.path(resDir, 'MT_percentage_summary.txt'), quote = F, sep = '\t', row.names = F, col.names = T)
  ## ---
  colnames(rRNAperSampComb)    <- c(colnames(rRNAperSampComb)[1], names(seuratObjList))
  write.table(x = rRNAperSampComb, file = file.path(resDir, 'rRNA_content_summary.txt'), quote = F, sep = '\t', row.names = F, col.names = T)
  ## ---
  colnames(noFilteredCellsSampComb)   <- c(colnames(noFilteredCellsSampComb)[1], names(seuratObjList))
  if (length(seuratQcProcessObjList)>1) {
    noFilteredCellsSampComb$Total       <- rowSums(noFilteredCellsSampComb[,-1])
  }
  write.table(x = noFilteredCellsSampComb, file = file.path(resDir, 'No_filtered_cells_summary.txt'), quote = F, sep = '\t', row.names = F, col.names = T)
  print('Step 2: END check mitochondrial content & normalization')
  print('---===---')
  ##--------------------------------------------------------------------------------------##
  return(list('countReadInOjb' = seuratObjList, 'qcProcessObj' = seuratQcProcessObjList, 'resDir' = resDir))
}
## ---
rRNAcontent <-  function(genomeSpecies) {
  if (genomeSpecies == 'human') return('^RP[SL]')
  if (genomeSpecies == 'mouse') return('^rp[sl]')
  if (genomeSpecies == 'human_mouse') return('^RP[SL]|^rp[sl]')
}
## ---
meata2list <- function(metadata) {
  cellrangerResList    <- list()
  for (i in 1:length(metadata$sample)) {
    cellrangerResList[[i]] <- as.character(metadata$path[i])
  }
  names(cellrangerResList) <- metadata$sample
  return(cellrangerResList)
}


##----------------------------------------------------------------------------------------

metadata <- read.csv(paste0("qc/qc_NS-DD-1s-DEC_metadata.txt"), sep='\t')

# shuffle samples
set.seed(123)
metadata <- metadata[sample(1:dim(metadata)[1], replace = F),]


qcRes <- processQC(metadata = metadata, 
          resDirName = paste0("qc/qc_NS-DD-1s-DEC_summary_FINAL"), 
          genomeSpecies = 'human', 
          minCells = 3, # defualt
          mtFiltering = T, 
          minFeatures = 100, 
          mtPerCutoff = 20,
          multiomics = T,
          extraFilter = F,
          )


lapply(names(qcRes$qcProcessObj), function(name) {
  file_path <- file.path(qcRes$resDir, paste0(name, ".rds"))
  saveRDS(qcRes$qcProcessObj[[name]], file = file_path)
})



