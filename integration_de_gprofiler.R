rm(list=ls())

# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("org.Mm.eg.db")
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources
library(org.Hs.eg.db)
# library(org.Mm.eg.db)
library(Seurat)
library(readxl)
library(writexl)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(stringr) # long text
library(openxlsx)
require(RColorBrewer)
require(ComplexHeatmap)
require(circlize)
require(digest)
require(cluster)
library(magick)
library(tidyverse)
library(patchwork) # wrap plots 
library(data.table) #fread


organism =  "hsapiens"

# names(table(Idents(rds)))

deResDir <- "integration_2/integration_2_leiden/"

# local
#### list of protein-coding genes ####
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
################# GO for DEGs (MAST) across major cell types, do "CD3_18h", "CD3_4h", "Ig_18h", "Ig_4h", "LPS_18h", "LPS_4h" vs "unstim_0h"################# 
#### Load DEGs for each contrast ----
excel.dir <- "results_wOrgClusterAnnotation_DEGs/stimuli_time_MAST/"
file.all <- list.files(paste0(deResDir,excel.dir))
files <- file.all[grepl("full.*SelClusters\\.xlsx$", file.all, ignore.case = TRUE)]
contrast <- gsub("expCondCompDeMarkers_|_full.*SelClusters\\.xlsx$", "", files)

# save each cell type to a tibble for each comparison 
for (i in 1:length(files)) {
  excel_file <- paste0(deResDir, excel.dir, files[i])
  sheet_names <- excel_sheets(excel_file)
  for (sheet_name in sheet_names) {
    assign(paste0(contrast[i], "_", sheet_name), read_excel(excel_file, sheet = sheet_name))
  }
}

#### function getUpSetInput ----
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


#### function go_up_down ----
go_up_down <- function(group = group, 
                       ct = ct,
                       tp = tp){
  # look into the top 10 most significant GOs per source
  go <- read_excel(paste0(deResDir,excel.dir, 'GOEnrichment/', group, "_", ct,'_GO_table.xlsx') , sheet = 'Sheet1')
  
  go_top10_GO_KEGG_TF <- go %>%
    filter(source %in% c('GO:MF', 'GO:CC', 'GO:BP', 'KEGG', 'TF')) %>%
    filter(grepl(tp, query)) %>%
    group_by(query, source) %>%
    arrange(p_value) %>%
    slice_head(n=10) %>%
    ungroup()
  
  if (any(go_top10_GO_KEGG_TF$significant == FALSE)) {
    print("There is at least one record where 'significant.' is FALSE.")
  } else {
    print("All records are significant.")
  }
  
  go_top10_GO_KEGG_TF$source <- factor(go_top10_GO_KEGG_TF$source,
                                       levels=c('GO:MF', 'GO:CC', 'GO:BP', 'KEGG', 'TF'))
  
  
  ylim <- ceiling(-log10(min(go_top10_GO_KEGG_TF$p_value)))
  # barchart
  p1 <- go_top10_GO_KEGG_TF %>%
    filter(query %in% paste0('unique_up_',group, "_", tp)) %>%
    ggplot(aes(x=factor(term_name, levels = term_name[order(-p_value)]), y=-log10(p_value)))+
    geom_bar(stat="identity", width=0.7, fill='#e63946')+
    theme_bw() +
    # theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
    # scale_y_reverse(limits = c(40, 0),breaks = seq(40, 0, by = -10)) +
    scale_y_reverse(limits = c(ylim, 0)) +
    # scale_y_continuous(limits = c(0, 150),breaks = seq(0, 150, by = 50)) +
    labs(x = NULL) +
    facet_grid(source~., scales="free", space="free") +
    # theme(strip.placement="outside") +
    coord_flip() +
    theme(strip.text.y = element_blank() ,
          strip.background = element_blank())
  
  # lollipop chart
  p1 <- go_top10_GO_KEGG_TF %>%
    filter(query %in% paste0('unique_up_',group, "_", tp)) %>%
    ggplot(aes(x=factor(term_name, levels = term_name[order(-p_value)]), y=-log10(p_value)))+
    geom_segment(aes(xend = factor(term_name, levels = term_name[order(-p_value)]), 
                     yend = 0), 
                 color = "#e63946", linewidth = 0.5) +  # Add the line segment
    geom_point(color = "#e63946", size = 2) +  # Add the point at the end of the line
    theme_bw() +
    # theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
    # scale_y_reverse(limits = c(40, 0),breaks = seq(40, 0, by = -10)) +
    scale_y_reverse(limits = c(ylim, 0)) +
    # scale_y_continuous(limits = c(0, 150),breaks = seq(0, 150, by = 50)) +
    labs(x = NULL) +
    facet_grid(source~., scales="free", space="free") +
    # theme(strip.placement="outside") +
    coord_flip() +
    theme(strip.text.y = element_blank() ,
          strip.background = element_blank()) 
  # + scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) # Place y-axis tick labels on the rCD3ht
  
  # barchart
  p2 <- go_top10_GO_KEGG_TF %>%
    filter(query %in% paste0('unique_down_',group, "_", tp)) %>%
    ggplot(aes(x=factor(term_name, levels = term_name[order(-p_value)]), y=-log10(p_value)))+
    geom_bar(stat="identity", width=0.7, fill='#a8dadc')+
    theme_bw() +
    # theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
    # scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
    scale_y_continuous(limits = c(0, ylim)) + 
    labs(x = "Term name") +
    facet_grid(source~., scales="free", space="free") +
    theme(strip.placement="outside") +
    coord_flip() +
    scale_x_discrete(position = "top")  # Place y-axis tick labels on the rCD3ht
  
  # lollipop chart
  p2 <- go_top10_GO_KEGG_TF %>%
    filter(query %in% paste0('unique_down_',group, "_", tp)) %>%
    ggplot(aes(x = factor(term_name, levels = term_name[order(-p_value)]), 
               y = -log10(p_value))) +
    geom_segment(aes(xend = factor(term_name, levels = term_name[order(-p_value)]), 
                     yend = 0), 
                 color = "#a8dadc", linewidth = 0.5) +  # Add the line segment
    geom_point(color = "#a8dadc", size = 2) +  # Add the point at the end of the line
    theme_bw() +
    labs(x = NULL) +
    scale_y_continuous(limits = c(0, ylim)) + 
    facet_grid(source ~ ., scales = "free", space = "free") +
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = "white", color = "black")) +
    coord_flip() +
    scale_x_discrete(#labels = function(x) str_wrap(x, width = 40),
      position = "top") # Place y-axis tick labels on the rCD3ht
  
  return(list(p1 = p1, p2 = p2, go_top10_GO_KEGG_TF= go_top10_GO_KEGG_TF))
  
}

go_up_down_db <- function(group = group, 
                          ct = ct,
                          conditions = conditions){
  go <- read_excel(paste0(deResDir,excel.dir, 'GOEnrichment/', group, "_", ct,'_GO_table.xlsx') , sheet = 'Sheet1')
  
  if (any(go$significant == FALSE)) {
    print("There is at least one record where 'significant.' is FALSE.")
  } else {
    print("All records are significant.")
  }
  
  # Loop over each condition and get top 10 terms per database (GO:MF, GO:CC, GO:BP, KEGG, TF)
  go_top10_GO_KEGG_TF <- list()
  
  go_top4_GO_KEGG_TF <- list()
  
  for (condition in conditions){
    # Filter and select top 10 terms for each database
    top_terms_10 <- go %>%
      filter(source %in% c('GO:MF', 'GO:CC', 'GO:BP', 'KEGG', 'TF')) %>%
      filter(grepl(condition, query)) %>%
      group_by(query, source) %>%
      arrange(p_value) %>%
      slice_head(n = 10) %>%
      ungroup()
    
    # Filter and select top 4 terms for each database
    top_terms_4 <- go %>%
      filter(source %in% c('GO:MF', 'GO:CC', 'GO:BP', 'KEGG', 'TF')) %>%
      filter(grepl(condition, query)) %>%
      group_by(query, source) %>%
      arrange(p_value) %>%
      slice_head(n = 4) %>%
      ungroup()
    
    # top_terms$source <- factor(top_terms$source,
    #                            levels=c('GO:MF', 'GO:CC', 'GO:BP', 'KEGG', 'TF'))
    
    # Store results in the list, each entry named by the condition
    go_top10_GO_KEGG_TF[[condition]] <- top_terms_10
    go_top4_GO_KEGG_TF[[condition]] <- top_terms_4
    
  }
  
  combined_go_top10_GO_KEGG_TF <- do.call(rbind, go_top10_GO_KEGG_TF)
  combined_go_top4_GO_KEGG_TF <- do.call(rbind, go_top4_GO_KEGG_TF)
  
  ylim <- ceiling(-log10(min(combined_go_top10_GO_KEGG_TF$p_value)))
  
  plots <-  list()
  
  color_values <- c(rep("#FE4088", 3), rep("#5F96FA", 3))

  names(color_values) <- conditions
  
  color_labels <- c(paste0("Shared up-regulated ", group),
                    paste0("Unique up-regulated ", group, "_4h"),
                    paste0("Unique up-regulated ", group, "_18h"),
                    paste0("Shared down-regulated ", group),
                    paste0("Unique down-regulated ", group, "_4h"),
                    paste0("Unique down-regulated ", group, "_18h"))
  
  for (db in unique(combined_go_top10_GO_KEGG_TF$source)){
    db_name<- str_replace_all(db, ":", "_")
    p <- combined_go_top10_GO_KEGG_TF %>%
      filter(grepl(db, source)) %>%
      group_by(query) %>%
      arrange(query, desc(p_value), .by_group = TRUE) %>%
      # Create a unique ordering for term_name within each query
      mutate(term_name_order = factor(term_name, levels = unique(term_name[order(-p_value)]))) %>%
      mutate(query = factor(query, levels = conditions)) %>%
      ggplot(aes(x = term_name_order, y = -log10(p_value))) +
      geom_segment(aes(xend = term_name_order, yend = 0, color = query), linewidth = 0.5) +  # Add the line segment
      geom_point(aes(color = query), size = 2) +  # Add the point at the end of the line
      scale_color_manual(values = color_values,
                         labels = color_labels) +
      theme_bw() +
      coord_flip() +  # Flip coordinates
      scale_y_reverse(limits = c(ylim, 0)) +
      labs(x = NULL, color = "Genes") +
      facet_grid(query ~ source, scales = "free", space = "free") +
      theme(strip.text.y = element_blank() ,
            strip.background = element_blank()) 
    # +
      # scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 40))
    
    # Save the plot to a PDF file
    ggsave(filename = paste0(deResDir, excel.dir, 'GOEnrichment/', group, "_", ct,"_", db_name, ".pdf"), plot = p, width = 8, height = 7)
    
    # Add the plot to the list
    plots[[db]] <- p
  }
  
  return(list(plots = plots, 
              combined_go_top10_GO_KEGG_TF = combined_go_top10_GO_KEGG_TF,
              combined_go_top4_GO_KEGG_TF = combined_go_top4_GO_KEGG_TF))
  
}

# --
#### B Ig ####
# Define the directory path

go_dir <- "GOEnrichment"

# Check if the directory exists
if (!dir.exists(paste0(deResDir, excel.dir, go_dir))) {
  # If it doesn't exist, create the directory
  dir.create(paste0(deResDir,excel.dir, go_dir))
  cat("Directory 'GOEnrichment' created.\n")
} else {
  cat("Directory 'GOEnrichment' already exists.\n")
}

celltype <- "B"
ct <- gsub("/| |-","",celltype) # replace any / space -
ct <- paste0("cluster_", ct)
upsetInput <- getUpSetInput(ct)
names(upsetInput$up) <- lapply(names(upsetInput$up), function(x){ return(paste0("up_",x)) })
names(upsetInput$dn) <- lapply(names(upsetInput$dn), function(x){ return(paste0("down_",x)) })
upsetInput.merge <- c(upsetInput$up, upsetInput$dn)

# for each treatment, do GO per cell type
for (group in c("Ig")) { # search the comparison that includes ig in upsetInput.merge
  upsetInput.group <- upsetInput.merge[grep(group, names(upsetInput.merge), value = TRUE)]
  
  # Extract the individual gene sets
  down_Ig_4h <- upsetInput.group[["down_Ig_4h/unstim_0h"]]
  down_Ig_18h <- upsetInput.group[["down_Ig_18h/unstim_0h"]]
  up_Ig_4h <- upsetInput.group[["up_Ig_4h/unstim_0h"]]
  up_Ig_18h <- upsetInput.group[["up_Ig_18h/unstim_0h"]]
  
  unique_down_Ig_4h <- setdiff(down_Ig_4h, down_Ig_18h)
  unique_down_Ig_18h <- setdiff(down_Ig_18h, down_Ig_4h)
  unique_up_Ig_4h <- setdiff(up_Ig_4h, up_Ig_18h)
  unique_up_Ig_18h <- setdiff(up_Ig_18h, up_Ig_4h)
  
  # Find overlapping gene sets
  overlap_down_Ig <- intersect(down_Ig_4h, down_Ig_18h)
  overlap_up_Ig <- intersect(up_Ig_4h, up_Ig_18h)
  
  # Add these unique sets back to upsetInput.group
  upsetInput.group[["unique_down_Ig_4h"]] <- unique_down_Ig_4h
  upsetInput.group[["unique_down_Ig_18h"]] <- unique_down_Ig_18h
  upsetInput.group[["unique_up_Ig_4h"]] <- unique_up_Ig_4h
  upsetInput.group[["unique_up_Ig_18h"]] <- unique_up_Ig_18h
  upsetInput.group[["overlap_down_Ig"]] <- overlap_down_Ig
  upsetInput.group[["overlap_up_Ig"]] <- overlap_up_Ig
  

  upsetInput.group.cnts <- data.frame(
    treatment = c("unique_up_Ig_4h", 
                  "unique_up_Ig_18h", 
                  "unique_down_Ig_4h", 
                  "unique_down_Ig_18h",
                  "overlap_down_Ig",
                  "overlap_up_Ig"),
    de_gene_cnts = c(length(unique_up_Ig_4h), 
                     length(unique_up_Ig_18h), 
                     length(unique_down_Ig_4h), 
                     length(unique_down_Ig_18h),
                     length(overlap_down_Ig),
                     length(overlap_up_Ig)))
  
  wb <- createWorkbook()
  addWorksheet(wb, "de_gene_cnts")
  writeData(wb, sheet = "de_gene_cnts", upsetInput.group.cnts)
  
  addWorksheet(wb, "unique_up_Ig_4h")
  writeData(wb, sheet = "unique_up_Ig_4h", data.frame(unique_up_Ig_4h))
  
  addWorksheet(wb, "overlap_up_Ig")
  writeData(wb, sheet = "overlap_up_Ig", data.frame(overlap_up_Ig))
  
  addWorksheet(wb, "unique_up_Ig_18h")
  writeData(wb, sheet = "unique_up_Ig_18h", data.frame(unique_up_Ig_18h))
  
  addWorksheet(wb, "unique_down_Ig_4h")
  writeData(wb, sheet = "unique_down_Ig_4h", data.frame(unique_down_Ig_4h))
  
  addWorksheet(wb, "unique_down_Ig_18h")
  writeData(wb, sheet = "unique_down_Ig_18h", data.frame(unique_down_Ig_18h))
  
  addWorksheet(wb, "overlap_down_Ig")
  writeData(wb, sheet = "overlap_down_Ig", data.frame(overlap_down_Ig))
  
  saveWorkbook(wb, file = paste0(deResDir, excel.dir, 'GOEnrichment/GO_', group, "_", ct, '_de_genes.xlsx'), 
               overwrite = TRUE)
  
  # write_xlsx(upsetInput.group.cnts, path = paste0(deResDir,excel.dir,'GOEnrichment/', group, "_", ct,'_de_genes_non_overlapping_cnts.xlsx'), col_names = T)
  
  # gost.res <- gost(upsetInput.group[grep("unique", names(upsetInput.group), value = TRUE)], organism = organism, correction_method = "fdr")
  
  gene_sets <- c(
    upsetInput.group[grep("unique", names(upsetInput.group), value = TRUE)],
    upsetInput.group[grep("overlap", names(upsetInput.group), value = TRUE)]
  )
  
  gost.res <- gost(gene_sets, organism = organism, correction_method = "fdr", evcodes=FALSE)
  
  gost.res$result$query <- factor(gost.res$result$query, 
                                  levels = c("overlap_up_Ig", "unique_up_Ig_4h", "unique_up_Ig_18h", 
                                             "overlap_down_Ig", "unique_down_Ig_4h", "unique_down_Ig_18h"))
  
  
  mygostplot <- gostplot(gost.res, interactive = F, capped = F)
  ggsave(paste0(deResDir,excel.dir,"GOEnrichment/", group, "_", ct, ".pdf"), #Ig_Bcells_.png
         mygostplot, height = 12, width = 8, units = "in")
  write_xlsx(as.data.frame(gost.res$result), path = paste0(deResDir,excel.dir,'GOEnrichment/', group, "_", ct,'_GO_table.xlsx'), col_names = T)
}

# # skip e.g. Ig_up_4h vs Ig_down_4h across 4 dbs on the same figure ----
# group <- "Ig"
# 
# # go %>%
# #   filter(source %in% c('GO:MF', 'GO:CC', 'GO:BP', 'KEGG', 'TF')) %>%
# #   filter(grepl(tp, query)) %>%
# #    filter(source %in% "KEGG") %>%
# #    print(n=40)
#  
# # for (tp in c("4h", "18h")){
# #   go_up_down(group = group, 
# #              ct = ct,
# #              tp = tp)
# # }
# 
# tp <- "4h"
# p1 <- go_up_down(group = group, 
#                  ct = ct,
#                  tp = tp)$p1
# p2 <- go_up_down(group = group, 
#                  ct = ct,
#                  tp = tp)$p2
# pdf(paste0(deResDir,excel.dir, 'GOEnrichment/', group, "_", tp, "_", ct,"_GO_KEGG_TF.pdf"),
#     width = 12,
#     height = 8)
# grid.arrange(p1, p2, ncol=2, widths=c(1,1))
# dev.off()
# 
# 
# 
# 
# tp <- "18h"
# p1 <- go_up_down(group = group, 
#                  ct = ct,
#                  tp = tp)$p1
# p2 <- go_up_down(group = group, 
#                  ct = ct,
#                  tp = tp)$p2
# pdf(paste0(deResDir,excel.dir, 'GOEnrichment/', group, "_", tp, "_", ct,"_GO_KEGG_TF.pdf"),
#     width = 12,
#     height = 8)
# grid.arrange(p1, p2, ncol=2, widths=c(1,1))
# dev.off()
# 
# 
#e.g. shared and unique up and down (6 conditions) in one db on the same figure ----
group <- "Ig"
# 6 conditions across 4 databases
conditions <- c("overlap_up_Ig", "unique_up_Ig_4h", "unique_up_Ig_18h", 
                "overlap_down_Ig", "unique_down_Ig_4h", "unique_down_Ig_18h")


plots <- go_up_down_db(group = group, ct = ct, conditions)$plots
combined_go_top10_GO_KEGG_TF1 <- go_up_down_db(group = group, ct = ct, conditions)$combined_go_top10_GO_KEGG_TF
# plots[[5]] +  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
#   theme(axis.text.y = element_text(size = 7))

combined_go_top4_GO_KEGG_TF1 <- go_up_down_db(group = group, ct = ct, conditions)$combined_go_top4_GO_KEGG_TF


#### CD4 T CD3  ####

# Define the directory path
go_dir <- "GOEnrichment"

# Check if the directory exists
if (!dir.exists(paste0(deResDir,excel.dir, go_dir))) {
  # If it doesn't exist, create the directory
  dir.create(paste0(deResDir,excel.dir, go_dir))
  cat("Directory 'GOEnrichment' created.\n")
} else {
  cat("Directory 'GOEnrichment' already exists.\n")
}

celltype <- "CD4 T"
ct <- gsub("/| |-","",celltype) # replace any / space -
ct <- paste0("cluster_", ct)
upsetInput <- getUpSetInput(ct)
names(upsetInput$up) <- lapply(names(upsetInput$up), function(x){ return(paste0("up_",x)) })
names(upsetInput$dn) <- lapply(names(upsetInput$dn), function(x){ return(paste0("down_",x)) })
upsetInput.merge <- c(upsetInput$up, upsetInput$dn)

# for each treatment, do GO per cell type
for (group in c("CD3")) { # search the comparison that includes CD3 in upsetInput.merge
  upsetInput.group <- upsetInput.merge[grep(group, names(upsetInput.merge), value = TRUE)]
  
  # Extract the individual gene sets
  down_CD3_4h <- upsetInput.group[["down_CD3_4h/unstim_0h"]]
  down_CD3_18h <- upsetInput.group[["down_CD3_18h/unstim_0h"]]
  up_CD3_4h <- upsetInput.group[["up_CD3_4h/unstim_0h"]]
  up_CD3_18h <- upsetInput.group[["up_CD3_18h/unstim_0h"]]
  
  unique_down_CD3_4h <- setdiff(down_CD3_4h, down_CD3_18h)
  unique_down_CD3_18h <- setdiff(down_CD3_18h, down_CD3_4h)
  unique_up_CD3_4h <- setdiff(up_CD3_4h, up_CD3_18h)
  unique_up_CD3_18h <- setdiff(up_CD3_18h, up_CD3_4h)
  
  # Find overlapping gene sets
  overlap_down_CD3 <- intersect(down_CD3_4h, down_CD3_18h)
  overlap_up_CD3 <- intersect(up_CD3_4h, up_CD3_18h)
  
  # Add these unique sets back to upsetInput.group
  upsetInput.group[["unique_down_CD3_4h"]] <- unique_down_CD3_4h
  upsetInput.group[["unique_down_CD3_18h"]] <- unique_down_CD3_18h
  upsetInput.group[["unique_up_CD3_4h"]] <- unique_up_CD3_4h
  upsetInput.group[["unique_up_CD3_18h"]] <- unique_up_CD3_18h
  upsetInput.group[["overlap_down_CD3"]] <- overlap_down_CD3
  upsetInput.group[["overlap_up_CD3"]] <- overlap_up_CD3
  
  
  upsetInput.group.cnts <- data.frame(
    treatment = c("unique_up_CD3_4h", 
                  "unique_up_CD3_18h", 
                  "unique_down_CD3_4h", 
                  "unique_down_CD3_18h",
                  "overlap_down_CD3",
                  "overlap_up_CD3"),
    de_gene_cnts = c(length(unique_up_CD3_4h), 
                     length(unique_up_CD3_18h), 
                     length(unique_down_CD3_4h), 
                     length(unique_down_CD3_18h),
                     length(overlap_down_CD3),
                     length(overlap_up_CD3)))
  
  wb <- createWorkbook()
  addWorksheet(wb, "de_gene_cnts")
  writeData(wb, sheet = "de_gene_cnts", upsetInput.group.cnts)
  
  addWorksheet(wb, "unique_up_CD3_4h")
  writeData(wb, sheet = "unique_up_CD3_4h", data.frame(unique_up_CD3_4h))
  
  addWorksheet(wb, "overlap_up_CD3")
  writeData(wb, sheet = "overlap_up_CD3", data.frame(overlap_up_CD3))
  
  addWorksheet(wb, "unique_up_CD3_18h")
  writeData(wb, sheet = "unique_up_CD3_18h", data.frame(unique_up_CD3_18h))
  
  addWorksheet(wb, "unique_down_CD3_4h")
  writeData(wb, sheet = "unique_down_CD3_4h", data.frame(unique_down_CD3_4h))
  
  addWorksheet(wb, "unique_down_CD3_18h")
  writeData(wb, sheet = "unique_down_CD3_18h", data.frame(unique_down_CD3_18h))
  
  addWorksheet(wb, "overlap_down_CD3")
  writeData(wb, sheet = "overlap_down_CD3", data.frame(overlap_down_CD3))
  
  saveWorkbook(wb, file = paste0(deResDir, excel.dir, 'GOEnrichment/GO_', group, "_", ct, '_de_genes.xlsx'), overwrite = TRUE)
  
  # write_xlsx(upsetInput.group.cnts, path = paste0(deResDir,excel.dir,'GOEnrichment/', group, "_", ct,'_de_genes_non_overlapping_cnts.xlsx'), col_names = T)
  
  # gost.res <- gost(upsetInput.group[grep("unique", names(upsetInput.group), value = TRUE)], organism = organism, correction_method = "fdr")
  
  gene_sets <- c(
    upsetInput.group[grep("unique", names(upsetInput.group), value = TRUE)],
    upsetInput.group[grep("overlap", names(upsetInput.group), value = TRUE)]
  )
  
  gost.res <- gost(gene_sets, organism = organism, correction_method = "fdr", evcodes=FALSE)
  
  gost.res$result$query <- factor(gost.res$result$query, 
                                  levels = c("overlap_up_CD3", "unique_up_CD3_4h", "unique_up_CD3_18h", 
                                             "overlap_down_CD3", "unique_down_CD3_4h", "unique_down_CD3_18h"))
  
  
  mygostplot <- gostplot(gost.res, interactive = F, capped = F)
  ggsave(paste0(deResDir,excel.dir,"GOEnrichment/", group, "_", ct, ".pdf"), #CD3_Bcells_.png
         mygostplot, height = 12, width = 8, units = "in")
  write_xlsx(as.data.frame(gost.res$result), path = paste0(deResDir,excel.dir,'GOEnrichment/', group, "_", ct,'_GO_table.xlsx'), col_names = T)
}



# # skip e.g. CD3_up_4h vs CD3_down_4h across 4 dbs on the same figure----
# group <- "CD3"
# 
# # go %>%
# #   filter(source %in% c('GO:MF', 'GO:CC', 'GO:BP', 'KEGG', 'TF')) %>%
# #   filter(grepl(tp, query)) %>%
# #    filter(source %in% "KEGG") %>%
# #    print(n=40)
# 
# 
# tp <- "4h"
# p1 <- go_up_down(group = group, 
#                  ct = ct,
#                  tp = tp)$p1
# p2 <- go_up_down(group = group, 
#                  ct = ct,
#                  tp = tp)$p2
# pdf(paste0(deResDir,excel.dir, 'GOEnrichment/', group, "_", tp, "_", ct,"_GO_KEGG_TF.pdf"),
#     width = 12,
#     height = 8)
# grid.arrange(p1, p2, ncol=2, widths=c(0.92,1))
# dev.off()
# 
# 
# 
# 
# tp <- "18h"
# p1 <- go_up_down(group = group, 
#                  ct = ct,
#                  tp = tp)$p1
# p2 <- go_up_down(group = group, 
#                  ct = ct,
#                  tp = tp)$p2
# pdf(paste0(deResDir,excel.dir, 'GOEnrichment/', group, "_", tp, "_", ct,"_GO_KEGG_TF.pdf"),
#     width = 12,
#     height = 8)
# grid.arrange(p1, p2, ncol=2, widths=c(1,0.88))
# 
# 
# dev.off()
# 
#e.g. shared and unique up and down (6 conditions) in one db on the same figure ----
group <- "CD3"
# 6 conditions across 4 databases
conditions <- c("overlap_up_CD3", "unique_up_CD3_4h", "unique_up_CD3_18h", 
                "overlap_down_CD3", "unique_down_CD3_4h", "unique_down_CD3_18h")


plots <- go_up_down_db(group = group, ct = ct, conditions)$plots
combined_go_top10_GO_KEGG_TF2 <- go_up_down_db(group = group, ct = ct, conditions)$combined_go_top10_GO_KEGG_TF
combined_go_top4_GO_KEGG_TF2 <- go_up_down_db(group = group, ct = ct, conditions)$combined_go_top4_GO_KEGG_TF

#### CD8 T CD3  ####

# Define the directory path
go_dir <- "GOEnrichment"

# Check if the directory exists
if (!dir.exists(paste0(deResDir,excel.dir, go_dir))) {
  # If it doesn't exist, create the directory
  dir.create(paste0(deResDir,excel.dir, go_dir))
  cat("Directory 'GOEnrichment' created.\n")
} else {
  cat("Directory 'GOEnrichment' already exists.\n")
}

celltype <- "CD8 T"
ct <- gsub("/| |-","",celltype) # replace any / space -
ct <- paste0("cluster_", ct)
upsetInput <- getUpSetInput(ct)
names(upsetInput$up) <- lapply(names(upsetInput$up), function(x){ return(paste0("up_",x)) })
names(upsetInput$dn) <- lapply(names(upsetInput$dn), function(x){ return(paste0("down_",x)) })
upsetInput.merge <- c(upsetInput$up, upsetInput$dn)

# for each treatment, do GO per cell type
for (group in c("CD3")) { # search the comparison that includes CD3 in upsetInput.merge
  upsetInput.group <- upsetInput.merge[grep(group, names(upsetInput.merge), value = TRUE)]
  
  # Extract the individual gene sets
  down_CD3_4h <- upsetInput.group[["down_CD3_4h/unstim_0h"]]
  down_CD3_18h <- upsetInput.group[["down_CD3_18h/unstim_0h"]]
  up_CD3_4h <- upsetInput.group[["up_CD3_4h/unstim_0h"]]
  up_CD3_18h <- upsetInput.group[["up_CD3_18h/unstim_0h"]]
  
  unique_down_CD3_4h <- setdiff(down_CD3_4h, down_CD3_18h)
  unique_down_CD3_18h <- setdiff(down_CD3_18h, down_CD3_4h)
  unique_up_CD3_4h <- setdiff(up_CD3_4h, up_CD3_18h)
  unique_up_CD3_18h <- setdiff(up_CD3_18h, up_CD3_4h)
  
  # Find overlapping gene sets
  overlap_down_CD3 <- intersect(down_CD3_4h, down_CD3_18h)
  overlap_up_CD3 <- intersect(up_CD3_4h, up_CD3_18h)
  
  # Add these unique sets back to upsetInput.group
  upsetInput.group[["unique_down_CD3_4h"]] <- unique_down_CD3_4h
  upsetInput.group[["unique_down_CD3_18h"]] <- unique_down_CD3_18h
  upsetInput.group[["unique_up_CD3_4h"]] <- unique_up_CD3_4h
  upsetInput.group[["unique_up_CD3_18h"]] <- unique_up_CD3_18h
  upsetInput.group[["overlap_down_CD3"]] <- overlap_down_CD3
  upsetInput.group[["overlap_up_CD3"]] <- overlap_up_CD3
  
  
  upsetInput.group.cnts <- data.frame(
    treatment = c("unique_up_CD3_4h", 
                  "unique_up_CD3_18h", 
                  "unique_down_CD3_4h", 
                  "unique_down_CD3_18h",
                  "overlap_down_CD3",
                  "overlap_up_CD3"),
    de_gene_cnts = c(length(unique_up_CD3_4h), 
                     length(unique_up_CD3_18h), 
                     length(unique_down_CD3_4h), 
                     length(unique_down_CD3_18h),
                     length(overlap_down_CD3),
                     length(overlap_up_CD3)))
  
  wb <- createWorkbook()
  addWorksheet(wb, "de_gene_cnts")
  writeData(wb, sheet = "de_gene_cnts", upsetInput.group.cnts)
  
  addWorksheet(wb, "unique_up_CD3_4h")
  writeData(wb, sheet = "unique_up_CD3_4h", data.frame(unique_up_CD3_4h))
  
  addWorksheet(wb, "overlap_up_CD3")
  writeData(wb, sheet = "overlap_up_CD3", data.frame(overlap_up_CD3))
  
  addWorksheet(wb, "unique_up_CD3_18h")
  writeData(wb, sheet = "unique_up_CD3_18h", data.frame(unique_up_CD3_18h))
  
  addWorksheet(wb, "unique_down_CD3_4h")
  writeData(wb, sheet = "unique_down_CD3_4h", data.frame(unique_down_CD3_4h))
  
  addWorksheet(wb, "unique_down_CD3_18h")
  writeData(wb, sheet = "unique_down_CD3_18h", data.frame(unique_down_CD3_18h))
  
  addWorksheet(wb, "overlap_down_CD3")
  writeData(wb, sheet = "overlap_down_CD3", data.frame(overlap_down_CD3))
  
  saveWorkbook(wb, file = paste0(deResDir, excel.dir, 'GOEnrichment/GO_', group, "_", ct, '_de_genes.xlsx'), overwrite = TRUE)
  
  # write_xlsx(upsetInput.group.cnts, path = paste0(deResDir,excel.dir,'GOEnrichment/', group, "_", ct,'_de_genes_non_overlapping_cnts.xlsx'), col_names = T)
  
  # gost.res <- gost(upsetInput.group[grep("unique", names(upsetInput.group), value = TRUE)], organism = organism, correction_method = "fdr")
  
  gene_sets <- c(
    upsetInput.group[grep("unique", names(upsetInput.group), value = TRUE)],
    upsetInput.group[grep("overlap", names(upsetInput.group), value = TRUE)]
  )
  
  gost.res <- gost(gene_sets, organism = organism, correction_method = "fdr", evcodes=FALSE)
  
  gost.res$result$query <- factor(gost.res$result$query, 
                                  levels = c("overlap_up_CD3", "unique_up_CD3_4h", "unique_up_CD3_18h", 
                                             "overlap_down_CD3", "unique_down_CD3_4h", "unique_down_CD3_18h"))
  
  
  mygostplot <- gostplot(gost.res, interactive = F, capped = F)
  ggsave(paste0(deResDir,excel.dir,"GOEnrichment/", group, "_", ct, ".pdf"), #CD3_Bcells_.png
         mygostplot, height = 12, width = 8, units = "in")
  write_xlsx(as.data.frame(gost.res$result), path = paste0(deResDir,excel.dir,'GOEnrichment/', group, "_", ct,'_GO_table.xlsx'), col_names = T)
}


# # skip e.g. CD3_up_4h vs CD3_down_4h across 4 dbs on the same figure----
# 
# group <- "CD3"
# 
# # go %>%
# #   filter(source %in% c('GO:MF', 'GO:CC', 'GO:BP', 'KEGG', 'TF')) %>%
# #   filter(grepl(tp, query)) %>%
# #    filter(source %in% "KEGG") %>%
# #    print(n=40)
# 
# 
# tp <- "4h"
# p1 <- go_up_down(group = group, 
#                  ct = ct,
#                  tp = tp)$p1
# p2 <- go_up_down(group = group, 
#                  ct = ct,
#                  tp = tp)$p2
# pdf(paste0(deResDir,excel.dir, 'GOEnrichment/', group, "_", tp, "_", ct,"_GO_KEGG_TF.pdf"),
#     width = 12,
#     height = 8)
# grid.arrange(p1, p2, ncol=2, widths=c(0.92,1))
# dev.off()
# 
# 
# 
# 
# tp <- "18h"
# p1 <- go_up_down(group = group, 
#                  ct = ct,
#                  tp = tp)$p1
# p2 <- go_up_down(group = group, 
#                  ct = ct,
#                  tp = tp)$p2
# pdf(paste0(deResDir,excel.dir, 'GOEnrichment/', group, "_", tp, "_", ct,"_GO_KEGG_TF.pdf"),
#     width = 12,
#     height = 8)
# grid.arrange(p1, p2, ncol=2, widths=c(1,0.92))
# 
# 
# dev.off()
# 
#e.g. shared and unique up and down (6 conditions) in one db on the same figure ----
group <- "CD3"
# 6 conditions across 4 databases
conditions <- c("overlap_up_CD3", "unique_up_CD3_4h", "unique_up_CD3_18h", 
                "overlap_down_CD3", "unique_down_CD3_4h", "unique_down_CD3_18h")


plots <- go_up_down_db(group = group, ct = ct, conditions)$plots
combined_go_top10_GO_KEGG_TF3 <- go_up_down_db(group = group, ct = ct, conditions)$combined_go_top10_GO_KEGG_TF
combined_go_top4_GO_KEGG_TF3 <- go_up_down_db(group = group, ct = ct, conditions)$combined_go_top4_GO_KEGG_TF

#### Mono/Mph LPS  ####

# Define the directory path
go_dir <- "GOEnrichment"

# Check if the directory exists
if (!dir.exists(paste0(deResDir,excel.dir, go_dir))) {
  # If it doesn't exist, create the directory
  dir.create(paste0(deResDir,excel.dir, go_dir))
  cat("Directory 'GOEnrichment' created.\n")
} else {
  cat("Directory 'GOEnrichment' already exists.\n")
}

celltype <- "Mono/Mph"
ct <- gsub("/| |-","",celltype) # replace any / space -
ct <- paste0("cluster_", ct)

upsetInput <- getUpSetInput(ct)
names(upsetInput$up) <- lapply(names(upsetInput$up), function(x){ return(paste0("up_",x)) })
names(upsetInput$dn) <- lapply(names(upsetInput$dn), function(x){ return(paste0("down_",x)) })
upsetInput.merge <- c(upsetInput$up, upsetInput$dn)

# for each treatment, do GO per cell type
for (group in c("LPS")) { # search the comparison that includes LPS in upsetInput.merge
  upsetInput.group <- upsetInput.merge[grep(group, names(upsetInput.merge), value = TRUE)]
  
  # Extract the individual gene sets
  down_LPS_4h <- upsetInput.group[["down_LPS_4h/unstim_0h"]]
  down_LPS_18h <- upsetInput.group[["down_LPS_18h/unstim_0h"]]
  up_LPS_4h <- upsetInput.group[["up_LPS_4h/unstim_0h"]]
  up_LPS_18h <- upsetInput.group[["up_LPS_18h/unstim_0h"]]
  
  unique_down_LPS_4h <- setdiff(down_LPS_4h, down_LPS_18h)
  unique_down_LPS_18h <- setdiff(down_LPS_18h, down_LPS_4h)
  unique_up_LPS_4h <- setdiff(up_LPS_4h, up_LPS_18h)
  unique_up_LPS_18h <- setdiff(up_LPS_18h, up_LPS_4h)
  
  # Find overlapping gene sets
  overlap_down_LPS <- intersect(down_LPS_4h, down_LPS_18h)
  overlap_up_LPS <- intersect(up_LPS_4h, up_LPS_18h)
  
  # Add these unique sets back to upsetInput.group
  upsetInput.group[["unique_down_LPS_4h"]] <- unique_down_LPS_4h
  upsetInput.group[["unique_down_LPS_18h"]] <- unique_down_LPS_18h
  upsetInput.group[["unique_up_LPS_4h"]] <- unique_up_LPS_4h
  upsetInput.group[["unique_up_LPS_18h"]] <- unique_up_LPS_18h
  upsetInput.group[["overlap_down_LPS"]] <- overlap_down_LPS
  upsetInput.group[["overlap_up_LPS"]] <- overlap_up_LPS
  
  
  upsetInput.group.cnts <- data.frame(
    treatment = c("unique_up_LPS_4h", 
                  "unique_up_LPS_18h", 
                  "unique_down_LPS_4h", 
                  "unique_down_LPS_18h",
                  "overlap_down_LPS",
                  "overlap_up_LPS"),
    de_gene_cnts = c(length(unique_up_LPS_4h), 
                     length(unique_up_LPS_18h), 
                     length(unique_down_LPS_4h), 
                     length(unique_down_LPS_18h),
                     length(overlap_down_LPS),
                     length(overlap_up_LPS)))
  
  wb <- createWorkbook()
  addWorksheet(wb, "de_gene_cnts")
  writeData(wb, sheet = "de_gene_cnts", upsetInput.group.cnts)
  
  addWorksheet(wb, "unique_up_LPS_4h")
  writeData(wb, sheet = "unique_up_LPS_4h", data.frame(unique_up_LPS_4h))
  
  addWorksheet(wb, "overlap_up_LPS")
  writeData(wb, sheet = "overlap_up_LPS", data.frame(overlap_up_LPS))
  
  addWorksheet(wb, "unique_up_LPS_18h")
  writeData(wb, sheet = "unique_up_LPS_18h", data.frame(unique_up_LPS_18h))
  
  addWorksheet(wb, "unique_down_LPS_4h")
  writeData(wb, sheet = "unique_down_LPS_4h", data.frame(unique_down_LPS_4h))
  
  addWorksheet(wb, "unique_down_LPS_18h")
  writeData(wb, sheet = "unique_down_LPS_18h", data.frame(unique_down_LPS_18h))
  
  addWorksheet(wb, "overlap_down_LPS")
  writeData(wb, sheet = "overlap_down_LPS", data.frame(overlap_down_LPS))
  
  saveWorkbook(wb, file = paste0(deResDir, excel.dir, 'GOEnrichment/GO_', group, "_", ct, '_de_genes.xlsx'), overwrite = TRUE)
  
  # write_xlsx(upsetInput.group.cnts, path = paste0(deResDir,excel.dir,'GOEnrichment/', group, "_", ct,'_de_genes_non_overlapping_cnts.xlsx'), col_names = T)
  
  # gost.res <- gost(upsetInput.group[grep("unique", names(upsetInput.group), value = TRUE)], organism = organism, correction_method = "fdr")
  
  gene_sets <- c(
    upsetInput.group[grep("unique", names(upsetInput.group), value = TRUE)],
    upsetInput.group[grep("overlap", names(upsetInput.group), value = TRUE)]
  )
  
  gost.res <- gost(gene_sets, organism = organism, correction_method = "fdr", evcodes=FALSE)
  
  gost.res$result$query <- factor(gost.res$result$query, 
                                  levels = c("overlap_up_LPS", "unique_up_LPS_4h", "unique_up_LPS_18h", 
                                             "overlap_down_LPS", "unique_down_LPS_4h", "unique_down_LPS_18h"))
  
  
  mygostplot <- gostplot(gost.res, interactive = F, capped = F)
  ggsave(paste0(deResDir,excel.dir,"GOEnrichment/", group, "_", ct, ".pdf"), #LPS_Bcells_.png
         mygostplot, height = 12, width = 8, units = "in")
  write_xlsx(as.data.frame(gost.res$result), path = paste0(deResDir,excel.dir,'GOEnrichment/', group, "_", ct,'_GO_table.xlsx'), col_names = T)

}


# # skip e.g. LPS_up_4h vs LPS_down_4h across 4 dbs on the same figure----
# 
# group <- "LPS"
# 
# # go %>%
# #   filter(source %in% c('GO:MF', 'GO:CC', 'GO:BP', 'KEGG', 'TF')) %>%
# #   filter(grepl(tp, query)) %>%
# #    filter(source %in% "KEGG") %>%
# #    print(n=40)
# 
# 
# tp <- "4h"
# go_top10_GO_KEGG_TF <- go_up_down(group = group, 
#            ct = ct,
#            tp = tp)$go_top10_GO_KEGG_TF
# ylim <- ceiling(-log10(min(go_top10_GO_KEGG_TF$p_value)))
# 
# 
# go_top10_GO_KEGG_TF %>%
#   filter(query %in% paste0('unique_up_', group, "_", tp)) %>%
#   group_by(term_name) %>%
#   filter(p_value == min(p_value)) %>%
#   ungroup()  # Ungroup after filtering
# 
# go_top10_GO_KEGG_TF %>%
#   filter(query %in% paste0('unique_up_', group, "_", tp)) %>%
#   group_by(term_name) %>%
#   filter(n() > 1) %>%
#   ungroup()
# 
# 
# # lollipop chart
# p1 <- go_top10_GO_KEGG_TF %>%
#   filter(query %in% paste0('unique_up_',group, "_", tp)) %>%
#   group_by(term_name) %>%
#   filter(p_value == min(p_value)) %>%
#   ungroup() %>% # when there are more than one record with the same term name, keep the lowest p value 
#   ggplot(aes(x=factor(term_name, levels = term_name[order(-p_value)]), y=-log10(p_value)))+
#   geom_segment(aes(xend = factor(term_name, levels = term_name[order(-p_value)]), 
#                    yend = 0), 
#                color = "#e63946", linewidth = 0.5) +  # Add the line segment
#   geom_point(color = "#e63946", size = 2) +  # Add the point at the end of the line
#   theme_bw() +
#   # theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
#   # scale_y_reverse(limits = c(40, 0),breaks = seq(40, 0, by = -10)) +
#   scale_y_reverse(limits = c(ylim, 0)) +
#   # scale_y_continuous(limits = c(0, 150),breaks = seq(0, 150, by = 50)) +
#   labs(x = NULL) +
#   facet_grid(source~., scales="free", space="free") +
#   # theme(strip.placement="outside") +
#   coord_flip() +
#   theme(strip.text.y = element_blank() ,
#         strip.background = element_blank()) 
# # + scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) # Place y-axis tick labels on the rCD3ht
# 
# 
# p2 <- go_up_down(group = group, 
#                  ct = ct,
#                  tp = tp)$p2
# pdf(paste0(deResDir,excel.dir, 'GOEnrichment/', group, "_", tp, "_", ct,"_GO_KEGG_TF.pdf"),
#     width = 12,
#     height = 8)
# grid.arrange(p1, p2, ncol=2, widths=c(0.89,1))
# dev.off()
# 
# 
# 
# 
# tp <- "18h"
# p1 <- go_up_down(group = group, 
#                  ct = ct,
#                  tp = tp)$p1
# p2 <- go_up_down(group = group, 
#                  ct = ct,
#                  tp = tp)$p2
# pdf(paste0(deResDir,excel.dir, 'GOEnrichment/', group, "_", tp, "_", ct,"_GO_KEGG_TF.pdf"),
#     width = 12,
#     height = 8)
# grid.arrange(p1, p2, ncol=2, widths=c(0.83,1))
# 
# 
# dev.off()
# 
# 
# 
#e.g. shared and unique up and down (6 conditions) in one db on the same figure ----
group <- "LPS"
# 6 conditions across 4 databases
conditions <- c("overlap_up_LPS", "unique_up_LPS_4h", "unique_up_LPS_18h", 
                "overlap_down_LPS", "unique_down_LPS_4h", "unique_down_LPS_18h")


plots <- go_up_down_db(group = group, ct = ct, conditions)$plots
combined_go_top10_GO_KEGG_TF4 <- go_up_down_db(group = group, ct = ct, conditions)$combined_go_top10_GO_KEGG_TF
combined_go_top4_GO_KEGG_TF4 <- go_up_down_db(group = group, ct = ct, conditions)$combined_go_top4_GO_KEGG_TF


##-- summary plots for the above ####
go_top4_GO_KEGG_TF <- list()
go_top4_GO_KEGG_TF[['B']] <- combined_go_top4_GO_KEGG_TF1
go_top4_GO_KEGG_TF[['CD4 T']] <- combined_go_top4_GO_KEGG_TF2
go_top4_GO_KEGG_TF[['CD8 T']] <- combined_go_top4_GO_KEGG_TF3
go_top4_GO_KEGG_TF[['Mono/Mph']] <- combined_go_top4_GO_KEGG_TF4

combined_go_top4_GO_KEGG_TF <- do.call(rbind, go_top4_GO_KEGG_TF)

conditions <- c("overlap_up", "unique_up_4h", "unique_up_18h",
                "overlap_down", "unique_down_4h", "unique_down_18h")

combined_go_top4_GO_KEGG_TF <- combined_go_top4_GO_KEGG_TF %>%
  rownames_to_column(var = "cell_type")  # Move rownames to a new column named 'new_column'

# Remove anything after the first '.' in the 'new_column'
combined_go_top4_GO_KEGG_TF$cell_type <- gsub("\\..*", "", combined_go_top4_GO_KEGG_TF$cell_type)
combined_go_top4_GO_KEGG_TF$query <- gsub("_CD3", "", combined_go_top4_GO_KEGG_TF$query)
combined_go_top4_GO_KEGG_TF$query <- gsub("_Ig", "", combined_go_top4_GO_KEGG_TF$query)
combined_go_top4_GO_KEGG_TF$query <- gsub("_LPS", "", combined_go_top4_GO_KEGG_TF$query)


combined_go_top4_GO_KEGG_TF <- combined_go_top4_GO_KEGG_TF %>%
  group_by(source, term_id, query) %>%
  mutate(cell_type_count = n_distinct(cell_type)) %>%
  ungroup()


# deal with the match class 1 in TF
combined_go_top4_GO_KEGG_TF <- combined_go_top4_GO_KEGG_TF %>%
  mutate(term_name = ifelse(source == "TF", gsub("; match class: 1", "", term_name), 
                            term_name))

combined_go_top4_GO_KEGG_TF <- combined_go_top4_GO_KEGG_TF %>%
  group_by(cell_type, query, source, term_name) %>%
  filter(!(source == "TF" & p_value != min(p_value))) %>%
  ungroup()


# ylim <- ceiling(-log10(min(combined_go_top4_GO_KEGG_TF$p_value)))

conditions <- c("overlap_up", "unique_up_4h", "unique_up_18h",
                "overlap_down", "unique_down_4h", "unique_down_18h")


color_values <- c(rep("#FE4088", 3), rep("#5F96FA", 3))

names(color_values) <- conditions

# color_labels <- c(paste0("Shared up-regulated"),
#                   paste0("Unique up-regulated 4h"),
#                   paste0("Unique up-regulated 18h"),
#                   paste0("Shared down-regulated"),
#                   paste0("Unique down-regulated 4h"),
#                   paste0("Unique down-regulated 18h"))

color_labels <- c(rep("Up-regulated", 3), rep("Down-regulated", 3))

custom_labels <- setNames(c('Shared', "Unique 4h", "Unique 18h",
                                        'Shared', "Unique 4h", "Unique 18h"), conditions)

combined_plots <- list()
for (db in c("GO:BP", "GO:MF", "TF")){
  # Loop through conditions
  plots <- list()  # To store individual plots
  ylim <- combined_go_top4_GO_KEGG_TF %>%
    group_by(source) %>%  # Group by the 'source' column
    summarize(min_pvalue_log = ceiling(-log10(min(p_value, na.rm = TRUE))))%>%
    filter(source == db) %>%
    pull(min_pvalue_log)
  for (i in seq_along(conditions)) {
    p <- combined_go_top4_GO_KEGG_TF %>%
      # Focus on a specific source and condition
      filter(grepl(db, source)) %>%
      filter(grepl(conditions[i], query)) %>%
      # Sort by cell_type_count, with the most shared term_ids first
      arrange(desc(cell_type_count), .by_group = TRUE) %>%
      # Reorder the term_name based on cell_type_count (for the y-axis order)
      mutate(term_name = factor(term_name, levels = unique(term_name[order(-cell_type_count)]))) %>%
      ggplot(aes(x = term_name, y = -log10(p_value))) +
      geom_segment(aes(xend = term_name, yend = 0, color = query), linewidth = 0.3) +  # Add the line segment
      geom_point(aes(color = query), size = 1.5) +  # Add the point at the end of the line
      scale_color_manual(values = color_values[i],
                         labels = color_labels[i]) +
      theme_bw() +
      coord_flip() +  # Flip coordinates
      scale_y_reverse(limits = c(ylim, 0)) +
      labs(x = NULL, color = "Genes") +
      facet_grid(query ~ cell_type, scales = "free", space = "free",
                 labeller = labeller(query = custom_labels)) +
      theme(
        axis.text.y = element_text(size = 6.5)  # Adjust size as needed
      )

    # Customize the strips and axis labels
    if (i != 1) {
      p <- p + theme(strip.text.x = element_blank(), # Remove strip for non-first plots
                     strip.background = element_blank(),
                     panel.grid.major = element_line(linewidth = 0.2, color = "lightgrey"),  # Major grid lines
                     panel.grid.minor = element_line(linewidth = 0.1, color = "lightgrey"))
    }
    if (i != length(conditions)) {
      p <- p + theme(axis.text.x = element_blank(), # Remove x-axis labels for non-last plots
                     axis.ticks.x = element_blank(),
                     axis.title.x = element_blank(),
                     strip.background = element_blank(),
                     panel.grid.major = element_line(linewidth = 0.2, color = "lightgrey"),  # Major grid lines
                     panel.grid.minor = element_line(linewidth = 0.1, color = "lightgrey"))
    }

    # Add each plot to the list
    plots[[i]] <- p
  }

  # Combine all plots using patchwork
  combined_plot <- wrap_plots(plots, ncol = 1)
  combined_plots[[db]] <- combined_plot
}


combined_plots[["GO:BP"]]
combined_plots[["GO:MF"]]
combined_plots[["TF"]]



# 9.45 by 26
# 9 by 27

pdf(paste0(deResDir, excel.dir, 
           'GOEnrichment/clsuter_B_CD4T_CD8T_MonoMph_TF.pdf'),
    width = 9, height = 10)
print(combined_plots[["TF"]])
dev.off()

pdf(paste0(deResDir, excel.dir, 
           'GOEnrichment/clsuter_B_CD4T_CD8T_MonoMph_GO_MF.pdf'),
    width = 8, height = 8)
print(combined_plots[["GO:MF"]])
dev.off()

pdf(paste0(deResDir, excel.dir, 
           'GOEnrichment/clsuter_B_CD4T_CD8T_MonoMph_GO_BP.pdf'),
    width = 8, height = 8)
print(combined_plots[["GO:BP"]])
dev.off()



wrap_plots(combined_plots, ncol = 1)


#---
#### some quick summary stats for database and GO terms ####

file_list <- list.files(path = paste0(deResDir,excel.dir, 'GOEnrichment'), pattern = "GO_table.xlsx$")

file_list <- file_list[c(1,2,3,4)]

results <- data.frame(source = character(),
                      term_id_count = character(),
                      condition = character(),
                       stringsAsFactors = FALSE)


for (file in file_list) {
  go_data <- read_excel(paste0(deResDir,excel.dir, 'GOEnrichment/',file))

  # count the number of terms per database per celltype
  result <- go_data %>%
    filter(significant == TRUE) %>%
    # select(term_id, source, query) %>%
    group_by(source) %>%
    summarize(term_id_count = n_distinct(term_id)) %>%
    ungroup()
  
  condition <- sub("_GO_table.xlsx$", "", file)
  result$condition <- condition
  
  results <- rbind(results, result)
  
  output_file <- paste0(sub("GO_table.xlsx$", "GO_summary_result.xlsx", file))
  
  # write_xlsx(result, paste0(deResDir,excel.dir, 'GOEnrichment/', output_file))

  # Print a message to indicate that the file has been processed
  cat("Processed:", file, "-> Saved:", output_file, "\n")
}


write_xlsx(results, paste0(deResDir,excel.dir, 'GOEnrichment/GO_cell_type_summary_results.xlsx'))



combined_results <- data.frame()

# Loop through the list of files
for (file in file_list) {
  
  # Read the current Excel file and select the necessary columns
  go_data <- read_excel(paste0(deResDir, excel.dir, 'GOEnrichment/', file)) %>%
    filter(significant == TRUE) %>%
    dplyr::select(term_id, source)
  
  # Append the current data frame to the combined data frame
  combined_results <- rbind(combined_results, go_data)
}

combined_results_sum <- combined_results %>%
  group_by(source) %>%
  summarize(term_id_count = n_distinct(term_id)) %>%
  ungroup()

# # A tibble: 11 × 2
# source term_id_count
# <chr>          <int>
#   1 CORUM            114
# 2 GO:BP           4209
# 3 GO:CC            825
# 4 GO:MF            701
# 5 HP              1258
# 6 HPA              543
# 7 KEGG             218
# 8 MIRNA            704
# 9 REAC             792
# 10 TF              3019
# 11 WP               348

write_xlsx(combined_results_sum, paste0(deResDir,excel.dir, 'GOEnrichment/GO_summary_results.xlsx'))



# the database

# Get available databases and their details\=
# test should be chnaged to a gene name 
databases <- gost(query = "test", organism = "hsapiens")$result$source %>% unique()

# Print the databases
print(databases)


#### look into IRF, E2F gene expresssion ####
# library(pheatmap)
library(ComplexHeatmap)
rds <- readRDS('integration_2/integration_2_leiden/RDS_Dir/integration_2_leiden_annot.rds')
DefaultAssay(rds) <- "RNA"

rds@meta.data <- rds@meta.data %>%
  # mutate(cluster_annotation_plot = str_replace_all(cluster_annotation, 
  #                                                  c("Alveolar Mph MT-positive" = "Alv Mph",
  #                                                    "EC general capillary" = "EC",
  #                                                    "Monocyte-derived Mph" = "Mono Mph",
  #                                                    "B cells" = "B",
  #                                                    "CD4 T cells" = "CD4 T",
  #                                                    "CD8 T cells" = "CD8 T",
  #                                                    "NK cells" = "NK",
  #                                                    "Migratory DC" = "DC")))


rds@meta.data$clt_ann_con <- paste(rds@meta.data$cluster_annotation_plot, rds@meta.data$expCond.stimuli.time, sep = " ")



# irf_genes <- Filter(function(x) grepl("^IRF", x), rownames(GetAssayData(object = rds, layer = "data")))
# [1] "IRF6"    "IRF2BP2" "IRF2"    "IRF1"    "IRF4"    "IRF5"    "IRF7"    "IRF9"    "IRF2BPL" "IRF8"    "IRF2BP1" "IRF3"   
# e2f_genes <- Filter(function(x) grepl("^E2F", x), rownames(GetAssayData(object = rds, layer = "data")))
# [1] "E2F2"     "E2F6"     "E2F3"     "E2F3-IT1" "E2F5"     "E2F7"     "E2F4"     "E2F1"     "E2F8"    

irf_genes <- c("IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF6", "IRF7", "IRF8", "IRF9")
e2f_genes <- c("E2F1", "E2F2", "E2F3", "E2F4", "E2F5", "E2F6", "E2F7", "E2F8")

to_keep <- c("B unstim_0h", "B Ig_4h", "B Ig_18h",
             "CD4 T unstim_0h", "CD4 T CD3_4h", "CD4 T CD3_18h",
             "CD8 T unstim_0h", "CD8 T CD3_4h", "CD8 T CD3_18h",
             "Mono Mph unstim_0h", "Mono Mph LPS_4h", "Mono Mph LPS_18h")
rds_sub <- subset(rds, cells = rownames(rds@meta.data[rds@meta.data$clt_ann_con %in% to_keep, ]))


rds_sub_scaled <- Seurat::ScaleData(object = rds_sub, do.scale = T, do.center = T, features = c(irf_genes, e2f_genes)) 


irf_e2f_exp <- as.data.frame(t(GetAssayData(object = rds_sub_scaled, layer = "scale.data")[c(irf_genes, e2f_genes), ]))


metadata <- rds_sub_scaled@meta.data[rownames(irf_e2f_exp), ]
irf_e2f_exp <- t(irf_e2f_exp)

irf_e2f_exp[irf_e2f_exp < -2.5] <- -2.5
irf_e2f_exp[irf_e2f_exp > 2.5] <- 2.5

# 
# 
# ordered_indices <- order(metadata$clt_ann_con)
# irf_e2f_exp_ordered <- irf_e2f_exp[ordered_indices, ]
# metadata_ordered <- metadata[ordered_indices, , drop = FALSE]
# 

# # Plot heatmap with row split by 'clt_ann' and without row names
# pheatmap(irf_e2f_exp_ordered, 
#          cluster_rows = FALSE, 
#          cluster_cols = FALSE, 
#          color = colorRampPalette(c("white", "blue"))(100),
#          main = "IRF and E2F Expression Heatmap",
#          annotation_row =  metadata_ordered["clt_ann_con"],  # Add row annotation
#          show_rownames = FALSE)  # Hide row names
# 
# pheatmap(irf_e2f_exp_ordered, 
#          cluster_rows = FALSE, 
#          cluster_cols = FALSE, 
#          color = colorRampPalette(c("white", "blue"))(100),
#          main = paste("Heatmap for", factor),
#          annotation_row = row_annotation_subset,  # Add row annotation
#          show_rownames = FALSE,
#          silent = TRUE) 
# 
# 

##-- complexheatmap ####
myCol <- colorRampPalette(c("#3B4CC0", "#FFFFFF", "#B40426"))(100)
# myCol <- colorRampPalette(c("#3B4CC0", "#89A0D0", "#FFFFFF", "#D08080", "#B40426"))(100)

# myCol <- colorRampPalette(brewer.pal(5, "RdYlBu"))(100)

myBreaks <- seq(-2.5, 2.5, length.out = 100)


col_ann <- data.frame(
  Celltypes = factor(metadata$cluster_annotation_plot,levels=c("B", "CD4 T", "CD8 T", "Mono Mph")),
  Conditions = factor(metadata$expCond.stimuli.time, 
                      levels=c('unstim_0h', 'Ig_4h', 'Ig_18h',
                               "CD3_4h", "CD3_18h", "LPS_4h", "LPS_18h")),
  stringsAsFactors = FALSE) #do not have to set factor 

col_colours <- list(
  Celltypes = c("B" = "#d0c1df",
            "CD4 T" = "#f2748f", 
            "CD8 T" = "#99c9a4", 
            "Mono Mph" = "#6d98ba"),
  Conditions = c('unstim_0h'= "#eff9f0", 
                 'Ig_4h'= "#ddc8c4", 
                 'Ig_18h'= "#6b4d57",
                 "CD3_4h"= "#ddc8c4", 
                 "CD3_18h"= "#6b4d57", 
                 "LPS_4h"= "#ddc8c4", 
                 "LPS_18h"= "#6b4d57"))

colAnn <- HeatmapAnnotation(
  df = col_ann,
  which = 'col', # 'col' (samples) or 'row' (gene) annotation?
  na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
  col = col_colours,
  annotation_height = 0.2,
  annotation_width = unit(1, 'cm'),
  simple_anno_size = unit(3, "mm"), # the height of the bar 
  gap = unit(1, 'mm'),
  annotation_legend_param = list(
    Celltypes = list(
      nrow = 4, # legend in 4 row2
      title = 'Celltypes',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 9),  # legend on the side
      labels_gp = gpar(fontsize = 9)),
    Conditions = list(
      nrow = 7, # legend in 4 row2
      title = 'Conditions',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 9),  # legend on the side
      labels_gp = gpar(fontsize = 9)))) # legend on the side


hmap <- Heatmap(irf_e2f_exp,
        column_split = factor(metadata$clt_ann,
                              levels=c("B unstim_0h", "B Ig_4h", "B Ig_18h",
                                     "CD4 T unstim_0h", "CD4 T CD3_4h", "CD4 T CD3_18h",
                                     "CD8 T unstim_0h", "CD8 T CD3_4h", "CD8 T CD3_18h",
                                     "Mono Mph unstim_0h", "Mono Mph LPS_4h", "Mono Mph LPS_18h")),
        # cluster_row_slices = FALSE,
        name = 'Scaled\nexpression',
        use_raster = TRUE,
        col = colorRamp2(myBreaks, myCol),
        # parameters for the colour-bar that represents gradient of expression
        heatmap_legend_param = list(
          at = c(-2.5, -2, -1, 0, 1, 2, 2.5),
          color_bar = 'continuous',
          legend_direction = 'vertical',
          legend_width = unit(6, 'cm'),
          legend_height = unit(4, 'cm'),
          title_position = 'topcenter',
          title_gp=gpar(fontsize = 9), # heatmap y axis title 
          labels_gp=gpar(fontsize = 9)), # heatmap gradient legend
        
        # row (gene) parameters
        cluster_rows = FALSE,
        show_row_dend = FALSE,
        row_title = 'Genes',
        row_title_side = 'left',
        row_title_gp = gpar(fontsize = 9),
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
        # left_annotation = rowAnn

        )

pdf(paste0(deResDir, excel.dir, "GOEnrichment/irf_e2f_complexheatmap.pdf"), height=2.85, width=7.49)
draw(hmap, 
     # + genelabels,
     heatmap_legend_side = 'right',
     annotation_legend_side = 'right')
dev.off()

##-- dotplot ####
rds@meta.data$cluster_annotation_expCond_stimuli_time <- paste0(
  rds@meta.data$cluster_annotation, 
  "|", 
  rds@meta.data$expCond.stimuli.time
)

to_keep <- c("B cells|unstim_0h", "B cells|Ig_4h", "B cells|Ig_18h",
             "CD4 T cells|unstim_0h", "CD4 T cells|CD3_4h", "CD4 T cells|CD3_18h",
             "CD8 T cells|unstim_0h", "CD8 T cells|CD3_4h", "CD8 T cells|CD3_18h",
             "Monocyte-derived Mph|unstim_0h", "Monocyte-derived Mph|LPS_4h", "Monocyte-derived Mph|LPS_18h")

rds_sub <- subset(rds, cells = rownames(rds@meta.data[rds@meta.data$cluster_annotation_expCond_stimuli_time %in% to_keep, ]))


dotplot <- DotPlot(rds_sub,
                   features = c(irf_genes, e2f_genes),
                   cols = c('#D3D3D3', '#CC0000'),
                   scale = T, scale.by = 'size',
                   dot.min = 0) + RotatedAxis()


irf_e2f_exp <- as.data.frame(t(GetAssayData(object = rds_sub, layer = "data")[c(irf_genes, e2f_genes), ]))

# Calculate the average expression per cluster
irf_e2f_avg_exp <- irf_e2f_exp %>%
  rownames_to_column(var = "cell") %>%
  inner_join(rds@meta.data %>% rownames_to_column(var = "cell"), by = "cell") %>%
  group_by(cluster_annotation_expCond_stimuli_time) %>%
  summarise(across(2:(ncol(irf_e2f_exp)+1), mean, na.rm = TRUE))

# Calculate the percentage of cells expressing each gene per cluster
irf_e2f_pct_exp <- irf_e2f_exp %>%
  rownames_to_column(var = "cell") %>%
  inner_join(rds@meta.data %>% rownames_to_column(var = "cell"), by = "cell") %>%
  group_by(cluster_annotation_expCond_stimuli_time) %>%
  summarise(across(2:(ncol(irf_e2f_exp)+1), ~mean(. > 0) * 100))

# Combine average expression and percentage expression data
irf_e2f_avg_pct_exp <- irf_e2f_avg_exp %>%
  pivot_longer(-cluster_annotation_expCond_stimuli_time, names_to = "gene", values_to = "expression") %>%
  inner_join(irf_e2f_pct_exp %>% pivot_longer(-cluster_annotation_expCond_stimuli_time, names_to = "gene", values_to = "pct.exp"), 
             by = c("cluster_annotation_expCond_stimuli_time", "gene"))

# Scale the expression data if needed
irf_e2f_avg_pct_exp <- irf_e2f_avg_pct_exp %>%
  group_by(gene) %>%
  mutate(expression = scale(expression))

irf_e2f_avg_pct_exp <- irf_e2f_avg_pct_exp %>%
  separate(
    cluster_annotation_expCond_stimuli_time, 
    into = c("cluster_annotation", "expCond_stimuli_time"), 
    sep = "\\|"
  )

irf_e2f_avg_pct_exp <- irf_e2f_avg_pct_exp %>%
  mutate(cluster_annotation_plot = case_when(
    cluster_annotation == "Alveolar Mph MT-positive" ~ "Alveolar Mph\nMT-positive",
    cluster_annotation == "EC general capillary" ~ "EC general\ncapillary",
    cluster_annotation == "Monocyte-derived Mph" ~ "Monocyte-\nderived Mph",
    TRUE ~ cluster_annotation  # Default case for other annotations
  ))

irf_e2f_avg_pct_exp <- irf_e2f_avg_pct_exp %>%
  mutate(pct.exp.shape = case_when(
    pct.exp < 25 ~ "<25%",
    pct.exp >= 25 & pct.exp < 50 ~ "25-50%",
    pct.exp >= 50 & pct.exp < 75 ~ "50-75%",
    pct.exp >= 75 ~ "75-100%"
  ))

irf_e2f_avg_pct_exp <- irf_e2f_avg_pct_exp %>%
  mutate(gene_cat = ifelse(grepl("^IRF", gene), "IRF", "E2F"))


irf_e2f_avg_pct_exp <- irf_e2f_avg_pct_exp %>%
  mutate(cluster_annotation_plot = str_replace_all(cluster_annotation_plot, 
                                                   c("Alveolar Mph\nMT-positive" = "Alv Mph",
                                                     "EC\ngeneral capillary" = "EC",
                                                     "Monocyte-\nderived Mph" = "Mono Mph",
                                                     "B cells" = "B",
                                                     "CD4 T cells" = "CD4 T",
                                                     "CD8 T cells" = "CD8 T",
                                                     "NK cells" = "NK",
                                                     "Migratory DC" = "DC")))

irf_e2f_avg_pct_exp <- irf_e2f_avg_pct_exp %>%
  mutate(expCond_stimuli_time = str_replace_all(expCond_stimuli_time, 
                                                   c("unstim_0h" = "Unstim_0h")))

pdf(paste0(deResDir, excel.dir, "GOEnrichment/irf_e2f_dotplot.pdf"), height=5.55, width=5.41)

ggplot(irf_e2f_avg_pct_exp, aes(x = factor(expCond_stimuli_time, 
                                           levels=c("Unstim_0h", "LPS_4h", "LPS_18h", "CD3_4h", "CD3_18h","Ig_4h", "Ig_18h")), 
                                   y = factor(gene, levels=unique(irf_e2f_avg_pct_exp$gene)),
                                   size = pct.exp, 
                                   fill = expression)) +
  # scale_shape_manual(values =  c("<25%" = 21, "25-50%" = 22, "50-75%" = 23, "75-100%" = 24)) +
  geom_point(shape=21) +
  # scale_color_manual(values =  c("<0.01" = "#0000ff", "Not significant" = "white", "NA" = "white")) +
  # scale_discrete_manual(aesthetics = "stroke", 
  #                       values =  c("<0.01" = 1, "Not significant" = 0, "NA" = 0), 
  #                       guide = "none") +
  scale_fill_gradientn(colors = colorRampPalette(c('#f9dbbd', '#ffa5ab', "#da627d", "#a53860", "#450920"))(100)) +
  facet_grid(gene_cat ~ cluster_annotation_plot, scales="free", space="free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "left",
        strip.background = element_rect(fill = "white", color = "black", size = 0.5),
        strip.text.y = element_blank(),
        strip.background.y = element_blank()) +
  labs(x = "Conditions",
       y = "Genes",
       size = "Percent expressed",
       fill = "Average expression")

dev.off()




################# GO for DEGs (MAST) | B unstim and ig | HLA-DQ2 vs None ################# 
# load functions getUpSetInput & go_up_down ----
#### B unstim & ig | HLA-DQ2 vs None ####
# Load DEGs for each contrast ----
excel.dir <- "results_wOrgClusterAnnotation_DEGs/HLA-DQA2|DQB2_ig_unstim_MAST/"
file.all <- list.files(paste0(deResDir,excel.dir))
files <- file.all[grepl("full.*SelClusters\\.xlsx$", file.all, ignore.case = TRUE)]
contrast <- gsub("expCondCompDeMarkers_|_full.*SelClusters\\.xlsx$", "", files)

# save each cell type to a tibble for each comparison 
for (i in 1:length(files)) {
  excel_file <- paste0(deResDir, excel.dir, files[i])
  sheet_names <- excel_sheets(excel_file)
  for (sheet_name in sheet_names) {
    assign(paste0(contrast[i], "_", sheet_name), read_excel(excel_file, sheet = sheet_name))
  }
}


go_dir <- "GOEnrichment"

# Check if the directory exists
if (!dir.exists(paste0(deResDir, excel.dir, go_dir))) {
  # If it doesn't exist, create the directory
  dir.create(paste0(deResDir,excel.dir, go_dir))
  cat("Directory 'GOEnrichment' created.\n")
} else {
  cat("Directory 'GOEnrichment' already exists.\n")
}

celltype <- "B"
ct <- gsub("/| |-","",celltype) # replace any / space -
ct <- paste0("cluster_", ct)
upsetInput <- getUpSetInput(ct)
names(upsetInput$up) <- lapply(names(upsetInput$up), function(x){ return(paste0("up_",x)) })
names(upsetInput$dn) <- lapply(names(upsetInput$dn), function(x){ return(paste0("down_",x)) })
upsetInput.merge <- c(upsetInput$up, upsetInput$dn)

gene_sets <- c(
  upsetInput.merge[grep("up", names(upsetInput.merge), value = TRUE)],
  upsetInput.merge[grep("down", names(upsetInput.merge), value = TRUE)]
)

gost.res <- gost(gene_sets, organism = organism, correction_method = "fdr", evcodes=FALSE)

gost.res$result$query <- factor(gost.res$result$query, 
                                levels = c("up_HLA/DQA2|DQB2/None", 
                                           "down_HLA/DQA2|DQB2/None"))


mygostplot <- gostplot(gost.res, interactive = F, capped = F)
ggsave(paste0(deResDir, excel.dir,"GOEnrichment/B_unstim_ig_HLA.pdf"),
       mygostplot, height = 4, width = 8, units = "in")
write_xlsx(as.data.frame(gost.res$result), 
           path = paste0(deResDir, excel.dir,
                         'GOEnrichment/B_unstim_ig_HLA_GO_table.xlsx'), 
           col_names = T)

#e.g. shared and unique up and down (2 conditions) in one db on the same figure ----
# 2 conditions across 4 databases
conditions <- c("up_HLA/DQA2|DQB2/None", 
                "down_HLA/DQA2|DQB2/None")

go <- read_excel(paste0(deResDir, excel.dir, 
                        'GOEnrichment/B_unstim_ig_HLA_GO_table.xlsx'), 
                 sheet = 'Sheet1')

if (any(go$significant == FALSE)) {
  print("There is at least one record where 'significant.' is FALSE.")
} else {
  print("All records are significant.")
}

# Loop over each condition and get top 10 terms per database (GO:MF, GO:CC, GO:BP, KEGG, TF)
go_top10_GO_KEGG_TF <- list()

go_top4_GO_KEGG_TF <- list()

for (condition in conditions){
  # Filter and select top 10 terms for each database
  top_terms_10 <- go %>%
    filter(source %in% c('GO:MF', 'GO:CC', 'GO:BP', 'KEGG', 'TF')) %>%
    filter(grepl(condition, query)) %>%
    group_by(query, source) %>%
    arrange(p_value) %>%
    slice_head(n = 10) %>%
    ungroup()
  
  # Filter and select top 4 terms for each database
  top_terms_4 <- go %>%
    filter(source %in% c('GO:MF', 'GO:CC', 'GO:BP', 'KEGG', 'TF')) %>%
    filter(grepl(condition, query)) %>%
    group_by(query, source) %>%
    arrange(p_value) %>%
    slice_head(n = 4) %>%
    ungroup()
  
  # top_terms$source <- factor(top_terms$source,
  #                            levels=c('GO:MF', 'GO:CC', 'GO:BP', 'KEGG', 'TF'))
  
  # Store results in the list, each entry named by the condition
  go_top10_GO_KEGG_TF[[condition]] <- top_terms_10
  go_top4_GO_KEGG_TF[[condition]] <- top_terms_4
  
}

combined_go_top10_GO_KEGG_TF <- do.call(rbind, go_top10_GO_KEGG_TF)
combined_go_top4_GO_KEGG_TF <- do.call(rbind, go_top4_GO_KEGG_TF)

ylim <- ceiling(-log10(min(combined_go_top10_GO_KEGG_TF$p_value)))

plots <-  list()

color_values <- c(rep("#FE4088", 1), rep("#5F96FA", 1))

names(color_values) <- conditions

color_labels <- c("Up-regulated", 
                  "Down-regulated")

for (db in unique(combined_go_top10_GO_KEGG_TF$source)){
  db_name<- str_replace_all(db, ":", "_")
  p <- combined_go_top10_GO_KEGG_TF %>%
    filter(grepl(db, source)) %>%
    group_by(query) %>%
    arrange(query, desc(p_value), .by_group = TRUE) %>%
    # Create a unique ordering for term_name within each query
    mutate(term_name_order = factor(term_name, levels = unique(term_name[order(-p_value)]))) %>%
    mutate(query = factor(query, levels = conditions)) %>%
    ggplot(aes(x = term_name_order, y = -log10(p_value))) +
    geom_segment(aes(xend = term_name_order, yend = 0, color = query), linewidth = 0.5) +  # Add the line segment
    geom_point(aes(color = query), size = 2) +  # Add the point at the end of the line
    scale_color_manual(values = color_values,
                       labels = color_labels) +
    theme_bw() +
    coord_flip() +  # Flip coordinates
    scale_y_reverse(limits = c(ylim, 0)) +
    labs(x = NULL, color = "Genes") +
    facet_grid(query ~ source, scales = "free", space = "free") +
    theme(strip.text.y = element_blank() ,
          strip.background = element_blank()) 
  # +
  # scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 40))
  
  # Save the plot to a PDF file
  ggsave(filename = paste0(deResDir, excel.dir, "GOEnrichment/B_unstim_ig_HLA_", db_name, ".pdf"), plot = p, width = 8, height = 7)
  
  # Add the plot to the list
  plots[[db]] <- p
}

plots 
combined_go_top10_GO_KEGG_TF
combined_go_top4_GO_KEGG_TF



##-- summary plots for the above ####
# go_top4_GO_KEGG_TF <- list()
# go_top4_GO_KEGG_TF[['HLA-DQA2/DQB2']] <- combined_go_top4_GO_KEGG_TF


conditions <- c("up_HLA/DQA2|DQB2/None", 
                "down_HLA/DQA2|DQB2/None")

combined_go_top4_GO_KEGG_TF <- combined_go_top4_GO_KEGG_TF %>%
  rownames_to_column(var = "cell_type")  # Move rownames to a new column named 'new_column'

combined_go_top4_GO_KEGG_TF$cell_type <- "B"


# deal with the match class 1 in TF
combined_go_top4_GO_KEGG_TF <- combined_go_top4_GO_KEGG_TF %>%
  mutate(term_name = ifelse(source == "TF", gsub("; match class: 1", "", term_name), 
                            term_name))

combined_go_top4_GO_KEGG_TF <- combined_go_top4_GO_KEGG_TF %>%
  group_by(cell_type, query, source, term_name) %>%
  filter(!(source == "TF" & p_value != min(p_value))) %>%
  ungroup()


color_values <- c(rep("#FE4088", 1), rep("#5F96FA", 1))

names(color_values) <- conditions


color_labels <- c(rep("Up-regulated", 1), rep("Down-regulated", 1))

custom_labels <- setNames(c("", 
                            ""), conditions)
cell_type_labels <- c("B" = "B Unstim and Ig")

combined_plots <- list()
for (db in c("GO:BP", "GO:MF", "TF")){
  # Loop through conditions
  plots <- list()  # To store individual plots
  ylim <- combined_go_top4_GO_KEGG_TF %>%
    group_by(source) %>%  # Group by the 'source' column
    summarize(min_pvalue_log = ceiling(-log10(min(p_value, na.rm = TRUE))))%>%
    filter(source == db) %>%
    pull(min_pvalue_log)
  for (i in seq_along(conditions)) {
    p <- combined_go_top4_GO_KEGG_TF %>%
      # Focus on a specific source and condition
      filter(grepl(db, source)) %>%
      filter(grepl(conditions[i], query, fixed = TRUE)) %>%
      ggplot(aes(x = term_name, y = -log10(p_value))) +
      geom_segment(aes(xend = term_name, yend = 0, color = query), linewidth = 0.3) +  # Add the line segment
      geom_point(aes(color = query), size = 1.5) +  # Add the point at the end of the line
      scale_color_manual(values = color_values[i],
                         labels = color_labels[i]) +
      theme_bw() +
      coord_flip() +  # Flip coordinates
      scale_y_reverse(limits = c(ylim, 0)) +
      labs(x = NULL, color = "HLA-DQA2/DQB2 expressed") +
      facet_grid(query ~ cell_type, 
                 scales = "free",
                 space = "free",
                 labeller = labeller(query = custom_labels,
                                     cell_type = cell_type_labels)) +
      theme(
        axis.text.x = element_text(size = 6.5)
        )
    
    # Customize the strips and axis labels
    if (i != 1) {
      p <- p + theme(strip.text.x = element_blank(), # Remove strip for non-first plots
                     strip.background = element_blank(),
                     panel.grid.major = element_line(linewidth = 0.2, color = "lightgrey"),  # Major grid lines
                     panel.grid.minor = element_line(linewidth = 0.1, color = "lightgrey"))
    }
    if (i != length(conditions)) {
      p <- p + theme(
                     axis.ticks.x = element_blank(),
                     axis.title.x = element_blank(),
                     strip.background = element_blank(),
                     panel.grid.major = element_line(linewidth = 0.2, color = "lightgrey"),  # Major grid lines
                     panel.grid.minor = element_line(linewidth = 0.1, color = "lightgrey"))
    }
    
    # Add each plot to the list
    plots[[i]] <- p
  }
  
  # Combine all plots using patchwork
  combined_plot <- wrap_plots(plots, ncol = 1)
  combined_plots[[db]] <- combined_plot
}



pdf(paste0(deResDir, excel.dir, 
           'GOEnrichment/B_unstim_ig_HLA_BP_MF_TF.pdf'),
    width = 7.07, height = 6.80)

wrap_plots(
  combined_plots[["GO:BP"]], 
  combined_plots[["GO:MF"]], 
  combined_plots[["TF"]], 
  ncol = 1
) + plot_layout(heights = c(1, 1, 1), 
                guides = "collect") 

dev.off()




#---

#### some quick summary stats for database and GO terms ####

file_list <- list.files(path = paste0(deResDir,excel.dir, 'GOEnrichment'), pattern = "GO_table.xlsx$")

file_list <- file_list[c(1)]

results <- data.frame(source = character(),
                      term_id_count = character(),
                      condition = character(),
                      stringsAsFactors = FALSE)


for (file in file_list) {
  go_data <- read_excel(paste0(deResDir,excel.dir, 'GOEnrichment/',file))
  
  # count the number of terms per database per celltype
  result <- go_data %>%
    filter(significant == TRUE) %>%
    # select(term_id, source, query) %>%
    group_by(source) %>%
    summarize(term_id_count = n_distinct(term_id)) %>%
    ungroup()
  
  condition <- sub("_GO_table.xlsx$", "", file)
  result$condition <- condition
  
  results <- rbind(results, result)
  
  output_file <- paste0(sub("GO_table.xlsx$", "GO_summary_result.xlsx", file))
  
  # write_xlsx(result, paste0(deResDir,excel.dir, 'GOEnrichment/', output_file))
  
  # Print a message to indicate that the file has been processed
  cat("Processed:", file, "-> Saved:", output_file, "\n")
}


write_xlsx(results, paste0(deResDir, excel.dir, 'GOEnrichment/GO_B_unstim_ig_HLA_summary_results.xlsx'))



# the database

# Get available databases and their details\=
# test should be changed to a gene name 
databases <- gost(query = "test", organism = "hsapiens")$result$source %>% unique()

# Print the databases
print(databases)


#### look into gene expresssion of TF ####
##-- dotplot ####
rds <- readRDS('integration_2/integration_2_leiden/RDS_Dir/integration_2_leiden_annot.rds')
DefaultAssay(rds) <- "RNA"

rds@meta.data$cluster_annotation3_expCond_stimuli_time <- paste0(
  rds@meta.data$cluster_annotation3, 
  "|", 
  rds@meta.data$expCond.stimuli.time
)

# irf_genes <- Filter(function(x) grepl("^IRF", x), rownames(GetAssayData(object = rds, layer = "data")))
# [1] "IRF6"    "IRF2BP2" "IRF2"    "IRF1"    "IRF4"    "IRF5"    "IRF7"    "IRF9"    "IRF2BPL" "IRF8"    "IRF2BP1" "IRF3"   
# e2f_genes <- Filter(function(x) grepl("^E2F", x), rownames(GetAssayData(object = rds, layer = "data")))
# [1] "E2F2"     "E2F6"     "E2F3"     "E2F3-IT1" "E2F5"     "E2F7"     "E2F4"     "E2F1"     "E2F8"    

tf_genes <- c("SREBP−1", "MAD5", "HEY2", "HES7", "ZF5", "SP1", "E2F4")
existing_genes <- intersect(unique(tf_genes), rownames(GetAssayData(object = rds, slot = "data")))
# [1] "SP1"  "E2F4" "HES7" "HEY2"

tf_genes <- c("SP1", "E2F4", "HES7", "HEY2")

intersect("SREBP1", rownames(GetAssayData(object = rds, slot = "data")))
Filter(function(x) grepl("^ZF", x), rownames(GetAssayData(object = rds, slot = "data")))

to_keep <- c("B|unstim_0h", "B|Ig_4h", "B|Ig_18h",
             "CD4 T|unstim_0h", "CD4 T|CD3_4h", "CD4 T|CD3_18h",
             "CD8 T|unstim_0h", "CD8 T|CD3_4h", "CD8 T|CD3_18h",
             "Mono/Mph|unstim_0h", "Mono/Mph|LPS_4h", "Mono/Mph|LPS_18h")

rds_sub <- subset(rds, cells = rownames(rds@meta.data[rds@meta.data$cluster_annotation3_expCond_stimuli_time %in% to_keep, ]))


dotplot <- DotPlot(rds_sub,
                   features = c(tf_genes),
                   cols = c('#D3D3D3', '#CC0000'),
                   scale = T, scale.by = 'size',
                   dot.min = 0) + RotatedAxis()


tf_exp <- as.data.frame(t(GetAssayData(object = rds_sub, layer = "data")[c(tf_genes), ]))

# Calculate the average expression per cluster
tf_avg_exp <- tf_exp %>%
  rownames_to_column(var = "cell") %>%
  inner_join(rds@meta.data %>% rownames_to_column(var = "cell"), by = "cell") %>%
  group_by(cluster_annotation3_expCond_stimuli_time) %>%
  summarise(across(2:(ncol(tf_exp)+1), mean, na.rm = TRUE))

# Calculate the percentage of cells expressing each gene per cluster
tf_pct_exp <- tf_exp %>%
  rownames_to_column(var = "cell") %>%
  inner_join(rds@meta.data %>% rownames_to_column(var = "cell"), by = "cell") %>%
  group_by(cluster_annotation3_expCond_stimuli_time) %>%
  summarise(across(2:(ncol(tf_exp)+1), ~mean(. > 0) * 100))

# Combine average expression and percentage expression data
tf_avg_pct_exp <- tf_avg_exp %>%
  pivot_longer(-cluster_annotation3_expCond_stimuli_time, names_to = "gene", 
               values_to = "expression") %>%
  inner_join(tf_pct_exp %>% pivot_longer(-cluster_annotation3_expCond_stimuli_time, 
                                         names_to = "gene", values_to = "pct.exp"), 
             by = c("cluster_annotation3_expCond_stimuli_time", "gene"))

# Scale the expression data if needed
tf_avg_pct_exp <- tf_avg_pct_exp %>%
  group_by(gene) %>%
  mutate(expression = scale(expression))

tf_avg_pct_exp <- tf_avg_pct_exp %>%
  separate(
    cluster_annotation3_expCond_stimuli_time, 
    into = c("cluster_annotation", "expCond_stimuli_time"), 
    sep = "\\|"
  )

tf_avg_pct_exp <- tf_avg_pct_exp %>%
  mutate(cluster_annotation_plot = case_when(
    # cluster_annotation == "Alveolar Mph MT-positive" ~ "Alveolar Mph\nMT-positive",
    # cluster_annotation == "EC general capillary" ~ "EC general\ncapillary",
    # cluster_annotation == "Monocyte-derived Mph" ~ "Monocyte-\nderived Mph",
    TRUE ~ cluster_annotation  # Default case for other annotations
  ))

tf_avg_pct_exp <- tf_avg_pct_exp %>%
  mutate(pct.exp.shape = case_when(
    pct.exp < 25 ~ "<25%",
    pct.exp >= 25 & pct.exp < 50 ~ "25-50%",
    pct.exp >= 50 & pct.exp < 75 ~ "50-75%",
    pct.exp >= 75 ~ "75-100%"
  ))

# tf_avg_pct_exp <-tf_avg_pct_exp %>%
#   mutate(gene_cat = ifelse(grepl("^IRF", gene), "IRF", "E2F"))

tf_avg_pct_exp <- tf_avg_pct_exp %>%
  mutate(gene_cat = case_when(
    gene == "HEY2" ~ "Up-regulated",
    gene == "HES7" ~ "Up-regulated",
    gene == "SP1" ~ "Dowm-regulated",
    gene == "E2F4" ~ "Dowm-regulated",
    TRUE ~ gene  # Default case for other annotations
  ))


# tf_avg_pct_exp <- tf_avg_pct_exp %>%
#   mutate(cluster_annotation_plot = str_replace_all(cluster_annotation_plot, 
#                                                    c("Alveolar Mph\nMT-positive" = "Alv Mph",
#                                                      "EC\ngeneral capillary" = "EC",
#                                                      "Monocyte-\nderived Mph" = "Mono Mph",
#                                                      "B cells" = "B",
#                                                      "CD4 T cells" = "CD4 T",
#                                                      "CD8 T cells" = "CD8 T",
#                                                      "NK cells" = "NK",
#                                                      "Migratory DC" = "DC")))

tf_avg_pct_exp <- tf_avg_pct_exp %>%
  mutate(expCond_stimuli_time = str_replace_all(expCond_stimuli_time, 
                                                c("unstim_0h" = "Unstim_0h")))

pdf(paste0(deResDir, excel.dir, "GOEnrichment/tf_dotplot.pdf"), height=3.03, width=5.10)

ggplot(tf_avg_pct_exp, aes(x = factor(expCond_stimuli_time, 
                                      levels=c("Unstim_0h", "LPS_4h", "LPS_18h", "CD3_4h", "CD3_18h","Ig_4h", "Ig_18h")), 
                                y = factor(gene, levels=unique(tf_avg_pct_exp$gene)),
                                size = pct.exp, 
                                fill = expression)) +
  # scale_shape_manual(values =  c("<25%" = 21, "25-50%" = 22, "50-75%" = 23, "75-100%" = 24)) +
  geom_point(shape=21) +
  # scale_color_manual(values =  c("<0.01" = "#0000ff", "Not significant" = "white", "NA" = "white")) +
  # scale_discrete_manual(aesthetics = "stroke", 
  #                       values =  c("<0.01" = 1, "Not significant" = 0, "NA" = 0), 
  #                       guide = "none") +
  scale_fill_gradientn(colors = colorRampPalette(c('#f9dbbd', '#ffa5ab', "#da627d", "#a53860", "#450920"))(100)) +
  facet_grid(factor(gene_cat,
                    levels=c("Down-regulated", "Up-regulated")) ~ 
               cluster_annotation_plot, scales="free", space="free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "right",
        strip.background = element_rect(fill = "white", color = "black", size = 0.5),
        strip.text.y = element_blank(),
        strip.background.y = element_blank()) +
  labs(x = "Conditions",
       y = "Genes",
       size = "Percent expressed",
       fill = "Scaled expression")

dev.off()



