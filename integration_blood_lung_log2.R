library(Seurat)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(data.table)
library(ggpubr) #cor
library(reshape2)
library(biomaRt)
library(ComplexHeatmap)
library(viridis)
library(tibble)
library(tidyr)
library(readxl)
library(circlize) #colorRamp2

seuratObj_annot <- readRDS('integration_2/integration_2_leiden/RDS_Dir/integration_2_leiden_annot.rds')

deResDir <- "integration_2/integration_2_leiden/"

blResDir <- "blood_lung/"

seuratObj_annot$expCond.celltype.stimuli.time <- paste0(Idents(seuratObj_annot), " ", seuratObj_annot$expCond.stimuli.time)

###################################################################################### 
###################################################################################### 
#### Below is the log 2 scale of gene expression, shared or complete DEGs ####
###################################################################################### 
###################################################################################### 
# 1.0 log2 change of DEGs from the bulk RNA-seq dataset  -----------------------------------
# orginally, we compared the avg exp -> gene (row names) by condition (column nams) matrix
bulk_log_exp <- read.table(paste0(deResDir, blResDir, "41588_2019_505_MOESM7_ESM"), header = T, sep = '\t')

# dim(bulk_log_exp)
# [1] 12011     8

bulk_log_exp[duplicated(bulk_log_exp), ]

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes <- c("ensembl_gene_id", "external_gene_name")  # Attributes to retrieve
gene_names <- getBM(attributes = attributes, filters = "ensembl_gene_id", values = bulk_log_exp$peak_id, mart = mart)
colnames(gene_names) <- c("peak_id", "peak_geneName")
gene_names <- gene_names %>%
  distinct()
bulk_log_exp <- merge(bulk_log_exp, gene_names, by = "peak_id", all.x = TRUE, all.y = FALSE)

sum(is.na(bulk_log_exp$peak_geneName))
dim(bulk_log_exp)

# q-value, was less than 0.01 and the absolute log2FC was greater than 1.
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6858557/
max(bulk_log_exp$adj.P.Val)
min(abs(bulk_log_exp$logFC))


# 2.0 log2 change from scRNA-seq dataset -----------------------------------
# list of protein-coding genes 
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

# get a specific set of deg results (gene name, log 2, p etc) from_getClusterExpCondDe_out_folder
get_degs <- function(excel.dir, compGroups, ident) {
  degs <- list()
  for (cp in compGroups) { #for each condition
    tryCatch({
      file.all <- list.files(excel.dir)
      file <- file.all[grepl(paste0(cp, "_adjSig_7SelClusters"), file.all, ignore.case = TRUE)]
      celltype <- gsub("\\W", "", ident)
      degs[[cp]] <- read_excel(paste0(excel.dir, file), sheet = paste0("cluster_",celltype))
    }, error = function(e) {
      message("Error occurred in compGroup ", cp, ": ", e$message)
    })
    next
  }
  return(degs)
}

y_groups <- list(
  "B" = c("Ig_4h-unstim_0h", "Ig_18h-unstim_0h"),
  "CD4 T" = c("CD3_4h-unstim_0h", "CD3_18h-unstim_0h"),
  "CD8 T" = c("CD3_4h-unstim_0h", "CD3_18h-unstim_0h"),
  "NK" = c("CD3_4h-unstim_0h", "CD3_18h-unstim_0h"),
  "Mono/Mph" = c("LPS_4h-unstim_0h", "LPS_18h-unstim_0h"))

y_stimuli <- c("Ig", "CD3", "CD3", "CD3", "LPS")

# get and merge degs for each y group
# find deg results in each names(y_groups) under the conditions listed in y_groups
y_genes <- lapply(names(y_groups), function(x) get_degs(
  excel.dir = paste0(deResDir, "results_wOrgClusterAnnotation_DEGs/stimuli_time_MAST/"),
  compGroups = y_groups[[x]],
  ident = x))

y_genes_groups <- lapply(seq_along(y_genes), function(i) {
  
  current_element <- y_genes[[i]]
  
  tibble_data_ls <- lapply(names(current_element), function(tibble_name) {
    tibble_data <- current_element[[tibble_name]]
    
    # add condition and cell type
    tibble_data %>%
      mutate(contrast = tibble_name, cell_type = names(y_groups)[i])
  })
  
  # Combine the tibbles in the current element into a single tibble
  bind_rows(tibble_data_ls)
})

y_genes_groups <- bind_rows(y_genes_groups)

colnames(y_genes_groups)[1] <- "geneName"

# protein-coding genes only 
y_genes_groups <- y_genes_groups %>%
  filter(geneName %in% pcgs)

y_genes_groups %>%
  group_by(contrast, cell_type) %>%
  dplyr::summarize(count = n(), .groups = 'drop')


# 3.0 shared DEGs -> calculate correlation -----------------------------------
plotscatter <- function(df2plot) {
  # myplot <- ggscatter(
  #   df2plot, x='Group1', y='Group2',
  #   add = "reg.line",
  #   add.params = list(color = "blue", fill = "lightgray"),
  #   conf.int = TRUE, conf.int.level = 0.9) +
  #   # coord_equal(ratio = 1) +
  #   scale_x_continuous(breaks=seq(0,20,1)) +
  #   scale_y_continuous(breaks=seq(0,20,1)) +
  #   labs(title = paste0("scatter ", x), x = x, y = y) +
  #   stat_regline_equation(label.y = 5.4) +
  #   stat_cor(method = "pearson", label.x.npc = "left", label.y = 5, r.accuracy = 0.01)
  # # ggsave(paste0(resDir, "/scatterPlot_wOrgClusterAnnotation/", n, ".jpg"), plot = myplot, width = 4, height = 4, units = "in")
  # # Extract R and p values
  # r_value <- ggplot2::ggplot_build(myplot)$data[[4]]$r
  # p_value <- ggplot2::ggplot_build(myplot)$data[[4]]$p
  
  lm_model <- lm(Group2 ~ Group1, data = df2plot)
  slope <- coef(lm_model)[2]  # Extract slope
  
  # Compute Pearson correlation
  cor_test <- cor.test(df2plot$Group1, df2plot$Group2, method = "pearson")
  r_value <- cor_test$estimate  # Extract Pearson's r
  p_value <- cor_test$p.value  
  return(list("r_value"=r_value, "p_value"=p_value, "slope"=slope))
}

y_genes_groups$cell_type_contrast <- paste(y_genes_groups$cell_type, y_genes_groups$contrast, sep = " | ")

correlation_values_df <- data.frame()

for (y in unique(y_genes_groups$cell_type_contrast)) {
  for (x in unique(bulk_log_exp$contrast)) {
    tryCatch({
      col1 <- y_genes_groups %>%
        filter(cell_type_contrast %in% y) %>%
        dplyr::select(geneName, avg_log2FC, cell_type_contrast) %>%
        column_to_rownames(var = "geneName")
      col2 <- bulk_log_exp %>%
        filter(contrast %in% x, 
               peak_geneName != "") %>% 
        dplyr::select(peak_geneName, logFC, contrast) %>%
        column_to_rownames(var = "peak_geneName")
      df2plot <- merge(col1, col2, by = "row.names", all = FALSE) %>%
        dplyr::select(Row.names, avg_log2FC, logFC)
      colnames(df2plot) <- c("rowname","Group1","Group2")
      scatter_values <- plotscatter(df2plot)
      new_row <- data.frame(blood = x, 
                            lung = y, 
                            r = scatter_values$r_value, 
                            p = scatter_values$p_value, 
                            slope = scatter_values$slope,
                            shared_genes = dim(df2plot)[1])
      correlation_values_df <- rbind(correlation_values_df, new_row)
      
    }, error = function(e) {
      message("Error occurred in compGroup ", x, " & ", y, ": ", e$message)
    })
    next
  }
}

write.table(correlation_values_df, paste0(deResDir, blResDir, "/bl_shared_log2_correlation_values.txt"), sep = "\t", row.name = FALSE, quote = FALSE)



# 3.0.1 shared DEGs -> pearson's r ->  plot the heatmap ---------------------------------------

# myCol <- colorRampPalette(c("#3B4CC0", "#FFFFFF", "#B40426"))(100)

correlation_values_df <- read.table(paste0(deResDir, blResDir, "/bl_shared_log2_correlation_values.txt"), sep = "\t", header = TRUE)

heatmap_mt <- correlation_values_df[,c(1,2,3)] %>% pivot_wider(names_from = lung, values_from = r) %>% column_to_rownames(var = "blood") %>% as.matrix()
mat_pval <- correlation_values_df[,c(1,2,4)] %>% pivot_wider(names_from = lung, values_from = p) %>% column_to_rownames(var = "blood") %>% as.matrix()
max(mat_pval, na.rm = TRUE)


ct_order <- c("B", "CD4 T","CD8 T","NK","Mono/Mph")
tm_order <- c("Ig","CD3","CD3","CD3","LPS")
timepoint <- c("4h", "18h")

colnames(heatmap_mt)[1] # follow this format ""B cells | Ig_4h-unstim_0h""

group_order <- c()
for (c in seq_along(ct_order)) {
  group_order <- c(group_order, paste0(ct_order[c], " | ", tm_order[c], "_", timepoint[1], "-unstim_0h"))
  group_order <- c(group_order, paste0(ct_order[c], " | ", tm_order[c], "_", timepoint[2], "-unstim_0h"))
}

group_order <- group_order[group_order %in% colnames(heatmap_mt)] 

# > group_order
# [1] "B cells | Ig_4h-unstim_0h"                "B cells | Ig_18h-unstim_0h"               "CD4 T cells | CD3_4h-unstim_0h"          
# [4] "CD4 T cells | CD3_18h-unstim_0h"          "CD8 T cells | CD3_4h-unstim_0h"           "CD8 T cells | CD3_18h-unstim_0h"         
# [7] "NK cells | CD3_4h-unstim_0h"              "NK cells | CD3_18h-unstim_0h"             "Monocyte-derived Mph | LPS_4h-unstim_0h" 
# [10] "Monocyte-derived Mph | LPS_18h-unstim_0h"


x_groups <- list(
  B = c("Bulk_B_S-Bulk_B_U", "Mem_B_S-Mem_B_U", "Naive_B_S-Naive_B_U"),
  
  CD4T = c("Effector_CD4pos_T_S-Effector_CD4pos_T_U", "Follicular_T_Helper_S-Follicular_T_Helper_U", 
           "Memory_Teffs_S-Memory_Teffs_U", "Naive_Teffs_S-Naive_Teffs_U", 
           "Th1_precursors_S-Th1_precursors_U", "Th17_precursors_S-Th17_precursors_U", 
           "Th2_precursors_S-Th2_precursors_U", "Memory_Tregs_S-Memory_Tregs_U", 
           "Naive_Tregs_S-Naive_Tregs_U", "Regulatory_T_S-Regulatory_T_U"), 
  
  CD8T = c("CD8pos_T_S-CD8pos_T_U", "Central_memory_CD8pos_T_S-Central_memory_CD8pos_T_U", 
           "Effector_memory_CD8pos_T_S-Effector_memory_CD8pos_T_U", "Naive_CD8_T_S-Naive_CD8_T_U",
           "Gamma_delta_T_S-Gamma_delta_T_U"),
  
  NK = c("Mature_NK_S-Mature_NK_U"),
  
  Myeloid = c("Monocytes_S-Monocytes_U")
)



us_order <- c(
  # B 
  "Bulk_B_S-Bulk_B_U", "Mem_B_S-Mem_B_U", "Naive_B_S-Naive_B_U",
  
  # CD4 T
  "Effector_CD4pos_T_S-Effector_CD4pos_T_U", "Follicular_T_Helper_S-Follicular_T_Helper_U", 
  "Memory_Teffs_S-Memory_Teffs_U", "Naive_Teffs_S-Naive_Teffs_U", 
  "Th1_precursors_S-Th1_precursors_U", "Th17_precursors_S-Th17_precursors_U", 
  "Th2_precursors_S-Th2_precursors_U", "Memory_Tregs_S-Memory_Tregs_U", 
  "Naive_Tregs_S-Naive_Tregs_U", "Regulatory_T_S-Regulatory_T_U", 
  
  # CD8 T
  "CD8pos_T_S-CD8pos_T_U", "Central_memory_CD8pos_T_S-Central_memory_CD8pos_T_U", 
  "Effector_memory_CD8pos_T_S-Effector_memory_CD8pos_T_U", "Naive_CD8_T_S-Naive_CD8_T_U", 
  "Gamma_delta_T_S-Gamma_delta_T_U",
  
  # NK
  "Mature_NK_S-Mature_NK_U", 
  
  # Mono Mph
  "Monocytes_S-Monocytes_U")

setdiff(rownames(heatmap_mt), us_order)

# match each rowname with a cell type 
# > x_split
# [1] "CD4T"    "CD4T"    "CD4T"    "CD8T"    "CD4T"    "B"       "B"       "B"       "CD8T"    "Myeloid" "CD4T"    "CD4T"    "CD4T"    "CD4T"    "CD8T"    "CD4T"    "NK"      "CD8T"    "CD4T"   
# [20] "NK"   
idx_list <- lapply(x_groups, function(x) which(rownames(heatmap_mt) %in% x))
x_split <- rep(NA, nrow(heatmap_mt))
for (n in names(idx_list)) {
  x_split[idx_list[[n]]] <- n
}

cell_types <- c(
  "Mono/Mph",
  "B", 
  "CD4 T", 
  "CD8 T", 
  "NK",
  "EC", 
  "DC")


selectedCol<-c( 'cornflowerblue', 
                'mediumvioletred',
                'wheat3', 
                'lightsalmon',  
                "#7CE3D8", 
                'yellowgreen',
                'seagreen')

names(selectedCol) <- cell_types
color.order <- c("mediumvioletred", 
                 'wheat3', 
                 'lightsalmon',  
                 "#7CE3D8", 
                 'cornflowerblue')

# significant only
hm <- Heatmap(heatmap_mt, name = "Pearson's r", na_col = "black",
              # col = colorRamp2(c(-0.1, 0, 0.1, 0.2, 0.3, 0.4), c("#60dbe8","white","#ffdd00","#ff9b54", "#EA4749", "#720026")),
              # col = colorRampPalette(c("#172896","#60dbe8","white","#ffdd00", "#EA4749"))(100),
              col = colorRamp2(c(-0.2,  0, 0.2, 0.4, 0.6, 0.8, 1),
                               c("#60dbe8", "white", 
                                 "#ffdd00", "#ff9b54", 
                                 "#EA4749", "#b2001d", 
                                 "#720026")),
              heatmap_legend_param = list(
                title = "Pearson's r",
                at = seq(-0.2, 1, 0.2),  # Controls tick marks on the legend
                labels = seq(-0.2, 1, 0.2)  # Labels match the tick marks
              ),
              # col = colorRampPalette(c("#3B4CC0", "#FFFFFF", "#B40426"))(100),
              # col = viridis(100),
              cluster_rows = F, 
              cluster_columns = F, 
              row_names_side = "left",
              column_names_side = "top", 
              row_order = us_order,
              column_order = group_order,
              # row_labels = gsub("^(.*) (.* .*?|unstim_0h)$", "\\2", rownames(heatmap_mt), perl = T), # "sample2 unstim_0h" to "unstim_0h"
              column_labels = gsub("unstim", "Unstim", gsub(".*\\|\\s*(.*)", "\\1", colnames(heatmap_mt), perl = TRUE)),
              row_names_gp = gpar(fontsize = 10),
              column_names_gp = gpar(fontsize = 10),
              row_split = factor(x_split, levels = names(x_groups)),
              column_split = factor(rep(ct_order, each = 2), levels = ct_order),
              # top_annotation = HeatmapAnnotation(Condition = factor(gsub("^(.*) \\| .*", "\\1", colnames(heatmap_mt)), levels = ct_order),
              #                                    col = list(Condition = selectedCol)),
              
              top_annotation = columnAnnotation(foo = anno_block(gp = gpar(fill = color.order),
                                                                 labels = c("B", "CD4 T", "CD8 T", "NK", "Mono/Mph"),
                                                                 labels_gp = gpar(col = "white", fontsize = 10))),
              
              left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = color.order),
                                                               labels = names(x_groups),
                                                               labels_gp = gpar(col = "white", fontsize = 10))),
              cell_fun = function(j, i, x, y, w, h, fill) {
                if (mat_pval[i, j]>=0.05){
                  grid.rect(x = x, y = y, width = w, height = h, gp = gpar(col = "lightgrey", fill = "lightgrey"))
                  # grid.text(sprintf("%.3f", mat_pval[i, j]), x, y, gp = gpar(fontsize = 3))
                }
              },
              row_title = NULL,
              column_title = NULL,
              row_gap = unit(2,"mm"),
              column_gap = unit(2,"mm")
)

pdf(paste0(deResDir, blResDir, "/bl_shared_log2_sig_correlation_values.pdf"), width = 5.74, height = 5.4)
draw(hm)
dev.off()



# significant & not significant
hm <- Heatmap(heatmap_mt, name = "Pearson's r", na_col = "black",
              # col = colorRamp2(c(-0.1, 0, 0.1, 0.2, 0.3, 0.4), c("#60dbe8","white","#ffdd00","#ff9b54", "#EA4749", "#720026")),
              # col = colorRampPalette(c("#172896","#60dbe8","white","#ffdd00", "#EA4749"))(100),
              col = colorRamp2(c(-0.2,  0, 0.2, 0.4, 0.6, 0.8, 1),
                               c("#60dbe8", "white", 
                                 "#ffdd00", "#ff9b54", 
                                 "#EA4749", "#b2001d", 
                                 "#720026")),
              heatmap_legend_param = list(
                title = "Pearson's r",
                at = seq(-0.2, 1, 0.2),  # Controls tick marks on the legend
                labels = seq(-0.2, 1, 0.2)  # Labels match the tick marks
              ),
              cluster_rows = F, 
              cluster_columns = F, 
              row_names_side = "left",
              column_names_side = "top", 
              row_order = us_order,
              column_order = group_order,
              # row_labels = gsub("^(.*) (.* .*?|unstim_0h)$", "\\2", rownames(heatmap_mt), perl = T), # "sample2 unstim_0h" to "unstim_0h"
              column_labels = gsub("unstim", "Unstim", gsub(".*\\|\\s*(.*)", "\\1", colnames(heatmap_mt), perl = TRUE)),
              row_names_gp = gpar(fontsize = 10),
              column_names_gp = gpar(fontsize = 10),
              row_split = factor(x_split, levels = names(x_groups)),
              column_split = factor(rep(ct_order, each = 2), levels = ct_order),
              # top_annotation = HeatmapAnnotation(Condition = factor(gsub("^(.*) \\| .*", "\\1", colnames(heatmap_mt)), levels = ct_order),
              #                                    col = list(Condition = selectedCol)),
              
              top_annotation = columnAnnotation(foo = anno_block(gp = gpar(fill = color.order),
                                                                 labels = c("B", "CD4 T", "CD8 T", "NK", "Mono/Mph"),
                                                                 labels_gp = gpar(col = "white", fontsize = 10))),
              
              left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = color.order),
                                                               labels = names(x_groups),
                                                               labels_gp = gpar(col = "white", fontsize = 10))),
              # cell_fun = function(j, i, x, y, w, h, fill) {
              #   if (mat_pval[i, j]>=0.05){
              #     grid.rect(x = x, y = y, width = w, height = h, gp = gpar(col = "lightgrey", fill = "lightgrey"))
              #     # grid.text(sprintf("%.3f", mat_pval[i, j]), x, y, gp = gpar(fontsize = 3))
              #   }
              # },
              row_title = NULL,
              column_title = NULL,
              row_gap = unit(2,"mm"),
              column_gap = unit(2,"mm")
)

pdf(paste0(deResDir, blResDir, "/bl_shared_log2_all_correlation_values.pdf"), width = 5.74, height = 5.4)
draw(hm)
dev.off()



# 3.0.2 shared DEGs -> slope -. plot the heatmap ---------------------------------------

# myCol <- colorRampPalette(c("#3B4CC0", "#FFFFFF", "#B40426"))(100)

correlation_values_df <- read.table(paste0(deResDir, blResDir, "/bl_shared_log2_correlation_values.txt"), sep = "\t", header = TRUE)

heatmap_mt <- correlation_values_df[,c(1,2,5)] %>% pivot_wider(names_from = lung, values_from = slope) %>% column_to_rownames(var = "blood") %>% as.matrix()
mat_pval <- correlation_values_df[,c(1,2,4)] %>% pivot_wider(names_from = lung, values_from = p) %>% column_to_rownames(var = "blood") %>% as.matrix()
max(mat_pval, na.rm = TRUE)


ct_order <- c("B", "CD4 T","CD8 T","NK","Mono/Mph")
tm_order <- c("Ig","CD3","CD3","CD3","LPS")
timepoint <- c("4h", "18h")

colnames(heatmap_mt)[1] # follow this format ""B cells | Ig_4h-unstim_0h""

group_order <- c()
for (c in seq_along(ct_order)) {
  group_order <- c(group_order, paste0(ct_order[c], " | ", tm_order[c], "_", timepoint[1], "-unstim_0h"))
  group_order <- c(group_order, paste0(ct_order[c], " | ", tm_order[c], "_", timepoint[2], "-unstim_0h"))
}

group_order <- group_order[group_order %in% colnames(heatmap_mt)] 

# > group_order
# [1] "B cells | Ig_4h-unstim_0h"                "B cells | Ig_18h-unstim_0h"               "CD4 T cells | CD3_4h-unstim_0h"          
# [4] "CD4 T cells | CD3_18h-unstim_0h"          "CD8 T cells | CD3_4h-unstim_0h"           "CD8 T cells | CD3_18h-unstim_0h"         
# [7] "NK cells | CD3_4h-unstim_0h"              "NK cells | CD3_18h-unstim_0h"             "Monocyte-derived Mph | LPS_4h-unstim_0h" 
# [10] "Monocyte-derived Mph | LPS_18h-unstim_0h"


x_groups <- list(
  B = c("Bulk_B_S-Bulk_B_U", "Mem_B_S-Mem_B_U", "Naive_B_S-Naive_B_U"),
  
  CD4T = c("Effector_CD4pos_T_S-Effector_CD4pos_T_U", "Follicular_T_Helper_S-Follicular_T_Helper_U", 
           "Memory_Teffs_S-Memory_Teffs_U", "Naive_Teffs_S-Naive_Teffs_U", 
           "Th1_precursors_S-Th1_precursors_U", "Th17_precursors_S-Th17_precursors_U", 
           "Th2_precursors_S-Th2_precursors_U", "Memory_Tregs_S-Memory_Tregs_U", 
           "Naive_Tregs_S-Naive_Tregs_U", "Regulatory_T_S-Regulatory_T_U"), 
  
  CD8T = c("CD8pos_T_S-CD8pos_T_U", "Central_memory_CD8pos_T_S-Central_memory_CD8pos_T_U", 
           "Effector_memory_CD8pos_T_S-Effector_memory_CD8pos_T_U", "Naive_CD8_T_S-Naive_CD8_T_U",
           "Gamma_delta_T_S-Gamma_delta_T_U"),
  
  NK = c("Mature_NK_S-Mature_NK_U"),
  
  Myeloid = c("Monocytes_S-Monocytes_U")
)



us_order <- c(
  # B 
  "Bulk_B_S-Bulk_B_U", "Mem_B_S-Mem_B_U", "Naive_B_S-Naive_B_U",
  
  # CD4 T
  "Effector_CD4pos_T_S-Effector_CD4pos_T_U", "Follicular_T_Helper_S-Follicular_T_Helper_U", 
  "Memory_Teffs_S-Memory_Teffs_U", "Naive_Teffs_S-Naive_Teffs_U", 
  "Th1_precursors_S-Th1_precursors_U", "Th17_precursors_S-Th17_precursors_U", 
  "Th2_precursors_S-Th2_precursors_U", "Memory_Tregs_S-Memory_Tregs_U", 
  "Naive_Tregs_S-Naive_Tregs_U", "Regulatory_T_S-Regulatory_T_U", 
  
  # CD8 T
  "CD8pos_T_S-CD8pos_T_U", "Central_memory_CD8pos_T_S-Central_memory_CD8pos_T_U", 
  "Effector_memory_CD8pos_T_S-Effector_memory_CD8pos_T_U", "Naive_CD8_T_S-Naive_CD8_T_U", 
  "Gamma_delta_T_S-Gamma_delta_T_U",
  
  # NK
  "Mature_NK_S-Mature_NK_U", 
  
  # Mono Mph
  "Monocytes_S-Monocytes_U")

setdiff(rownames(heatmap_mt), us_order)

# match each rowname with a cell type 
# > x_split
# [1] "CD4T"    "CD4T"    "CD4T"    "CD8T"    "CD4T"    "B"       "B"       "B"       "CD8T"    "Myeloid" "CD4T"    "CD4T"    "CD4T"    "CD4T"    "CD8T"    "CD4T"    "NK"      "CD8T"    "CD4T"   
# [20] "NK"   
idx_list <- lapply(x_groups, function(x) which(rownames(heatmap_mt) %in% x))
x_split <- rep(NA, nrow(heatmap_mt))
for (n in names(idx_list)) {
  x_split[idx_list[[n]]] <- n
}

cell_types <- c(
  "Mono/Mph",
  "B", 
  "CD4 T", 
  "CD8 T", 
  "NK",
  "EC", 
  "DC")


selectedCol<-c( 'cornflowerblue', 
                'mediumvioletred',
                'wheat3', 
                'lightsalmon',  
                "#7CE3D8", 
                'yellowgreen',
                'seagreen')

names(selectedCol) <- cell_types
color.order <- c("mediumvioletred", 
                 'wheat3', 
                 'lightsalmon',  
                 "#7CE3D8", 
                 'cornflowerblue')


min(heatmap_mt)
max(heatmap_mt)


# significant & not significant correlations -> plot slope to see if it is biased
hm <- Heatmap(heatmap_mt, name = "Slope", na_col = "black",
              # col = colorRamp2(c(-0.1, 0, 0.1, 0.2, 0.3, 0.4), c("#60dbe8","white","#ffdd00","#ff9b54", "#EA4749", "#720026")),
              # col = colorRampPalette(c("#172896","#60dbe8","white","#ffdd00", "#EA4749"))(100),
              col = colorRamp2(c(-2, -1, 0, 2, 4, 6, 8),
                               c("#172896", "#60dbe8", "white", "#ffdd00","#ff9b54", "#EA4749", "#720026")),
              # heatmap_legend_param = list(
              #   title = "Slope",
              #   at = seq(-0.1, 0.4, 0.1),  # Controls tick marks on the legend
              #   labels = seq(-0.1, 0.4, 0.1)  # Labels match the tick marks
              # ),
              cluster_rows = F, 
              cluster_columns = F, 
              row_names_side = "left",
              column_names_side = "top", 
              row_order = us_order,
              column_order = group_order,
              # row_labels = gsub("^(.*) (.* .*?|unstim_0h)$", "\\2", rownames(heatmap_mt), perl = T), # "sample2 unstim_0h" to "unstim_0h"
              column_labels = gsub("unstim", "Unstim", gsub(".*\\|\\s*(.*)", "\\1", colnames(heatmap_mt), perl = TRUE)),
              row_names_gp = gpar(fontsize = 10),
              column_names_gp = gpar(fontsize = 10),
              row_split = factor(x_split, levels = names(x_groups)),
              column_split = factor(rep(ct_order, each = 2), levels = ct_order),
              # top_annotation = HeatmapAnnotation(Condition = factor(gsub("^(.*) \\| .*", "\\1", colnames(heatmap_mt)), levels = ct_order),
              #                                    col = list(Condition = selectedCol)),
              
              top_annotation = columnAnnotation(foo = anno_block(gp = gpar(fill = color.order),
                                                                 labels = c("B", "CD4 T", "CD8 T", "NK", "Mono/Mph"),
                                                                 labels_gp = gpar(col = "white", fontsize = 10))),
              
              left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = color.order),
                                                               labels = names(x_groups),
                                                               labels_gp = gpar(col = "white", fontsize = 10))),
              # cell_fun = function(j, i, x, y, w, h, fill) {
              #   if (mat_pval[i, j]>=0.05){
              #     grid.rect(x = x, y = y, width = w, height = h, gp = gpar(col = "lightgrey", fill = "lightgrey"))
              #     # grid.text(sprintf("%.3f", mat_pval[i, j]), x, y, gp = gpar(fontsize = 3))
              #   }
              # },
              row_title = NULL,
              column_title = NULL,
              row_gap = unit(2,"mm"),
              column_gap = unit(2,"mm")
)

pdf(paste0(deResDir, blResDir, "/bl_shared_log2_all_slope_values.pdf"), width = 5.74, height = 5.4)
draw(hm)
dev.off()




# 3.1 complete DEGs -> calculate correlation -----------------------------------
plotscatter <- function(df2plot) {
  # myplot <- ggscatter(
  #   df2plot, x='Group1', y='Group2',
  #   add = "reg.line",
  #   add.params = list(color = "blue", fill = "lightgray"),
  #   conf.int = TRUE, conf.int.level = 0.9) +
  #   # coord_equal(ratio = 1) +
  #   scale_x_continuous(breaks=seq(0,20,1)) +
  #   scale_y_continuous(breaks=seq(0,20,1)) +
  #   labs(title = paste0("scatter ", x), x = x, y = y) +
  #   stat_regline_equation(label.y = 5.4) +
  #   stat_cor(method = "pearson", label.x.npc = "left", label.y = 5, r.accuracy = 0.01)
  # # ggsave(paste0(resDir, "/scatterPlot_wOrgClusterAnnotation/", n, ".jpg"), plot = myplot, width = 4, height = 4, units = "in")
  # # Extract R and p values
  # r_value <- ggplot2::ggplot_build(myplot)$data[[4]]$r
  # p_value <- ggplot2::ggplot_build(myplot)$data[[4]]$p
  
  lm_model <- lm(Group2 ~ Group1, data = df2plot)
  slope <- coef(lm_model)[2]  # Extract slope
  
  # Compute Pearson correlation
  cor_test <- cor.test(df2plot$Group1, df2plot$Group2, method = "pearson")
  r_value <- cor_test$estimate  # Extract Pearson's r
  p_value <- cor_test$p.value  
  return(list("r_value"=r_value, "p_value"=p_value, "slope"=slope))
}

y_genes_groups$cell_type_contrast <- paste(y_genes_groups$cell_type, y_genes_groups$contrast, sep = " | ")

correlation_values_df <- data.frame()

for (y in unique(y_genes_groups$cell_type_contrast)) {
  for (x in unique(bulk_log_exp$contrast)) {
    tryCatch({
      col1 <- y_genes_groups %>%
        filter(cell_type_contrast %in% y) %>%
        dplyr::select(geneName, avg_log2FC, cell_type_contrast) %>%
        column_to_rownames(var = "geneName")
      col2 <- bulk_log_exp %>%
        filter(contrast %in% x, 
               peak_geneName != "") %>% 
        dplyr::select(peak_geneName, logFC, contrast) %>%
        column_to_rownames(var = "peak_geneName")
      df2plot <- merge(col1, col2, by = "row.names", all = TRUE) %>%
        dplyr::select(Row.names, avg_log2FC, logFC)
      
      df2plot[is.na(df2plot)] <- 0
      
      colnames(df2plot) <- c("rowname","Group1","Group2")
      scatter_values <- plotscatter(df2plot)
      
      new_row <- data.frame(blood = x, 
                            lung = y, 
                            r = scatter_values$r_value, 
                            p = scatter_values$p_value, 
                            slope = scatter_values$slope, 
                            blood_degs = dim(col2)[1],
                            lung_degs = dim(col1)[1]
                            )
      correlation_values_df <- rbind(correlation_values_df, new_row)
      
    }, error = function(e) {
      message("Error occurred in compGroup ", x, " & ", y, ": ", e$message)
    })
    next
  }
}

write.table(correlation_values_df, paste0(deResDir, blResDir, "/bl_complete_log2_correlation_values.txt"), sep = "\t", row.name = FALSE, quote = FALSE)



# 3.1.1 complete DEGs -> pearson's r ->  plot the heatmap ---------------------------------------

# myCol <- colorRampPalette(c("#3B4CC0", "#FFFFFF", "#B40426"))(100)

correlation_values_df <- read.table(paste0(deResDir, blResDir, "/bl_complete_log2_correlation_values.txt"), sep = "\t", header = TRUE)

heatmap_mt <- correlation_values_df[,c(1,2,3)] %>% pivot_wider(names_from = lung, values_from = r) %>% column_to_rownames(var = "blood") %>% as.matrix()
mat_pval <- correlation_values_df[,c(1,2,4)] %>% pivot_wider(names_from = lung, values_from = p) %>% column_to_rownames(var = "blood") %>% as.matrix()
max(mat_pval, na.rm = TRUE)


ct_order <- c("B", "CD4 T","CD8 T","NK","Mono/Mph")
tm_order <- c("Ig","CD3","CD3","CD3","LPS")
timepoint <- c("4h", "18h")

colnames(heatmap_mt)[1] # follow this format ""B cells | Ig_4h-unstim_0h""

group_order <- c()
for (c in seq_along(ct_order)) {
  group_order <- c(group_order, paste0(ct_order[c], " | ", tm_order[c], "_", timepoint[1], "-unstim_0h"))
  group_order <- c(group_order, paste0(ct_order[c], " | ", tm_order[c], "_", timepoint[2], "-unstim_0h"))
}

group_order <- group_order[group_order %in% colnames(heatmap_mt)] 

# > group_order
# [1] "B cells | Ig_4h-unstim_0h"                "B cells | Ig_18h-unstim_0h"               "CD4 T cells | CD3_4h-unstim_0h"          
# [4] "CD4 T cells | CD3_18h-unstim_0h"          "CD8 T cells | CD3_4h-unstim_0h"           "CD8 T cells | CD3_18h-unstim_0h"         
# [7] "NK cells | CD3_4h-unstim_0h"              "NK cells | CD3_18h-unstim_0h"             "Monocyte-derived Mph | LPS_4h-unstim_0h" 
# [10] "Monocyte-derived Mph | LPS_18h-unstim_0h"


x_groups <- list(
  B = c("Bulk_B_S-Bulk_B_U", "Mem_B_S-Mem_B_U", "Naive_B_S-Naive_B_U"),
  
  CD4T = c("Effector_CD4pos_T_S-Effector_CD4pos_T_U", "Follicular_T_Helper_S-Follicular_T_Helper_U", 
           "Memory_Teffs_S-Memory_Teffs_U", "Naive_Teffs_S-Naive_Teffs_U", 
           "Th1_precursors_S-Th1_precursors_U", "Th17_precursors_S-Th17_precursors_U", 
           "Th2_precursors_S-Th2_precursors_U", "Memory_Tregs_S-Memory_Tregs_U", 
           "Naive_Tregs_S-Naive_Tregs_U", "Regulatory_T_S-Regulatory_T_U"), 
  
  CD8T = c("CD8pos_T_S-CD8pos_T_U", "Central_memory_CD8pos_T_S-Central_memory_CD8pos_T_U", 
           "Effector_memory_CD8pos_T_S-Effector_memory_CD8pos_T_U", "Naive_CD8_T_S-Naive_CD8_T_U",
           "Gamma_delta_T_S-Gamma_delta_T_U"),
  
  NK = c("Mature_NK_S-Mature_NK_U"),
  
  Myeloid = c("Monocytes_S-Monocytes_U")
)



us_order <- c(
  # B 
  "Bulk_B_S-Bulk_B_U", "Mem_B_S-Mem_B_U", "Naive_B_S-Naive_B_U",
  
  # CD4 T
  "Effector_CD4pos_T_S-Effector_CD4pos_T_U", "Follicular_T_Helper_S-Follicular_T_Helper_U", 
  "Memory_Teffs_S-Memory_Teffs_U", "Naive_Teffs_S-Naive_Teffs_U", 
  "Th1_precursors_S-Th1_precursors_U", "Th17_precursors_S-Th17_precursors_U", 
  "Th2_precursors_S-Th2_precursors_U", "Memory_Tregs_S-Memory_Tregs_U", 
  "Naive_Tregs_S-Naive_Tregs_U", "Regulatory_T_S-Regulatory_T_U", 
  
  # CD8 T
  "CD8pos_T_S-CD8pos_T_U", "Central_memory_CD8pos_T_S-Central_memory_CD8pos_T_U", 
  "Effector_memory_CD8pos_T_S-Effector_memory_CD8pos_T_U", "Naive_CD8_T_S-Naive_CD8_T_U", 
  "Gamma_delta_T_S-Gamma_delta_T_U",
  
  # NK
  "Mature_NK_S-Mature_NK_U", 
  
  # Mono Mph
  "Monocytes_S-Monocytes_U")

setdiff(rownames(heatmap_mt), us_order)

# match each rowname with a cell type 
# > x_split
# [1] "CD4T"    "CD4T"    "CD4T"    "CD8T"    "CD4T"    "B"       "B"       "B"       "CD8T"    "Myeloid" "CD4T"    "CD4T"    "CD4T"    "CD4T"    "CD8T"    "CD4T"    "NK"      "CD8T"    "CD4T"   
# [20] "NK"   
idx_list <- lapply(x_groups, function(x) which(rownames(heatmap_mt) %in% x))
x_split <- rep(NA, nrow(heatmap_mt))
for (n in names(idx_list)) {
  x_split[idx_list[[n]]] <- n
}

cell_types <- c(
  "Mono/Mph",
  "B", 
  "CD4 T", 
  "CD8 T", 
  "NK",
  "EC", 
  "DC")


selectedCol<-c( 'cornflowerblue', 
                'mediumvioletred',
                'wheat3', 
                'lightsalmon',  
                "#7CE3D8", 
                'yellowgreen',
                'seagreen')

names(selectedCol) <- cell_types
color.order <- c("mediumvioletred", 
                 'wheat3', 
                 'lightsalmon',  
                 "#7CE3D8", 
                 'cornflowerblue')

# significant only
hm <- Heatmap(heatmap_mt, name = "Pearson's r", na_col = "black",
              # col = colorRamp2(c(-0.1, 0, 0.1, 0.2, 0.3, 0.4), c("#60dbe8","white","#ffdd00","#ff9b54", "#EA4749", "#720026")),
              # col = colorRampPalette(c("#172896","#60dbe8","white","#ffdd00", "#EA4749"))(100),
              col = colorRamp2(c(-0.1, -0.05, 0, 0.1, 0.2, 0.3, 0.4),
                               c("#172896", "#60dbe8", "white", "#ffdd00","#ff9b54", "#EA4749", "#720026")),
              heatmap_legend_param = list(
                title = "Pearson's r",
                at = seq(-0.1, 0.4, 0.1),  # Controls tick marks on the legend
                labels = seq(-0.1, 0.4, 0.1)  # Labels match the tick marks
              ),
              # col = colorRampPalette(c("#3B4CC0", "#FFFFFF", "#B40426"))(100),
              # col = viridis(100),
              cluster_rows = F, 
              cluster_columns = F, 
              row_names_side = "left",
              column_names_side = "top", 
              row_order = us_order,
              column_order = group_order,
              # row_labels = gsub("^(.*) (.* .*?|unstim_0h)$", "\\2", rownames(heatmap_mt), perl = T), # "sample2 unstim_0h" to "unstim_0h"
              column_labels = gsub("unstim", "Unstim", gsub(".*\\|\\s*(.*)", "\\1", colnames(heatmap_mt), perl = TRUE)),
              row_names_gp = gpar(fontsize = 10),
              column_names_gp = gpar(fontsize = 10),
              row_split = factor(x_split, levels = names(x_groups)),
              column_split = factor(rep(ct_order, each = 2), levels = ct_order),
              # top_annotation = HeatmapAnnotation(Condition = factor(gsub("^(.*) \\| .*", "\\1", colnames(heatmap_mt)), levels = ct_order),
              #                                    col = list(Condition = selectedCol)),
              
              top_annotation = columnAnnotation(foo = anno_block(gp = gpar(fill = color.order),
                                                                 labels = c("B", "CD4 T", "CD8 T", "NK", "Mono/Mph"),
                                                                 labels_gp = gpar(col = "white", fontsize = 10))),
              
              left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = color.order),
                                                               labels = names(x_groups),
                                                               labels_gp = gpar(col = "white", fontsize = 10))),
              cell_fun = function(j, i, x, y, w, h, fill) {
                if (mat_pval[i, j]>=0.05){
                  grid.rect(x = x, y = y, width = w, height = h, gp = gpar(col = "lightgrey", fill = "lightgrey"))
                  # grid.text(sprintf("%.3f", mat_pval[i, j]), x, y, gp = gpar(fontsize = 3))
                }
              },
              row_title = NULL,
              column_title = NULL,
              row_gap = unit(2,"mm"),
              column_gap = unit(2,"mm")
)

pdf(paste0(deResDir, blResDir, "/bl_complete_log2_sig_correlation_values.pdf"), width = 5.74, height = 5.4)
draw(hm)
dev.off()



# significant & not significant
hm <- Heatmap(heatmap_mt, name = "Pearson's r", na_col = "black",
              # col = colorRamp2(c(-0.1, 0, 0.1, 0.2, 0.3, 0.4), c("#60dbe8","white","#ffdd00","#ff9b54", "#EA4749", "#720026")),
              # col = colorRampPalette(c("#172896","#60dbe8","white","#ffdd00", "#EA4749"))(100),
              col = colorRamp2(c(-0.1, -0.05, 0, 0.1, 0.2, 0.3, 0.4),
                               c("#172896", "#60dbe8", "white", "#ffdd00","#ff9b54", "#EA4749", "#720026")),
              heatmap_legend_param = list(
                title = "Pearson's r",
                at = seq(-0.1, 0.4, 0.1),  # Controls tick marks on the legend
                labels = seq(-0.1, 0.4, 0.1)  # Labels match the tick marks
              ),
              cluster_rows = F, 
              cluster_columns = F, 
              row_names_side = "left",
              column_names_side = "top", 
              row_order = us_order,
              column_order = group_order,
              # row_labels = gsub("^(.*) (.* .*?|unstim_0h)$", "\\2", rownames(heatmap_mt), perl = T), # "sample2 unstim_0h" to "unstim_0h"
              column_labels = gsub("unstim", "Unstim", gsub(".*\\|\\s*(.*)", "\\1", colnames(heatmap_mt), perl = TRUE)),
              row_names_gp = gpar(fontsize = 10),
              column_names_gp = gpar(fontsize = 10),
              row_split = factor(x_split, levels = names(x_groups)),
              column_split = factor(rep(ct_order, each = 2), levels = ct_order),
              # top_annotation = HeatmapAnnotation(Condition = factor(gsub("^(.*) \\| .*", "\\1", colnames(heatmap_mt)), levels = ct_order),
              #                                    col = list(Condition = selectedCol)),
              
              top_annotation = columnAnnotation(foo = anno_block(gp = gpar(fill = color.order),
                                                                 labels = c("B", "CD4 T", "CD8 T", "NK", "Mono/Mph"),
                                                                 labels_gp = gpar(col = "white", fontsize = 10))),
              
              left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = color.order),
                                                               labels = names(x_groups),
                                                               labels_gp = gpar(col = "white", fontsize = 10))),
              # cell_fun = function(j, i, x, y, w, h, fill) {
              #   if (mat_pval[i, j]>=0.05){
              #     grid.rect(x = x, y = y, width = w, height = h, gp = gpar(col = "lightgrey", fill = "lightgrey"))
              #     # grid.text(sprintf("%.3f", mat_pval[i, j]), x, y, gp = gpar(fontsize = 3))
              #   }
              # },
              row_title = NULL,
              column_title = NULL,
              row_gap = unit(2,"mm"),
              column_gap = unit(2,"mm")
)

pdf(paste0(deResDir, blResDir, "/bl_complete_log2_all_correlation_values.pdf"), width = 5.74, height = 5.4)
draw(hm)
dev.off()



# 3.1.2 complete DEGs -> slope -. plot the heatmap ---------------------------------------

# myCol <- colorRampPalette(c("#3B4CC0", "#FFFFFF", "#B40426"))(100)

correlation_values_df <- read.table(paste0(deResDir, blResDir, "/bl_complete_log2_correlation_values.txt"), sep = "\t", header = TRUE)

heatmap_mt <- correlation_values_df[,c(1,2,5)] %>% pivot_wider(names_from = lung, values_from = slope) %>% column_to_rownames(var = "blood") %>% as.matrix()
mat_pval <- correlation_values_df[,c(1,2,4)] %>% pivot_wider(names_from = lung, values_from = p) %>% column_to_rownames(var = "blood") %>% as.matrix()
max(mat_pval, na.rm = TRUE)


ct_order <- c("B", "CD4 T","CD8 T","NK","Mono/Mph")
tm_order <- c("Ig","CD3","CD3","CD3","LPS")
timepoint <- c("4h", "18h")

colnames(heatmap_mt)[1] # follow this format ""B cells | Ig_4h-unstim_0h""

group_order <- c()
for (c in seq_along(ct_order)) {
  group_order <- c(group_order, paste0(ct_order[c], " | ", tm_order[c], "_", timepoint[1], "-unstim_0h"))
  group_order <- c(group_order, paste0(ct_order[c], " | ", tm_order[c], "_", timepoint[2], "-unstim_0h"))
}

group_order <- group_order[group_order %in% colnames(heatmap_mt)] 

# > group_order
# [1] "B cells | Ig_4h-unstim_0h"                "B cells | Ig_18h-unstim_0h"               "CD4 T cells | CD3_4h-unstim_0h"          
# [4] "CD4 T cells | CD3_18h-unstim_0h"          "CD8 T cells | CD3_4h-unstim_0h"           "CD8 T cells | CD3_18h-unstim_0h"         
# [7] "NK cells | CD3_4h-unstim_0h"              "NK cells | CD3_18h-unstim_0h"             "Monocyte-derived Mph | LPS_4h-unstim_0h" 
# [10] "Monocyte-derived Mph | LPS_18h-unstim_0h"


x_groups <- list(
  B = c("Bulk_B_S-Bulk_B_U", "Mem_B_S-Mem_B_U", "Naive_B_S-Naive_B_U"),
  
  CD4T = c("Effector_CD4pos_T_S-Effector_CD4pos_T_U", "Follicular_T_Helper_S-Follicular_T_Helper_U", 
           "Memory_Teffs_S-Memory_Teffs_U", "Naive_Teffs_S-Naive_Teffs_U", 
           "Th1_precursors_S-Th1_precursors_U", "Th17_precursors_S-Th17_precursors_U", 
           "Th2_precursors_S-Th2_precursors_U", "Memory_Tregs_S-Memory_Tregs_U", 
           "Naive_Tregs_S-Naive_Tregs_U", "Regulatory_T_S-Regulatory_T_U"), 
  
  CD8T = c("CD8pos_T_S-CD8pos_T_U", "Central_memory_CD8pos_T_S-Central_memory_CD8pos_T_U", 
           "Effector_memory_CD8pos_T_S-Effector_memory_CD8pos_T_U", "Naive_CD8_T_S-Naive_CD8_T_U",
           "Gamma_delta_T_S-Gamma_delta_T_U"),
  
  NK = c("Mature_NK_S-Mature_NK_U"),
  
  Myeloid = c("Monocytes_S-Monocytes_U")
)



us_order <- c(
  # B 
  "Bulk_B_S-Bulk_B_U", "Mem_B_S-Mem_B_U", "Naive_B_S-Naive_B_U",
  
  # CD4 T
  "Effector_CD4pos_T_S-Effector_CD4pos_T_U", "Follicular_T_Helper_S-Follicular_T_Helper_U", 
  "Memory_Teffs_S-Memory_Teffs_U", "Naive_Teffs_S-Naive_Teffs_U", 
  "Th1_precursors_S-Th1_precursors_U", "Th17_precursors_S-Th17_precursors_U", 
  "Th2_precursors_S-Th2_precursors_U", "Memory_Tregs_S-Memory_Tregs_U", 
  "Naive_Tregs_S-Naive_Tregs_U", "Regulatory_T_S-Regulatory_T_U", 
  
  # CD8 T
  "CD8pos_T_S-CD8pos_T_U", "Central_memory_CD8pos_T_S-Central_memory_CD8pos_T_U", 
  "Effector_memory_CD8pos_T_S-Effector_memory_CD8pos_T_U", "Naive_CD8_T_S-Naive_CD8_T_U", 
  "Gamma_delta_T_S-Gamma_delta_T_U",
  
  # NK
  "Mature_NK_S-Mature_NK_U", 
  
  # Mono Mph
  "Monocytes_S-Monocytes_U")

setdiff(rownames(heatmap_mt), us_order)

# match each rowname with a cell type 
# > x_split
# [1] "CD4T"    "CD4T"    "CD4T"    "CD8T"    "CD4T"    "B"       "B"       "B"       "CD8T"    "Myeloid" "CD4T"    "CD4T"    "CD4T"    "CD4T"    "CD8T"    "CD4T"    "NK"      "CD8T"    "CD4T"   
# [20] "NK"   
idx_list <- lapply(x_groups, function(x) which(rownames(heatmap_mt) %in% x))
x_split <- rep(NA, nrow(heatmap_mt))
for (n in names(idx_list)) {
  x_split[idx_list[[n]]] <- n
}

cell_types <- c(
  "Mono/Mph",
  "B", 
  "CD4 T", 
  "CD8 T", 
  "NK",
  "EC", 
  "DC")


selectedCol<-c( 'cornflowerblue', 
                'mediumvioletred',
                'wheat3', 
                'lightsalmon',  
                "#7CE3D8", 
                'yellowgreen',
                'seagreen')

names(selectedCol) <- cell_types
color.order <- c("mediumvioletred", 
                 'wheat3', 
                 'lightsalmon',  
                 "#7CE3D8", 
                 'cornflowerblue')


min(heatmap_mt)
max(heatmap_mt)


# significant & not significant correlations -> plot slope to see if it is biased
hm <- Heatmap(heatmap_mt, name = "Slope", na_col = "black",
              # col = colorRamp2(c(-0.1, 0, 0.1, 0.2, 0.3, 0.4), c("#60dbe8","white","#ffdd00","#ff9b54", "#EA4749", "#720026")),
              col = colorRampPalette(c("#172896","#60dbe8","white","#ffdd00", "#EA4749"))(100),
              # col = colorRamp2(c(-0.1, -0.05, 0, 0.1, 0.2, 0.3, 0.4),
              #                  c("#172896", "#60dbe8", "white", "#ffdd00","#ff9b54", "#EA4749", "#720026")),
              # heatmap_legend_param = list(
              #   title = "Slope",
              #   at = seq(-0.1, 0.4, 0.1),  # Controls tick marks on the legend
              #   labels = seq(-0.1, 0.4, 0.1)  # Labels match the tick marks
              # ),
              cluster_rows = F, 
              cluster_columns = F, 
              row_names_side = "left",
              column_names_side = "top", 
              row_order = us_order,
              column_order = group_order,
              # row_labels = gsub("^(.*) (.* .*?|unstim_0h)$", "\\2", rownames(heatmap_mt), perl = T), # "sample2 unstim_0h" to "unstim_0h"
              column_labels = gsub("unstim", "Unstim", gsub(".*\\|\\s*(.*)", "\\1", colnames(heatmap_mt), perl = TRUE)),
              row_names_gp = gpar(fontsize = 10),
              column_names_gp = gpar(fontsize = 10),
              row_split = factor(x_split, levels = names(x_groups)),
              column_split = factor(rep(ct_order, each = 2), levels = ct_order),
              # top_annotation = HeatmapAnnotation(Condition = factor(gsub("^(.*) \\| .*", "\\1", colnames(heatmap_mt)), levels = ct_order),
              #                                    col = list(Condition = selectedCol)),
              
              top_annotation = columnAnnotation(foo = anno_block(gp = gpar(fill = color.order),
                                                                 labels = c("B", "CD4 T", "CD8 T", "NK", "Mono/Mph"),
                                                                 labels_gp = gpar(col = "white", fontsize = 10))),
              
              left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = color.order),
                                                               labels = names(x_groups),
                                                               labels_gp = gpar(col = "white", fontsize = 10))),
              # cell_fun = function(j, i, x, y, w, h, fill) {
              #   if (mat_pval[i, j]>=0.05){
              #     grid.rect(x = x, y = y, width = w, height = h, gp = gpar(col = "lightgrey", fill = "lightgrey"))
              #     # grid.text(sprintf("%.3f", mat_pval[i, j]), x, y, gp = gpar(fontsize = 3))
              #   }
              # },
              row_title = NULL,
              column_title = NULL,
              row_gap = unit(2,"mm"),
              column_gap = unit(2,"mm")
)

pdf(paste0(deResDir, blResDir, "/bl_complete_log2_all_slope_values.pdf"), width = 5.74, height = 5.4)
draw(hm)
dev.off()




# 3.2 focus on blood DEGs -> find if there is any in lung -> calculate correlation -----------------------------------
plotscatter <- function(df2plot) {
  # myplot <- ggscatter(
  #   df2plot, x='Group1', y='Group2',
  #   add = "reg.line",
  #   add.params = list(color = "blue", fill = "lightgray"),
  #   conf.int = TRUE, conf.int.level = 0.9) +
  #   # coord_equal(ratio = 1) +
  #   scale_x_continuous(breaks=seq(0,20,1)) +
  #   scale_y_continuous(breaks=seq(0,20,1)) +
  #   labs(title = paste0("scatter ", x), x = x, y = y) +
  #   stat_regline_equation(label.y = 5.4) +
  #   stat_cor(method = "pearson", label.x.npc = "left", label.y = 5, r.accuracy = 0.01)
  # # ggsave(paste0(resDir, "/scatterPlot_wOrgClusterAnnotation/", n, ".jpg"), plot = myplot, width = 4, height = 4, units = "in")
  # # Extract R and p values
  # r_value <- ggplot2::ggplot_build(myplot)$data[[4]]$r
  # p_value <- ggplot2::ggplot_build(myplot)$data[[4]]$p
  
  lm_model <- lm(Group2 ~ Group1, data = df2plot)
  slope <- coef(lm_model)[2]  # Extract slope
  
  # Compute Pearson correlation
  cor_test <- cor.test(df2plot$Group1, df2plot$Group2, method = "pearson")
  r_value <- cor_test$estimate  # Extract Pearson's r
  p_value <- cor_test$p.value  
  return(list("r_value"=r_value, "p_value"=p_value, "slope"=slope))
}

y_genes_groups$cell_type_contrast <- paste(y_genes_groups$cell_type, y_genes_groups$contrast, sep = " | ")

correlation_values_df <- data.frame()

for (y in unique(y_genes_groups$cell_type_contrast)) {
  for (x in unique(bulk_log_exp$contrast)) {
    tryCatch({
      col1 <- y_genes_groups %>%
        filter(cell_type_contrast %in% y) %>%
        dplyr::select(geneName, avg_log2FC, cell_type_contrast) %>%
        column_to_rownames(var = "geneName")
      col2 <- bulk_log_exp %>%
        filter(contrast %in% x, 
               peak_geneName != "") %>% 
        dplyr::select(peak_geneName, logFC, contrast) %>%
        column_to_rownames(var = "peak_geneName")
      df2plot <- merge(col1, col2, by = "row.names", all.y = TRUE) %>%
        dplyr::select(Row.names, avg_log2FC, logFC)
      
      df2plot[is.na(df2plot)] <- 0
      
      colnames(df2plot) <- c("rowname","Group1","Group2")
      scatter_values <- plotscatter(df2plot)
      
      new_row <- data.frame(blood = x, 
                            lung = y, 
                            r = scatter_values$r_value, 
                            p = scatter_values$p_value, 
                            slope = scatter_values$slope, 
                            blood_degs = dim(col2)[1],
                            lung_degs = dim(col1)[1]
      )
      correlation_values_df <- rbind(correlation_values_df, new_row)
      
    }, error = function(e) {
      message("Error occurred in compGroup ", x, " & ", y, ": ", e$message)
    })
    next
  }
}

write.table(correlation_values_df, paste0(deResDir, blResDir, "/bl_log2_correlation_values.txt"), sep = "\t", row.name = FALSE, quote = FALSE)



# 3.0.1 shared DEGs -> pearson's r ->  plot the heatmap ---------------------------------------

# myCol <- colorRampPalette(c("#3B4CC0", "#FFFFFF", "#B40426"))(100)

correlation_values_df <- read.table(paste0(deResDir, blResDir, "/bl_log2_correlation_values.txt"), sep = "\t", header = TRUE)

heatmap_mt <- correlation_values_df[,c(1,2,3)] %>% pivot_wider(names_from = lung, values_from = r) %>% column_to_rownames(var = "blood") %>% as.matrix()
mat_pval <- correlation_values_df[,c(1,2,4)] %>% pivot_wider(names_from = lung, values_from = p) %>% column_to_rownames(var = "blood") %>% as.matrix()
max(mat_pval, na.rm = TRUE)


ct_order <- c("B", "CD4 T","CD8 T","NK","Mono/Mph")
tm_order <- c("Ig","CD3","CD3","CD3","LPS")
timepoint <- c("4h", "18h")

colnames(heatmap_mt)[1] # follow this format ""B cells | Ig_4h-unstim_0h""

group_order <- c()
for (c in seq_along(ct_order)) {
  group_order <- c(group_order, paste0(ct_order[c], " | ", tm_order[c], "_", timepoint[1], "-unstim_0h"))
  group_order <- c(group_order, paste0(ct_order[c], " | ", tm_order[c], "_", timepoint[2], "-unstim_0h"))
}

group_order <- group_order[group_order %in% colnames(heatmap_mt)] 

# > group_order
# [1] "B cells | Ig_4h-unstim_0h"                "B cells | Ig_18h-unstim_0h"               "CD4 T cells | CD3_4h-unstim_0h"          
# [4] "CD4 T cells | CD3_18h-unstim_0h"          "CD8 T cells | CD3_4h-unstim_0h"           "CD8 T cells | CD3_18h-unstim_0h"         
# [7] "NK cells | CD3_4h-unstim_0h"              "NK cells | CD3_18h-unstim_0h"             "Monocyte-derived Mph | LPS_4h-unstim_0h" 
# [10] "Monocyte-derived Mph | LPS_18h-unstim_0h"


x_groups <- list(
  B = c("Bulk_B_S-Bulk_B_U", "Mem_B_S-Mem_B_U", "Naive_B_S-Naive_B_U"),
  
  CD4T = c("Effector_CD4pos_T_S-Effector_CD4pos_T_U", "Follicular_T_Helper_S-Follicular_T_Helper_U", 
           "Memory_Teffs_S-Memory_Teffs_U", "Naive_Teffs_S-Naive_Teffs_U", 
           "Th1_precursors_S-Th1_precursors_U", "Th17_precursors_S-Th17_precursors_U", 
           "Th2_precursors_S-Th2_precursors_U", "Memory_Tregs_S-Memory_Tregs_U", 
           "Naive_Tregs_S-Naive_Tregs_U", "Regulatory_T_S-Regulatory_T_U"), 
  
  CD8T = c("CD8pos_T_S-CD8pos_T_U", "Central_memory_CD8pos_T_S-Central_memory_CD8pos_T_U", 
           "Effector_memory_CD8pos_T_S-Effector_memory_CD8pos_T_U", "Naive_CD8_T_S-Naive_CD8_T_U",
           "Gamma_delta_T_S-Gamma_delta_T_U"),
  
  NK = c("Mature_NK_S-Mature_NK_U"),
  
  Myeloid = c("Monocytes_S-Monocytes_U")
)



us_order <- c(
  # B 
  "Bulk_B_S-Bulk_B_U", "Mem_B_S-Mem_B_U", "Naive_B_S-Naive_B_U",
  
  # CD4 T
  "Effector_CD4pos_T_S-Effector_CD4pos_T_U", "Follicular_T_Helper_S-Follicular_T_Helper_U", 
  "Memory_Teffs_S-Memory_Teffs_U", "Naive_Teffs_S-Naive_Teffs_U", 
  "Th1_precursors_S-Th1_precursors_U", "Th17_precursors_S-Th17_precursors_U", 
  "Th2_precursors_S-Th2_precursors_U", "Memory_Tregs_S-Memory_Tregs_U", 
  "Naive_Tregs_S-Naive_Tregs_U", "Regulatory_T_S-Regulatory_T_U", 
  
  # CD8 T
  "CD8pos_T_S-CD8pos_T_U", "Central_memory_CD8pos_T_S-Central_memory_CD8pos_T_U", 
  "Effector_memory_CD8pos_T_S-Effector_memory_CD8pos_T_U", "Naive_CD8_T_S-Naive_CD8_T_U", 
  "Gamma_delta_T_S-Gamma_delta_T_U",
  
  # NK
  "Mature_NK_S-Mature_NK_U", 
  
  # Mono Mph
  "Monocytes_S-Monocytes_U")

setdiff(rownames(heatmap_mt), us_order)

# match each rowname with a cell type 
# > x_split
# [1] "CD4T"    "CD4T"    "CD4T"    "CD8T"    "CD4T"    "B"       "B"       "B"       "CD8T"    "Myeloid" "CD4T"    "CD4T"    "CD4T"    "CD4T"    "CD8T"    "CD4T"    "NK"      "CD8T"    "CD4T"   
# [20] "NK"   
idx_list <- lapply(x_groups, function(x) which(rownames(heatmap_mt) %in% x))
x_split <- rep(NA, nrow(heatmap_mt))
for (n in names(idx_list)) {
  x_split[idx_list[[n]]] <- n
}

cell_types <- c(
  "Mono/Mph",
  "B", 
  "CD4 T", 
  "CD8 T", 
  "NK",
  "EC", 
  "DC")


selectedCol<-c( 'cornflowerblue', 
                'mediumvioletred',
                'wheat3', 
                'lightsalmon',  
                "#7CE3D8", 
                'yellowgreen',
                'seagreen')

names(selectedCol) <- cell_types
color.order <- c("mediumvioletred", 
                 'wheat3', 
                 'lightsalmon',  
                 "#7CE3D8", 
                 'cornflowerblue')

# significant only
hm <- Heatmap(heatmap_mt, name = "Pearson's r", na_col = "black",
              # col = colorRamp2(c(-0.1, 0, 0.1, 0.2, 0.3, 0.4), c("#60dbe8","white","#ffdd00","#ff9b54", "#EA4749", "#720026")),
              # col = colorRampPalette(c("#172896","#60dbe8","white","#ffdd00", "#EA4749"))(100),
              col = colorRamp2(c(-0.2, 0, 0.2, 0.4, 0.6),
                               c("#60dbe8", 
                                 "white", 
                                 "#ffdd00", 
                                 "#ff9b54", 
                                 # "#EA4749", 
                                 "#b2001d")),
              heatmap_legend_param = list(
                title = "Pearson's r",
                at = seq(-0.2, 0.6, 0.2),  # Controls tick marks on the legend
                labels = seq(-0.2, 0.6, 0.2)  # Labels match the tick marks
              ),
              # col = colorRampPalette(c("#3B4CC0", "#FFFFFF", "#B40426"))(100),
              # col = viridis(100),
              cluster_rows = F, 
              cluster_columns = F, 
              row_names_side = "left",
              column_names_side = "top", 
              row_order = us_order,
              column_order = group_order,
              # row_labels = gsub("^(.*) (.* .*?|unstim_0h)$", "\\2", rownames(heatmap_mt), perl = T), # "sample2 unstim_0h" to "unstim_0h"
              column_labels = gsub("unstim", "Unstim", gsub(".*\\|\\s*(.*)", "\\1", colnames(heatmap_mt), perl = TRUE)),
              row_names_gp = gpar(fontsize = 10),
              column_names_gp = gpar(fontsize = 10),
              row_split = factor(x_split, levels = names(x_groups)),
              column_split = factor(rep(ct_order, each = 2), levels = ct_order),
              # top_annotation = HeatmapAnnotation(Condition = factor(gsub("^(.*) \\| .*", "\\1", colnames(heatmap_mt)), levels = ct_order),
              #                                    col = list(Condition = selectedCol)),
              
              top_annotation = columnAnnotation(foo = anno_block(gp = gpar(fill = color.order),
                                                                 labels = c("B", "CD4 T", "CD8 T", "NK", "Mono/Mph"),
                                                                 labels_gp = gpar(col = "white", fontsize = 10))),
              
              left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = color.order),
                                                               labels = names(x_groups),
                                                               labels_gp = gpar(col = "white", fontsize = 10))),
              cell_fun = function(j, i, x, y, w, h, fill) {
                if (mat_pval[i, j]>=0.05){
                  grid.rect(x = x, y = y, width = w, height = h, gp = gpar(col = "lightgrey", fill = "lightgrey"))
                  # grid.text(sprintf("%.3f", mat_pval[i, j]), x, y, gp = gpar(fontsize = 3))
                }
              },
              row_title = NULL,
              column_title = NULL,
              row_gap = unit(2,"mm"),
              column_gap = unit(2,"mm")
)

pdf(paste0(deResDir, blResDir, "/bl_log2_sig_correlation_values.pdf"), width = 5.74, height = 5.4)
draw(hm)
dev.off()



# significant & not significant
hm <- Heatmap(heatmap_mt, name = "Pearson's r", na_col = "black",
              # col = colorRamp2(c(-0.1, 0, 0.1, 0.2, 0.3, 0.4), c("#60dbe8","white","#ffdd00","#ff9b54", "#EA4749", "#720026")),
              # col = colorRampPalette(c("#172896","#60dbe8","white","#ffdd00", "#EA4749"))(100),
              col = colorRamp2(c(-0.2, 0, 0.2, 0.4, 0.6),
                               c("#60dbe8", 
                                 "white", 
                                 "#ffdd00", 
                                 "#ff9b54",
                                 # "#EA4749",
                                 "#b2001d"
                                 )),
              heatmap_legend_param = list(
                title = "Pearson's r",
                at = seq(-0.2, 0.6, 0.2),  # Controls tick marks on the legend
                labels = seq(-0.2, 0.6, 0.2)  # Labels match the tick marks
              ),
              cluster_rows = F, 
              cluster_columns = F, 
              row_names_side = "left",
              column_names_side = "top", 
              row_order = us_order,
              column_order = group_order,
              # row_labels = gsub("^(.*) (.* .*?|unstim_0h)$", "\\2", rownames(heatmap_mt), perl = T), # "sample2 unstim_0h" to "unstim_0h"
              column_labels = gsub("unstim", "Unstim", gsub(".*\\|\\s*(.*)", "\\1", colnames(heatmap_mt), perl = TRUE)),
              row_names_gp = gpar(fontsize = 10),
              column_names_gp = gpar(fontsize = 10),
              row_split = factor(x_split, levels = names(x_groups)),
              column_split = factor(rep(ct_order, each = 2), levels = ct_order),
              # top_annotation = HeatmapAnnotation(Condition = factor(gsub("^(.*) \\| .*", "\\1", colnames(heatmap_mt)), levels = ct_order),
              #                                    col = list(Condition = selectedCol)),
              
              top_annotation = columnAnnotation(foo = anno_block(gp = gpar(fill = color.order),
                                                                 labels = c("B", "CD4 T", "CD8 T", "NK", "Mono/Mph"),
                                                                 labels_gp = gpar(col = "white", fontsize = 10))),
              
              left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = color.order),
                                                               labels = names(x_groups),
                                                               labels_gp = gpar(col = "white", fontsize = 10))),
              # cell_fun = function(j, i, x, y, w, h, fill) {
              #   if (mat_pval[i, j]>=0.05){
              #     grid.rect(x = x, y = y, width = w, height = h, gp = gpar(col = "lightgrey", fill = "lightgrey"))
              #     # grid.text(sprintf("%.3f", mat_pval[i, j]), x, y, gp = gpar(fontsize = 3))
              #   }
              # },
              row_title = NULL,
              column_title = NULL,
              row_gap = unit(2,"mm"),
              column_gap = unit(2,"mm")
)

pdf(paste0(deResDir, blResDir, "/bl_log2_all_correlation_values.pdf"), width = 5.74, height = 5.4)
draw(hm)
dev.off()



# 3.0.2 shared DEGs -> slope -. plot the heatmap ---------------------------------------

# myCol <- colorRampPalette(c("#3B4CC0", "#FFFFFF", "#B40426"))(100)

correlation_values_df <- read.table(paste0(deResDir, blResDir, "/bl_shared_log2_correlation_values.txt"), sep = "\t", header = TRUE)

heatmap_mt <- correlation_values_df[,c(1,2,5)] %>% pivot_wider(names_from = lung, values_from = slope) %>% column_to_rownames(var = "blood") %>% as.matrix()
mat_pval <- correlation_values_df[,c(1,2,4)] %>% pivot_wider(names_from = lung, values_from = p) %>% column_to_rownames(var = "blood") %>% as.matrix()
max(mat_pval, na.rm = TRUE)


ct_order <- c("B", "CD4 T","CD8 T","NK","Mono/Mph")
tm_order <- c("Ig","CD3","CD3","CD3","LPS")
timepoint <- c("4h", "18h")

colnames(heatmap_mt)[1] # follow this format ""B cells | Ig_4h-unstim_0h""

group_order <- c()
for (c in seq_along(ct_order)) {
  group_order <- c(group_order, paste0(ct_order[c], " | ", tm_order[c], "_", timepoint[1], "-unstim_0h"))
  group_order <- c(group_order, paste0(ct_order[c], " | ", tm_order[c], "_", timepoint[2], "-unstim_0h"))
}

group_order <- group_order[group_order %in% colnames(heatmap_mt)] 

# > group_order
# [1] "B cells | Ig_4h-unstim_0h"                "B cells | Ig_18h-unstim_0h"               "CD4 T cells | CD3_4h-unstim_0h"          
# [4] "CD4 T cells | CD3_18h-unstim_0h"          "CD8 T cells | CD3_4h-unstim_0h"           "CD8 T cells | CD3_18h-unstim_0h"         
# [7] "NK cells | CD3_4h-unstim_0h"              "NK cells | CD3_18h-unstim_0h"             "Monocyte-derived Mph | LPS_4h-unstim_0h" 
# [10] "Monocyte-derived Mph | LPS_18h-unstim_0h"


x_groups <- list(
  B = c("Bulk_B_S-Bulk_B_U", "Mem_B_S-Mem_B_U", "Naive_B_S-Naive_B_U"),
  
  CD4T = c("Effector_CD4pos_T_S-Effector_CD4pos_T_U", "Follicular_T_Helper_S-Follicular_T_Helper_U", 
           "Memory_Teffs_S-Memory_Teffs_U", "Naive_Teffs_S-Naive_Teffs_U", 
           "Th1_precursors_S-Th1_precursors_U", "Th17_precursors_S-Th17_precursors_U", 
           "Th2_precursors_S-Th2_precursors_U", "Memory_Tregs_S-Memory_Tregs_U", 
           "Naive_Tregs_S-Naive_Tregs_U", "Regulatory_T_S-Regulatory_T_U"), 
  
  CD8T = c("CD8pos_T_S-CD8pos_T_U", "Central_memory_CD8pos_T_S-Central_memory_CD8pos_T_U", 
           "Effector_memory_CD8pos_T_S-Effector_memory_CD8pos_T_U", "Naive_CD8_T_S-Naive_CD8_T_U",
           "Gamma_delta_T_S-Gamma_delta_T_U"),
  
  NK = c("Mature_NK_S-Mature_NK_U"),
  
  Myeloid = c("Monocytes_S-Monocytes_U")
)



us_order <- c(
  # B 
  "Bulk_B_S-Bulk_B_U", "Mem_B_S-Mem_B_U", "Naive_B_S-Naive_B_U",
  
  # CD4 T
  "Effector_CD4pos_T_S-Effector_CD4pos_T_U", "Follicular_T_Helper_S-Follicular_T_Helper_U", 
  "Memory_Teffs_S-Memory_Teffs_U", "Naive_Teffs_S-Naive_Teffs_U", 
  "Th1_precursors_S-Th1_precursors_U", "Th17_precursors_S-Th17_precursors_U", 
  "Th2_precursors_S-Th2_precursors_U", "Memory_Tregs_S-Memory_Tregs_U", 
  "Naive_Tregs_S-Naive_Tregs_U", "Regulatory_T_S-Regulatory_T_U", 
  
  # CD8 T
  "CD8pos_T_S-CD8pos_T_U", "Central_memory_CD8pos_T_S-Central_memory_CD8pos_T_U", 
  "Effector_memory_CD8pos_T_S-Effector_memory_CD8pos_T_U", "Naive_CD8_T_S-Naive_CD8_T_U", 
  "Gamma_delta_T_S-Gamma_delta_T_U",
  
  # NK
  "Mature_NK_S-Mature_NK_U", 
  
  # Mono Mph
  "Monocytes_S-Monocytes_U")

setdiff(rownames(heatmap_mt), us_order)

# match each rowname with a cell type 
# > x_split
# [1] "CD4T"    "CD4T"    "CD4T"    "CD8T"    "CD4T"    "B"       "B"       "B"       "CD8T"    "Myeloid" "CD4T"    "CD4T"    "CD4T"    "CD4T"    "CD8T"    "CD4T"    "NK"      "CD8T"    "CD4T"   
# [20] "NK"   
idx_list <- lapply(x_groups, function(x) which(rownames(heatmap_mt) %in% x))
x_split <- rep(NA, nrow(heatmap_mt))
for (n in names(idx_list)) {
  x_split[idx_list[[n]]] <- n
}

cell_types <- c(
  "Mono/Mph",
  "B", 
  "CD4 T", 
  "CD8 T", 
  "NK",
  "EC", 
  "DC")


selectedCol<-c( 'cornflowerblue', 
                'mediumvioletred',
                'wheat3', 
                'lightsalmon',  
                "#7CE3D8", 
                'yellowgreen',
                'seagreen')

names(selectedCol) <- cell_types
color.order <- c("mediumvioletred", 
                 'wheat3', 
                 'lightsalmon',  
                 "#7CE3D8", 
                 'cornflowerblue')


min(heatmap_mt)
max(heatmap_mt)


# significant & not significant correlations -> plot slope to see if it is biased
hm <- Heatmap(heatmap_mt, name = "Slope", na_col = "black",
              # col = colorRamp2(c(-0.1, 0, 0.1, 0.2, 0.3, 0.4), c("#60dbe8","white","#ffdd00","#ff9b54", "#EA4749", "#720026")),
              # col = colorRampPalette(c("#172896","#60dbe8","white","#ffdd00", "#EA4749"))(100),
              col = colorRamp2(c(-2, -1, 0, 2, 4, 6, 8),
                               c("#172896", "#60dbe8", "white", "#ffdd00","#ff9b54", "#EA4749", "#720026")),
              # heatmap_legend_param = list(
              #   title = "Slope",
              #   at = seq(-0.1, 0.4, 0.1),  # Controls tick marks on the legend
              #   labels = seq(-0.1, 0.4, 0.1)  # Labels match the tick marks
              # ),
              cluster_rows = F, 
              cluster_columns = F, 
              row_names_side = "left",
              column_names_side = "top", 
              row_order = us_order,
              column_order = group_order,
              # row_labels = gsub("^(.*) (.* .*?|unstim_0h)$", "\\2", rownames(heatmap_mt), perl = T), # "sample2 unstim_0h" to "unstim_0h"
              column_labels = gsub("unstim", "Unstim", gsub(".*\\|\\s*(.*)", "\\1", colnames(heatmap_mt), perl = TRUE)),
              row_names_gp = gpar(fontsize = 10),
              column_names_gp = gpar(fontsize = 10),
              row_split = factor(x_split, levels = names(x_groups)),
              column_split = factor(rep(ct_order, each = 2), levels = ct_order),
              # top_annotation = HeatmapAnnotation(Condition = factor(gsub("^(.*) \\| .*", "\\1", colnames(heatmap_mt)), levels = ct_order),
              #                                    col = list(Condition = selectedCol)),
              
              top_annotation = columnAnnotation(foo = anno_block(gp = gpar(fill = color.order),
                                                                 labels = c("B", "CD4 T", "CD8 T", "NK", "Mono/Mph"),
                                                                 labels_gp = gpar(col = "white", fontsize = 10))),
              
              left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = color.order),
                                                               labels = names(x_groups),
                                                               labels_gp = gpar(col = "white", fontsize = 10))),
              # cell_fun = function(j, i, x, y, w, h, fill) {
              #   if (mat_pval[i, j]>=0.05){
              #     grid.rect(x = x, y = y, width = w, height = h, gp = gpar(col = "lightgrey", fill = "lightgrey"))
              #     # grid.text(sprintf("%.3f", mat_pval[i, j]), x, y, gp = gpar(fontsize = 3))
              #   }
              # },
              row_title = NULL,
              column_title = NULL,
              row_gap = unit(2,"mm"),
              column_gap = unit(2,"mm")
)

pdf(paste0(deResDir, blResDir, "/bl_shared_log2_all_slope_values.pdf"), width = 5.74, height = 5.4)
draw(hm)
dev.off()




## -----------------
## -----------------
## -----------------
## -----------------
## -----------------
## -----------------
## -----------------
## -----------------
###################################################################################### 
###################################################################################### 
#### Below is the orginal version, not considering log 2 scale of gene expression ####
###################################################################################### 
###################################################################################### 
## -----------------

# 1.0 prepare input for bulk RNAseq dataset -----------------------------------

bulk_exp <- read.table(paste0(deResDir, blResDir, "GSE118165_RNA_gene_abundance.txt"), header = T, sep = '\t')

bulk_exp <- setDT(bulk_exp, keep.rownames = "ID")[]

bulk_exp <- as.data.frame(bulk_exp)  # Convert back to data frame

dim(bulk_exp)

# map gene names to ensembl IDs in the gene abundance table (gene of NAs & duplicated genes were removed)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes <- c("ensembl_gene_id", "external_gene_name")  # Attributes to retrieve
gene_names <- getBM(attributes = attributes, filters = "ensembl_gene_id", values = bulk_exp$ID, mart = mart)
colnames(gene_names) <- c("ID", "geneName")
bulk_exp <- merge(bulk_exp, gene_names, by = "ID")
bulk_exp <- bulk_exp[bulk_exp$geneName != "", ]
bulk_exp <- bulk_exp[!duplicated(bulk_exp$geneName),]
rownames(bulk_exp) <- bulk_exp[["geneName"]]
idx2rm <- which(colnames(bulk_exp) %in% c("ID","geneName"))
bulk_exp <- bulk_exp[,-idx2rm]

library(xlsx)
library(readxl)
library(openxlsx)
library(writexl)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(patchwork)
library(grid)
library(gtable)
library(ggpubr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(GetoptLong)
library(preprocessCore)

# visualize exp matrix in box plots
long_df <- bulk_exp %>% rownames_to_column(var = "Gene") %>% melt(id.vars = "Gene")
# bp.before <- ggplot(long_df, aes(x = variable, y = value)) + 
#   geom_boxplot() +
#   scale_y_log10() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   labs(x = "Sample", y = "Expression Level", title = "Gene Expression Across Samples")
# 
# quantile normalize
bulk_exp_matrix <- as.matrix(bulk_exp)
normalized_matrix <- normalize.quantiles(bulk_exp_matrix)
rownames(normalized_matrix) <- rownames(bulk_exp)
colnames(normalized_matrix) <- colnames(bulk_exp)
bulk_exp_norm <- as.data.frame(normalized_matrix)
long_df <- bulk_exp_norm %>% rownames_to_column(var = "Gene") %>% melt(id.vars = "Gene")
# bp.after <- ggplot(long_df, aes(x = variable, y = value)) + 
#   geom_boxplot() +
#   scale_y_log10() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   labs(x = "Sample", y = "Expression Level", title = "Gene Expression Across Samples")
# bp.comp <- bp.before/bp.after
# ggsave(paste0(resDir, "/exp_norm_boxplot.png"), bp.comp, units = "in", width = 20, height = 10)

# calculate avg. gene expression for U/S respectively
grp_df <- data.frame(
  colname = colnames(bulk_exp_norm),
  group = gsub("^.*?\\.", "", colnames(bulk_exp_norm)),
  celltype = gsub("^.*?\\.(.*?)\\..*", "\\1", colnames(bulk_exp_norm)),
  var = gsub("^.*\\.", "", colnames(bulk_exp_norm))
)

# > View(bulk_exp_norm)
# > which(grp_df$group == g)
# [1] 161
# > g
# [1] "Immature_NK.U"
# > colnames(bulk_exp_norm)[161]
# [1] "X1010.Immature_NK.U"

grp_avg_exp <- list()
for (g in unique(grp_df$group)) {
  grp_idx <- which(grp_df$group == g)
  if (is.null(dim(bulk_exp_norm[, grp_idx]))) {
    grp_avg_exp[[g]] <- bulk_exp_norm[, grp_idx]
  }else{
    grp_avg_exp[[g]] <- rowMeans(bulk_exp_norm[, grp_idx])
  }
}
grp_avg_exp_df <- data.frame(grp_avg_exp)
long_df <- grp_avg_exp_df %>% rownames_to_column(var = "Gene") %>% melt(id.vars = "Gene")
# bp.avg <- ggplot(long_df, aes(x = variable, y = value)) + 
#   geom_boxplot() +
#   scale_y_log10() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   labs(x = "Sample", y = "Expression Level", title = "Gene Expression Across Samples")
# ggsave(paste0(resDir, "/exp_avg_boxplot.png"), bp.avg, units = "in", width = 10, height = 6)


# 2.0 calculate correlation between our dataset and GSE118165 for DE gene unions -------------------------------------



# get a specific set of deg results (gene name, log 2, p etc) from_getClusterExpCondDe_out_folder
get_degs <- function(excel.dir, compGroups, ident) {
  degs <- list()
  for (cp in compGroups) { #for each condition
    tryCatch({
      file.all <- list.files(excel.dir)
      file <- file.all[grepl(paste0(cp, "_adjSig_allClusters"), file.all, ignore.case = TRUE)]
      celltype <- gsub("\\W", "", ident)
      degs[[cp]] <- read_excel(paste0(excel.dir, file), sheet = paste0("cluster_",celltype))
    }, error = function(e) {
      message("Error occurred in compGroup ", n, ": ", e$message)
    })
    next
  }
  return(degs)
}

plotscatter <- function(df2plot) {
  myplot <- ggscatter(
    df2plot, x='Group1', y='Group2',
    add = "reg.line",
    add.params = list(color = "blue", fill = "lightgray"),
    conf.int = TRUE, conf.int.level = 0.9) +
    # coord_equal(ratio = 1) +
    scale_x_continuous(breaks=seq(0,20,1)) +
    scale_y_continuous(breaks=seq(0,20,1)) +
    labs(title = paste0("scatter ", n), x = n, y = m) +
    stat_regline_equation(label.y = 5.4) +
    stat_cor(method = "pearson", label.x.npc = "left", label.y = 5, r.accuracy = 0.01)
  # ggsave(paste0(resDir, "/scatterPlot_wOrgClusterAnnotation/", n, ".jpg"), plot = myplot, width = 4, height = 4, units = "in")
  # Extract R and p values
  r_value <- ggplot2::ggplot_build(myplot)$data[[4]]$r
  p_value <- ggplot2::ggplot_build(myplot)$data[[4]]$p
  return(list("r_value"=r_value, "p_value"=p_value))
}

# list groups/samples of interest for comparison
x_groups <- list(
  B = c("Bulk_B.U", "Bulk_B.S", "Mem_B.U", "Mem_B.S", 
        "Naive_B.U", "Naive_B.S", "Plasmablasts.U"),
  
  CD4T = c("Effector_CD4pos_T.U", "Effector_CD4pos_T.S", "Follicular_T_Helper.U", "Follicular_T_Helper.S",
           "Memory_Teffs.U", "Memory_Teffs.S", "Naive_Teffs.U", "Naive_Teffs.S",
           "Th1_precursors.U", "Th1_precursors.S", "Th17_precursors.U", "Th17_precursors.S",
           "Th2_precursors.U", "Th2_precursors.S", "Memory_Tregs.U", "Memory_Tregs.S", 
           "Naive_Tregs.U", "Naive_Tregs.S", "Regulatory_T.U", "Regulatory_T.S"), 
  
  CD8T = c("CD8pos_T.U", "CD8pos_T.S", "Central_memory_CD8pos_T.U", 
           "Central_memory_CD8pos_T.S", "Effector_memory_CD8pos_T.U",
           "Effector_memory_CD8pos_T.S", "Naive_CD8_T.U", "Naive_CD8_T.S"),
  
  NK = c("Gamma_delta_T.U", "Gamma_delta_T.S", "Mature_NK.U",
    "Mature_NK.S", "Memory_NK.U", "Immature_NK.U"),
  
  Myeloid = c("Monocytes.U", "Monocytes.S", "Myeloid_DCs.U", "pDCs.U")
)


y_groups <- list(
  "B" = c("Ig_4h-unstim_0h", "Ig_18h-unstim_0h"),
  "CD4 T" = c("CD3_4h-unstim_0h", "CD3_18h-unstim_0h"),
  "CD8 T" = c("CD3_4h-unstim_0h", "CD3_18h-unstim_0h"),
  "NK" = c("CD3_4h-unstim_0h", "CD3_18h-unstim_0h"),
  "Mono/Mph" = c("LPS_4h-unstim_0h", "LPS_18h-unstim_0h"))

y_stimuli <- c("Ig", "CD3", "CD3", "CD3", "LPS", "LPS")

# get and merge degs for each y group
# find deg results in each names(y_groups) under the conditions listed in y_groups
y_genes <- lapply(names(y_groups), function(x) get_degs(
  excel.dir = paste0(deResDir, "results_wOrgClusterAnnotation_DEGs/stimuli_time_MAST/"),
  compGroups = y_groups[[x]],
  ident = x))


# y_genes_union represents the unique, non-duplicated list of genes identified as 
# differentially expressed across the two conditions in y_groups for celltype
# essentially list of deg unions per cell type in y groups
y_genes_union <- lapply(y_genes, function(x) {
  tb <- do.call(rbind, x)
  return(unique(tb[[1]]))
})

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

# get avg exp for y genes union
subsetOn <- "expCond.stimuli.time"
timepoint <- c("4h", "18h")
# > seq_along(y_genes_union)
# [1] 1 2 3 4 5 6
y_genes_union_avgExp_stimuli_timepoint1 <- lapply(seq_along(y_genes_union), function(x) {
  n <- names(y_groups)[[x]] #  "B cells"
  l <- names(table(seuratObj_annot[[subsetOn]])) # "CD3_18h"   "CD3_4h"    "Ig_18h"    "Ig_4h"     "LPS_18h"   "LPS_4h"    "unstim_0h"
  l <- l[grep(timepoint[1], l)] # "CD3_4h" "Ig_4h"  "LPS_4h"
  mydata <- getScatterPlot(
    rds = seuratObj_annot, 
    subsetOn = subsetOn, # "expCond.stimuli.time"
    subset = l[grep(y_stimuli[x], l)], # "Ig_4h"
    expCondCheck = "idents", 
    selectedGroups = list(n, n), # ("B cells", "B cells")
    subsetGenes = y_genes_union[[x]], # select degs
    interactive = F
  )
  return(mydata)
})



y_genes_union_avgExp_stimuli_timepoint2 <- lapply(seq_along(y_genes_union), function(x) {
  n <- names(y_groups)[[x]]
  l <- names(table(seuratObj_annot[[subsetOn]]))
  l <- l[grep(timepoint[2], l)]
  mydata <- getScatterPlot(
    rds = seuratObj_annot, 
    subsetOn = subsetOn, 
    subset = l[grep(y_stimuli[x], l)],
    expCondCheck = "idents", 
    selectedGroups = list(n, n),
    subsetGenes = y_genes_union[[x]],
    interactive = F
  )
  return(mydata)
})


y_genes_union_avgExp_unstim <- lapply(seq_along(y_genes_union), function(x) {
  n <- names(y_groups)[[x]]
  mydata <- getScatterPlot(
    rds = seuratObj_annot, 
    subsetOn = subsetOn, 
    subset = "unstim_0h",
    expCondCheck = "idents", 
    selectedGroups = list(n, n),
    subsetGenes = y_genes_union[[x]],
    interactive = F)
  return(mydata)
})



# calculate correlation values
correlation_values_df <- data.frame()
for (x in seq_along(y_genes_union_avgExp_stimuli_timepoint1)) {
  n <- names(y_groups)[[x]] # "B cells"
  for (m in colnames(grp_avg_exp_df)) {
    tryCatch({
      col1 <- data.frame(y_genes_union_avgExp_stimuli_timepoint1[[x]][,1], row.names = rownames(y_genes_union_avgExp_stimuli_timepoint1[[x]])) # gene names row, one column being the average expression 
      col2 <- data.frame(grp_avg_exp_df[[m]], row.names = rownames(grp_avg_exp_df))
      df2plot <- merge(col1, col2, by = "row.names", all = FALSE) # inner join, three columns 
      colnames(df2plot) <- c("rowname","Group1","Group2")
      scatter_values <- plotscatter(df2plot)
      new_row <- data.frame(x = paste0(n, " ", y_stimuli[x], " ", timepoint[1]), y = m, R = scatter_values$r_value, p = scatter_values$p_value)
      correlation_values_df <- rbind(correlation_values_df, new_row)
    }, error = function(e) {
      message("Error occurred in compGroup ", n, " & ", m, ": ", e$message)
    })
    next
  }
}

for (x in seq_along(y_genes_union_avgExp_stimuli_timepoint2)) {
  n <- names(y_groups)[[x]]
  for (m in colnames(grp_avg_exp_df)) {
    tryCatch({
      col1 <- data.frame(y_genes_union_avgExp_stimuli_timepoint2[[x]][,1], row.names = rownames(y_genes_union_avgExp_stimuli_timepoint2[[x]]))
      col2 <- data.frame(grp_avg_exp_df[[m]], row.names = rownames(grp_avg_exp_df))
      df2plot <- merge(col1, col2, by = "row.names", all = FALSE)
      colnames(df2plot) <- c("rowname","Group1","Group2")
      scatter_values <- plotscatter(df2plot)
      new_row <- data.frame(x = paste0(n, " ", y_stimuli[x], " ", timepoint[2]), y = m, R = scatter_values$r_value, p = scatter_values$p_value)
      correlation_values_df <- rbind(correlation_values_df, new_row)
    }, error = function(e) {
      message("Error occurred in compGroup ", n, " & ", m, ": ", e$message)
    })
    next
  }
}
for (x in seq_along(y_genes_union_avgExp_unstim)) {
  n <- names(y_groups)[[x]]
  for (m in colnames(grp_avg_exp_df)) {
    tryCatch({
      col1 <- data.frame(y_genes_union_avgExp_unstim[[x]][,1], row.names = rownames(y_genes_union_avgExp_unstim[[x]]))
      col2 <- data.frame(grp_avg_exp_df[[m]], row.names = rownames(grp_avg_exp_df))
      df2plot <- merge(col1, col2, by = "row.names", all = FALSE)
      colnames(df2plot) <- c("rowname","Group1","Group2")
      scatter_values <- plotscatter(df2plot)
      new_row <- data.frame(x = paste0(n, " unstim_0h"), y = m, R = scatter_values$r_value, p = scatter_values$p_value)
      correlation_values_df <- rbind(correlation_values_df, new_row)
    }, error = function(e) {
      message("Error occurred in compGroup ", n, " & ", m, ": ", e$message)
    })
    next
  }
}

# try row.names = False
write.table(correlation_values_df, paste0(deResDir, blResDir, "/correlation_values_bl_de_union.txt"), sep = "\t", quote = FALSE)




# 3.0 plot corrplot ------------------------------------------------------------
corr_df <- read.table(paste0(deResDir, blResDir, "/correlation_values_bl_de_union.txt"), sep = "\t")
colnames(corr_df)[3] <- "r"
colnames(corr_df)[4] <- "p_value"

corr_df$r <- as.numeric(corr_df$r)
corr_df$p_value <- as.numeric(corr_df$p_value)

unique(corr_df$x)

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

cell_types <- c("B cells", 
                "CD4 T cells", 
                "CD8 T cells", 
                "NK cells",
                "Monocyte-derived Mph",
                "Alveolar Mph MT-positive")
stimuli.time <- c("unstim_0h", 
                  "LPS_4h", "LPS_18h",
                  "CD3_4h", "CD3_18h",
                  "Ig_4h", "Ig_18h")      
# Initialize an empty vector to store the ordered combinations
order_vector1 <- c("B cells unstim_0h", "B cells Ig 4h", "B cells Ig 18h",
                   "CD4 T cells unstim_0h", "CD4 T cells CD3 4h", "CD4 T cells CD3 18h",    
                   "CD8 T cells unstim_0h", "CD8 T cells CD3 4h", "CD8 T cells CD3 18h",
                   "NK cells unstim_0h", "NK CD3 4h", "NK CD3 18h",
                   "Monocyte-derived Mph unstim_0h", "Monocyte-derived Mph LPS 4h", "Monocyte-derived Mph LPS 18h",       
                   "Alveolar Mph MT-positive unstim_0h", "Alveolar Mph MT-positive LPS 4h", "Alveolar Mph MT-positive LPS 18h")

unique(corr_df$y)
 
        
corr_matrix_ordered <- corr_matrix[order_vector, order_vector]

library(corrplot)

pdf("integration/integration_1/integration_1_louvain/cell_clustering/correlation.pdf", width = 10.73, height = 6.59)
# Plotting the matrix
corrplot(corr_matrix, type = "lower", order = "original", 
         tl.col = "black", 
         tl.srt = 45, 
         tl.cex = 0.6,
         # tl.pos="n",
         # col = colorRampPalette(c(rep("#62a1db", 10),
         #                          "white", "#b01111"))(400),
         method = "circle"
         # ,
         # col.lim = c(0.64, 1)
         )

dev.off()


# 4.0 plot heatmap ------------------------------------------------------------
correlation_values_df <- read.table(paste0(deResDir, blResDir, "/correlation_values_bl_de_union.txt"), sep = "\t", header = TRUE)

correlation_values_df <- correlation_values_df[!(correlation_values_df$x %in% "Alveolar Mph MT-positive unstim_0h"), ]
correlation_values_df <- correlation_values_df[!(correlation_values_df$x %in% "Alveolar Mph MT-positive LPS 4h"), ]
correlation_values_df <- correlation_values_df[!(correlation_values_df$x %in% "Alveolar Mph MT-positive LPS 18h"), ]

heatmap_mt <- correlation_values_df[,-4] %>% pivot_wider(names_from = y, values_from = R) %>% column_to_rownames(var = "x") %>% as.matrix()
mat_pval <- correlation_values_df[,-3] %>% pivot_wider(names_from = y, values_from = p) %>% column_to_rownames(var = "x") %>% as.matrix()
max(mat_pval, na.rm = TRUE)



ct_order <- c("B","CD4 T","CD8 T","NK","Mono/Mph")
tm_order <- c("Ig","CD3","CD3","CD3","LPS")

group_order <- c()
for (c in seq_along(ct_order)) {
  group_order <- c(group_order, paste0(ct_order[c], " unstim_0h"))
  group_order <- c(group_order, paste0(ct_order[c], " ", tm_order[c], " ", timepoint[1]))
  group_order <- c(group_order, paste0(ct_order[c], " ", tm_order[c], " ", timepoint[2]))
}
group_order <- group_order[group_order %in% rownames(heatmap_mt)] # "B cells unstim_0h", "B cells Ig 4h" etc

us_order_list <- lapply(x_groups, function (x) {
  us <- gsub(".*\\.(S|U)$", "\\1", x)
  us_order_idx <- c(which(us=="U"),which(us=="S"))
  us_order <- x[us_order_idx]
})
us_order <- unlist(us_order_list, use.names = F) 

# > us_order
# [1] "Bulk_B.U"                   "Mem_B.U"                    "Naive_B.U"                  "Plasmablasts.U"             "Bulk_B.S"                   "Mem_B.S"                   
# [7] "Naive_B.S"                  "Effector_CD4pos_T.U"        "Follicular_T_Helper.U"      "Memory_Teffs.U"             "Naive_Teffs.U"              "Th1_precursors.U"          
# [13] "Th17_precursors.U"          "Th2_precursors.U"           "Memory_Tregs.U"             "Naive_Tregs.U"              "Regulatory_T.U"             "Effector_CD4pos_T.S"       
# [19] "Follicular_T_Helper.S"      "Memory_Teffs.S"             "Naive_Teffs.S"              "Th1_precursors.S"           "Th17_precursors.S"          "Th2_precursors.S"          
# [25] "Memory_Tregs.S"             "Naive_Tregs.S"              "Regulatory_T.S"             "CD8pos_T.U"                 "Central_memory_CD8pos_T.U"  "Effector_memory_CD8pos_T.U"
# [31] "Naive_CD8_T.U"              "CD8pos_T.S"                 "Central_memory_CD8pos_T.S"  "Effector_memory_CD8pos_T.S" "Naive_CD8_T.S"              "Gamma_delta_T.U"           
# [37] "Mature_NK.U"                "Memory_NK.U"                "Immature_NK.U"              "Gamma_delta_T.S"            "Mature_NK.S"                "Monocytes.U"               
# [43] "Myeloid_DCs.U"              "pDCs.U"                     "Monocytes.S"    

# match each column in heatmap_mt with blood cell type
idx_list <- lapply(x_groups, function(x) which(colnames(heatmap_mt) %in% x))
x_split <- rep(NA, ncol(heatmap_mt))
for (n in names(idx_list)) {
  x_split[idx_list[[n]]] <- n
}




cell_types <- c(
                "Monocyte-derived Mph",
                "B cells", 
                "CD4 T cells", 
                "CD8 T cells", 
                "NK cells",
                "EC general capillary", 
                "Migratory DC")



selectedCol <- c("#6d98ba", "#d4c2fc", "#ff758f", "#98c9a3", "#ffb4a2", "#ddb892", "#b2f7ef")

names(selectedCol) <- cell_types
color.order <- c("#d4c2fc", "#ff758f", "#98c9a3", "#ffb4a2", "#6d98ba")

hm <- Heatmap(heatmap_mt, name = "Pearson's r", na_col = "black",
              # col = colorRamp2(c(-0.1, 0, 0.2, 0.4, 0.6, 0.7), c("#60dbe8","white","#ffdd00","#ff9b54", "#EA4749", "#720026")),
              col = colorRampPalette(c("#172896","#60dbe8","white","#ffdd00","#EA4749"))(100),
              # col = viridis(100),
              cluster_rows = F, 
              cluster_columns = F, 
              row_names_side = "left",
              column_names_side = "top", 
              row_order = group_order,
              column_order = us_order,
              row_labels = gsub("^(.*) (.* .*?|unstim_0h)$", "\\2", rownames(heatmap_mt), perl = T),
              column_labels = gsub("\\.[US]$", "", colnames(heatmap_mt), perl = T),
              row_names_gp = gpar(fontsize = 10),
              column_names_gp = gpar(fontsize = 10),
              row_split = factor(gsub("^(.*) (.* .*?|unstim_0h)$", "\\1", rownames(heatmap_mt)), levels = ct_order),
              column_split = factor(x_split, levels = names(x_groups)),
              left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = color.order),
                                                               labels = ct_order,
                                                               labels_gp = gpar(col = "white", fontsize = 10))),
              top_annotation = HeatmapAnnotation(Condition = factor(gsub(".*\\.(S|U)$", "\\1", colnames(heatmap_mt)), levels = c("U","S")),
                                                 col = list(Condition = c("U"="#172896", "S"="#EA4749"))),
              cell_fun = function(j, i, x, y, w, h, fill) {
                if (mat_pval[i, j]>=0.05){
                  grid.rect(x = x, y = y, width = w, height = h, gp = gpar(col = "grey", fill = "white"))
                  # grid.text(sprintf("%.3f", mat_pval[i, j]), x, y, gp = gpar(fontsize = 3)) 
                }
              },
              row_title = NULL,
              column_title = NULL,
              row_gap = unit(2,"mm"),
              column_gap = unit(2,"mm"))

pdf(paste0(deResDir, blResDir, "/correlation_values_bl_de_union.pdf"), width = 10.46, height = 4.46)
draw(hm)
dev.off()


# 5.0 plot heatmap per single-cell and bulk RNA sample group ----------------------------------
g <- list(
  "B cells"="B",
  "CD4 T cells"="CD4T",
  "CD8 T cells"="CD8T",
  "NK cells"="NK",
  "Monocyte-derived Mph"="Myeloid")


global_min <- min(heatmap_mt, na.rm = TRUE)
global_max <- max(heatmap_mt, na.rm = TRUE)
color_scale <- colorRamp2(c(0.2, 0.3, 0.4, 0.5, 0.6), c("#172896","#60dbe8","white","#ffdd00","#EA4749"))

for (c in names(y_groups)[-6]) {
  heatmap_mt.g <- heatmap_mt[grep(c, rownames(heatmap_mt)),colnames(heatmap_mt) %in% x_groups[[g[[c]]]]]
  mat_pval.g <- mat_pval[grep(c, rownames(mat_pval)),colnames(mat_pval) %in% x_groups[[g[[c]]]]]
  group_order.g <- group_order[group_order %in% rownames(heatmap_mt.g)]
  us_order.g <- us_order[us_order %in% colnames(heatmap_mt.g)]
  hm <- Heatmap(heatmap_mt.g, name = "R", na_col = "black",
                # col = colorRamp2(c(min(heatmap_mt.g), max(heatmap_mt.g)), c("white","#EA4749")),
                # col =  colorRampPalette(c("#172896","#60dbe8","white","#ffdd00","#EA4749"))(100), 
                col = color_scale,
                cluster_rows = F, 
                cluster_columns = F, 
                row_names_side = "left",
                column_names_side = "top", 
                row_order = group_order.g,
                column_order = us_order.g,
                row_labels = rownames(heatmap_mt.g),
                column_labels = gsub("\\.[US]$", "", colnames(heatmap_mt.g), perl = T),
                row_names_gp = gpar(fontsize = 10),
                column_names_gp = gpar(fontsize = 10),
                top_annotation = HeatmapAnnotation(Condition = factor(gsub(".*\\.(S|U)$", "\\1", colnames(heatmap_mt.g)), levels = c("U","S")),
                                                   col = list(Condition = c("U"="#172896", "S"="#EA4749"))),
                cell_fun = function(j, i, x, y, w, h, fill) {
                  if (mat_pval.g[i, j]>=0.05){
                    grid.rect(x = x, y = y, width = w, height = h, gp = gpar(col = "grey", fill = "white"))
                    # grid.text(sprintf("%.3f", mat_pval.g[i, j]), x, y, gp = gpar(fontsize = 3)) 
                  }
                },
                row_title = NULL,
                column_title = NULL,
                row_gap = unit(2,"mm"),
                column_gap = unit(2,"mm")
                )
  pdf(paste0(deResDir, blResDir, "/correlation_heatmap__bl_de_union_", sub("\\W", "_", c), ".pdf"), width = 5, height = 3)
  draw(hm)
  dev.off()
}


# 6.0 pull up some r values from 2.0 ####
corr <- read.table(paste0(deResDir, blResDir, "/correlation_values_bl_de_union.txt"), sep = "\t", header=1)

corr %>%
  filter(x == "B cells unstim_0h") %>%
  arrange(desc(R)) 

corr %>%
  filter(x == "B cells unstim_0h" & grepl("B\\.", y)) %>% 
  arrange(desc(R))  

corr %>%
  filter(grepl("B cells Ig", x) & grepl("B\\.", y)) %>% 
  arrange(desc(R))  

corr %>%
  filter(grepl("B cells", x) & (y %in% x_groups[["B"]])) %>% 
  arrange(desc(x), desc(R)) 

corr %>%
  filter(grepl("CD4 T cells", x) & (y %in% x_groups[["CD4T"]])) %>% 
  arrange(desc(x), desc(R))  


corr %>%
  filter(grepl("CD8 T cells", x) & (y %in% x_groups[["CD8T"]])) %>% 
  arrange(desc(x), desc(R))  

corr %>%
  filter(grepl("NK cells", x) & (y %in% x_groups[["NK"]])) %>% 
  arrange(desc(x), desc(R))  



corr %>%
  filter(grepl("Monocyte", x) & (y %in% x_groups[["Myeloid"]])) %>% 
  arrange(desc(x), desc(R))  




# # skip 7.0 PCoA to compare bulk_exp_norm from 1.0 with seuratObj_annot ------------------------
# 
# library(vegan) # This is vegan 2.6-8
# 
# sc_matrix <- as.matrix(GetAssayData(seuratObj_annot, layer = "data"))
# shared_genes <- intersect(rownames(bulk_exp_norm), rownames(sc_matrix))
# 
# bulk_exp_norm_sub <- bulk_exp_norm[shared_genes, ]
# sc_matrix_sub <- sc_matrix[shared_genes, ]
# 
# bulk_sc_exp_matrix <- cbind(bulk_exp_norm_sub, sc_matrix_sub)
# 
# saveRDS(bulk_sc_exp_matrix, file = "integration/integration_1/integration_1_louvain/blood_lung/bulk_sc_exp_matrix.rds")
# 
# # Reopen the RDS file on randi | DO NOT NEED TO 
# # scp /Users/funingtian/Dropbox/UChicago/CRI_ftian1/scRNAseq/integration/integration_1/integration_1_louvain/blood_lung/bulk_sc_exp_matrix.rds ftian1@randi.cri.uchicago.edu:/gpfs/data/schoettler-lab/ftian1/robi/integration_1/integration_1_louvain/blood_lung
# bulk_sc_exp_matrix <- readRDS("integration/integration_1/integration_1_louvain/blood_lung/bulk_sc_exp_matrix.rds")
# # bulk_sc_exp_matrix <- readRDS("integration_1/integration_1_louvain/blood_lung/bulk_sc_exp_matrix.rds")
# 
# 
# bulk_sc_pca <- prcomp(t(bulk_sc_exp_matrix), scale. = TRUE)
# 
# # View PCA results
# summary(bulk_sc_pca)
# 
# saveRDS(bulk_sc_pca, file = "integration/integration_1/integration_1_louvain/blood_lung/bulk_sc_pca.rds")
# 
# 
# # Get PCA scores
# bulk_sc_pca <- readRDS("integration/integration_1/integration_1_louvain/blood_lung/bulk_sc_pca.rds")
#   
# bulk_sc_pca_scores <- bulk_sc_pca$x
# 
# # Create a data frame with PCA results
# # Assuming you have a vector indicating the source (bulk or single-cell)
# source <- c(rep("Blood", ncol(bulk_exp_norm_sub)), rep("Lung", ncol(sc_matrix_sub)))
# pca_df <- data.frame(PC1 = bulk_sc_pca_scores[, 1], PC2 = bulk_sc_pca_scores[, 2], Source = source)
# 
# 
# x_groups <- list(
#   B = c(
#     "Bulk_B.U",
#     "Bulk_B.S",
#     "Mem_B.U",
#     "Mem_B.S",
#     "Naive_B.U",
#     "Naive_B.S",
#     "Plasmablasts.U"
#   ),
#   CD4T = c(
#     "Effector_CD4pos_T.U",
#     "Effector_CD4pos_T.S",
#     "Follicular_T_Helper.U",
#     "Follicular_T_Helper.S",
#     "Memory_Teffs.U",
#     "Memory_Teffs.S",
#     "Naive_Teffs.U",
#     "Naive_Teffs.S",
#     "Th1_precursors.U",
#     "Th1_precursors.S",
#     "Th17_precursors.U",
#     "Th17_precursors.S",
#     "Th2_precursors.U",
#     "Th2_precursors.S",
#     "Memory_Tregs.U",
#     "Memory_Tregs.S",
#     "Naive_Tregs.U",
#     "Naive_Tregs.S",
#     "Regulatory_T.U",
#     "Regulatory_T.S"
#   ),
#   CD8T = c(
#     "CD8pos_T.U",
#     "CD8pos_T.S",
#     "Central_memory_CD8pos_T.U",
#     "Central_memory_CD8pos_T.S",
#     "Effector_memory_CD8pos_T.U",
#     "Effector_memory_CD8pos_T.S",
#     "Naive_CD8_T.U",
#     "Naive_CD8_T.S"
#   ),
#   NK = c(
#     "Gamma_delta_T.U",
#     "Gamma_delta_T.S",
#     "Mature_NK.U",
#     "Mature_NK.S",
#     "Memory_NK.U",
#     "Immature_NK.U"
#   ),
#   Myeloid = c("Monocytes.U", "Monocytes.S", "Myeloid_DCs.U", "pDCs.U")
# )
# 
# # 
# # y_groups <- list(
# #   "B cells" = c("Ig_4h-unstim_0h", "Ig_18h-unstim_0h"),
# #   "CD4 T cells" = c("CD3_4h-unstim_0h", "CD3_18h-unstim_0h"),
# #   "CD8 T cells" = c("CD3_4h-unstim_0h", "CD3_18h-unstim_0h"),
# #   "NK cells" = c("CD3_4h-unstim_0h", "CD3_18h-unstim_0h"),
# #   "Monocyte-derived Mph" = c("LPS_4h-unstim_0h", "LPS_18h-unstim_0h"),
# #   "Alveolar Mph MT-positive" = c("LPS_4h-unstim_0h", "LPS_18h-unstim_0h"))
# # 
# # y_stimuli <- c("Ig", "CD3", "CD3", "CD3", "LPS", "LPS")
# # 
# 
# for (i in 1:166) {
#   row_name <- strsplit(rownames(pca_df)[i], "\\.")[[1]]
#   
#   if (length(row_name) > 1) {
#     # Get the part after the first period
#     name <- paste(row_name[2:length(row_name)], collapse = ".") 
#     for (group_name in names(x_groups)) {
#       if (name %in% x_groups[[group_name]]) {
#         pca_df$expCond.celltype.stimuli.time[i] <- paste(group_name, name) # Combine group name and row name
#         pca_df$expCond.celltype[i] <- group_name
#         break # Exit the loop once a match is found
#       }
#     }
#   }
# }
# 
# 
# pca_df$expCond.celltype.stimuli.time[167:nrow(pca_df)] <- seuratObj_annot@meta.data[rownames(pca_df)[167:nrow(pca_df)], "expCond.celltype.stimuli.time"]
# pca_df$expCond.celltype[167:nrow(pca_df)] <- as.character(seuratObj_annot@meta.data[rownames(pca_df)[167:nrow(pca_df)], "cluster_annotation"])
# 
# 
# pca_df <- pca_df %>%
#   mutate(expCond.celltype.stimuli.time = str_replace_all(expCond.celltype.stimuli.time, 
#                                                    c("Alveolar Mph MT-positive" = "Alv Mph",
#                                                      "EC general capillary" = "EC",
#                                                      "Monocyte-derived Mph" = "Mono Mph",
#                                                      "B cells" = "B",
#                                                      "CD4 T cells" = "CD4 T",
#                                                      "CD8 T cells" = "CD8 T",
#                                                      "NK cells" = "NK",
#                                                      "Migratory DC" = "DC")),
#          expCond.celltype = str_replace_all(expCond.celltype, 
#                                                          c("Alveolar Mph MT-positive" = "Alv Mph",
#                                                            "EC general capillary" = "EC",
#                                                            "Monocyte-derived Mph" = "Mono Mph",
#                                                            "B cells" = "B",
#                                                            "CD4 T cells" = "CD4 T",
#                                                            "CD8 T cells" = "CD8 T",
#                                                            "CD4T" = "CD4 T",
#                                                            "CD8T" = "CD8 T",
#                                                            "NK cells" = "NK",
#                                                            "Migratory DC" = "DC")))
# 
# # Visualize PCA results
# pca_df %>%
#   filter(!expCond.celltype %in% "EC",
#          !expCond.celltype %in% "DC",
#          !expCond.celltype %in% "Alv Mph") %>%
#   ggplot(aes(x = PC1, y = PC2, color = expCond.celltype, shape = Source)) +
#   geom_point(size=1.8) +
#   labs(title = "PCA of Bulk & scRNA-seq Data", 
#        x = "Principal Component 1", 
#        y = "Principal Component 2", 
#        color = "Cell Type",
#        shape = "Bulk/scRNA-seq") +
#   scale_color_manual(values = c("B" = "#d4c2fc", 
#                                 "CD4 T" = "#ff758f",
#                                 "CD8 T" = "#98c9a3",
#                                 "NK" = "#ffb4a2",
#                                 "Mono Mph" = "#6d98ba",
#                                 "Myeloid" = "#6d98ba")) +  # Customize colors
#   theme_bw()
# 
# seuratObj_annot@meta.data %>%
#   filter(percent.mt >=10)%>%
#   dim()
# 
# seuratObj_annot@meta.data %>%
#   dim()
# 
# # 8.0 