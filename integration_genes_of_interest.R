rm(list = ls())


library(Seurat)
library(readxl)
library(writexl)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(gtable)
library(grid)
library(gridExtra)
library(stringr)   
library(patchwork)
library(VennDiagram)


input.dir <- "integration/integration_1/integration_1_louvain/genes_of_interest/"
out.dir <- "integration_2/integration_2_leiden/genes_of_interest/"

# Read the Excel file starting from line 3
genes <- read_excel("supplementary_tables_20240708[38].xlsx", sheet = "S8", skip = 3)
genes <- genes[c(1,2,4,5)]
colnames(genes) <- c("COA_genes", "COA_scores", "AOA_genes", "AOA_scores")

genes_0903 <- read.table(paste0(input.dir, "Dot_plot_cell_type_gene_list_09032024.txt"), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(genes_0903) <- c("genes")

rds <- readRDS('integration_2/integration_2_leiden/RDS_Dir/integration_2_leiden_annot.rds')

DefaultAssay(rds) <- "RNA"

#### skip dotplot genes_0903 celltype treatment ####
rds@meta.data$cluster_annotation2_expCond_stimuli_time <- paste0(
  rds@meta.data$cluster_annotation2, 
  "|", 
  rds@meta.data$expCond.stimuli.time
)


# Extract data for the specified features
missing_genes <- setdiff(unique(genes_0903$genes), rownames(GetAssayData(object = rds, slot = "data")))

# Display missing genes, if any
if(length(missing_genes) > 0) {
  print("These genes are missing:")
  print(missing_genes)
} else {
  print("All genes are present.")
}

# all genes are present

existing_genes <- intersect(unique(genes_0903$genes), rownames(GetAssayData(object = rds, slot = "data")))

dotplot <- DotPlot(rds,
        features = existing_genes,
                   cols = c('#D3D3D3', '#CC0000'),
                   scale = T, scale.by = 'size',
                   dot.min = 0) + RotatedAxis()


genes_0903_exp <- as.data.frame(t(GetAssayData(object = rds, layer = "data")[existing_genes, ]))

# Calculate the average expression per cluster
genes_0903_avg_exp <- genes_0903_exp %>%
  rownames_to_column(var = "cell") %>%
  inner_join(rds@meta.data %>% rownames_to_column(var = "cell"), by = "cell") %>%
  group_by(cluster_annotation_expCond_stimuli_time) %>%
  summarise(across(2:(ncol(genes_0903_exp)+1), mean, na.rm = TRUE))

# Calculate the percentage of cells expressing each gene per cluster
genes_0903_pct_exp <- genes_0903_exp %>%
  rownames_to_column(var = "cell") %>%
  inner_join(rds@meta.data %>% rownames_to_column(var = "cell"), by = "cell") %>%
  group_by(cluster_annotation_expCond_stimuli_time) %>%
  summarise(across(2:(ncol(genes_0903_exp)+1), ~mean(. > 0) * 100))

# Combine average expression and percentage expression data
genes_0903_avg_pct_exp <- genes_0903_avg_exp %>%
  pivot_longer(-cluster_annotation_expCond_stimuli_time, names_to = "gene", values_to = "expression") %>%
  inner_join(genes_0903_pct_exp %>% pivot_longer(-cluster_annotation_expCond_stimuli_time, names_to = "gene", values_to = "pct.exp"), 
             by = c("cluster_annotation_expCond_stimuli_time", "gene"))

# Scale the expression data if needed
genes_0903_avg_pct_exp <- genes_0903_avg_pct_exp %>%
  group_by(gene) %>%
  mutate(expression = scale(expression))

genes_0903_avg_pct_exp <- genes_0903_avg_pct_exp %>%
  separate(
    cluster_annotation_expCond_stimuli_time, 
    into = c("cluster_annotation", "expCond_stimuli_time"), 
    sep = "\\|"
  )

genes_0903_avg_pct_exp <- genes_0903_avg_pct_exp %>%
  mutate(cluster_annotation_plot = case_when(
    cluster_annotation == "Alveolar Mph MT-positive" ~ "Alveolar Mph\nMT-positive",
    cluster_annotation == "EC general capillary" ~ "EC general\ncapillary",
    cluster_annotation == "Monocyte-derived Mph" ~ "Monocyte-\nderived Mph",
    TRUE ~ cluster_annotation  # Default case for other annotations
  ))

genes_0903_avg_pct_exp <- genes_0903_avg_pct_exp %>%
  mutate(pct.exp.shape = case_when(
    pct.exp < 25 ~ "<25%",
    pct.exp >= 25 & pct.exp < 50 ~ "25-50%",
    pct.exp >= 50 & pct.exp < 75 ~ "50-75%",
    pct.exp >= 75 ~ "75-100%"
  ))


ggplot(genes_0903_avg_pct_exp, aes(x = factor(expCond_stimuli_time, 
                                              levels=c("unstim_0h", "LPS_4h", "LPS_18h", "CD3_4h", "CD3_18h","Ig_4h", "Ig_18h")), 
                                   y = gene,
                                   shape = pct.exp.shape, 
                                   fill = expression)) +
  scale_shape_manual(values =  c("<25%" = 21, "25-50%" = 22, "50-75%" = 23, "75-100%" = 24)) +
  geom_point(size=2.7) +
  scale_shape_manual(values =  c("<25%" = 21, "25-50%" = 22, "50-75%" = 23, "75-100%" = 24)) +
  # scale_color_manual(values =  c("<0.01" = "#0000ff", "Not significant" = "white", "NA" = "white")) +
  # scale_discrete_manual(aesthetics = "stroke", 
  #                       values =  c("<0.01" = 1, "Not significant" = 0, "NA" = 0), 
  #                       guide = "none") +
  scale_fill_gradientn(colors = colorRampPalette(c('#ffffff', '#e0fbfc', "#98c1d9", "#3d5a80"))(100)) +
  facet_grid(. ~ cluster_annotation_plot, scales="free", space="free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "left",
        strip.background = element_rect(fill = "white", color = "black", size = 1),
        strip.text.y = element_blank(),
        strip.background.y = element_blank()) +
  labs(x = "Celltype, stimuli and time",
       y = "Genes",
       shape = "Percent expressed",
       fill = "Scaled average expression")












pdf(paste0(out.dir,"/feature_plot_AOA_genes.pdf"), width = 12, height = 7)
print(FeaturePlot_AOA)
dev.off()



#### dotplot genes_0903 celltype  ####


rds@meta.data$cluster_annotation3_expCond_stimuli_time <- paste0(
  rds@meta.data$cluster_annotation3, 
  "|", 
  rds@meta.data$expCond.stimuli.time
)


# Extract data for the specified features
missing_genes <- setdiff(unique(genes_0903$genes), rownames(GetAssayData(object = rds, layer = "data")))

# Display missing genes, if any
if(length(missing_genes) > 0) {
  print("These genes are missing:")
  print(missing_genes)
} else {
  print("All genes are present.")
}

# all genes are present

existing_genes <- intersect(unique(genes_0903$genes), rownames(GetAssayData(object = rds, slot = "data")))
# existing_genes <- unique(genes_0903$genes)
dotplot <- DotPlot(rds,
                   features = existing_genes,
                   cols = c('#D3D3D3', '#CC0000'),
                   scale = T, scale.by = 'size',
                   dot.min = 0) + RotatedAxis()


genes_0903_exp <- as.data.frame(t(GetAssayData(object = rds, layer = "data")[existing_genes, ]))

# Calculate the average expression per cluster
genes_0903_avg_exp <- genes_0903_exp %>%
  rownames_to_column(var = "cell") %>%
  inner_join(rds@meta.data %>% rownames_to_column(var = "cell"), by = "cell") %>%
  group_by(cluster_annotation3) %>%
  summarise(across(2:(ncol(genes_0903_exp)+1), mean, na.rm = TRUE))

# Calculate the percentage of cells expressing each gene per cluster
genes_0903_pct_exp <- genes_0903_exp %>%
  rownames_to_column(var = "cell") %>%
  inner_join(rds@meta.data %>% rownames_to_column(var = "cell"), by = "cell") %>%
  group_by(cluster_annotation3) %>%
  summarise(across(2:(ncol(genes_0903_exp)+1), ~mean(. > 0) * 100))

# Combine average expression and percentage expression data
genes_0903_avg_pct_exp <- genes_0903_avg_exp %>%
  pivot_longer(-cluster_annotation3, names_to = "gene", values_to = "expression") %>%
  inner_join(genes_0903_pct_exp %>% pivot_longer(-cluster_annotation3, names_to = "gene", values_to = "pct.exp"), 
             by = c("cluster_annotation3", "gene")) %>%
  filter(cluster_annotation3 %in% c( "B", "CD4 T", 
                                     "CD8 T", "NK",
                                     "Mono/Mph", "EC", "DC"))

# Scale the expression data if needed
genes_0903_avg_pct_exp <- genes_0903_avg_pct_exp %>%
  group_by(gene) %>%
  mutate(expression = scale(expression))


genes_0903_avg_pct_exp <- genes_0903_avg_pct_exp %>%
  mutate(pct.exp.shape = case_when(
    pct.exp < 25 ~ "<25%",
    pct.exp >= 25 & pct.exp < 50 ~ "25-50%",
    pct.exp >= 50 & pct.exp < 75 ~ "50-75%",
    pct.exp >= 75 ~ "75-100%"
  ))

genes_0903_avg_pct_exp_filtered <- genes_0903_avg_pct_exp %>%
  filter(pct.exp > 0)


pdf(paste0(out.dir,"/genes_0903_dotplot.pdf"), width = 12.47, height = 2.36)
ggplot(genes_0903_avg_pct_exp, aes(y = factor(cluster_annotation3,levels=c("B", 
                                                                          "CD4 T", 
                                                                          "CD8 T", 
                                                                          "NK",
                                                                          "Mono/Mph",
                                                                          "EC", 
                                                                          "DC")), 
                                   x = factor(gene, levels=unique(gene)),
                                   size = pct.exp, 
                                   fill = expression)) +
  # scale_shape_manual(values =  c("<25%" = 21, "25-50%" = 22, "50-75%" = 23, "75-100%" = 24)) +
  geom_point(shape=21) +
  # scale_color_manual(values =  c("<0.01" = "#0000ff", "Not significant" = "white", "NA" = "white")) +
  # scale_discrete_manual(aesthetics = "stroke", 
  #                       values =  c("<0.01" = 1, "Not significant" = 0, "NA" = 0), 
  #                       guide = "none") +
  scale_fill_gradientn(colors = colorRampPalette(c('#f9dbbd', '#ffa5ab', "#da627d", "#a53860", "#450920"))(100)) +
  # facet_grid(. ~ cluster_annotation, scales="free", space="free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "right",
        legend.box = "horizontal",
        legend.title = element_text(hjust = 0.5, vjust = 0, size = 10),  # Center-align legend titles
        strip.background = element_rect(fill = "white", color = "black", size = 1),
        strip.text.y = element_blank(),
        strip.background.y = element_blank()) +
  labs(y = "Celltype",
       x = "Genes",
       size = "Percent\nexpressed",
       fill = "Scaled\naverage\nexpression") 

dev.off()




#############################################
#############################################
#### feature plots umap for genes of interest ####
#############################################
#############################################
##-- 1.0 HLA-DQA2 HLA_DQB2 ####
names(table(Idents(rds)))

rds_sub <- subset(rds, subset = expCond.stimuli.time %in% c('unstim_0h', 'Ig_4h', 'Ig_18h'))

rds_sub$expCond.stimuli.time <- factor(rds_sub$expCond.stimuli.time,
                                       levels = c('unstim_0h', 'Ig_4h', 'Ig_18h'))

levels(rds_sub$expCond.stimuli.time)[levels(rds_sub$expCond.stimuli.time) == "unstim_0h"] <- "Unstim_0h"

# FeaturePlot(
#   rds,
#   features = c('HLA-DQA2'),
#   cols = c("lightgray", "red"))
# 
# FeaturePlot(
#   rds,
#   features = c('HLA-DQB2'),
#   cols = c("lightgray", "red"))

p1 <- FeaturePlot(rds_sub, features=c('HLA-DQA2'), 
            split.by = "expCond.stimuli.time", 
            combine=TRUE,
            keep.scale = 'all',
            order = TRUE,
            cols = c("#f4f1de", "#3c096c")            
            ) & theme(legend.position = "right") &
  theme_bw() &
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.y.right = element_blank(),  # Remove ticks on the right side
    axis.text.y.right = element_blank(),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    axis.text = element_text(size = 10),   # Adjust size as needed
    axis.title = element_text(size = 12)  # Adjust size as needed
  )&
  labs(x = "UMAP 1", y = "UMAP 2")


p2 <- FeaturePlot(rds_sub, features=c('HLA-DQB2'), 
            split.by = "expCond.stimuli.time", 
            combine=TRUE,
            keep.scale = 'all',
            order = TRUE,
            cols = c("#f4f1de", "#3c096c"),
) & theme(legend.position = "right") &
  theme_bw() &
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.y.right = element_blank(),  # Remove ticks on the right side
    axis.text.y.right = element_blank(),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    axis.text = element_text(size = 10),   # Adjust size as needed
    axis.title = element_text(size = 12)  # Adjust size as needed
  )&
  labs(x = "UMAP 1", y = "UMAP 2")


pdf(paste0(out.dir,"/feature_plot_HLA-DQAB2.pdf"), width = 9.14, height = 4.7)

p1/p2

dev.off()


##-- 2.0 selected asthma genes ####
names(table(Idents(rds)))

FeaturePlot_AOA <- FeaturePlot(rds, features = genes$AOA_genes[1:11])
FeaturePlot_COA1 <- FeaturePlot(rds, features = genes$COA_genes[1:12])
FeaturePlot_COA2 <- FeaturePlot(rds, features = genes$COA_genes[13:24])
FeaturePlot_COA3 <- FeaturePlot(rds, features = genes$COA_genes[25:36])

pdf(paste0(out.dir,"/feature_plot_AOA_genes.pdf"), width = 12, height = 7)
print(FeaturePlot_AOA)
dev.off()

pdf(paste0(out.dir,"/feature_plot_COA1_genes.pdf"), width = 12, height = 7)
print(FeaturePlot_COA1)
dev.off()

pdf(paste0(out.dir,"/feature_plot_COA2_genes.pdf"), width = 12, height = 7)
print(FeaturePlot_COA2)
dev.off()

pdf(paste0(out.dir,"/feature_plot_COA3_genes.pdf"), width = 12, height = 7)
print(FeaturePlot_COA3)
dev.off()



##-- 3.0 th2 gene markers ####
names(table(Idents(rds)))

th2_genes <- c(
  "IL4", "IL5", "IL13", 
  "CCR3", "CCR4", "CCR6", 
  "PTGDR2", "TCF7", "LEF1", 
  "SELL", "KLF2", "CD109", 
  "GPR15", "CD27", "SLC2A3")

missing_genes <- setdiff(th2_genes, rownames(GetAssayData(object = rds, layer = "data")))

# Display missing genes, if any
if(length(missing_genes) > 0) {
  print("These genes are missing:")
  print(missing_genes)
} else {
  print("All genes are present.")
}

FeaturePlot_th2 <- FeaturePlot(rds, features = th2_genes)

pdf(paste0(out.dir,"/feature_plot_th2_gene_markers.pdf"), width = 19, height = 17)
print(FeaturePlot_th2)
dev.off()









##-- 3.0 CISH CD3 CD4 CD28 with a focus on CD4 T ####
names(table(Idents(rds)))

cish_genes <- c(
  "CISH", "CD3", "CD4", "CD28"
  )

missing_genes <- setdiff(cish_genes, rownames(GetAssayData(object = rds, layer = "data")))

# Display missing genes, if any
if(length(missing_genes) > 0) {
  print("These genes are missing:")
  print(missing_genes)
} else {
  print("All genes are present.")
}

# [1] "These genes are missing:"
# [1] "CD3"

cish_genes <- setdiff(cish_genes, missing_genes)

FeaturePlot_cish <- FeaturePlot(rds, features = cish_genes)

pdf(paste0(out.dir,"/feature_plot_cish_genes_markers.pdf"), width = 19, height = 4)
print(FeaturePlot_cish)
dev.off()










#### top risk genes dot plot genes of interest vs cell types vs treatments ####
# DotPlot <- DotPlot(rds, 
#         features = genes$AOA_genes[1:11], 
#         cols = c('#D3D3D3', '#CC0000'), 
#         scale = T, scale.by = 'size', 
#         dot.min = 0.01) + RotatedAxis()



gene_scores <- genes %>%
  pivot_longer(cols = c(COA_genes, AOA_genes), 
               names_to = "gene_type", 
               values_to = "gene") %>%
  pivot_longer(cols = c(COA_scores, AOA_scores), 
               names_to = "score_type", 
               values_to = "score") %>%
  mutate(gene_type = str_remove(gene_type, "_genes"),
         score_type = str_remove(score_type, "_scores")) %>%
  filter(gene_type == score_type) %>%
  drop_na()

gene_scores <- gene_scores[, !names(gene_scores) %in% "score_type"]

# Create a category based on COA and/or AOA
genes_in_coa <- gene_scores %>% filter(gene_type == "COA") %>% pull(gene)
genes_in_aoa <- gene_scores %>% filter(gene_type == "AOA") %>% pull(gene)

gene_scores <- gene_scores %>%
  mutate(risk = case_when(
    gene %in% genes_in_coa & gene %in% genes_in_aoa ~ "COA_AOA",
    gene %in% genes_in_coa ~ "COA",
    gene %in% genes_in_aoa ~ "AOA",
    TRUE ~ "Other" # For genes not in either list
  ))

gene_scores <- gene_scores %>%
  mutate(position = ifelse(risk == "COA_AOA", "dodge", "stack"))

Seurat::DefaultAssay(rds) <- "RNA"
rds@meta.data$cluster_annotation3_expCond_stimuli_time <- paste0(
  rds@meta.data$cluster_annotation3, 
  "|", 
  rds@meta.data$expCond.stimuli.time
)

# Extract data for the specified features
COA_AOA_exp <- as.data.frame(t(GetAssayData(object = rds, layer = "data")[unique(c(genes$COA_genes, genes$AOA_genes[1:11])), ]))


# Calculate some stats
COA_AOA_exp_meta <- COA_AOA_exp %>%
  rownames_to_column(var = "cell") %>%
  inner_join(rds@meta.data %>% rownames_to_column(var = "cell"), by = "cell")

any(is.na(COA_AOA_exp_meta)) # any nas? 

COA_AOA_exp_meta <- COA_AOA_exp_meta %>%
  mutate(cluster_annotation_plot = case_when(
    # cluster_annotation == "Alveolar Mph MT-positive" ~ "Alveolar Mph\nMT-positive",
    # cluster_annotation == "EC general capillary" ~ "EC general\ncapillary",
    # cluster_annotation == "Monocyte-derived Mph" ~ "Monocyte-\nderived Mph",
    TRUE ~ cluster_annotation3  # Default case for other annotations
  ))

p_values <- data.frame(gene = character(),
                       expCond.stimuli.time = character(),
                       cell_type = character(),
                       p_value = numeric(),
                       stringsAsFactors = FALSE)

cell_types <-  c( "Mono/Mph", "B", "CD4 T", "CD8 T", "NK", "EC", "DC")
baseline <- "unstim_0h"

for (gene in colnames(COA_AOA_exp)) {
  for (cell_type in cell_types) {
    # print(cell_type)
    df1 <- COA_AOA_exp_meta %>% filter(cluster_annotation_plot == cell_type)
    if (baseline %in% df1$expCond.stimuli.time) {
      for (treatment in unique(df1$expCond.stimuli.time)[2:7]) { # exclude the baseline
        df2 <- df1 %>% filter(expCond.stimuli.time %in% c(baseline, treatment))
        # print(sprintf("Gene %s | Treatment comparision %s %s", gene, baseline, treatment))
        res <- wilcox.test(
          df2[df2$expCond.stimuli.time == baseline, gene],
          df2[df2$expCond.stimuli.time == treatment, gene],
          exact = FALSE
        )
        # Store the p-value if it is not NA
        if (is.na(res$p.value)) {
          p_values <- rbind(p_values, data.frame(gene = gene, 
                                                 expCond.stimuli.time = treatment, 
                                                 cell_type = cell_type,
                                                 p_value = "NA"))
        }else{
          p_values <- rbind(p_values, data.frame(gene = gene, 
                                                 expCond.stimuli.time = treatment, 
                                                 cell_type = cell_type,
                                                 p_value = res$p.value))
        }
      }
    }
  }
}


p_values <- p_values %>%
  rename(expCond_stimuli_time = expCond.stimuli.time)

p_values <- p_values %>%
  rename(cluster_annotation_plot = cell_type)


p_values <- p_values %>%
  mutate(p_value = as.numeric(p_value),   # Ensure p_value column is numeric
         p_value_sig = case_when(
           p_value < 0.01 ~ "<0.01",      # Mark significant p-values
           TRUE ~ "Not significant"       # Mark non-significant p-values
         ))



# Calculate the average expression per cluster
COA_AOA_avg_exp <- COA_AOA_exp %>%
  rownames_to_column(var = "cell") %>%
  inner_join(rds@meta.data %>% 
               rownames_to_column(var = "cell") %>%
               filter(cluster_annotation3 %in% cell_types), 
             by = "cell") %>%
  group_by(cluster_annotation3_expCond_stimuli_time) %>%
  summarise(across(2:42, mean, na.rm = TRUE))

# Calculate the percentage of cells expressing each gene per cluster
COA_AOA_pct_exp <- COA_AOA_exp %>%
  rownames_to_column(var = "cell") %>%
  inner_join(rds@meta.data %>% 
               rownames_to_column(var = "cell") %>%
               filter(cluster_annotation3 %in% cell_types), 
             by = "cell") %>%
  group_by(cluster_annotation3_expCond_stimuli_time) %>%
  summarise(across(2:42, ~mean(. > 0) * 100))

# Combine average expression and percentage expression data
COA_AOA_avg_pct_exp <- COA_AOA_avg_exp %>%
  pivot_longer(-cluster_annotation3_expCond_stimuli_time, names_to = "gene", values_to = "expression") %>%
  inner_join(COA_AOA_pct_exp %>% pivot_longer(-cluster_annotation3_expCond_stimuli_time, names_to = "gene", values_to = "pct.exp"), 
             by = c("cluster_annotation3_expCond_stimuli_time", "gene"))

# Scale the expression data if needed
COA_AOA_avg_pct_exp <- COA_AOA_avg_pct_exp %>%
  group_by(gene) %>%
  mutate(expression = scale(expression))

COA_AOA_avg_pct_exp <- COA_AOA_avg_pct_exp %>%
  mutate(risk = case_when(
    gene %in% genes$AOA_genes[1:11] & gene %in% genes$COA_genes ~ "COA_AOA", # If in both lists
    gene %in% genes$AOA_genes[1:11] & !(gene %in% genes$COA_genes) ~ "AOA", # If only in AOA list
    gene %in% genes$COA_genes & !(gene %in% genes$AOA_genes[1:11]) ~ "COA", # If only in COA list
    TRUE ~ "Other"  # For genes not in either list
  ))

COA_AOA_avg_pct_exp <- COA_AOA_avg_pct_exp %>%
  separate(
    cluster_annotation3_expCond_stimuli_time, 
    into = c("cluster_annotation3", "expCond_stimuli_time"), 
    sep = "\\|"
  )

COA_AOA_avg_pct_exp <- COA_AOA_avg_pct_exp %>%
  mutate(pct.exp.shape = case_when(
    pct.exp < 25 ~ "<25%",
    pct.exp >= 25 & pct.exp < 50 ~ "25-50%",
    pct.exp >= 50 & pct.exp < 75 ~ "50-75%",
    pct.exp >= 75 ~ "75-100%"
  ))

COA_AOA_avg_pct_exp <- COA_AOA_avg_pct_exp %>%
  mutate(cluster_annotation_plot = case_when(
    # cluster_annotation == "Alveolar Mph MT-positive" ~ "Alveolar Mph\nMT-positive",
    # cluster_annotation == "EC general capillary" ~ "EC general\ncapillary",
    # cluster_annotation == "Monocyte-derived Mph" ~ "Monocyte-\nderived Mph",
    TRUE ~ cluster_annotation3  # Default case for other annotations
  ))


# add p alues
COA_AOA_avg_pct_exp_p_values <- COA_AOA_avg_pct_exp %>%
  left_join(p_values, by = c("cluster_annotation_plot", 
                              "expCond_stimuli_time", 
                              "gene"))

COA_AOA_avg_pct_exp_p_values <- COA_AOA_avg_pct_exp_p_values %>%
  mutate(p_value_sig = ifelse(is.na(p_value_sig), "NA", p_value_sig))

# COA_AOA_avg_pct_exp_p_values <- COA_AOA_avg_pct_exp_p_values %>%
#   mutate(cluster_annotation_plot = str_replace_all(cluster_annotation_plot, 
#                                                    c("Alveolar Mph\nMT-positive" = "Alv Mph",
#                                                      "EC general\ncapillary" = "EC",
#                                                      "Monocyte-\nderived Mph" = "Mono Mph",
#                                                      "B cells" = "B",
#                                                      "CD4 T cells" = "CD4 T",
#                                                      "CD8 T cells" = "CD8 T",
#                                                      "NK cells" = "NK",
#                                                      "Migratory DC" = "DC")))

any(is.na(COA_AOA_avg_pct_exp_p_values)) # would be yes this time due to unstim vvs unsim

# Create the custom dot plot
selectedCol <- DiscretePalette(n = length(unique(COA_AOA_avg_pct_exp$cluster_annotation3)), palette = 'alphabet')

pdf(paste0(out.dir,"/dot_plot_COA_AOA_celltype_treament_pvalues.pdf"), width = 13, height = 7)

p1 <- ggplot(COA_AOA_avg_pct_exp_p_values, aes(x = factor(expCond_stimuli_time, levels=c("unstim_0h", 
                                                                          "LPS_4h", "LPS_18h",
                                                                          "CD3_4h", "CD3_18h",
                                                                          "Ig_4h", "Ig_18h")), 
                                y = factor(gene, levels=rev(c(unique(c(genes$COA_genes, genes$AOA_genes[1:11]))))),
                                shape = pct.exp.shape, 
                                fill = expression, 
                                color = factor(p_value_sig, levels=c("<0.01", "Not significant", "NA")),
                                stroke = p_value_sig)) +
  geom_point(size=2.7) +
  scale_shape_manual(values =  c("<25%" = 21, "25-50%" = 22, "50-75%" = 23, "75-100%" = 24)) +
  scale_color_manual(values =  c("<0.01" = "#0000ff", "Not significant" = "white", "NA" = "white")) +
  scale_discrete_manual(aesthetics = "stroke", 
                        values =  c("<0.01" = 1, "Not significant" = 0, "NA" = 0), 
                        guide = "none") +
  scale_fill_gradientn(colors = colorRampPalette(c('white', '#eae2b7', '#ffba08',
                                                   "#dc2f02","#9d0208",
                                                   "#370617", "#03071e"))(100),limits = c(-3.0, 7.5)) +
  facet_grid(risk ~ 
               factor(cluster_annotation_plot, 
                      levels = cell_types), 
             scales="free", space="free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "left",
        strip.background = element_rect(fill = "white", color = "black", linewidth = 1),
        strip.text.y = element_blank(),
        strip.background.y = element_blank()) +
  labs(x = "Celltype, stimuli and time", 
       y = "Genes",
       shape = "Percent expressed",
       fill = "Average expression",
       color = expression(italic("p") * " value"))


p2 <- ggplot(gene_scores, aes(x = score, 
           y = factor(gene, levels=rev(c(unique(c(genes$COA_genes, genes$AOA_genes[1:11]))))), 
           color = factor(gene_type, levels=c("COA", "AOA")),
           group = factor(gene_type, levels=c("COA", "AOA")))) +
  # geom_line(linewidth=0.5) +
  geom_point(size=2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 1),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  scale_color_manual(values = c("#faa4bd", "#533b4d")) +
  scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, by = 1)) +
  labs(x = "Risk score", y = "Genes", Color = "Category", color = "Risk type") +
  facet_grid(risk ~ ., scales = "free", space = "free")
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

p1 + p2 + plot_layout(ncol = 2, widths = c(10, 1)) 


dev.off()


gene_scores %>%
  filter(gene == "BACH2")

genes %>%
  filter(COA_genes == "BACH2")



#### ALL risk genes dot plot genes of interest vs cell types vs treatments ####
# Read the Excel file starting from line 3
genes <- read_excel("AdditionalFile1_SupplementaryTables_GM.xlsx", sheet = "S8", skip = 3)
genes <- genes[c(1,2,4,5)]
# different col names
colnames(genes) <- c("AOA_genes", "AOA_scores", "COA_genes", "COA_scores")

gene_scores <- genes %>%
  pivot_longer(cols = c(COA_genes, AOA_genes), 
               names_to = "gene_type", 
               values_to = "gene") %>%
  pivot_longer(cols = c(COA_scores, AOA_scores), 
               names_to = "score_type", 
               values_to = "score") %>%
  mutate(gene_type = str_remove(gene_type, "_genes"),
         score_type = str_remove(score_type, "_scores")) %>%
  filter(gene_type == score_type) %>%
  drop_na()

gene_scores <- gene_scores[, !names(gene_scores) %in% "score_type"]

# Create a category based on COA and/or AOA
genes_in_coa <- gene_scores %>% filter(gene_type == "COA") %>% pull(gene)
genes_in_aoa <- gene_scores %>% filter(gene_type == "AOA") %>% pull(gene)

gene_scores <- gene_scores %>%
  mutate(risk = case_when(
    gene %in% genes_in_coa & gene %in% genes_in_aoa ~ "COA_AOA",
    gene %in% genes_in_coa ~ "COA",
    gene %in% genes_in_aoa ~ "AOA",
    TRUE ~ "Other" # For genes not in either list
  ))

gene_scores <- gene_scores %>%
  mutate(position = ifelse(risk == "COA_AOA", "dodge", "stack"))

Seurat::DefaultAssay(rds) <- "RNA"
rds@meta.data$cluster_annotation3_expCond_stimuli_time <- paste0(
  rds@meta.data$cluster_annotation3, 
  "|", 
  rds@meta.data$expCond.stimuli.time
)

# Extract data for the specified features
missing_genes <- setdiff(unique(gene_scores$gene), rownames(GetAssayData(object = rds, layer = "data")))

# Display missing genes, if any
if(length(missing_genes) > 0) {
  print("These genes are missing:")
  print(missing_genes)
} else {
  print("All genes are present.")
}

# [1] "SOWAHA"   "TPD52L3"  "AQP2"     "KRTAP3-3" "KRTAP4-1" "FLJ20373" "BLTP1"    "OPTC"     "GPR182"   "PRM2"    
# [11] "AQP5"

existing_genes <- intersect(unique(gene_scores$gene), rownames(GetAssayData(object = rds, layer = "data")))

COA_AOA_exp <- as.data.frame(t(GetAssayData(object = rds, layer = "data")[existing_genes, ]))


# Calculate some stats
COA_AOA_exp_meta <- COA_AOA_exp %>%
  rownames_to_column(var = "cell") %>%
  inner_join(rds@meta.data %>% rownames_to_column(var = "cell"), by = "cell")

any(is.na(COA_AOA_exp_meta)) # any nas? 

COA_AOA_exp_meta <- COA_AOA_exp_meta %>%
  mutate(cluster_annotation_plot = case_when(
    # cluster_annotation == "Alveolar Mph MT-positive" ~ "Alveolar Mph\nMT-positive",
    # cluster_annotation == "EC general capillary" ~ "EC general\ncapillary",
    # cluster_annotation == "Monocyte-derived Mph" ~ "Monocyte-\nderived Mph",
    TRUE ~ cluster_annotation3  # Default case for other annotations
  ))

p_values <- data.frame(gene = character(),
                       expCond.stimuli.time = character(),
                       cell_type = character(),
                       p_value = numeric(),
                       stringsAsFactors = FALSE)

baseline <- "unstim_0h"
cell_types <-  c("Mono/Mph", "B", "CD4 T", "CD8 T", "NK", "EC", "DC")

for (gene in colnames(COA_AOA_exp)) {
  for (cell_type in cell_types) {
    # print(cell_type)
    df1 <- COA_AOA_exp_meta %>% filter(cluster_annotation_plot == cell_type)
    if (baseline %in% df1$expCond.stimuli.time) {
      for (treatment in unique(df1$expCond.stimuli.time)[2:7]) { # exclude the baseline
        df2 <- df1 %>% filter(expCond.stimuli.time %in% c(baseline, treatment))
        # print(sprintf("Gene %s | Treatment comparision %s %s", gene, baseline, treatment))
        res <- wilcox.test(
          df2[df2$expCond.stimuli.time == baseline, gene],
          df2[df2$expCond.stimuli.time == treatment, gene],
          exact = FALSE
        )
        # Store the p-value if it is not NA
        if (is.na(res$p.value)) {
          p_values <- rbind(p_values, data.frame(gene = gene, 
                                                 expCond.stimuli.time = treatment, 
                                                 cell_type = cell_type,
                                                 p_value = "NA"))
        }else{
          p_values <- rbind(p_values, data.frame(gene = gene, 
                                                 expCond.stimuli.time = treatment, 
                                                 cell_type = cell_type,
                                                 p_value = res$p.value))
        }
      }
    }
  }
}


p_values <- p_values %>%
  rename(expCond_stimuli_time = expCond.stimuli.time)

p_values <- p_values %>%
  rename(cluster_annotation_plot = cell_type)

p_values <- p_values %>%
  mutate(p_value = as.numeric(p_value),   # Ensure p_value column is numeric
         p_value_sig = case_when(
           p_value < 0.01 ~ "<0.01",      # Mark significant p-values
           TRUE ~ "Not significant"       # Mark non-significant p-values
         ))

# Calculate the average expression per cluster
COA_AOA_avg_exp <- COA_AOA_exp %>%
  rownames_to_column(var = "cell") %>%
  inner_join(rds@meta.data %>% 
               rownames_to_column(var = "cell") %>%
               filter(cluster_annotation3 %in% cell_types), 
             by = "cell") %>%
  group_by(cluster_annotation3_expCond_stimuli_time) %>%
  summarise(across(2:(ncol(COA_AOA_exp) + 1), ~ mean(.x, na.rm = TRUE)))


# Calculate the percentage of cells expressing each gene per cluster
COA_AOA_pct_exp <- COA_AOA_exp %>%
  rownames_to_column(var = "cell") %>%
  inner_join(rds@meta.data %>% 
               rownames_to_column(var = "cell") %>%
               filter(cluster_annotation3 %in% cell_types), 
             by = "cell") %>%
  group_by(cluster_annotation3_expCond_stimuli_time) %>%
  summarise(across(2:(ncol(COA_AOA_exp)+1), ~mean(. > 0) * 100))

# Combine average expression and percentage expression data
COA_AOA_avg_pct_exp <- COA_AOA_avg_exp %>%
  pivot_longer(-cluster_annotation3_expCond_stimuli_time, names_to = "gene", values_to = "expression") %>%
  inner_join(COA_AOA_pct_exp %>% pivot_longer(-cluster_annotation3_expCond_stimuli_time, names_to = "gene", values_to = "pct.exp"), 
             by = c("cluster_annotation3_expCond_stimuli_time", "gene"))

# Scale the expression data if needed
COA_AOA_avg_pct_exp <- COA_AOA_avg_pct_exp %>%
  group_by(gene) %>%
  mutate(expression = scale(expression))

COA_AOA_avg_pct_exp <- COA_AOA_avg_pct_exp %>%
  mutate(risk = case_when(
    gene %in% (genes %>% filter(!is.na(AOA_genes)))$AOA_genes & gene %in% genes$COA_genes ~ "COA_AOA", # If in both lists
    gene %in% (genes %>% filter(!is.na(AOA_genes)))$AOA_genes & !(gene %in% genes$COA_genes) ~ "AOA", # If only in AOA list
    gene %in% genes$COA_genes & !(gene %in% (genes %>% filter(!is.na(AOA_genes)))$AOA_genes) ~ "COA", # If only in COA list
    TRUE ~ "Other"  # For genes not in either list
  ))

COA_AOA_avg_pct_exp <- COA_AOA_avg_pct_exp %>%
  separate(
    cluster_annotation3_expCond_stimuli_time, 
    into = c("cluster_annotation3", "expCond_stimuli_time"), 
    sep = "\\|"
  )

COA_AOA_avg_pct_exp <- COA_AOA_avg_pct_exp %>%
  mutate(pct.exp.shape = case_when(
    pct.exp < 25 ~ "<25%",
    pct.exp >= 25 & pct.exp < 50 ~ "25-50%",
    pct.exp >= 50 & pct.exp < 75 ~ "50-75%",
    pct.exp >= 75 ~ "75-100%"
  ))


COA_AOA_avg_pct_exp <- COA_AOA_avg_pct_exp %>%
  mutate(cluster_annotation_plot = case_when(
    # cluster_annotation == "Alveolar Mph MT-positive" ~ "Alveolar Mph\nMT-positive",
    # cluster_annotation == "EC general capillary" ~ "EC general\ncapillary",
    # cluster_annotation == "Monocyte-derived Mph" ~ "Monocyte-\nderived Mph",
    TRUE ~ cluster_annotation3  # Default case for other annotations
  ))

# add p alues
COA_AOA_avg_pct_exp_p_values <- COA_AOA_avg_pct_exp %>%
  left_join(p_values, by = c("cluster_annotation_plot", 
                             "expCond_stimuli_time", 
                             "gene"))

COA_AOA_avg_pct_exp_p_values <- COA_AOA_avg_pct_exp_p_values %>%
  mutate(p_value_sig = ifelse(is.na(p_value_sig), "NA", p_value_sig))

# COA_AOA_avg_pct_exp_p_values <- COA_AOA_avg_pct_exp_p_values %>%
#   mutate(cluster_annotation_plot = str_replace_all(cluster_annotation_plot, 
#                                                    c("Alveolar Mph\nMT-positive" = "Alv Mph",
#                                                      "EC general\ncapillary" = "EC",
#                                                      "Monocyte-\nderived Mph" = "Mono Mph",
#                                                      "B cells" = "B",
#                                                      "CD4 T cells" = "CD4 T",
#                                                      "CD8 T cells" = "CD8 T",
#                                                      "NK cells" = "NK",
#                                                      "Migratory DC" = "DC")))

any(is.na(COA_AOA_avg_pct_exp_p_values)) # would be yes this time due to unstim vs unstim

#### -- Create the custom dot plot 212 ####
pdf(paste0(out.dir,"/dot_plot_COA_AOA_celltype_treament_pvalues_212.pdf"), width = 13, height = 24)

# adjust the min and max exp
p1 <- ggplot(COA_AOA_avg_pct_exp_p_values,
             aes(x = factor(expCond_stimuli_time, levels=c("unstim_0h", 
                                                           "LPS_4h", "LPS_18h",
                                                           "CD3_4h", "CD3_18h",
                                                           "Ig_4h", "Ig_18h")), 
                                               y = factor(gene, 
                                                          levels=rev(c(unique(c(genes$COA_genes, (genes %>% filter(!is.na(AOA_genes)))$AOA_genes))))),
                                               shape = pct.exp.shape, 
                                               fill = expression, 
                                               color = factor(p_value_sig, levels=c("<0.01", "Not significant", "NA")),
                                               stroke = p_value_sig)) +
  geom_point(size=2.7) +
  scale_shape_manual(values =  c("<25%" = 21, "25-50%" = 22, "50-75%" = 23, "75-100%" = 24)) +
  scale_color_manual(values =  c("<0.01" = "#0000ff", "Not significant" = "white", "NA" = "white")) +
  scale_discrete_manual(aesthetics = "stroke", 
                        values =  c("<0.01" = 1, "Not significant" = 0, "NA" = 0), 
                        guide = "none") +
  scale_fill_gradientn(colors = colorRampPalette(c('white', '#eae2b7', '#ffba08',
                                                   "#dc2f02","#9d0208",
                                                   "#370617", "#03071e"))(100),limits = c(-3.5, 7.5)) +
  facet_grid(risk ~ 
               factor(cluster_annotation_plot,
                      levels = cell_types), scales="free", space="free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "left",
        strip.background = element_rect(fill = "white", color = "black", size = 1),
        axis.text.y = element_text(size = 9),
        strip.text.y = element_blank(),
        strip.background.y = element_blank()) +
  labs(x = "Celltype, stimuli and time", 
       y = "Genes",
       shape = "Percent expressed",
       fill = "Average expression",
       color = expression(italic("p") * " value"))


p2 <- ggplot(gene_scores, aes(x = score, 
                              y = factor(gene, levels=rev(c(unique(c(genes$COA_genes, (genes %>% filter(!is.na(AOA_genes)))$AOA_genes))))), 
                              color = factor(gene_type, levels=c("COA", "AOA")),
                              group = factor(gene_type, levels=c("COA", "AOA")))) +
  # geom_line(linewidth=0.5) +
  geom_point(size=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.background = element_rect(fill = "white", color = "black", size = 1),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  scale_color_manual(values = c("#faa4bd", "#533b4d")) +
  scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, by = 1)) +
  labs(x = "Risk score", y = "Genes", Color = "Category", color = "Risk type") +
  facet_grid(risk ~ ., scales = "free", space = "free")
# theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

p1 + p2 + plot_layout(ncol = 2, widths = c(10, 1)) 


dev.off()


gene_scores %>%
  filter(gene == "BACH2")

genes %>%
  filter(COA_genes == "BACH2")


#### -- Create the custom dot plot 212 additional filters ####
selected_genes <- COA_AOA_avg_pct_exp_p_values%>%
  filter(pct.exp >= 50) %>% 
  group_by(gene) %>%  
  filter(any(p_value_sig < 0.01)) %>% 
  ungroup() %>%
  filter(expression[,1] >= 2.5) %>% 
  select(gene)

COA_AOA_avg_pct_exp_p_values_filtered  <- COA_AOA_avg_pct_exp_p_values %>%
  filter(gene %in% selected_genes$gene)

  
pdf(paste0(out.dir,"/dot_plot_COA_AOA_cellType_treament_pvalues_212_pct.exp_50_exp2.5_sig.pdf"), width = 13, height = 8)

p1 <- ggplot(COA_AOA_avg_pct_exp_p_values_filtered,
             aes(x = factor(expCond_stimuli_time, levels=c("unstim_0h", 
                                                           "LPS_4h", "LPS_18h",
                                                           "CD3_4h", "CD3_18h",
                                                           "Ig_4h", "Ig_18h")), 
                 y = factor(gene, levels=rev(c(unique(c(genes$COA_genes, (genes %>% filter(!is.na(AOA_genes)))$AOA_genes))))),
                 shape = pct.exp.shape, 
                 fill = expression, 
                 color = factor(p_value_sig, levels=c("<0.01", "Not significant", "NA")),
                 stroke = p_value_sig)) +
  geom_point(size=2.7) +
  scale_shape_manual(values =  c("<25%" = 21, "25-50%" = 22, "50-75%" = 23, "75-100%" = 24)) +
  scale_color_manual(values =  c("<0.01" = "#0000ff", "Not significant" = "white", "NA" = "white")) +
  scale_discrete_manual(aesthetics = "stroke", 
                        values =  c("<0.01" = 1, "Not significant" = 0, "NA" = 0), 
                        guide = "none") +
  scale_fill_gradientn(colors = colorRampPalette(c('white', '#eae2b7', '#ffba08',
                                                   "#dc2f02","#9d0208",
                                                   "#370617", "#03071e"))(100),limits = c(-3.5, 7.5)) +
  facet_grid(risk ~ cluster_annotation_plot, scales="free", space="free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "left",
        strip.background = element_rect(fill = "white", color = "black", size = 1),
        # axis.text.y = element_text(size = 9),
        strip.text.y = element_blank(),
        strip.background.y = element_blank(),
        plot.title = element_text(size = 12, face = "bold")) +
  labs(x = "Celltype, stimuli and time", 
       y = "Genes",
       shape = "Percent expressed",
       fill = "Average expression",
       color = expression(italic("p") * " value"),
       title = "Highly expressed (both avg exp ad % exp) genes in at least one cell type with p-value < 0.01")


p2 <- ggplot(gene_scores %>% filter(gene %in% unique(COA_AOA_avg_pct_exp_p_values_filtered$gene)), aes(x = score, 
                              y = factor(gene, levels=rev(c(unique(c(genes$COA_genes, (genes %>% filter(!is.na(AOA_genes)))$AOA_genes))))), 
                              color = factor(gene_type, levels=c("COA", "AOA")),
                              group = factor(gene_type, levels=c("COA", "AOA")))) +
  # geom_line(linewidth=0.5) +
  geom_point(size=2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.background = element_rect(fill = "white", color = "black", size = 1),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  scale_color_manual(values = c("#faa4bd", "#533b4d")) +
  scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, by = 1)) +
  labs(x = "Risk score", y = "Genes", Color = "Category", color = "Risk type") +
  facet_grid(risk ~ ., scales = "free", space = "free")
# theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

p1 + p2 + plot_layout(ncol = 2, widths = c(10, 1)) 


dev.off()

#### -- venn diagram ####
COA_AOA_avg_pct_exp_p_values %>%
  group_by(risk) %>%
  summarise(gene_count = n_distinct(gene))


aoa_cnts1 <- COA_AOA_avg_pct_exp_p_values %>%
  group_by(risk) %>%
  summarise(gene_count = n_distinct(gene)) %>% 
  filter(risk == "AOA") %>% pull(gene_count)
coa_cnts1 <- COA_AOA_avg_pct_exp_p_values %>%
  group_by(risk) %>%
  summarise(gene_count = n_distinct(gene)) %>% 
  filter(risk == "COA") %>% pull(gene_count)
aoa_coa_cnts1 <- COA_AOA_avg_pct_exp_p_values %>%
  group_by(risk) %>%
  summarise(gene_count = n_distinct(gene)) %>% 
  filter(risk == "COA_AOA") %>% pull(gene_count)  

aoa_cnts2 <- gene_scores %>%
  group_by(risk) %>%
  summarise(gene_count = n_distinct(gene)) %>% 
  filter(risk == "AOA") %>% pull(gene_count)
coa_cnts2 <- gene_scores %>%
  group_by(risk) %>%
  summarise(gene_count = n_distinct(gene)) %>% 
  filter(risk == "COA") %>% pull(gene_count)
aoa_coa_cnts2 <- gene_scores %>%
  group_by(risk) %>%
  summarise(gene_count = n_distinct(gene)) %>% 
  filter(risk == "COA_AOA") %>% pull(gene_count)  

venn1 <- draw.pairwise.venn(
  area1 = aoa_cnts1 + aoa_coa_cnts1,
  area2 = coa_cnts1 + aoa_coa_cnts1,
  cross.area = aoa_coa_cnts1,
  category = c("AOA", "COA"),
  fill = c("#533b4d", "#faa4bd"),
  lty = "blank",
  alpha = 0.5,
  cat.cex = 1.2,  # Category label size
  cex = 1.2,  # Number size
  cat.pos = c(-10, 10),  # Category position
  cat.dist = 0.05  # Category distance from the circles
)

venn2 <- draw.pairwise.venn(
  area1 = aoa_cnts2 + aoa_coa_cnts2,
  area2 = coa_cnts2 + aoa_coa_cnts2,
  cross.area = aoa_coa_cnts2,
  category = c("AOA", "COA"),
  fill = c("#533b4d", "#faa4bd"),
  lty = "blank",
  alpha = 0.5,
  cat.cex = 1.2,  # Category label size
  cex = 1.2,  # Number size
  cat.pos = c(-10, 10),  # Category position
  cat.dist = 0.05  # Category distance from the circles
  )

pdf(paste0(out.dir,"/venn_COA_AOA_212_201.pdf"), width = 6.13, height = 2.39)

grid.arrange(venn2, venn1, ncol = 2)

dev.off()

#### -- How many asthmatic genes are DEGs from DE analysis? ####
deResDir <- "integration/integration_1/integration_1_louvain/"
# Load DEGs for each contrast 
excel.dir <- "results_wOrgClusterAnnotation_DEGs/stimuli_time_MAST/"
file.all <- list.files(paste0(deResDir,excel.dir))
files <- file.all[grepl("full_allClusters.xlsx", file.all, ignore.case = TRUE)]
contrast <- gsub("expCondCompDeMarkers_|_full_allClusters.xlsx", "", files)

# save each cell type to a tibble for each comparison 
for (i in 1:length(files)) {
  excel_file <- paste0(deResDir, excel.dir, files[i])
  sheet_names <- excel_sheets(excel_file)
  for (sheet_name in sheet_names) {
    assign(paste0(contrast[i], "_", sheet_name), read_excel(excel_file, sheet = sheet_name))
  }
}


getUpSetInput <- function(degtbl = "AllCellTypes$", logFC.cutoff = 0.25, fdr.cutoff = 0.05) {
  grp <- ls(pattern = degtbl, envir = parent.frame())
  
  # filter by threshold and get a list of deduplicated gene names for each 
  sig.degs.list <- lapply(grp, function(x) {
    degs <- get(x)
    sig.degs      <- degs %>% dplyr::filter(p_val_adj<=fdr.cutoff & abs(avg_log2FC)>=logFC.cutoff) %>% dplyr::arrange(desc(avg_log2FC))
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


# regardless up or down regualted, put degs together
celltypes <- c("Alveolar Mph MT-positive", 
              "EC general capillary", 
              "Monocyte-derived Mph",
              "B cells",
              "CD4 T cells",
              "CD8 T cells",
              "NK cells",
              "Migratory DC")



# length(upsetInput$up$`CD3_18h/unstim_0h`) + length(upsetInput$dn$`CD3_18h/unstim_0h`)
# length(upsetInput$both$`CD3_18h/unstim_0h`)


# Create an empty data frame to store results
res <- data.frame(
  gene = character(),
  condition = character(),
  celltype = character(),
  stringsAsFactors = FALSE
)

# Loop over each cell type
for (celltype in celltypes) {
  # Clean the cell type name
  ct <- gsub("/| |-","",celltype)
  
  # Get the UpSet input for the cell type
  upsetInput <- getUpSetInput(ct)
  
  # Loop over the conditions in the upsetInput$both list
  for (condition in names(upsetInput$both)) {
    # Get the list of genes for the condition
    genes <- upsetInput$both[[condition]]
    
    # Check if there are any genes to add
    if (length(genes) > 0) {
      # Create a temporary data frame with the genes, condition, and cell type
      tmp <- data.frame(
        gene = genes,
        condition = condition,
        celltype = celltype,
        stringsAsFactors = FALSE
      )
      
      # Append the temporary data frame to the results data frame
      res <- rbind(res, tmp)
    } else {
      # Print the celltype and condition if no genes are present
      print(paste("No genes for celltype:", ct, "condition:", condition))
    }
  }
}

# Display the final results data frame
print(res)

# spot check one line
res %>%
  filter(condition %in% "Ig_18h/unstim_0h",
         celltype%in% "Mono/Mph") %>%
  dim()

degs_COA_AOA <- inner_join(res, gene_scores, 
                           by = "gene", relationship = "many-to-many")


length((degs_COA_AOA$gene))
length(unique(degs_COA_AOA$gene))

cat("AOA_COA genes that are DEGs: ", length(unique(degs_COA_AOA$gene)), "\n")
#AOA_COA genes that are DEGs:  136 

# Asthmatic genes that are DEGs & DEGs in one cell type
degs_AOA_COA_multiple_celltypes <- degs_COA_AOA %>%
  group_by(gene) %>%
  summarize(celltype_count = n_distinct(celltype)) %>%
  filter(celltype_count > 1) 

cat("Number of genes present in more than one cell type: ", nrow(degs_AOA_COA_multiple_celltypes), "\n")
# Number of genes present in more than one cell type:  99 

degs_AOA_COA_multiple_conditions <- degs_COA_AOA %>%
  group_by(celltype, gene) %>%
  summarize(condition_count = n_distinct(condition)) %>%
  filter(condition_count > 1) 

cat("Number of genes present in more than one condition within the same cell type: ", length(unique(degs_AOA_COA_multiple_conditions$gene)), "\n")
# Number of genes present in more than one condition within the same cell type:  112 


#### -- How many asthmatic genes are DEGs from the plot itself #### 

degs_COA_AOA <- COA_AOA_avg_pct_exp_p_values[c('cluster_annotation', 
                                               "expCond_stimuli_time",
                                               'gene',
                                               "cluster_annotation_plot",
                                               "p_value",
                                               "p_value_sig")] %>%
  filter(p_value_sig %in% "<0.01")


cat("AOA_COA genes that are DEGs: ", length(unique(degs_COA_AOA$gene)), "\n")
#AOA_COA genes that are DEGs:  192 

# Asthmatic genes that are DEGs & DEGs in one cell type
degs_AOA_COA_multiple_celltypes <- degs_COA_AOA %>%
  group_by(gene) %>%
  summarize(celltype_count = n_distinct(cluster_annotation_plot)) %>%
  filter(celltype_count > 1) 

cat("Number of genes present in more than one cell type: ", nrow(degs_AOA_COA_multiple_celltypes), "\n")
# Number of genes present in more than one cell type:  173 


# run the HLA below

degs_hla <- hla_avg_pct_exp_p_values[c('cluster_annotation', 
                                               "expCond_stimuli_time",
                                               'gene',
                                               "cluster_annotation_plot",
                                               "p_value",
                                               "p_value_sig")] %>%
  filter(p_value_sig %in% "<0.01")


cat("HLA genes that are DEGs: ", length(unique(degs_hla$gene)), "\n")
#AOA_COA genes that are DEGs:  24 

# Asthmatic genes that are DEGs & DEGs in one cell type
degs_hla_multiple_celltypes <- degs_hla %>%
  group_by(gene) %>%
  summarize(celltype_count = n_distinct(cluster_annotation_plot)) %>%
  filter(celltype_count > 1) 

cat("Number of genes present in more than one cell type: ", nrow(degs_hla_multiple_celltypes), "\n")
# Number of genes present in more than one cell type:  24


# together
degs_COA_AOA_hla <- rbind(degs_COA_AOA, degs_hla)

cat("HLA genes that are DEGs: ", length(unique(degs_COA_AOA_hla$gene)), "\n")
#AOA_COA genes that are DEGs:  216 

# Asthmatic genes that are DEGs & DEGs in one cell type
degs_COA_AOA_hla_multiple_celltypes <- degs_COA_AOA_hla %>%
  group_by(gene) %>%
  summarize(celltype_count = n_distinct(cluster_annotation_plot)) %>%
  filter(celltype_count > 1) 

cat("Number of genes present in more than one cell type: ", nrow(degs_COA_AOA_hla_multiple_celltypes), "\n")
# Number of genes present in more than one cell type:  197


COA_AOA_hla <- rbind(COA_AOA_avg_pct_exp_p_values[c('cluster_annotation', 
                                              "expCond_stimuli_time",
                                              'gene',
                                              "cluster_annotation_plot",
                                              "p_value",
                                              "p_value_sig")],
                     hla_avg_pct_exp_p_values[c('cluster_annotation', 
                                                "expCond_stimuli_time",
                                                'gene',
                                                "cluster_annotation_plot",
                                                "p_value",
                                                "p_value_sig")])
length(unique(COA_AOA_hla$gene))   #225                 
length(unique(COA_AOA_avg_pct_exp_p_values$gene)) #201
length(unique(hla_avg_pct_exp_p_values$gene)) #24



#### HLA genes ####
genes_hla <- read_excel(paste0(input.dir, "HLA_Genes.xlsx"), col_names = FALSE, sheet = 1)
colnames(genes_hla) <- c("gene")
  
  
Seurat::DefaultAssay(rds) <- "RNA"
rds@meta.data$cluster_annotation_expCond_stimuli_time <- paste0(
  rds@meta.data$cluster_annotation, 
  "|", 
  rds@meta.data$expCond.stimuli.time
)

# Extract data for the specified features
missing_genes <- setdiff(unique(genes_hla$gene), rownames(GetAssayData(object = rds, layer = "data")))

# Display missing genes, if any
if(length(missing_genes) > 0) {
  print("These genes are missing:")
  print(missing_genes)
} else {
  print("All genes are present.")
}

# [1] "These genes are missing:"
# [1] "HLA-DRB3" "HLA-DRB4"

# all genes that start with HLA
# Filter(function(x) grepl("^HLA", x), rownames(GetAssayData(object = rds, layer = "data")))

existing_genes <- intersect(unique(genes_hla$gene), rownames(GetAssayData(object = rds, layer = "data")))

hla_exp <- as.data.frame(t(GetAssayData(object = rds, layer = "data")[existing_genes, ]))


# Calculate some stats
hla_exp_meta <- hla_exp %>%
  rownames_to_column(var = "cell") %>%
  inner_join(rds@meta.data %>% rownames_to_column(var = "cell"), by = "cell")

any(is.na(hla_exp_meta)) # any nas? 

hla_exp_meta <- hla_exp_meta %>%
  mutate(cluster_annotation_plot = case_when(
    cluster_annotation == "Alveolar Mph MT-positive" ~ "Alveolar Mph\nMT-positive",
    cluster_annotation == "EC general capillary" ~ "EC general\ncapillary",
    cluster_annotation == "Monocyte-derived Mph" ~ "Monocyte-\nderived Mph",
    TRUE ~ cluster_annotation  # Default case for other annotations
  ))

p_values <- data.frame(gene = character(),
                       expCond.stimuli.time = character(),
                       cell_type = character(),
                       p_value = numeric(),
                       stringsAsFactors = FALSE)

baseline <- "unstim_0h"


for (gene in colnames(hla_exp)) {
  for (cell_type in unique(hla_exp_meta$cluster_annotation_plot)) {
    # print(cell_type)
    df1 <- hla_exp_meta %>% filter(cluster_annotation_plot == cell_type)
    if (baseline %in% df1$expCond.stimuli.time) {
      for (treatment in unique(df1$expCond.stimuli.time)[2:7]) { # exclude the baseline
        df2 <- df1 %>% filter(expCond.stimuli.time %in% c(baseline, treatment))
        # print(sprintf("Gene %s | Treatment comparision %s %s", gene, baseline, treatment))
        res <- wilcox.test(df2[df2$expCond.stimuli.time == baseline, gene],
          df2[df2$expCond.stimuli.time == treatment, gene],
          exact = FALSE)
        if (is.na(res$p.value)) {
          p_values <- rbind(p_values, data.frame(gene = gene, 
                                                 expCond.stimuli.time = treatment, 
                                                 cell_type = cell_type,
                                                 p_value = "NA"))
        }else{
          p_values <- rbind(p_values, data.frame(gene = gene, 
                                                 expCond.stimuli.time = treatment, 
                                                 cell_type = cell_type,
                                                 p_value = res$p.value))
        }
      }
    }
  }
}




p_values <- p_values %>%
  rename(expCond_stimuli_time = expCond.stimuli.time)

p_values <- p_values %>%
  rename(cluster_annotation_plot = cell_type)

p_values <- p_values %>%
  mutate(p_value = as.numeric(p_value),   # Ensure p_value column is numeric
         p_value_sig = case_when(
           p_value < 0.01 ~ "<0.01",      # Mark significant p-values
           TRUE ~ "Not significant"       # Mark non-significant p-values
         ))


# Calculate the average expression per cluster
hla_avg_exp <- hla_exp %>%
  rownames_to_column(var = "cell") %>%
  inner_join(rds@meta.data %>% rownames_to_column(var = "cell"), by = "cell") %>%
  group_by(cluster_annotation_expCond_stimuli_time) %>%
  summarise(across(2:(ncol(hla_exp) + 1), ~ mean(.x, na.rm = TRUE)))


# Calculate the percentage of cells expressing each gene per cluster
hla_pct_exp <- hla_exp %>%
  rownames_to_column(var = "cell") %>%
  inner_join(rds@meta.data %>% rownames_to_column(var = "cell"), by = "cell") %>%
  group_by(cluster_annotation_expCond_stimuli_time) %>%
  summarise(across(2:(ncol(hla_exp)+1), ~mean(. > 0) * 100))

# Combine average expression and percentage expression data
hla_avg_pct_exp <- hla_avg_exp %>%
  pivot_longer(-cluster_annotation_expCond_stimuli_time, names_to = "gene", values_to = "expression") %>%
  inner_join(hla_pct_exp %>% pivot_longer(-cluster_annotation_expCond_stimuli_time, names_to = "gene", values_to = "pct.exp"), 
             by = c("cluster_annotation_expCond_stimuli_time", "gene"))

# Scale the expression data if needed
hla_avg_pct_exp <- hla_avg_pct_exp %>%
  group_by(gene) %>%
  mutate(expression = scale(expression))

hla_avg_pct_exp <- hla_avg_pct_exp %>%
  separate(
    cluster_annotation_expCond_stimuli_time, 
    into = c("cluster_annotation", "expCond_stimuli_time"), 
    sep = "\\|"
  )

hla_avg_pct_exp <- hla_avg_pct_exp %>%
  mutate(pct.exp.shape = case_when(
    pct.exp < 25 ~ "<25%",
    pct.exp >= 25 & pct.exp < 50 ~ "25-50%",
    pct.exp >= 50 & pct.exp < 75 ~ "50-75%",
    pct.exp >= 75 ~ "75-100%"
  ))


hla_avg_pct_exp <- hla_avg_pct_exp %>%
  mutate(cluster_annotation_plot = case_when(
    cluster_annotation == "Alveolar Mph MT-positive" ~ "Alveolar Mph\nMT-positive",
    cluster_annotation == "EC general capillary" ~ "EC general\ncapillary",
    cluster_annotation == "Monocyte-derived Mph" ~ "Monocyte-\nderived Mph",
    TRUE ~ cluster_annotation  # Default case for other annotations
  ))

# add p alues
hla_avg_pct_exp_p_values <- hla_avg_pct_exp %>%
  left_join(p_values, by = c("cluster_annotation_plot", 
                             "expCond_stimuli_time", 
                             "gene"))

hla_avg_pct_exp_p_values <- hla_avg_pct_exp_p_values %>%
  mutate(p_value_sig = ifelse(is.na(p_value_sig), "NA", p_value_sig))

hla_avg_pct_exp_p_values <- hla_avg_pct_exp_p_values %>%
  mutate(cluster_annotation_plot = str_replace_all(cluster_annotation_plot, 
                                                   c("Alveolar Mph\nMT-positive" = "Alv Mph",
                                                     "EC general\ncapillary" = "EC",
                                                     "Monocyte-\nderived Mph" = "Mono Mph",
                                                     "B cells" = "B",
                                                     "CD4 T cells" = "CD4 T",
                                                     "CD8 T cells" = "CD8 T",
                                                     "NK cells" = "NK",
                                                     "Migratory DC" = "DC")))

any(is.na(hla_avg_pct_exp_p_values)) # would be yes this time due to unstim vs unstim



#### -- Create the custom dot plot ####
# adjust the min and max exp

pdf(paste0(out.dir,"/dot_plot_hla_pvalues.pdf"), width = 10, height = 4.97)

ggplot(hla_avg_pct_exp_p_values,
             aes(x = factor(expCond_stimuli_time, levels=c("unstim_0h", 
                                                           "LPS_4h", "LPS_18h",
                                                           "CD3_4h", "CD3_18h",
                                                           "Ig_4h", "Ig_18h")), 
                 y = factor(gene, levels=rev(genes_hla$gene)),
                 shape = pct.exp.shape, 
                 fill = expression, 
                 color = factor(p_value_sig, levels=c("<0.01", "Not significant", "NA")),
                 stroke = p_value_sig)) +
  geom_point(size=2.7) +
  scale_shape_manual(values =  c("<25%" = 21, "25-50%" = 22, "50-75%" = 23, "75-100%" = 24)) +
  scale_color_manual(values =  c("<0.01" = "#0000ff", "Not significant" = "white", "NA" = "white")) +
  scale_discrete_manual(aesthetics = "stroke", 
                        values =  c("<0.01" = 1, "Not significant" = 0, "NA" = 0), 
                        guide = "none") +
  scale_fill_gradientn(colors = colorRampPalette(c('white', '#eae2b7', '#ffba08',
                                                   "#dc2f02","#9d0208",
                                                   "#370617", "#03071e"))(100),limits = c(-3.5, 7.5)) +
  facet_grid(. ~ cluster_annotation_plot, scales="free", space="free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "left",
        strip.background = element_rect(fill = "white", color = "black", linewidth = 1),
        axis.text.y = element_text(size = 9),
        strip.text.y = element_blank(),
        strip.background.y = element_blank()) +
  labs(x = "Celltype, stimuli and time", 
       y = "Genes",
       shape = "Percent expressed",
       fill = "Average expression",
       color = expression(italic("p") * " value"))

dev.off()




