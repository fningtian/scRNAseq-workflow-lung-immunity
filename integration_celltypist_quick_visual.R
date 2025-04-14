# setwd("/Users/funingtian/Dropbox/UChicago/CRI_ftian1/scRNAseq")
setwd("/gpfs/data/schoettler-lab/ftian1/robi/integration_2")

.libPaths("/ess/home/home1/ftian1/R/x86_64-pc-linux-gnu-library/4.2")
.libPaths("/gpfs/data/icelake-apps/software/gcc-12.1.0/R/4.2.1/lib64/R/library")

library(Seurat)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(data.table)
library(patchwork)


#### run on randi | 1.0 extract metadata from seurat objects ####
seuratObj <- readRDS('integration_2_leiden/RDS_Dir/integration_2_leiden.rds')

# export metadata
write.csv(as.matrix(seuratObj@meta.data), 
          file = 'integration_2_leiden/metadata/integration_2_leiden_meta.data.csv', quote = F)

# export the data layer
write.csv(as.matrix(seuratObj@assays$RNA$counts), file = "integration_2_leiden/metadata/integration_2_leiden_rna_counts.csv", quote = FALSE)

#### run on randi | 2.0 celltypist ####

#### 3.0 analyze celltypist predictions ####

setwd("/Users/funingtian/Dropbox/UChicago/CRI_ftian1/scRNAseq")

library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


pred <- read.csv("integration_2/integration_2_leiden/metadata/integration_2_leiden_meta.data.celltypist.csv")


# stacked bar chart showing the percentage of each predicted_labels for each seurat_clusters
# if a label accounts for less than 5% of a seurat_cluster, make it as "Underrepresented"
pred_clts <- pred %>%
  group_by(seurat_clusters, predicted_labels) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  mutate(predicted_labels = ifelse(percentage < 5, "Underrepresented", predicted_labels)) %>%
  group_by(seurat_clusters, predicted_labels) %>%
  summarise(
    count = sum(count),  # Combine counts for "Underrepresented"
    percentage = sum(percentage),  # Sum percentages for "Underrepresented"
    .groups = "drop"
  )

# Create the stacked bar chart
p1 <- ggplot(pred_clts, aes(x = factor(seurat_clusters), y = percentage, fill = predicted_labels)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    x = "Seurat clusters",
    y = "Percentage",
    fill = "Cell types",
  ) +
  scale_fill_manual(values = c(
    "#d8f8ff",  # Green
    "#BEAED4",  # Purple
    "#FDC086",  # Orange
    "#386CB0",  # Blue
    "#F0027F",  # Pink
    "#BF5B17",  # Brown
    "#1B9E77",  # Teal
    "#D95F02",  # Dark Orange
    "#7570B3",  # Violet
    '#74a287',  # Light Green
    "#E6AB02",  # Mustard Yellow
    "grey"   # Earthy Brown
  )) + 
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1)
  )

p2 <- ggplot(pred_clts, aes(x = factor(seurat_clusters), y = count, fill = predicted_labels)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    x = "Seurat clusters",
    y = "Frequency",
    fill = "Cell types",
  ) +
  scale_fill_manual(values = c(
    "#d8f8ff",  # Green
    "#BEAED4",  # Purple
    "#FDC086",  # Orange
    "#386CB0",  # Blue
    "#F0027F",  # Pink
    "#BF5B17",  # Brown
    "#1B9E77",  # Teal
    "#D95F02",  # Dark Orange
    "#7570B3",  # Violet
    '#74a287',  # Light Green
    "#E6AB02",  # Mustard Yellow
    "grey"   # Earthy Brown
  )) + 
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5))

pdf(paste0("integration_2/integration_2_leiden/metadata/celltypist_seurat_clts.pdf"), width = 8.11, height = 5.09)

p1 / p2 + plot_layout(guides = "collect") &
  theme(legend.position = "right")
dev.off()

pred_clts %>%
  filter(seurat_clusters %in% 15)



#### 4.0 add the metadata ####
setwd("/Users/funingtian/Dropbox/UChicago/CRI_ftian1/scRNAseq")

seuratObj <- readRDS('integration_2/integration_2_leiden/RDS_Dir/integration_2_leiden.rds')
metadata <- read.csv("integration_2/integration_2_leiden/metadata/integration_2_leiden_meta.data.celltypist.csv")


duplicate_cells_meta <- which(duplicated(metadata$Unnamed..0))
if (length(duplicate_cells_meta) > 0) {
  print("Duplicate cell names found in metadata:")
  print(metadata$Unnamed..0[duplicate_cells_meta])
}

rownames(metadata) <- metadata$Unnamed..0
metadata$Unnamed..0 <- NULL

# update the metadata
seuratObj@meta.data <- metadata


seuratObj_annot <- RenameIdents(seuratObj,
                     ## ---
                     `2`  = 'B',
                     `3`  = 'B',
                     `13` = 'B',
                     ## ---
                     `9`  = 'NK',
                     ## ---
                     `6`  = 'CD4 T',
                     `12`  = 'CD4 T',
                     `15` = "CD4 T",
                     `5` = 'CD8 T',
                     ## ---
                     `16` = 'DC',
                     ## ---
                     `1`  = 'Mono Mph',
                     `4`  = 'Mono Mph',
                     `7`  = 'Mono Mph',
                     `8`  = 'Mono Mph',
                     `10`  = 'Mono Mph',
                     `11` = 'Mono Mph',
                     `14` = 'Mono Mph',
                     `17`  = 'Mono Mph'
)

renamed_idents <- Idents(seuratObj_annot)
renamed_idents_df <- data.frame(cell = names(renamed_idents), renamed_identity = renamed_idents)

# add the renamed identities to meta.data
seuratObj_annot@meta.data$cluster_annotation <- renamed_idents_df$renamed_identity[match(rownames(seuratObj_annot@meta.data), renamed_idents_df$cell)]

seuratObj_annot@meta.data <- seuratObj_annot@meta.data %>%
  mutate(across(
    c(expCond.time, expCond.stimuli.time, expCond.donor.stimuli.time),
    ~ gsub("RPMI|CTS", "0h", .)
  ))

saveRDS(seuratObj_annot, 'integration_2/integration_2_leiden/RDS_Dir/integration_2_leiden_annot.rds')


#### 5.0 customize tsne and umap visualizations ####
# load the saved intergrated seurat object 

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
workdir <- "integration_2/integration_2_leiden"

# umap clusters with and without labeling and split by donors 
umapCluster <- DimPlot(seuratObj_annot, reduction = "umap", cols = selectedCol, label = T, label.size = 6, repel = T) + 
  labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')
umapClusterNolabel <- DimPlot(seuratObj_annot, reduction = "umap", cols = selectedCol, label = F, repel = T) + 
  labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')
umapSplit <- DimPlot(seuratObj_annot, reduction = "umap", cols = selectedCol, label = T, 
                     label.size = 4, repel = T, split.by = 'expCond.donor') +
  labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')


pdf(file = paste0(workdir, "/umap_plot_noLabel_integrate_orgAnnotation.pdf"), 
    width = 5.7, height = 7)
print(umapClusterNolabel + theme1wLegend + guides(colour = guide_legend(nrow=4, byrow=TRUE, override.aes = list(size=6))))
dev.off()

pdf(file = paste0(workdir, "/umap_plot_wLabel_integrate_orgAnnotation.pdf"), 
    width = 5.7, height = 7)
print(umapCluster + theme1wLegend + guides(colour = guide_legend(nrow=4, byrow=TRUE, override.aes = list(size=6))) )
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

