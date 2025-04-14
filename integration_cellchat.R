.libPaths(c("/ess/home/home1/ftian1/R/x86_64-pc-linux-gnu-library/4.2", 
            "/gpfs/data/icelake-apps/software/gcc-12.1.0/R/4.2.1/lib64/R/library"))


library(Seurat)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
# devtools::install_github("jinworks/CellChat")
library(igraph) # if on randi 
library(CellChat)
library(patchwork)
library(tidyr) #function split
library(ComplexHeatmap)
library(purrr) # stacked violin plot 
library(readxl)
library(forcats) #fct_reorder
library(stringr) # str detect
library(cowplot)
library(tibble) # column to rownames

options(stringsAsFactors = FALSE)


# tutorial to follow 
# https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html


#### DONE 1.0 run on randi | recreate the cellchat object based on de genes ####

dir = "/gpfs/data/schoettler-lab/ftian1/robi/integration_2/integration_2_leiden/results_wOrgClusterAnnotation_DEGs/stimuli_time_MAST"

files <- list.files(path = dir, pattern = "_adjSig_7SelClusters.xlsx$", full.names = TRUE)

# Function to read the first column of each sheet in a file
read_first_column_from_sheets <- function(file) {
  sheet_names <- excel_sheets(file) # Get sheet names
  all_genes <- map(sheet_names, function(sheet) {
    read_excel(file, sheet = sheet) %>%
      select(1) %>%
      pull()
  }) %>% 
    unlist() # Combine all vectors into one
  return(all_genes)
}

# Read and combine all gene lists from all files, then get unique genes
all_genes <- files %>%
  map(read_first_column_from_sheets) %>%
  unlist() %>%
  unique()
# Display the unique gene list
# print(all_genes)

gtf_path <- "/gpfs/data/schoettler-lab/software/cellranger/refdata-gex-GRCh38-2024-A/genes/protein_coding_genes.gtf"
gtf <- fread(gtf_path, sep = "\t", header = FALSE)

colnames(gtf) <- c("seqname", "source", "feature", "start", "end", 
                   "score", "strand", "frame", "attribute")

gtf[, gene_id := sub('.*gene_id "([^"]+)".*', '\\1', attribute)]
gtf[, gene_name := sub('.*gene_name "([^"]+)".*', '\\1', attribute)]
gtf[, gene_type := sub('.*gene_type "([^"]+)".*', '\\1', attribute)]

# gtf <- gtf[feature == "gene" & gtf == "protein_coding"]

pcgs <- unique(gtf$gene_name)


print(paste0("Number of DEGs: ", length(all_genes)))

all_genes <- intersect(all_genes, pcgs)

print(paste0("Number of protein-coding DEGs: ", length(all_genes)))


seuratObj_annot <- readRDS('integration_2/integration_2_leiden/RDS_Dir/integration_2_leiden_annot.rds')


seuratObj_annot <- subset(seuratObj_annot, subset = cluster_annotation3 %in% c("B", "CD4 T", "CD8 T", "NK", "Mono/Mph", "EC", "DC"))
seuratObj_annot$cluster_annotation3 <- factor(
  seuratObj_annot$cluster_annotation3,
  levels = c("B", "CD4 T", "CD8 T", "NK", "Mono/Mph", "EC", "DC")
)

seuratObj_annot@meta.data$cluster_annotation3.expCond.stimuli.time <- paste(seuratObj_annot@meta.data$cluster_annotation3, 
                                                                           seuratObj_annot@meta.data$expCond.stimuli.time, sep = "|")

# data.input <- seuratObj_annot[["RNA"]]$data # normalized data matrix
# For Seurat version >= “5.0.0”, get the normalized data via `seurat_object[["RNA"]]$data`
# labels <- Idents(seuratObj_annot)
# meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellchat <- createCellChat(object = seuratObj_annot, group.by = "cluster_annotation3.expCond.stimuli.time", assay = "RNA")


CellChatDB <- CellChatDB.human

# showDatabaseCategory(CellChatDB)
# 
# dplyr::glimpse(CellChatDB$interaction)
# table(CellChatDB$interaction$version)
# table(CellChatDB$interaction$annotation)

# CellChatDB.use <- subsetDB(CellChatDB) ##-- need to look into this 
# simply use the default CellChatDB. We do not suggest to use it in this way because CellChatDB v2 
# includes "Non-protein Signaling" (i.e., metabolic and synaptic signaling). 

# set the used database in the object
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling", "Cell-Cell Contact"))
# table(CellChatDB.use$interaction$annotation)
# table(CellChatDB.use$interaction$annotation)

# set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
#futurhttps://rpubs.com/marswh/boxplote::plan("multisession", workers = 4) # do parallel

features <- intersect(all_genes, row.names(cellchat@data.signaling))

print(length(all_genes)) #
print(length(row.names(cellchat@data.signaling))) #
print(length(features)) #


# cellchat <- identifyOverExpressedGenes(cellchat) # skip
features.name = "features"
features.sig <- features
cellchat@var.features[[features.name]] <- unique(features.sig)
features.name <- paste0(features.name, ".info")
cellchat@var.features[[features.name]] <- unique(features.sig)


cellchat <- identifyOverExpressedInteractions(cellchat)

# inference of cell-cell communication network
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# compute network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
# cellchat <- computeNetSimilarity(cellchat, type = "functional")
# cellchat <- netEmbedding(cellchat, type = "functional",umap.method = "uwot")
# cellchat <- netClustering(cellchat, type = "functional") # this will be run locally

saveRDS(cellchat, file="integration_2/integration_2_leiden/RDS_Dir/integration_2_leiden_annot_cellchat_de.rds")






# #### SKIP 2.0 run on local server | subsetting cellchat object ####
# cellchat <- readRDS('integration/integration_1/integration_1_louvain/RDS_Dir/integration_1_louvain_annot_cellchat_de.rds')
# 
# tmp <- cellchat
# 
# cell_types_reset <- c("Mono/Mph",
#                       "B", 
#                       "CD4 T", 
#                       "CD8 T", 
#                       "NK",
#                       "EC", 
#                       "DC")
# 
# 
# cell_types_reset_level <- c(rev(grep(cell_types_reset[1], levels(tmp@idents), value = T)),
#                             rev(grep(cell_types_reset[2], levels(tmp@idents), value = T)),
#                             rev(grep(cell_types_reset[3], levels(tmp@idents), value = T)),
#                             rev(grep(cell_types_reset[4], levels(tmp@idents), value = T)),
#                             rev(grep(cell_types_reset[5], levels(tmp@idents), value = T)),
#                             rev(grep(cell_types_reset[6], levels(tmp@idents), value = T)),
#                             rev(grep(cell_types_reset[7], levels(tmp@idents), value = T))
#                             )
# 
# levels(tmp@idents) <- cell_types_reset_level # this does not update net and netP
# 
# whole_color_use_doughnut <- c(
#   colorRampPalette(colors = rev(c("#eae2b7", "#f77f00")))(length(grep(cell_types[1], unique(tmp@idents), value = T))),
#   colorRampPalette(colors = rev(c("#C6DBEF", "#075a84")))(length(grep(cell_types[2], unique(tmp@idents), value = T))), 
#   colorRampPalette(colors = rev(c("#dfd3e8", "#9163cb")))(length(grep(cell_types[3], unique(tmp@idents), value = T))), 
#   colorRampPalette(colors = rev(c("#ffd6e0", "#ff4d6d")))(length(grep(cell_types[4], unique(tmp@idents), value = T))), 
#   colorRampPalette(colors = rev(c("#d7e0ab", "#679436")))(length(grep(cell_types[5], unique(tmp@idents), value = T))), 
#   colorRampPalette(colors = rev(c("#fed9b7", "#f07167")))(length(grep(cell_types[6], unique(tmp@idents), value = T))),
#   colorRampPalette(colors = rev(c("#c2b5a7", "#603808")))(length(grep(cell_types[7], unique(tmp@idents), value = T))),
#   colorRampPalette(colors = rev(c("#b2f7ef", "#07beb8")))(length(grep(cell_types[8], unique(tmp@idents), value = T)))
# )
# 
# whole_color_use2_doughnut <- c(
#   colorRampPalette(colors = c("#eae2b7", "#f77f00"))(length(grep(cell_types[1], unique(tmp@idents), value = T))),
#   colorRampPalette(colors = c("#C6DBEF", "#075a84"))(length(grep(cell_types[2], unique(tmp@idents), value = T))), 
#   colorRampPalette(colors = c("#dfd3e8", "#9163cb"))(length(grep(cell_types[3], unique(tmp@idents), value = T))), 
#   colorRampPalette(colors = c("#ffd6e0", "#ff4d6d"))(length(grep(cell_types[4], unique(tmp@idents), value = T))), 
#   colorRampPalette(colors = c("#d7e0ab", "#679436"))(length(grep(cell_types[5], unique(tmp@idents), value = T))), 
#   colorRampPalette(colors = c("#fed9b7", "#f07167"))(length(grep(cell_types[6], unique(tmp@idents), value = T))),
#   colorRampPalette(colors = c("#c2b5a7", "#603808"))(length(grep(cell_types[7], unique(tmp@idents), value = T))),
#   colorRampPalette(colors = c("#b2f7ef", "#07beb8"))(length(grep(cell_types[8], unique(tmp@idents), value = T)))
# )
# 
# 
# # name the colors based on the cell types
# names(whole_color_use_doughnut) <- unlist(lapply(cell_types_reset, function(ct) grep(ct, levels(tmp@idents), value = T)))
# names(whole_color_use2_doughnut) <- unlist(lapply(cell_types_reset, function(ct) grep(ct, levels(tmp@idents), value = T)))
# 
# groupSize <- as.numeric(table(tmp@idents))
# 
# counts <- tmp@meta %>% 
#   group_by(cluster_annotation.expCond.stimuli.time) %>% 
#   count() 
# 
# pdf(paste0("./integration/integration_1/integration_1_louvain/cellchat_de/doughnut_cell_cnts.pdf"), 
#     width = 9, height = 6)
# 
# # quick doughnut plot for the %
# ggplot(counts, aes(x = 2, y = n, fill = factor(cluster_annotation.expCond.stimuli.time, levels=cell_types_reset_level))) +
#   geom_bar(stat = "identity", width = 1, color = "white") +
#   coord_polar(theta = "y") +    # Make it circular (polar coordinates)
#   xlim(0.5, 2.5) +              # Create space for the hole in the middle
#   theme_void() +                # Remove all background, gridlines, and axes
#   # theme(legend.position = "none") +  # Position the legend
#   labs(fill = "Cluster Annotation") +  # Label the legend
#   scale_fill_manual(values = whole_color_use2_doughnut) +
#   guides(fill = guide_legend(title = "Cluster Annotation")) +
#   geom_text(aes(label = ifelse(n > 3000, n, "")), position = position_stack(vjust = 0.5))  # Add counts as labels
# 
# dev.off()
# 
# exclude_cells <- c("Alveolar Mph MT-positive|LPS_4h",
#                    "Alveolar Mph MT-positive|LPS_18h")
# 
# cellchat <- subsetCellChat(cellchat, idents.use = unique(cellchat@idents)[!unique(cellchat@idents) %in% exclude_cells])
# # saveRDS(cellchat, file="integration/integration_1/integration_1_louvain/RDS_Dir/integration_1_louvain_annot_cellchat_de_subset.rds")
# 
# 
# # cellchat <- readRDS('integration/integration_1/integration_1_louvain/RDS_Dir/integration_1_louvain_annot_cellchat_de_subset.rds')
# 
# # inference of cell-cell communication network
# # cellchat <- computeCommunProb(cellchat, type = "triMean")
# # cellchat <- filterCommunication(cellchat, min.cells = 10)
# # cellchat <- computeCommunProbPathway(cellchat)
# # cellchat <- aggregateNet(cellchat)
# 
# # compute network centrality scores
# # cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
# # cellchat <- computeNetSimilarity(cellchat, type = "functional")
# # cellchat <- netEmbedding(cellchat, type = "functional",umap.method = "uwot")
# # cellchat <- netClustering(cellchat, type = "functional")
# 
# saveRDS(cellchat, file="integration/integration_1/integration_1_louvain/RDS_Dir/integration_1_louvain_annot_cellchat_de_subset.rds")
# 
#### 3.0 additional functions on cellchat before making plots  ####
cellchat <- readRDS("integration_2/integration_2_leiden/RDS_Dir/integration_2_leiden_annot_cellchat_de.rds")


# compute network centrality scores
# cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
# cellchat <- computeNetSimilarity(cellchat, type = "functional")
# cellchat <- netEmbedding(cellchat, type = "functional",umap.method = "uwot")
# cellchat <- netClustering(cellchat, type = "functional")

# cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "net")


#### 4.0 plot set up ####
# follow the order of levels

# name the colors based on the cell types


cell_types <- c(
                "B", 
                "CD4 T", 
                "CD8 T", 
                "DC", 
                "EC", 
                "Mono/Mph", 
                "NK")


whole_color_use <- c(
  colorRampPalette(colors = rev(c("#faace0", "#771155")))(length(grep(cell_types[1], unique(cellchat@idents), value = T))),
  colorRampPalette(colors = rev(c("#c2b5a7", "#603808")))(length(grep(cell_types[3], unique(cellchat@idents), value = T))), 
  colorRampPalette(colors = rev(c("#fed9b7", "#f07167")))(length(grep(cell_types[4], unique(cellchat@idents), value = T))), 
  colorRampPalette(colors = rev(c("#a9fccf", "#194a2f")))(length(grep(cell_types[4], unique(cellchat@idents), value = T))), 
  colorRampPalette(colors = rev(c("#d7e0ab", "#679436")))(length(grep(cell_types[5], unique(cellchat@idents), value = T))), 
  colorRampPalette(colors = rev(c("#d4e3ff", "#43639c")))(length(grep(cell_types[6], unique(cellchat@idents), value = T))),
  colorRampPalette(colors = rev(c("#b2f7ef", "#07beb8")))(length(grep(cell_types[7], unique(cellchat@idents), value = T)))
)



# name the colors based on the cell types
names(whole_color_use) <- unlist(lapply(cell_types, function(ct) grep(ct, levels(cellchat@idents), value = T)))




groupSize <- as.numeric(table(cellchat@idents))



#### plot 1 interaction weights/strength ####
par(mfrow = c(1,1), xpd=TRUE)

pdf(paste0("./integration_2/integration_2_leiden/cellchat_de/netVisual_circle_all.pdf"), 
    width = 10, height = 6)
netVisual_circle(cellchat@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 vertex.size.max = 5,
                 title.name = "Interaction weights/strength",
                 arrow.size = 0.05, 
                 color.use = whole_color_use, 
                 vertex.label.cex = 0.3)
dev.off()


mat <- cellchat@net$weight

pdf(paste0("./integration_2/integration_2_leiden/cellchat_de/netVisual_circle_per_celltype.pdf"), 
    width = 10, height = 6)
par(mfrow = c(2,4), xpd=TRUE)


for (j in (0:6)){
  if (j==0){
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    tmp = mat2
    tmp[c(1:7), ] <- mat[c(1:7), ]
    netVisual_circle(tmp, 
                     vertex.weight = groupSize, 
                     weight.scale = T, 
                     edge.weight.max = max(mat), 
                     title.name = cell_types[1],
                     color.use = whole_color_use,
                     vertex.label.cex = 0.2)
  }else{
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[c((1+j*7):(7+j*7)), ] <- mat[c((1+j*7):(7+j*7)), ]
    netVisual_circle(mat2, 
                     vertex.weight = groupSize,
                     weight.scale = T, 
                     edge.weight.max = max(mat),
                     title.name = cell_types[j+1],
                     color.use = whole_color_use,
                     vertex.label.cex = 0.2)
  }
}
dev.off()



for (j in (0:6)) {
  cell_type <- gsub("/", "_", cell_types[j+1])
  pdf(paste0("./integration_2/integration_2_leiden/cellchat_de/netVisual_circle_",
             cell_type, ".pdf"), width = 10, height = 6)
  par(mfrow = c(2,4), xpd=TRUE)
  for (i in rev((1+j*7):(7+j*7))) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, 
                     vertex.weight = groupSize, 
                     weight.scale = T, 
                     edge.weight.max = max(mat), 
                     title.name = rownames(mat)[i],
                     color.use = whole_color_use,
                     vertex.label.cex = 0.2)
  }
  dev.off()
}













#### plot 2 incoming/outgoing strength ####
# View(netAnalysis_signalingRole_scatter)
object <- cellchat
centr <- slot(object, "netP")$centr
outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
dimnames(outgoing) <- list(levels(object@idents), names(centr))
dimnames(incoming) <- dimnames(outgoing)
x.measure = "outdeg"
y.measure = "indeg"

for (i in 1:length(centr)) {
  outgoing[, i] <- centr[[i]][[x.measure]]
  incoming[, i] <- centr[[i]][[y.measure]]
}

signaling = NULL
if (is.null(signaling)) {
  message("Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways")
} else {
  message("Signaling role analysis on the cell-cell communication network from user's input")
  signaling <- signaling[signaling %in% object@netP$pathways]
  if (length(signaling) == 0) {
    stop("There is no significant communication for the input signaling. All the significant signaling are shown in `object@netP$pathways`")
  }
  outgoing <- outgoing[, signaling, drop = FALSE]
  incoming <- incoming[, signaling, drop = FALSE]
}


outgoing.cells <- rowSums(outgoing)
incoming.cells <- rowSums(incoming)
num.link <- aggregateNet(object, signaling = signaling, 
                         return.object = FALSE, remove.isolate = FALSE)$count
num.link <- rowSums(num.link) + colSums(num.link) - diag(num.link)
df <- data.frame(x = outgoing.cells, y = incoming.cells, 
                 labels = names(incoming.cells), Count = num.link)
df$labels <- factor(df$labels, levels = names(incoming.cells))

df <- df %>% separate(labels, into = c("celltype", "stimili.time"), sep = "\\|", remove = FALSE)

pdf(paste0("./integration_2/integration_2_leiden/cellchat_de/netAnalysis_signalingRole_scatter_per_celltype.pdf"), 
    width = 10, height = 5.5)
ggplot(data = df, aes(x, y)) + 
  geom_point(aes(size = Count, colour = labels, fill = labels)) +
  CellChat_theme_opts() +
  theme(text = element_text(size = 10), legend.key.height = grid::unit(0.15, "in")) + 
  labs(title = NULL, x = "Outgoing interaction strength", y= "Incoming interaction strength") + 
  theme(plot.title = element_text(size = 10, face = "plain")) + 
  theme(axis.line.x = element_line(linewidth = 0.25),
        axis.line.y = element_line(linewidth = 0.25)) +
  scale_fill_manual(values = ggplot2::alpha(whole_color_use, alpha = 0.6), drop = FALSE) + 
  guides(fill = "none") +
  scale_size_continuous(range = c(2, 6)) +
  scale_colour_manual(values = whole_color_use, drop = FALSE) + 
  guides(colour = "none") +
  ggrepel::geom_text_repel(mapping = aes(label = stimili.time, colour = labels),
                           size = 2, show.legend = F, 
                           segment.size = 0.2, segment.alpha = 0.5) +
  theme(legend.position = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        legend.position = "right",
        strip.background = element_rect(fill = "white", color = "black", size = 1),
        axis.text.y = element_text(size = 9),
        strip.text.y = element_blank(),
        strip.placement = "inside", 
        strip.background.x = element_blank()) +
  xlim(c(0,45)) +
  ylim(c(0,45)) +
  geom_abline(slope = 1, color = "lightgrey", linetype = "dashed") +
  facet_wrap(. ~ celltype, scales = "free", ncol = 4, nrow = 2)
dev.off()




#### plot 3 contribution of signal in incoming/outgoing strength per celltype ####
#### netP ####
# run the line below
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

pdf(paste0("./integration_2/integration_2_leiden/cellchat_de/netAnalysis_signalingRole_heatmap_all_netP.pdf"), 
    width = 7, height = 8)




netAnalysis_signalingRole_heatmap_top_reorder <- function (object, signaling = NULL, 
                                                           pattern = c("outgoing", "incoming", "all"), 
                                                           slot.name = "netP", color.use = NULL, 
                                                           color.heatmap = "BuGn", 
                                                           title = NULL, 
                                                           width = 10, 
                                                           height = 8, 
                                                           # new feature, loop through each column and keep the top n ligand names with highest strength
                                                           top = NULL, 
                                                           cellchat_ident_level = NULL, 
                                                           font.size = 6, font.size.title = 10, 
                                                           cluster.rows = FALSE, 
                                                           cluster.cols = FALSE) {
  pattern <- match.arg(pattern)
  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
  }
  centr <- slot(object, slot.name)$centr
  outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(object@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  
  for (i in 1:length(centr)) {
    outgoing[, i] <- centr[[i]]$outdeg
    incoming[, i] <- centr[[i]]$indeg
  }
  if (pattern == "outgoing") {
    mat <- t(outgoing)
    if (!is.null(cellchat_ident_level)){
      mat <- mat[, cellchat_ident_level]
    }
    if (!is.null(top)){
      top_values <- lapply(1:ncol(mat), function(col_idx) {
        sort(mat[, col_idx], decreasing = TRUE)[1:top]
      })
      rows_with_top_values <- unique(unlist(lapply(1:ncol(mat), function(col_idx) {
        top_values <- top_values[[col_idx]]
        which(mat[, col_idx] %in% top_values)
      })))
      mat <- mat[rows_with_top_values, ]
    }
    legend.name <- "Outgoing"
  }
  else if (pattern == "incoming") {
    mat <- t(incoming)
    if (!is.null(cellchat_ident_level)){
      mat <- mat[, cellchat_ident_level]
    }
    if (!is.null(top)){
      top_values <- lapply(1:ncol(mat), function(col_idx) {
        sort(mat[, col_idx], decreasing = TRUE)[1:top]
      })
      rows_with_top_values <- unique(unlist(lapply(1:ncol(mat), function(col_idx) {
        top_values <- top_values[[col_idx]]
        which(mat[, col_idx] %in% top_values)
      })))
      mat <- mat[rows_with_top_values, ]
    }
    legend.name <- "Incoming"
  }
  else if (pattern == "all") {
    mat <- t(outgoing + incoming)
    if (!is.null(cellchat_ident_level)){
      mat <- mat[, cellchat_ident_level]
    }
    if (!is.null(top)){
      top_values <- lapply(1:ncol(mat), function(col_idx) {
        sort(mat[, col_idx], decreasing = TRUE)[1:top]
      })
      rows_with_top_values <- unique(unlist(lapply(1:ncol(mat), function(col_idx) {
        top_values <- top_values[[col_idx]]
        which(mat[, col_idx] %in% top_values)
      })))
      mat <- mat[rows_with_top_values, ]
    }
    legend.name <- "Overall"
  }
  if (is.null(title)) {
    title <- paste0(legend.name, " signaling patterns")
  }
  else {
    title <- paste0(paste0(legend.name, " signaling patterns"), 
                    " - ", title)
  }
  if (!is.null(signaling)) {
    mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
    mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
    idx <- match(rownames(mat1), signaling)
    mat[idx[!is.na(idx)], ] <- mat1
    dimnames(mat) <- list(signaling, colnames(mat1))
  }
  mat.ori <- mat
  mat <- sweep(mat, 1L, apply(mat, 1, max), "/", check.margin = FALSE)
  mat[mat == 0] <- NA
  if (is.null(color.use)) {
    color.use <- scPalette(length(colnames(mat)))
  }
  color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, 
                                                                             name = color.heatmap))))(100)
  df <- data.frame(group = colnames(mat)) 
  rownames(df) <- colnames(mat) 
  names(color.use) <- colnames(mat) 
  
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
                                      which = "column", show_legend = FALSE, show_annotation_name = FALSE, 
                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(mat.ori), 
                                                  border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
                          show_annotation_name = FALSE)
  pSum <- rowSums(mat.ori)
  pSum.original <- pSum
  pSum <- -1/log(pSum)
  pSum[is.na(pSum)] <- 0
  idx1 <- which(is.infinite(pSum) | pSum < 0)
  if (length(idx1) > 0) {
    values.assign <- seq(max(pSum) * 1.1, max(pSum) * 1.5, 
                         length.out = length(idx1))
    position <- sort(pSum.original[idx1], index.return = TRUE)$ix
    pSum[idx1] <- values.assign[match(1:length(idx1), position)]
  }
  ha1 = rowAnnotation(Strength = anno_barplot(pSum, border = FALSE), 
                      show_annotation_name = FALSE)
  if (min(mat, na.rm = T) == max(mat, na.rm = T)) {
    legend.break <- max(mat, na.rm = T)
  }
  else {
    legend.break <- c(round(min(mat, na.rm = T), digits = 1), 
                      round(max(mat, na.rm = T), digits = 1))
  }
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", 
                name = "Relative strength", bottom_annotation = col_annotation, 
                top_annotation = ha2, right_annotation = ha1, cluster_rows = cluster.rows, 
                cluster_columns = cluster.rows, row_names_side = "left", 
                row_names_rot = 0, row_names_gp = gpar(fontsize = font.size), 
                column_names_gp = gpar(fontsize = font.size), 
                width = unit(width, "cm"), height = unit(height, "cm"), 
                column_title = title, 
                column_title_gp = gpar(fontsize = font.size.title), 
                column_names_rot = 90, 
                heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"), 
                                            title_position = "leftcenter-rot", 
                                            border = NA, 
                                            at = legend.break, 
                                            legend_height = unit(20, "mm"),
                                            labels_gp = gpar(fontsize = 8), grid_width = unit(2,"mm")))
  return(list(plot=ht1, mat=mat, df=df, mat.ori=mat.ori))
}


mat <- netAnalysis_signalingRole_heatmap_top_reorder(cellchat, slot.name = "netP",
                                                     font.size = 7,
                                                     top = 5, # top 5 pathwyas per condition 
                                                     height = 20,
                                                     width = 12, 
                                                     pattern = "all",
                                                     color.heatmap = "YlGnBu",
                                                     color.use = whole_color_use)$mat




# df <- netAnalysis_signalingRole_heatmap_reorder(cellchat, slot.name = "netP",
#                                                  font.size = 7,
#                                                 top = 50,
#                                                  height = 20,
#                                                  width = 12,
#                                                  pattern = "all",
#                                                  color.heatmap = "YlGnBu",
#                                                  color.use = whole_color_use2)$df

mat.ori <- netAnalysis_signalingRole_heatmap_top_reorder(cellchat, slot.name = "netP",
                                                         font.size = 7,
                                                         top = 5,
                                                         height = 20,
                                                         width = 12, 
                                                         pattern = "all",
                                                         color.heatmap = "YlGnBu",
                                                         color.use = whole_color_use)$mat.ori


cell_types_reset <- c("Mono/Mph",
                      "B", 
                      "CD4 T", 
                      "CD8 T", 
                      "NK",
                      "EC", 
                      "DC")


cell_types_reset_level <- c(rev(grep(cell_types_reset[1], levels(cellchat@idents), value = T)),
                            rev(grep(cell_types_reset[2], levels(cellchat@idents), value = T)),
                            rev(grep(cell_types_reset[3], levels(cellchat@idents), value = T)),
                            rev(grep(cell_types_reset[4], levels(cellchat@idents), value = T)),
                            rev(grep(cell_types_reset[5], levels(cellchat@idents), value = T)),
                            rev(grep(cell_types_reset[6], levels(cellchat@idents), value = T)),
                            rev(grep(cell_types_reset[7], levels(cellchat@idents), value = T)))

whole_color_use_reset <- c(
  colorRampPalette(colors = c("#d4e3ff", "#43639c"))(length(grep(cell_types_reset[1], unique(cellchat@idents), value = T))),
  colorRampPalette(colors = c("#faace0", "#771155"))(length(grep(cell_types_reset[2], unique(cellchat@idents), value = T))), 
  colorRampPalette(colors = c("#c2b5a7", "#603808"))(length(grep(cell_types_reset[3], unique(cellchat@idents), value = T))), 
  colorRampPalette(colors = c("#fed9b7", "#f07167"))(length(grep(cell_types_reset[4], unique(cellchat@idents), value = T))), 
  colorRampPalette(colors = c("#b2f7ef", "#07beb8"))(length(grep(cell_types_reset[5], unique(cellchat@idents), value = T))), 
  colorRampPalette(colors = c("#d7e0ab", "#679436"))(length(grep(cell_types_reset[6], unique(cellchat@idents), value = T))),
  colorRampPalette(colors = c("#a9fccf", "#194a2f"))(length(grep(cell_types_reset[7], unique(cellchat@idents), value = T)))
)

color.heatmap <- 'YlGnBu'
color.use <- whole_color_use_reset

color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, 
                                                                           name = color.heatmap))))(100)
font.size = 7
height = 10
width = 10
font.size.title = 10
cluster.rows = FALSE
cluster.cols = FALSE
title = NULL
mat.ori <- mat.ori[, cell_types_reset_level]
mat <- mat[, cell_types_reset_level]

pSum <- rowSums(mat.ori)
pSum.original <- pSum
pSum <- -1/log(pSum)
pSum[is.na(pSum)] <- 0
idx1 <- which(is.infinite(pSum) | pSum < 0)
if (length(idx1) > 0) {
  values.assign <- seq(max(pSum) * 1.1, max(pSum) * 1.5, 
                       length.out = length(idx1))
  position <- sort(pSum.original[idx1], index.return = TRUE)$ix
  pSum[idx1] <- values.assign[match(1:length(idx1), position)]
}
pSum <- sort(pSum, decreasing = TRUE)

mat.ori <- mat.ori[names(pSum), ]
mat <- mat[names(pSum), ]


df <- data.frame(group = colnames(mat))
rownames(df) <- colnames(mat)
names(color.use) <- colnames(mat)

col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
                                    which = "column", show_legend = FALSE, show_annotation_name = FALSE, 
                                    simple_anno_size = grid::unit(0.2, "cm"))
ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(mat.ori), 
                                                border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
                        show_annotation_name = FALSE)


# barplot for row annotation
ha1 = rowAnnotation(Strength = anno_barplot(pSum, border = FALSE, gp = gpar(fill = "white")), 
                    show_annotation_name = FALSE)

# ha1 = rowAnnotation(Strength = anno_points(pSum, border = FALSE, gp = gpar(color = "grey")), 
#                     show_annotation_name = FALSE)





if (min(mat, na.rm = T) == max(mat, na.rm = T)) {
  legend.break <- max(mat, na.rm = T)
}else {
  legend.break <- c(round(min(mat, na.rm = T), digits = 1), 
                    round(max(mat, na.rm = T), digits = 1))
}
Heatmap(mat, col = color.heatmap.use, na_col = "white", 
        name = "Relative strength", bottom_annotation = col_annotation, 
        top_annotation = ha2, right_annotation = ha1, 
        cluster_rows = cluster.rows, 
        cluster_columns = cluster.rows, 
        row_names_side = "left", 
        row_names_rot = 0, row_names_gp = gpar(fontsize = font.size), 
        column_names_gp = gpar(fontsize = font.size), 
        # width = unit(width, "cm"), height = unit(height, "cm"), 
        column_title = title, 
        column_title_gp = gpar(fontsize = font.size.title), 
        column_names_rot = 90, 
        heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"), 
                                    title_position = "leftcenter-rot", 
                                    border = NA, 
                                    at = legend.break, 
                                    legend_height = unit(20, "mm"),
                                    labels_gp = gpar(fontsize = 8), grid_width = unit(2,"mm")))



dev.off()



h3 <- netAnalysis_signalingRole_heatmap(cellchat, slot.name = "net",
                                        font.size = 7,
                                        height = 20,
                                        width = 12,
                                        pattern = "all",
                                        color.heatmap = "YlGnBu",
                                        color.use = whole_color_use)
h3



#### net ####
# run the line below
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "net")

pdf(paste0("./integration_2/integration_2_leiden/cellchat_de/netAnalysis_signalingRole_heatmap_all_net.pdf"), 
    width = 7, height = 10)




netAnalysis_signalingRole_heatmap_top_reorder <- function (object, signaling = NULL, 
                                                       pattern = c("outgoing", "incoming", "all"), 
                                                       slot.name = "netP", color.use = NULL, 
                                                       color.heatmap = "BuGn", 
                                                       title = NULL, 
                                                       width = 10, 
                                                       height = 8, 
                                                       # new feature, loop through each column and keep the top n ligand names with highest strength
                                                       top = NULL, 
                                                       cellchat_ident_level = NULL, 
                                                       font.size = 6, font.size.title = 10, 
                                                       cluster.rows = FALSE, 
                                                       cluster.cols = FALSE) {
  pattern <- match.arg(pattern)
  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
  }
  centr <- slot(object, slot.name)$centr
  outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(object@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  
  for (i in 1:length(centr)) {
    outgoing[, i] <- centr[[i]]$outdeg
    incoming[, i] <- centr[[i]]$indeg
  }
  if (pattern == "outgoing") {
    mat <- t(outgoing)
    if (!is.null(cellchat_ident_level)){
      mat <- mat[, cellchat_ident_level]
    }
    if (!is.null(top)){
      top_values <- lapply(1:ncol(mat), function(col_idx) {
        sort(mat[, col_idx], decreasing = TRUE)[1:top]
      })
      rows_with_top_values <- unique(unlist(lapply(1:ncol(mat), function(col_idx) {
        top_values <- top_values[[col_idx]]
        which(mat[, col_idx] %in% top_values)
      })))
      mat <- mat[rows_with_top_values, ]
    }
    legend.name <- "Outgoing"
  }
  else if (pattern == "incoming") {
    mat <- t(incoming)
    if (!is.null(cellchat_ident_level)){
      mat <- mat[, cellchat_ident_level]
    }
    if (!is.null(top)){
      top_values <- lapply(1:ncol(mat), function(col_idx) {
        sort(mat[, col_idx], decreasing = TRUE)[1:top]
      })
      rows_with_top_values <- unique(unlist(lapply(1:ncol(mat), function(col_idx) {
        top_values <- top_values[[col_idx]]
        which(mat[, col_idx] %in% top_values)
      })))
      mat <- mat[rows_with_top_values, ]
    }
    legend.name <- "Incoming"
  }
  else if (pattern == "all") {
    mat <- t(outgoing + incoming)
    if (!is.null(cellchat_ident_level)){
      mat <- mat[, cellchat_ident_level]
    }
    if (!is.null(top)){
      top_values <- lapply(1:ncol(mat), function(col_idx) {
        sort(mat[, col_idx], decreasing = TRUE)[1:top]
      })
      rows_with_top_values <- unique(unlist(lapply(1:ncol(mat), function(col_idx) {
        top_values <- top_values[[col_idx]]
        which(mat[, col_idx] %in% top_values)
      })))
      mat <- mat[rows_with_top_values, ]
    }
    legend.name <- "Overall"
  }
  if (is.null(title)) {
    title <- paste0(legend.name, " signaling patterns")
  }
  else {
    title <- paste0(paste0(legend.name, " signaling patterns"), 
                    " - ", title)
  }
  if (!is.null(signaling)) {
    mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
    mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
    idx <- match(rownames(mat1), signaling)
    mat[idx[!is.na(idx)], ] <- mat1
    dimnames(mat) <- list(signaling, colnames(mat1))
  }
  mat.ori <- mat
  mat <- sweep(mat, 1L, apply(mat, 1, max), "/", check.margin = FALSE)
  mat[mat == 0] <- NA
  if (is.null(color.use)) {
    color.use <- scPalette(length(colnames(mat)))
  }
  color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, 
                                                                             name = color.heatmap))))(100)
  df <- data.frame(group = colnames(mat)) 
  rownames(df) <- colnames(mat) 
  names(color.use) <- colnames(mat) 
  
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
                                      which = "column", show_legend = FALSE, show_annotation_name = FALSE, 
                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(mat.ori), 
                                                  border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
                          show_annotation_name = FALSE)
  pSum <- rowSums(mat.ori)
  pSum.original <- pSum
  pSum <- -1/log(pSum)
  pSum[is.na(pSum)] <- 0
  idx1 <- which(is.infinite(pSum) | pSum < 0)
  if (length(idx1) > 0) {
    values.assign <- seq(max(pSum) * 1.1, max(pSum) * 1.5, 
                         length.out = length(idx1))
    position <- sort(pSum.original[idx1], index.return = TRUE)$ix
    pSum[idx1] <- values.assign[match(1:length(idx1), position)]
  }
  ha1 = rowAnnotation(Strength = anno_barplot(pSum, border = FALSE), 
                      show_annotation_name = FALSE)
  if (min(mat, na.rm = T) == max(mat, na.rm = T)) {
    legend.break <- max(mat, na.rm = T)
  }
  else {
    legend.break <- c(round(min(mat, na.rm = T), digits = 1), 
                      round(max(mat, na.rm = T), digits = 1))
  }
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", 
                name = "Relative strength", bottom_annotation = col_annotation, 
                top_annotation = ha2, right_annotation = ha1, cluster_rows = cluster.rows, 
                cluster_columns = cluster.rows, row_names_side = "left", 
                row_names_rot = 0, row_names_gp = gpar(fontsize = font.size), 
                column_names_gp = gpar(fontsize = font.size), 
                width = unit(width, "cm"), height = unit(height, "cm"), 
                column_title = title, 
                column_title_gp = gpar(fontsize = font.size.title), 
                column_names_rot = 90, 
                heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"), 
                                            title_position = "leftcenter-rot", 
                                            border = NA, 
                                            at = legend.break, 
                                            legend_height = unit(20, "mm"),
                                            labels_gp = gpar(fontsize = 8), grid_width = unit(2,"mm")))
  return(list(plot=ht1, mat=mat, df=df, mat.ori=mat.ori))
}


mat <- netAnalysis_signalingRole_heatmap_top_reorder(cellchat, slot.name = "net",
                                                 font.size = 7,
                                                 top = 5,
                                                 height = 20,
                                                 width = 12, 
                                                 pattern = "all",
                                                 color.heatmap = "YlGnBu",
                                                 color.use = whole_color_use)$mat

# df <- netAnalysis_signalingRole_heatmap_reorder(cellchat, slot.name = "net",
#                                                  font.size = 7,
#                                                 top = 50,
#                                                  height = 20,
#                                                  width = 12,
#                                                  pattern = "all",
#                                                  color.heatmap = "YlGnBu",
#                                                  color.use = whole_color_use2)$df

mat.ori <- netAnalysis_signalingRole_heatmap_top_reorder(cellchat, slot.name = "net",
                                                     font.size = 7,
                                                     top = 5,
                                                     height = 20,
                                                     width = 12, 
                                                     pattern = "all",
                                                     color.heatmap = "YlGnBu",
                                                     color.use = whole_color_use)$mat.ori



cell_types_reset <- c("Mono/Mph",
                      "B", 
                      "CD4 T", 
                      "CD8 T", 
                      "NK",
                      "EC", 
                      "DC")


cell_types_reset_level <- c(rev(grep(cell_types_reset[1], levels(cellchat@idents), value = T)),
                            rev(grep(cell_types_reset[2], levels(cellchat@idents), value = T)),
                            rev(grep(cell_types_reset[3], levels(cellchat@idents), value = T)),
                            rev(grep(cell_types_reset[4], levels(cellchat@idents), value = T)),
                            rev(grep(cell_types_reset[5], levels(cellchat@idents), value = T)),
                            rev(grep(cell_types_reset[6], levels(cellchat@idents), value = T)),
                            rev(grep(cell_types_reset[7], levels(cellchat@idents), value = T))
                            )

whole_color_use_reset <- c(
  colorRampPalette(colors = c("#d4e3ff", "#43639c"))(length(grep(cell_types_reset[1], unique(cellchat@idents), value = T))),
  colorRampPalette(colors = c("#faace0", "#771155"))(length(grep(cell_types_reset[2], unique(cellchat@idents), value = T))), 
  colorRampPalette(colors = c("#c2b5a7", "#603808"))(length(grep(cell_types_reset[3], unique(cellchat@idents), value = T))), 
  colorRampPalette(colors = c("#fed9b7", "#f07167"))(length(grep(cell_types_reset[4], unique(cellchat@idents), value = T))), 
  colorRampPalette(colors = c("#b2f7ef", "#07beb8"))(length(grep(cell_types_reset[5], unique(cellchat@idents), value = T))), 
  colorRampPalette(colors = c("#d7e0ab", "#679436"))(length(grep(cell_types_reset[6], unique(cellchat@idents), value = T))),
  colorRampPalette(colors = c("#a9fccf", "#194a2f"))(length(grep(cell_types_reset[7], unique(cellchat@idents), value = T)))
)

color.heatmap <- 'YlGnBu'
color.use <- whole_color_use_reset


color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, 
                                                                           name = color.heatmap))))(100)
font.size = 6
height = 10
width = 10
font.size.title = 10
cluster.rows = FALSE
cluster.cols = FALSE
title = NULL
mat.ori <- mat.ori[, cell_types_reset_level]
mat <- mat[, cell_types_reset_level]

pSum <- rowSums(mat.ori)
pSum.original <- pSum
pSum <- -1/log(pSum)
pSum[is.na(pSum)] <- 0
idx1 <- which(is.infinite(pSum) | pSum < 0)
if (length(idx1) > 0) {
  values.assign <- seq(max(pSum) * 1.1, max(pSum) * 1.5, 
                       length.out = length(idx1))
  position <- sort(pSum.original[idx1], index.return = TRUE)$ix
  pSum[idx1] <- values.assign[match(1:length(idx1), position)]
}
pSum <- sort(pSum, decreasing = TRUE)

# cd39_index <- which(names(pSum) == "CD39")
# 
# # Subset pSum to keep only up to "CD39"
# pSum_filtered <- pSum[1:cd39_index]


mat.ori <- mat.ori[names(pSum), ]
mat <- mat[names(pSum), ]


df <- data.frame(group = colnames(mat))
rownames(df) <- colnames(mat)
names(color.use) <- colnames(mat)

col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
                                    which = "column", show_legend = FALSE, show_annotation_name = FALSE, 
                                    simple_anno_size = grid::unit(0.2, "cm"))
ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(mat.ori), 
                                                border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
                        show_annotation_name = FALSE)


# barplot for row annotation
ha1 = rowAnnotation(Strength = anno_barplot(pSum, border = FALSE, gp = gpar(fill = "white")), 
                    show_annotation_name = FALSE)

# ha1 = rowAnnotation(Strength = anno_points(pSum, border = FALSE, gp = gpar(color = "grey")), 
#                     show_annotation_name = FALSE)





if (min(mat, na.rm = T) == max(mat, na.rm = T)) {
  legend.break <- max(mat, na.rm = T)
}else {
  legend.break <- c(round(min(mat, na.rm = T), digits = 1), 
                    round(max(mat, na.rm = T), digits = 1))
}
Heatmap(mat, col = color.heatmap.use, na_col = "white", 
        name = "Relative strength", bottom_annotation = col_annotation, 
        top_annotation = ha2, right_annotation = ha1, 
        cluster_rows = cluster.rows, 
        cluster_columns = cluster.rows, 
        row_names_side = "left", 
        row_names_rot = 0, row_names_gp = gpar(fontsize = font.size), 
        column_names_gp = gpar(fontsize = font.size), 
        # width = unit(width, "cm"), height = unit(height, "cm"), 
        column_title = title, 
        column_title_gp = gpar(fontsize = font.size.title), 
        column_names_rot = 90, 
        heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"), 
                                    title_position = "leftcenter-rot", 
                                    border = NA, 
                                    at = legend.break, 
                                    legend_height = unit(20, "mm"),
                                    labels_gp = gpar(fontsize = 8), grid_width = unit(2,"mm")))



dev.off()



h3 <- netAnalysis_signalingRole_heatmap(cellchat, slot.name = "net",
                                        font.size = 7,
                                        height = 20,
                                        width = 12,
                                        pattern = "all",
                                        color.heatmap = "YlGnBu",
                                        color.use = whole_color_use)
h3






##--------------------------------------------
##--------------------------------------------
##--------------------------------------------
##--------------------------------------------
##--------------------------------------------
##--------------------------------------------



# 5.0 new figures  -----------------------------------------------------------------
cellchat <- readRDS("integration_2/integration_2_leiden/RDS_Dir/integration_1_louvain_annot_cellchat_de_subset.rds")

cell_types_reset <- c("B", 
                      "CD4 T", 
                      "CD8 T", 
                      "NK",
                      "Mono/ Mph",
                      "DC",
                      "EC")

cell_types_reset_abbr <- c("B", 
                           "CD4 T", 
                           "CD8 T", 
                           "NK", 
                           "Mono/Mph", 
                           "DC", 
                           "EC")

cell_types_reset_level <- c(rev(grep(cell_types_reset[1], levels(cellchat@idents), value = T)),
                            rev(grep(cell_types_reset[2], levels(cellchat@idents), value = T)),
                            rev(grep(cell_types_reset[3], levels(cellchat@idents), value = T)),
                            rev(grep(cell_types_reset[4], levels(cellchat@idents), value = T)),
                            rev(grep(cell_types_reset[5], levels(cellchat@idents), value = T)),
                            rev(grep(cell_types_reset[6], levels(cellchat@idents), value = T)),
                            rev(grep(cell_types_reset[7], levels(cellchat@idents), value = T)),
                            rev(grep(cell_types_reset[8], levels(cellchat@idents), value = T)))

cell_types_reset_abbr_level <- cell_types_reset_level
cell_types_reset_abbr_level <- gsub(" cells", "", cell_types_reset_abbr_level)
cell_types_reset_abbr_level <- gsub("Monocyte-derived Mph", "Mono Mph", cell_types_reset_abbr_level)
cell_types_reset_abbr_level <- gsub("Alveolar Mph MT-positive", "Alv Mph", cell_types_reset_abbr_level)
cell_types_reset_abbr_level <- gsub("Migratory DC", "DC", cell_types_reset_abbr_level)
cell_types_reset_abbr_level <- gsub("EC general capillary", "EC", cell_types_reset_abbr_level)
cell_types_reset_abbr_level <- gsub("unstim", "Unstim", cell_types_reset_abbr_level)
cell_types_reset_abbr_level <- gsub("\\|", " | ", cell_types_reset_abbr_level)


## a. number of interactions ----

cell_types_reset_level_cnts <- factor(
  c(rep("B", 7), 
    rep("CD4 T", 7),
    rep("CD8 T", 7),
    rep("NK", 7),
    rep("Mono Mph", 7),
    rep("DC", 7),
    rep("EC", 7)),
  levels = c("B", "CD4 T", "CD8 T", "NK", "Mono Mph", "DC", "EC")
)



# cell_type_colors <- c("#f77f00", "#075a84", "#9163cb", "#ff4d6d", "#679436", "#f07167", "#603808", "#07beb8")

cell_type_colors <- c(
  "B" = "#d4c2fc",     
  "CD4 T" = "#ff758f",   
  "CD8 T" = "#98c9a3",   
  "NK" = "#ffb4a2",  
  "Mono Mph" = "#6d98ba", 
  "DC" = "#b2f7ef",      
  "EC" = "#ddb892" 
)


mat <- cellchat@net$count[cell_types_reset_level, cell_types_reset_level]

for (i in 1:length(rownames(mat))) {
  # rownames
  rownames(mat)[i] <- gsub(" cells", "", rownames(mat)[i])
  rownames(mat)[i] <- gsub("Monocyte-derived Mph", "Mono Mph", rownames(mat)[i])
  rownames(mat)[i] <- gsub("Alveolar Mph MT-positive", "Alv Mph", rownames(mat)[i])
  rownames(mat)[i] <- gsub("Migratory DC", "DC", rownames(mat)[i])
  rownames(mat)[i] <- gsub("EC general capillary", "EC", rownames(mat)[i])
  rownames(mat)[i] <- gsub("unstim", "Unstim", rownames(mat)[i])
  # colnames
  colnames(mat)[i] <- gsub(" cells", "", colnames(mat)[i])
  colnames(mat)[i] <- gsub("Monocyte-derived Mph", "Mono Mph", colnames(mat)[i])
  colnames(mat)[i] <- gsub("Alveolar Mph MT-positive", "Alv Mph", colnames(mat)[i])
  colnames(mat)[i] <- gsub("Migratory DC", "DC", colnames(mat)[i])
  colnames(mat)[i] <- gsub("EC general capillary", "EC", colnames(mat)[i])
  colnames(mat)[i] <- gsub("unstim", "Unstim", colnames(mat)[i])
}

total_target <- colSums(mat)
total_source <- rowSums(mat)
bar_ylim <- c(0, max(total_source, total_target))

# # export source data
# mat %>% 
#   as_tibble(rownames = "source") %>% 
#   save_table("source_data/figure_3a", "Figure 3a")

# could set up a function
# plot_n_interactions similar to the format below 
# https://github.com/csbg/neuroblastoma/blob/main/plot_figure_3_S3.R
BASE_TEXT_SIZE_PT = 5.5
TITLE_TEXT_SIZE_PT = 11

pdf(paste0("./integration/integration_1/integration_1_louvain/cellchat_de/plot_n_interactions.pdf"),
    width = 7, height = 6.6)
    
# make plot
Heatmap(
  mat,
  col = RColorBrewer::brewer.pal(9, "PuBuGn"),
  name = "Number of\ninteractions",
  
  cluster_rows = FALSE, 
  row_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
  row_names_side = "left",
  row_title = "Source (ligand)", 
  row_title_gp = gpar(fontsize = TITLE_TEXT_SIZE_PT), 
  
  cluster_columns = FALSE, 
  column_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
  column_title = "Target (receptor)",
  column_title_side = "bottom",
  column_title_gp = gpar(fontsize = TITLE_TEXT_SIZE_PT), 
  
  # width = unit(20, "mm"),
  # height = unit(20, "mm"),
  
  top_annotation = HeatmapAnnotation(
    count_bar = anno_barplot(
      total_target,
      ylim = bar_ylim,
      border = FALSE, 
      axis = TRUE, # set FALSE if y-axis is not needed
      gp = gpar(fill = "gray70", col = "gray70")
    ),
    count_text = anno_text(
      total_target, 
      gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    ), 
    simple_anno_size_adjust = TRUE,
    show_annotation_name = FALSE
  ), 
  
  right_annotation = rowAnnotation(
    count_text = anno_text(
      total_source, 
      gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    ), 
    count_bar = anno_barplot(
      total_source,
      ylim = bar_ylim,
      border = FALSE, 
      axis = TRUE, # set FALSE if y-axis is not needed
      gp = gpar(fill = "gray70", col = "gray70")
    ),
    simple_anno_size_adjust = TRUE,
    show_annotation_name = FALSE
  ), 
  
  left_annotation = rowAnnotation(
    cell_type = anno_block(
      gp = gpar(fill = cell_type_colors),
      labels = names(cell_type_colors),
      labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      width = unit(4, "mm")
      )
  ),
  
  
  bottom_annotation = columnAnnotation(
    cell_type = anno_block(
      gp = gpar(fill = cell_type_colors),
      labels = names(cell_type_colors),
      labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      height = unit(4, "mm")
    )
  ),

  heatmap_legend_param = list(
    at = c(min(mat), max(mat)),
    title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
    title_position = "topleft", 
    border = NA, 
    legend_height = unit(20, "mm"),
    labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
    grid_width = unit(2, "mm")
  ),
  
  # Split rows and columns by cell type
  row_split = cell_types_reset_level_cnts,
  column_split = cell_types_reset_level_cnts,
)

dev.off()


# (p <- plot_n_interactions())
# ggsave_publication("n_interactions", plot = p, width = 7, height = 6)

## b. top ligand-receptor pair/pwy per cell type  ----


sig_data <-
  cellchat %>% 
  subsetCommunication(thresh = NA) %>% 
  as_tibble() %>% 
  group_by(interaction_name_2) %>% 
  mutate(prob.norm = prob / max(prob)) %>% 
  ungroup() 
# %>% 
#   mutate(
#     source = factor(source, c("NB", "pDC", "M", "B", "T", "NK", "SC", "E")),
#     target = factor(target, c("NB", "pDC", "M", "B", "T", "NK", "SC", "E")),
#   )

signif = 0.05

vis_data_list <- list()

# by cell type and condition
for (cell_type_condition in unique(cell_types_reset_level)) {
  vis_data <- 
    sig_data %>%
    filter(
      pathway_name %in% cellchat@netP$pathways, 
      source == cell_type_condition,  # Using `==` to filter for the specific cell type and condition
      target != cell_type_condition,
      pval <= signif
    ) %>%
    group_by(interaction = interaction_name_2) %>%  # Group by cell_type_condition and interaction
    summarise(
      prob.norm = sum(prob.norm),
      count = n(),
      pathway_name = first(pathway_name),
      .groups = 'drop'
    ) %>%
    ungroup() %>%
    mutate(prob.avg = prob.norm / count) %>% 
    mutate(interaction = fct_reorder(interaction, prob.norm))  %>% # Order interactions by prob.norm
    mutate(cell_type_condition = cell_type_condition)
  
  vis_data_list[[cell_type_condition]] <- vis_data  # Store each vis_data in a list
}

# by cell type 
for (cell_type in cell_types_reset) {
  vis_data <- 
    sig_data %>%
    filter(
      pathway_name %in% cellchat@netP$pathways, 
      source %in% rev(grep(cell_type, levels(cellchat@idents), value = T)),
      !target %in% rev(grep(cell_type, levels(cellchat@idents), value = T)),
      pval <= signif
    ) %>%
    group_by(interaction = interaction_name_2) %>%  # Group by cell_type_condition and interaction
    summarise(
      prob.norm = sum(prob.norm),
      count = n(),
      pathway_name = first(pathway_name),
      .groups = 'drop'
    ) %>%
    ungroup() %>%
    mutate(prob.avg = prob.norm / count) %>% 
    mutate(interaction = fct_reorder(interaction, prob.norm))  %>% # Order interactions by prob.norm
    mutate(cell_type_condition = cell_type)
  
  vis_data_list[[cell_type]] <- vis_data  # Store each vis_data in a list
}

# Combine all vis_data into a single data frame
combined_vis_data <- bind_rows(vis_data_list)
combined_vis_data$cell_type_condition_abbr <- combined_vis_data$cell_type_condition
combined_vis_data$condition <- combined_vis_data$cell_type_condition

combined_vis_data <- combined_vis_data %>%
  mutate(
    cell_type_condition_abbr = gsub(" cells", "", cell_type_condition_abbr),
    cell_type_condition_abbr = gsub("Monocyte-derived Mph", "Mono Mph", cell_type_condition_abbr),
    cell_type_condition_abbr = gsub("Alveolar Mph MT-positive", "Alv Mph", cell_type_condition_abbr),
    cell_type_condition_abbr = gsub("Migratory DC", "DC", cell_type_condition_abbr),
    cell_type_condition_abbr = gsub("EC general capillary", "EC", cell_type_condition_abbr),
    cell_type_condition_abbr = gsub("unstim", "Unstim", cell_type_condition_abbr)
  )

combined_vis_data <- combined_vis_data %>%
  mutate(condition = ifelse(str_detect(cell_type_condition_abbr, "\\|"), 
                            sub(".*\\|", "", cell_type_condition_abbr), 
                            "All"))

# Make plot

# B, CD4_T, CD8_T, NK, Mono_Mph, Alv_Mph, DC, EC

closest_higher_multiple_of_5 <- function(x) {
  return(ceiling(x / 5) * 5)
}


plot_contribution_celltype_condition <- function(i=i){

  BASE_TEXT_SIZE_MM = 1.8
  BASE_LINE_SIZE = 0.3
  BASE_TEXT_SIZE_PT_strip = 8
  BASE_TEXT_SIZE_PT_x = 6
  BASE_TEXT_SIZE_PT_y = 5
  
  tmp <- combined_vis_data %>%
    filter(cell_type_condition %in% rev(grep(cell_types_reset[i], levels(cellchat@idents), value = T)))
  
  
  num <- tmp %>%
    filter(!condition %in% "All") %>%
    summarise(max_prob_norm = max(prob.norm, na.rm = TRUE)) %>%
    pull(max_prob_norm)
  
  closest_num <- closest_higher_multiple_of_5(num)
  
  p <- ggplot(tmp, aes(interaction, prob.norm)) +
    geom_col(fill = "#ffd500", width = 0.8) +
    geom_hline(yintercept = 0, 
               color = "grey92",
               linewidth = BASE_LINE_SIZE
    ) +
    geom_text(
      aes(y = (closest_num/2), label = pathway_name),
      size = BASE_TEXT_SIZE_MM,
      fontface = "italic"
    ) +
    xlab("Ligand-Receptor Pair") +
    ylab(paste0("Total outgoing communication score from ", cell_types_reset_abbr[i])) +
    coord_flip() +
    scale_y_continuous(limits = c(0, closest_num)) +
    facet_grid(~ factor(condition, levels = unique(condition)),
               space = "free", scales = "free") + 
    theme(
      panel.background = element_rect(fill = "white"),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = BASE_TEXT_SIZE_PT_strip),
      axis.title.x = element_text(hjust = 1),
      # axis.ticks.y = element_blank(),
      panel.grid.major.x = element_line(
        color = "grey92",
        linewidth = BASE_LINE_SIZE
      ),
      axis.text.y = element_text(size = BASE_TEXT_SIZE_PT_y),  
      axis.text.x = element_text(size = BASE_TEXT_SIZE_PT_x), 
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = BASE_LINE_SIZE),  # Adjust border thickness
      axis.ticks = element_line(linewidth = BASE_LINE_SIZE) 
    )
  
  
  return(plot = p)
  
}


plot_contribution_celltype <- function(i=i){
  
  BASE_TEXT_SIZE_MM = 1.8
  BASE_LINE_SIZE = 0.3
  BASE_TEXT_SIZE_PT_strip = 8
  BASE_TEXT_SIZE_PT_x = 6
  BASE_TEXT_SIZE_PT_y = 5
  
  tmp <- combined_vis_data %>%
    filter(cell_type_condition_abbr == cell_types_reset_abbr[i]) %>%
    mutate(interaction = fct_reorder(interaction, prob.norm, .desc = FALSE))
  
  # num <- tmp %>%
  #   summarise(max_prob_norm = max(prob.norm, na.rm = TRUE)) %>%
  #   pull(max_prob_norm)
  
  num <- combined_vis_data %>%
    filter(condition == "All") %>%
    summarise(max_prob_norm = max(prob.norm, na.rm = TRUE)) %>%
    pull(max_prob_norm)
  
  closest_num <- closest_higher_multiple_of_5(num)
  
  p <- ggplot(tmp, aes(interaction, prob.norm)) +
    geom_col(fill = "#ffd500", width = 0.8) +
    geom_hline(yintercept = 0, 
               color = "grey92",
               linewidth = BASE_LINE_SIZE
    ) +
    geom_text(
      aes(y = ceiling(closest_num / 2) , label = pathway_name),
      size = BASE_TEXT_SIZE_MM,
      fontface = "italic"
    ) +
    xlab("Ligand-Receptor Pair") +
    ylab(paste0("Total outgoing communication score from ", cell_types_reset_abbr[i])) +
    coord_flip() +
    scale_y_continuous(limits = c(0, closest_num)) +
    facet_grid(~ factor(condition, levels = unique(condition)),
               space = "free", scales = "free") + 
    theme(
      panel.background = element_rect(fill = "white"),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = BASE_TEXT_SIZE_PT_strip),
      axis.title.x = element_text(hjust = 1),
      # axis.ticks.y = element_blank(),
      panel.grid.major.x = element_line(
        color = "grey92",
        linewidth = BASE_LINE_SIZE
      ),
      axis.text.y = element_text(size = BASE_TEXT_SIZE_PT_y),  
      axis.text.x = element_text(size = BASE_TEXT_SIZE_PT_x), 
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = BASE_LINE_SIZE),  # Adjust border thickness
      axis.ticks = element_line(linewidth = BASE_LINE_SIZE) 
    )
  
  
  return(plot = p)
  
}

for (i in seq(1, 8)) {
  
  # Set the width based on the value of i
  height <- if (i %in% c(seq(1,4), 8)) {
    6.54  # Width for i = 1, 2, 3, 4, 8
  } else {
    8.96  # Width for i = 5, 6, 7
  }
  
  pdf(paste0("./integration_2/integration_2_leiden/cellchat_de/plot_contribution_celltype_condition_", 
             cell_types_reset_abbr[i], ".pdf"),
      width = 8.96,
      height = height)
  
  # Create the plots
  plot1 <- plot_contribution_celltype(i=i)
  plot2 <- plot_contribution_celltype_condition(i=i)
  
  # Combine the plots and print
  print(plot_grid(plot1, plot2, ncol = 2, rel_widths = c(1, 3)))
  
  dev.off()
}


## c. dot plot top ligand-receptor pair/pwy per cell type and condition  ----

highest_records <- combined_vis_data %>%
  filter(condition != "All") %>% 
  group_by(cell_type_condition_abbr) %>% 
  filter(prob.norm == max(prob.norm, na.rm = TRUE)) %>% 
  ungroup() 

source_type <-  rev(grep(cell_types_reset[1], levels(cellchat@idents), value = T))
source_type <- c("B cells|LPS_4h", 
                 "B cells|LPS_18h", 
                 "B cells|Ig_4h", 
                 "B cells|Ig_18h")
interactions <- c("LGALS9 - CD45",
                  "IL7 - (IL7R+IL2RG)",
                  "MIF - (CD74+CD44)")
filtered_interactions <- highest_records %>%
  filter(cell_type_condition %in% source_type)



signif = 0.05
#prepare data
vis_data <- 
  sig_data %>% 
  filter(
    interaction_name_2 %in% interactions, 
    source %in% source_type,
    pval <= signif
  ) %>%
  mutate(
    pathway_name = factor(pathway_name),
    # source = fct_recode(source, "from NB (source) to" = "NB"),
  )


vis_data <- vis_data %>%
  mutate(
    source = gsub(" cells", "", source),
    source = gsub("Monocyte-derived Mph", "Mono Mph", source),
    source = gsub("Alveolar Mph MT-positive", "Alv Mph", source),
    source = gsub("Migratory DC", "DC", source),
    source = gsub("EC general capillary", "EC", source),
    source = gsub("unstim", "Unstim", source),
    source = gsub("\\|", " | ", source),
    
    
    target = gsub(" cells", "", target),
    target = gsub("Monocyte-derived Mph", "Mono Mph", target),
    target = gsub("Alveolar Mph MT-positive", "Alv Mph", target),
    target = gsub("Migratory DC", "DC", target),
    target = gsub("EC general capillary", "EC", target),
    target = gsub("unstim", "Unstim", target),
    target = gsub("\\|", " | ", target),
    
    
    source = factor(source, levels = cell_types_reset_abbr_level),
    target = factor(target, levels = cell_types_reset_abbr_level)
  )


# export source data
# vis_data %>% 
#   select(source, target, pathway = pathway_name,
#          ligand_receptor_pair = interaction_name_2,
#          relative_comm_score = prob.norm) %>% 
#   save_table("source_data/figure_3c", "Figure 3c")

library(latex2exp)

BASE_TEXT_SIZE_MM = 1.8
BASE_LINE_SIZE = 0.2
BASE_TEXT_SIZE_PT_strip = 8
BASE_TITLE_SIZE_PT = 9
BASE_TEXT_SIZE_PT = 7



# make plot  
ggplot(vis_data, aes(target, interaction_name_2)) +
  geom_point(
    aes(color = prob.norm, size = pmin(1, -log10(pval)))
  ) + 
  scale_radius(
    TeX("-log_{10} p-Value_{cap.}"),
    range = c(1, 5),
    guide = "none"
  ) +
  scale_color_distiller(
    "Relative \ncommunication score",
    type = "seq",
    palette = "YlOrBr",
    direction = 1
  ) +
  xlab("Cell type (target)") +
  ylab("Ligand-receptor pair\nand pathway") +
  facet_grid(
    vars(pathway_name),
    vars(source),
    scales = "free",
    space = "free",
    switch = "both"
  ) + 
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                               size = BASE_TEXT_SIZE_PT, 
                               margin = margin(t = 1, b = 1, r = 1, l = 1, unit = "mm")),
    axis.text.y = element_text(size = BASE_TEXT_SIZE_PT),
    axis.ticks.length.x = unit(2, "mm"),
    axis.ticks.x = element_line(
      arrow = arrow(length = unit(1,"mm"), ends = "first", type = "closed")
    ),
    axis.title = element_text(size = BASE_TITLE_SIZE_PT),
    strip.background.y = element_rect(
      fill = "white"
      # ,
      # size = BASE_LINE_SIZE # thickness of the border of the strip
    ),
    strip.text = element_text(size = BASE_TEXT_SIZE_PT_strip),
    strip.background.x = element_rect(
      fill = "white"
      # ,
      # size = BASE_LINE_SIZE # thickness of the border of the strip
    ),
    # strip.switch.pad.grid = unit(0, "pt"), # distance of the grid with the strip 
    legend.text = element_text(size = BASE_TEXT_SIZE_PT),
    legend.title = element_text(size = BASE_TEXT_SIZE_PT),
    legend.key.width = unit(3.5, "mm"),
    legend.key.height = unit(3.5, "mm"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = BASE_LINE_SIZE),
    axis.ticks = element_line(linewidth = BASE_LINE_SIZE) 
   ) 


## d. sending, receiving ----

# cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
unique((highest_records$interaction))
centralities <- cellchat@netP$centr[["MIF_CD74_CD44"]]

mat <- 
  centralities %>% 
  as_tibble() %>% 
  select(
    sender = outdeg,
    receiver = indeg,
    mediator = flowbet,
    influencer = info
  ) %>%
  mutate(
    across(everything(), ~. / max(.)),
    cell_type = names(centralities$outdeg) %>% 
      factor(levels = cell_types_reset_level)
  ) %>%
  arrange(cell_type) %>% 
  column_to_rownames("cell_type") %>% 
  as.matrix() %>% 
  t()

for (i in 1:length(colnames(mat))) {
  # # rownames
  # rownames(mat)[i] <- gsub(" cells", "", rownames(mat)[i])
  # rownames(mat)[i] <- gsub("Monocyte-derived Mph", "Mono Mph", rownames(mat)[i])
  # rownames(mat)[i] <- gsub("Alveolar Mph MT-positive", "Alv Mph", rownames(mat)[i])
  # rownames(mat)[i] <- gsub("Migratory DC", "DC", rownames(mat)[i])
  # rownames(mat)[i] <- gsub("EC general capillary", "EC", rownames(mat)[i])
  # rownames(mat)[i] <- gsub("unstim", "Unstim", rownames(mat)[i])
  # colnames
  colnames(mat)[i] <- gsub(" cells", "", colnames(mat)[i])
  colnames(mat)[i] <- gsub("Monocyte-derived Mph", "Mono Mph", colnames(mat)[i])
  colnames(mat)[i] <- gsub("Alveolar Mph MT-positive", "Alv Mph", colnames(mat)[i])
  colnames(mat)[i] <- gsub("Migratory DC", "DC", colnames(mat)[i])
  colnames(mat)[i] <- gsub("EC general capillary", "EC", colnames(mat)[i])
  colnames(mat)[i] <- gsub("unstim", "Unstim", colnames(mat)[i])
}

# # export source data
# mat %>% 
#   as_tibble(rownames = "centrality_measure") %>% 
#   save_table(source_data_filename, source_data_sheet)

# make plot
Heatmap(
  mat, 
  col = RColorBrewer::brewer.pal(9, "PuBu"),
  name = "importance",
  
  cluster_rows = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
  
  cluster_columns = FALSE, 
  column_names_side = "top",
  column_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
  column_title = str_glue("cell type"),
  column_title_side = "top",
  column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
  
  # width = unit(20, "mm"),
  # height = unit(10, "mm"), 
  
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT, fontface = "plain"), 
    title_position = "leftcenter-rot", 
    border = NA, 
    at = c(0, 1), 
    legend_height = unit(10, "mm"), 
    labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
    grid_width = unit(2, "mm")
  )
)


plot_centrality <- function(pathway,
                            row_names_side = "left",
                            source_data_filename = NULL,
                            source_data_sheet = NULL) {
  #prepare data
  
}

(p <- plot_centrality(
  "MIF",
  source_data_filename = "source_data/figure_3d",
  source_data_sheet = "Figure 3d"
))

