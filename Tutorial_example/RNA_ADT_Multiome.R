setwd('/home/rstudio/Projects/Rare_cell/Tutorial_example/')

library(Seurat)
library(RareQ)
library(dplyr)
library(ggplot2)

# Read example data
# ## The data and preprocessing steps can be found at https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis
bm = readRDS('data/bmcite.RDS')  # Example data from data folder

DefaultAssay(bm) <- 'RNA'
bm <- NormalizeData(bm) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

DefaultAssay(bm) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the
VariableFeatures(bm) <- rownames(bm[["ADT"]])
bm <- NormalizeData(bm, normalization.method = 'CLR', margin = 2) %>%
  ScaleData() %>% RunPCA(reduction.name = 'apca')



# Identify multimodal neighbors. These will be stored in the neighbors slot,
# and can be accessed using bm[['weighted.nn']]
# The WNN graph can be accessed at bm[["wknn"]],
# and the SNN graph used for clustering at bm[["wsnn"]]
# Cell-specific modality weights can be accessed at bm$RNA.weight
bm <- FindMultiModalNeighbors(
  bm, reduction.list = list("pca", "apca"),
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
)

bm <- RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

## The spca reduction will be used to calculate neighborhood for rare cell detection
DefaultAssay(bm) = 'RNA'
bm <- RunSPCA(bm, assay = 'RNA', graph = 'wsnn')
bm <- FindNeighbors(object = bm,
                    k.param = 20,
                    compute.SNN = F,
                    prune.SNN = 0,
                    reduction = "spca",
                    dims = 1:50,
                    force.recalc = F, return.neighbor = T)

# Use RareQ to derive both major and rare cell clusters
cluster = FindRare(sc_object = bm)
bm$cluster = cluster


# Label cluster according to cluster size
cluster.cnt <- sort(table(bm$cluster))
bm$cluster_sort = factor(as.character(bm$cluster), levels=names(cluster.cnt), labels = 1:length(cluster.cnt), ordered = T)

# Visualization
cols <- c("#532C8A","#c19f70","#f9decf","#c9a997","#B51D8D","#9e6762","#3F84AA","#F397C0",
          "#C594BF","#DFCDE4","#eda450","#635547","#C72228","#EF4E22","#f77b59","#989898",
          "#7F6874","#8870ad","#65A83E","#EF5A9D","#647a4f","#FBBE92","#354E23","#139992",
          "#C3C388","#8EC792","#0F4A9C","#8DB5CE","#1A1A1A","#FACB12","#C9EBFB","#DABE99",
          "#ed8f84","#005579","#CDE088","#BBDCA8","#F6BFCB"
)

getPalette = colorRampPalette(cols[c(2,3,1,5,6,7,8,9,11,12,13,14,16,18,19,20,22,23,24,27,28,29,30,31,34)])

DimPlot(bm, group.by='cluster_sort', reduction = 'wnn.umap') +
  scale_color_manual(values = getPalette(length(unique(bm$cluster_sort))))

