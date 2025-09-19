setwd('/home/rstudio/Projects/Rare_cell/Tutorial_example/')

library(Seurat)
library(RareQ)
library(ggplot2)

# Read example data
obj = readRDS('data/Jurkat.RDS')  # Example data from data folder
counts = obj@assays$RNA@counts

# Preprocessing scRNA-seq data
sc_object <- CreateSeuratObject(count=counts, project = "sc_object", min.cells = 3)
sc_object$percent.mt <- PercentageFeatureSet(sc_object, pattern = "^MT-")
sc_object <- subset(sc_object, percent.mt<20)
sc_object <- NormalizeData(sc_object)
sc_object <- FindVariableFeatures(sc_object, nfeatures=2000)
sc_object <- ScaleData(sc_object)
sc_object <- RunPCA(sc_object, npcs=50)
sc_object <- RunUMAP(sc_object, dims=1:50)
sc_object <- FindNeighbors(object = sc_object,
                           k.param = 20,
                           compute.SNN = F,
                           prune.SNN = 0,
                           reduction = "pca",
                           dims = 1:50,
                           force.recalc = F,
                           return.neighbor = T)

# Use RareQ to derive both major and rare cell clusters
cluster <- FindRare(sc_object)
table(cluster)

sc_object$cluster = cluster

# Label cluster according to cluster size
cluster.cnt <- sort(table(sc_object$cluster))
sc_object$cluster_sort = factor(as.character(sc_object$cluster), levels=names(cluster.cnt), labels = 1:length(cluster.cnt), ordered = T)

# Visualization
cols <- c("#532C8A","#c19f70","#f9decf","#c9a997","#B51D8D","#9e6762","#3F84AA","#F397C0",
          "#C594BF","#DFCDE4","#eda450","#635547","#C72228","#EF4E22","#f77b59","#989898",
          "#7F6874","#8870ad","#65A83E","#EF5A9D","#647a4f","#FBBE92","#354E23","#139992",
          "#C3C388","#8EC792","#0F4A9C","#8DB5CE","#1A1A1A","#FACB12","#C9EBFB","#DABE99",
          "#ed8f84","#005579","#CDE088","#BBDCA8","#F6BFCB"
)

getPalette = colorRampPalette(cols[c(2,3,1,5,6,7,8,9,11,12,13,14,16,18,19,20,22,23,24,27,28,29,30,31,34)])

DimPlot(sc_object, group.by='cluster_sort') +
  scale_color_manual(values = getPalette(length(unique(sc_object$cluster_sort))))


