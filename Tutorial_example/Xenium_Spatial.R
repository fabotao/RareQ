setwd('/home/rstudio/Projects/Rare_cell/Tutorial_example/')

library(Seurat)
library(RareQ)
library(dplyr)
library(ggplot2)


dat = readRDS('data/Xenium_Mouse_brain.RDS')
count = dat@assays$RNA@counts
meta.info = dat@meta.data

sc_object <- CreateSeuratObject(count=count, project = "sc_object", min.cells = 3)
sc_object <- NormalizeData(sc_object, scale.factor = 80) %>% ScaleData()  #
sc_object <- RunPCA(sc_object, features = dimnames(sc_object)[[1]], npcs = 50)
sc_object <- FindNeighbors(object = sc_object,
                           k.param = 20,
                           compute.SNN = F,
                           prune.SNN = 0,
                           reduction = "pca",
                           dims = 1:20,
                           force.recalc = F, return.neighbor = T)

# Use RareQ to derive both major and rare cell clusters
cluster = FindRare(sc_object = sc_object)
sc_object$cluster = cluster
sc_object$X = meta.info$X
sc_object$Y = meta.info$Y


# Visualization
cols <- c("#532C8A","#c19f70","#f9decf","#c9a997","#B51D8D","#9e6762","#3F84AA","#F397C0",
          "#C594BF","#DFCDE4","#eda450","#635547","#C72228","#EF4E22","#f77b59","#989898",
          "#7F6874","#8870ad","#65A83E","#EF5A9D","#647a4f","#FBBE92","#354E23","#139992",
          "#C3C388","#8EC792","#0F4A9C","#8DB5CE","#1A1A1A","#FACB12","#C9EBFB","#DABE99",
          "#ed8f84","#005579","#CDE088","#BBDCA8","#F6BFCB"
)

getPalette = colorRampPalette(cols[c(2,3,1,5,6,7,8,9,11,12,13,14,16,18,19,20,22,23,24,27,28,29,30,31,34)])

cluster.df = sc_object@meta.data
p.cluster <- ggplot(data = cluster.df) + geom_point(aes(x = X, y = Y, color=factor(cluster)), size=0.001) +
  theme_bw() + theme_minimal() +
  theme(panel.grid = element_blank(), legend.position = 'none',
        axis.text = element_blank(),
        axis.title = element_blank()) +
  scale_color_manual(values = getPalette(length(unique(cluster.df$cluster))))
p.cluster

# CA2 may have the same color with CA3 in plot, user can check it with the reference cluster 82936 in dat.
cnt = table(cluster[dat$cluster==82936])
cnt
# Plot CA2 region
plot(cluster.df$X, cluster.df$Y, cex=0.1, col=ifelse(cluster.df$cluster==names(cnt)[which.max(cnt)],'red','grey'))


