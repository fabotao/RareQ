setwd('/home/rstudio/Projects/Rare_cell/Tutorial_example/')

library(Seurat)
library(Signac)
library(RareQ)
library(dplyr)
library(ggplot2)
library(GenomicRanges)

# Read example data
dat = readRDS('data/Mouse_gdT.RDS')  # Example data from data folder


### RNA
sc_object <- CreateSeuratObject(count=dat$RNA, project = "sc_object", min.cells = 3)

### ATAC
ATAC <- dat$Peaks
n_peaks <- dim(ATAC)[1]

seqnames = sapply(strsplit(dimnames(ATAC)[[1]], ':', fixed = T), '[', 1)
range.info = sapply(strsplit(dimnames(ATAC)[[1]], ':', fixed = T), '[', 2)

#
random_gr <- GRanges(
  seqnames = seqnames,
  ranges = IRanges(
    start = as.integer(sapply(strsplit(range.info, '-', fixed = T), '[', 1)),
    end = as.integer(sapply(strsplit(range.info, '-', fixed = T), '[', 2))
  ),
  strand = "*",
  score = rnorm(n_peaks, mean = 10, sd = 2)                         # 随机分数
)

# 创建 ChromatinAssay
chrom_assay <- CreateChromatinAssay(
  counts = ATAC,
  ranges = random_gr,
  genome = "mm10",
  assay = "Peak"
)


sc_object[["ATAC"]] <- chrom_assay

# RNA analysis
DefaultAssay(sc_object) <- "RNA"
sc_object <- NormalizeData(sc_object) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(sc_object) <- "ATAC"
sc_object <- RunTFIDF(sc_object)
sc_object <- FindTopFeatures(sc_object, min.cutoff = 'q0')
sc_object <- RunSVD(sc_object)

## Integration with WNN
sc_object <- FindMultiModalNeighbors(sc_object, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
sc_object <- RunUMAP(sc_object, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

## The spca reduction will be used to calculate neighborhood for rare cell detection
sc_object <- RunSPCA(sc_object, assay = 'RNA', graph = 'wsnn')
sc_object <- FindNeighbors(object = sc_object,
                      k.param = 20,
                      compute.SNN = F,
                      prune.SNN = 0,
                      reduction = "spca",
                      dims = 1:50,
                      force.recalc = F, return.neighbor = T)


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

DimPlot(sc_object, group.by='cluster_sort', reduction = 'wnn.umap') +
  scale_color_manual(values = getPalette(length(unique(sc_object$cluster_sort))))
