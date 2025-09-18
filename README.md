# RareQ
RareQ (R package) is a network-propagation-based algorithm for rare cell type identification, that demonstrates superior performance and higher scalability to existing algorithms.


## Required R modules
```R
    R >= 4.0.0
    Seurat >= 4.0.2
    Signac >= 1.9.0  # For preprocessing scATAC-seq data
```

## Installation
```R
  # Install in R with devtools
  library(devtools)
  install_github('fabotao/RareQ')
```

## Usage
For scRNA-seq, CITE-seq data and cell-segmented spatial data, count matrix is preprocessed by Seurat package

```R
  library(RareQ)
  library(Seurat) 
  
  # Read example data
  obj = readRDS('data/Jurkat.RDS)  # Example data from data folder
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
  DimPlot(sc_object, group.by='cluster')
  
```



## Citation


