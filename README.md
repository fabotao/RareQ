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
  obj = readRDS('data/Jurkat.obj)  # Example data from data folder
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
For scATAC-seq data, count matrix is preprocessed by Signac and Seurat 

```R
  library(RareQ)
  library(Seurat) 
  library(Signac)
  
  # Preprocessing scATAC-seq data from peak count matrix
  ATAC <- read_mtx('peak_cell.mtx')
  
  set.seed(123)
  n_peaks <- dim(ATAC$X)[1]  
  
  # Generate random genomic range, User can also use default genomic range using standard procedure.
  random_gr <- GRanges(
    seqnames = paste0("chr", sample(1:22, n_peaks, replace = TRUE)),  
    ranges = IRanges(
      start = sample(1:1e6, n_peaks),                                 
      width = sample(100:500, n_peaks, replace = TRUE)                
    ),
    strand = "*",
    score = rnorm(n_peaks, mean = 10, sd = 2)                       
  )
  
  # Create ChromatinAssay
  chrom_assay <- CreateChromatinAssay(
    counts = ATAC$X,          
    ranges = random_gr,       
    genome = "mm10",        
    assay = "Peak"          
  )
  
  # Create Seurat object
  seurat_atac <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "Peak",       
    project = "ATAC_Project"
  )
  
  seurat_atac <- RunTFIDF(seurat_atac)         
  seurat_atac <- FindTopFeatures(seurat_atac, min.cutoff = 10)  
  seurat_atac <- RunSVD(seurat_atac)            
  
  seurat_atac <- RunUMAP(seurat_atac, reduction = "lsi", dims = 2:30)
  seurat_atac <- FindNeighbors(seurat_atac, reduction = "lsi", dims = 2:30,
                               k.param = 20, compute.SNN = F,return.neighbor = T)
  
  cluster = FindRare(sc_object = seurat_atac, assay = 'Peak')

```



## Citation


