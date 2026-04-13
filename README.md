# RareQ
RareQ is an R package for identifying rare cell populations in single-cell and cell-segmented spatial datasets. It uses Q-value-guided network propagation and is designed to be accurate, scalable, and robust.

## Dependencies
```r
R >= 4.0.0
Seurat >= 4.0.2
Signac >= 1.9.0  # required only for scATAC-seq preprocessing
```

## Installation
```r
library(devtools)
install_github("fabotao/RareQ")
```

## Quick start (scRNA-seq)
The [Jurkat example dataset](https://github.com/fabotao/RareQ/blob/main/data/Jurkat.RDS) is included for demonstration.

```r
library(RareQ)
library(Seurat)

# Load example data
obj <- readRDS("Jurkat.RDS")
counts <- obj@assays$RNA@counts

# Preprocess scRNA-seq data
sc_object <- CreateSeuratObject(counts = counts, project = "sc_object", min.cells = 3)
sc_object <- NormalizeData(sc_object)
sc_object <- FindVariableFeatures(sc_object, nfeatures = 2000)
sc_object <- ScaleData(sc_object)
sc_object <- RunPCA(sc_object, npcs = 50)
sc_object <- RunUMAP(sc_object, dims = 1:50)
sc_object <- FindNeighbors(
  object = sc_object,
  k.param = 20,
  compute.SNN = FALSE,
  prune.SNN = 0,
  reduction = "pca",
  dims = 1:50,
  force.recalc = FALSE,
  return.neighbor = TRUE
)

# Identify major and rare clusters
cluster <- FindRare(sc_object)
table(cluster)

sc_object$cluster <- cluster
DimPlot(sc_object, group.by = "cluster")
```

## Optional: consensus clustering with `ConsensusRare`
`FindRare` is deterministic and generally robust. If you want additional stability checks, `ConsensusRare` runs `FindRare` repeatedly on shuffled cell orders and aggregates the results through consensus clustering.

> **Note**: `ConsensusRare` is slower and more memory-intensive than `FindRare`, especially for large datasets.

```r
library(RareQ)
library(Seurat)

# Load example data
obj <- readRDS("Jurkat.RDS")
counts <- obj@assays$RNA@counts

# Preprocess scRNA-seq data
sc_object <- CreateSeuratObject(counts = counts, project = "sc_object", min.cells = 3)
sc_object <- NormalizeData(sc_object)
sc_object <- FindVariableFeatures(sc_object, nfeatures = 2000)
sc_object <- ScaleData(sc_object)
sc_object <- RunPCA(sc_object, npcs = 50)
sc_object <- RunUMAP(sc_object, dims = 1:50)

# Run consensus clustering
cluster <- ConsensusRare(
  sc_object,
  assay = "RNA",
  reduction = "pca",
  dims = 1:50,
  k.param = 20,
  k = 6,
  Q_cut = 0.6,
  ratio = 0.2,
  reps = 30
)

table(cluster)
sc_object$cluster <- cluster
DimPlot(sc_object, group.by = "cluster")
```

## Tutorials
Tutorial HTML notebooks are available in the `Tutorials/` directory:

1. [scRNA_analysis](https://xiaolab-xjtu.github.io/RareQ/Tutorials/scRNA_analysis.html): scRNA-seq analysis using Jurkat data.
2. [scRNA_scATAC_analysis](https://xiaolab-xjtu.github.io/RareQ/Tutorials/scRNA_scATAC_analysis.html): joint scRNA-seq/scATAC-seq multiome analysis.
3. [scRNA_ADT_analysis](https://xiaolab-xjtu.github.io/RareQ/Tutorials/scRNA_ADT_analysis.html): CITE-seq analysis with RNA and ADT modalities.
4. [Xenium_spatial_analysis](https://xiaolab-xjtu.github.io/RareQ/Tutorials/Xenium_spatial_analysis.html): cell-segmented Xenium spatial analysis.

Optional tutorial for consensus workflow:

- [Consensus_Tutorial](https://xiaolab-xjtu.github.io/RareQ/Tutorials/Consensus_Tutorial.html): `ConsensusRare` on Jurkat scRNA-seq data.

Related tutorial datasets:

- [Tutorial example datasets](https://zenodo.org/records/17190972/files/Tutorial_example.rar?download=1)

## Simulation resources
Simulation scripts are in `Simulation/`.

- [Simulation code walkthrough](https://xiaolab-xjtu.github.io/RareQ/Simulation/Simulated_PBMC1_2_3.html)
- [Simulation datasets](https://zenodo.org/records/17190972/files/simulation.rar?download=1)

## Citation
If you use RareQ in your work, please cite the.

## License
Copyright © 2026 XiaoLab@XJTU.
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
