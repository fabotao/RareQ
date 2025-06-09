#'
#' @title Create a pseudo-Seurat object to accommodate emerging data modalities in rare cell identification.
#'
#' @description Use provided low-dimensional embeddings processed from any single cell data.
#'
#' @usage CreateObject(
#'   embed.matrix,
#'   assay.name='RNA',
#'   dims = NULL
#' )
#'
#' @param embed.matrix Low-dimensional cell embeddings like PCA, SVD ....
#' @param assay.name The default assay name used when creating the pseudo-Seurat object.
#' @param dims Dimensions used in finding nearest neighbors. Such as 1:50, 2:30,...
#' @return pseudo-Seurat object that can be used by FindRare to infer rare cell populations.
#' @export
#' @examples
#'
CreateObject <- function(embed.matrix, assay.name='RNA', dims=NULL){
  shape = dim(embed.matrix)

  N = max(shape)  # Assume the cells are much larger than the reduced-dimension
  p = min(shape)

  if(shape[1] > shape[2]){
    embeddings = embed.matrix
  }else{
    embeddings = t(embed.matrix)
  }
  if(is.null(dimnames(embeddings)[[1]])){
    cell.id <- paste0('cell', 1:N)
  }else{
    cell.id <- dimnames(embeddings)[[1]]
  }
  count <- matrix(rbinom(N * p, size=50, prob=0.1), ncol=N, nrow=p)
  dimnames(count) <- list(paste0('feature', 1:p), cell.id)

  reduction.key <- "PC_"
  rownames(x = embeddings) <- dimnames(embeddings)[[1]]
  colnames(x = embeddings) <- paste0(reduction.key, 1:dim(embeddings)[2])
  reduction.data <- SeuratObject::CreateDimReducObject(
    embeddings = embeddings,
    #loadings = NULL,
    assay = assay.name,
    stdev = 0,
    key = 'PC_',
    #misc = 0
  )

  sc_object <- CreateSeuratObject(count=count, project = "sc_object", min.cells = 1)
  sc_object[['pca']] <- reduction.data
  if(is.null(dims)){
    sc_object <- FindNeighbors(sc_object, reduction = 'pca', dims=1:p, k.param = 20, return.neighbor = T)
  }else{
    sc_object <- FindNeighbors(sc_object, reduction = 'pca', dims=dims, k.param = 20, return.neighbor = T)
  }
  return(sc_object)
}







#'
#' @title Compute the Q values for each cell.
#'
#' @description Compute the Q values for each cell for users to choose rare cell clusters.
#'
#' @usage CreateObject(
#'   embed.matrix,
#'   assay.name='RNA',
#'   dims = NULL
#' )
#'
#' @param embed.matrix Low-dimensional cell embeddings like PCA, SVD ....
#' @param assay.name The default assay name used when creating the pseudo-Seurat object.
#' @param dims Dimensions used in finding nearest neighbors. Such as 1:50, 2:30,...
#' @return pseudo-Seurat object that can be used by FindRare to infer rare cell populations.
#' @export
#' @examples
#'
ComputeQ <- function(sc_object, assay='RNA'){
  assay.all <- Assays(sc_object)
  if(!assay %in% assay.all){
    stop(paste0("The ", assay," assay does not exist. Please choose from ", assay.all))
  }
  nn.slot <- paste0(assay,'.nn')
  if(!nn.slot %in% names(sc_object@neighbors)){
    stop(paste0("The ", nn.slot," slot does not exist. Please run FindNeighbors with return.neighbor = TRUE"))
  }
  knn.matrix <- sc_object@neighbors[[nn.slot]]@nn.idx
  knn.dist <- sc_object@neighbors[[nn.slot]]@nn.dist
  N = dim(knn.matrix)[1]
  cell.cpt <- apply(knn.matrix[,1:6],1,function(x){.get.compact(knn.matrix = knn.matrix, x)})
  return(cell.cpt)
}

