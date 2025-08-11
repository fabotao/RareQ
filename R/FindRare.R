#'
#' @title Identify rare cell populations in single cell genomic data and cell-segmented spatial data.
#'
#' @description RareQ identifies rare cell types via mutual connectivity (Q)-based network propagation.
#'
#' @usage FindRare(
#'   sc_object,
#'   assay = "RNA",
#'   k = 6,
#'   Q_cut = 0.6,
#'   ratio = 0.2,
#' )
#'
#' @param sc_object A seurat object containing assays, reduction, meta data etc.
#' @param assay The Assay object of the input Seurat object for inferring rare cell types. Default is 'RNA'. User can adjust according to data modality.
#' @param k The size of neighborhood in computing Q values for rare cell detection. k should be between 5 and the total nearest neighbors (e.g., 20) obtained in processing step. The default value is 6.
#' @param Q_cut The default Q value for retaining rare clusters. Q should be between 0.5 and 1. The default value is 0.6.
#' @param ratio The threshold for merging clusters. When the proportion of a cluster's connections with its neighboring clusters among its total connections is larger than this threshold, merge the two cell groups. The default value is 0.2.
#' @param max_iter The maximum iterations for propagation if convergence is not reached. Default value is 100.
#' @return Cluster assignments of cells including major and rare cell types.
#' @export
#' @examples
#' ## Starts from counts
#' library(RareQ)
#' library(Seurat)
#'
#' ## Data prepocessing
#' sc_object <- CreateSeuratObject(count=counts, project = "sc_object", min.cells = 3)
#' sc_object$percent.mt <- PercentageFeatureSet(sc_object, pattern = "^MT-")
#' sc_object <- subset(sc_object, percent.mt<20)
#' sc_object <- NormalizeData(sc_object)
#' sc_object <- FindVariableFeatures(sc_object, nfeatures=2000)
#' sc_object <- ScaleData(sc_object, do.center = T)
#' sc_object <- RunPCA(sc_object, npcs=50)
#' sc_object <- FindNeighbors(object = sc_object,
#'                            k.param = 20,
#'                            compute.SNN = F,
#'                            prune.SNN = 0,
#'                            reduction = "pca",
#'                            dims = 1:50,
#'                            force.recalc = F, return.neighbor = T)
#' ## Find rare cell types
#' clusters <- FindRare(sc_object)
#'

FindRare <- function(sc_object, assay='RNA', k=6, Q_cut=0.6, ratio=0.2, max_iter=100){
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
  k.param = dim(knn.matrix)[2]  # The number of nearest neighbors obtained in data preprocessing.
  N = dim(knn.matrix)[1]
  cell.cpt <- apply(knn.matrix[,1:k],1,function(x){.get.compact(knn.matrix = knn.matrix, x, k.param=k.param)})
  clu0 <- 1:N
  change = T
  iter = 0
  while(change){
    clu <- apply(knn.matrix[,1:k],1,function(x){
      cpt.x <- cell.cpt[x]
      return(clu0[min(x[which(cpt.x == max(cpt.x))])])
    })
    iter = iter + 1
    if(all(clu==clu0) | iter > max_iter){change = F}
    clu0 = clu
  }
  knn.vote <- function(cluster){
    change = T
    cluster0 <- c(cluster)
    iter = 0
    while(change){
      cluster.adj <- knn_vote(V=cluster0, M=knn.matrix[,1:k], steps = 1)
      iter = iter + 1
      if(all(cluster0 == cluster.adj) | iter > max_iter){change = F}
      cluster0 = cluster.adj
    }
    return(cluster0)
  }
  cluster0 <- knn.vote(clu0)
  cluster.connectivity <- tapply(1:N, cluster0, function(x){.get.compact(knn.matrix = knn.matrix, x, k.param=k.param)})
  cluster.count <- table(cluster0)
  bad.communities <- names(cluster.connectivity)[cluster.connectivity < 0.5 | (cluster.count < 5)] #
  good.cl <- setdiff(names(cluster.connectivity), c(bad.communities))
  good.id <- which(cluster0 %in% good.cl)
  change = T
  iter = 0
  while(change){
    cluster.adj <- knn_vote(V=cluster0, M=knn.matrix[,], steps = 1)
    cluster0.old <- cluster0
    cluster0 <- ifelse(1:N %in% good.id, cluster0, cluster.adj) # Keep the good cells fixed
    iter = iter + 1
    if(all(cluster0 == cluster0.old) | iter > max_iter * 2){change=F} #
  }
  cell.count <- table(knn.matrix)
  cluster.edge <- tapply(1:N, cluster0, function(x){
    len <- min(c(k.param, length(x)))
    return(sum(cell.count[(x)])/(length(x) * len))
  })
  cluster.connectivity <- tapply(1:N, cluster0, function(x){.get.compact(knn.matrix = knn.matrix, x, k.param=k.param)})
  cluster.count <- table(cluster0)
  bad.cluster <- as.integer(names(cluster.edge)[cluster.edge >= 2 | cluster.connectivity < 1 | cluster.count < 8])
  if(length(bad.cluster)==0){
    cluster.good <- (cluster0)
  }else{
    combine = T
    while(combine){
      bad.cluster.old = bad.cluster
      cluster.id <- tapply(1:N, cluster0, function(x){x})
      new.cluster <- unlist(lapply(bad.cluster, function(x){
        id.x = cluster.id[[as.character(x)]]
        len = length(id.x)
        neib.id <- c(knn.matrix[id.x,1:min(c(len, k.param))])
        neib.cell <- unique(neib.id)
        cluster.neib.count = table(cluster0[neib.id])
        cl.oth <- names(cluster.neib.count)[names(cluster.neib.count) != as.character(x)]
        if(length(cl.oth)>0){
          cand.cl <- (cl.oth)[which.max(cluster.neib.count[cl.oth])]
          self.edge <- cluster.edge[[as.character(x)]]
          cand.connectivity = cluster.connectivity[[cand.cl]]
          self.connectivity = cluster.connectivity[[as.character(x)]]
          comb.connectivity = .get.compact(knn.matrix = knn.matrix, c(cluster.id[[cand.cl]], id.x), k.param=k.param)
          neib.cluster.id <- cluster.id[[cand.cl]]
          neib.cell.cand <- neib.cell[cluster0[neib.cell] == as.integer(cand.cl)]
          if(is.na(self.connectivity) | is.na(self.edge)){
            return(as.integer(cand.cl))
            }
          if((self.connectivity < 0.6 & mean(cell.cpt[id.x] < Q_cut)) | self.edge >= 2){
            return(as.integer(cand.cl))
          }else if(self.connectivity < 0.65 & (comb.connectivity - cand.connectivity >= 0.01) & (comb.connectivity - self.connectivity >= 0.01) & mean(cell.cpt[id.x] < Q_cut)){
            return(as.integer(cand.cl))
          }else if((comb.connectivity - cand.connectivity >= 0.05) & (comb.connectivity - self.connectivity >= 0.05) & mean(cell.cpt[id.x] < Q_cut)){
            return(as.integer(cand.cl))
          }else{
            if(len < 8){
              delta.dist <- (apply(knn.dist[id.x,,drop=F],1,diff))
              ord.delta.dist <- apply(delta.dist, 2, function(x){order(x, decreasing = T)[2]})
              if(sum(ord.delta.dist == len) == len & self.edge < 1.1){
                return(x)
              }else{
                return(as.integer(cand.cl))
              }
            }else if(len <= (c(N * 0.01))){
              if(comb.connectivity > cand.connectivity &
                 comb.connectivity > self.connectivity &
                 cluster.count[[cand.cl]] <= (c(N * 0.01)) &
                 cluster.neib.count[cand.cl]/cluster.neib.count[as.character(x)] >= ratio ){ # & mean(cell.cpt[id.x] < 0.6)
                return(as.integer(cand.cl))
              }else if(comb.connectivity > cand.connectivity &
                       comb.connectivity > self.connectivity &
                       cluster.count[[cand.cl]] > (c(N * 0.01)) &
                       cluster.neib.count[cand.cl]/cluster.neib.count[as.character(x)] >= min(c(ratio + 0.1, 1)) & mean(cell.cpt[id.x]) < Q_cut){
                return(as.integer(cand.cl))
              }else{
                return(x)
              }
            }else if(len > (c(N * 0.01)) & len <= N * 0.03){
              if(self.connectivity < 0.8 &
                 comb.connectivity > cand.connectivity &
                 comb.connectivity > self.connectivity &
                 mean(cell.cpt[id.x]) < Q_cut){
                return(as.integer(cand.cl))
              }else if(comb.connectivity > cand.connectivity &
                       comb.connectivity > self.connectivity &
                       cluster.count[[cand.cl]] > (c(N * 0.01)) &
                       cluster.neib.count[cand.cl]/cluster.neib.count[as.character(x)] >= ratio &
                       mean(cell.cpt[id.x]) < Q_cut){ #
                return(as.integer(cand.cl))
              }else{
                return(x)
              }
            }else{
              if(self.connectivity < 0.85 &
                 comb.connectivity > cand.connectivity &
                 comb.connectivity > self.connectivity &
                 cluster.neib.count[cand.cl]/cluster.neib.count[as.character(x)] >= ratio){
                return(as.integer(cand.cl))
              }else if(comb.connectivity > cand.connectivity &
                       comb.connectivity > self.connectivity &
                       cluster.count[[cand.cl]] > (c(N * 0.01)) &
                       cluster.neib.count[cand.cl]/cluster.neib.count[as.character(x)] >= ratio){
                return(as.integer(cand.cl))
              } else{
                return(x)
              }
            }
          }
        }else{
          if(len < 8){
            delta.dist <- (apply(knn.dist[id.x,,drop=F],1,diff))
            cluster.neib.count = table(cluster0[c(knn.matrix[id.x,1:(len + 1)])])
            cl.oth <- names(cluster.neib.count)[names(cluster.neib.count) != as.character(x)]
            cand.cl <- (cl.oth)[which.max(cluster.neib.count[cl.oth])]
            ord.delta.dist <- apply(delta.dist, 2, function(x){order(x, decreasing = T)[2]})
            if(sum(ord.delta.dist == len) == len &
               .get.compact(knn.matrix = knn.matrix, cluster.id[[cand.cl]], k.param=k.param) >= .get.compact(knn.matrix = knn.matrix, c(cluster.id[[cand.cl]], id.x), k.param=k.param)){
              return(x)
            }else{
              return(as.integer((cl.oth)[which.max(cluster.neib.count[cl.oth])]))
            }
          }else{
            return(x)
          }
        }
      }))
      cluster0.new = cluster0
      for(m in 1:length(bad.cluster)){
        bad.cl = bad.cluster[m]
        new.cl = new.cluster[m]
        if(new.cl %in% bad.cluster){
          if(new.cluster[which(bad.cluster == new.cl)]==bad.cl){
            cluster0.new[cluster.id[[as.character(bad.cl)]]] = min(c(bad.cl, new.cl))
          }else{
            cluster0.new[cluster.id[[as.character(bad.cl)]]] = new.cl
          }
        }else{
          cluster0.new[cluster.id[[as.character(bad.cl)]]] = new.cl
        }
      }
      cluster0 = cluster0.new
      cluster.edge <- tapply(1:N, cluster0, function(x){
        len <- min(c(k.param, length(x)))
        return(sum(cell.count[(x)])/(length(x) * len))
      })
      cluster.count <- table(cluster0)
      cluster.connectivity <- tapply(1:N, cluster0, function(x){.get.compact(knn.matrix = knn.matrix, x, k.param=k.param)})
      bad.cluster <- as.integer(names(cluster.edge)[cluster.edge >= 2 | cluster.connectivity < 1 | cluster.count < 8])
      if(length(bad.cluster)==0 | setequal(bad.cluster, bad.cluster.old)) combine = F
    }
    cluster.good <- (cluster0)
  }
  cluster.id <- tapply(1:N, cluster.good, function(x){x})
  cluster.good.temp <- cluster.good
  cluster.connectivity <- tapply(1:N, cluster.good, function(x){.get.compact(knn.matrix = knn.matrix, x, k.param=k.param)})
  refined.cl <- (sapply(cluster.id, function(x){
    len <- length(x)
    x.cl = cluster.good[x[1]]
    if(len >= N * 0.01){
      return(x.cl)
    }else{
      neib.id <- c(knn.matrix[x,1:min(c(len, k.param))])
      neib.cell <- unique(neib.id)
      cluster.neib.count = table(cluster.good[neib.id])
      cluster.cell.count = table(cluster.good[neib.cell])
      ratio <- cluster.neib.count/cluster.cell.count
      cl.oth <- names(ratio)[names(ratio) != as.character(x.cl)]
      if(length(cl.oth) > 0){
        cand.cl <- (cl.oth)[which.max(ratio[cl.oth])]
        cand.connectivity = cluster.connectivity[[cand.cl]]
        self.connectivity = cluster.connectivity[[as.character(x.cl)]]
        comb.connectivity = .get.compact(knn.matrix = knn.matrix, c(cluster.id[[cand.cl]], x), k.param=k.param)
        neib.cell.cand <- neib.cell[cluster.good[neib.cell] == as.integer(cand.cl)]
        neib.cluster.id <- cluster.id[[cand.cl]]
        if(ratio[as.character(x.cl)] * 0.4 <= ratio[cand.cl] &
           comb.connectivity > self.connectivity & comb.connectivity > cand.connectivity & mean(cell.cpt[neib.cell.cand])>=0.5){
          if(.get.compact(knn.matrix = knn.matrix, union(x, neib.cell.cand), k.param=k.param) > .get.compact(knn.matrix = knn.matrix, c(x), k.param=k.param) &
             .get.compact(knn.matrix = knn.matrix, setdiff(neib.cluster.id, neib.cell.cand), k.param=k.param) > .get.compact(knn.matrix = knn.matrix, neib.cluster.id, k.param=k.param)){
            return(c(x.cl, neib.cell.cand))
          }else{
            return(as.integer(cand.cl))
          }
        }else{
          return(x.cl)
        }
      }else{
        return(x.cl)
      }
    }
  }))
  cluster.good.new = cluster.good.temp
  cluster.id <- tapply(1:N, cluster.good.temp, function(x){x})
  for(cl in names(cluster.id)){
    init.cl = as.integer(cl)
    new.cl = refined.cl[[cl]]
    if(length(new.cl)==1){
      if(init.cl!=new.cl){
        cluster.good.new[cluster.id[[as.character(init.cl)]]] = new.cl
      }
    }else{
      cluster.good.new[new.cl[2:length(new.cl)]] <- new.cl[1]
    }
  }
  return(cluster.good.new)
}


