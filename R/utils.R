
.get.compact <- function(knn.matrix, id, k.param=20){
  len <- length(id)
  if(len < k.param){
    snn.id <- c(knn.matrix[id,1:len])
    len.base <- max(c(len, 5))
    return(sum(snn.id %in% id)/(len.base * len.base))
  }else {
    snn.id <- c(knn.matrix[id,1:k.param])
    len.base <- k.param
    return(sum(snn.id %in% id)/(len.base * len))
  }
}

