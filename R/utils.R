
.get.compact <- function(knn.matrix, id){
  len <- length(id)
  if(len < 20){
    snn.id <- c(knn.matrix[id,1:len])
    len.base <- max(c(len, 5))
    return(sum(snn.id %in% id)/(len.base * len.base))
  }else {
    snn.id <- c(knn.matrix[id,1:20])
    len.base <- 20
    return(sum(snn.id %in% id)/(len.base * len))
  }
}

