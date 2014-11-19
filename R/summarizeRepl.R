summarizeRepl <- function(data,type=median){

  annot <- unique(rownames(data))
  dat <- matrix(NA, ncol=dim(data)[2], nrow=length(annot))

  for (i in 1:length(annot)) {
    id <- which(rownames(data) == annot[i])
    dat[i,] <- apply(data[id, ], 2, type, na.rm=T)
  }

  colnames(dat) <- colnames(data)
  rownames(dat) <- annot

  return(dat)
}

