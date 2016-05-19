#' Rearrange matrix based on binary partition vector
#'
#' @param X numeric array or matrix
#' @param partition binary partition vector with 1 representing block 1
#' and 0 representing block 2
#' @param fix.first don't change the order of the first column
#'
#' @return numeric matrix with a minimal row sum variance
#'
#' @export
#'
#' @author Kris Boudt, \email{kris.boudt@@vub.ac.be}
#' @author Steven Vanduffel, \email{steven.vanduffel@@vub.ac.be}
#' @author Kristof Verbeken, \email{kristof.verbeken@@vub.ac.be}
rearrangepartition <- function(X, partition) {
  if (!is.matrix(X)) {
    stop("matrix argument is a numeric matrix")
  }

  if (!is.vector(partition)) {
    stop("partition is not a vector")
  }

  if (!all(partition %in% c(T, F))) {
    stop("partition vector has elements that are not binary")
  }

  if (length(partition) != ncol(X)) {
    stop("Partition incompatible with matrix.")
  }

  # make sure the first block is the smaller one
  if (sum(partition) > length(partition)/2) {
    partition=!partition
  }
  
  if (sum(partition)>0) # if the partition is nondegenerate, rearrange
  {
    block1 <- X[, partition]
    block2 <- X[,!partition]
    # rearrange block 1
    block1 <- rearrange(block1, block2)
    X[, partition] <- block1
  }
  
  return(X)
}
