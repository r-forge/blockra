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
rearrangepartition <- function(X, partition, fix.first = TRUE) {
  if (!is.matrix(X)) {
    stop("matrix argument is a numeric matrix")
  }

  if (!is.vector(partition)) {
    stop("partition is not a vector")
  }

  if (length(partition[!partition %in% c(0, 1)]) > 0) {
    stop("partition vector has elements that are not binary")
  }

  if (length(partition) != ncol(X)) {
    stop("Partition incompatible with matrix.")
  }

  # Either fix the order of the first partition
  # OR make sure the first block is the largest
  if ((fix.first && partition[1] == 1) ||
      (!fix.first && sum(partition) > length(partition) / 2)) {
    partition[partition == 1] <- 2
    partition[partition == 0] <- 1
    partition[partition == 2] <- 0
  }

  block1 <- X[, partition == 1]
  block2 <- X[, partition == 0]

  # rearrange block 1
  block1.ra <- rearrange(block1, block2)
  X[, partition == 1] <- block1.ra

  return(X)
}
