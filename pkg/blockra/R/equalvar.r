#' Find two blocks with equal variance
#'
#' The matrix is divided in two blocks with an equal variance using the greedy
#' implementation of the partition problem on the covariances between the columns
#' and the total row sums.
#'
#' The objective is to find two block with an variance that is as close as
#' possible to the variance of the other block without significantly impacting
#' the performance.
#'
#' The first step is to calculate the covariances between all columns and
#' the total row sums of the matrix. Using the greedy implementation of the
#' partition problem, we find two groups with variances that are as close as
#' possible to eachother.
#'
#' @param X numeric array or matrix
#'
#' @return numeric binary vector with the partition
#'
#' @export
#'
#' @keywords internal
#'
#' @references \url{Lhttp://en.wikipedia.org/wiki/Partition_problem#The_greedy_algorithm}
#'
#' @author Kris Boudt, \email{kris.boudt@@vub.ac.be}
#' @author Steven Vanduffel, \email{steven.vanduffel@@vub.ac.be}
#' @author Kristof Verbeken, \email{kristof.verbeken@@vub.ac.be}
equalvar <- function(X) {
  swept       <- X - tcrossprod(rep(1, nrow(X)), colMeans(X))
  covariances <- crossprod(swept, rowSums(swept)) / nrow(swept - 1)
  covariances <- sort(covariances, index.return = TRUE)
  partition   <- rep(0, length(covariances$x))

  block1 <- 0
  block2 <- 0

  for (element in 1 : length(covariances$x)) {
    value <- covariances$x[element]
    index <- covariances$ix[element]

    if ((value > 0 && block1 < block2) ||
        (value < 0 && block1 > block2)) {
      block1 <- block1 + value
      partition[index] <- 1
    } else {
      block2 <- block2 + value
    }
  }

  return(partition)
}
