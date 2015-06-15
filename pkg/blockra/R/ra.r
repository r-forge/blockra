#' Rearrangement Algorithm
#'
#' Function which performs the traditional RA
#'
#' @param X numeric array or matrix
#' @param epsilon target variance of row sums is epsilon multiplied by the mean of the matrix variances
#' @param shuffle randomly permute each column of the matrix before rearrangement
#'
#' @return numeric matrix with a minimal row sum variance
#'
#' @export
#'
#' @seealso \code{\link{blockra}} for the block rearrangement algorithm
#' @seealso \code{\link{brave}} for the block rearrangement variance equalizer
#'
#' @references \url{LINK TO RA PAPER}
#'
#' @author Kris Boudt, \email{kris.boudt@@vub.ac.be}
#' @author Steven Vanduffel, \email{steven.vanduffel@@vub.ac.be}
#' @author Kristof Verbeken, \email{kristof.verbeken@@vub.ac.be}
ra <- function(X, epsilon = 0.1, shuffle = TRUE, fix.first = TRUE) {
  if (shuffle) X <- shufflematrix(X, fix.first)

  var.new   <- var(rowSums(X))
  var.old   <- 2 * var.new
  target    <- epsilon * mean(apply(X, 2, var))

  while ((var.new > target) && (var.new < var.old)) {

    for (col in 1 : ncol(X)) {
      current         <- sort(X[, col])
      other.sums      <- rowSums(X[ , -col])
      other.sums.rank <- rank(-other.sums, ties.method = "random")
      X[, col]   <- current[other.sums.rank]
    }

    var.old <- var.new
    var.new <- var(rowSums(X))
  }

  return(X)
}
