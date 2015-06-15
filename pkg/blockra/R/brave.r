#' Block ReArrangement Variance Equalizer
#'
#' @param X numeric array or matrix
#' @param epsilon target variance of row sums is epsilon multiplied by the mean of the matrix variances
#' @param shuffle randomly permute each column of the matrix before rearrangement
#'
#' @return numeric matrix with a minimal row sum variance
#'
#' @export
#'
#' @seealso \code{\link{ra}} for the rearrangement algorithm
#' @seealso \code{\link{blockra}} for the block rearrangement algorithm
#'
#' @references \url{LINK TO RA PAPER}
#'
#' @author Kris Boudt, \email{kris.boudt@@vub.ac.be}
#' @author Steven Vanduffel, \email{steven.vanduffel@@vub.ac.be}
#' @author Kristof Verbeken, \email{kristof.verbeken@@vub.ac.be}
brave <- function(X, epsilon = 0.1, shuffle = TRUE, fix.first = TRUE) {
  if (shuffle) X <- shufflematrix(X, fix.first)

  var.new   <- var(rowSums(X))
  var.old   <- 2 * var.new
  target    <- epsilon * mean(apply(X, 2, var))

  while ((var.new > target) && (var.new < var.old)) {
    partition <- equalvar(X)
    X         <- rearrangepartition(X, partition, fix.first)
    var.old   <- var.new
    var.new   <- var(rowSums(X))
  }

  return(X)
}
