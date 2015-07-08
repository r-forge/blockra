#' Block rearrangement algorithm
#'
#' Function which performs the block ra of Bernard and McLeish
#'
#' @param X numeric array or matrix
#' @param epsilon target variance of row sums is epsilon multiplied by the mean of the matrix variances
#' @param shuffle randomly permute each column of the matrix before rearrangement
#' @param fix.first don't change the order of the first column
#' @param obj objective function that is minimized, default is variance
#'
#' @return numeric matrix with a minimal row sum variance
#'
#' @export
#'
#' @seealso The \code{\link{ra}} for the rearrangement algorithm
#'
#' @references \url{LINK TO RA PAPER}
#'
#' @author Kris Boudt, \email{kris.boudt@@vub.ac.be}
#' @author Steven Vanduffel, \email{steven.vanduffel@@vub.ac.be}
#' @author Kristof Verbeken, \email{kristof.verbeken@@vub.ac.be}
bra <- function(X, epsilon = 0.1, shuffle = TRUE, fix.first = TRUE, obj = var) {
  itermax = 1e3
  if (shuffle) 
    X <- shufflematrix(X, fix.first)
  obj.new <- obj(rowSums(X))
  obj.old <- 2 * obj.new
  obj.target <- epsilon * mean(apply(X, 2, obj))
  citer <- 0
  while ((obj.new > obj.target ) & (citer<itermax) ) {
    partition <- sample(0:1, ncol(X), replace = TRUE)
    X <- rearrangepartition(X, partition, fix.first)
    obj.old <- obj.new
    obj.new <- obj(rowSums(X))
  }
  return(X)
}
