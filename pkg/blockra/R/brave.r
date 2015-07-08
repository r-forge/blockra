#' Block ReArrangement Variance Equalizer
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
#' @seealso \code{\link{ra}} for the rearrangement algorithm
#' @seealso \code{\link{blockra}} for the block rearrangement algorithm
#'
#' @references \url{LINK TO RA PAPER}
#'
#' @author Kris Boudt, \email{kris.boudt@@vub.ac.be}
#' @author Steven Vanduffel, \email{steven.vanduffel@@vub.ac.be}
#' @author Kristof Verbeken, \email{kristof.verbeken@@vub.ac.be}
brave <- function(X, epsilon = 0.1, shuffle = TRUE, fix.first = TRUE, obj = var) {
  if (shuffle) 
    X <- shufflematrix(X, fix.first)
  obj.new <- obj(rowSums(X))
  obj.old <- 2 * obj.new
  obj.target <- epsilon * mean(apply(X, 2, obj))
  while ((obj.new > obj.target) ) {
    partition <- equalvar(X)
    X <- rearrangepartition(X, partition, fix.first)
    obj.old <- obj.new
    obj.new <- obj(rowSums(X))
    if( obj.new == obj.old ){
      partition <- sample(0:1, ncol(X), replace = TRUE)
      X <- rearrangepartition(X, partition, fix.first) 
      obj.old <- obj.new
      obj.new <- obj(rowSums(X))
    }
  }
  return(X)
}
