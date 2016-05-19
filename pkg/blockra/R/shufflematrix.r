#' Shuffle matrix
#'
#' Permute randomly each column of the matrix, except for the first column
#'
#' @param X numeric array or matrix
#' @param fix.first don't change the order of the first column
#'
#' @return randomly permuted numeric matrix
#'
#' @export
#'
#' @seealso The \code{\link{ra}} for the rearrangement algorithm
#'
#' @seealso \code{\link{ra}} for the rearrangement algorithm
#' @seealso \code{\link{blockra}} for the block rearrangement algorithm
#'
#' @references \url{http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2634669}
#'
#' @author Kris Boudt, \email{kris.boudt@@vub.ac.be}
#' @author Edgars Jakobsons, \email{edgars.jakobsons@math.ethz.ch}
#' @author Steven Vanduffel, \email{steven.vanduffel@@vub.ac.be}
#' @author Kristof Verbeken, \email{kristof.verbeken@@vub.ac.be}
#' 
shufflematrix <- function(X) {
  return(apply(X, 2, sample))
}
