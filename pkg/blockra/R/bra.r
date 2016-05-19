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
#' @return list of rearranged matrix, rowsums (descending), iterations, current objective, iterations of objective, whether converged
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
#' 
bra <- function(X, f=var, shuffle=T, maxiter=1e3, stalliter=ncol(X), abs.tol=0, rel.tol=0, f.target=-Inf) {
  d<-ncol(X)
  if (shuffle) 
    X <- shufflematrix(X)
  fiter<-rep(NA,maxiter)
  niter <- 1
  fval <- f(rowSums(X))
  fprev<-Inf
  fiter[1]<-fval
  converged<-FALSE
  
  while ((fval>f.target)&&(niter<maxiter)&&(converged==F)) {
    niter <- niter+1
    partition <- sample(c(T,F), d, replace = TRUE)
    X <- rearrangepartition(X,partition)
    fval<-f(rowSums(X))
    fiter[niter]<-fval
    
    if (niter>stalliter)
      fprev<-fiter[niter-stalliter]
    converged<-(fprev-fval)<=max(abs.tol,rel.tol*abs(fval))
  }
  
  # sort according to decreasing rowsums
  S<-rowSums(X)
  ord<-order(S,decreasing=T)
  S<-S[ord]
  X<-X[ord,]
  # prepare output
  out<-{}
  out$X<-X
  out$S<-S
  out$niter<-niter
  out$fval<-fval
  fiter<-fiter[1:niter]
  out$fiter<-fiter
  out$converged<-converged
  return(out)
}