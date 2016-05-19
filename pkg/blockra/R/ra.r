#' Rearrangement Algorithm
#'
#' Function which performs the traditional RA
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
#' @seealso \code{\link{blockra}} for the block rearrangement algorithm
#' @seealso \code{\link{brave}} for the block rearrangement variance equalizer
#'
#' @references \url{LINK TO RA PAPER}
#'
#' @author Kris Boudt, \email{kris.boudt@@vub.ac.be}
#' @author Steven Vanduffel, \email{steven.vanduffel@@vub.ac.be}
#' @author Kristof Verbeken, \email{kristof.verbeken@@vub.ac.be}
ra <- function(X, f=var, shuffle=T, maxiter=1e3, stalliter=ncol(X), abs.tol=0, rel.tol=0, f.target=-Inf) {
  d<-ncol(X)
  if (shuffle) 
    X <- shufflematrix(X)
  fiter<-rep(NA,maxiter)
  niter <- 1
  fval <- f(rowSums(X))
  fprev<-Inf
  fiter[1]<-fval
  converged<-FALSE
  
  col<-1
  
  while ((fval>f.target)&&(niter<maxiter)&&(converged==F)) {
    niter <- niter+1
    
    current <- sort(X[,col]) # ascending
    other.sums <- rowSums(X[,-col])
    other.sums.rank <- rank(-other.sums, ties.method = "first") # descending
    X[, col] <- current[other.sums.rank]
    
    fval<-f(rowSums(X))
    fiter[niter]<-fval
    
    if(col==d){col<-1}else{col<-col+1 }
    
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