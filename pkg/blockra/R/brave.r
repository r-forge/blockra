#' Block ReArrangement Variance Equalizer
#'
#' @param f function of the rowsums to be minimized
#' @param shuffle randomly permute each column of the matrix before rearrangement
#' @param maxiter number of maximum iterations
#' @param stalliter convergence if no improvement after stalliter iteration
#' @param abs.tol abs convergence crit
#' @param rel.tol relative convergence criterion
#' @param f.target value of the objective function at which absolute converge is reached
#' @param version version of the brave use
#' @param method partitioning method used
#' @param nochange.tol switch to random block condition
#'
#' @return list of rearranged matrix, rowsums (descending), iterations, current objective, iterations of objective, whether converged
#'
#' @export
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
brave <- function(X, f=var, shuffle=F, maxiter=1e3, stalliter=maxiter, abs.tol=0, rel.tol=0, f.target=-Inf, version=1,nochange.tol=0,method="KK") {
  d<-ncol(X)
  if (shuffle) 
    X <- shufflematrix(X)
  fiter<-rep(NA,maxiter)
  blocksizeiter <-rep(NA,maxiter) #kb
  niter <- 1
  fval <- f(rowSums(X))
  fprev<-Inf
  fiter[1]<-fval
  partition.rand<-logical(maxiter)
  partition.prev<-logical(d)
  RA.col<-0
  RA.switch<-FALSE
  converged<-FALSE
  varprev<-Inf
  
  while ((fval>f.target)&&(niter<maxiter)&&(converged==F)) {
    niter <- niter+1
    partition <- equalvar(X,partition.prev,method)
    blocksizeiter[niter] <- min(sum(partition),d-sum(partition)) #kb
    X <- rearrangepartition(X,partition)
    fval<-f(rowSums(X))
    varcurr<-var(rowSums(X))
    fiter[niter]<-fval
    # if the objective is unchanged, select a random partition
    if((varprev-varcurr<=nochange.tol*varcurr)&&(niter<maxiter)){
      niter <- niter+1
      partition.rand[niter]<-TRUE
      switch(version,
             { #1 random block
               partition <- sample(c(T,F), d, replace = TRUE)
               blocksizeiter[niter] <- min(sum(partition),d-sum(partition)) #kb
             },
             { #2 random single column
               partition<-logical(d)
               partition[sample(1:d,1)]<-TRUE
               blocksizeiter[niter] <- 1
             },
             { #3 next column in sequence
               RA.col<-RA.col+1
               if(RA.col>d){RA.col<-1}
               partition<-logical(d)
               partition[RA.col]<-TRUE
               blocksizeiter[niter] <- 1
             },
             { #4 permanently switch to RA
               RA.switch<-TRUE
               partition<-logical(d)
             }
      )
      if (!RA.switch){
      X <- rearrangepartition(X, partition) 
      fval<-f(rowSums(X))
      fiter[niter]<-fval
      }
    }
    varprev<-var(rowSums(X))
    partition.prev<-partition
    if (niter>stalliter) {fprev<-fiter[niter-stalliter]}
    converged<-(RA.switch||((fprev-fval)<=max(abs.tol,rel.tol*abs(fval))))
  }
  
  if (RA.switch){
    ra.out<-ra(X, f=f, shuffle=shuffle, maxiter=maxiter-niter+1, stalliter=stalliter, abs.tol=abs.tol, rel.tol=rel.tol, f.target=f.target)
    niter.tot<-niter-1+ra.out$niter
    fiter[niter:niter.tot]<-ra.out$fiter
    partition.rand[niter:niter.tot]<-TRUE
    blocksizeiter[niter:niter.tot]<-1
    niter<-niter.tot
    fval<-ra.out$fval
    converged<-ra.out$converged
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
  out$fiter<-fiter[1:niter]
  out$rand<-partition.rand[1:niter]
  out$converged<-converged
  out$blocksizeiter <- blocksizeiter[1:niter] #kb
  return(out)
}
