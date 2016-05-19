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
#' @references \url{http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2634669}
#'
#' @author Kris Boudt, \email{kris.boudt@@vub.ac.be}
#' @author Edgars Jakobsons, \email{edgars.jakobsons@math.ethz.ch}
#' @author Steven Vanduffel, \email{steven.vanduffel@@vub.ac.be}
#' @author Kristof Verbeken, \email{kristof.verbeken@@vub.ac.be}
#' 
#' to test: d<-10; X<-matrix(runif(d^2,max=10),ncol=d); partition<-sample(c(T,F),d,replace=T); print(equalvar(X,partition,"greedy"));print(equalvar(X,partition,"KK"));
equalvar <- function(X,partition.prev=logical(ncol(X)),method="KK") {
  
  d<-ncol(X)
  if (!((length(partition.prev)==d)&&is.logical(partition.prev)&&is.vector(partition.prev)))
    stop("previous partition is not a logical vector of length compatible with X")
  
  equalize<-switch(method,greedy=greedy,KK=KK,rapartition=rapartition,stop("Partitioning method not available!"))
  partition<-logical(d)
  
  demean <- X - tcrossprod(rep(1, nrow(X)), colMeans(X))
  covars <- crossprod(demean, rowSums(demean)) / (nrow(demean)-1)
  
  # make sure the first block is the smaller one
  if (sum(partition.prev) > d/2) {
    partition.prev=!partition.prev
  }  
  
  partA<-equalize(covars[partition.prev])
  partB<-equalize(c(partA$diff,covars[!partition.prev])) # we add the difference of previously assigned elements as the first number
  
  if (all(!partition.prev)){ # in case the (smaller) forst block is empty, take the partition of the larger block
    partition <- partB$partition
  }
  else{
    
    if (partB$partition[1]){  # combine the partitions: if the difference belongs to the F block, then invert partA$partition
      partition[partition.prev] <- partA$partition
    }
    else{
      partition[partition.prev] <- !partA$partition
    }
    partition[!partition.prev] <- partB$partition[-1]
    
  }
  #cat(method,t(covars)%*%(2*partition-1),"\n")
  
  if(all(partition.prev==partition)){warning("Same partition as previously!")}
  
  return(partition)
}
