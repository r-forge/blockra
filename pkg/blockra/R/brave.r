brave <- function(X, epsilon = 0.1) {
  variance.new <- var(rowSums(X)) 
  variance.old <- 2 * variance.new
  iterations <- 0
  time.start <- proc.time()
  
  target <- epsilon * mean(CalculateVariances(X))
  
  while (variance.new > target) {
    partition <- EqualBlocks(X)
    # p <<- profr(partition <- EqualBlocks(X))    
    
    X <- RearrangeSelection(X, partition)
    
    iterations <- iterations + 1
    variance.old <- variance.new
    variance.new <- var(rowSums(X))
    
  }
  
  time.elapsed <- proc.time() - time.start
  
  cat("\r\n---------")
  cat("\r\nBRAVE (new)")
  cat("\r\n---------")
  cat("\r\nIterations: ", iterations)
  cat("\r\nTime: ", time.elapsed[3])
  cat("\r\n---------")
  
  return(iterations)
}

EqualBlocks <- function(X) {  
  # Create X of 2 rows
  # First row has initial order of X
  # Second row has covariances associated to that order
#   covariances <- X(1 : ncol(X), 
#                         ncol = ncol(X), 
#                         nrow = 10,
#                         byrow = TRUE)
  
  # Calculate covariances with total row sums for each column
  swept <- X - tcrossprod(rep(1, nrow(X)), colMeans(X))
  # covariances[2, ] <- crossprod(swept, rowSums(swept)) /(nrow(swept-1))
  covariances <- crossprod(swept, rowSums(swept)) /(nrow(swept-1))
  
  partition <- FindPartition(covariances)
  
  return(partition)
}


FindPartition <- function(covariances) {

  
  block1 <- 0
  block2 <- 0
  covariances <- sort(covariances, index.return = TRUE)
  partition <- rep(0, length(covariances$x))
  
  for (element in 1 : length(covariances$x)) {
    value <- covariances$x[element]
    
    if (value > 0) {
      if (block1 < block2) {
        block1 <- block1 + value
        partition[covariances$ix[element]] <- 1
      } else {
        block2 <- block2 + value
      }      
    } else {
      if (block1 > block2) {
        block1 <- block1 + value
        partition[covariances$ix[element]] <- 1
      } else {
        block2 <- block2 + value
      }      
    }
  }
  
  return(partition)
}




CalculateVariances <- function(X) {
  variances <- numeric(ncol(X))
  
  for (column in 1 : ncol(X)) {
    currentCol <- X[ , column]
    currentVar <- var(currentCol)
    variances[column] <- currentVar
  }
  
  return(variances)
}