BRA2 <- function(X, epsilon = 0.1) {
  
  variance.new <- var(rowSums(X)) 
  variance.old <- 2 * variance.new
  iterations <- 0
  time.start <- proc.time()
  
  target <- epsilon * mean(CalculateVariances(X))
  
  while ((variance.new > target) && (variance.new < variance.old)) {
    
    
    partition <- sample(0 : 1, ncol(X), replace = TRUE)
    X <- RearrangeSelection(X, partition)

    iterations <- iterations + 1
    variance.old <- variance.new
    variance.new <- var(rowSums(X))
  }
  
  time.elapsed <- proc.time() - time.start
  
  cat("\r\n---------")
  cat("\r\nBRA2")
  cat("\r\n---------")
  cat("\r\nIterations: ", iterations)
  cat("\r\nTime: ", time.elapsed[3])
  cat("\r\n---------")
  
  return(X)
}


RearrangeSelection <- function(X, selection) {
  # Rearrange the X based on a selection
  #
  # Args:
  #   X: X to be rearranged
  #   selection: Vector with selection 
  #
  # Returns:
  #   The rearranged X
  
  # Get submatrices based on selection
  X.first <- X[ , selection > 0]
  X.second <- X[ , selection <= 0]    
  

  X <- Rearrange(X.first, X.second)

  
  return(X)
}

Rearrange <- function(block1, block2) {
  # Inversely sort a small X against a large X  
  #
  # Args:
  #   block1: small partition
  #   block2: large partition
  #
  # Returns:
  #   The rearranged X
  
  
  # Make sure smallp is smaller than largep
  if(ncol(block2) < ncol(block1)) {
    Rearrange(block2, block1)
  }
  
  rank.block2 <- rank(RowSums(block2))
  rank.target <- (nrow(block1) + 1) - rank.block2
  index.initial <- sort(RowSums(block1), index.return = TRUE)$ix
  
  block1 <- block1[index.initial[rank.target], ]
  
  # Bind both partition together
  rearranged <- cbind(block1, block2)
  
  return(rearranged)
}

RowSums <- function(X) {
  if (is.X(X))
    return(rowSums(X))
  
  return(X)
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