#
# Rearrangement Algorithm
# 

RA <- function(X, epsilon = 0.1) {
  # Apply Rearrangement Algorithm until variance of row sums
  # is smaller than epsilon
  #
  # Args:
  #   X: Original X to be rearranged
  #   epsilon: Error margin
  #
  # Returns:
  #   The rearranged X
  
  variance.new <- var(rowSums(X)) 
  variance.old <- 2 * variance.new
  iterations <- 0
  time.start <- proc.time()
  
  target <- epsilon * mean(CalculateVariances(X))
  
  while (variance.new > target) {
    for (col in 1 : ncol(X)) {
      # Store current column
      current <- X[, col]
      
      # Sort the current column
      current <- sort(current)
      
      # Take sums of other columbs
      other.sums <- rowSums(X[ , -col])
      
      # Rank other sums
      other.sums.rank <- rank(-other.sums, ties.method = "random")
      
      # Save current col
      X[, col] <- current[other.sums.rank]
    }
    
    
    
    iterations <- iterations + 1
    variance.old <- variance.new
    variance.new <- var(rowSums(X))
    
    print(iterations)
    print(variance.old)
    print(variance.new)
    
  }
  
  time.elapsed <- proc.time() - time.start
  
  cat("\r\n---------")
  cat("\r\nRA")
  cat("\r\n---------")
  cat("\r\nIterations: ", iterations * ncol(X))
  cat("\r\nTime: ", time.elapsed[3])
  cat("\r\n---------")
  
  return(X)
}

