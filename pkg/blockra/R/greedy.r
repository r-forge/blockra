# To test:
# nums<-round(runif(10,min=-100,max=100)); greed<-greedy(nums); kk<-KK(nums); print(c(greed,kk)); print(c(nums%*%(2*greed$partition-1),nums%*%(2*kk$partition-1)))

greedy <- function(X) {
  
  n<-length(X)
  partition<-rep(NA,n) # i-th entry contains the block index {T,F} for i-th input number
  sumA<-0
  sumB<-0

  if (n==0){
    out<-{}
    out$diff<-numeric(0)
    out$partition<-logical(0)
    return(out)
  }
  
  Sidx <- sort(abs(X), decreasing = TRUE,index=TRUE)$ix
  S<-X[Sidx]

  for (i in 1:n){
    if (xor(S[i] > 0,sumA>sumB)){
      sumA<-sumA+S[i]
      partition[Sidx[i]]<-TRUE
    }
    else{
      sumB<-sumB+S[i]
      partition[Sidx[i]]<-FALSE
    }
    #print(c(sumA,sumB))
  }
  if (sumA<sumB){
    partition<-!partition
  }
  
  # prepare output
  out<-{}
  out$diff <- abs(sumA-sumB)
  out$partition <- partition
  
  return(out)
}