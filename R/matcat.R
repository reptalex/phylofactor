matcat <- function(A,B){
  if (typeof(B) != 'double' || typeof(A) != 'double'){stop('input must be two matrices')}
  
  p1 <- dim(A)[2]
  p2 <- dim(B)[2]
  
  m1 <- dim(A)[1]
  m2 <- dim(B)[1]
  
  C <- A
  #First, we'll bind the overlaps
  dum <- matrix(0,nrow=m1,ncol=p2)
  ix1 <- rownames(A)
  ix2 <- rownames(B)
  overlap <- intersect(ix1,ix2)
  rownames(dum) <- ix1
  colnames(dum) <- colnames(B)
  ind1 <- match(ix2[ix2 %in% overlap],ix1) 
  ind2 <- match(ix1[ix1 %in% overlap],ix2)
  dum[ind1,] = B[ind2,]
  C <- cbind(C,dum)
  
  #Then, we'll bind the non-overlaps
  D <- B[!(ix2 %in% ix1),]
  dum <- matrix(0,nrow=dim(D)[1],ncol=p1)
  colnames(dum)<- colnames(A)
  rownames(dum) <- rownames(D)
  D <- cbind(dum,D)
  C <- rbind(C,D)
  return(C)
}