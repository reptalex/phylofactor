#' Function for aligning rows and concatenating two matrices.
#' @export
#' @param A matrix
#' @param B matrix
#' @return Matrix where upper-left block is matrix A, upper-right block is the abundances in B of parts in A. Lower left block is abundances in A of parts in B not found in A, and lower-right block is abundances in B of parts not found in A.
#' @examples
#' A <- matrix(rnorm(12),ncol=4)
#' rownames(A) <- 1:3
#' B <- matrix(rnorm(36),ncol=6)
#' rownames(B) <- sample(6)
#' matcat(A,B)

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
