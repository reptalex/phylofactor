################# phyloreg ############################
phyloreg <- function(y,x,frmla,choice,...){
  #input: dependent variable y, independent variable x, formula 'frmla' and 'choice' to be used to
  #compare different phylogenetic partitions
  ##performs individual regression.
  M <- data.frame(y,x)
  names(M) <- c('Data','X')
  return(glm(frmla,data = M,...))
}################# phyloreg ############################
phyloreg <- function(y,x,frmla,choice,...){
  #input: dependent variable y, independent variable x, formula 'frmla' and 'choice' to be used to
  #compare different phylogenetic partitions
  ##performs individual regression.
  M <- data.frame(y,x)
  names(M) <- c('Data','X')
  return(glm(frmla,data = M,...))
}