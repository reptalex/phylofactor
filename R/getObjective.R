#' objective function for \code{\link{gpf}}
#' @export
#' @param grp list containing two disjoint lists of species, such as thouse output from \code{\link{getGroups}}
#' @param tree phylo class object
#' @param DF data table containing Species, observations (N), and sample
#' @param X meta-data containing variables in formula

getObjective <- function(grp,tree,DF,X,SampleFrame,size=1,frmla){
  ###################### define & merge factor ############
  # ff <- factorFrame(grp,tree$tip.label)
  # r <- unique(ff[G=='R']$Size)
  # s <- unique(ff[G=='S']$Size)
  # data.table::setkey(DF,Species)
  # DF2 <- DF[ff]
  # # DF2 <- data.table:::`[.data.table`(DF,ff)
  # setkey(DF2,G,sample)
  # 
  # ###################### marginal-stable aggregatoin ######
  # ## insert code to add zeros here:
  # DF2 <- DF2[,list(Successes=sum(N),
  #                  Failures=size*(unique(Size))-sum(N),
  #                  Size=unique(Size)),by=list(G,sample)]
  # setkey(DF2,sample)
  # 
  # ##### Inputting zeros
  #   setkey(SampleFrame,sample,G)
  #   setkey(DF2,sample,G)
  #   DF2 <- DF2[SampleFrame]
  #   ix <- is.na(DF2$Successes)
  #   DF2[ix & G=='R']$Size <- r
  #   DF2[ix & G=='S']$Size <- s
  #   DF2[ix & G=='R']$Failures <- size*r
  #   DF2[ix & G=='S']$Failures <- size*s
  #   DF2[ix]$Successes <- 0
  # 
  # setkey(X,sample)
  # DF2 <- DF2[X]
  DF2 <- mAggregation(DF,grp,tree,X,size,SampleFrame)
  r <- unique(DF2[G=='R',Size])
  s <- unique(DF2[G=='S',Size])
  ss <- tryCatch(anova(glm(frmla,data=DF2,
                           family=binomial),
                       test='Chisq'),error=function(e) 0)
  

  if ('anova' %in% class(ss)){
    nms <- grepl(':G',rownames(ss))
    objective <- sum(ss$Deviance[nms])
  } else {
    warning(paste('Model failed for a group. Objective set to 0 for phyloGroups of sizes',r,'and',s,sep=' '))
    objective <- 0
  }
  return(objective)
}