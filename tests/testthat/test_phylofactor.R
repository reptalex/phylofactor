test_phylofactor <- function(){
  # data("FTmicrobiome")
  # tree <- FTmicrobiome$tree
  # Taxonomy <- FTmicrobiome$taxonomy
  # tree <- ape::drop.tip(tree,setdiff(tree$tip.label,sample(tree$tip.label,20)))
  # 
  # Taxonomy <- Taxonomy[match(tree$tip.label,Taxonomy[,1]),]
  # X <- as.factor(c(rep(0,5),rep(1,5)))
  # 
  # ### Simulate data ###
  # Factornodes <- c(37,27)
  # Factoredges <- sapply(Factornodes,FUN=function(n,tree) which(tree$edge[,2]==n),tree=tree)
  # edgelabels(c('PF 1','PF 2'),edge=Factoredges,cex=2,bg='red')
  # sigClades <- phangorn::Descendants(tree,Factornodes,type='tips')
  # 
  # Data <- matrix(rlnorm(20*10,meanlog = 8,sdlog = .5),nrow=20)
  # rownames(Data) <- tree$tip.label
  # colnames(Data) <- X
  # Data[sigClades[[1]],X==0] <- Data[sigClades[[1]],X==0]*8
  # Data[sigClades[[2]],X==1] <- Data[sigClades[[2]],X==1]*9
  # Data <- t(clo(t(Data)))
  # Bins <- bins(G=sigClades,set=1:20)
  # 
  # 
  # ### PhylOFactor ###
  # PF <- PhyloFactor(Data,tree,X,nfactors=2)
  # PF$bins
  # check_tr all(PF$bins %in% Bins)
  # all.equal(PF,PhyloFactor(Data,tree,X,nfactors=2,ncores=2))
   expect_that(1+1,equals(3))
}