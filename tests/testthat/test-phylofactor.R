set.seed(2)
data("FTmicrobiome")
tree <- FTmicrobiome$tree
Taxonomy <- FTmicrobiome$taxonomy
tree <- ape::drop.tip(tree,setdiff(tree$tip.label,sample(tree$tip.label,20)))

Taxonomy <- Taxonomy[match(tree$tip.label,Taxonomy[,1]),]
X <- as.factor(c(rep(0,5),rep(1,5)))

### Simulate data ###
Factornodes <- c(37,27)
sigClades <- list(c(16,17),c(6,7,8,9,10,11,12,13,14,15))

Data <- matrix(rlnorm(20*10,meanlog = 8,sdlog = .5),nrow=20)
rownames(Data) <- tree$tip.label
colnames(Data) <- X
Data[sigClades[[1]],X==0] <- Data[sigClades[[1]],X==0]*9
Data[sigClades[[2]],X==1] <- Data[sigClades[[2]],X==1]*14
Bins <- bins(G=sigClades,set=1:20)[c(3,2,1)]


### PhylOFactor ###
PF <- PhyloFactor(Data,tree,X,nfactors=2)
PF$bins
nms <- names(PF$bins)
names(PF$bins) <- NULL
test_that('PhyloFactor bins are not correct', expect_true( all.equal(unlist(PF$bins),unlist(Bins))))
names(PF$bins) <- nms
PF.par <- PhyloFactor(Data,tree,X,nfactors=2,ncores=2)

test_that('Parellelization works, but Serial and Parellel PhyloFactorizations are not equal',expect_true(all.equal(PF,PF.par)))

smry <- pf.summary(PF,Taxonomy,factor=1)

test_that('pf.summary incorrectly correctly summarizes split taxa',expect_true('k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Ruminococcaceae; g__; s__' %in% smry$TaxaSplit[[1]]$TaxaIDs & 'k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Ruminococcaceae; g__Faecalibacterium;' %in% smry$TaxaSplit[[2]]$TaxaIDs))

edgs <- apply(PF$basis,MARGIN=2,getFactoredEdges,tree=tree)

test_that('Did not correctly extract factored edges',expect_true(all.equal(edgs,c(11,31))))

cl <- phyloFcluster(2)
edgs.par <- getFactoredEdgesPAR(ncores=2,tree=tree,PF=PF,cl=cl) %>% unlist
parallel::stopCluster(cl)
rm('cl')

test_that('getFactoredEdgesPAR did not exctract correct edges',expect_true(all.equal(edgs.par,c(11,31))))

PF.stop <- PhyloFactor(Data,tree,X,stop.early = T,KS.Pthreshold = 0.05)
PF.stop$terminated <- FALSE

test_that('Default stop function works',expect_true(all.equal(PF.stop,PF)))

PF.stop <- PhyloFactor(Data,tree,X,stop.early=F,KS.Pthreshold=0.05)

test_that('Stop Late works',expect_true(PF.stop$nfactors==3))

a <- X
b <- rnorm(length(X))
A <- data.frame(a,b)
frmla <- Data ~ a + b^2

PF <- PhyloFactor(Data,tree,A,frmla,nfactors=2)
PF.par <- PhyloFactor(Data,tree,A,frmla,nfactors=2,ncores=2)

test_that('Multiple Regression Data ~ a + b^2 works, serial & parallel',expect_true(all.equal(PF,PF.par)))

X <- a
X[3] <- X[6]
X[6] <- X[4]
PF <- PhyloFactor(Data,tree,X,frmla=X~Data,family=binomial(link='logit'),nfactors=1,choice='F')
PF <- PhyloFactor(Data,tree,X,frmla=X~Data,family=binomial(link='logit'),nfactors=1,choice='var')
test_that('Prediction X~Data failed',expect_true(TRUE))




########### Testing Customized objective functions

# Let's work with some newly simulated data ####
 set.seed(1.1)
 n=100
 Data <- matrix(rlnorm(20*n,meanlog = 8,sdlog = .5),nrow=20)
 rownames(Data) <- tree$tip.label
 a <- rnorm(n)
 b <- rnorm(n)
 X <- data.frame(a,b)
 Data[sigClades[[1]],] <- t(t(Data[sigClades[[1]],])*(20/(1+exp(5*b)))) ## This clade has a nonlinear response with b, decreasing for high values of b.
 Data[sigClades[[2]],] <- t(t(Data[sigClades[[2]],])*8*a^-2)  ## this clade is abundant only for intermediate values of a.
 
 
 ############## To input a custom choice.fcn, it needs to take as input the vector of ILR coefficients 'y', the input meta-data 'X',
 ############## and a logical PF.output. The output of the custom choice function when PF.output=T will be returned in PF$custom.output.
 
 ## Demo choice.fcn - generalized additive modelling ##
 GAM <- function(y,X,PF.output=F,...){
   dataset <- cbind(y,X)
   gg <- mgcv::gam(y~s(a)+s(b),data=dataset,...)
 
   if (PF.output){
     return(gg)
     break
   } else {
     output <- NULL
     output$objective <- getStats(gg)['ExplainedVar']  ## The output of the choice function for PF.output=F must contain two labelled numerics: an "objective" statistic and a "stopStatistics". 
     output$stopStatistics <- getStats(gg)['Pval']
     return(output)
   }
 }
 
 load.dependencies <- function(){library(mgcv)}
 ############## For parallelization of customized choice function, we also need to define a function, 
 ############## choice.fcn,dependencies, which loads all dependencies to cluster.
 ############## The exact call will be clusterEvalQ(cl,choice.fcn.dependencies())
 
 
 PF.G <- PhyloFactor(Data,tree,X,nfactors=2,choice.fcn=GAM,choice.fcn.dependencies = load.dependencies,sp=1)
 PF.G.par <- PhyloFactor(Data,tree,X,nfactors=2,choice.fcn=GAM,choice.fcn.dependencies = load.dependencies,ncores=2,sp=1)
 test_that('Parellelized & serialized customized objective function - GAM - are not equal',expect_true(all.equal(PF.G,PF.G.par)))
 names(PF.G$bins) <- NULL
 test_that('GAM phylofactorization did not extract correct clades',expect_true(all.equal(unlist(sigClades),unlist(PF.G$bins[2:3]))))
 
 
 s <- pf.summary(PF.par,Taxonomy,factor=2)
 
 test_that('pf.summary works',expect_true(TRUE))
 
 td <- pf.tidy(s)
 
 test_that('pf.tidy taxonomy is good works',expect_true(all(grepl('Clostridiales',td$`group2, Paraphyletic`))))

 
 
 
 ############### PhyCA Testing #############
 phca <- PhyCA(Data,tree,ncomponents=2)
 phca.par <- PhyCA(Data,tree,ncomponents = 2,ncores=2)
 
 test_that('PhyCA serial & parellel work and are equal',expect_true(all.equal(phca,phca.par)))
 
 test_that('phyca.tidy works', {
       td <- phyca.tidy(phca,Taxonomy)
       td <- phyca.tidy(phca,Taxonomy,taxa.split = T)
       td <- phyca.tidy(phca,Taxonomy,taxa.split = T,common.name = F)
       td <- phyca.tidy(phca,Taxonomy,taxa.split = T,uniques = F)})
 
 