context('Checking core PhyloFactor functionality and summary functions')

set.seed(2)
data("FTmicrobiome",package = 'phylofactor')
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


### PhyloFactor ###
PF <- PhyloFactor(Data,tree,X,nfactors=2)
PF$bins
nms <- names(PF$bins)
names(PF$bins) <- NULL
test_that('PhyloFactor bins are not correct', expect_true( all.equal(unlist(PF$bins),unlist(Bins))))
names(PF$bins) <- nms
PF.par <- PhyloFactor(Data,tree,X,nfactors=2,ncores=2)

test_that('Parellelization works, but Serial and Parellel PhyloFactorizations are not equal',expect_true(all.equal(PF,PF.par)))

smry <- pf.summary(PF,Taxonomy,factor=1)

test_that('pf.summary incorrectly correctly summarizes split taxa',
          expect_true('k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Ruminococcaceae; g__; s__' %in% smry$TaxaSplit[[1]]$TaxaIDs & 'k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Ruminococcaceae; g__Faecalibacterium' %in% smry$TaxaSplit[[2]]$TaxaIDs))


options(warn=-1)
PF.stop <- PhyloFactor(Data,tree,X,stop.early = T,KS.Pthreshold = 0.05)
PF.stop$terminated <- FALSE
options(warn=0)

test_that('Default stop function works',expect_true(all.equal(PF.stop,PF)))
options(warn=-1)
PF.stop <- PhyloFactor(Data,tree,X,stop.early=F,KS.Pthreshold=0.05)
options(warn=0)

test_that('Stop Late works',expect_true(PF.stop$nfactors==3))

###### non-default methods: method=='max.var'
pf.var <- PhyloFactor(Data,tree,method='max.var')
test_that('method="max.var" works',expect_true(pf.var$nfactors==19))
pf.var <- PhyloFactor(Data,tree,method='max.var',ncores=2,nfactors=2)
test_that('method="max.var" works in parallel',expect_true(pf.var$nfactors==2))

###### method='gam'
X <- data.frame('a'=rnorm(10))
pf.gam <- PhyloFactor(Data,tree,X=X,method='gam',frmla=Data~s(a),nfactors=2,sp=1)
test_that('method="gam" works',expect_true(pf.gam$nfactors==2))
pf.gam <- PhyloFactor(Data,tree,X=X,method='gam',frmla=Data~s(a),nfactors=2,ncores=2,sp=1)
test_that('method="gam" works in parallel',expect_true(pf.gam$nfactors==2))