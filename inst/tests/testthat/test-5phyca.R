context('Checking PhyCA and phyca.tidy')

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



############### PhyCA Testing #############
invisible(capture.output(phca <- PhyCA(Data,tree,ncomponents=2)))
test_that('PhyCA captures correct bins',
          expect_true(all(all.equal(phca$bins[[3]],sigClades[[1]]),all.equal(phca$bins[[2]],sigClades[[2]]))))
invisible(capture.output(suppressWarnings(phca <- PhyCA(Data[sample(nrow(Data)),],tree,ncomponents=2))))
test_that('PhyCA works with scrambled data',
          expect_true(all(all.equal(phca$bins[[3]],sigClades[[1]]),all.equal(phca$bins[[2]],sigClades[[2]]))))
invisible(capture.output(phca.par <- PhyCA(Data,tree,ncomponents = 2,ncores=2)))

test_that('PhyCA serial & parellel work and are equal',expect_true(all.equal(phca,phca.par)))

test_that('phyca.tidy works', {
  capture.output( {td <- phyca.tidy(phca,Taxonomy)
  td <- phyca.tidy(phca,Taxonomy,taxa.split = T)
  td <- phyca.tidy(phca,Taxonomy,taxa.split = T,common.name = F)
  td <- phyca.tidy(phca,Taxonomy,taxa.split = T,uniques = F)})})