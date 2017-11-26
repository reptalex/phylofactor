context('Checking getFactoredEdges and getFactoredEdgesPAR')

set.seed(2)
data("FTmicrobiome",package='phylofactor')
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

PF <- PhyloFactor(Data,tree,X,nfactors=2)

edgs <- apply(PF$basis,MARGIN=2,getFactoredEdges,tree=tree)

test_that('Did not correctly extract factored edges from phylofactor basis',expect_true(all.equal(edgs,c(11,31))))

edgs.par <- getFactoredEdgesPAR(ncores=2,tree=tree,PF=PF) %>% unlist

test_that('getFactoredEdgesPAR can extract edges from phylofactor object',expect_true(exists('edgs.par')))
test_that('getFactoredEdgesPAR extracts correct edges from phylofactor object',expect_true(all.equal(edgs.par,c(11,31))))
