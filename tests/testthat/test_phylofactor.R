set.seed(2)
data("FTmicrobiome")
tree <- FTmicrobiome$tree
Taxonomy <- FTmicrobiome$taxonomy
tree <- ape::drop.tip(tree,setdiff(tree$tip.label,sample(tree$tip.label,20)))

Taxonomy <- Taxonomy[match(tree$tip.label,Taxonomy[,1]),]
X <- as.factor(c(rep(0,5),rep(1,5)))

### Simulate data ###
Factornodes <- c(37,27)
sigClades <- phangorn::Descendants(tree,Factornodes,type='tips')

Data <- matrix(rlnorm(20*10,meanlog = 8,sdlog = .5),nrow=20)
rownames(Data) <- tree$tip.label
colnames(Data) <- X
Data[sigClades[[1]],X==0] <- Data[sigClades[[1]],X==0]*8
Data[sigClades[[2]],X==1] <- Data[sigClades[[2]],X==1]*9
Data <- t(clo(t(Data)))
Bins <- bins(G=sigClades,set=1:20)


### PhylOFactor ###
PF <- PhyloFactor(Data,tree,X,nfactors=2)
PF$bins

test_that('PhyloFactor bins are correct', expect_true( all(PF$bins %in% Bins )))
            
PF.par <- PhyloFactor(Data,tree,X,nfactors=2,ncores=2)

test_that('Serial and Parellel are equal',expect_true(all.equal(PF,PF.par)))

smry <- pf.summary(PF,Taxonomy,factor=1)

test_that('pf.summary correctly summarizes split taxa',expect_true('k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Ruminococcaceae; g__; s__' %in% smry$TaxaSplit[[1]]$TaxaIDs & 'k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Ruminococcaceae; g__Faecalibacterium;' %in% smry$TaxaSplit[[2]]$TaxaIDs))


