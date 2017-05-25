context('Checking PhyloFactor multiple regression functionality')

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
a <- X
b <- rnorm(length(X))
A <- data.frame(a,b)
frmla <- Data ~ a + b^2

PF <- PhyloFactor(Data,tree,A,frmla,nfactors=2)
PF.par <- PhyloFactor(Data,tree,A,frmla,nfactors=2,ncores=2)

test_that('Multiple Regression Data ~ a + b^2 works, serial & parallel',expect_true(all.equal(PF,PF.par)))


s <- pf.summary(PF.par,Taxonomy,factor=2)

test_that('pf.summary works for multiple regression',expect_true(TRUE))

td <- pf.tidy(s)

test_that('pf.tidy taxonomy is good for multiple regression',expect_true(all(grepl('Clostridiales',td$`group2, Paraphyletic`))))
