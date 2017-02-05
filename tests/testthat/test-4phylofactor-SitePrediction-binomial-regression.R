context('Checking PhyloFactor meta-data prediction via logit regression')

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
X[3] <- X[6]
X[6] <- X[4]
PF <- PhyloFactor(Data,tree,X,frmla=X~Data,family=binomial(link='logit'),nfactors=1,choice='F')
PF <- PhyloFactor(Data,tree,X,frmla=X~Data,family=binomial(link='logit'),nfactors=1,choice='var')
test_that('Prediction X~Data failed',expect_true(TRUE))

s <- pf.summary(PF,Taxonomy,1)
test_that('pf.summary works with X~Data regression',expect_true(exists('s')))
td <- pf.tidy(s)
test_that('pf.tidy works with X~Data regression',expect_true(exists('td')))


#### Multiple Regression

a <- rnorm(length(X))
A <- data.frame('a'=a,'Sample_Site'=X)
rm('PF')

PF <- PhyloFactor(Data,tree,X=A,frmla=Sample_Site~a+Data,family=binomial(link='logit'),nfactors=1,choice='var')

test_that('X~Data+a multiple regression works',expect_true(exists('PF')))