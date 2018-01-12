context('Checking gpf')

set.seed(1)
m <- 50
n <- 200
tree <- ape::rtree(m)
X <- data.table::data.table('y'=rnorm(n),
               'z'=rnorm(n,sd=0.5),
               'Sample'=sapply(1:n,FUN=function(s) paste('Sample',s)),
               key='Sample')
#we'll partition by 'y'.
binom.size=3
clade <- phangorn::Descendants(tree,75,'tips')[[1]]
clade2 <- phangorn::Descendants(tree,53,'tips')[[1]]

######## presence/absence dataset with affected clade #######
## most species have higher P{present} with y
ilogit <- function(eta) 1/(1+exp(-eta))
eta <- X$z+X$y
p <- ilogit(eta)
M <- matrix(rbinom(m*n,binom.size,rep(p,times=m)),nrow=m,ncol=n,byrow=T)
rownames(M) <- tree$tip.label
colnames(M) <- X$Sample

#### the first clade decreases with y ####
eta1 <- X$z-X$y
p1 <- ilogit(eta1)
for (species in clade){
  M[species,] <- rbinom(n,binom.size,p1)
}
#### the second clade weakly decreases with y ####
eta2 <- X$z-.3*X$y
p2 <- ilogit(eta2)
for (species in clade2){
  M[species,] <- rbinom(n,binom.size,p2)
}

####################### to partition on y, must have phylo* #########
invisible(capture.output(pf <- gpf(M,tree,X,frmla.phylo=cbind(Successes,Failures)~z+phylo*y,nfactors=2,
         binom.size=binom.size,family=binomial(link='logit'),
         PartitioningVariables='y',algorithm='mStable')))
test_that('algorithm "mStable" works',expect_true(all.equal(pf$groups[[1]][[1]],clade)))

### partition vector of data controlling for sample effort ###
set.seed(1)
effort <- rnorm(50)
eta <- effort-3
eta[clade] <- eta[clade]+6
eta[clade2] <- eta[clade2]+8
Data <- data.table::data.table('Species'=tree$tip.label,effort,Z=rbinom(50,1,ilogit(eta)),'Sample'=1)

invisible(capture.output(pf <- gpf(Data,tree,frmla.phylo=Z~effort+phylo,nfactors=2,algorithm='phylo',family=binomial)))
test_that('algorithm "phylo" works', expect_true(all.equal(pf$groups[[1]][[1]],clade) & all.equal(pf$groups[[2]][[1]],clade2)))



DF <- matrix.to.phyloframe(M,data.name='Successes')
DF[,Failures:=binom.size-Successes]
data.table::setkey(DF,Sample)
DF <- data.table:::`[.data.table`(DF,X)
invisible(capture.output(pf <- gpf(DF,tree,frmla=cbind(Successes,Failures)~z+y,
                               PartitioningVariables='y',
                               algorithm='CoefContrast',
                               family=binomial,
                               nfactors=2)))
test_that('algorithm "CoefContrast" works', expect_true(all.equal(pf$groups[[1]][[1]],clade) & all.equal(pf$groups[[2]][[1]],clade2)))
