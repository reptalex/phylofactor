context('Checking gpf')

set.seed(1)
m <- 50
n <- 50
tree <- ape::rtree(m)
MetaData <- data.table::data.table('y'=rnorm(n),
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
eta <- MetaData$z+MetaData$y
p <- ilogit(eta)
M <- matrix(rbinom(m*n,binom.size,rep(p,times=m)),nrow=m,ncol=n,byrow=T)
rownames(M) <- tree$tip.label
colnames(M) <- MetaData$Sample

#### the first clade decreases with y ####
eta1 <- MetaData$z-1.5*MetaData$y
p1 <- ilogit(eta1)
for (species in clade){
  M[species,] <- rbinom(n,binom.size,p1)
}
#### the second clade weakly decreases with y ####
eta2 <- MetaData$z-.7*MetaData$y
p2 <- ilogit(eta2)
for (species in clade2){
  M[species,] <- rbinom(n,binom.size,p2)
}

Data <- list('Successes'=M,'Failures'=binom.size-M)

####################### to partition on y, must have phylo* #########
invisible(capture.output(pf <- gpf(Data,tree,MetaData=MetaData,
                                   frmla.phylo=cbind(Successes,Failures)~z+phylo*y,nfactors=2,
                                   family=binomial(link='logit'),
                                   PartitioningVariables='y',algorithm='mStable')))
test_that('algorithm "mStable"',expect_true(all.equal(pf$groups[[1]][[1]],clade)))

### partition vector of data controlling for sample effort ###
set.seed(1)
effort <- rnorm(50)
eta <- effort-3
eta[clade] <- eta[clade]+6
eta[clade2] <- eta[clade2]+8
Data <- data.table::data.table('Species'=tree$tip.label,effort,Z=rbinom(50,1,ilogit(eta)),'Sample'=1)

invisible(capture.output(pf <- gpf(Data,tree,frmla.phylo=Z~effort+phylo,nfactors=2,algorithm='phylo',family=binomial)))
test_that('algorithm "phylo"', expect_true(all.equal(pf$groups[[1]][[1]],clade) & all.equal(pf$groups[[2]][[1]],clade2)))

w <- 1:nrow(Data)
invisible(capture.output(pf <- gpf(Data,tree,frmla.phylo=Z~effort+phylo,nfactors=1,algorithm='phylo',family=binomial,weights=w)))
test_that('algorithm "phylo" with weights', expect_true(all(pf$models[[1]]$prior.weights==w)))


w <- 1:nrow(Data)
invisible(capture.output(pf <- gpf(Data,tree,frmla.phylo=Z~s(effort)+phylo,nfactors=1,algorithm='phylo',model.fcn=mgcv::gam,family=binomial,weights=w)))
test_that('algorithm "phylo" with gam and weights', expect_true(all(pf$models[[1]]$prior.weights==w)))


DF <- matrix.to.phyloframe(M,data.name='Successes')
DF[,Failures:=binom.size-Successes]
data.table::setkey(DF,Sample)
DF <- data.table:::`[.data.table`(DF,MetaData)
invisible(capture.output(pf <- gpf(DF,tree,frmla=cbind(Successes,Failures)~z+y,
                               PartitioningVariables='y',
                               algorithm='CoefContrast',
                               family=binomial,
                               nfactors=2)))
test_that('algorithm "CoefContrast"', expect_true(all.equal(pf$groups[[1]][[1]],clade) & all.equal(pf$groups[[2]][[1]],clade2)))



### test gam-functionality
set.seed(1)
m <- 50
n <- 20
tree <- rtree(m)
MetaData <- data.table('y'=rnorm(n),
                       'z'=rnorm(n,sd=0.5),
                       'Sample'=sapply(1:n,FUN=function(s) paste('Sample',s)),
                       key='Sample')
#we'll partition by 'y'.
binom.size=10
clade <- phangorn::Descendants(tree,75,'tips')[[1]]
clade2 <- phangorn::Descendants(tree,53,'tips')[[1]]

## most species have higher P{present} with y
eta <- .5*MetaData$z+MetaData$y
p <- ilogit(eta)
M <- matrix(rbinom(m*n,binom.size,rep(p,times=m)),nrow=m,ncol=n,byrow=T)
rownames(M) <- tree$tip.label
colnames(M) <- MetaData$Sample

#### the first clade decreases with y ####
eta1 <- .5*MetaData$z+2*MetaData$y^2
p1 <- ilogit(eta1)
for (species in clade){
  M[species,] <- rbinom(n,binom.size,p1)
}
#### the second clade weakly decreases with y ####
eta2 <- .5*MetaData$z-4*MetaData$y
p2 <- ilogit(eta2)
for (species in clade2){
  M[species,] <- rbinom(n,binom.size,p2)
}

DF <- matrix.to.phyloframe(M,data.name='Successes')
DF[,Failures:=binom.size-Successes]
setkey(DF,Sample)
DF <- DF[MetaData]


model.z <- glm(cbind(Successes,Failures)~z,family=binomial,data=DF)
DF[,z.fit:=coef(model.z)['z']*z]

invisible(capture.output(pf.gam <- gpf(DF,tree,frmla.phylo=cbind(Successes,Failures)~offset(z.fit)+s(y,by=phylo),
          PartitioningVariables='y',family=binomial,nfactors=2,ncores=2,
          model.fcn = mgcv::gam,algorithm = 'phylo')))

testthat::test_that('model.fcn=mgcv::gam, algorithm "phylo" works and produces the correct factors', 
          testthat::expect_true(all.equal(pf.gam$groups[[1]][[1]],clade2) & all.equal(pf.gam$groups[[2]][[1]],clade)))



samples <- DF$Sample
###################################
DF$Sample <- NULL
##serial
invisible(capture.output(pf <- gpf(DF,tree,frmla.phylo=cbind(Successes,Failures)~offset(z.fit)+s(y,by=phylo),
          PartitioningVariables='y',family=binomial,
          nfactors=2,model.fcn = mgcv::gam,algorithm = 'phylo')))

testthat::test_that('Binomial gam with offset and smoothing',
                    testthat::expect_true(all.equal(pf$groups[[1]][[1]],clade2) & all.equal(pf$groups[[2]][[1]],clade)))

##parallel
invisible(capture.output(pf <- gpf(DF,tree,frmla.phylo=cbind(Successes,Failures)~offset(z.fit)+s(y,by=phylo),
          PartitioningVariables='y',family=binomial,
          nfactors=2,ncores=2,model.fcn = mgcv::gam,algorithm = 'phylo')))
testthat::test_that('Parallelized Binomial gam with offset and smoothing',
                    testthat::expect_true(all.equal(pf$groups[[1]][[1]],clade2) & all.equal(pf$groups[[2]][[1]],clade)))



## test parallelized, mix-algorithm, glm with offset
invisible(capture.output(pf <- gpf(DF,tree,frmla=cbind(Successes,Failures)~offset(z.fit)+y,
          frmla.phylo=cbind(Successes,Failures)~offset(z.fit)+y*phylo,
          PartitioningVariables='y',family=binomial,nfactors=2,ncores=2)))
testthat::test_that('Parallelized binomial regression, algorithm=mix with offset and smoothing',
                    testthat::expect_true(all.equal(pf$groups[[1]][[1]],clade2) & all.equal(pf$groups[[2]][[1]],clade)))


##test binomial mStable
DF$Sample <- samples
invisible(capture.output(pf <- gpf(DF,tree,frmla.phylo=cbind(Successes,Failures)~offset(z.fit)+y*phylo,
          PartitioningVariables='y',family=binomial,
          nfactors=2,algorithm = 'mStable')))
testthat::test_that('mStable binomial regression with offset',
                    testthat::expect_true(all.equal(pf$groups[[1]][[1]],clade2) & all.equal(pf$groups[[2]][[1]],clade)))


##test mStable with gam
invisible(capture.output(pf <- gpf(DF,tree,frmla.phylo=cbind(Successes,Failures)~offset(z.fit)+s(y,by=phylo),
          PartitioningVariables='y',family=binomial,
          nfactors=2,model.fcn = mgcv::gam,algorithm = 'mStable')))
testthat::test_that('mStable binomial gam with offset',
                    testthat::expect_true(all.equal(pf$groups[[1]][[1]],clade2) & all.equal(pf$groups[[2]][[1]],clade)))


#### testing other algorithms:
## CoefContrast with poisson
invisible(capture.output(pf <- gpf(DF,tree,frmla=Successes~offset(z.fit)+y,
          PartitioningVariables = 'y',family=poisson,nfactors=2,
          algorithm='CoefContrast')))
testthat::test_that('CoefContrast poisson regression with offset',
                    testthat::expect_true(all.equal(pf$groups[[1]][[1]],clade2) & all.equal(pf$groups[[2]][[1]],clade)))

## CoefContrast with binomial
invisible(capture.output(pf <- gpf(DF,tree,frmla=cbind(Successes,Failures)~offset(z.fit)+y,
          PartitioningVariables = 'y',family=binomial,nfactors=2,
          algorithm='CoefContrast')))
testthat::test_that('CoefContrast binomial regression with offset',
                    testthat::expect_true(all.equal(pf$groups[[1]][[1]],clade2) & all.equal(pf$groups[[2]][[1]],clade)))


## mix with binomial
invisible(capture.output(pf <- gpf(DF,tree,frmla=cbind(Successes,Failures)~offset(z.fit)+y,
          frmla.phylo=cbind(Successes,Failures)~offset(z.fit)+y*phylo,
          PartitioningVariables = 'y',family=binomial,nfactors=2)))
testthat::test_that('mix algorithm binomial regression with offset',
                    testthat::expect_true(all.equal(pf$groups[[1]][[1]],clade2) & all.equal(pf$groups[[2]][[1]],clade)))


### try mStable with input matrix:
MetaData <- DF[,c('Sample','y','z','z.fit')]
MetaData <- MetaData[!duplicated(MetaData)]
Mats <- list('Successes'=M,'Failures'=binom.size-M)
invisible(capture.output(pf <- gpf(Mats,tree,MetaData=MetaData,
          frmla.phylo=cbind(Successes,Failures)~offset(z.fit)+phylo*y,
          algorithm='mStable',nfactors=2,family=binomial)))
testthat::test_that('mStable binomial regression with offset and input matrices+MetaData',
                    testthat::expect_true(all.equal(pf$groups[[1]][[1]],clade2) & all.equal(pf$groups[[2]][[1]],clade)))


### check error for improper input of one matrix for binomial
invisible(capture.output(err <- tryCatch(pf <- gpf(Mats[[1]],tree,MetaData=MetaData,
          frmla.phylo=cbind(Successes,Failures)~offset(z.fit)+phylo*y,
          algorithm='mStable',nfactors=2,family=binomial),error=function(e) e)))
testthat::test_that('Input only one matrix for binomial mStable returns appropriate error',
                    testthat::expect_true(grepl('for family=binomial. Must be either data frame with Successes',as.character(err))))



## minimal input: frmla.phylo defines partitioning variables as all those with phylo* or *phylo.
invisible(capture.output(pf <- gpf(DF,tree,frmla.phylo=Successes~offset(z.fit)+phylo*y,
          nfactors=2,family=poisson)))
testthat::test_that('Minimial input - data tree, frmla.phylo - works with family=poisson',
                    testthat::expect_true(all.equal(pf$groups[[1]][[1]],clade2) & all.equal(pf$groups[[2]][[1]],clade)))
testthat::test_that('Minimal input Partitioning Variables are correct',
                    testthat::expect_true(pf$PartitioningVariables=='y'))

### single vector of data
set.seed(1)
effort <- .5*rnorm(m)
eta <- 3+0.5*effort
eta[clade] <- eta[clade]+1
eta[clade2] <- eta[clade2]-2
DF <- data.frame('Species'=tree$tip.label,'N'=rpois(m,exp(eta)),'effort'=effort)
fit <- glm(N~effort,family=poisson,data=DF)
DF$effort.fit <- coef(fit)['effort']*effort


invisible(capture.output(pf <- gpf(DF,tree,frmla.phylo = N~phylo,algorithm='phylo',nfactors=2)))
testthat::test_that('gpf can partition single vectors of data without "Sample" column',
                    testthat::expect_true(all.equal(pf$groups[[1]][[1]],clade)&all.equal(pf$groups[[2]][[1]],clade2)))

#offset effort with a global effect of effort on N
invisible(capture.output(pf <- gpf(DF,tree,frmla.phylo = N~offset(effort.fit)+phylo,algorithm='phylo',nfactors=2)))
testthat::test_that('gpf can partition single vectors of data with offset',
                    testthat::expect_true(length(pf$models[[1]]$offset)==m))