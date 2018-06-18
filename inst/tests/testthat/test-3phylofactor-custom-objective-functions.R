context('Checking PhyloFactor with customized objective functions GAM and HUTCH')
options(warn=-1)
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

load.dependencies <- 'library(mgcv)'

invisible(capture.output(suppressWarnings(PF.G <- PhyloFactor(Data,tree,X,nfactors=2,choice.fcn=GAM,cluster.depends = load.dependencies,sp=c(1,1)))))
invisible(capture.output(suppressWarnings(PF.G.par <- PhyloFactor(Data,tree,X,nfactors=2,choice.fcn=GAM,cluster.depends = load.dependencies,ncores=2,sp=c(1,1)))))
test_that('Parellelized & serialized customized objective function - GAM - are not equal',expect_true(all.equal(PF.G,PF.G.par)))
names(PF.G$bins) <- NULL
test_that('GAM phylofactorization did not extract correct clades',expect_true(all.equal(unlist(sigClades),unlist(PF.G$bins[2:3]))))

PF.G$custom.output

invisible(capture.output(suppressWarnings(pf.gg <- PhyloFactor(Data,tree,X,nfactors=2,frmla=Data~s(a)+s(b),method='gam',sp=c(1,1)))))

test_that('built-in method="gam" is equal to customized GAM',
          expect_equal(lapply(PF.G$custom.output,coef),lapply(pf.gg$custom.output,coef)))
 ######################## Finding Hutchisonian Niches #####################################
 ### Example of how to use PhyloFactor to identify Gaussian-shapped Hutchinsonian niches ###
 set.seed(1)
 n=1000
 A <- 20
 mu=-1
 sigma=0.9
 Data <- matrix(rlnorm(20*n,meanlog = 8,sdlog = .5),nrow=20)
 rownames(Data) <- tree$tip.label
 X <- rnorm(n)
 Data[sigClades[[1]],] <- t(t(Data[sigClades[[1]],])*A*exp(-(((X-mu)^2)/(2*sigma^2))))
 Data <- t(clo(t(Data)))

 frmla=Data~X+I(X^2)
 invisible(capture.output(PF.Gaus <- PhyloFactor(Data,tree,frmla=frmla,X,nfactors=1,ncores=7)))

 test_that('Correctly identified Clade with Gaussian-shaped Niche',expect_true(all.equal(sigClades[[1]],PF.Gaus$bins[[2]])))
 options(warn=0)
 