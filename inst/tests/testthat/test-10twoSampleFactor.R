context('Checking core PhyloFactor functionality and summary functions')

set.seed(2)
data("FTmicrobiome",package = 'phylofactor')
tree <- FTmicrobiome$tree
tree <- ape::drop.tip(tree,setdiff(tree$tip.label,sample(tree$tip.label,20)))

sigClades <- list(c(16,17),c(6,7,8,9,10,11,12,13,14,15))

Z <- rep(0,20)
Z[sigClades[[1]]] <- 1
Z[sigClades[[2]]] <- 1

test_that('twoSampleFactor with method="Fisher" works',
          expect_error(invisible(capture.output(pf <- twoSampleFactor(Z,tree,2,method='Fisher')))))

invisible(capture.output(pf <- twoSampleFactor(Z,tree,2,method='Fisher')))

test_that('twoSampleFactor with method="Fisher" captures correct clades',
          expect_true(all.equal(sigClades[[2]],pf$bins[[2]]) & all.equal(sigClades[[1]],pf$bins[[3]])))

set.seed(1)
D <- 300
tree <- rtree(D)
n1 <- 477
n2 <- 332
c1 <- phangorn::Descendants(tree,n1,'tips')[[1]]
c2 <- phangorn::Descendants(tree,n2,'tips')[[1]]

Z <- rnorm(D)
Z[c1] <- Z[c1]+1
Z[c2] <- Z[c2]-2


test_that('twoSampleFactor works in parallel',
          expect_error(invisible(capture.output(pf <- twoSampleFactor(Z,tree,2,ncores=2))),NA))

invisible(capture.output(pf <- twoSampleFactor(Z,tree,2,ncores=2)))

test_that('twoSampleFactor captures correct clades',
          expect_true(all(sapply(pf$bins,length)==c(126,50,124))))