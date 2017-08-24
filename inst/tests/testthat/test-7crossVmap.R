context('Checking crossVmap')


set.seed(1)
tree <- ape::rtree(7)
c1 <- tree$tip.label[1:6]
c2 <- setdiff(tree$tip.label,tree$tip.label[c(3,4,5)])

grp1 <- tree$tip.label[c(1,2)]
grp2 <- tree$tip.label[c(5,6)]
Grps <- list(grp1,grp2)

test_that('CrossVmap works when ignoring interruptions & no interruptions present',
  expect_true(all.equal(crossVmap(Grps,tree,c1,c2),
                        list('9'=c('t7','t2'),'13'='t3')))
)

grp1 <- tree$tip.label[c(3,4)]
grp2 <- tree$tip.label[c(5)]
Grps <- list(grp1,grp2)
c1 <- c(grp1,grp2,tree$tip.label[6])
c2 <- setdiff(tree$tip.label,tree$tip.label[4])

##correct output should be [[t3]],[[t5]]
test_that('crossVmap works when ignoring interruptions & interruptions are present',
          expect_true(all.equal(crossVmap(Grps,tree,c1,c2),
                                list('11'='t5','5'='t6'))))

test_that('crossVmap works when not ignoring interruptions',
          expect_true(all.equal(crossVmap(Grps,tree,c1,c2,ignore.interruptions = FALSE),
                                list('11'='t5','10'=c('t7','t2'),'12'='t1','5'='t6'))))
