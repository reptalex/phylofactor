## phylofactor: R package for phylogentic factorization of biological data

This package provides functions to extract and visualize the phylogenetic factors in biological datasets. To get started, download the packsage using \code{devtools}

```{r cars}
devtools::install_github('reptalex/phylofactor')
? Phylofactor
```

For more details, a tutorial is available online [here](http://media.wix.com/ugd/0119a1_f9ff04d5ada440829cc2d942b8b9f928.pdf "Phylofactor Tutorial").  If you have any questions, please feel free to contact me!

``` {r echo=FALSE}
library(phylofactor)
data('FTmicrobiome')
PF <- FTmicrobiome$PF
binPhyloPlot(PF,factor=3,edge.width=2)
```