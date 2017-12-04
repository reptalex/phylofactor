## phylofactor: R package for phylogentic factorization of biological data

This package provides functions to extract and visualize the phylogenetic factors in biological datasets. To get started, first you need to download two packages from Bioconductor:

```{r install}
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
biocLite("ggtree")
```
Now, you're ready to download the \code{phylofactor} packsage using \code{devtools}:

```{r install}
devtools::install_github('reptalex/phylofactor')
? PhyloFactor
```

For detailed information, a tutorial is available online [here](http://media.wix.com/ugd/0119a1_951b32eb0abe4f228f0d6fd4ae11a0e8.pdf "Phylofactor Tutorial").  If you have any questions, please feel free to post issues and contact me!

```{r plot}
library(phylofactor)
data('FTmicrobiome')
pf <- FTmicrobiome$PF
pf.tree(pf,factors=1:3)
```