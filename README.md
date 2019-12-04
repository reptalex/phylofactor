## phylofactor: R package for phylogentic factorization of biological data


This package provides functions to extract and visualize the phylogenetic factors in biological datasets. To get started, download R and Rtools. In R, you first need to download two packages from Bioconductor:

```{r install}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Biostrings","ggtree"))
```
Now, you're ready to download the \code{phylofactor} packsage using \code{devtools}:

```{r install}
devtools::install_github('reptalex/phylofactor')
```
The help functions in this package contain a decent amount of detail for how to perform phylofactorization.

```{r PhyloFactor}
? PhyloFactor
? gpf
```

For more detailed information, a tutorial is available online [here](https://docs.wixstatic.com/ugd/0119a1_099ae20df8424af9a38585dcebc0d45a.pdf "Phylofactor Tutorial").  If you have any questions, please feel free to post issues and contact me!
