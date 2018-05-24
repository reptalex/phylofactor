#' Fecal Tongue Microbiome Dataset
#' 
#' A dataset containing an OTU table, tree, sample-site meta-data, taxonomy and
#' phylofactor objects with 100 factors each.
#' 
#' @format List containing 6 elements:
#' \describe{
#'    \item{OTUTable}{OTUTable containing 2713 OTUs and 40 samples}
#'    \item{tree}{phylo object containing the OTUs in \code{rownames(OTUTable)}}
#'    \item{X}{sample-site vector corresponding to columns of \code{OTUTable}}
#'    \item{taxonomy}{list containing mapping of OTU ids in columns of \code{OTUTable}
#'                    to greengenes taxonomy (must be semicolon-separated strings)}
#'    \item{PF}{default phylofactorization of trimmed \code{OTUTable} with \code{choice='var'}}
#'    \item{PF.stat}{phylofactorization of trimmed \code{OTUTable} with \code{choice='F'}}
#' }
"FTmicrobiome"