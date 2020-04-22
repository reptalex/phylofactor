#' Summarize phylofactor objects
#' @export
#' @param PF phylofactor object
#' @param taxonomy taxonomy (optional) for taxonomic trimming. First column contains tree tip-labels, second column has semicolon-delimited taxonomic character strings.
#' @param factor Integer (optional). If NULL, will use print.phylofactor
#' @param taxon.trimming Chracter string. How to trim taxonomic groups. Must be either "sup" (group by shortest unique prefix for every taxon in group distinguishing it from taxa in complementary group), "taxon" for taxon levels in second column of \code{taxonomy}, or "species" for no taxonomic binning.
#' @param output.signal logical - whether or not to output signal. Will attempt to call model.fcn or customized choice.fcn to estimate signal (and thus relative importance) for each taxonomic group from \code{taxon.trimming}.
#' @return phylofactor.summary object containing data, group summary and, with input taxonomy, taxonomic break-down of signal
summary.phylofactor <- function(PF,taxonomy=NULL,factor=NULL,taxon.trimming='sup',output.signal=T){
  ### must output the following:
  # (0) factors[n,] description - if null, summarize whole object w/ biggest, #signif at holm cutoff, adjusted P-vals etc.
  # (1) model summary
  # (2) data for this factor
        # - PhyloFactor & PhyCA: transform.fcn --> contrast.fcn and OTUTable
        # - gpf: data frame input to model
        # - twoSample: NA
        # - ALL: fitted.values
  # (3) taxa table: (optional sup, full-taxon-binning,or species-level data)
        # - #species, taxa, mean (sorted), signal (tbd, possibly sorted)
  # (4) taxa.split
  # (5) species.lists
  if (is.null(factor)){
    print(PF)
  } else {
    if (factor>PF$nfactors){
      stop('factor input exceeds number of factors in phylofactor object')
    }
    if (is.null(taxonomy)){
      taxon.trimming <- 'species'
    } else {
      if (!all(class(taxonomy)=='data.frame')){
        taxonomy <- as.data.frame(taxonomy)
      }
      if (taxon.trimming %in% c('sup','taxon')){
        brks <- grepl('\\[',taxonomy[,2]) | grepl('\\]',taxonomy[,2])
        if (any(brks)){
          taxonomy[brks,2] <- gsub('\\[','~',taxonomy[brks,2])
          taxonomy[brks,2] <- gsub('\\]','~',taxonomy[brks,2])
        }
        curls <- grepl('\\{',taxonomy[,2]) | grepl('\\}',taxonomy[,2])
        if (any(curls)){
          taxonomy[curls,2] <- gsub('\\{','~',taxonomy[curls,2])
          taxonomy[curls,2] <- gsub('\\}','~',taxonomy[curls,2])
        }
        parens <- grepl('\\(',taxonomy[,2]) | grepl('\\)',taxonomy[,2])
        if (any(parens)){
          taxonomy[parens,2] <- gsub('\\(','~',taxonomy[parens,2])
          taxonomy[parens,2] <- gsub('\\)','~',taxonomy[parens,2])
        }
      }
    }
    if (!taxon.trimming %in% c('sup','taxon','species')){
      stop('unknown input taxon.trimming. Must be either "sup", "taxon", or "species".')
    }


    output <- NULL
    output$group.summary <- PF$factors[factor,]

    if (factor==1){
      Grps <- getPhyloGroups(PF$tree)
    } else {
      Grps <- pf.getPhyloGroups(PF,factor)
    }

    ### model.summary
    if (PF$phylofactor.fcn == 'PhyloFactor'){
      if (!is.null(PF$custom.output)){
        output$model.summary <- tryCatch(summary(PF$custom.output[[factor]]),
                 error=function(e) 'could not obtain summary for customized output')
      } else {
        output$model.summary <- summary(PF$models[[factor]])
      }
    } else if (PF$phylofactor.fcn=='gpf'){
      if (PF$algorithm!='CoefContrast'){
        output$model.summary <- summary(PF$models[[factor]])
      } else {
        if (all(PF$PartitioningVariables %in% colnames(PF$coefficient.matrix))){
          B <- t(PF$coefficient.matrix[,PF$PartitioningVariables,drop=F]/PF$coefficient.SE[,PF$PartitioningVariables,drop=F]) %*% PF$basis[,factor]
        } else {
          pvs_not_found <- PF$PartitioningVariables[!PF$PartitioningVariables %in% colnames(PF$coefficient.matrix)]
          pvs_found <- setdiff(PF$PartitioningVariables,pvs_not_found)
          ix <- sapply(pvs_not_found,FUN=function(a,b) grepl(a,b),b=colnames(PF$coefficient.matrix)) %>% apply(MARGIN=1,any)
          ix <- ix | colnames(PF$coefficient.matrix) %in% pvs_found
          if (!any(ix)){
            stop('Could neither match nor grep any PartitioningVariables in the coefficients of model. Try running stats::coefficients on your input model.fcn for a single species to determine the appropriate names for PartitioningVariables.')
          }
          B <- t(PF$coefficient.matrix[,ix,drop=F]/PF$coefficient.SE[,ix,drop=F]) %*% PF$basis[,factor]
        }
        output$model.summary <- data.frame('PartitioningVariable'=PF$PartitioningVariables,'CoefContrast'=B)
      }
    } else if (PF$phylofactor.fcn=='PhyCA'){
      output$model.summary <- data.frame('Percent_Explained_Variance'=PF$factors[factor,'ExpVar'])
    } else {
      output$model.summary <- data.frame('method'=PF$method,
                                         'Test_Statistic'=PF$objective[factor],
                                         'Pval'=PF$pvals[factor])
    }

    ### data ###
    if (PF$phylofactor.fcn %in% c('PhyloFactor','PhyCA')){
      y <- PF$contrast.fcn(PF$groups[[factor]],PF$transform.fcn(PF$Data))
      if (is.null(PF$X)){
        output$data <- data.frame('Data'=y)
      } else {
        if (!all(class(PF$X)=='data.frame')){
          output$data <- data.frame('Data'=y,'X'=PF$X)
        } else {
          output$data <- cbind(data.frame('Data'=y),PF$X)
        }
      }
      if (!is.null(PF$models)){
        output$data$fitted.values <- PF$models[[factor]]$fitted.values
      }
    } else if (PF$phylofactor.fcn == 'gpf'){
      if (PF$algorithm!='CoefContrast'){
        output$data <- PF$models[[factor]]$data
        output$data$fitted.values <- NA
        output$data$fitted.values[setdiff(1:nrow(output$data),PF$models[[factor]]$na.action)] <- PF$models[[factor]]$fitted.values
      } else {
        output$data <- list('Coefficient_matrix'=PF$coefficient.matrix,
                            'Coefficient_SE'=PF$coefficient.SE,
                            'Contrast_Vector'=PF$basis[,factor,drop=F],
                            'Partitioning_Variables'=PF$PartitioningVariables)
      }
    } else {
      output$data <- list('Group1'=PF$Data[PF$groups[[factor]][[1]]],
                          'Group2'=PF$Data[PF$groups[[factor]][[2]]])
    }

    ### species.list
    output$species.list <- pf.groupsTospecies(PF)[[factor]]

    ### taxa.split
    if (!is.null(taxonomy)){
      output$taxa.split <- pf.taxa(PF,taxonomy,factor)
    } else {
      output$taxa.split <- NULL
    }


    ### taxon.table ##################################################
    if (taxon.trimming=='species'){
      spp1 <- output$species.list[[1]]
      spp2 <- output$species.list[[2]]
      group1.taxa <- spp1
      group2.taxa <- spp2
      species.assignment1 <- spp1
      species.assignment2 <- spp2

    } else {
      group1.taxonomy <- taxonomy[match(output$species.list[[1]],as.character(taxonomy[,1])),]
      group2.taxonomy <- taxonomy[match(output$species.list[[2]],as.character(taxonomy[,1])),]
      if (taxon.trimming=='sup'){
        group1.taxa <- output$taxa.split[[1]]
        group2.taxa <- output$taxa.split[[2]]
      } else {
        group1.taxa <- unique(group1.taxonomy[,2])
        group2.taxa <- unique(group2.taxonomy[,2])
      }
      group1.taxa <- group1.taxa[order(nchar(group1.taxa),decreasing = T)]
      group2.taxa <- group2.taxa[order(nchar(group2.taxa),decreasing = T)]

      species.assignment1 <- vector(mode='list',length=length(group1.taxa))
      names(species.assignment1) <- group1.taxa
      species.assignment2 <- vector(mode='list',length=length(group2.taxa))
      names(species.assignment2) <- group2.taxa

      if (length(group1.taxa)==1){
        tx <- 1
      } else {
        tx <- sapply(group1.taxa,FUN=function(tax,taxa) grepl(tax,taxa),group1.taxonomy[,2]) %>%
          apply(MARGIN=1,FUN=function(g) min(which(g)))
      }
      for (i in 1:length(group1.taxa)){
        species.assignment1[[i]] <- output$species.list[[1]][tx==i]
      }

      if (length(group2.taxa)==1){
        tx <- 1
      } else {
        tx <- sapply(group2.taxa,FUN=function(tax,taxa) grepl(tax,taxa),group2.taxonomy[,2]) %>%
                    apply(MARGIN=1,FUN=function(g) min(which(g)))
      }
      for (i in 1:length(group2.taxa)){
        species.assignment2[[i]] <- output$species.list[[2]][tx==i]
      }

      output$taxon.species.assignments <- list('group1'=species.assignment1,
                                               'group2'=species.assignment2)
    }


    Tbl1 <- data.frame('Taxon'=group1.taxa,'nSpecies'=sapply(species.assignment1,length))
    Tbl2 <- data.frame('Taxon'=group2.taxa,'nSpecies'=sapply(species.assignment2,length))


    if (PF$phylofactor.fcn %in% c('PhyloFactor','PhyCA')){
      transform.means1 <- sapply(species.assignment1,FUN=function(ix,PF) mean(PF$transform.fcn(PF$Data[ix,,drop=F]),na.rm=T),PF)
      transform.vars1 <- sapply(species.assignment1,FUN=function(ix,PF) mean(apply(PF$transform.fcn(PF$Data[ix,,drop=F]),1,var,na.rm=T),na.rm=T),PF)
      transform.means2 <- sapply(species.assignment2,FUN=function(ix,PF) mean(PF$transform.fcn(PF$Data[ix,,drop=F]),na.rm=T),PF)
      transform.vars2 <- sapply(species.assignment2,FUN=function(ix,PF) mean(apply(PF$transform.fcn(PF$Data[ix,,drop=F]),1,var,na.rm=T),na.rm=T),PF)
      raw.means1 <- sapply(species.assignment1,FUN=function(ix,PF) mean(PF$Data[ix,,drop=F],na.rm=T),PF)
      raw.vars1 <- sapply(species.assignment1,FUN=function(ix,PF) mean(apply(PF$Data[ix,,drop=F],1,var,na.rm=T),na.rm=T),PF)
      raw.means2 <- sapply(species.assignment2,FUN=function(ix,PF) mean(PF$Data[ix,,drop=F],na.rm=T),PF)
      raw.vars2 <- sapply(species.assignment2,FUN=function(ix,PF) mean(apply(PF$Data[ix,,drop=F],1,var,na.rm=T),na.rm=T),PF)
      raw.sum1 <- mapply(raw.means1,FUN = function(a,b) a*b,b=sapply(species.assignment1,length))
      raw.sum2 <- mapply(raw.means2,FUN = function(a,b) a*b,b=sapply(species.assignment2,length))

      Stats1 <- data.frame('raw.sum'=raw.sum1,'raw.mean'=raw.means1,'raw.var'=raw.vars1,
                           'transformed.mean'=transform.means1,'transformed.var'=transform.vars1,row.names = NULL)
      Stats2 <- data.frame('raw.sum'=raw.sum2,'raw.mean'=raw.means2,'raw.var'=raw.vars2,
                         'transformed.mean'=transform.means2,'transformed.var'=transform.vars2,row.names=NULL)

      Tbl1 <- cbind(Tbl1,Stats1)
      Tbl2 <- cbind(Tbl2,Stats2)
      Tbl1 <- Tbl1[order(Tbl1$raw.sum,decreasing = T),]
      Tbl2 <- Tbl2[order(Tbl2$raw.sum,decreasing = T),]

    }

    if (output.signal){
      if (PF$phylofactor.fcn=='gpf' & !PF$algorithm=='mStable'){
        if (!key(PF$Data)=='Species'){
          setkey(PF$Data,Species)
        }
      }
      Grps1 <- lapply(species.assignment1,FUN=function(spp,tree) match(spp,tree$tip.label),PF$tree)
      Grps2 <- lapply(species.assignment2,FUN=function(spp,tree) match(spp,tree$tip.label),PF$tree)
      Grps1 <- lapply(Grps1,FUN=function(g1,g2) list(g1,g2),g2=unlist(PF$groups[[factor]][2]))
      Grps2 <- lapply(Grps2,FUN=function(g1,g2) list(g2,g1),g2=unlist(PF$groups[[factor]][1]))

      Tbl1$signal <- sapply(Grps1,getSignal,PF)
      Tbl2$signal <- sapply(Grps2,getSignal,PF)
      Tbl1 <- Tbl1[order(Tbl1$signal,decreasing = T),]
      Tbl2 <- Tbl2[order(Tbl2$signal,decreasing = T),]
    }
    rownames(Tbl1) <- NULL
    rownames(Tbl2) <- NULL
    output$signal.table <- list('Group1'=Tbl1,'Group2'=Tbl2)

    if (PF$phylofactor.fcn %in% c('PhyloFactor','gpf')){
      if (!is.null(PF$models)){
        formula <- Reduce(paste,deparse(PF$models[[1]]$formula))
      } else {
        if (PF$phylofactor.fcn=='PhyloFactor'){
          if (PF$choice=='custom'){
            
            if (is.null(tryCatch(PF$custom.output[[1]]$formula,error=function(e) NULL))){
              formula <- 'Customized'
            } else {
              formula <- paste('
                               Formula                   : ',Reduce(paste,deparse(PF$custom.output[[1]]$formula)),sep='')
            }
          } else if (PF$method=='max.var'){
              formula <- NULL
          }
        } else {
            if (PF$algorithm=='CoefContrast'){
              formula <- Reduce(paste,deparse(PF$species.models[[1]]$formula))
            }
        }
      }
    } else {
      formula <- NULL
    }
    output$info <- list('phylofactor.fcn'=PF$phylofactor.fcn,
                        'method'=PF$method,'choice'=PF$choice,
                        'algorithm'=PF$algorithm,'formula'=formula,
                        'factor'=factor,'nEdges'=length(Grps))
    class(output) <- 'phylofactor.summary'
    return(output)
  }
}