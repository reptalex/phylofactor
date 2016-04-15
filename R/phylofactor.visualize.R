#' Phylofactor Ordination-Visualization
#'
#' @param PF Phylofactor object
#' @param dimension Either 2 or 3 - the dimension of the ordination-visualization plot
#' @param default.colors Logical indicating if samples should be colored using default settings. If default.colors=F, user must provide appropriate color arguments for plot or scatterplot3d, and no legend will be made.
#' @param Legend Logical indicating if legend should be included. Must also have default.colors=T
#' @param lx x location for legend. If \code{dimension} = 2, x must be a numeric. If \code{dimension} = 3, x must be an scatterplot3d()$xyz.convert list
#' @param ly y location for legend, only applicable if \code{dimension}=2
#' @param colorbar.name The name of the colorbar for quantitative independent variable, \code{PF$X}. If independent variable is not a factor, \code{phylofactor.visualization} assumes they are quantitative and will make a colorbar on the sorted independent variable.
#' @param colorbar.ticks Number of ticks in colorbar
#' @param ... Additional arguments. If \code{dimension==2}, these are arguments for \code{\link{plot}}. Otherwise, they are additional arguments for \code{\link{scatterplot3d}}
#' @examples
#' data(FTmicrobiome)
#' phylofactor.visualize(FTmicrobiome$PF)
#' phylofactor.visualize(FTmicrobiome$PF,dimension=3)
#'
#' pf <- FTmicrobiome$PF
#' pf$X <- 1:40
#' phylofactor.visualize(pf,dimension=3,colorbar.name='sample number')

phylofactor.visualize <- function(PF,dimension=2,default.colors=T,Legend=T,lx=NULL,ly=NULL,colorbar.name=NULL,colorbar.ticks=5,...){

  if (!(dimension %in% c(2,3))){stop('dimension must be either 2 or 3')}

  PROJ <- phylofactor.projection(PF,nfactors=dimension)

  if (default.colors){
    if (is.factor(PF$X)){
      xfactor=T
    ## this is my attempt to guess colors of factors
      n <- length(unique(PF$X))
      cols <- rainbow(n)[match(PF$X,sort(unique(PF$X)))]
      clg <- cols[match(sort(unique(PF$X)),PF$X)]
    } else {
      xfactor=F
      n <- length(unique(PF$X))
      colfunc <- colorRampPalette(c('blue','red','yellow'))(n)
      cols <- colfunc[match(PF$X,sort(unique(PF$X)))]
      clg <- cols[match(sort(unique(PF$X)),PF$X)]

    }
  }

  if (exists('xlab')==F){xlab='PF 1'}
  if (exists('ylab')==F){ylab='PF 2'}
  if (exists('zlab')==F){zlab='PF 3'}
  if (exists('pch')==F){pch=18}
  if (exists('cex')==F){cex=2}

  par(mfrow=c(1,1))
  par(mar=c(5,4,4,1))
  if (exists('xlim') || exists('ylim') || exists('zlim')==FALSE){
    lms <- c(min(PROJ),max(PROJ))
    if(dimension==2){
      if (default.colors){
        if (xfactor){
          plot(PROJ[1,],PROJ[2,],xlab=xlab,ylab=ylab,col=cols,pch=pch,cex=cex,...)
        } else {
          #including colorbar legend is more complicated...

          layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
          plot(PROJ[1,],PROJ[2,],xlab=xlab,ylab=ylab,col=cols,pch=pch,cex=cex,...)

          legend_image <- as.raster(matrix(colfunc, ncol=1))
          if (is.null(colorbar.name)){
            colorbar.name='Legend'
          }
          if (is.null(colorbar.seq)){
            if (is.null(colorbar.ticks)){
              colorbar.ticks = 5
            }
            colorbar.seq=seq(min(PF$X),max(PF$X),l=colorbar.ticks)
          }

          plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend')
          text(x=1.5, y = seq(0,1,l=colorbar.ticks), labels = colorbar.seq)
          rasterImage(legend_image, 0, 0, 1,1)
        }
      } else {
        plot(PROJ[1,],PROJ[2,],xlab=xlab,ylab=ylab,pch=pch,cex=cex,...)
      }
    } else {
      if (default.colors){
        if (xfactor){
        s3d <- scatterplot3d(PROJ[1,],PROJ[2,],PROJ[3,],xlab=xlab,ylab=ylab,zlab=zlab,color = cols,pch=pch,cex.symbols = cex,...)
        } else {
          layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
          s3d <- scatterplot3d(PROJ[1,],PROJ[2,],PROJ[3,],xlab=xlab,ylab=ylab,zlab=zlab,color = cols,pch=pch,cex.symbols = cex)
          legend_image <- as.raster(matrix(colfunc, ncol=1))
          if (is.null(colorbar.name)){
            colorbar.name='Legend'
          }
          if (is.null(colorbar.seq)){
            if (is.null(colorbar.ticks)){
              colorbar.ticks = 5
            }
            colorbar.seq=seq(min(PF$X),max(PF$X),l=colorbar.ticks)
          }

          plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend')
          text(x=1.5, y = seq(0,1,l=colorbar.ticks), labels = colorbar.seq)
          rasterImage(legend_image, 0, 0, 1,1)

        }
      } else {
        s3d <- scatterplot3d(PROJ[1,],PROJ[2,],PROJ[3,],xlab=xlab,ylab=ylab,zlab=zlab,pch=pch,cex.symbols = cex,...)
      }
    }
  }

  if (Legend==T && default.colors && xfactor){
    lbls <- PF$X %>% unique %>% as.list %>% sapply(.,toString,simplify=T)
    if (is.null(lx)){
      if (dimension==2){
        lx <- (max(PROJ[1,])+min(PROJ[1,]))*.2
      } else {
        lx <- max(PROJ[1,])*.8
        lx <- s3d$xyz.convert(x=min(PROJ[1,])-.2*(max(PROJ[1,])-min(PROJ[1,])),y=min(PROJ[2,]),z=max(PROJ[3,])*2)
      }
    }
    if (is.null(ly)){
      if (dimension==2){
        ly <- (max(PROJ[2,])+min(PROJ[2,]))*.8
      }
    }

    legend(x=lx,y=ly,pch=pch,cex=cex,col=clg,legend=lbls)
  }

}
