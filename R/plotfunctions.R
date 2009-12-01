## Plot functions

dyebias.rgplot <- function
(data, slide,
 iGSDBs,
 dyebias.percentile=5,
 application.subset=TRUE,
 output=NULL,
 xlim = c(log2(50),log2(50000)),
 ylim = c(log2(50),log2(50000)),
 xticks = c(100,1000,10000,10000),
 yticks = c(100,1000,10000,10000),
 pch = 19,
 cex = 0.3,
 cex.lab = 1.4,
  ...) { 
  file.device <- .set.output(output)
  
  
  dyebias <- .merge.dyebias(data=data, iGSDBs=iGSDBs)
  
  subset <- if(length(application.subset)==1) T else application.subset[,slide]

  lR <- ( maA(data) + maM(data)/2 )
  lG <-  maA(data) - maM(data)/2 

  plot(lR[subset, slide] ~ lG[subset,slide],
       pch=pch, cex=cex, cex.lab=cex.lab,
       xlab="Green", ylab="Red",
       axes=FALSE,
       xlim=xlim, ylim=ylim, ...)

  axis(1, xticks, main="Green", at= log2(xticks))
  axis(2, yticks, main="Red", at= log2(yticks))
  box()
  abline(log2(2),1,col="grey")
  abline(0,1,col="black")
  abline(-log2(2),1,col="grey")
  
  rg.subsets <- .rg.subsets(dyebias, subset, dyebias.percentile)
  green.subset <- rg.subsets$green
  red.subset <- rg.subsets$red
  
  points(lR[green.subset, slide] ~ lG[green.subset, slide], pch=pch, cex=cex, col="green")
  points(lR[red.subset, slide] ~ lG[red.subset, slide], pch=pch, cex=cex, col="red")

  ##  points(lR[,slide] ~ lG[,slide], pch=pch, cex=cex) ## overplot the rest for maximum effect :-)

  if(file.device) {
    dev.off()
  }
} ##dyebias.rgplot

dyebias.maplot <- function
(data, slide,
 iGSDBs, dyebias.percentile=5,
 application.subset=TRUE,
 output=NULL,
 xlim = c(6,16),
 ylim = c(-2,2),
 pch = 19,
 cex = 0.3,
 cex.lab = 1.4,
 ...) { 
  ### as rgplot.it, but now giving an MA plot
  file.device <- .set.output(output)
  

  dyebias <- .merge.dyebias(data=data,iGSDBs=iGSDBs)
  
  subset <- if(length(application.subset)==1) T else application.subset[,slide]

  M=maM(data)
  A=maA(data)
  
  plot(x=A[subset,slide],y=M[subset, slide], 
       pch=pch, cex=cex, cex.lab=cex.lab,
       xlab="A", ylab="M",
       xlim=xlim, ylim=ylim
       , ...)

  box()
  abline(0,0,col="black")
  abline(1,0,col="grey")
  abline(-1,0,col="grey")
  
  rg.subsets <- .rg.subsets(dyebias, subset, dyebias.percentile)
  green.subset <- rg.subsets$green
  red.subset <- rg.subsets$red
  
  points(M[green.subset, slide] ~ A[green.subset, slide], pch=pch, cex=cex, col="green")
  points(M[red.subset, slide] ~ A[red.subset, slide], pch=pch, cex=cex, col="red")

  ##  points(lR[,slide] ~ lG[,slide], pch=pch, cex=cex) ## overplot the rest for maximum effect :-)

  if(file.device) {
    dev.off()
  }

} ##dyebias.maplot

dyebias.boxplot <- function(data, iGSDBs, dyebias.percentile=5, application.subset=TRUE,
                       order, output=NULL,
                       ylim=c(-4,4),
                       ...) { 
  file.device <- .set.output(output)

  dyebias <- .merge.dyebias(data=data,iGSDBs=iGSDBs)
  
  slide.bias <- .slide.bias(data, dyebias, dyebias.percentile, application.subset=application.subset)
  reds <- slide.bias[["red"]]
  greens <- slide.bias[["green"]]
    
  ## order the slides by increasing dyebias:
  if (is.null(order)) { 
    order <- .slide.bias.order(slide.bias)
  } else {                              #re-use previous ordering
    if (is.logical(order) && !order) {
      order <- 1:maNsamples(data)       # don't order. NOTE maNsamples is wrong name
    }
  }

  greens.ordered <- list()
  reds.ordered <- list()

  for (i in 1:maNsamples(data)) {       # NOTE: maNsamples is wrong name
    greens.ordered[[i]] <- greens[[ order[i] ]]
      reds.ordered[[i]] <- reds [[ order[i] ]]
  }
  names(greens.ordered) <- names(greens)[order]
  names(reds.ordered) <- names(reds)[order]
  
  par(pin=c(5.2,5.2), mar=c(4.2,4.2,2.1,2.1))

  boxplot(greens.ordered, col="green", outpch=NA, whiskcol="black", whisklty=1,ylim=ylim,
          cex=0.6, cex.lab=1.4, ylab="slide bias (M)", xlab="slide", na.rm=TRUE, ...)

  abline(h=0,col="black")

  boxplot(reds.ordered, col="red", outpch=NA, whiskcol="black", whisklty=1, na.rm=TRUE,  add=TRUE, ...)

  if(file.device) {
    dev.off()
  }

  return(order)                         #so the next one can use the same ordering
}                                       # dyebias.boxplot

dyebias.trendplot <- function(data,
                              iGSDBs,
                              dyebias.percentile=5,
                              application.subset=TRUE,
                              n.bins=20,
                              order, output=NULL,
                              ylim=c(-1,1),
                              cex=0.3,
                              lty=1,
                              lwd=1,
                              type="median",
                              main="dye bias trend plot",
                              xlab="slide bias rank",
                              ylab="M",
                              sub=NULL, #title, default should suffice
                              ...) { 

  plot=list(median=T, mean=T, gene=T, median.gene=T, worst.gene=T) # what to plot, in the end
  
  if(is.null(plot[[type]])) {
    stop("plotfunctions.R: type '", type ,"' not recognized" , call.=T)
  }
  
  all.bins <- seq(0,1, by=1/n.bins)
  all.quantiles <- quantile(iGSDBs$dyebias, all.bins, na.rm=TRUE)

  file.device <- .set.output(output)

  dyebias <- .merge.dyebias(data=data,iGSDBs=iGSDBs)

  slide.bias <- .slide.bias(data, dyebias, dyebias.percentile, application.subset=application.subset)
    
  ## order the slides by increasing dyebias:
  if (is.null(order)) { 
    order <- .slide.bias.order(slide.bias)
  }

  x <- as.array(1:maNsamples(data))     # maNsamples is wrong name
  y <- apply(x, 1, function(i){maM(data)[,order[i]]})

### do the averaging etc. per bin:
  m <- matrix(0, nrow=n.bins, ncol=maNsamples(data)) # maNsamples(data) is wrong name ...
  plot$median <- m
  plot$mean <- m
  plot$gene <- m                        #meaning: one closest to bin-median
  plot$median.gene <- m
  plot$worst.gene <- m

  for (i in 1:n.bins) {
    min <- all.quantiles[i]
    max <- all.quantiles[i+1]

    slice <- y[ min <= dyebias$dyebias  & dyebias$dyebias < max, ]

    plot$median[i,] <- apply( slice , 2, median, na.rm=TRUE)
    plot$mean[i,] <- apply( slice , 2, mean, na.rm=TRUE)

    ## following takes time, so only do if needed:
    if(length(grep("gene", type))) {
      ## best gene has lowest sum of squares with median, across all slides:
      sum.sq <- apply( as.array(1:nrow(slice)), 1, function(row) { sum((plot$median[i,] - slice[row,])^2)})
      ranking <- order(sum.sq)
      plot$gene[i,] <- slice[ ranking[1], ]

      # median gene (in terms of distance to the median of the bins)
      plot$median.gene[i,] <- slice[ ranking[round(length(ranking)/2)], ]
      
      ## worst gene has lowest correlation with median, across all slides:
      cor <- apply( as.array(1:nrow(slice)), 1, function(row) { cor(plot$median[i,], slice[row,])})
      plot$worst.gene[i,] <- slice[ order(cor)[1], ]
    }

  }                                     # n.bins

  mid <- length(all.quantiles[ all.quantiles <= 0])
  if(mid ==0) {
    mid <- 1                   # happens when all iGSDB > 0. In that case, just draw
  }  

  red <- mid:n.bins
  
  if (is.null(sub)) {
    sprintf("trend (type '%s') of dye bias per slide, binned by dye bias (%d bins)",
            type, n.bins)
  }

  ## first plot all greens, then overplot with red (if needed)
  ## reason is that dimension becomes funny if there's is one green line.
  matplot(x, t((plot[[type]])[,]), pch=19, type="l", col="green", cex=cex,
          xlab=xlab, ylab=ylab, ylim=ylim, main=main,
          sub=sub,
          axes=FALSE,
          lty=lty, lwd=lwd, ...)

  axis(side=1, at=1:maNsamples(data),labels=1:maNsamples(data))
  axis(side=2)

  if(mid < n.bins) {                   # i.e. unless all iGSDBs < 0
    matlines(x, t((plot[[type]])[red,]), pch=19, type="l", col="red",
             cex=cex, lty=lty,lwd=lwd)
  }
  abline(h=0,col="black")

  if(file.device) {
    dev.off()
  }

  return(order)                         #so the next one can use the same ordering
}                                       #dyebias.trendplot

dyebias.monotonicity <- function(data,
                                 iGSDBs,
                                 dyebias.percentile=5,
                                 order=NULL) { 

  slopes <- .monotonicity(data,iGSDBs, dyebias.percentile, order)



  test <- .mann.kendall(slopes)
  test$order <- order            # so the next one can use the same ordering 
  return(test)

}                                       # dyebias.monotonicity

dyebias.monotonicityplot <- function(data,
                                     iGSDBs,
                                     dyebias.percentile=5,
                                     order=NULL,
                                     output=NULL,
                                     pch=19, cex=0.3,
                                     cex.lab=1.4,
                                     ylim=c(-0.2, 0.2),
                                     xlab="rank", ylab="slope",
                                     sub=NULL,
                                     ...
                                     ) { 

  file.device <- .set.output(output)

  slopes <- .monotonicity(data,iGSDBs, dyebias.percentile, order)

  if(is.null(sub)) {
    test <- .mann.kendall(slopes)
    sub=sprintf("tau=%.3g, p=%.3g", test$estimate, test$p.value) 
  }
  
  n <- length(slopes)
  green <- floor(n*dyebias.percentile*0.01)
  red <- floor(n*(1-dyebias.percentile*0.01))

  black <- (green+1):(red-1)
  green <- 1:green
  red <- red:n

  plot(x=black, y=slopes[black],
       xlab=xlab, ylab=ylab,cex.lab=cex.lab, pch=pch, cex=cex,ylim=ylim, col="black", sub=sub, ...)

  points(x=green, y=slopes[green], pch=pch, cex=cex, col="green")
  points(x=red,   y=slopes[red],   pch=pch,cex=cex, col="red")
  
  if(file.device) {
    dev.off()
  }
  return(order)
}                                       #dyebias.monotonicityplot


## ===== support functions ================================================


.set.output <- function(output) {
  ## set output device based on file name extension (if any). Return T if
  ## it's a file (so you know if a call to dev.off() is needed afterwards)
  
  if(is.null(output))  {                #reuse existing output channel, or start one
    return(FALSE)
  }

  if(output=="")  {                     #reuse existing channel, or start one
    return(FALSE)
  }

  if (output=="X11") {                   # meaning: new X11 display
    X11()
    return(F)
  }

  if (output=="windows") {               # meaning: new Windows display
    windows()
    return(F)
  }

  if (output=="quartz") {                # meaning: new Quartz display
    quartz()
    return(F)
  }

  if (length(grep("\\.pdf$", output)) > 0) {
    pdf(output)
    return(T)
  }

  if (length(grep("\\.png$", output)) > 0) {
    png(output)
    return(T)
  }

  if ( (length(grep("\\.eps$", output)) > 0)
      || (length(grep("\\.ps$", output)) > 0)) {
    postscript(output)
    return(T)
  }
  stop("plotfunctions.R: .set.output: Unknown extension. Only know .pdf .png .eps .ps \n", call. =TRUE)
}                                       #.set.output

.merge.dyebias <- function(data, iGSDBs) {
  ## put dyebias in right shape: one value per spot, in order of the data
  .check.required.columns(frame=iGSDBs,
                          wanted.columns=c("reporterId", "dyebias", "A"),
                          error.suffix=".merge.dyebias")
  
  ## merge dyebias into it, for selecting the percentiles
  gnames <- data.frame(reporterId=maLabels(maGnames(data)))
  gnames$dyebias <- NULL # may already be part of it, gives dyebias.{x,y} name after merge
  gnames$A <- NULL

  gnames$rank <- 1:length(gnames[[1]])       #for sorting back
  dyebias <- merge(gnames, iGSDBs, by="reporterId", all.x=TRUE)
  dyebias <- dyebias[ order(dyebias$rank), ]
  dyebias[ is.na(dyebias$dyebias), "dyebias"] <- 0.0 ## genes for which we don't have an estimate get zero bias
  dyebias[ is.na(dyebias$A), "A"] <- 666 ## an evil A :-)
  dyebias
}                                       # .merge.dyebias

.rg.subsets <- function(dyebias, subset, dyebias.percentile) {
  # return list(red,green) of indexes of the top dyebias.percentile red and green biased genes within subset
  large.dyebias.cutoff.green <- quantile( dyebias[ subset , "dyebias"], probs=0.01*dyebias.percentile, na.rm=TRUE) 
  large.dyebias.cutoff.red <- quantile( dyebias[ subset , "dyebias"], probs= (1- 0.01*dyebias.percentile), na.rm=TRUE) 
  
  large.dyebias.green.ids <- dyebias[ (subset & dyebias$dyebias < large.dyebias.cutoff.green), "reporterId"]
  large.dyebias.red.ids <- dyebias[ (subset & dyebias$dyebias > large.dyebias.cutoff.red), "reporterId"]
  
  list(green=dyebias$reporterId %in% large.dyebias.green.ids, red=dyebias$reporterId %in% large.dyebias.red.ids)
}                                       #.rg.subsets


.slide.bias <- function(data, dyebias, dyebias.percentile, application.subset=TRUE) {
  ## return a list(red, green) with the M's of the top dyebias.percentile red
  ## and green biases genes (this is used to order the slides by slide bias)

  reds <- list()
  greens <- list()

  for (slide in 1:maNsamples(data))  {
    the.subset <- if (length(application.subset)==1) T else application.subset[,slide]

    rg.subsets <- .rg.subsets(dyebias, the.subset, dyebias.percentile)
    green.subset <- rg.subsets$green
    red.subset <- rg.subsets$red
    
    slide.name <- maLabels(maTargets(data))[slide]

    greens[[slide.name]] <- maM(data)[green.subset, slide]
    reds[[slide.name]] <- maM(data)[red.subset, slide]
  }

  list(red=reds, green=greens)
}                                       # .slide.bias

.slide.bias.order <- function(slide.bias) {
  ## return ordering for slides by bias
  order( sapply(slide.bias$red, median, na.rm=TRUE)
        - sapply(slide.bias$green, median, na.rm=TRUE))
}                                       # .slide.bias.order


.monotonicity <- function(data,
                          iGSDBs,
                          dyebias.percentile=5,
                          order=NULL) { 
  ## alternatives would be Page's trend test or Jonckheere-Terpstra test

  ## only meaning full for uncorrected data, really
  dyebias <- .merge.dyebias(data=data,iGSDBs=iGSDBs)

  slide.bias <- .slide.bias(data, dyebias, dyebias.percentile,
                            application.subset=TRUE)
    
  ## order the slides by increasing dyebias:
  if (is.null(order)) { 
    order <- .slide.bias.order(slide.bias)
  }

  x <- as.array(1:maNsamples(data))     # maNsamples is wrong name
  y <- apply(x, 1, function(i){maM(data)[order(dyebias$dyebias) ,order[i]]})
  ## y is now one gene per row, rows sorted by iGSDB, columns sorted by slide bias
  
  # insist on at least 4 proper values:
  y <- y [apply(y, 1, function(row){sum(!is.na(row))})>=4, ]

  slope <- function(data) {
    lsfit(x=1:length(data), y=data, intercept=TRUE)$coefficients["X"]
  }

  slopes <- apply(y, 1, function(row){slope(row)})
  ## alternative is to use mann.kendall for the 'slope' (as well as for tau), like
  ##   slopes <- apply(y, 1, function(row){mann.kendall(row)$estimate})
  ## results not as good though

  slopes
}                                       # .monotonicity

.mann.kendall <- function(x) {
  cor.test(1:length(x), x, method="kendall")
}
