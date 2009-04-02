###
### dye bias normalization according to Margaritis et al. 2009
###
###
### Main routines. All functions starting with "dyebias." are exported
###
###

dyebias.estimate.iGSDBs <- function
(data.norm,
 is.balanced=TRUE,
 reference="ref",
 verbose=FALSE) {

  averaging.function <- mean               # to use when is.balanced==TRUE

  if (verbose) {
    print("Estimating the intrinsic gene-specific dyebiases ...")
  }

  here=", in dyebias.R:dyebias.estimate.iGSDBs"
  
  if (class(data.norm) != "marrayNorm") {
    stop("Need marrayNorm object", here, call.= TRUE)
  }

  data.norm.sorted <- data.norm[order(maLabels(maGnames(data.norm))), ]
  
  n.spots <- maNspots(data.norm.sorted)
  ## n.slides <- maNtargets(data.norm.sorted)  ## not available in marray ...
  n.slides <- maNsamples(data.norm.sorted) ## wrong name! we have 5 hybs, so 10 samples!

  
  ## for code using reporter.filter as a filter on the oligos for which to find iGSDBs, see  before revision 2890

  if (verbose) {
    print(sprintf("Found %d slides containing %d spots", n.slides, n.spots))
  }

  ## average the replicate spots. (For version not averaging the replicate
  ## spots, see prior to rev. 2885)

  reporter.ids <- maLabels(maGnames(data.norm.sorted))
  if(length(reporter.ids)==0) {
    stop("No reporter IDs found (has maLabels been set?))",here, call. = TRUE)
  }

  unique.ids <- unique(reporter.ids)
  n.unique <- length(unique.ids)

  if (n.unique == length(reporter.ids)) {
    Mavg <- maM(data.norm.sorted)
    Aavg <- maA(data.norm.sorted)
  } else { ## need to average the spots (which takes time, so only do if needed)
    if(verbose) {
      print("Averaging replicate spots ...")
    }
    
    Mavg <- matrix(0, n.unique, n.slides)
    Aavg <- matrix(0, n.unique, n.slides)

    for (i in 1:n.slides) { 
      Mavg[, i] <- .average.replicates( data.frame(id=reporter.ids, value=maM(data.norm.sorted)[,i], stringsAsFactors=FALSE))
      Aavg[, i] <- .average.replicates( data.frame(id=reporter.ids, value=maA(data.norm.sorted)[,i], stringsAsFactors=FALSE))
    }

    if(verbose) {
      print("Done averaging replicate spots")
    }
  }

  if (is.balanced) {                    #we can just take average; much faster
    ### note that the reference argument is not used
    if (verbose) { print("Balanced design, simply averaging") }
    return(data.frame(reporterId=unique.ids,
                      dyebias=apply(Mavg, 1, averaging.function, na.rm=TRUE),
                      A=apply(Aavg,1, averaging.function, na.rm=TRUE),
                      stringsAsFactors=FALSE
                      ))
  }

  ### else: do full limma model
  library(limma)
  targets <- maInfo(maTargets(data.norm.sorted))
  .check.required.columns(targets,
                          wanted.columns=c("Cy3","Cy5"),
                          error.suffix=here)

  design <- .set.design(targets, reference, verbose)
  if(verbose) {
    print("Design\n");
    print(design);
  }

  data.limma <- new("MAList", list(M=Mavg, A=Aavg, weights=matrix(1, n.unique, n.slides),
                                   printer=list(ndups=1, spacing=1),
                                   genes=data.frame(reporterId=unique.ids, Control="gene"),
                                   targets=targets))
  
  if (verbose) { print("fitting linear model ...")}
  fit <- lmFit(data.limma, design)
  efit <- eBayes(fit)
  results <- .limma.to.dataframe(efit)
  
  if(verbose) {
    print("intrinsic GSDB's:\n")
    print(summary(results$Coef.Dye))    # name 'Dye' comes from design
  }
  
  if(verbose) {
    print("Done estimating the intrinsic gene-specific dyebiases.")
  }
  
  final <- data.frame(reporterId = results$Genes.reporterId, dyebias = results$Coef.Dye, A=results$A,
                      stringsAsFactors=FALSE)
  return(final)
}                                       # dyebias.estimate.iGSDBs

dyebias.apply.correction <-  function
(data.norm,
 iGSDBs,
 estimator.subset=TRUE, 
 application.subset=TRUE, 
 dyebias.percentile=5,
 minmaxA.perc=25,
 minA.abs=NULL,                         
 maxA.abs=NULL,
 verbose=FALSE
 ) {  
  here <- ", in dyebias.R:dyebias.apply.correction"

  
#  if (verbose) {
#    print("Invoking dye bias correction. Please cite Margaritis et al. 2008 etc.\n")
#  }

  .check.required.columns(frame=iGSDBs,
                          wanted.columns=c("reporterId", "dyebias", "A"),
                          error.suffix=here)
  
  if(any(duplicated(iGSDBs$reporterId))) {
    stop("table of intrinsic GSDB's contains duplicates, which confuses me greatly", here)
  }

  if(length(application.subset)==1) {
    application.subset <- maM(data.norm)
    application.subset[,] <- TRUE
  } else if ( length(application.subset) == maNspots(data.norm) )   {
    application.subset <- rep(application.subset, maNsamples(data.norm))
    dim(application.subset) <- dim(maM(data.norm))
  } else if ( !all(dim(application.subset) == dim(maM(data.norm)) ))  {
    stop("wrong dimensions for application.subset argument")
  }

  ##  reporter.info <- maInfo(maGnames(data.norm))
  reporter.info <- data.frame(reporterId=maLabels(maGnames(data.norm)))

  ## put (replace) the dyebiases into reporter.info. If we don't have an iGSDB, make it zero
  reporter.info$dyebias <- NULL
  reporter.info$rank <- 1:length(reporter.info[[1]])
  reporter.info <- merge(reporter.info, iGSDBs, by="reporterId", all.x=TRUE)
  reporter.info [ is.na(reporter.info$dyebias), "dyebias"] <- 0.00
  reporter.info <- reporter.info[order(reporter.info$rank),] # sort back to original order
  reporter.info$rank <- NULL


  estimators <- .find.estimators(reporter.info=reporter.info,
                                 estimator.subset=estimator.subset,
                                 dyebias.percentile=dyebias.percentile,
                                 minA.abs=minA.abs, 
                                 maxA.abs=maxA.abs,
                                 minmaxA.perc=minmaxA.perc, 
                                 verbose=verbose)
  n.slides <- length(maInfo(maTargets(data.norm))[[1]])
##  slide.names <- maLabels(maTargets(data.norm))

  n.spots <- maNspots(data.norm)

  Mcor <- matrix(0, n.spots, n.slides)
  Rcor <- matrix(0, n.spots, n.slides)
  Gcor <- matrix(0, n.spots, n.slides)

  summary <- data.frame(slide="", file="",
                        green.bias="", red.bias="",
                        green.correction="", red.correction="",
                        avg.correction="",
                        var.ratio="", reduction.perc="", p.value="")

  summary.names  <- names(summary) 

  summary  <- summary[ summary$slide != "",] #empty it

  for (i in 1:n.slides) {
    file <- maLabels(maTargets(data.norm))[i]

    if (verbose) {
      print(sprintf("Dyebias correcting slide %d (%s)...\n",i, file))
    }
    
    ## actual dye effects of the bottom/top 5% percentile reporters on this slide:
    slide.bias.green = median( maM(data.norm)[ estimators$green.subset, i], na.rm=TRUE)
    slide.bias.red = median( maM(data.norm)[ estimators$red.subset, i], na.rm=TRUE)
    slide.correction.green = slide.bias.green / estimators$green.effect
    slide.correction.red = slide.bias.red / estimators$red.effect

    slide.correction.avg <- mean( c(slide.correction.red, slide.correction.green ))

    if(verbose) { 
      print(sprintf("Estimated slide bias and correction\n: green bias: %g, red bias: %g, green correction: %g, red correction: %g\nTotal correction %g\n",
                    slide.bias.green, slide.bias.red,
                    slide.correction.green, slide.correction.red,
                    slide.correction.avg))
    }

    excluded <- !application.subset[,i]
    if( length(excluded) > 0) {
      if(verbose) {
        print(sprintf("Excluding %d out of %d spots from dyebias correction\n",
                      sum(excluded), n.spots))
      }
    }

    ## finally, apply the correction:
    Mcor[ !excluded, i ] <- maM(data.norm[ !excluded, i]) -
      as.matrix(reporter.info[ !excluded, "dyebias"] * slide.correction.avg)
    Mcor[  excluded, i ] <- maM(data.norm[excluded,i])

    var.ratio.included <- (sd(Mcor[!excluded,i], na.rm=TRUE) / sd(maM(data.norm)[!excluded,i], na.rm=TRUE))^2
    p.value <- (var.test(Mcor[!excluded, i], maM(data.norm)[!excluded, i]))$p.value
    ## Levene test (more robust with non-Normal distributions) always gives even lower
    ## p-values (!)
    reduction.perc <- 100*(1.0-var.ratio.included)

    if(verbose) {
      print(sprintf("variance ratio (i.e. corrected/uncorrected): %.3g\nOveral reduction in variance is %.1f %% (p=%.3g)\n",
                    var.ratio.included, reduction.perc, p.value))
    }

    summary <- rbind(summary,
                     c(as.character(i),
                       as.character(file),
                       as.character(slide.bias.green),
                       as.character(slide.bias.red),
                       as.character(slide.correction.green),
                       as.character(slide.correction.red),
                       as.character(slide.correction.avg),
                       as.character(var.ratio.included),
                       as.character(reduction.perc),
                       as.character(p.value)))
    
    Rcor[,i] <- 2^(maA(data.norm[,i]) + 0.5*Mcor[,i])
    Gcor[,i] <- 2^(maA(data.norm[,i]) - 0.5*Mcor[,i])

    if (verbose) {
      print(sprintf("Done correcting slide %d (%s)\n",i, file))
    }
  } ## i in 1:nslides

  names(summary) <- summary.names
  
  data.dyecorr <- data.norm

  maM(data.dyecorr) <- Mcor
  maA(data.dyecorr) <- maA(data.norm)
  
  if( .have.umcu.version() ) {          # our own marray version has "maR<-" and "maG<-"
    maR(data.dyecorr) <- Rcor           # functions, which will set the corresponding 
    maG(data.dyecorr) <- Gcor           # slots (absent in the 'real' marray)
  } 

  ##  maLabels(maTargets(data.dyecorr)) <- slide.names

  if(verbose) {
    print("Done dyebias-correcting the slides. Summary of the variance-reduction percentages:\n")
    print(summary(as.numeric(summary$reduction.perc)))
  }
  
  return(list(data.corrected=data.dyecorr, estimators=estimators, summary=summary))
}                                      # dyebias.apply.correction

dyebias.application.subset <- function
(data.raw=NULL,
 min.SNR=1.5,
 use.background=FALSE,
 maxA=15
 ) { 
  here <- ", in utils.R:dyebias.application.weights"

  n.slides <- length(maInfo(maTargets(data.raw))[[1]])

  weights <- matrix(1.0, maNspots(data.raw), maNsamples(data.raw) )

  for (i in 1:n.slides) {

    if(use.background) {

      if  (  all(  dim(maRb(data.raw)) == 0 ) 
           || all(  dim(maGb(data.raw)) == 0 ) ) {
        stop("use.background==TRUE but no weights are available", here)
      }
    
      unmeasurable <- ( ( maRf(data.raw)[,i] / maRb(data.raw)[,i] < min.SNR)
                       | (maGf(data.raw)[,i] / maGb(data.raw)[,i] < min.SNR))
    } else {                            #estimate it from min. foreground:
      minRf <- min(maRf(data.raw)[,i])
      minGf <- min(maGf(data.raw)[,i])
      
      unmeasurable <- (( maRf(data.raw[,i]) / minRf < min.SNR)
                       | (maGf(data.raw[,i]) / minGf < min.SNR))
    }

    saturated <- (maRf(data.raw[,i]) > 2^maxA | maGf(data.raw[,i]) > 2^maxA )
    
    excluded <- unmeasurable | saturated
    
    weights[excluded, i] = 0
  }
  return(weights)
}                           #dyebias.application.subset

## --- support routines, not exported

.check.required.columns <- function
## simple checking function
(frame, wanted.columns, error.suffix) {
  for (col in wanted.columns) {
    if (is.null(frame[, col])) {
      stop(sprintf("column '%s' is missing. Columns found: %s.%s",
                   col, paste(names(frame),collapse=" "), error.suffix), call. = TRUE)
    }
  }
}                           #.check.required.columns

.find.estimators <- function
## returns a list( {green,red}.{ids,cutoff,subset,effect} of 
## reporters (and effects) having the strongest intrinsic gene-specific dye bias.
## These reporters and values are used to estimate the slide bias.
(reporter.info,
 estimator.subset=TRUE,          # subset of reporters to consider as estimators
 dyebias.percentile=5, # The percentile of iGSDB's to use for estimating the slide bias.
 minA.abs=NULL,                         # genes with average A lower than minA.abs or 
 maxA.abs=NULL,                         # higher than maxAbs are never used as slide bias
 minmaxA.perc=25,                       # estimators. If minmaxA.perc is given, a percentile
 verbose=FALSE                              # is used.
 ) {
  here <- ", in: dyebias.R: .find.estimators"

  if(is.null(estimator.subset)) {
    estimator.subset <- TRUE;                        #i.e.all
  } else { 
    
    if( (length(estimator.subset) > 1 ) &&
       length(reporter.info[[1]]) != length(estimator.subset)) {
      stop("expected estimator.subset and reporter.info to be of equal length", here, call. =TRUE)
    }
  }

  if ( !is.null(minA.abs) ) { 
    minA <- minA.abs
  } else {
    minA <- quantile( reporter.info[ estimator.subset , "A"], probs=0.01*minmaxA.perc, na.rm=TRUE )
  }

  if ( !is.null(maxA.abs) ) { 
    maxA <- maxA.abs
  } else {
    maxA <- quantile( reporter.info[ estimator.subset , "A"], probs=(1-0.01*minmaxA.perc), na.rm=TRUE )
  }
  
  if(verbose) {
    print(sprintf("Average log2(signal) for estimator genes: min=%g, max=%g\nExcluding those outside %g < A < %g\n",
                  min(reporter.info[,"A"]), max(reporter.info[,"A"]),
                  minA, maxA))
  }

  rightrange <- (reporter.info$A > minA & reporter.info$A < maxA)
  subset <- ( estimator.subset & rightrange)
  
  large.dyebias.cutoff.green = quantile( reporter.info[ subset , "dyebias"], probs=0.01*dyebias.percentile, na.rm=TRUE ) 
  large.dyebias.cutoff.red = quantile( reporter.info[ subset , "dyebias"], probs= (1- 0.01*dyebias.percentile), na.rm=TRUE) 

  large.dyebias.green.ids = reporter.info[ (subset & reporter.info$dyebias < large.dyebias.cutoff.green), "reporterId"]
  large.dyebias.red.ids = reporter.info[ (subset & reporter.info$dyebias > large.dyebias.cutoff.red), "reporterId"]

  green.subset = reporter.info$reporterId %in% large.dyebias.green.ids
  red.subset = reporter.info$reporterId %in% large.dyebias.red.ids

  if (sum(green.subset) < 10 ) {
    stop("could not find enough green estimator genes", here, call. = TRUE)
  }
  if (sum(red.subset) < 10 ) {
    stop("could not find enough red estimator genes", here, call. = TRUE)
  }

  green.effect = median(reporter.info[ green.subset , "dyebias"])
  red.effect = median(reporter.info[ red.subset, "dyebias"])

  strongest=list(
    green.ids = large.dyebias.green.ids
    , green.cutoff = large.dyebias.cutoff.green
    , green.subset = green.subset
    , green.effect = green.effect
    , red.ids = large.dyebias.red.ids
    , red.cutoff = large.dyebias.cutoff.red
    , red.subset = red.subset
    , red.effect = red.effect
    )

  if(verbose) {
    print(sprintf("Top %g %% reporters (%d spots; %d of them candidates, %d of them in the right range, %d both) used to estimate the slide bias:\n\tgreen:\tcutoff: %g, nprobes: %d (%d spots), median: %g\n\tred:\tcutoff: %g, nprobes: %d (%d spots), median: %g\n",
                  dyebias.percentile, length(reporter.info[[1]]), sum(estimator.subset), sum(rightrange), sum(subset), 
                  strongest$green.cutoff, length(strongest$green.ids), sum(strongest$green.subset), strongest$green.effect, 
                  strongest$red.cutoff,   length(strongest$red.ids),   sum(strongest$red.subset),   strongest$red.effect
                  ))
  }
  return( strongest)
}                                     # .find.estimators

.average.replicates <- function
## average replicate spots within one slide.
(frame,                                 #has to have columns "id" and "value"
 averaging.function=mean) {             #(median is much slower!)
  here <- "dyebias.R: .average.replicates"

  avg.data <- tapply(frame[, "value"], frame[, "id"], averaging.function, na.rm=TRUE)
  ids = rownames(avg.data)
  if (! all(order(ids) ==  1:length(ids))) { 
    stop("averaging mixed up id's", here, call. = TRUE)
  }

  ## mean of a factor level with only NA's becomes NaN, which is bit useless; make them NA's:
  avg.data[ is.nan(avg.data)] <- NA
  return(avg.data)
}                                       #average.replicates

.limma.to.dataframe <- function(fit, results = NULL, digits = 3, adjust = "none") {

  my.is <- function(thing, expected.class){      # is() from methods gives load errors ...
    if(is.null(thing)){return(FALSE);}
    class <- attributes(thing)$class
    if(is.null(class)){return(FALSE)}
    return (class==expected.class)
  }

  if (!my.is(fit, "MArrayLM"))
    stop("fit should be an MArrayLM object")

  if (!is.null(results) && !my.is(results, "TestResults"))
    stop("results should be a TestResults object")

  if (is.null(fit$t) || is.null(fit$p.value))
    fit <- eBayes(fit)
  
  p.value <- as.matrix(fit$p.value)
  for (j in 1:ncol(p.value)) p.value[, j] <- p.adjust(p.value[,j], method = adjust)

  rn <- function(x, digits) {
    if (is.null(x))
      NULL
    else
      round(x, digits=digits)
  }

  ### note: weird way of doing this, but needed to preserve the structure
  tab <- list()
  tab$A <- rn(fit$Amean, digits = digits - 1)

  ## following three 'columns'  (i.e. Coef, t and p.value) have as many 'sub-columns'
  ## as there are effects. So the "mutVSwt" effect will end up in 
  ## tab$Coef.mutVSwt, tab$t.mutVSwt and tab$p.value.mutVSwt, etc.
  tab$Coef <- rn(fit$coefficients, digits = digits)
  tab$t <- rn(fit$t, digits = digits - 1)
  tab$p.value <- rn(p.value, digits = digits + 2)

  ## HOWEVER, if there's just one effect (say mutVSwt), there is an R-(or limma) bug that results in
  ## completely wrong column names (i.e.
  ##   A mutVSwt mutVSwt.1 mutVSwt.2 F F.p.value Genes.Name Genes.rid
  ## rather than
  ##   A p.value.mutVSwt t.mutVSwt p.value.mutVSwt F F.p.value Genes.Name Genes.rid
  ## i.e. they seem to be overwritten #$%^&*). Solve this by redoing the 
  ## assignment explicitly for this case:
  if (ncol(fit$design) == 1) {
    tab$Coef =NULL
    tab$t =NULL
    tab$p.value =NULL

    effect.name=colnames(fit$design)

    colname=paste("Coef",effect.name, sep=".")
    tab[[ colname ]]= as.vector(rn(fit$coefficients, digits = digits))

    colname=paste("t",effect.name, sep=".")
    tab[[ colname ]]= as.vector(rn(fit$t, digits = digits))

    colname=paste("p.value",effect.name, sep=".")
    tab[[ colname ]]= as.vector(rn(fit$p.value, digits = digits))
  }

  ### following is single columns
  tab$F <- rn(fit$F, digits = digits - 1)
  tab$F.p.value <- rn(fit$F.p.value, digits = digits + 2)

  tab$Res <- unclass(results)

  ### Genes consist of sub-columns Name, rid and Controls:
  tab$Genes <- fit$genes
  
  return(data.frame(tab, check.names=FALSE, stringsAsFactor=FALSE));
} ## .limma.to.dataframe

.set.design <- function (targets, reference, verbose) {
  ## finds design to be used
  here=", in dyebias.R: set.design"
  n.slides <- length(targets[[1]])

  if (length(reference)==1) {
    if(verbose) { print("Found common-reference design");} 
    design <- modelMatrix(targets, ref=reference)
    return(cbind(design, Dye=1))
  } else {
    if(verbose) { print("Found set-of-common-references design");} 

    ## see if all refs and non-refs (which have must have ref as a prefix)
    ## sum up properly:
    cy.any <- c(targets$Cy3, targets$Cy5)
    list <- apply(as.array(paste("^",reference, sep="")), 1,
                  function(r){grep(r,cy.any)})
    

    if ( do.call("sum", args=lapply(list, length))  != 2*n.slides) {
      stop("Problem with naming of non-reference samples: must all start with reference name",
           here, call.= TRUE)
    }

    ## (other errors are caught by modelMatrix)
    
    sub.matrices <- lapply( as.list(reference), 
                           function(ref){modelMatrix( targets[grep(sprintf("^%s",ref), targets$Cy3),], ref=ref)}) 
    design <- do.call("blockDiag", args=sub.matrices)
    design <- cbind(design, "Dye" = rep(1,n.slides))
    if(sum(design) != n.slides) {
      stop("Problem with design: expected it to sum up to number of slides", here, call.= TRUE)
    }
    return(design[order(rownames(design)), ]) 
  }
}                                       #set.design

.have.umcu.version <- function() {
  ## tell if we're running our own modified version of marray or not
  ## (not a very precise test, but never mind)
  return(package.version(pkg="marray")=="1.5.8")
}
