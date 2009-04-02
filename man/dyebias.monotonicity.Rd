\name{dyebias.monotonicity}
\alias{dyebias.monotonicity}

\title{
  Calculate the degree of monotonicity of the dye bias across the slides.
}

\description{
  If you order genes by their iGSDB, and hybridizations by slide
  bias, the graphs of each gene should form a 'fan' out of the origin
  (see also dyebias.trendplot). This function gives measure of the
  extent to which this is true.
}

\usage{
  dyebias.monotonicity(data,
                       iGSDBs, dyebias.percentile = 5,
                       order = NULL)
}

\arguments{

  \item{data}{The \code{marrayNorm} to determine the degree of
  monotonicity for
    
  }

  \item{iGSDBs}{A data frame with intrinsic gene-specific dye biases,
    the same as that used in \code{\link{dyebias.apply.correction}},
    probably returned by \cr
    \code{\link{dyebias.estimate.iGSDBs}}; see
    there for documentation.
    
  }

  \item{dyebias.percentile}{
    The percentile of intrinsic gene specific dye biases (iGSDBs) for
    which to highlight the reporters. Default should suffice in almost
    all cases.
  }

  \item{order}{
    If \code{order==NULL}, the slides are sorted by
    increasing slide bias prior to boxplotting. This is typically done
    for data that is not yet dye bias corrected.  (This order is also
    returned in the object returned). If \code{order!=NULL}, the slides are put
    into this order first.  This is typically done for a dye
    bias-corrected data set, using the order of the uncorrected set.
  }
}

\details{

  The total dye bias appears to be the product of iGSDB and slide
  bias. In other words, it is monotonous (always increasing or always
  decreasing), both with respect to the intrinsic gene specific dye bias
  and with respect to the slide bias. This function orders genes by
  their iGSDB and the slides by slide bias. Subsequently a linear
  regression of each gene is done, with x being the slide bias rank, (not
  the slide bias itself), and y being the M. The slopes of each linear
  regression line should be an increasing array of values, representing
  the 'fan' of lines. The degree to which this array is increasing is
  tested using the Mann-Kendall test, and is returned. In the case of
  uncorrected data, tau is generally larger than 0.3. After correction,
  tau should be close to zero.

}

\value{

  A \code{dyebias.monotonicity} uses \code{cor.test}, which returns
  \code{htest} object. To this list an extra element, \code{order}, is
  added, which indicates the ordering of the data set by slide bias.
  The degree of monotonicity is indicated by the \code{estimate}
  element; its signficance by the \code{p.value} element.
}

\note{

  This function can take very long to compute, since it calculates
  regressions for each gene.

}


\references{
   Margaritis, T., Lijnzaad, P., van~Leenen, D., Bouwmeester, D.,
  Kemmeren, P., van~Hooff, S.R and Holstege, F.C.P. (2009).
  Adaptable gene-specific dye bias correction for two-channel DNA microarrays.
  \emph{Molecular Systems Biology, submitted}
}

\author{Philip Lijnzaad \email{p.lijnzaad@umcutrecht.nl} }

\seealso{
  \code{\link{dyebias.trendplot}},
  \code{\link{dyebias.monotonicityplot}}
}

\examples{
  \dontshow{
     options(stringsAsFactors = FALSE)

     library(dyebias)
     library(dyebiasexamples)
     data(data.raw)
     data(data.norm)

     ### obtain estimate for the iGSDBs:
     iGSDBs.estimated <- dyebias.estimate.iGSDBs(data.norm,
                                                 is.balanced=TRUE,
                                                 verbose=FALSE)

     ### choose the estimators and which spots to correct:
     estimator.subset <- dyebias.umcu.proper.estimators(maInfo(maGnames(data.norm)))

     application.subset <- maW(data.norm) == 1 &
                   dyebias.application.subset(data.raw=data.raw, use.background=TRUE)

     ### do the correction:
     correction <- dyebias.apply.correction(data.norm=data.norm,
                                            iGSDBs = iGSDBs.estimated,
                                            estimator.subset=estimator.subset,
                                            application.subset = application.subset,
                                            verbose=FALSE)
  }                                     % dontshow

  cat("monotonicity before correction")
  monotonicity <- dyebias.monotonicity(data=data.norm,
                                       iGSDBs=iGSDBs.estimated, 
                                       order=NULL)
  monotonicity

  cat("monotonicity after correction")

  dyebias.monotonicity(data=correction$data.corrected,
                                       iGSDBs=iGSDBs.estimated, 
                                       order= monotonicity$order)


}

\keyword{misc}
\keyword{htest}
