\name{dyebias.trendplot}
\alias{dyebias.trendplot}

\title{Creates a trend-plot of all reporters, binned by dye bias,
  with all slides ordered by slide bias. 
}

\description{

  The aim of this routine is to show the monotonicity of the total dye
  bias in the (uncorrected) data set. This is to judge whether the total
  dye bias of one reporter in one hybridization indeed behaves as the
  product of an intrinsic gene specific dye bias (iGSDB) and a slide
  specific factor (the slide bias), which is at the heart of the GASSCO method.

  Showing the total dye bias of all reporters is too overwhelming,
  therefore the medians of the total dye bias after binning by intrinsic
  gene specific dye bias (as given in \code{dyebias$dyebias}) are
  plotted.
}
 

\usage{
dyebias.trendplot(data, iGSDBs, dyebias.percentile=5,
                  application.subset=TRUE, n.bins=20, order,
                  output=NULL, ylim=c(-1,1), cex=0.3, lty=1, lwd=1,
                  type="median", main="dye bias trend plot",
                  xlab="slide bias rank", ylab="M", sub=NULL, \dots)
}

\arguments{

  \item{data}{The \code{marrayNorm} to trendplot. 
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

  \item{application.subset}{
    The set of reporters that was eligible for dye bias correction; same
    argument as for \code{\link{dyebias.apply.correction}}.
  }

  \item{n.bins}{The number of bins into which to classify the reporters,
  based on their intrinsic gene-specific dye bias. The median of each
  bin is plotted.}

  \item{type}{What to print for each bin and hybridization. Valid values are:
    \itemize{
      \item{median}{Plot the bin-median: the
        median of all genes in that bin and hybridization}
      \item{mean}{Plot the mean of all genes in that bin and hybridization}
      \item{gene}{Plot the single gene that has the lowest squared distance to the
        bin-median over \emph{all} hybridizations}
      \item{median.gene}{Plot the single gene that has the median
        squared distance to the bin-medians over \emph{all} hybridizations}
      \item{worst.gene}{Plot the single gene that correlates
        worst with the bin-median over \emph{all} hybridizations}
    }
  }

  \item{order}{
    If \code{order==NULL}, the slides are sorted by
    increasing slide bias prior to boxplotting. This is typically done
    for data that is not yet dye bias corrected.  This order is also
    returned as a value. If \code{order!=NULL}, the slides are put into
    this order before trendplotting.  This is typically done for a dye
    bias-corrected data set, using the order of the uncorrected set. }

  \item{output}{ Specifies the output. If \code{NULL}, the existing
    output device is used; if \code{output} is one of  \code{"X11",
      "windows", "quartz"}, 
    a new X11 (Unix)/windows (Windows)/quartz (Mac) device is created.
    If \code{output} is a string ending in one of \code{".pdf", ".png",
    ".eps", ".ps"} is given, a file of that name and type is created and
    closed afterwards.}

  \item{ylim, lty, lwd, main, sub, cex, xlab, ylab}{As for \code{matplot()}}

  \item{\dots}{Other arguments are passed on to \code{matplot()}. }
}

\value{
  The order obtained, for use in a later call to this same function.
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

  ## show trend plots of uncorrected and corrected next to each other:
  ylim <- c(-0.6, 0.6)

  layout(matrix(1:2, nrow=1,ncol=2))

  order <- dyebias.trendplot(data=data.norm, 
                        iGSDBs=iGSDBs.estimated, # from e.g. dyebias.estimate.iGSDBs
                        order=NULL,              # i.e., order by increasing slide bias
                        output=NULL,
                        main="before correction",
                        ylim=ylim
                        )

  order <- dyebias.trendplot(data=correction$data.corrected, # from dyebias.apply.correction
                        iGSDBs=iGSDBs.estimated,
                        order=order,             # order by the original slide bias
                        output=NULL,
                        main="after correction",
                        ylim=ylim
                       ) 

}

\author{Philip Lijnzaad \email{p.lijnzaad@umcutrecht.nl} }

\seealso{
  \code{\link{dyebias.estimate.iGSDBs}},
  \code{\link{dyebias.apply.correction}},
  \code{\link{dyebias.rgplot}},
  \code{\link{dyebias.maplot}},
  \code{\link{dyebias.monotonicity}}
  \code{\link{dyebias.monotonicityplot}}
}


\references{
   Margaritis, T., Lijnzaad, P., van Leenen, D., Bouwmeester, D.,
  Kemmeren, P., van Hooff, S.R and Holstege, F.C.P. (2009).
  Adaptable gene-specific dye bias correction for two-channel DNA microarrays.
  \emph{Molecular Systems Biology}, 5:266, 2009. doi: 10.1038/msb.2009.21.

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}                   % silly, but one keyword is compulsary
\keyword{misc}                    % silly, but one keyword is compulsary
