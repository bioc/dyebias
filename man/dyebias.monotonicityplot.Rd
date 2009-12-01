\name{dyebias.monotonicityplot}
\alias{dyebias.monotonicityplot}

\title{
  Show the degree of monotonicity of the dye bias across the slides.
}
\description{

  If you order genes by their iGSDB, and hybridizations by slide bias,
  the graphs of each gene should form a 'fan' out of the origin (see
  also dyebias.trendplot). This function plots the regression slope of
  each gene, ordered by iGSDB and slide bias.  If the uncorrected total
  dye bias is indeed monotonous, an increasing trend should be visible.

}

  
\usage{
  dyebias.monotonicityplot(data,
                           iGSDBs,
                           dyebias.percentile = 5,
                           order = NULL,
                           output = NULL,
                           pch = 19, cex = 0.3, cex.lab = 1.4,
                           ylim = c(-0.2, 0.2),
                           xlab = "rank", ylab = "slope",
                           sub = NULL,
                           ...)
}

%- maybe also 'usage' for other objects documented here.
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

  \item{output}{ Specifies the output. If \code{NULL}, the existing
    output device is used; if \code{output} is one of  \code{"X11",
      "windows", "quartz"}, 
    a new X11 (Unix)/windows (Windows)/quartz (Mac) device is created.
    If \code{output} is a string ending in one of \code{".pdf", ".png",
    ".eps", ".ps"} is given, a file of that name and type is created and
    closed afterwards.}

  \item{pch, cex, cex.lab, ylim, xlab, ylab}{
    As for \code{plot()}
  }

  \item{sub}{The subtitle. If \code{sub==NULL}, a string giving the tau and
    p-value of the Mann-Kendall test is used; otherwise as for \code{plot()}
  }

  \item{\dots}{Other arguments are passed on to \code{plot()}
  }
}


\details{
  The total dye bias appears to be the product of iGSDB and slide
  bias. In other words, it is monotonous (always increasing or always
  decreasing), both with respect to the intrinsic gene specific dye bias
  and with respect to the slide bias. This function orders genes by
  their iGSDB and the slides by slide bias. Subsequently a linear
  regression of each gene is done, with x being the slide bias rank (not
  the slide bias itself), and y being the M. The slopes of each linear
  regression line should be an increasing array of values, representing
  the 'fan' of lines. The array of slopes is plotted (versus the
  rank). Generally, a clear trend is visible for uncorrected
  hybridizations, and the trend has disappeared after dye bias correction.
  
}

\value{
  The order of the slide bias is returned, for use in plotting the
  behaviour of the regression slopes in the corrected data set.
}

\note{

  This function takes very long to compute, since it calculates
  regressions for each gene.

}


\author{Philip Lijnzaad \email{p.lijnzaad@umcutrecht.nl} }

\seealso{
  \code{\link{dyebias.monotonicity}},
  \code{\link{dyebias.trendplot}}
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

  layout(matrix(1:2, nrow=1,ncol=2))

  order <- dyebias.monotonicityplot(data=data.norm, 
                        iGSDBs=iGSDBs.estimated, # from e.g. dyebias.estimate.iGSDBs
                        order=NULL,              # i.e., order by increasing slide bias
                        output=NULL,
                        main="before correction"
                        )

  order <- dyebias.monotonicityplot(data=correction$data.corrected,
                        iGSDBs=iGSDBs.estimated,
                        order=order,             # order by the original slide bias
                        output=NULL,
                        main="after correction"
                       ) 

}



\references{
   Margaritis, T., Lijnzaad, P., van~Leenen, D., Bouwmeester, D.,
  Kemmeren, P., van~Hooff, S.R and Holstege, F.C.P. (2009).
  Adaptable gene-specific dye bias correction for two-channel DNA microarrays.
  \emph{Molecular Systems Biology, submitted}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}                   % silly, but one keyword is compulsary
\keyword{misc}                    % silly, but one keyword is compulsary
