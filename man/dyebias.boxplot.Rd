\name{dyebias.boxplot}
\alias{dyebias.boxplot}

\title{Creates boxplots of the reporters with the strongest dye bias}

\description{

  The aim of this routine is to show the magnitude of the dye bias
  across the data set, as well as the extent to which the GASSCO method
  could get rid of it.  Typically, two boxplots would be shown, one
  before, one after dye bias correction. For esthetic reasons, the
  boxplots are usually ordered by the overal slide bias of the
  uncorrected data set. See also Margaritis et~al. (2009), Fig.~1 and 3.

}

\synopsis{

dyebias.boxplot <- function(data, iGSDBs, dyebias.percentile=5, application.subset=TRUE,
                            order, output=NULL,
                            ylim=c(-4,4),
                            \dots)
}

\arguments{

  \item{data}{The \code{marrayNorm} object to boxplot. 
  }

  \item{iGSDBs}{A data frame with intrinsic gene-specific dye biases,
    the same as that used in \code{\link{dyebias.apply.correction}},
    probably returned by \cr
    \code{\link{dyebias.estimate.iGSDBs}}; see
    there for documentation.
  }

  \item{dyebias.percentile}{
    The percentile of intrinsic gene specific dye biases (iGSDBs) for
    which to highlight the reporters. 
  }

  \item{application.subset}{
    The set of reporters that was eligible for dye bias correction; same
    argument as for \code{\link{dyebias.apply.correction}}.
  }

  \item{order}{

    If \code{order==FALSE}, no ordering of slides prior to boxplotting
    takes place.  If \code{order==NULL}, the slides are sorted by
    increasing slide bias prior to boxplotting. This is typically done
    for data that is not yet dye bias corrected.  This order is also
    returned as a value. If an \code{order!=NULL}, the slides are put
    this order before boxplotting.  This is typically done for a dye
    bias-corrected data set, using the order of the uncorrected set.
    (See also Fig. 3 in the paper). }

  \item{output}{ Specifies the output. If \code{NULL}, the existing
    output device is used; if \code{output} is one of  \code{"X11",
      "windows", "quartz"}, 
    a new X11 (Unix)/windows (Windows)/quartz (Mac) device is created.
    If \code{output} is a string ending in one of \code{".pdf", ".png",
    ".eps", ".ps"} is given, a file of that name and type is created and
    closed afterwards.}

  \item{ylim}{As for \code{boxplot()}}
  
  \item{\dots}{
    Other arguments (such as \code{main}, etc.)
    are passed on to \code{boxplot()}.
  }
}

\value{
  The order obtained, for use in a later call to this same function.
}

\seealso{
  \code{\link{dyebias.estimate.iGSDBs}},
  \code{\link{dyebias.apply.correction}},
  \code{\link{dyebias.rgplot}},
  \code{\link{dyebias.maplot}},
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

  ylim <- c(-1, 1)

  layout(matrix(1:2, nrow=1,ncol=2))

  order <- dyebias.boxplot(data=data.norm, 
                        iGSDBs=iGSDBs.estimated,  # from e.g. dyebias.estimate.iGSDBs
                        order=NULL,               # i.e., order by increasing slide bias
                        output=NULL,
                        main="before correction",
                        ylim=ylim)

  order <- dyebias.boxplot(data=correction$data.corrected, # from dyebias.apply.correction
                        iGSDBs=iGSDBs.estimated,
                        order=order,              # order by the original slide bias
                        output=NULL,
                        main="after correction",
                        ylim=ylim
                        )

}

\author{Philip Lijnzaad \email{p.lijnzaad@umcutrecht.nl} }


\references{
  Margaritis, T., Lijnzaad, P., van~Leenen, D., Bouwmeester, D.,
  Kemmeren, P., van~Hooff, S.R and Holstege, F.C.P. (2009)
  Adaptable gene-specific dye bias correction for two-channel DNA microarrays.
  \emph{Molecular Systems Biology, submitted}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}                   % silly, but one keyword is compulsary
\keyword{misc}                    % silly, but one keyword is compulsary
