\name{dyebias.application.subset}
\alias{dyebias.application.subset}

\title{Return a subset of reporters that can be dye bias-corrected
  reasonably well}

\description{

  Convenience function returning a subset of reporters that can be
  expected to be corrected reasonably well. Often, the logical AND of
  this set and that of \code{maW(data.norm) == 1.0} is used.  The
  resulting subset is passed as the \code{application.subset}-argument
  to \code{\link{dyebias.apply.correction}}.
}

\arguments{

  \item{data.raw}{
    A \code{marrayRaw} object whose normalized data is to be dye bias-corrected.
  }

  \item{min.SNR}{
    The minimum signal to noise ratio to require. It is loosely defined
    here as the foreground over the background signal. The background
    signal may not be real; see below.
  }

  \item{use.background}{
    Logical indicating whether or not to use the background signals \cr 
    \code{maRb(data.raw)} and \code{maGb(data.raw)}. If the \code{data.raw}
    object does not have them, specify \code{use.background=FALSE}. This
    will use the smallest foreground of all reporters instead of the
    real backgrounds. 
  }
  
  \item{maxA}{
    The maximum signal that is still allowed.
  }
  
}

\value{

  A matrix of logicals with the same dimensions as those of
  \code{maRf{data.raw}} is returned.
}

\examples{
  ## First load data and estimate the iGSDBs
  ## (see dyebias.estimate.iGSDBs)

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
  }                                    % dontshow

  ### choose the estimators and which spots to correct:
  estimator.subset <- dyebias.umcu.proper.estimators(maInfo(maGnames(data.norm)))

  ### choose which genes to dye bias correct. Typically, this is based
  ### both on flagged spots and intensity
  application.subset <- maW(data.norm) == 1 &
          dyebias.application.subset(data.raw=data.raw, use.background=TRUE)

  summary(application.subset)

}

\details{

  This routine requires an \code{marrayRaw} object since only that
  contains the background intensities. If you only have normalized data,
  use something like

  \preformatted{
    bg <- matrix(0.5, nrow=maNspots(data.norm), ncol=maNsamples(data.norm))
    data.raw <- new("marrayRaw", maRf=maR(data.norm), maGf=maG(data.norm),
                                 maRb=bg, maGb=bg, maW=maW(data.norm))
  }
}

\seealso{\code{\link{dyebias.apply.correction}}}

\author{Philip Lijnzaad \code{p.lijnzaad@umcutrecht.nl}} 

\references{
  Margaritis, T., Lijnzaad, P., van Leenen, D., Bouwmeester, D.,
  Kemmeren, P., van Hooff, S.R and Holstege, F.C.P. (2009)
  Adaptable gene-specific dye bias correction for two-channel DNA microarrays.
  \emph{Molecular Systems Biology}, 5:266, 2009. doi: 10.1038/msb.2009.21.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}                          % silly, but one keyword is compulsary
