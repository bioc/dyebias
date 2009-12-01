\name{dyebias.rgplot}
\alias{dyebias.rgplot}
\alias{dyebias.maplot}

\title{Produce scatterplots of the hybridization, with strongest dye
  biases highlighted. }

\description{

  Plots the $log_2(R)$ vs. $log_2(G)$ (or alternatively $M$ vs. $A$)
  signal of one slide, highlighting the reporters with the strongest red
  and green dye bias. Two lines indicate two-fold change. See also
  Margaritis et al. (2009), Fig. 1

}

\synopsis{

  dyebias.rgplot(data, slide, iGSDBs, dyebias.percentile=5,
                 application.subset=TRUE, output=NULL,
                 xlim = c(log2(50),log2(50000)),
                 ylim = c(log2(50),log2(50000)),
                 xticks = c(100,1000,10000,10000),
                 yticks = c(100,1000,10000,10000),
                 pch = 19, cex = 0.3, cex.lab = 1.4,
                 \dots)

  dyebias.maplot(data, slide, iGSDBs, dyebias.percentile=5,
                 application.subset=TRUE, output=NULL,
                 xlim = c(6,16), ylim = c(-2,2),
                 pch = 19, cex = 0.3,  cex.lab = 1.4, \dots)
}

\arguments{

  \item{data}{The \code{marrayNorm} object to plot one slide of.  }

  \item{slide}{The index of the slide to plot; must be > 1, and < \code{maNsamples(data)}}

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

  \item{output}{ Specifies the output. If \code{NULL}, the existing
    output device is used; if \code{output} is one of  \code{"X11",
      "windows", "quartz"}, 
    a new X11 (Unix)/windows (Windows)/quartz (Mac) device is created.
    If \code{output} is a string ending in one of \code{".pdf", ".png",
    ".eps", ".ps"} is given, a file of that name and type is created and
    closed afterwards.}

  \item{xlim,ylim, xticks, yticks,pch,cex,cex.lab}{
    Graphical parameters; see \code{par()}
   }
  
  \item{\dots}{Other arguments (such as \code{main} etc.) are passed on to \code{plot()}.
  }
}

\value{
 None.
}


\note{The highlighted spots are \emph{all spots} with an iGSDB that lies
  in the top- or bottom- \code{dyebias.percentile} of iGSDBS. That is, not just
  the estimator genes are highlighted.}

\author{Philip Lijnzaad \email{p.lijnzaad@umcutrecht.nl} }

\seealso{
  \code{\link{dyebias.estimate.iGSDBs}},
  \code{\link{dyebias.apply.correction}},
  \code{\link{dyebias.rgplot}},
  \code{\link{dyebias.maplot}},
  \code{\link{dyebias.boxplot}},
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

  application.subset <- maW(data.norm) == 1 & dyebias.application.subset(data.raw=data.raw, use.background=TRUE)

  ### do the correction:
  correction <- dyebias.apply.correction(data.norm=data.norm,
                                         iGSDBs = iGSDBs.estimated,
                                         estimator.subset=estimator.subset,
                                         application.subset = application.subset,
                                         verbose=FALSE)

  }                                     % dontshow

  ## show both an RG-plot and an MA-plot of the uncorrected data and the
  ## corrected data next to each other. 

  slide <- 3                               # or any other other, of course

  layout(matrix(1:4, nrow=2,ncol=2, byrow=TRUE))

  dyebias.rgplot(data=data.norm,
                 slide=slide,
                 iGSDBs=iGSDBs.estimated,   # from dyebias.estimate.iGSDBs
                 main=sprintf("RG-plot, uncorrected, slide \%d", slide),
                 output=NULL)

  dyebias.rgplot(data=correction$data.corrected,
                 slide=slide,
                 iGSDBs=iGSDBs.estimated,
                 main=sprintf("RG-plot, corrected, slide \%d", slide),
                 output=NULL)


  dyebias.maplot(data=data.norm,
                 slide=slide,
                 iGSDBs=iGSDBs.estimated,
                 main=sprintf("MA-plot, uncorrected, slide \%d",slide),
                 output=NULL)

  dyebias.maplot(data=correction$data.corrected,
                 slide=slide,
                 iGSDBs=iGSDBs.estimated,
                 main=sprintf("MA-plot, corrected, slide \%d",slide),
                 output=NULL)
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
