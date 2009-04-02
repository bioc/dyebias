\name{dyebias.apply.correction}
\alias{dyebias.apply.correction}

\title{Perform dye bias correction using the GASSCO method}

\description{Corrects the gene- and slide specific dye bias in a data set, using
  the GASSCO method by Margaritis et~al.
}

\synopsis{dyebias.apply.correction(data.norm,
                         iGSDBs,
                         estimator.subset=TRUE, 
                         application.subset=TRUE, 
                         dyebias.percentile=5,
                         minmaxA.perc=25,
                         minA.abs=NULL,                         
                         maxA.abs=NULL,
                         verbose=FALSE)
}

\arguments{

  \item{data.norm}{
    A \code{marrayNorm} object containing the data whose dye bias should
    be corrected. This object must be a complete \code{marrayNorm} object. In particular,
    \code{maLabels(maGnames(data.norm))} should be set and indicate the
    identities of the spots. Spots with the same ID should contain the same
    oligo or cDNA sequence, and will receive the same dye bias correction.
  }

  \item{iGSDBs}{ 
    A data frame with the intrinsic gene specific dye bias per reporter
    (i.e., oligo or cDNA).  The data frame would typically have come
    from a call to \code{\link{dyebias.estimate.iGSDBs}}, but this is 
    not necessary; other estimates can also be used.

    The data frame must have (at least) the following columns:
    \itemize{
      \item{reporterId}{The name of the reporter. This must match the IDs
        in %% cannot use \cr for some reason here; hence the new paragraph.
        
        \code{maLabels(maGnames(data.norm))}
      }
      
      \item{dyebias}{ An estimate of the dye bias }
      
      \item{A}{
        
        The average expression value \eqn{A} of this reporter. (\eqn{A =
          (log_2(R)+log_2(G))/2 = (log_2(Cy5)+log_2(Cy3))/2 }). The
        \eqn{A}-value is used to base exclusions on. If you don't have it,
        you can use any value (but realize that the \code{minmaxA.perc,
          minA.abs, maxA.abs} arguments are still applied).
      }
    }                                   % itemize

    The order of the rows in this data frame is irrelevant. There must
    be no rows with duplicate \code{reporterId} in this frame. 
    
    For any reporter in \code{data.norm} that is not in the
    \code{iGSDBs} data frame, an iSGDB of 0.00 is used, i.e. data from
    such reporters is not dye bias-corrected.

  }
      
  \item{estimator.subset}{

    An index indicating which reporters are fit
    to be used as estimators of the slide bias. This set of reporters is
    used throughout the whole data set. Reporters that are typically
    excluded are those corresponding to parasitic DNA elements or
    mitochondrial genes.
  }

  \item{application.subset}{
    An index indicating which values must be dye
    bias-corrected. It should be either a vector with as many values as
    spots, or a matrix  of the same dimensions as
    \code{maM(data.norm)}. In former case, the selected spots on all
    slides with be dye bias-corrected; in the latter, selected spots on
    selected slides will corrected.

    Often it is prudent not to dye bias-correct measurements that are
    close to the detection limit or close to signal saturation.  A
    convenience function for this is provided; see \cr
    \code{\link{dyebias.application.subset}}.
  }

  \item{dyebias.percentile}{
    The slide bias estimation uses a small subset of reporters having
    the strongest green or red iGSDB, as specified by this
    percentile. The default should suffice in practically all cases.
  }

  \item{minmaxA.perc}{

    To obtain a robust estimate of the slide bias, the range of the
    average expression \eqn{A} is trimmed by \code{minmaxA.perc} percent
    on both sides; only reporters lying inside this trimmed range are
    considered as estimators of the slide bias. The default value is 25,
    meaning that top \code{dyebias.percentile} red- and green-biased
    spots within the the middle two average expression quartiles are
    used. This should suffice in practically all cases.

  }

  \item{minA.abs}{

    If specified, reporters with an average expression
    (\eqn{A}) lower than this value are never considered as estimators
    of the slide bias. If not specified, reporters with an
    \eqn{A}-percentile < \code{minmaxA.perc} are not considered.

  }

  \item{maxA.abs}{ If specified, reporters with an average expression
  (\eqn{A}) greater than this are never considered as estimators of the
  slide bias. If not specified, reporters with an \eqn{A}-percentile <
  \code{100-minmaxA.perc} are not considered.}


  \item{verbose}{ Logical speficying whether to be verbose or not }
}


\value{

  The data returned is a list wit the following elements
  
  \item{data.corrected}{A \code{marrayNorm} object of the same 'shape' as
    the input \code{data.norm}, but with corrected \eqn{M} values.
  }

  \item{estimators}{Another list, containing the details of the
    reporters that were used to obtain an estimate of the slide bias. 
    The contents of the \code{estimators} list are:
    \itemize{
      \item{green.ids}{The IDs of the reporters having the strongest green effect.}
      \item{green.cutoff}{All reporters in green.ids have an iGSDB below this value.}
      \item{green.subset}{An index into the reporters having the strongest green effect.}
      \item{red.ids}{The IDs of the reporters having the strongest red effect.}
      \item{red.cutoff}{All reporters in green.ids have an iGSDB above this value.}
      \item{red.subset}{An index into the reporters having the strongest
        red effect.}
    }
  }

  \item{summary}{
    A data frame summarizing the correction process per slide. It
    consist of the following columns:
    \itemize{
      \item{slide}{The slide number}
      \item{file}{Which file it came from}
      \item{green.bias}{The green bias of this slide}
      \item{red.bias}{The red bias of this slide}
      \item{green.correction}{The correction based on only the
        green bias of this slide}
      \item{red.correction}{The correction based on only the
        red bias of this slide}
      \item{avg.correction}{The total correction factor of this slide. This
        is in fact the slide bias}
      \item{var.ratio}{The ratio of the variance of \eqn{M} after and before the
        correction. The smaller this number, the smaller the variance of \eqn{M}
        around the mean has become, providing a measure of the success of the
        dye bias correction. Only data points that were in the
        \code{application.subset} are considered.}
      \item{reduction.perc}{As \code{var.ratio}, but expressed as a percentage.
        The larger this value, the greater the correction.}
      \item{p.value}{The p-value for the signficance of the reduction in
        variance (\eqn{F}-test; \eqn{H_0}: variances before and after correction are
        identical)}
    }
  }
}

\details{

  This function corrects the gene-specific dye bias of two-colour
  microarrays using the GASSCO method. This method is general, robust
  and fast, and is based on the observation that the total bias per gene
  is the product of a slide-specific factor (strongly related to the
  labeling percentage) and an intrinsic gene-specific factor (iGSDB),
  which is strongly related to the probe sequence.

  The slide bias is estimated from the total bias of the
  \code{dyebias.percentile} percentage of reporters having the strongest
  iGSDB. The iGSDBs can be estimated with \cr
  \code{\link{dyebias.estimate.iGSDBs}}.

  If the signal of certain oligos is too weak, or in contrast, tends to
  be saturated, they are no good estimator of the slide bias.
  Therefore, only reporters with an average expression level \eqn{A}
  that is not too extreme are allowed to be slide bias estimators. (This
  is the reason for the \code{A}-column in the \code{iGSDBs} data
  frame).

  Full control over which reporters to allow as slide bias estimators is
  given by the arguments \code{minmaxA.perc}, \code{minA.abs}, and
  \code{maxA.abs}; see there for details. To not exclude any reporter
  (e.g., when \eqn{A} is not available and therefore artificially set),
  you can use \code{minA.abs= -Inf} and \code{maxA.abs = Inf}.

  For further details concerning the method, see the \code{dyebias}
  vignette and the publication. If your research benefits from using this
  package, we kindly request that you cite this work.
}

\author{Philip Lijnzaad \email{p.lijnzaad@umcutrecht.nl} }

\note{
  Note that the input data should be normalized, and that the dye swaps
  should \strong{not} have been swapped back (if needed, this can of
  course be done afterwards).
}

\seealso{
  \code{\link{dyebias.estimate.iGSDBs}},
  \code{\link{dyebias.application.subset}},
  \code{\link{dyebias.rgplot}},
  \code{\link{dyebias.maplot}},
  \code{\link{dyebias.boxplot}},
  \code{\link{dyebias.trendplot}}
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

  ### choose which genes to dye bias correct:
  application.subset <- (maW(data.norm) == 1 &
               dyebias.application.subset(data.raw=data.raw, use.background=TRUE))

  ### do the correction:
  correction <- dyebias.apply.correction(data.norm=data.norm,
                                         iGSDBs = iGSDBs.estimated,
                                         estimator.subset=estimator.subset,
                                         application.subset = application.subset,
                                         verbose=FALSE)
  
  \dontrun{
     edit(correction$summary)
  }

  ## give overview:
  correction$summary[,c("slide", "file", "reduction.perc", "p.value")]

  ## and summary:
  summary(as.numeric(correction$summary[, "reduction.perc"]))
}

\references{
  Margaritis, T., Lijnzaad, P., van~Leenen, D., Bouwmeester, D.,
  Kemmeren, P., van~Hooff, S.R and Holstege, F.C.P. (2009).
  Adaptable gene-specific dye bias correction for two-channel DNA microarrays.
  \emph{Molecular Systems Biology, submitted}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}                          % silly, but one keyword is compulsary
