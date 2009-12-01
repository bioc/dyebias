\name{dyebias.estimate.iGSDBs}

\alias{dyebias.estimate.iGSDBs}

\title{Estimate intrinsic gene specific dye biases (part of the GASSCO method)}

\description{
  Obtain estimates for the instrinsic gene-specific dye bias (iGSDB)
  using a set of normalized data, as part of the GASSCO method.
}

\synopsis{
  dyebias.estimate.iGSDBs(data.norm, is.balanced=TRUE, reference="ref", verbose=FALSE)
}

\arguments{

  \item{data.norm}{

    A \code{marrayNorm} object containing the data for estimating the
    dye bias. This object is supposed to be complete. In particular, \cr
    \code{maLabels(maGnames(data.norm))} must be set and must indicate the
    identities of the reporter sequence (i.e., oligo or cDNA
    sequence) of each spot. This helps identify replicate spots, which
    are averaged as part of the estimation. 

    If the data is unbalanced (so \code{is.balanced} is \code{FALSE}), \cr
    \code{maInfo(maTargets(data.norm))} is also required, and should
    contain at least two attributes: \code{Cy5} and \code{Cy3}. Both
    should indicate the factor value for the respective channel.}

  \item{is.balanced}{
    Logical indicating whether the data set represents a balanced design
    (which is in fact the most common case). A design is balanced if all
    factor values are present an equal number of times in both the
    forward and reverse dye orientations. A self-self design is by
    definition balanced (even if the number of slides is uneven). If
    \code{is.balanced} is \code{TRUE}, the iGSDB estimate is
    obtained by simply averaging, per reporter, all \eqn{M} values
    (and the value of the \code{reference} argument is ignored).

    If \code{is.balanced==FALSE}, the design is inferred from the
    \code{reference} argument, and subsequently the \code{limma}
    package is used to model the dye effect. This is typically done for
    an unbalanced data set, but there is no harm in setting
    \code{is.balanced=FALSE} for a design that by itself is already
    balanced. If there are no missing values in the data, the results of
    using the simple average and the limma procedure are identical
    (although LIMMA takes longer to compute the iGSDBs). If the data set
    contains many missing data points (NA's), the limma estimates differ
    slightly from the simple averaged estimates (although it is not
    clear which ones are better).
  }

  \item{reference}{ If the design is a single common reference,
    \code{reference} should be this common reference. If the design
    consists of a set of common reference designs, \code{reference}
    should be a vector listing all the common references, and the name
    of the factor value that is not the common reference should have its
    own common reference as a prefix. E.g., if two mutant strains
    \code{mutA} and \code{mutB} were assayed, each against a different
    reference \code{ref1} and \code{ref2}, the \code{reference}-argument
    would be \code{c("ref1", "ref2")}, and the \code{Cy3} and \code{Cy5}
    attributes of \code{maInfo(maTargets(data.norm))} must contain
    values from \code{"ref1:mutA", "ref2:mutA", "ref1:mutB",
      "ref2:mutB"}. (The colon is not important; the prefix is).
   }

  \item{verbose}{
    Logical, indicating wether or not to be verbose
  }
}

\value{

  A data frame is returned with as many rows as there are reporters
  (replicate spots have been averaged), and the following columns: \cr

  \item{reporterId}{The name of the reporter}
  \item{dyebias}{The intrinsic gene-specific dye bias (iGSDB) of this reporter}
  \item{A}{The average expression level of this reporter in the given data set}
  \item{p.value}{The \emph{p}-value for the \code{dyebias} ($H_0$:
          \code{dyebias} = 0). All \emph{p}-value are set to NA if
  they were not estimated (i.e., if limma was not run because
  \code{is.balanced} was TRUE)}

  This data frame is typically used as input to \code{\link{dyebias.apply.correction}}.
}

\examples{

  \dontshow{
     options(stringsAsFactors = FALSE)

     library(dyebias)
     library(dyebiasexamples)
     data(data.raw)
     data(data.norm)
  }                                     % dontshow

  iGSDBs.estimated <- dyebias.estimate.iGSDBs(data.norm,
                                             is.balanced=TRUE,
                                             verbose=FALSE)
  summary(iGSDBs.estimated)

 \dontrun{
    hist(iGSDBs.estimated$dyebias, breaks=50)
  }
}


\details{

  This function implements the first step of the GASSCO method:
  estimating the so-called intrinsic gene specific dye biases, or
  briefly iGSDB. They can be estimated from a (preferably large) data
  set containing either self-self experiments, or dye-swapped slides.

  The assumption underlying this approach is that with self-selfs, or
  with pairs of dye swaps, the only effect that can lead to systematic
  changes between Cy5 and Cy3, is in fact the dye effect.

  There are two cases to distinguish, the balanced case, and the
  unbalanced case. In the balanced case, the iGSDB estimate is simply
  the average \eqn{M} (where \cr \eqn{M = log_2(R/G) = log_2(Cy5/Cy3)}) over all
  slides.  A set of slides is balanced if all factor values are present
  in as many dye-swapped as non-dye-swapped slides.  A set of self-self
  slides is in fact a degenerate form of this, and is therefore also
  balanced.

  In the unbalanced case, one could omit slides until the data set is
  balanced. However, this is wasteful as we can use linear modelling to
  obtain estimates. We use the limma package for this (Smyth, 2005). The only
  unbalanced designs currently supported are a common reference design,
  and a set of common reference designs.

  There are is no weights or subset argument to this function; the 
  estimation is done for all reporters found. If there are replicate
  spots, they are averaged prior to the estimation (the reason being
  that we are not interested in p-values for the estimate)

  Having obtained the iGSDB estimates, the corrections can be applied
  to either to the hybridizations given by the \code{data.norm} argument,
  or to a different set of slides that is thought to have very similar
  iGSDBs. Applying the corrections is done with \cr \code{\link{dyebias.apply.correction}}. 
}

\seealso{dyebias.apply.correction}




\note{
  Note that the input data should be normalized, and that the dye swaps
  should \strong{not} have been swapped back.  After all, we're
  interested in the difference of Cy5 over Cy3, \strong{not} the
  difference of experiment over reference.
}

\author{Philip Lijnzaad \email{p.lijnzaad@umcutrecht.nl}} 

\references{

  Margaritis, T., Lijnzaad, P., van Leenen, D., Bouwmeester, D.,
  Kemmeren, P., van Hooff, S.R and Holstege, F.C.P. (2009)
  Adaptable gene-specific dye bias correction for two-channel DNA microarrays.
  \emph{Molecular Systems Biology}, 5:266, 2009. doi: 10.1038/msb.2009.21.

  Dudoit, S. and Yang, Y.H. (2002)  Bioconductor R packages for
  exploratory analysis and normalization of cDNA microarray data. In:
  Parmigiani, G., Garrett, E.S. , Irizarry, R.A., and Zeger, S.L. (eds.)
  \emph{The Analysis of Gene Expression Data: Methods and Software},
  Springer, New York.

%%   Smyth, G.K. (2004) Linear models and empirical Bayes methods for
%%   assessing differential expression in microarray
%%   experiments. \emph{Statistical Applications in Genetics and Molecular
%%     Biology}, \bold{3}, article 3.

   Smyth, G.K. (2005) Limma: linear models for microarray data.
   In: Gentleman, R., Carey, V., Dudoit, S., Irizarry, R. and
   Huber, W. (eds). \emph{Bioinformatics and Computational Biology
   Solutions using R and Bioconductor}, Springer, New York.

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}                          % silly, but one keyword is compulsary
