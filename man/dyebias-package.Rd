\name{dyebias-package}
\alias{dyebias-package}
\alias{dyebias}

\docType{package}

\title{ 
  The GASSCO method for correcting for slide-dependent gene-specific dye bias 
}

\description{
 Many two-colour hybridizations suffer from a dye bias that is both gene-specific and slide-specific. The former depends on the content of the nucleotide used for labeling; the latter depends on the labeling percentage. The slide-dependency was hitherto not recognized, and made addressing the artefact impossible.  Given a reasonable number of dye-swapped pairs of hybridizations, or of same vs. same hybridizations, both the gene- and slide-biases can be estimated and corrected using the GASSCO method (Margaritis et al., Mol. Sys. Biol. 5:266 (2009), doi:10.1038/msb.2009.21)
}

\details{
  \tabular{ll}{
    Package: \tab dyebias\cr
    Type: \tab Package\cr
    Version: \tab  1.7.1\cr
    Date: \tab 26 May 2010\cr
    Licence: \tab GPL-3\cr
  }
}

\author{ 
  Philip Lijnzaad and Thanasis Margaritis\cr
  Philip Lijnzaad <p.lijnzaad@umcutrecht.nl> (Maintainer). 
}

\seealso{
  
  \code{\link{dyebias.application.subset}},
  \code{\link{dyebias.apply.correction}},
  \code{\link{dyebias.boxplot}},
  \code{\link{dyebias.estimate.iGSDBs}},
  \code{\link{dyebias.monotonicity}},
  \code{\link{dyebias.monotonicityplot}},
  \code{\link{dyebias.rgplot}},
  \code{\link{dyebias.trendplot}}.
}

\keyword{package}
