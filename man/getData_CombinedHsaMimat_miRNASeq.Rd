\name{getData_CombinedHsaMimat_miRNASeq}
\alias{getData_CombinedHsaMimat_miRNASeq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
getData_CombinedHsaMimat_miRNASeq
}
\description{
Retrieves a list of values for all samples/diseases for the given combined HSA and MIMAT id.
}
\usage{
getData_CombinedHsaMimat_miRNASeq(theGeneEq, theZipFile="/rsrch1/bcb/batcheffects/GENE_REPORT/GeneSurvey.zip", theUseDeltaFlag = FALSE, theRemoveDupFlag = TRUE, theVerboseFlag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theGeneEq}{ A vector of combined HSA and MIMAT id names as found in the getNames_ function for this data type. This consists of ids similar to "hsa-mir-663.MIMAT0003326" and "hsa-mir-663b.MIMAT0005867". This means a single HSA format id can have multiple entries (MIMAT ids). }
  \item{theZipFile}{ Default: For DQS Servers. If you are on a server where the GENE_REPORT directory is on a different path, you can pass in a value to override the default of "/rsrch1/bcb/batcheffects/GENE_REPORT/GeneSurvey.zip". }
  \item{theUseDeltaFlag}{ Default: FALSE. A flag indicating to use the differences in value between tumor and normal samples. }
  \item{theRemoveDupFlag}{ Default: FALSE. Remove the duplicate barcodes from the same control samples used in different diseases. }
  \item{theVerboseFlag}{ Default: FALSE. TRUE means write all output, which can be very verbose. }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a matrix of values, with colnames being the sample ids (barcodes) and rownames the HSA/MIMAT ids. Samples and row ids are in sorted order, not requested order.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#getData_CombinedHsaMimat_miRNASeq(getNames_CombinedHsaMimat_miRNASeq()[1])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
