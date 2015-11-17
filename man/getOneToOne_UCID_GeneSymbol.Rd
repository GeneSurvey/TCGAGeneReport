\name{getOneToOne_UCID_GeneSymbol}
\alias{getOneToOne_UCID_GeneSymbol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getOneToOne_UCID_GeneSymbol}
\description{Given one or more gene symbols, give the list of corresponding UCSC ids}
\usage{
getOneToOne_UCID_GeneSymbol(theGeneSymbol, theIdDataDir = "/rsrch1/bcb/batcheffects/GENE_REPORT/iddata/downloads", theVerboseFlag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theGeneSymbol}{One or a vector of gene symbols from the getOneToOne_GeneSymbol_List. Note: these values may not correspond to values in provided data. }
  \item{theIdDataDir}{ Default: For DQS Servers. If you are on a server where the GENE_REPORT directory is on a different path, you can pass in a value to override the default of "/rsrch1/bcb/batcheffects/GENE_REPORT/iddata/downloads". NOTE: this is a subdirectory of GENE_REPORT. }
  \item{theVerboseFlag}{ Default: FALSE. TRUE means write all output, which can be very verbose. }
}
\details{
Given one or more gene symbols, give the list of corresponding UCSC ids.
Note: these values may not correspond to values in provided data.
}
\value{
A vector of values with the original request vector as names().
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (theGeneSymbol, theIdDataDir = "/rsrch1/bcb/batcheffects/GENE_REPORT/iddata/downloads", 
    theVerboseFlag = FALSE) 
{
    getOneToOne_Name_internal(theGeneSymbol, theIdDataDir, "getOneToOne_UCID_GeneSymbol", 
        theVerboseFlag = theVerboseFlag)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
