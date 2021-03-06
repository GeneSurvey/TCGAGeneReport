\name{plot_GeneSymbol_Meth27}
\alias{plot_GeneSymbol_Meth27}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot_GeneSymbol_Meth27
}
\description{
Plots stripcharts and replicate diagrams for the given arguments.
Also plots heatmaps for probes.
Uses the gene symbol in the argument to look up probes and plots for those probes.
}
\usage{
plot_GeneSymbol_Meth27(theGene, theOutputDir, theZipFile="/rsrch1/bcb/batcheffects/GENE_REPORT/GeneSurvey.zip", theVerboseFlag = FALSE, theUseDeltaFlag = FALSE, theReplicateFlag=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theGene}{ The gene symbol to plot which corresponds to an id from the appropriate getNames_ function and maps to one or more probes. }
  \item{theOutputDir}{ Directory to which to write files. Directory will be created if needed/possible. }
  \item{theZipFile}{ Default: For DQS Servers. If you are on a server where the GENE_REPORT directory is on a different path, you can pass in a value to override the default of "/rsrch1/bcb/batcheffects/GENE_REPORT/GeneSurvey.zip". }
  \item{theVerboseFlag}{ Default: FALSE. TRUE means write all output, which can be very verbose. }
  \item{theUseDeltaFlag}{ Default: FALSE. A flag indicating to use the differences in value between tumor and normal samples. }
  \item{theReplicateFlag}{ Default: FALSE. A flag indicating to write the replicate plots. }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of files output.
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
function (theGene, theOutputDir, theZipFile="/rsrch1/bcb/batcheffects/GENE_REPORT/GeneSurvey.zip", 
    theVerboseFlag = FALSE, theUseDeltaFlag = FALSE) 
{
    plot_fromGene_internal(theGene, theOutputDir, theZipFile, 
        theMapGeneToProbesFunction = getNames_ProbeFromGeneSymbol_Meth27, 
        theReadProbeFunction = getMetadata_Probe_Meth27, theReadGeneFunction = getData_Probe_Meth27, 
        theDataType = "Meth27", theDataTypeLabel = "Beta Value", 
        theVerboseFlag = theVerboseFlag, theUseDeltaFlag = theUseDeltaFlag)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
