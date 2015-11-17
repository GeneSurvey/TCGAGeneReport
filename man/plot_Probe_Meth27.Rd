\name{plot_Probe_Meth27}
\alias{plot_Probe_Meth27}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot_Probe_Meth27
}
\description{
Plots stripcharts and replicate diagrams for the given arguments.
}
\usage{
plot_Probe_Meth27(theGeneEq, theOutputDir, theGeneReportDir = "/rsrch1/bcb/batcheffects/GENE_REPORT", theVerboseFlag = FALSE, theUseDeltaFlag = FALSE, theReplicateFlag=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theGeneEq}{ The probe to plot which corresponds to an id from the appropriate getNames_ function. }
  \item{theOutputDir}{ Directory to which to write files. Directory will be created if needed/possible. }
  \item{theGeneReportDir}{ Default: For DQS Servers. If you are on a server where the GENE_REPORT directory is on a different path, you can pass in a value to override the default of "/rsrch1/bcb/batcheffects/GENE_REPORT". }
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
function (theGeneEq, theOutputDir, theGeneReportDir = "/rsrch1/bcb/batcheffects/GENE_REPORT", 
    theVerboseFlag = FALSE, theUseDeltaFlag = FALSE) 
{
    myProbe <- getMetadata_Probe_Meth27(theGeneEq, file.path(theGeneReportDir, 
        "data"))
    myGeneData <- NULL
    myBarcodeDiseases <- NULL
    myBarcodeSampleType <- NULL
    myTag <- ""
    if (TRUE == theUseDeltaFlag) {
        myTag <- "Delta"
    }
    myGeneData <- getData_Probe_Meth27(theGeneEq, theGeneReportDir, 
        theVerboseFlag = theVerboseFlag, theUseDeltaFlag = theUseDeltaFlag)
    myBarcodeDiseases <- getMetadataPop_BarcodeDisease(file.path(theGeneReportDir, 
        "data"), theVerboseFlag = theVerboseFlag)
    myBarcodeSampleType <- getMetadataPop_BarcodeSamplecode(file.path(theGeneReportDir, 
        "data"), theVerboseFlag = theVerboseFlag)
    geneEqPre <- paste(padChromosomeName(myProbe@mChromosome), 
        "-", myProbe@mProbeLocation, "-", sep = "")
    plotGenericOutput(theGeneEq, theOutputDir, myGeneData, myBarcodeDiseases, 
        myBarcodeSampleType, "Meth27", "Beta Value", theVerboseFlag = theVerboseFlag, 
        theTag = myTag, theGeneEqPre = geneEqPre)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
