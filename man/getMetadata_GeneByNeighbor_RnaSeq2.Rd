\name{getMetadata_GeneByNeighbor_RnaSeq2}
\alias{getMetadata_GeneByNeighbor_RnaSeq2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getMetadata_GeneByNeighbor_RnaSeq}
\description{
Get the metadata about the gene(s) within a given range.
}
\usage{
getMetadata_GeneByNeighbor_RnaSeq2(theStartPosition, theStopPosition, theChromosome, theStrand, theDataDir = "/rsrch1/bcb/batcheffects/GENE_REPORT/data", theVerboseFlag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theStartPosition}{ Inclusive position to start search. }
  \item{theStopPosition}{ Inclusive position to stop search. }
  \item{theChromosome}{ Chromosome (as string) to search, such as "X", "Y", "1", or "21" - NULL skips. }
  \item{theStrand}{Strand as a string "+" or "-" - NULL skips}
  \item{theDataDir}{ Default: For DQS Servers. If you are on a server where the GENE_REPORT directory is on a different path, you can pass in a value to override the default of "/rsrch1/bcb/batcheffects/GENE_REPORT/data". NOTE: this is a subdirectory of GENE_REPORT. }
  \item{theVerboseFlag}{ Default: FALSE. TRUE means write all output, which can be very verbose. }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The return value is a list of objects with gene data. The following functions work to pull data out of the object.
  \item{get.GeneSymbol}{ \link{get.GeneSymbol} }
  \item{get.GeneId}{ \link{get.GeneId} }
  \item{get.VersionIndex}{ \link{get.VersionIndex} }
  \item{get.LocationStart}{ \link{get.LocationStart} }
  \item{get.LocationEnd}{ \link{get.LocationEnd} }
  \item{get.Chromosome}{ \link{get.Chromosome} }
  \item{get.Strand}{ \link{get.Strand} }
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
#print(getMetadata_Gene_RnaSeq("?|100130426"))
#print(getMetadata_Gene_RnaSeq2("?|26823"))
#print(getMetadata_Gene_RnaSeq("A1BG"))
#print(getMetadata_Gene_RnaSeq2("tAKR"))
#print(getMetadata_Gene_RnaSeq("foo"))
#print(getMetadata_Gene_RnaSeq2("SNORD116-19"))
#print(getMetadata_Gene_RnaSeq("CT45A4"))
#print(getMetadata_GeneByNeighbor_RnaSeq(53700000, 53800000, "12", "-"))
#print(getMetadata_GeneByNeighbor_RnaSeq2(53700000, 53800000, "12", "+"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
