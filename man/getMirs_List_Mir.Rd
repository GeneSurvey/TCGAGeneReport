\name{getMirs_List_Mir}
\alias{getMirs_List_Mir}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getMirs_List_Mir}
\description{Get a list of valid mir ids (this is just for the getMir functions--not for the mirRNASeq data.)}
\usage{
getMirs_List_Mir(theZipFile="/rsrch1/bcb/batcheffects/GENE_REPORT/GeneSurvey.zip", theVerboseFlag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theZipFile}{ Default: For DQS Servers. If you are on a server where the GENE_REPORT directory is on a different path, you can pass in a value to override the default of "/rsrch1/bcb/batcheffects/GENE_REPORT/GeneSurvey.zip". NOTE: this is a subdirectory of GENE_REPORT. }
  \item{theVerboseFlag}{ Default: FALSE. TRUE means write all output, which can be very verbose. }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Return a vector of string ids
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
