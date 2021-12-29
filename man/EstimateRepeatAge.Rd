\name{EstimateRepeatAge}
\alias{EstimateRepeatAge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Estimate Age of Repeat Subfamily Insertion into Genome
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Uses repeat annotation, genomic peak intervals, and genomic substitution rate to calculate evolutionary age of repeat subfamilies in accessible and inaccessible genomic regions.
}
\usage{
EstimateRepeatAge(repeatMasterFile, peakFile, substRate)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{repeatMasterFile}{
%%     ~~Describe \code{repeatMasterFile} here~~
data frame containing repeat information, generated via \link{read_rm}.
}
  \item{peakFile}{
%%     ~~Describe \code{peakFile} here~~
data frame containing peak intervals and annotation.
}
  \item{substRate}{
%%     ~~Describe \code{substRate} here~~
numeric indicating nucleotide substitution rate for genome of interest.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A data frame with columns corresponding to repeat subfamily, estimated age, and peak association status.
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
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
