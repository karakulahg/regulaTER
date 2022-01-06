\name{MakePlots}
\alias{MakePlots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Create Exploratory Scatter and Bar Plots for Repeat Regions Enriched in Accessible Genomic Peaks
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Creates automated plots of significantly repeat regions enriched in accessible
genomic peaks identified by the TEAR function, using ggplot2 functions.
}
\usage{
MakePlots(enrichTEARResult, outDir, width = 4, height = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{enrichTEARResult}{
%%     ~~Describe \code{enrichPARsResult} here~~
a list of length 3, with the names "RepeatName", "RepeatType", and "RepeatFamily",
and data frames generated via \link{binom.test}, as a result of the \link{TEAR}
function.
}
  \item{outDir}{
%%     ~~Describe \code{outDir} here~~
character vector identifying target output directory of generated plots.
}
  \item{width}{
%%     ~~Describe \code{width} here~~
numeric identifying the width (x-axis length) of the resulting image files.
}
  \item{height}{
%%     ~~Describe \code{height} here~~
numeric identifying the height (y-axis length) of the resulting image files.
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
None. Generated output is placed in the directory identified by the outDir parameter.
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
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
