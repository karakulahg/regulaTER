\name{EnrichPARs}
\alias{EnrichPARs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculate Enrichment of Repeats Associated with Genomic Peak Regions
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculates the p-values of peak and repeat associations using negative binomial distribution and multi testing correction, based on genome shuffling of peak regions.
}
\usage{
EnrichPARs(inputPeakFile, pathList, numberOfShuffle = 1, repeatMaskerFile, format, minoverlap = 0L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputPeakFile}{
%%     ~~Describe \code{inputPeakFile} here~~
GRanges object of peak information in one of narrowPeak, broadPeak, or summits file formats, annotated via ChIPseeker and read using the readPeakFile function, with a minimum of seqnames, ranges, strand required columns, as well as summit and annotation columns.
}
  \item{pathList}{
%%     ~~Describe \code{pathList} here~~
a list of length 7, with the names "Promoter", "Exon", "Intron", "5UTR", "3UTR", "Downstream", and "genomeSizePath", and character values declaring the file paths of BED files identifying genomic regions of respective annotations and chromosome sizes. Refer to \link{ShufflePeaks} for further information.
}
  \item{numberOfShuffle}{
%%     ~~Describe \code{numberOfShuffle} here~~
integer determining the number of shuffle operations to be performed for calculation of peak - repeat association enrichment.
}
  \item{repeatMaskerFile}{
%%     ~~Describe \code{repeatMaskerFile} here~~
data frame containing repeat information, generated via \link{read_rm}.
}
  \item{format}{
%%     ~~Describe \code{format} here~~
character vector identifying peak file format, one of "narrow", "broad", or "summit".
}
  \item{alternative}{
%%     ~~Describe \code{alternative} here~~
character vector identifying what the null hypothesis of the binomial test is. One of "greater", "less", or "two-sided", or the initial letter of such.
}
  \item{minobserved}{
%%     ~~Describe \code{minoverlap} here~~
integer identifying minimum number of observed peak - repeat associations to be included in the output.
}
  \item{minoverlap}{
%%     ~~Describe \code{minoverlap} here~~
integer identifying minimum number of overlapping bases between "broad" peak and repeat intervals.
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
A list with three components:
\item{RepeatName}{Results of the \link{binom.test} as applied to repeat subfamilies, using FDR adjustment of p-values.}
\item{RepeatFamily}{Results of the \link{binom.test} as applied to repeat families, using FDR adjustment of p-values.}
\item{RepeatType}{Results of the \link{binom.test} as applied to repeat clases, using FDR adjustment of p-values.}
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
