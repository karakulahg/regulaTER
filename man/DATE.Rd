\name{DATE}
\alias{DATE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
DEG Associated Transposable Elements - Identify Enrichment of Repeat Elements in Proximity of Differentially Expressed Genes
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculates the p-values of repeat regions associated with peaks and differentially expressed gene (DEG) associations using negative binomial distribution and multi testing correction, compared to random sampling of non-peak associated repeats in proximity of DEGs.
}
\usage{
DATE(enrichTEARResult, peaks, rmsk, genes, alternative = "greater", minobserved = 10, numberOfShuffle = 100, distance = 1e+05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{enrichTEARResult}{
%%     ~~Describe \code{enrichPARsResult} here~~
a list of length 3, returned by the TEAR function.
}
  \item{peaks}{
%%     ~~Describe \code{peaks} here~~
data frame containing peak intervals and annotation.
}
  \item{rmsk}{
%%     ~~Describe \code{rmsk} here~~
data frame containing repeat information, generated via \link{read_rm}.
}
  \item{genes}{
%%     ~~Describe \code{genes} here~~
data frame with six columns - geneID: character vector identifying Ensembl IDs of differentially expressed genes, or another suitable unique interval identifier; geneName: character vector identifying official gene symbols of differentially expressed genes; seqnames: character vector identifying chromosome names of gene intervals; start: integer vector identifying starting nucleotides of gene intervals; end: integer vector identifying end nucleotides of gene intervals; strand: character vector identifying strand information of genes.
}
  \item{alternative}{
%%     ~~Describe \code{alternative} here~~
character vector identifying what the null hypothesis of the binomial test is. One of "greater", "less", or "two-sided", or the initial letter of such.
}
  \item{minobserved}{
%%     ~~Describe \code{minoverlap} here~~
integer identifying minimum number of observed peak - repeat associations to be included in the output.
}
  \item{numberOfShuffle}{
%%     ~~Describe \code{numberOfShuffle} here~~
integer determining the number of shuffle operations to be performed for calculation of repeat - gene association enrichment.
}
  \item{distance}{
%%     ~~Describe \code{distance} here~~
integer determining the distance in number of nucleotides, within which repeats are considered associated with DEGs.
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
A data frame of binomial test results, multitest corrected using FDR adjustment, identifying repeat subfamilies enriched near DEGs.
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
