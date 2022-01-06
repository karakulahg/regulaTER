\name{ShufflePeaks}
\alias{ShufflePeaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Shuffle Peak Ranges in Genomic Ranges of Shared Gene Context

}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Using a genomic peak file annotated via ChIPseeker  and a list of files identifying complement genomic regions corresponding to ChIPseeker annotation classes (promoter, 5' UTR, 3' UTR, exon, intron, downstream, and intergenic), shuffles the peak ranges in the genome to identify likelihood of peak - repeat relationships. For details, see \link{EnrichPARs}.
}
\usage{
ShufflePeaks(peakFile, pathList, seed = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peakFile}{
%%     ~~Describe \code{peakFile} here~~
GRanges object of peak information in one of narrowPeak, broadPeak, or summits file formats, annotated via ChIPseeker and read using the readPeakFile function, with a minimum of seqnames, ranges, strand required columns, as well as summit and annotation columns.
}
  \item{pathList}{
%%     ~~Describe \code{pathList} here~~
a list of length 7, with the names "Promoter", "Exon", "Intron", "5UTR", "3UTR", "Downstream", and "genomeSizePath", and character values declaring the file paths of BED files identifying complement genomic regions of respective annotations and chromosome sizes. Refer to details for further information.
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
integer used as seed value to obtain shuffle results. It is 0 by default and will result in unreproducible shuffled regions. To ensure results can be reproduced, replace by any positive integer value.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
ChIPseeker annotations are divided into seven major classes. Promoter regions are identified as up to 3 KB upstream of transcription start site, while downstream regions are identified as up to 3 KB downstream of transcription end site. Intergenic regions, while not explicitly defined for this function, are identifed as any genomic region that does not fit in with any other annotation class. These BED files can be obtained via \href{https://genome.ucsc.edu/index.html}{UCSC Genome Browser} using the Table Browser tool, or similar tools for other annotations of the user's choosing.

In the current version of regulaTER, the complement intervals of the regions must be provided for each object in the list. They use the same format as the obtained BED files, but cover all genomic regions not found in the category. The user can generate these files using another genome arithmetic tool, such as bedtools complement (-i <BED> -g <genome>), available on the [bedtools Suite](https://bedtools.readthedocs.io/en/latest/index.html).

genomeSizePath points to a file with two columns, corresponding to chromosome name and size, delimited by tabs.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A GRanges object.
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
