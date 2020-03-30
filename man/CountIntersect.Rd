\name{CountIntersect}
\alias{CountIntersect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Find and Count Repeat Elements Overlapping with Genomic Intervals

}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Using a narrowPeak, broadPeak, or summits file obtained via MACS 2 and annotated via ChIPseeker, and a RepeatMasker file, aggregates repeat types, families and classes intersecting with the intervals in the peak file.
}
\usage{
CountIntersect(repeatMaskerFile, inputPeakFile, format, minoverlap = 0L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{repeatMaskerFile}{
%%     ~~Describe \code{repeatMaskerFile} here~~
data frame containing repeat information, generated via \link{read_rm}.
}
  \item{inputPeakFile}{
%%     ~~Describe \code{inputPeakFile} here~~
data frame containing peak intervals and annotation.
}
  \item{format}{
%%     ~~Describe \code{format} here~~
character vector identifying peak file format, one of "narrow", "broad", or "summit".
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
A list with three components.
  \item{RepeatName}{A tibble with two columns, respectively showing the repeat name and the count of instances intersecting with input peaks.}
  \item{RepeatFamily}{A tibble with two columns, respectively showing the repeat family and the count of instances intersecting with input peaks.}
  \item{RepeatClass}{A tibble with two columns, respectively showing the repeat class and the count of instances intersecting with input peaks.}
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
function (repeatMaskerFile, inputPeakFile, format, minoverlap = 0L)
{
    gr.input <- inputPeakFile
    rmsk <- FormattingRM(repeatMaskerFile)
    library(GenomicRanges)
    gr.rmsk <- with(rmsk, GenomicRanges::GRanges(chr, IRanges(start,
        end), strand = strand))
    values(gr.rmsk) <- DataFrame(RepeatName = rmsk$repeat_name,
        RepeatFamily = rmsk$repeat_family, RepeatType = rmsk$repeat_type)
    if (format == "narrow") {
        if (ncol(elementMetadata(gr.input)) != 0) {
            GenomicRanges::start(gr.input) <- GenomicRanges::start(gr.input) +
                gr.input$peak
            GenomicRanges::end(gr.input) <- GenomicRanges::start(gr.input) +
                1
        }
    }
    if (format == "braod") {
        m <- GenomicRanges::findOverlaps(gr.rmsk, gr.input, ignore.strand = TRUE,
            minoverlap = minoverlap)
    }
    else {
        m <- GenomicRanges::findOverlaps(gr.rmsk, gr.input, ignore.strand = TRUE)
    }
    gr.rmsk.matched <- gr.rmsk[queryHits(m)]
    mcols(gr.rmsk.matched) <- cbind.data.frame(mcols(gr.rmsk.matched),
        mcols(gr.input[subjectHits(m)]))
    df.rmsk.matched <- as.data.frame(gr.rmsk.matched)
    diff.rName <- setdiff(rmsk$repeat_name, df.rmsk.matched$RepeatName)
    diff.rFamily <- setdiff(rmsk$repeat_family, df.rmsk.matched$RepeatFamily)
    diff.rType <- setdiff(rmsk$repeat_type, df.rmsk.matched$RepeatType)
    library(dplyr)
    x <- df.rmsk.matched \%>\% count(RepeatName)
    colnames(x) <- c("RepeatName", "nRepeatName")
    x <- rbind(x, data.frame(RepeatName = diff.rName, nRepeatName = rep.int(0,
        length(diff.rName))))
    x1 <- df.rmsk.matched \%>\% count(RepeatFamily)
    colnames(x1) <- c("RepeatFamily", "nRepeatFamily")
    x1 <- rbind(x1, data.frame(RepeatFamily = diff.rFamily, nRepeatFamily = rep.int(0,
        length(diff.rFamily))))
    x2 <- df.rmsk.matched \%>\% count(RepeatType)
    colnames(x2) <- c("RepeatType", "nRepeatType")
    x2 <- rbind(x2, data.frame(RepeatType = diff.rType, nRepeatType = rep.int(0,
        length(diff.rType))))
    all.counts <- list(x, x1, x2)
    return(all.counts)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
