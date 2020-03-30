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
data frame containing peak intervals and annotation.
}
  \item{pathList}{
%%     ~~Describe \code{pathList} here~~
a list of length 7, with the names "Promoter", "Exon", "Intron", "5UTR", "3UTR", "Downstream", and "genomeSizePath", and character values declaring the file paths of BED files identifying genomic regions of respective annotations and chromosome sizes. Refer to \link{EnrichPARs} for further information.
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (inputPeakFile, pathList, numberOfShuffle = 1, repeatMaskerFile,
    format, minoverlap = 0L)
{
    gr.input <- MakeGrangeObj(inputPeakFile = inputPeakFile)
    observe.counts <- CountIntersect(repeatMaskerFile, gr.input,
        format, minoverlap)
    gr <- ShufflePeaks(inputPeakFile, pathList)
    expected.counts <- CountIntersect(repeatMaskerFile, gr, format,
        minoverlap)
    if (numberOfShuffle > 1) {
        Rname <- as.data.frame(expected.counts[[1]])
        Rfamily <- as.data.frame(expected.counts[[2]])
        Rtype <- as.data.frame(expected.counts[[3]])
        for (i in 1:numberOfShuffle) {
            tmp <- ShufflePeaks(inputPeakFile, pathList)
            tmp.counts <- CountIntersect(repeatMaskerFile, tmp,
                format, minoverlap)
            tmp.Rname <- as.data.frame(tmp.counts[[1]])
            colnames(tmp.Rname) <- c("RepeatName", i)
            tmp.Rfamily <- as.data.frame(tmp.counts[[2]])
            colnames(tmp.Rfamily) <- c("RepeatFamily", i)
            tmp.Rtype <- as.data.frame(tmp.counts[[3]])
            colnames(tmp.Rtype) <- c("RepeatType", i)
            Rname <- merge(Rname, tmp.Rname, by = "RepeatName")
            Rfamily <- merge(Rfamily, tmp.Rfamily, by = "RepeatFamily")
            Rtype <- merge(Rtype, tmp.Rtype, by = "RepeatType")
        }
        Rname$Mean <- round(rowMeans(Rname[, c(2:ncol(Rname))]))
        Rfamily$Mean <- round(rowMeans(Rfamily[, c(2:ncol(Rfamily))]))
        Rtype$Mean <- round(rowMeans(Rtype[, c(2:ncol(Rtype))]))
        expected.counts <- list(Rname, Rfamily, Rtype)
    }
    rmsk <- FormattingRM(rmsk)
    rmsk.counts <- CountRM(rmsk)
    all.RepeatName <- merge(as.data.frame(rmsk.counts[[1]]),
        as.data.frame(observe.counts[[1]]), by = "RepeatName")
    all.RepeatName <- merge(all.RepeatName, as.data.frame(expected.counts[[1]][,
        c(1, ncol(expected.counts[[1]]))]), by = "RepeatName")
    colnames(all.RepeatName) <- c("RepeatName", "rmsk", "observe",
        "expected")
    all.RepeatFamily <- merge(as.data.frame(rmsk.counts[[2]]),
        as.data.frame(observe.counts[[2]]), by = "RepeatFamily")
    all.RepeatFamily <- merge(all.RepeatFamily, as.data.frame(expected.counts[[2]][,
        c(1, ncol(expected.counts[[1]]))]), by = "RepeatFamily")
    colnames(all.RepeatFamily) <- c("RepeatFamily", "rmsk", "observe",
        "expected")
    all.RepeatType <- merge(as.data.frame(rmsk.counts[[3]]),
        as.data.frame(observe.counts[[3]]), by = "RepeatType")
    all.RepeatType <- merge(all.RepeatType, as.data.frame(expected.counts[[3]][,
        c(1, ncol(expected.counts[[1]]))]), by = "RepeatType")
    colnames(all.RepeatType) <- c("RepeatType", "rmsk", "observe",
        "expected")
    test <- function(x, p, n) {
        binom.test(x, p, n)
    }
    b.rName <- mapply(test, all.RepeatName$observe, all.RepeatName$rmsk,
        (all.RepeatName$expected/all.RepeatName$rmsk))
    all.RepeatName$p.value <- do.call(rbind, b.rName["p.value",
        ])
    all.RepeatName$p.adjust.value <- p.adjust(all.RepeatName$p.value,
        method = "fdr", n = length(all.RepeatName$p.value))
    b.rFamily <- mapply(test, all.RepeatFamily$observe, all.RepeatFamily$rmsk,
        (all.RepeatFamily$expected/all.RepeatFamily$rmsk))
    all.RepeatFamily$p.value <- do.call(rbind, b.rFamily["p.value",
        ])
    all.RepeatFamily$p.adjust.value <- p.adjust(all.RepeatFamily$p.value,
        method = "fdr", n = length(all.RepeatFamily$p.value))
    b.rType <- mapply(test, all.RepeatType$observe, all.RepeatType$rmsk,
        (all.RepeatType$expected/all.RepeatType$rmsk))
    all.RepeatType$p.value <- do.call(rbind, b.rType["p.value",
        ])
    all.RepeatType$p.adjust.value <- p.adjust(all.RepeatType$p.value,
        method = "fdr", n = length(all.RepeatType$p.value))
    binom.test.results <- list(RepeatName = subset(all.RepeatName,
        p.value < 0.001 & observe > expected), RepeatFamily = subset(all.RepeatFamily,
        p.value < 0.001 & observe > expected), RepeatType = subset(all.RepeatType,
        p.value < 0.001 & observe > expected))
    return(binom.test.results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
