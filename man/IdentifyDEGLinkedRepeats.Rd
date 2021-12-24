\name{IdentifyDEGLinkedRepeats}
\alias{IdentifyDEGLinkedRepeats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Identify Enrichment of Repeat Elements in Proximity of Differentially Expressed Genes
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculates the p-values of peak associated repeat (PAR) and differentially expressed gene (DEG) associations using negative binomial distribution and multi testing correction, compared to random sampling of non-peak associated repeats in proximity of DEGs.
}
\usage{
IdentifyDEGLinkedRepeats(enrichPARsResult, peaks, rmsk, genes, numberOfShuffle = 100, distance = 1e+05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{enrichPARsResult}{
%%     ~~Describe \code{enrichPARsResult} here~~
a list of length 3, returned by the EnrichPARs function.
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
  \item{alternative}
%%     ~~Describe \code{alternative} here~~
character vector identifying what the null hypothesis of the binomial test is. One of "greater", "less", or "two-sided", or the initial letter of such.
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (enrichPARsResult, peaks, rmsk, genes, numberOfShuffle = 100,
    distance = 1e+05)
{
    count <- enrichPARsResult
    input.file <- peaks
    repeatList <- count[[1]][which(count[[1]]["observe"] >= 2),
        c("RepeatName", "observe")]
    gr.input <- MakeGrangeObj(inputPeakFile = input.file)
    rmsk <- FormattingRM(rmsk)
    overlapped <- GetOverlap(rmsk = rmsk, gr.input = gr.input,
        format = "narrow")
    repeats.overlapped <- subset(overlapped, RepeatName \%in\%
        repeatList$RepeatName)
    d <- distanceToNearest(x = overlapped, subject = MakeGrangeObj(genes))
    m <- d[which(elementMetadata(d)$distance < distance), ]
    result.overlapped <- overlapped[queryHits(m)]
    observed.counts <- CountElements(result.overlapped, rmsk)[[1]]
    gr.rmsk <- MakeGrangeObj(rmsk)
    gr.nonPAR <- gr.rmsk[!gr.rmsk \%over\% overlapped, ]
    last.nonPAR <- data.frame()
    isfirsttime <- "true"
    if (nrow(repeatList) != 0) {
        for (i in 1:nrow(repeatList)) {
            part.nonPAR <- gr.nonPAR[sample(length(gr.nonPAR),
                repeatList$observe[i]), ]
            d <- distanceToNearest(x = part.nonPAR, subject = MakeGrangeObj(genes))
            m <- d[which(elementMetadata(d)$distance < distance),
                ]
            part.nonPAR <- part.nonPAR[queryHits(m)]
            if (length(part.nonPAR) != 0) {
                elementMetadata(part.nonPAR)$target <- as.character(repeatList$RepeatName[i])
                if (isfirsttime == "true") {
                  last.nonPAR <- part.nonPAR
                  isfirsttime <- "false"
                }
                else {
                  last.nonPAR <- append(last.nonPAR, part.nonPAR)
                }
            }
            else {
                isfirsttime == "true"
                last.nonPAR <- part.nonPAR
            }
        }
    }
    df.nonPAR <- as.data.frame(last.nonPAR)
    diff.rName <- setdiff(rmsk$repeat_name, df.nonPAR$target)
    library(dplyr)
    x <- data.frame()
    if (nrow(df.nonPAR) != 0) {
        x <- df.nonPAR \%>\% count(target)
        colnames(x) <- c("RepeatName", "nRepeatName")
    }
    expected.counts <- rbind(x, data.frame(RepeatName = diff.rName,
        nRepeatName = rep.int(0, length(diff.rName))))
    if (numberOfShuffle > 1) {
        for (i in 1:(numberOfShuffle)) {
            last.nonPAR <- data.frame()
            isfirsttime <- "true"
            for (j in 1:nrow(repeatList)) {
                tmp.nonPAR <- gr.nonPAR[sample(length(gr.nonPAR),
                  repeatList$observe[j]), ]
                d <- distanceToNearest(x = tmp.nonPAR, subject = MakeGrangeObj(genes))
                m <- d[which(elementMetadata(d)$distance < distance),
                  ]
                tmp.nonPAR <- tmp.nonPAR[queryHits(m)]
                if (length(tmp.nonPAR) != 0) {
                  elementMetadata(tmp.nonPAR)$target <- as.character(repeatList$RepeatName[j])
                  if (isfirsttime == "true") {
                    last.nonPAR <- tmp.nonPAR
                    isfirsttime <- "false"
                  }
                  else {
                    last.nonPAR <- append(last.nonPAR, tmp.nonPAR)
                  }
                }
                else {
                  isfirsttime == "true"
                  last.nonPAR <- tmp.nonPAR
                }
            }
            df.nonPAR <- as.data.frame(last.nonPAR)
            diff.rName <- setdiff(rmsk$repeat_name, df.nonPAR$target)
            x <- data.frame()
            if (nrow(df.nonPAR) != 0) {
                x <- df.nonPAR \%>\% count(target)
                colnames(x) <- c("RepeatName", "nRepeatName")
            }
            tmp <- rbind(x, data.frame(RepeatName = diff.rName,
                nRepeatName = rep.int(0, length(diff.rName))))
            expected.counts <- merge(expected.counts, tmp, by = "RepeatName")
        }
        expected.counts$Mean <- round(rowMeans(expected.counts[,
            c(2:ncol(expected.counts))]))
    }
    else {
        names(expected.counts)[2] <- "Mean"
    }
    rmsk.counts <- CountRM(rmsk)
    all.RepeatName <- merge(as.data.frame(rmsk.counts[[1]]),
        as.data.frame(observed.counts), by = "RepeatName")
    all.RepeatName <- merge(all.RepeatName, expected.counts[,
        c("RepeatName", "Mean")], by = "RepeatName")
    colnames(all.RepeatName) <- c("RepeatName", "rmsk", "observe",
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
    return(all.RepeatName)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
