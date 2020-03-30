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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (repeatMasterFile, peakFile, substRate)
{
    rmsk <- FormattingRM(repeatMasterFile, goal = "calcAge")
    gr.rmsk <- MakeGrangeObj(rmsk)
    gr.peak <- MakeGrangeObj(peakFile)
    overlapped <- GetOverlap(rmsk = rmsk, gr.input = gr.peak,
        format = "narrow", goal = "calcAge")
    df.PAR <- as.data.frame(overlapped)
    gr.nonPAR <- gr.rmsk[!gr.rmsk \%over\% overlapped, ]
    df.nonPAR <- as.data.frame(gr.nonPAR)
    x1 <- data.frame()
    if (nrow(df.PAR) != 0) {
        x1 <- df.PAR[, c("RepeatName", "perc_div")] \%>\% group_by(RepeatName) \%>\%
            summarise(Age = (mean(perc_div)))
    }
    x1$Type <- "PAR"
    x2 <- data.frame()
    if (nrow(df.nonPAR) != 0) {
        x2 <- df.nonPAR[, c("repeat_name", "perc_div")] \%>\% group_by(repeat_name) \%>\%
            summarise(Age = (mean(perc_div)))
    }
    x2$Type <- "nonPAR"
    colnames(x2)[1] <- "RepeatName"
    x <- rbind(x1, x2)
    x$Age <- (x$Age * 10)/substRate
    return(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
