\name{FindMotifs}
\alias{FindMotifs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Find Enriched Genomic Motifs in Repeats Enriched in Promoters using HOMER
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Uses HOMER to identify enriched genomic motifs, such as transcription factor binding sites, in provided
list of DEG promoter enriched PARs, and places the results into user provided directory.
}
\usage{
FindMotifs(df, repeatMaskerFile, outDir, homerPath)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
%%     ~~Describe \code{df} here~~
data frame of binomial test results, multitest corrected using FDR adjustment, identifying enriched
repeat subfamilies, obtained using EnrichPARs or IdentiftDEFLinkedRepeats functions.
}
  \item{repeatMaskerFile}{
%%     ~~Describe \code{repeatMaskerFile} here~~
data frame containing repeat information, generated via \link{read_rm}.
}
  \item{outDir}{
%%     ~~Describe \code{outDir} here~~
character vector identifying target output directory of motif enrichment analysis results.
}
  \item{homerPath}{
%%     ~~Describe \code{homerPath} here~~
character vector identifying full path to HOMER binary directory.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
A local installation of HOMER genomic sequencing peak analysis software is required for
use of the FindMotifs function. HOMER software can be obtained via the instructions on the
[HOMER website](http://homer.ucsd.edu/homer/).

Heinz S, Benner C, Spann N, Bertolino E et al. Simple Combinations of Lineage-Determining Transcription Factors Prime cis-Regulatory Elements Required for Macrophage and B Cell Identities. Mol Cell 2010 May 28;38(4):576-589. PMID: 20513432

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
None. Generated output is placed in the directory identified by the outDir parameter. For details
of the HOMER software output, please see \href{http://homer.ucsd.edu/homer/}{HOMER website}.
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
function (df, repeatMaskerFile, outDir, homerPath)
{
    rmsk <- FormattingRM(repeatMaskerFile)
    binom.test.results <- df
    name <- rmsk[which(rmsk$repeat_name \%in\% as.vector(binom.test.results$RepeatName$RepeatName)),
        ]
    name$ID <- row.names(name)
    df.name <- name[, c(1, 2, 3, 8, 5, 4)]
    write.table(df.name, file = paste0(outDir, "/repeatName.bed"),
        quote = F, sep = "\t", row.names = F, col.names = F)
    library(marge)
    options(homer_path = homerPath)
    options(stringsAsFactors = F)
    if (nrow(df.name) > 0) {
        for (i in 1:10) {
            dir <- paste0(outDir, "/margeOutput/asRepeatName/",
                df.name[i, ]$repeat_name)
            print(dir)
            dir.create(dir, recursive = T)
            find_motifs_genome(df.name[i, ], dir, "hg38", overwrite = T)
        }
    }
  }
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
