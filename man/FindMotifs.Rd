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
list of repeats enriched in accessible DEG promoters, and places the results into user provided directory.
}
\usage{
FindMotifs(df, repeatMaskerFile, peak, genes, distance, genome, outDir, homerPath, type, numberOfMotifs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
%%     ~~Describe \code{df} here~~
data frame of binomial test results, multitest corrected using FDR adjustment, identifying enriched
repeat subfamilies, obtained using EnrichPARs or IdentifyDEGLinkedRepeats functions.
}
  \item{repeatMaskerFile}{
%%     ~~Describe \code{repeatMaskerFile} here~~
data frame containing repeat information, generated via \link{read_rm}.
}
  \item{peak}{
%%     ~~Describe \code{peak} here~~
data frame containing peak intervals and annotation.
}
  \item{genes}{
%%     ~~Describe \code{genes} here~~
data frame with six columns - geneID: character vector identifying Ensembl IDs of differentially expressed genes,
or another suitable unique interval identifier; geneName: character vector identifying official gene symbols of
differentially expressed genes; seqnames: character vector identifying chromosome names of gene intervals;
start: integer vector identifying starting nucleotides of gene intervals; end: integer vector identifying
end nucleotides of gene intervals; strand: character vector identifying strand information of genes. Should
be the same object as used for IdentifyDEGLinkedRepeats if type is "linkedRepeats".
}
  \item{distance}{
%%     ~~Describe \code{distance} here~~
integer determining the distance in number of nucleotides, within which repeats are considered associated with DEGs.
Shouldbe the same distance as used for IdentifyDEGLinkedRepeats if type is "linkedRepeats".
}
  \item{genome}{
%%     ~~Describe \code{genome} here~~
character vector identifying the genome build used for HOMER operations. Usually an organism common name or a genome build name.
For details, see the HOMER configuration manual.
}
  \item{outDir}{
%%     ~~Describe \code{outDir} here~~
character vector identifying target output directory of motif enrichment analysis results.
}
  \item{homerPath}{
%%     ~~Describe \code{homerPath} here~~
character vector identifying full path to HOMER binary directory.
}
  \item{type}{
%%     ~~Describe \code{type} here~~
character vector identifying whether the df input argument is generated via EnrichPARs or IdentifyDEGLinkedRepeats. One of "enrichPeak"
or "linkedRepeats", for the two functions, respectively. This argument has no default and must be supplied for the function to work.
}
  \item{numberOfMotifs}{
%%     ~~Describe \code{type} here~~
boolean value identifying whether only the top ten repeats (TRUE) or all repeats identified in the list (FALSE) should be included in the
motif discovery analysis.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
A local installation of HOMER genomic sequencing peak analysis software is required for
use of the FindMotifs function. HOMER software can be obtained via the instructions on the
[HOMER website](http://homer.ucsd.edu/homer/).

Heinz S, Benner C, Spann N, Bertolino E et al. Simple Combinations of Lineage-Determining Transcription Factors Prime cis-Regulatory Elements Required for Macrophage and B Cell Identities. Mol Cell 2010 May 28;38(4):576-589. PMID: 20513432

Users should note that the motif discovery for a single repeat subfamily can take a
considerable amount of time, and if "numberOfMotifs = FALSE", it is recommended that users only
provide a subset of their input for the "df" argument to analyze selected repeats of interest.
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
