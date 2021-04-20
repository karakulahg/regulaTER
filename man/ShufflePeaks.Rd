\name{ShufflePeaks}
\alias{ShufflePeaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Shuffle Peak Ranges in Genomic Ranges of Shared Gene Context

}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Using a genomic peak file annotated via ChIPseeker  and a list of files identifying genomic regions corresponding to ChIPseeker annotation classes (promoter, 5' UTR, 3' UTR, exon, intron, downstream, and intergenic), shuffles the peak ranges in the genome to identify likelihood of peak - repeat relationships. For details, see \link{EnrichPARs}.
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
a list of length 7, with the names "Promoter", "Exon", "Intron", "5UTR", "3UTR", "Downstream", and "genomeSizePath", and character values declaring the file paths of BED files identifying genomic regions of respective annotations and chromosome sizes. Refer to details for further information.}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
integer used as seed value to obtain shuffle results. It is 0 by default and will result in unreproducible shuffled regions. To ensure results can be reproduced, replace by any positive integer value.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
ChIPseeker annotations are divided into seven major classes. Promoter regions are identified as up to 3 KB upstream of transcription start site, while downstream regions are identified as up to 3 KB downstream of transcription end site. Intergenic regions, while not explicitly defined for this function, are identifed as any genomic region that does not fit in with any other annotation class. These BED files can be obtained via \href{https://genome.ucsc.edu/index.html}{UCSC Genome Browser} using the Table Browser tool, or similar tools for other annotations of the user's choosing.

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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (peakFile, pathList, seed = 0)
{
    gr.input <- MakeGrangeObj(inputPeakFile = peakFile)
    gr.Promoter <- gr.input[grepl("Promoter", gr.input$annotation),
        ]
    gr.Exon <- gr.input[grepl("Exon", gr.input$annotation), ]
    gr.Intron <- gr.input[grepl("Intron", gr.input$annotation),
        ]
    gr.5UTR <- gr.input[grepl("5' UTR", gr.input$annotation),
        ]
    gr.3UTR <- gr.input[grepl("3' UTR", gr.input$annotation),
        ]
    gr.Intergenic <- gr.input[grepl("Intergenic", gr.input$annotation),
        ]
    gr.Downstream <- gr.Intron <- gr.input[grepl("Downstream",
        gr.input$annotation), ]
    write.table(as.data.frame(gr.Promoter), file = "Promoter.bed",
        quote = F, sep = "\t", row.names = F, col.names = F)
    system(paste(" bedtools shuffle -i Promoter.bed -g ", pathList$genomeSizePath,
        " -incl ", pathList$Promoter, " > shuffled.promoters "))
    sh.promoter <- read.csv("shuffled.promoters", header = F,
        sep = "\t")
    write.table(as.data.frame(gr.Exon), file = "Exon.bed", quote = F,
        sep = "\t", row.names = F, col.names = F)
    system(paste(" bedtools shuffle -i Exon.bed -g ", pathList$genomeSizePath,
        " -incl ", pathList$Exon, " > shuffled.exons "))
    sh.exon <- read.csv("shuffled.exons", header = F, sep = "\t")
    write.table(as.data.frame(gr.Intron), file = "Intron.bed",
        quote = F, sep = "\t", row.names = F, col.names = F)
    system(paste(" bedtools shuffle -i Intron.bed -g ", pathList$genomeSizePath,
        " -incl ", pathList$Intron, " > shuffled.introns "))
    sh.intron <- read.csv("shuffled.introns", header = F, sep = "\t")
    write.table(as.data.frame(gr.5UTR), file = "FiveUTR.bed",
        quote = F, sep = "\t", row.names = F, col.names = F)
    system(paste(" bedtools shuffle -i FiveUTR.bed -g ", pathList$genomeSizePath,
        " -incl ", pathList$`5UTR`, " > shuffled.5UTR "))
    sh.5UTR <- read.csv("shuffled.5UTR", header = F, sep = "\t")
    write.table(as.data.frame(gr.3UTR), file = "ThreeUTR.bed",
        quote = F, sep = "\t", row.names = F, col.names = F)
    system(paste(" bedtools shuffle -i ThreeUTR.bed -g ", pathList$genomeSizePath,
        " -incl ", pathList$`3UTR`, " > shuffled.3UTR "))
    sh.3UTR <- read.csv("shuffled.3UTR", header = F, sep = "\t")
    write.table(as.data.frame(gr.Downstream), file = "Downstream.bed",
        quote = F, sep = "\t", row.names = F, col.names = F)
    system(paste(" bedtools shuffle -i Downstream.bed -g ", pathList$genomeSizePath,
        " -incl ", pathList$Downstream, " > shuffled.downstream "))
    sh.downstream <- read.csv("shuffled.downstream", header = F,
        sep = "\t")
    write.table(as.data.frame(gr.Intergenic), file = "Intergenic.bed",
        quote = F, sep = "\t", row.names = F, col.names = F)
    system(paste(" bedtools shuffle -i Intergenic.bed -g ", pathList$genomeSizePath,
        " -excl ", pathList$Promoter, " -excl ", pathList$Exon,
        " -excl ", pathList$Intron, " -excl ", pathList$`5UTR`,
        " -excl ", pathList$`3UTR`, " -excl ", pathList$Downstream,
        " > shuffled.intergenic "))
    sh.intergenic <- read.csv("shuffled.intergenic", header = F,
        sep = "\t")
    all.shuffeledAnnots <- rbind(sh.promoter, sh.exon, sh.intron,
        sh.5UTR, sh.3UTR, sh.intergenic, sh.downstream)
    colnames(all.shuffeledAnnots) <- colnames(peakFile)
    gr <- MakeGrangeObj(all.shuffeledAnnots)
    return(gr)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
