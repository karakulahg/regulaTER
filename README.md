#### ToolX - Packages : ####
# library(biomartr)
# library(GenomicRanges)
# library(dplyr)
# library(valr)

#### Tutorial Steps ####

# rmsk<-biomartr::read_rm("../test/hg19.fa.out") # repeat annotation
# input.file <- read.csv("../test/Annotated-Sorted_09_018U_004MSoton_MCF7-ZEB1_unind_Pol2_hs_i91_peaks.narrowPeak.tsv", header=TRUE, stringsAsFactors=FALSE, sep = "\t")

# gr.peak <- MakeGrangeObj(inputPeakFile = input.file)

# peak.repeat.counts <- CountIntersect(repeatMaskerFile = rmsk, inputPeakFile = gr.peak)

# shuffled.peak.repeat.counts <- MakeShuffle(inputPeakFile = gr.peak, genomeSizePath = "../test/hg19.chrom.sizes", numberOfShuffle = 3, repeatMaskerFile = rmsk)


