#### ToolX - Packages : ####
# library(biomartr)
# library(GenomicRanges)
# library(dplyr)
# library(valr)
# library(ToolX)

#### Tutorial Steps ####

# rmsk<-biomartr::read_rm("../test/hg19.fa.out") # repeat annotation
# input.file <- read.csv("../test/Annotated-Sorted_09_018U_004MSoton_MCF7-ZEB1_unind_Pol2_hs_i91_peaks.narrowPeak.tsv", header=TRUE, stringsAsFactors=FALSE, sep = "\t")

# gr.peak <- MakeGrangeObj(inputPeakFile = input.file)

##### peak.repeat.counts <- CountIntersect(repeatMaskerFile = rmsk, inputPeakFile = gr.peak)

# pathList <- list("Promoter" = "/home/nazmiye/Desktop/BIP/test/hg19.Promoter",
#      "Exon" = "/home/nazmiye/Desktop/BIP/test/hg19.Exons",
#      "Intron" = "/home/nazmiye/Desktop/BIP/test/hg19.Introns",
#      "5UTR" = "/home/nazmiye/Desktop/BIP/test/hg19.5Prime",
#      "3UTR" = "/home/nazmiye/Desktop/BIP/test/hg19.3Prime",
#      "Downstream" = "/home/nazmiye/Desktop/BIP/test/hg19.Downstream",
#      "genomeSizePath" = "/home/nazmiye/Desktop/BIP/test/hg19.chrom.sizes"
#   )


# test <- EnrichPARs(inputPeakFile = input.file, ShuffledPeak = gr.sh, pathList = pathList, numberOfShuffle = 2, repeatMaskerFile = rmsk)




