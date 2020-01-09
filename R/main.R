
#### countIntersect description ####
#
#  Takes repeatMaster file and annotated narrowPeak file as input
#  Calculates number of repeats with overlapping peaks based on given category
#  (to add: category member filtering, input peak file type, multiple file type support, overlap expectation for broadPeak)
#
#### countIntersect function ####

MakeGrangeObj <- function(inputPeakFile){


  inputPeakFile <- as.data.frame(inputPeakFile)
  hit <- inputPeakFile$strand == "."
  inputPeakFile$strand[hit] <- "*"

  #### converts input peak data frame into GRanges object ####

  library(GenomicRanges)
  gr.input <-
    with(inputPeakFile,
         GenomicRanges::GRanges(seqnames, IRanges(start, end), strand = strand))
  if(ncol(inputPeakFile) > 6){
    values(gr.input) <- inputPeakFile[, 6:ncol(inputPeakFile)]
  }

  return(gr.input)
}


CountIntersect <-
  function(repeatMaskerFile,
           inputPeakFile) {

    gr.input <- inputPeakFile
    rmsk <- repeatMaskerFile

    #### rearranges repeatMasker file for function convenience ####

    library(biomartr)
    last <- as.data.frame(stringr::str_split_fixed(rmsk$matching_class, "/", 2))
    rmsk <-
      data.frame(
        "chr" = rmsk$qry_id,
        "start" = rmsk$qry_start,
        "end" = rmsk$qry_end,
        "strand" = rmsk$matching_repeat,
        "repeat_name" = rmsk$repeat_id,
        "repeat_type" = last$V1,
        "repeat_family" = last$V2
      )
    rmsk$strand <- as.character(rmsk$strand)
    rmsk$strand <- replace(rmsk$strand, rmsk$strand == "C", "-")

    #### converts repeatMasker data frame into GRanges object ####

    library(GenomicRanges)
    gr.rmsk <- with(rmsk, GenomicRanges::GRanges(chr, IRanges(start, end), strand = strand))
    values(gr.rmsk) <-
      DataFrame(
        RepeatName = rmsk$repeat_name,
        RepeatFamily = rmsk$repeat_family,
        RepeatType = rmsk$repeat_type
      )


    #### converts gr.input range into single nucleotide at summit location ####
    if(ncol(elementMetadata(gr.input))!= 0){
      GenomicRanges::start(gr.input) <- GenomicRanges::start(gr.input) + gr.input$peak
      GenomicRanges::end(gr.input) <- GenomicRanges::start(gr.input) + 1
    }

    #### finds repeat ranges with overlapping summits ####

    m <- GenomicRanges::findOverlaps(gr.rmsk, gr.input, ignore.strand = TRUE)
    gr.rmsk.matched <- gr.rmsk[queryHits(m)]

    #### adds the metadata from gr2 to GRanges of intersecting peaks ####

    mcols(gr.rmsk.matched) <- cbind.data.frame(mcols(gr.rmsk.matched),
                                               mcols(gr.input[subjectHits(m)]))

    #### converts matched GRanges object into data frame and counts instances of given column name ####

    df.rmsk.matched <- as.data.frame(gr.rmsk.matched)
    library(dplyr)
    x <- df.rmsk.matched %>% count(RepeatName)
    colnames(x) <- c("RepeatName", "nRepeatName")
    x.1 <- df.rmsk.matched %>% count(RepeatFamily)
    colnames(x.1) <- c("RepeatFamily","nRepeatFamily")
    x.2 <- df.rmsk.matched %>% count(RepeatType)
    colnames(x.2) <- c("RepeatType", "nRepeatType")

    all.counts <- list(x,x.1,x.2)

    ####

    return(all.counts)

  }



#### get shuffle genome interval with using bedr shuffle function ####

MakeShuffle<-function(inputPeakFile,genomeSizePath, numberOfShuffle=1,repeatMaskerFile ){

  gr.input <- inputPeakFile

  library(valr)
  genome <- read_genome(genomeSizePath)
  gr<-bed_shuffle(gr.input, genome)
  colnames(gr)<-c("seqnames","start","end","strand")
  gr <- MakeGrangeObj(gr)
  counts <- CountIntersect(repeatMaskerFile, gr)
  if(numberOfShuffle > 1){

    Rname <- as.data.frame(counts[[1]])
    Rfamily <- as.data.frame(counts[[2]])
    Rtype <- as.data.frame(counts[[3]])

    for(i in 1:numberOfShuffle){

      tmp <- bed_shuffle(gr.input, genome)
      colnames(tmp)<-c("seqnames","start","end","strand")
      tmp <- MakeGrangeObj(tmp)
      tmp.counts <- CountIntersect(repeatMaskerFile, tmp)

      Rname <- merge(Rname, as.data.frame(tmp.counts[[1]]),by = "RepeatName")
      Rfamily <- merge(Rfamily, as.data.frame(tmp.counts[[2]]),by = "RepeatFamily")
      Rtype <- merge(Rtype, as.data.frame(tmp.counts[[3]]),by = "RepeatType")

    }

    Rname$Mean <- round(rowMeans(Rname[,c(2:ncol(Rname))]))
    Rfamily$Mean <- round(rowMeans(Rfamily[,c(2:ncol(Rfamily))]))
    Rtype$Mean <- round(rowMeans(Rtype[,c(2:ncol(Rtype))]))

    counts <- list(Rname,Rfamily,Rtype)

  }


  return(counts)
}

  ####






