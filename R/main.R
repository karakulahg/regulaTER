
#### countIntersect description ####
#
#  Takes repeatMaster file and annotated narrowPeak file as input
#  Calculates number of repeats with overlapping peaks based on given category
#  (to add: category member filtering, input peak file type, multiple file type support, overlap expectation for broadPeak)
#
#### countIntersect function ####

CountIntersect <-
  function(repeatMaskerFile,
           inputPeakFile,
           colName,
           path = NULL) {
    rmsk <- repeatMaskerFile
    input.file <- inputPeakFile

    #### rearranges repeatMasker file for function convenience ####

    library(biomartr)
    last <-
      as.data.frame(stringr::str_split_fixed(rmsk$matching_class, "/", 2))
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
    gr.rmsk <-
      with(rmsk, GenomicRanges::GRanges(chr, IRanges(start, end), strand = strand))
    values(gr.rmsk) <-
      DataFrame(
        RepeatName = rmsk$repeat_name,
        RepeatFamily = rmsk$repeat_family,
        RepeatType = rmsk$repeat_type
      )
    gr.rmsk

    #### converts input peak data frame into GRanges object ####

    library(GenomicRanges)
    gr.input <-
      with(input.file,
           GenomicRanges::GRanges(seqnames, IRanges(start, end), strand = strand))
    values(gr.input) <- input.file[, 6:ncol(input.file)]
    gr.input

    #### converts gr.input range into single nucleotide at summit location ####

    GenomicRanges::start(gr.input) <-
      GenomicRanges::start(gr.input) + gr.input$peak
    GenomicRanges::end(gr.input) <- GenomicRanges::start(gr.input) + 1

    #### finds repeat ranges with overlapping summits ####

    m <-
      GenomicRanges::findOverlaps(gr.rmsk, gr.input, ignore.strand = TRUE)
    gr.rmsk.matched <- gr.rmsk[queryHits(m)]

    #### adds the metadata from gr2 to GRanges of intersecting peaks ####

    mcols(gr.rmsk.matched) <- cbind.data.frame(mcols(gr.rmsk.matched),
                                               mcols(gr.input[subjectHits(m)]))

    #### converts matched GRanges object into data frame and counts instances of given column name ####

    df.rmsk.matched <- as.data.frame(gr.rmsk.matched)
    library(dplyr)
    x <- df.rmsk.matched %>%
      count(df.rmsk.matched[, colName])

    #### writes output counts into .tsv file ####

    if (is.null(path)) {
      write.table(
        x,
        paste0(colName, "-count.tsv"),
        col.names = T,
        row.names = F,
        sep = "\t",
        quote = F
      )
    } else {
      write.table(
        x,
        path,
        col.names = T,
        row.names = F,
        sep = "\t",
        quote = F
      )
    }

  ####

  }
