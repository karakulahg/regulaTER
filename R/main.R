

source("R/plotScripts.R")

MakeGrangeObj <- function(inputPeakFile) {
  options(warn = -1)
  inputPeakFile <- as.data.frame(inputPeakFile)
  hit <- inputPeakFile$strand == "."
  inputPeakFile$strand[hit] <- "*"

  #### converts input peak data frame into GRanges object ####

  library(GenomicRanges)
  gr.input <-
    with(inputPeakFile,
         GenomicRanges::GRanges(seqnames, IRanges(start, end), strand = strand))
  if (ncol(inputPeakFile) > 6) {
    elementMetadata(gr.input) <- inputPeakFile[, 5:ncol(inputPeakFile)]
  }
  if (is.element("strand", colnames(elementMetadata(gr.input)))) {
    elementMetadata(gr.input) <-
      elementMetadata(gr.input)[, -which(colnames(elementMetadata(gr.input)) == "strand")]
  }

  return(gr.input)
}

#### Re arranges repeat Masker file for function convenience ####
FormattingRM <- function(rmsk, goal = "base") {
  options(warn = -1)
  library(biomartr)
  last <-
    as.data.frame(stringr::str_split_fixed(rmsk$matching_class, "/", 2))
  new.rmsk <-
    data.frame(
      "seqnames" = rmsk$qry_id,
      "start" = rmsk$qry_start,
      "end" = rmsk$qry_end,
      "strand" = rmsk$matching_repeat,
      "repeat_name" = rmsk$repeat_id,
      "repeat_type" = last$V1,
      "repeat_family" = last$V2
    )
  new.rmsk$strand <- as.character(new.rmsk$strand)
  new.rmsk$strand <-
    replace(new.rmsk$strand, new.rmsk$strand == "C", "-")

  if (goal == "calcAge") {
    new.rmsk$perc_div <- as.numeric(rmsk$perc_div)
  }

  hit <- new.rmsk$repeat_family == ""
  new.rmsk <- new.rmsk[!hit, ]

  hit2 <- new.rmsk$repeat_type == "Unknown"
  new.rmsk <- new.rmsk[!hit2, ]

  return(new.rmsk)
}

CountRM <- function(rmsk) {
  options(warn = -1)
  library(dplyr)
  x <- rmsk %>% count(repeat_name)
  colnames(x) <- c("RepeatName", "nRepeatName")

  x1 <- rmsk %>% count(repeat_type)
  colnames(x1) <- c("RepeatType", "nRepeatType")

  x2 <- rmsk %>% count(repeat_family)
  colnames(x2) <- c("RepeatFamily", "nRepeatFamily")

  rmsk.counts <- list(x, x2, x1)

  return(rmsk.counts)
}

GetOverlap <-
  function(rmsk,
           gr.input,
           format,
           minoverlap = 0L,
           goal = "base") {
    options(warn = -1)
    #### converts repeat Masker data frame into GRanges object ####

    library(GenomicRanges)
    gr.rmsk <-
      with(rmsk,
           GenomicRanges::GRanges(seqnames, IRanges(start, end), strand = strand))
    elementMetadata(gr.rmsk) <-
      DataFrame(
        RepeatName = rmsk$repeat_name,
        RepeatFamily = rmsk$repeat_family,
        RepeatType = rmsk$repeat_type
      )
    if (goal == "calcAge") {
      elementMetadata(gr.rmsk) <-
        DataFrame(elementMetadata(gr.rmsk), perc_div = as.numeric(rmsk$perc_div))
    }

    #### converts gr.input range into single nucleotide at summit location ####

    if (format == "narrow") {
      if (ncol(elementMetadata(gr.input)) != 0) {
        gr.temp <- gr.input
        newend <-
          GenomicRanges::start(gr.temp) + as.data.frame(gr.temp)[, "summit"] + 1
        newstart <-
          GenomicRanges::start(gr.temp) + as.data.frame(gr.temp)[, "summit"]
        gr.new <-
          with(
            as.data.frame(gr.temp),
            GenomicRanges::GRanges(
              seqnames,
              IRanges(start = as.integer(newstart), end = as.integer(newend)),
              strand = strand
            )
          )
        elementMetadata(gr.new) <- elementMetadata(gr.input)
        gr.input <- gr.new
      }

    }

    #### finds repeat ranges with overlapping summits ####

    if (format == "braod") {
      m <-
        GenomicRanges::findOverlaps(gr.rmsk,
                                    gr.input,
                                    ignore.strand = TRUE,
                                    minoverlap = minoverlap)
    } else{
      m <-
        GenomicRanges::findOverlaps(gr.rmsk, gr.input, ignore.strand = TRUE)
    }
    gr.rmsk.matched <- gr.rmsk[queryHits(m)]

    #### adds the metadata from gr2 to GRanges of intersecting peaks ####

    mcols(gr.rmsk.matched) <- cbind.data.frame(mcols(gr.rmsk.matched),
                                               mcols(gr.input[subjectHits(m)]))
    return(gr.rmsk.matched)
  }

CountElements <- function(gr.rmsk.matched, rmsk) {
  options(warn = -1)
  #### converts matched GRanges object into data frame and counts instances of given column name ####

  df.rmsk.matched <- as.data.frame(gr.rmsk.matched)
  diff.rName <-
    setdiff(rmsk$repeat_name, df.rmsk.matched$RepeatName)
  diff.rFamily <-
    setdiff(rmsk$repeat_family, df.rmsk.matched$RepeatFamily)
  diff.rType <-
    setdiff(rmsk$repeat_type, df.rmsk.matched$RepeatType)

  library(dplyr)
  x <- df.rmsk.matched %>% count(RepeatName)
  colnames(x) <- c("RepeatName", "nRepeatName")
  x <-
    rbind(x, data.frame(RepeatName = diff.rName, nRepeatName = rep.int(0, length(diff.rName))))

  x1 <- df.rmsk.matched %>% count(RepeatFamily)
  colnames(x1) <- c("RepeatFamily", "nRepeatFamily")
  x1 <-
    rbind(x1,
          data.frame(
            RepeatFamily = diff.rFamily,
            nRepeatFamily = rep.int(0, length(diff.rFamily))
          ))

  x2 <- df.rmsk.matched %>% count(RepeatType)
  colnames(x2) <- c("RepeatType", "nRepeatType")
  x2 <-
    rbind(x2,
          data.frame(RepeatType = diff.rType, nRepeatType = rep.int(0, length(diff.rType))))

  all.counts <- list(x, x1, x2)

  ####
  return(all.counts)
}

#### count Intersect description ####
#
#  Takes repeat Master file and annotated narrow Peak file as input
#  Calculates number of repeats with overlapping peaks based on given category
#  (to add: category member filtering, input peak file type, multiple file type support, overlap expectation for broad Peak)
#
#### count Intersect function ####
CountIntersect <-
  function(repeatMaskerFile,
           inputPeakFile,
           format,
           minoverlap = 0L) {
    options(warn = -1)
    gr.input <- inputPeakFile
    rmsk <- FormattingRM(repeatMaskerFile)


    gr.rmsk.matched <- GetOverlap(rmsk,
                                  gr.input, format, minoverlap)

    all.counts <- CountElements(gr.rmsk.matched, rmsk)

    return(all.counts)

  }

ShufflePeaks <- function(peakFile, pathList, seed = 0) {
  options(warn = -1)
  #### get shuffle genome interval with using bedr shuffle function ####

  gr.input <- peakFile
  write.table(
    gr.input,
    file = "before.shuffle.input.peak.bed",
    quote = F,
    sep = "\t",
    row.names = F,
    col.names = T
  )

  gr.Promoter <- gr.input[grepl('Promoter', gr.input$annotation), ]
  gr.Exon <- gr.input[grepl('Exon', gr.input$annotation), ]
  gr.Intron <- gr.input[grepl('Intron', gr.input$annotation), ]
  gr.5UTR <- gr.input[grepl("5' UTR", gr.input$annotation), ]
  gr.3UTR <- gr.input[grepl("3' UTR", gr.input$annotation), ]
  gr.Intergenic <-
    gr.input[grepl('Intergenic', gr.input$annotation), ]
  gr.Downstream <-
    gr.input[grepl('Downstream', gr.input$annotation), ]

  write.table(
    gr.Promoter,
    file = "Promoter.bed",
    quote = F,
    sep = "\t",
    row.names = F,
    col.names = F
  )
  system(
    paste(
      " bedtools shuffle -i Promoter.bed -g ",
      pathList$genomeSizePath,
      " -excl ",
      pathList$Promoter,
      " > shuffled.promoters "
    )
  )
  sh.promoter <-
    read.csv("shuffled.promoters", header = F, sep = "\t")

  write.table(
    as.data.frame(gr.Exon),
    file = "Exon.bed",
    quote = F,
    sep = "\t",
    row.names = F,
    col.names = F
  )
  system(
    paste(
      " bedtools shuffle -i Exon.bed -g ",
      pathList$genomeSizePath,
      " -excl ",
      pathList$Exon,
      " > shuffled.exons "
    )
  )
  sh.exon <- read.csv("shuffled.exons", header = F, sep = "\t")

  write.table(
    as.data.frame(gr.Intron),
    file = "Intron.bed",
    quote = F,
    sep = "\t",
    row.names = F,
    col.names = F
  )
  system(
    paste(
      " bedtools shuffle -i Intron.bed -g ",
      pathList$genomeSizePath,
      " -excl ",
      pathList$Intron,
      " > shuffled.introns "
    )
  )
  sh.intron <- read.csv("shuffled.introns", header = F, sep = "\t")

  write.table(
    as.data.frame(gr.5UTR),
    file = "FiveUTR.bed",
    quote = F,
    sep = "\t",
    row.names = F,
    col.names = F
  )
  system(
    paste(
      " bedtools shuffle -i FiveUTR.bed -g ",
      pathList$genomeSizePath,
      " -excl ",
      pathList$`5UTR`,
      " > shuffled.FiveUTR "
    )
  )
  sh.5UTR <- read.csv("shuffled.FiveUTR", header = F, sep = "\t")

  write.table(
    as.data.frame(gr.3UTR),
    file = "ThreeUTR.bed",
    quote = F,
    sep = "\t",
    row.names = F,
    col.names = F
  )
  system(
    paste(
      " bedtools shuffle -i ThreeUTR.bed -g ",
      pathList$genomeSizePath,
      " -excl ",
      pathList$`3UTR`,
      " > shuffled.ThreeUTR "
    )
  )
  sh.3UTR <- read.csv("shuffled.ThreeUTR", header = F, sep = "\t")

  write.table(
    as.data.frame(gr.Downstream),
    file = "Downstream.bed",
    quote = F,
    sep = "\t",
    row.names = F,
    col.names = F
  )
  system(
    paste(
      " bedtools shuffle -i Downstream.bed -g ",
      pathList$genomeSizePath,
      " -excl ",
      pathList$Downstream,
      " > shuffled.downstream "
    )
  )
  sh.downstream <-
    read.csv("shuffled.downstream", header = F, sep = "\t")

  write.table(
    as.data.frame(gr.Intergenic),
    file = "Intergenic.bed",
    quote = F,
    sep = "\t",
    row.names = F,
    col.names = F
  )
  system(
    paste(
      " bedtools shuffle -i Intergenic.bed -g ",
      pathList$genomeSizePath,
      " -excl ",
      pathList$Promoter,
      " -excl ",
      pathList$Exon,
      " -excl ",
      pathList$Intron,
      " -excl ",
      pathList$`5UTR`,
      " -excl ",
      pathList$`3UTR`,
      " -excl ",
      pathList$Downstream,
      " > shuffled.intergenic "
    )
  )
  sh.intergenic <-
    read.csv("shuffled.intergenic", header = F, sep = "\t")

  all.shuffeledAnnots <-
    rbind(sh.promoter,
          sh.exon,
          sh.intron,
          sh.5UTR,
          sh.3UTR,
          sh.intergenic,
          sh.downstream)
  colnames(all.shuffeledAnnots) <- colnames(as.data.frame(peakFile))
  gr <- MakeGrangeObj(all.shuffeledAnnots)

  return(gr)

}

EnrichPARs <-
  function(inputPeakFile,
           pathList,
           numberOfShuffle = 1,
           repeatMaskerFile,
           format,
           minoverlap = 0L,
           outdir) {
    options(warn = -1)
    # getting input Peak File as grange objects
    gr.input <- inputPeakFile
    # to write peak file as bed file
    write.table(
      as.data.frame(gr.input),
      file = "before.observed.input.peak.bed",
      quote = F,
      sep = "\t",
      row.names = F,
      col.names = T
    )
    # to count matched elements between repeat maskers and peaks
    observe.counts <-
      CountIntersect(repeatMaskerFile, gr.input, format, minoverlap)
    # to call shuffle function
    gr <- ShufflePeaks(inputPeakFile, pathList)
    # to count matched elements between repeat masker and shuffled peaks
    expected.counts <-
      CountIntersect(repeatMaskerFile, gr, format, minoverlap)

    if (numberOfShuffle > 1) {
      # to split different data frames of listed data frame for expected count list
      Rname <- as.data.frame(expected.counts[[1]])
      Rfamily <- as.data.frame(expected.counts[[2]])
      Rtype <- as.data.frame(expected.counts[[3]])

      for (i in 1:numberOfShuffle) {
        start <- Sys.time()
        # to do shuffle each time
        tmp <-  ShufflePeaks(inputPeakFile, pathList)
        end <- Sys.time()
        print(paste("shuffle ", i, " - sh time :", (end - start)))

        start <- Sys.time()
        #to calculate count for each shuffelled peak and repeats
        tmp.counts <-
          CountIntersect(repeatMaskerFile, tmp, format, minoverlap)
        end <- Sys.time()
        print(paste("Countintersect time :", (end - start)))

        # to separate data frames in new counted list
        tmp.Rname <- as.data.frame(tmp.counts[[1]])
        colnames(tmp.Rname) <- c("RepeatName", i)
        tmp.Rfamily <- as.data.frame(tmp.counts[[2]])
        colnames(tmp.Rfamily) <- c("RepeatFamily", i)
        tmp.Rtype <- as.data.frame(tmp.counts[[3]])
        colnames(tmp.Rtype) <- c("RepeatType", i)

        # merge types of counts with new calculated ones
        Rname <- merge(Rname, tmp.Rname, by = "RepeatName")
        Rfamily <- merge(Rfamily, tmp.Rfamily, by = "RepeatFamily")
        Rtype <- merge(Rtype, tmp.Rtype, by = "RepeatType")

      }
      # to calculate means for each data frames and add new column
      Rname$Mean <- round(rowMeans(Rname[, c(2:ncol(Rname))]))
      Rfamily$Mean <- round(rowMeans(Rfamily[, c(2:ncol(Rfamily))]))
      Rtype$Mean <- round(rowMeans(Rtype[, c(2:ncol(Rtype))]))

      # to calculate True means for each data frames and add new column
      Rname$TrueMean <- rowMeans(Rname[, c(2:ncol(Rname))])
      Rfamily$TrueMean <- rowMeans(Rfamily[, c(2:ncol(Rfamily))])
      Rtype$TrueMean <- rowMeans(Rtype[, c(2:ncol(Rtype))])

      # to combine all manipulated data frames for counts results of types as list
      expected.counts <- list(Rname, Rfamily, Rtype)

    }

    # to manipulate getting raw repeat masker file as more regular format
    rmsk <- FormattingRM(repeatMaskerFile)
    # to calculate counts of types in repeat masker
    rmsk.counts <- CountRM(rmsk)

    # to combine count results of each same types for repeat and observed and expected counts

    all.RepeatName <-
      merge(as.data.frame(rmsk.counts[[1]]),
            as.data.frame(observe.counts[[1]]),
            by = "RepeatName")
    all.RepeatName <-
      merge(all.RepeatName, as.data.frame(expected.counts[[1]][, c("RepeatName", "Mean", "TrueMean")]), by = "RepeatName")
    colnames(all.RepeatName) <-
      c("RepeatName", "rmsk", "observed", "expected", "TrueMean")


    all.RepeatFamily <-
      merge(as.data.frame(rmsk.counts[[2]]),
            as.data.frame(observe.counts[[2]]),
            by = "RepeatFamily")
    all.RepeatFamily <-
      merge(all.RepeatFamily, as.data.frame(expected.counts[[2]][, c("RepeatFamily", "Mean", "TrueMean")]), by = "RepeatFamily")
    colnames(all.RepeatFamily) <-
      c("RepeatFamily", "rmsk", "observed", "expected", "TrueMean")

    all.RepeatType <-
      merge(as.data.frame(rmsk.counts[[3]]),
            as.data.frame(observe.counts[[3]]),
            by = "RepeatType")
    all.RepeatType <-
      merge(all.RepeatType, as.data.frame(expected.counts[[3]][, c("RepeatType", "Mean", "TrueMean")]), by = "RepeatType")
    colnames(all.RepeatType) <-
      c("RepeatType", "rmsk", "observed", "expected", "TrueMean")

    # calculate binomial test function
    test <-
      function(x, p, n) {
        binom.test(x, p, n, alternative = "greater", conf.level = 0.95)
      }

    # apply binomial test for each results

    b.rName <-
      mapply(
        test,
        all.RepeatName$observed,
        all.RepeatName$rmsk,
        (all.RepeatName$expected / all.RepeatName$rmsk)
      )
    all.RepeatName$p.value <- do.call(rbind, b.rName["p.value", ])
    all.RepeatName$p.adjust.value <-
      p.adjust(
        all.RepeatName$p.value,
        method = "fdr",
        n = length(all.RepeatName$p.value)
      )
    all.RepeatName$obsOnTrueMean <-
      all.RepeatName$observed / all.RepeatName$TrueMean
    all.RepeatName$p.value[all.RepeatName$observed < 11] <- NA
    all.RepeatName$p.adjust.value[all.RepeatName$observed < 11] <- NA

    b.rFamily <-
      mapply(
        test,
        all.RepeatFamily$observed,
        all.RepeatFamily$rmsk,
        (all.RepeatFamily$expected / all.RepeatFamily$rmsk)
      )
    all.RepeatFamily$p.value <- do.call(rbind, b.rFamily["p.value", ])
    all.RepeatFamily$p.adjust.value <-
      p.adjust(
        all.RepeatFamily$p.value,
        method = "fdr",
        n = length(all.RepeatFamily$p.value)
      )
    all.RepeatFamily$obsOnTrueMean <-
      all.RepeatFamily$observed / all.RepeatFamily$TrueMean
    all.RepeatFamily$p.value[all.RepeatFamily$observed < 11] <- NA
    all.RepeatFamily$p.adjust.value[all.RepeatFamily$observed < 11] <-
      NA

    b.rType <-
      mapply(
        test,
        all.RepeatType$observed,
        all.RepeatType$rmsk,
        (all.RepeatType$expected / all.RepeatType$rmsk)
      )
    all.RepeatType$p.value <- do.call(rbind, b.rType["p.value", ])
    all.RepeatType$p.adjust.value <-
      p.adjust(
        all.RepeatType$p.value,
        method = "fdr",
        n = length(all.RepeatType$p.value)
      )
    all.RepeatType$obsOnTrueMean <-
      all.RepeatType$observed / all.RepeatType$TrueMean
    all.RepeatType$p.value[all.RepeatType$observed < 11] <- NA
    all.RepeatType$p.adjust.value[all.RepeatType$observed < 11] <- NA

    # to sort results as p-adjust values
    all.RepeatName <-
      all.RepeatName[order(all.RepeatName$p.adjust.value), ]
    all.RepeatFamily <-
      all.RepeatFamily[order(all.RepeatFamily$p.adjust.value), ]
    all.RepeatType <-
      all.RepeatType[order(all.RepeatType$p.adjust.value), ]

    # write each type results as separately to output directory
    write.csv(
      all.RepeatName,
      paste0(
        outdir,
        "_",
        numberOfShuffle,
        "Final_Shuffle_beforeSubset_ALLRepeatName.csv"
      ),
      row.names = F,
      quote = F
    )
    write.csv(
      all.RepeatFamily,
      paste0(
        outdir,
        "_",
        numberOfShuffle,
        "Final_Shuffle_beforeSubset_ALLRepeatFamily.csv"
      ),
      row.names = F,
      quote = F
    )
    write.csv(
      all.RepeatType,
      paste0(
        outdir,
        "_",
        numberOfShuffle,
        "Final_Shuffle_beforeSubset_ALLRepeatType.csv"
      ),
      row.names = F,
      quote = F
    )

    # to collect results to a list and returned
    list <-
      list(
        "RepeatName" = all.RepeatName,
        "RepeatFamily" = all.RepeatFamily,
        "RepeatType" = all.RepeatType
      )
    # remove bed shuffle sub log files
    system("rm shuffled.* *.bed")
    return(list)


  }


calculate_background_par <-
  function(df, repeatMaskerFile, peakFile) {
    options(warn = -1)
    df <- df[which(df$p.adjust.value <= 0.05 & df$observed >= 10), ]

    rmsk <- FormattingRM(repeatMaskerFile)

    matched.rmsk <-
      rmsk[which(rmsk$repeat_name %in% as.vector(df$RepeatName)), ]

    matched.rmsk$ID <-
      row.names(matched.rmsk) # for homer unique identifier

    gr.rmsk.par <- GetOverlap(matched.rmsk,
                              peakFile, format = "narrow")

    gr.matched.rmsk <- MakeGrangeObj(matched.rmsk)

    ## extract gr.rmsk not overlapping with gr.overlapped granges
    gr.nonPAR <- gr.matched.rmsk[!gr.matched.rmsk %over% gr.rmsk.par, ]


    list <-
      list("Query" = as.data.frame(gr.rmsk.par)[, c(1, 2, 3, 4, 5, 6, 7, 8, 11)],
           "Background" = as.data.frame(gr.nonPAR))

    return(list)

  }


calculate_background_linkedRepeats <-
  function(df,
           repeatMaskerFile,
           peakFile,
           genes,
           distance) {
    options(warn = -1)
    df <- df[which(df$p.adjust.value <= 0.05 & df$observed >= 10), ]

    rmsk <- FormattingRM(repeatMaskerFile)

    matched.rmsk <-
      rmsk[which(rmsk$repeat_name %in% as.vector(df$RepeatName)), ]

    matched.rmsk$ID <-
      row.names(matched.rmsk) # for homer unique identifier

    gr.rmsk.par <- GetOverlap(matched.rmsk,
                              peakFile, format = "narrow")

    gr.matched.rmsk <- MakeGrangeObj(matched.rmsk)


    d <-
      distanceToNearest(x = gr.matched.rmsk, subject = MakeGrangeObj(genes))
    m <- d[which(elementMetadata(d)$distance < distance),]
    b <- d[which(elementMetadata(d)$distance >= distance),]
    queries <- gr.matched.rmsk[queryHits(m)]
    backgrounds <- gr.matched.rmsk[queryHits(b)]

    list <-
      list("Query" = as.data.frame(queries),
           "Background" = as.data.frame(backgrounds))

    return(list)

  }


FindMotifs <-
  function(df,
           repeatMaskerFile,
           peak,
           genes,
           distance,
           genome,
           outDir,
           homerPath,
           type) {
    options(warn = -1)
    library(marge)
    options("homer_path" = homerPath)
    options(stringsAsFactors = F)

    if (type == "enrichPeak") {
      list <- calculate_background_par(df, repeatMaskerFile, peak)
      queries <- list$Query
      backgrounds <- list$Background
      Rnames <-  unique(queries$RepeatName)[1:10]

      if (length(Rnames) > 0) {
        for (name in Rnames) {
          dir <- paste0(outDir, "/margeOutput/asRepeatName/", name)
          print(dir)
          dir.create(dir, recursive = T)
          find_motifs_genome(queries[which(queries$RepeatName == name), ],
                             dir,
                             genome,
                             overwrite = T,
                             background = backgrounds[which(backgrounds$repeat_name == name), ])
        }
      }
    } else if (type == "linkedRepeats") {
      list <-
        calculate_background_linkedRepeats(df, repeatMaskerFile, peak, genes, distance)
      queries <- list$Query
      backgrounds <- list$Background
      Rnames <-  unique(queries$repeat_name)[1:10]

      if (length(Rnames) > 0) {
        for (name in Rnames) {
          dir <- paste0(outDir, "/margeOutput/asRepeatName/", name)
          print(dir)
          dir.create(dir, recursive = T)
          find_motifs_genome(queries[which(queries$RepeatName == name), ],
                             dir,
                             genome,
                             overwrite = T,
                             background = backgrounds[which(backgrounds$repeat_name == name), ])
        }
      }
    } else{
      return(print("please give correct type paremeter!"))
    }



  }

IdentifyDEGLinkedRepeats <-
  function(enrichPARsResult,
           peaks,
           rmsk,
           genes,
           numberOfShuffle = 100,
           distance = 100000) {
    options(warn = -1)
    count <- enrichPARsResult
    write.csv(count[[1]],
              paste0("count.csv"),
              row.names = F,
              quote = F)
    # getting enriched counts equals and more than 2 for observed repeat counts
    repeatList <-
      count[[1]][which(count[[1]]["observed"] >= 2), c("RepeatName", "observed")]
    gr.input <- peaks
    write.csv(
      as.data.frame(gr.input),
      paste0("gr.input.csv"),
      row.names = F,
      quote = F
    )
    # regulate repeat masker file
    rmsk <- FormattingRM(rmsk)
    # overlapped repeats and peaks
    overlapped <-
      GetOverlap(rmsk = rmsk,
                 gr.input = gr.input,
                 format = "narrow")
    # get overlapped has repeat names matched in enrich results
    repeats.overlapped <-
      subset(overlapped, RepeatName %in% repeatList$RepeatName)
    write.csv(
      repeats.overlapped,
      paste0("repeats.overlapped.csv"),
      row.names = F,
      quote = F
    )
    # find overlapped results nearest to given genes
    d <-
      distanceToNearest(x = repeats.overlapped, subject = MakeGrangeObj(genes))
    write.csv(as.data.frame(d),
              paste0("d.csv"),
              row.names = F,
              quote = F)
    m <- d[which(elementMetadata(d)$distance < distance),]
    result.overlapped <- repeats.overlapped[queryHits(m)]

    # calculate repeats counts related genes
    observed.counts <- CountElements(result.overlapped, rmsk)[[1]]

    gr.rmsk <- MakeGrangeObj(rmsk)

    ## extract gr.rmsk not overlapping with gr.overlapped granges
    gr.nonPAR <- gr.rmsk[gr.rmsk %outside% overlapped, ]
    last.nonPAR <- data.frame()
    isfirsttime <- "true"
    write.csv(
      as.data.frame(gr.nonPAR),
      paste0("gr.nonPAR.csv"),
      row.names = F,
      quote = F
    )

    if (nrow(repeatList) != 0) {
      for (i in 1:nrow(repeatList)) {
        # to get match
        temp.nonPAR <-
          gr.nonPAR[gr.nonPAR$repeat_name == repeatList$RepeatName[i]]
        part.nonPAR <-
          temp.nonPAR[sample(length(temp.nonPAR), repeatList$observed[i]),]
        d <-
          distanceToNearest(x = part.nonPAR, subject = MakeGrangeObj(genes))
        m <- d[which(elementMetadata(d)$distance < distance),]
        part.nonPAR <- part.nonPAR[queryHits(m)]

        if (length(part.nonPAR) != 0) {
          elementMetadata(part.nonPAR)$target <-
            as.character(repeatList$RepeatName[i])
          if (isfirsttime == "true") {
            last.nonPAR <- part.nonPAR
            isfirsttime <- "false"
          } else{
            last.nonPAR <- append(last.nonPAR, part.nonPAR)
          }
        } else if (length(part.nonPAR) == 0 && isfirsttime == "true") {
          last.nonPAR <- part.nonPAR
        }
      }

    }


    #### converts matched GRanges object into data frame and counts instances of given column name ####

    df.nonPAR <- as.data.frame(last.nonPAR)
    write.csv(df.nonPAR,
              paste0("df1.csv"),
              row.names = F,
              quote = F)
    diff.rName <- setdiff(rmsk$repeat_name, df.nonPAR$target)

    library(dplyr)
    x <- data.frame()
    if (nrow(df.nonPAR) != 0) {
      x <- df.nonPAR %>% count(target)
      colnames(x) <- c("RepeatName", "nRepeatName")
    }

    expected.counts <-
      rbind(x, data.frame(RepeatName = diff.rName, nRepeatName = rep.int(0, length(diff.rName)))) #????
    write.csv(expected.counts,
              paste0("ec1.csv"),
              row.names = F,
              quote = F)

    if (numberOfShuffle > 1) {
      for (i in 1:(numberOfShuffle)) {
        last.nonPAR <- data.frame()
        isfirsttime <- "true"
        for (j in 1:nrow(repeatList)) {
          temp.nonPAR <-
            gr.nonPAR[gr.nonPAR$repeat_name == repeatList$RepeatName[j]]
          tmp.nonPAR <-
            temp.nonPAR[sample(length(temp.nonPAR), repeatList$observed[j]),]
          d <-
            distanceToNearest(x = tmp.nonPAR, subject = MakeGrangeObj(genes))
          m <- d[which(elementMetadata(d)$distance < distance),]
          tmp.nonPAR <- tmp.nonPAR[queryHits(m)]
          if (length(tmp.nonPAR) != 0) {
            elementMetadata(tmp.nonPAR)$target <-
              as.character(repeatList$RepeatName[j])
            if (isfirsttime == "true") {
              last.nonPAR <- tmp.nonPAR
              isfirsttime <- "false"
            } else{
              last.nonPAR <- append(last.nonPAR, tmp.nonPAR)
            }
          } else if (length(part.nonPAR) == 0 &&
                     isfirsttime == "true") {
            last.nonPAR <- part.nonPAR
          }
        }

      #### converts matched GRanges object into data frame and counts instances of given column name ####

      df.nonPAR <- as.data.frame(last.nonPAR)
      write.csv(df.nonPAR,
                paste0("df2.csv"),
                row.names = F,
                quote = F)

      diff.rName <- setdiff(rmsk$repeat_name, df.nonPAR$target)
      x <- data.frame()
      if (nrow(df.nonPAR) != 0) {
        x <- df.nonPAR %>% count(target)
        colnames(x) <- c("RepeatName", "nRepeatName")
      }
      tmp <-
        rbind(x,
              data.frame(
                RepeatName = diff.rName,
                nRepeatName = rep.int(0, length(diff.rName))
              )) ###???
      expected.counts <-
        merge(expected.counts, tmp, by = "RepeatName")

    }
    expected.counts$Mean <-
      round(rowMeans(expected.counts[, c(2:ncol(expected.counts))]))
    write.csv(expected.counts,
              paste0("ec2.csv"),
              row.names = F,
              quote = F)

  } else{
    names(expected.counts)[2] <- "Mean"
  }

rmsk.counts <- CountRM(rmsk)

all.RepeatName <-
  merge(as.data.frame(rmsk.counts[[1]]),
        as.data.frame(observed.counts),
        by = "RepeatName")
write.csv(all.RepeatName,
          paste0("ar1.csv"),
          row.names = F,
          quote = F)

all.RepeatName <-
  merge(all.RepeatName, expected.counts[, c("RepeatName", "Mean")], by = "RepeatName")
colnames(all.RepeatName) <-
  c("RepeatName", "rmsk", "observed", "expected")
write.csv(all.RepeatName,
          paste0("ar2.csv"),
          row.names = F,
          quote = F)


test <- function(x, p, n) {
  binom.test(x, p, n)
}

b.rName <-
  mapply(
    test,
    all.RepeatName$observed,
    all.RepeatName$rmsk,
    (all.RepeatName$expected / all.RepeatName$rmsk)
  )
all.RepeatName$p.value <- do.call(rbind, b.rName["p.value", ])
all.RepeatName$p.adjust.value <-
  p.adjust(all.RepeatName$p.value,
           method = "fdr",
           n = length(all.RepeatName$p.value))

all.RepeatName <-
  all.RepeatName[order(all.RepeatName$p.adjust.value), ]

return(all.RepeatName)
}

EstimateRepeatAge <-
  function(repeatMasterFile, peakFile, substRate) {
    options(warn = -1)
    rmsk <- FormattingRM(repeatMasterFile, goal = "calcAge")
    gr.rmsk <- MakeGrangeObj(rmsk)
    gr.peak <- peakFile

    overlapped <-
      GetOverlap(
        rmsk = rmsk,
        gr.input = gr.peak,
        format = "narrow",
        goal = "calcAge"
      )

    df.PAR <- as.data.frame(overlapped)

    ## extract gr.rmsk not overlapping with gr.overlapped granges
    gr.nonPAR <- gr.rmsk[!gr.rmsk %over% overlapped, ]
    df.nonPAR <- as.data.frame(gr.nonPAR)

    ## calculate mean of perc_div for each repeat
    x1 <- data.frame()
    if (nrow(df.PAR) != 0) {
      x1 <- df.PAR[, c("RepeatName", "perc_div")] %>%
        group_by(RepeatName) %>%
        summarise(Age = (mean(perc_div)))
    }
    x1$Type <- "PAR"

    x2 <- data.frame()
    if (nrow(df.nonPAR) != 0) {
      x2 <- df.nonPAR[, c("repeat_name", "perc_div")] %>%
        group_by(repeat_name) %>%
        summarise(Age = (mean(perc_div)))
    }
    x2$Type <- "nonPAR"
    colnames(x2)[1] <- "RepeatName"

    x <- rbind(x1, x2)

    ## calculate age of each repeat when annotated on genome
    x$Age <- (x$Age * 10) / substRate

    return(x)

  }

getInterval <- function(input, dataset) {
  options(warn = -1)
  library(biomaRt)
  ensembl = biomaRt::useEnsembl(biomart = "ensembl",
                                dataset = dataset,
                                verbose = FALSE)
  df <- scan(input, character())

  gene_ids <-
    getBM(
      attributes = c("ensembl_gene_id", "external_gene_name"),
      filters = c("external_gene_name"),
      values = df,
      mart = ensembl,
      verbose = T,
      useCache = FALSE
    )

  data <-
    getBM(
      attributes = c(
        "ensembl_gene_id",
        "external_gene_name",
        "chromosome_name",
        "start_position",
        "end_position",
        "strand"
      ),
      filters = c("ensembl_gene_id"),
      # listFilters(ensembl)
      values = gene_ids$ensembl_gene_id,
      mart = ensembl,
      verbose = FALSE,
      useCache = FALSE
    )

  names(data) <-
    c("geneID", "geneName", "seqnames", "start", "end", "strand")

  hit <- which(data$strand == "1")
  data$strand[hit] <- "+"

  hit <- which(data$strand == "-1")
  data$strand[hit] <- "-"

  data$seqnames <-
    sapply(1:nrow(data), function(x)
      gsub("^\\d$", paste0("chr", data$seqnames[x]), data$seqnames[x]))

  return(data)
}

MakePlots <- function(enrichPARsResult, outDir) {
  options(warn = -1)
  MakePlotsRepeatName(enrichPARsResult$RepeatName, outDir)
  MakePlotsRepeatType(enrichPARsResult$RepeatType, outDir)
  MakePlotsRepeatFamily(enrichPARsResult$RepeatFamily, outDir)


}
