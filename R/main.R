
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

#### rearranges repeatMasker file for function convenience ####
FormattingRM <- function(rmsk){

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

  hit <- rmsk$repeat_family == ""
  rmsk <- rmsk[!hit,]

  hit2 <- rmsk$repeat_type == "Unknown"
  rmsk <- rmsk[!hit2,]

  return(rmsk)
}

CountRM <- function(rmsk){

  library(dplyr)
  x <- rmsk %>% count(repeat_name)
  colnames(x) <- c("RepeatName", "nRepeatName")

  x1 <- rmsk %>% count(repeat_type)
  colnames(x1) <- c("RepeatType", "nRepeatType")

  x2 <- rmsk %>% count(repeat_family)
  colnames(x2) <- c("RepeatFamily","nRepeatFamily")

  rmsk.counts <- list(x,x2,x1)

  return(rmsk.counts)
}

CountIntersect <-
  function(repeatMaskerFile,
           inputPeakFile) {

    gr.input <- inputPeakFile
    rmsk <- FormattingRM(repeatMaskerFile)

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
    diff.rName <- setdiff(rmsk$repeat_name, df.rmsk.matched$RepeatName)
    diff.rFamily <- setdiff(rmsk$repeat_family, df.rmsk.matched$RepeatFamily)
    diff.rType <- setdiff(rmsk$repeat_type, df.rmsk.matched$RepeatType)

    library(dplyr)
    x <- df.rmsk.matched %>% count(RepeatName)
    colnames(x) <- c("RepeatName", "nRepeatName")
    x <- rbind(x, data.frame(RepeatName = diff.rName, nRepeatName = rep.int(0, length(diff.rName))))

    x1 <- df.rmsk.matched %>% count(RepeatFamily)
    colnames(x1) <- c("RepeatFamily","nRepeatFamily")
    x1 <- rbind(x1, data.frame(RepeatFamily = diff.rFamily, nRepeatFamily = rep.int(0, length(diff.rFamily))))

    x2 <- df.rmsk.matched %>% count(RepeatType)
    colnames(x2) <- c("RepeatType", "nRepeatType")
    x2 <- rbind(x2, data.frame(RepeatType = diff.rType, nRepeatType = rep.int(0, length(diff.rType))))

    all.counts <- list(x,x1,x2)

    ####

    return(all.counts)

  }


pathList <- list("Promoter" = "/home/nazmiye/Desktop/BIP/test/hg19.Promoter",
     "Exon" = "/home/nazmiye/Desktop/BIP/test/hg19.Exons",
     "Intron" = "/home/nazmiye/Desktop/BIP/test/hg19.Introns",
     "5UTR" = "/home/nazmiye/Desktop/BIP/test/hg19.5Prime",
     "3UTR" = "/home/nazmiye/Desktop/BIP/test/hg19.3Prime",
     "Downstream" = "/home/nazmiye/Desktop/BIP/test/hg19.Downstream",
     "genomeSizePath" = "/home/nazmiye/Desktop/BIP/test/hg19.chrom.sizes"
  )


ShufflePeaks <- function(peakFile, pathList, seed = 0){

  gr.input <- inputPeakFile

  gr.Promoter <- gr.input[grepl('Promoter', gr.input$annotation),]
  gr.Exon <- gr.input[grepl('Exon', gr.input$annotation),]
  gr.Intron <- gr.input[grepl('Intron', gr.input$annotation),]
  gr.5UTR <- gr.input[grepl("5' UTR", gr.input$annotation),]
  gr.3UTR <- gr.input[grepl("3' UTR", gr.input$annotation),]
  gr.Intergenic <- gr.input[grepl('Intergenic', gr.input$annotation),]
  gr.Downstream <- gr.Intron <- gr.input[grepl('Downstream', gr.input$annotation),]

#   library(valr)
#   genome <- read_genome(genomeSizePath)

  write.table(as.data.frame(gr.Promoter), file="inc.Promoter.bed", quote=F, sep="\t", row.names=F, col.names=F)
  system(paste(" bedtools shuffle -i ",pathList$Promoter," -g ",pathList$genomeSizePath," -incl inc.Promoter.bed > shuffled.promoters "))

  sh.Promoter <- bed_shuffle(gr.Promoter, genome, incl = tbl.promoter)
  colnames(sh.Promoter)<-c("seqnames","start","end","strand")

  sh.Exon <- bed_shuffle(gr.Exon, genome)
  colnames(sh.Exon)<-c("seqnames","start","end","strand")

  sh.Intron <- bed_shuffle(gr.Intron, genome)
  colnames(sh.Intron)<-c("seqnames","start","end","strand")

  sh.5UTR <- bed_shuffle(gr.5UTR, genome)
  colnames(sh.5UTR)<-c("seqnames","start","end","strand")

  sh.3UTR <- bed_shuffle(gr.3UTR, genome)
  colnames(sh.3UTR)<-c("seqnames","start","end","strand")

  sh.Intergenic <- bed_shuffle(gr.Intergenic, genome)
  colnames(sh.Intergenic)<-c("seqnames","start","end","strand")

  sh.Downstream <- bed_shuffle(gr.Downstream, genome)
  colnames(sh.Downstream)<-c("seqnames","start","end","strand")




  gr <- MakeGrangeObj(gr)

  return(gr)

}




#### get shuffle genome interval with using bedr shuffle function ####

EnrichPARs <- function(inputPeakFile, ShuffledPeak, genomeSizePath, numberOfShuffle=1, repeatMaskerFile ){

  gr.input <- inputPeakFile
  observe.counts <- CountIntersect(repeatMaskerFile, gr.input)

  gr <- ShuffledPeak
  expected.counts <- CountIntersect(repeatMaskerFile, gr)

  if(numberOfShuffle > 1){

    Rname <- as.data.frame(expected.counts[[1]])
    Rfamily <- as.data.frame(expected.counts[[2]])
    Rtype <- as.data.frame(expected.counts[[3]])

    for(i in 1:numberOfShuffle){

      tmp <- bed_shuffle(gr.input, genome)
      colnames(tmp)<-c("seqnames","start","end","strand")
      tmp <- MakeGrangeObj(tmp)
      tmp.counts <- CountIntersect(repeatMaskerFile, tmp)

      tmp.Rname <- as.data.frame(tmp.counts[[1]])
      colnames(tmp.Rname) <- c("RepeatName",i)
      tmp.Rfamily <- as.data.frame(tmp.counts[[2]])
      colnames(tmp.Rfamily) <- c("RepeatFamily",i)
      tmp.Rtype <- as.data.frame(tmp.counts[[3]])
      colnames(tmp.Rtype) <- c("RepeatType",i)


      Rname <- merge(Rname, tmp.Rname, by = "RepeatName")
      Rfamily <- merge(Rfamily, tmp.Rfamily, by = "RepeatFamily")
      Rtype <- merge(Rtype, tmp.Rtype, by = "RepeatType")

    }

    Rname$Mean <- round(rowMeans(Rname[,c(2:ncol(Rname))]))
    Rfamily$Mean <- round(rowMeans(Rfamily[,c(2:ncol(Rfamily))]))
    Rtype$Mean <- round(rowMeans(Rtype[,c(2:ncol(Rtype))]))

    expected.counts <- list(Rname,Rfamily,Rtype)

  }

  rmsk <- FormattingRM(rmsk)
  rmsk.counts <- CountRM(rmsk)


  all.RepeatName <- merge(as.data.frame(rmsk.counts[[1]]),as.data.frame(observe.counts[[1]]), by = "RepeatName")
  all.RepeatName <- merge(all.RepeatName, as.data.frame(expected.counts[[1]][,c(1,6)]), by = "RepeatName")
  colnames(all.RepeatName) <- c("RepeatName","rmsk","observe","expected")


  all.RepeatFamily <- merge(as.data.frame(rmsk.counts[[2]]),as.data.frame(observe.counts[[2]]), by = "RepeatFamily")
  all.RepeatFamily <- merge(all.RepeatFamily, as.data.frame(expected.counts[[2]][,c(1,6)]), by = "RepeatFamily")
  colnames(all.RepeatFamily) <- c("RepeatFamily","rmsk","observe","expected")

  all.RepeatType <- merge(as.data.frame(rmsk.counts[[3]]),as.data.frame(observe.counts[[3]]), by = "RepeatType")
  all.RepeatType <- merge(all.RepeatType, as.data.frame(expected.counts[[3]][,c(1,6)]), by = "RepeatType")
  colnames(all.RepeatType) <- c("RepeatType","rmsk","observe","expected")


  test <- function(x, p, n){binom.test(x, p, n)}

  b.rName <- mapply(test, all.RepeatName$observe, all.RepeatName$rmsk, (all.RepeatName$expected/all.RepeatName$rmsk))
  all.RepeatName$p.value <- do.call(rbind, b.rName["p.value",])

  b.rFamily <- mapply(test,all.RepeatFamily$observe, all.RepeatFamily$rmsk, (all.RepeatFamily$expected/all.RepeatFamily$rmsk))
  all.RepeatFamily$p.value <- do.call(rbind, b.rFamily["p.value",])

  b.rType <- mapply(test, all.RepeatType$observe, all.RepeatType$rmsk, (all.RepeatType$expected/all.RepeatType$rmsk))
  all.RepeatType$p.value <- do.call(rbind, b.rType["p.value",])


  binom.test.results <- list("RepeatName" = subset(all.RepeatName, p.value < 1e-03 & observe > expected), "RepeatFamily" = subset(all.RepeatFamily, p.value < 1e-03 & observe > expected), "RepeatType" = subset(all.RepeatType, p.value < 1e-03 & observe > expected))

  return(binom.test.results)

}


  ####


homer <- function(df ,repeatMaskerFİle, class = "name", outDir){

  rmsk <- repeatMaskerFile
  df <- binom.test.results

  name <- rmsk[which(rmsk$repeat_name %in% as.vector(binom.test.results$RepeatName$RepeatName)),]
  family <- rmsk[which(rmsk$repeat_family %in% as.vector(binom.test.results$RepeatFamily$RepeatFamily)),]
  type <- rmsk[which(rmsk$repeat_type %in% as.vector(binom.test.results$RepeatType$RepeatType)),]

#   all.annot <- rbind(name,family,type)
#   all.annot$ID <- row.names(all.annot)
#
#   df <- all.annot[,c(1,2,3,8,5,4)]

  name$ID <- row.names(name)
  df.name <- name[,c(1,2,3,8,5,4)]

  family$ID <- row.names(family)
  df.family <- family[,c(1,2,3,8,5,4)]

  type$ID <- row.names(type)
  df.type <- type[,c(1,2,3,8,5,4)]


  write.table(df.name, file = paste0(outDir,"/repeatName.bed"), quote=F, sep="\t", row.names=F, col.names=F,)
  write.table(df.family, file = paste0(outDir,"/repeatFamily.bed"), quote=F, sep="\t", row.names=F, col.names=F)
  write.table(df.type, file = paste0(outDir,"/repeatType.bed"), quote=F, sep="\t", row.names=F, col.names=F)



  if(nrow(df.name) > 0){

    dir.create(path = paste0(outDir,"/Name"))
    system(paste("findMotifsGenome.pl ",paste0(outDir,"/Name/repeatName.bed"), "hg19 . -size 100 "))

  }

  if(nrow(df.family) > 0){

    dir.create(path = paste0(outDir,"/Family"))
    system(paste("findMotifsGenome.pl ",paste0(outDir,"/Family/repeatFamily.bed"), "hg19 . -size 100 "))

  }

  if(nrow(df.type) > 0){

    dir.create(path = paste0(outDir,"/Type"))
    system(paste("findMotifsGenome.pl ",paste0(outDir,"/Type/repeatType.bed"), "hg19 . -size 100 "))

  }






}








