
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
           inputPeakFile, format, minoverlap=0L) {

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

    if(format == "narrow"){

      if(ncol(elementMetadata(gr.input))!= 0){
        GenomicRanges::start(gr.input) <- GenomicRanges::start(gr.input) + gr.input$peak
        GenomicRanges::end(gr.input) <- GenomicRanges::start(gr.input) + 1
      }

    }

    #### finds repeat ranges with overlapping summits ####

    if(format == "braod"){
      m <- GenomicRanges::findOverlaps(gr.rmsk, gr.input, ignore.strand = TRUE, minoverlap = minoverlap)
    }else{
      m <- GenomicRanges::findOverlaps(gr.rmsk, gr.input, ignore.strand = TRUE)
    }
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



ShufflePeaks <- function(peakFile, pathList, seed = 0){

  gr.input <- MakeGrangeObj(inputPeakFile = peakFile)

  gr.Promoter <- gr.input[grepl('Promoter', gr.input$annotation),]
  gr.Exon <- gr.input[grepl('Exon', gr.input$annotation),]
  gr.Intron <- gr.input[grepl('Intron', gr.input$annotation),]
  gr.5UTR <- gr.input[grepl("5' UTR", gr.input$annotation),]
  gr.3UTR <- gr.input[grepl("3' UTR", gr.input$annotation),]
  gr.Intergenic <- gr.input[grepl('Intergenic', gr.input$annotation),]
  gr.Downstream <- gr.Intron <- gr.input[grepl('Downstream', gr.input$annotation),]


  write.table(as.data.frame(gr.Promoter), file="Promoter.bed", quote=F, sep="\t", row.names=F, col.names=F)
  system(paste(" bedtools shuffle -i Promoter.bed -g ",pathList$genomeSizePath," -incl ",pathList$Promoter," > shuffled.promoters "))
  sh.promoter <- read.csv("shuffled.promoters", header = F, sep = "\t")

  write.table(as.data.frame(gr.Exon), file="Exon.bed", quote=F, sep="\t", row.names=F, col.names=F)
  system(paste(" bedtools shuffle -i Exon.bed -g ",pathList$genomeSizePath," -incl ",pathList$Exon," > shuffled.exons "))
  sh.exon <- read.csv("shuffled.exons", header = F, sep = "\t")

  write.table(as.data.frame(gr.Intron), file="Intron.bed", quote=F, sep="\t", row.names=F, col.names=F)
  system(paste(" bedtools shuffle -i Intron.bed -g ",pathList$genomeSizePath," -incl ",pathList$Intron," > shuffled.introns "))
  sh.intron <- read.csv("shuffled.introns", header = F, sep = "\t")

  write.table(as.data.frame(gr.5UTR), file="FiveUTR.bed", quote=F, sep="\t", row.names=F, col.names=F)
  system(paste(" bedtools shuffle -i FiveUTR.bed -g ",pathList$genomeSizePath," -incl ",pathList$`5UTR`," > shuffled.5UTR "))
  sh.5UTR <- read.csv("shuffled.5UTR", header = F, sep = "\t")

  write.table(as.data.frame(gr.3UTR), file="ThreeUTR.bed", quote=F, sep="\t", row.names=F, col.names=F)
  system(paste(" bedtools shuffle -i ThreeUTR.bed -g ",pathList$genomeSizePath," -incl ",pathList$`3UTR`," > shuffled.3UTR "))
  sh.3UTR <- read.csv("shuffled.3UTR", header = F, sep = "\t")

  write.table(as.data.frame(gr.Downstream), file="Downstream.bed", quote=F, sep="\t", row.names=F, col.names=F)
  system(paste(" bedtools shuffle -i Downstream.bed -g ",pathList$genomeSizePath," -incl ",pathList$Downstream," > shuffled.downstream "))
  sh.downstream <- read.csv("shuffled.downstream", header = F, sep = "\t")

  write.table(as.data.frame(gr.Intergenic), file="Intergenic.bed", quote=F, sep="\t", row.names=F, col.names=F)
  system(paste(" bedtools shuffle -i Intergenic.bed -g ",pathList$genomeSizePath," -excl ",pathList$Promoter,
               " -excl ",pathList$Exon," -excl ",pathList$Intron," -excl ",pathList$`5UTR`," -excl ",pathList$`3UTR`,
               " -excl ",pathList$Downstream," > shuffled.intergenic "))
  sh.intergenic <- read.csv("shuffled.intergenic", header = F, sep = "\t")


  all.shuffeledAnnots <- rbind(sh.promoter, sh.exon, sh.intron, sh.5UTR, sh.3UTR, sh.intergenic, sh.downstream)
  colnames(all.shuffeledAnnots) <- colnames(peakFile)
  gr <- MakeGrangeObj(all.shuffeledAnnots)

  return(gr)

}




#### get shuffle genome interval with using bedr shuffle function ####

EnrichPARs <- function(inputPeakFile, pathList, numberOfShuffle=1, repeatMaskerFile, format, minoverlap=0L){

  gr.input <- MakeGrangeObj(inputPeakFile = inputPeakFile)
  observe.counts <- CountIntersect(repeatMaskerFile, gr.input,format, minoverlap)

  gr <- ShufflePeaks(inputPeakFile, pathList)
  expected.counts <- CountIntersect(repeatMaskerFile, gr, format, minoverlap)

  if(numberOfShuffle > 1){

    Rname <- as.data.frame(expected.counts[[1]])
    Rfamily <- as.data.frame(expected.counts[[2]])
    Rtype <- as.data.frame(expected.counts[[3]])

    for(i in 1:numberOfShuffle){

      tmp <-  ShufflePeaks(inputPeakFile, pathList)
      tmp.counts <- CountIntersect(repeatMaskerFile, tmp, format, minoverlap)

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
  all.RepeatName <- merge(all.RepeatName, as.data.frame(expected.counts[[1]][,c(1,ncol(expected.counts[[1]]))]), by = "RepeatName")
  colnames(all.RepeatName) <- c("RepeatName","rmsk","observe","expected")


  all.RepeatFamily <- merge(as.data.frame(rmsk.counts[[2]]),as.data.frame(observe.counts[[2]]), by = "RepeatFamily")
  all.RepeatFamily <- merge(all.RepeatFamily, as.data.frame(expected.counts[[2]][,c(1,ncol(expected.counts[[1]]))]), by = "RepeatFamily")
  colnames(all.RepeatFamily) <- c("RepeatFamily","rmsk","observe","expected")

  all.RepeatType <- merge(as.data.frame(rmsk.counts[[3]]),as.data.frame(observe.counts[[3]]), by = "RepeatType")
  all.RepeatType <- merge(all.RepeatType, as.data.frame(expected.counts[[3]][,c(1,ncol(expected.counts[[1]]))]), by = "RepeatType")
  colnames(all.RepeatType) <- c("RepeatType","rmsk","observe","expected")


  test <- function(x, p, n){binom.test(x, p, n)}

  b.rName <- mapply(test, all.RepeatName$observe, all.RepeatName$rmsk, (all.RepeatName$expected/all.RepeatName$rmsk))
  all.RepeatName$p.value <- do.call(rbind, b.rName["p.value",])
  all.RepeatName$p.adjust.value <- p.adjust(all.RepeatName$p.value, method = "fdr", n = length(all.RepeatName$p.value))

  b.rFamily <- mapply(test,all.RepeatFamily$observe, all.RepeatFamily$rmsk, (all.RepeatFamily$expected/all.RepeatFamily$rmsk))
  all.RepeatFamily$p.value <- do.call(rbind, b.rFamily["p.value",])
  all.RepeatFamily$p.adjust.value <- p.adjust(all.RepeatFamily$p.value, method = "fdr", n = length(all.RepeatFamily$p.value))

  b.rType <- mapply(test, all.RepeatType$observe, all.RepeatType$rmsk, (all.RepeatType$expected/all.RepeatType$rmsk))
  all.RepeatType$p.value <- do.call(rbind, b.rType["p.value",])
  all.RepeatType$p.adjust.value <- p.adjust(all.RepeatType$p.value, method = "fdr", n = length(all.RepeatType$p.value))


  binom.test.results <- list("RepeatName" = subset(all.RepeatName, p.value < 1e-03 & observe > expected), "RepeatFamily" = subset(all.RepeatFamily, p.value < 1e-03 & observe > expected), "RepeatType" = subset(all.RepeatType, p.value < 1e-03 & observe > expected))

  return(binom.test.results)

}


####


FindMotifs <- function(df ,repeatMaskerFile, outDir, homerPath){

  rmsk <- FormattingRM(repeatMaskerFile)
  binom.test.results <- df

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

  library(marge)
  options("homer_path" = homerPath)
  options(stringsAsFactors = F)

  if(nrow(df.name) > 0){

    find_motifs_genome(df.name, paste0(outDir,"/margeOutput/asRepeatName"), "hg19", overwrite = T)

  }

  if(nrow(df.family) > 0){

    find_motifs_genome(df.family, paste0(outDir,"/margeOutput/asRepeatFamily"), "hg19", overwrite = T)

  }

  if(nrow(df.type) > 0){

    find_motifs_genome(df.type, paste0(outDir,"/margeOutput/asRepeatType"), "hg19", overwrite = T)

  }


}







