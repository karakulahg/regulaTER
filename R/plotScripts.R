
options(stringsAsFactors = F)
library(ggplot2)


MakePlotsRepeatName <- function(repeat.names, outDir, width, height){

  #### plotting script ####

  # removes unplottable values

  repeat.names.sub <-
    repeat.names[!is.na(repeat.names$p.adjust.value), ]

  # adds fraction of contributing peaks

  repeat.names.sub$contributingRatio = repeat.names.sub$observed / repeat.names.sub$rmsk

  # adds significance labels

  repeat.names.sub$is.significant <- FALSE
  repeat.names.sub$is.significant[which(repeat.names.sub$p.adjust.value <= 0.05)] <-
    TRUE

  # plots the p.value vs ratio plot

  if (nrow(repeat.names.sub) != 0){

    ggplot(repeat.names.sub,
           aes(obsOnTrueMean, -log10(p.adjust.value), colour = is.significant)) + geom_point() +
      theme_light() + geom_text(aes(label = ifelse(is.significant == TRUE, as.character(RepeatName), '')), hjust = 0.5, vjust = -.7) +
      labs(x = "Fold Enrichment", y = "- log 10 (FDR)", colour = "Significance")

    # ggplot(repeat.names.sub,
    #        aes(obsOnTrueMean, -log10(p.adjust.value), color = is.significant)) + geom_point(size = 7) +
    #   theme_classic() +
    #   labs(x = "Fold Enrichment", y = "- log 10 (FDR)", colour = "Significance") +
    #   scale_color_manual(values=c("#253494", "#a50f15")) + theme(
    #     text = element_text( color = "grey20"))


    ggsave(paste0(outDir,"/Fig1_RepeatName.pdf"), width = width, height = height)

    # subsets significant repeats
    repeat.names.sub[repeat.names.sub$is.significant == TRUE,] -> repeat.names.sub.2

    if (nrow(repeat.names.sub.2) != 0){

      # finds coefficient between y-axes

      coeff = 1 / (max(repeat.names.sub.2$observed) / max(repeat.names.sub.2$contributingRatio))


      ggplot(repeat.names.sub[repeat.names.sub$is.significant == TRUE,], aes(x = RepeatName, group = 1)) +
        geom_col(aes(y = observed)) + geom_line(aes(y = (observed / rmsk) / coeff)) + theme_light() + labs(x = "Repeat Name") +
        scale_y_continuous(name = "Observed Repeats", sec.axis = sec_axis( ~ . * coeff, name = "Ratio of Observed Repeats")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.7)) + theme(axis.text=element_text(size=6))


      # ggplot(repeat.names.sub[repeat.names.sub$is.significant == TRUE,], aes(x = RepeatName, group = 1)) +
      #   geom_col(fill = "#253494",aes(y = observed)) + geom_line(color = "#737373",aes(y = (observed / rmsk) / coeff)) + theme_classic() + labs(x = "Repeat Name") +
      #   scale_y_continuous(name = "Observed Repeats", sec.axis = sec_axis( ~ . * coeff, name = "Ratio of Observed Repeats")) +
      #   theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.7)) + theme(axis.text=element_text(size=6))


      # plots the peak contribution plot
      ggsave(paste0(outDir,"/Fig2_RepeatName.pdf"), width = width, height = height)

    }
  }
}

MakePlotsRepeatType <- function(repeat.types, outDir, width, height){

  #### plotting script ####

  # removes unplottable values

  repeat.types.sub <-
    repeat.types[!is.na(repeat.types$p.adjust.value), ]

  # adds fraction of contributing peaks

  repeat.types.sub$contributingRatio = repeat.types.sub$observed / repeat.types.sub$rmsk

  # adds significance labels

  repeat.types.sub$is.significant <- FALSE
  repeat.types.sub$is.significant[which(repeat.types.sub$p.adjust.value <= 0.05)] <-
    TRUE

  # plots the p.value vs ratio plot

  if (nrow(repeat.types.sub) != 0){

    ggplot(repeat.types.sub,
           aes(obsOnTrueMean, -log10(p.adjust.value), colour = is.significant)) + geom_point() + theme_light() +
      geom_text(aes(label =ifelse(is.significant == TRUE, as.character(RepeatType), '')), hjust = 0.5, vjust = -.7) +
      labs(x = "Fold Enrichment", y = "- log 10 (FDR)", colour = "Significance")


    # ggplot(repeat.types.sub,
    #        aes(obsOnTrueMean, -log10(p.adjust.value), color = is.significant)) + geom_point(size = 7) +
    #   theme_classic() +
    #   labs(x = "Fold Enrichment", y = "- log 10 (FDR)", colour = "Significance") +
    #   scale_color_manual(values=c("#253494", "#a50f15")) + theme(
    #     text = element_text( color = "grey20"))

    ggsave(paste0(outDir,"/Fig1_RepeatType.pdf"), width = width, height = height)

    # subsets significant repeats
    repeat.types.sub[repeat.types.sub$is.significant == TRUE,] -> repeat.types.sub.2

    if (nrow(repeat.types.sub.2) != 0){

      # finds coefficient between y-axes

      coeff = 1 / (max(repeat.types.sub.2$observed) / max(repeat.types.sub.2$contributingRatio))

      # plots the peak contribution plot

      ggplot(repeat.types.sub[repeat.types.sub$is.significant == TRUE,], aes(x = RepeatType, group = 1)) + geom_col(aes(y = observed)) +
        geom_line(aes(y = (observed / rmsk) / coeff)) + theme_light() + labs(x = "Repeat Name") +
        scale_y_continuous(name = "Observed Repeats", sec.axis = sec_axis( ~ . * coeff, name = "Ratio of Observed Repeats")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.7)) + theme(axis.text=element_text(size=6))

      # ggplot(repeat.types.sub.2[repeat.types.sub.2$is.significant == TRUE,], aes(x = RepeatType, group = 1)) +
      #   geom_col(fill = "#253494",aes(y = observed)) + geom_line(color = "#737373",aes(y = (observed / rmsk) / coeff)) + theme_classic() + labs(x = "Repeat Type") +
      #   scale_y_continuous(name = "Observed Repeats", sec.axis = sec_axis( ~ . * coeff, name = "Ratio of Observed Repeats")) +
      #   theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.7)) + theme(axis.text=element_text(size=6))



      ggsave(paste0(outDir,"/Fig2_RepeatType.pdf"), width = width, height = height)

    }
  }
}

MakePlotsRepeatFamily <- function(repeat.families, outDir, width, height){

  #### plotting script ####

  # removes unplottable values

  repeat.families.sub <-
    repeat.families[!is.na(repeat.families$p.adjust.value), ]

  # adds fraction of contributing peaks

  repeat.families.sub$contributingRatio = repeat.families.sub$observed / repeat.families.sub$rmsk

  # adds significance labels

  repeat.families.sub$is.significant <- FALSE
  repeat.families.sub$is.significant[which(repeat.families.sub$p.adjust.value <= 0.05)] <-
    TRUE

  # plots the p.value vs ratio plot

  if (nrow(repeat.families.sub) != 0){

    ggplot(repeat.families.sub,
           aes(obsOnTrueMean, -log10(p.adjust.value), colour = is.significant)) + geom_point() + theme_light() +
      geom_text(aes(label = ifelse( is.significant == TRUE, as.character(RepeatFamily), '')), hjust = 0.5, vjust = -.7) +
      labs(x = "Fold Enrichment", y = "- log 10 (FDR)", colour = "Significance")



    # ggplot(repeat.families.sub,
    #        aes(obsOnTrueMean, -log10(p.adjust.value), color = is.significant)) + geom_point(size = 7) +
    #   theme_classic() +
    #   labs(x = "Fold Enrichment", y = "- log 10 (FDR)", colour = "Significance") +
    #   scale_color_manual(values=c("#253494", "#a50f15")) + theme(
    #     text = element_text( color = "grey20"))

    ggsave(paste0(outDir,"/Fig1_RepeatFamily.pdf"), width = width, height = height)

    # subsets significant repeats
    repeat.families.sub[repeat.families.sub$is.significant == TRUE,] -> repeat.families.sub.2

    if (nrow(repeat.families.sub.2) != 0){

      # finds coefficient between y-axes

      coeff = 1 / (max(repeat.families.sub.2$observed) / max(repeat.families.sub.2$contributingRatio))

      # plots the peak contribution plot

      ggplot(repeat.families.sub[repeat.families.sub$is.significant == TRUE,], aes(x = RepeatFamily, group = 1)) + geom_col(aes(y = observed)) +
        geom_line(aes(y = (observed / rmsk) / coeff)) + theme_light() + labs(x = "Repeat Name") +
        scale_y_continuous(name = "Observed Repeats",sec.axis = sec_axis( ~ . * coeff, name = "Ratio of Observed Repeats")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.7)) + theme(axis.text=element_text(size=6))


      # ggplot(repeat.families.sub.2[repeat.families.sub.2$is.significant == TRUE,], aes(x = RepeatFamily, group = 1)) +
      #   geom_col(fill = "#253494",aes(y = observed)) + geom_line(color = "#737373",aes(y = (observed / rmsk) / coeff)) + theme_classic() + labs(x = "Repeat Family") +
      #   scale_y_continuous(name = "Observed Repeats", sec.axis = sec_axis( ~ . * coeff, name = "Ratio of Observed Repeats")) +
      #   theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.7)) + theme(axis.text=element_text(size=6))




      ggsave(paste0(outDir,"/Fig2_RepeatFamily.pdf"), width = width, height = height)
    }
  }
}

