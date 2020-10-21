options(stringsAsFactors = F)

library(ggplot2)

Filtering <- function(df,stage){
  df[df$observed < 10 |
       df$expected == 0, "p.value"] <- NA
  df[df$observed < 10 |
       df$expected == 0, "p.adjust.value"] <- NA

  df$log.p <- -log10(df$p.adjust.value)

  df$isSignificant <-
    ifelse(df$p.adjust.value <= 0.05 &
             df$observed >= 10,
           "True",
           "False")

  df$stage <- stage
  return(df)
}


#### Fig 1 ####

MakeFigure1 <- function(df, outDir, type){
  order(df$log.p, decreasing = T)[1:15] -> top15.df
  if(type == "RepeatName"){
    ggplot(df[!is.na(df$p.adjust.value),],
           aes(x = obsOnTrueMean, y = log.p)) +
      geom_point(aes(color = isSignificant)) + theme_light() +
      geom_text(aes(label =ifelse(p.adjust.value <= 0.05, RepeatName, "")), hjust = -0.2, vjust = -0.2) +
      coord_cartesian(xlim = c(0, max(df[!is.na(df$p.adjust.value), "obsOnTrueMean"]) * 1.1),
                      ylim = c(0, max(df[!is.na(df$p.adjust.value), "log.p"]) * 1.1))

    ggplot(df[top15.df,],
           aes(x = obsOnTrueMean, y = log.p)) +
      geom_point(aes(color = isSignificant)) + theme_light() +
      geom_text(aes(label =ifelse(p.adjust.value <= 0.05, RepeatName, "")),
                hjust = -0.2, vjust = -0.2) +
      coord_cartesian(xlim = c(0, max(df[!is.na(df$p.adjust.value), "obsOnTrueMean"]) * 1.1),
                      ylim = c(0, max(df[!is.na(df$p.adjust.value), "log.p"]) * 1.1))

    ggsave(paste0(outDir,"/Fig1.pdf"), width = 4, height = 4)

  }else if(type == "RepeatType"){
    ggplot(df[!is.na(df$p.adjust.value),],
           aes(x = obsOnTrueMean, y = log.p)) +
      geom_point(aes(color = isSignificant)) + theme_light() +
      geom_text(aes(label =ifelse(p.adjust.value <= 0.05, RepeatType, "")), hjust = -0.2, vjust = -0.2) +
      coord_cartesian(xlim = c(0, max(df[!is.na(df$p.adjust.value), "obsOnTrueMean"]) * 1.1),
                      ylim = c(0, max(df[!is.na(df$p.adjust.value), "log.p"]) * 1.1))

    ggplot(df[top15.df,],
           aes(x = obsOnTrueMean, y = log.p)) +
      geom_point(aes(color = isSignificant)) + theme_light() +
      geom_text(aes(label =ifelse(p.adjust.value <= 0.05, RepeatType, "")),
                hjust = -0.2, vjust = -0.2) +
      coord_cartesian(xlim = c(0, max(df[!is.na(df$p.adjust.value), "obsOnTrueMean"]) * 1.1),
                      ylim = c(0, max(df[!is.na(df$p.adjust.value), "log.p"]) * 1.1))

    ggsave(paste0(outDir,"/Fig1.pdf"), width = 4, height = 4)

  }else if(type == "RepeatFamily"){
    ggplot(df[!is.na(df$p.adjust.value),],
           aes(x = obsOnTrueMean, y = log.p)) +
      geom_point(aes(color = isSignificant)) + theme_light() +
      geom_text(aes(label =ifelse(p.adjust.value <= 0.05, RepeatFamily, "")), hjust = -0.2, vjust = -0.2) +
      coord_cartesian(xlim = c(0, max(df[!is.na(df$p.adjust.value), "obsOnTrueMean"]) * 1.1),
                      ylim = c(0, max(df[!is.na(df$p.adjust.value), "log.p"]) * 1.1))

    ggplot(df[top15.df,],
           aes(x = obsOnTrueMean, y = log.p)) +
      geom_point(aes(color = isSignificant)) + theme_light() +
      geom_text(aes(label =ifelse(p.adjust.value <= 0.05, RepeatFamily, "")),
                hjust = -0.2, vjust = -0.2) +
      coord_cartesian(xlim = c(0, max(df[!is.na(df$p.adjust.value), "obsOnTrueMean"]) * 1.1),
                      ylim = c(0, max(df[!is.na(df$p.adjust.value), "log.p"]) * 1.1))

    ggsave(paste0(outDir,"/Fig1.pdf"), width = 4, height = 4)

  }


}


#### Fig 2 ####
MakeFigure2 <- function(df, outDir, type){

  order(df$log.p, decreasing = T)[1:15] -> top15.df

  conversion.coef.fig2 <-
    max(df[top15.df, ]$observed) / (max((
      df[top15.df, ]$observed / df[top15.df, ]$rmsk
    ) * 100))

  if(type == "RepeatName"){
    ggplot(df[top15.df, ], aes(x = RepeatName)) +
      geom_line(data = df[top15.df, ],
                mapping = aes(
                  x = RepeatName,
                  y = (observed / rmsk) * 100 * conversion.coef.fig2,group = 1),inherit.aes = F)+
      geom_col(aes(y = observed, alpha = 0.5), show.legend = F) +
      theme_light() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
      scale_y_continuous(
        name = "Number of Observed PARs",
        sec.axis = sec_axis(~ . / conversion.coef.fig2, name = "Percentage of Observed PARs in Repeat Category")
      )


    ggsave(paste0(outDir,"/Fig2.pdf"), width = 4, height = 4)

  }else if(type == "RepeatType"){
    ggplot(df[top15.df, ], aes(x = RepeatType)) +
      geom_line(data = df[top15.df, ],
                mapping = aes(
                  x = RepeatType,
                  y = (observed / rmsk) * 100 * conversion.coef.fig2,group = 1),inherit.aes = F)+
      geom_col(aes(y = observed, alpha = 0.5), show.legend = F) +
      theme_light() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
      scale_y_continuous(
        name = "Number of Observed PARs",
        sec.axis = sec_axis(~ . / conversion.coef.fig2, name = "Percentage of Observed PARs in Repeat Category")
      )


    ggsave(paste0(outDir,"/Fig2.pdf"), width = 4, height = 4)

  }else if(type == "RepeatFamily"){
    ggplot(df[top15.df, ], aes(x = RepeatFamily)) +
      geom_line(data = df[top15.df, ],
                mapping = aes(
                  x = RepeatFamily,
                  y = (observed / rmsk) * 100 * conversion.coef.fig2,group = 1),inherit.aes = F)+
      geom_col(aes(y = observed, alpha = 0.5), show.legend = F) +
      theme_light() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
      scale_y_continuous(
        name = "Number of Observed PARs",
        sec.axis = sec_axis(~ . / conversion.coef.fig2, name = "Percentage of Observed PARs in Repeat Category")
      )


    ggsave(paste0(outDir,"/Fig2.pdf"), width = 4, height = 4)

  }

}
#### Fig 3 ####

MakeFigure3 <- function(df, outDir, type){


  order(df$log.p, decreasing = T)[1:15] -> top15.df

  conversion.coef.fig3 <-
    max(df[top15.df, ]$observed) / max(df[top15.df, ]$log.p)

  if(type == "RepeatName"){
    ggplot(df[top15.df, ], aes(x = RepeatName)) + geom_line(
      data = df[top15.df, ],
      mapping = aes(
        x = RepeatName,
        y = log.p * conversion.coef.fig3,
        group = 1
      ),
      inherit.aes = F
    ) + geom_col(aes(y = observed, alpha = 0.5), show.legend = F) + theme_light() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + scale_y_continuous(
      name = "Number of Observed PARs",
      sec.axis = sec_axis(~ . / conversion.coef.fig3, name = "-log10(Adjusted P Value)")
    )

    ggsave(paste0(outDir,"/Fig3.pdf"), width = 4, height = 4)

  }else if(type == "RepeatType"){
    ggplot(df[top15.df, ], aes(x = RepeatType)) + geom_line(
      data = df[top15.df, ],
      mapping = aes(
        x = RepeatType,
        y = log.p * conversion.coef.fig3,
        group = 1
      ),
      inherit.aes = F
    ) + geom_col(aes(y = observed, alpha = 0.5), show.legend = F) + theme_light() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + scale_y_continuous(
      name = "Number of Observed PARs",
      sec.axis = sec_axis(~ . / conversion.coef.fig3, name = "-log10(Adjusted P Value)")
    )

    ggsave(paste0(outDir,"/Fig3.pdf"), width = 4, height = 4)

  }else if(type == "RepeatFamily"){
    ggplot(df[top15.df, ], aes(x = RepeatFamily)) + geom_line(
      data = df[top15.df, ],
      mapping = aes(
        x = RepeatFamily,
        y = log.p * conversion.coef.fig3,
        group = 1
      ),
      inherit.aes = F
    ) + geom_col(aes(y = observed, alpha = 0.5), show.legend = F) + theme_light() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + scale_y_continuous(
      name = "Number of Observed PARs",
      sec.axis = sec_axis(~ . / conversion.coef.fig3, name = "-log10(Adjusted P Value)")
    )

    ggsave(paste0(outDir,"/Fig3.pdf"), width = 4, height = 4)

  }

}

#### Fig 4 ####
MakeFigure4 <- function(df, outDir, type){

  order(df$log.p, decreasing = T)[1:15] -> top15.df

  conversion.coef.fig4 <-
    max(df[top15.df, ]$obsOnTrueMean) / max(df[top15.df, ]$log.p)

  if(type == "RepeatName"){
    ggplot(df[top15.df, ], aes(x = RepeatName)) + geom_line(
      data = df[top15.df, ],
      mapping = aes(
        x = RepeatName,
        y = log.p * conversion.coef.fig4,
        group = 1
      ),
      inherit.aes = F
    ) + geom_col(aes(y = obsOnTrueMean, alpha = 0.5), show.legend = F) + theme_light() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + scale_y_continuous(
      name = "PAR Fold Enrichment",
      sec.axis = sec_axis(~ . / conversion.coef.fig4, name = "-log10(Adjusted P Value)")
    )

    ggsave(paste0(outDir,"/Fig4.pdf"), width = 4, height = 4)

  }else if(type == "RepeatType"){
    ggplot(df[top15.df, ], aes(x = RepeatType)) + geom_line(
      data = df[top15.df, ],
      mapping = aes(
        x = RepeatType,
        y = log.p * conversion.coef.fig4,
        group = 1
      ),
      inherit.aes = F
    ) + geom_col(aes(y = obsOnTrueMean, alpha = 0.5), show.legend = F) + theme_light() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + scale_y_continuous(
      name = "PAR Fold Enrichment",
      sec.axis = sec_axis(~ . / conversion.coef.fig4, name = "-log10(Adjusted P Value)")
    )

    ggsave(paste0(outDir,"/Fig4.pdf"), width = 4, height = 4)

  }else if(type == "RepeatFamily"){
    ggplot(df[top15.df, ], aes(x = RepeatFamily)) + geom_line(
      data = df[top15.df, ],
      mapping = aes(
        x = RepeatFamily,
        y = log.p * conversion.coef.fig4,
        group = 1
      ),
      inherit.aes = F
    ) + geom_col(aes(y = obsOnTrueMean, alpha = 0.5), show.legend = F) + theme_light() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + scale_y_continuous(
      name = "PAR Fold Enrichment",
      sec.axis = sec_axis(~ . / conversion.coef.fig4, name = "-log10(Adjusted P Value)")
    )

    ggsave(paste0(outDir,"/Fig4.pdf"), width = 4, height = 4)

  }

}
