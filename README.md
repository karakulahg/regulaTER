# ToolX: An R package for ..
This repo is currently under review. Citation information will be provided as soon as our work is accepted. 
### What is this package used for? 


#### What are the dependencies for ToolX ?
1. [R](https://www.r-project.org/) version should be version 3.5+
2. While using r programming, we suggest you to use [Rstudio](https://www.rstudio.com/products/rstudio/download/) which is the R statistical computing environment to use and understand functions TEffectR well.
3. [Bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html) is required on your local computer.
4. [devtools](https://cran.r-project.org/web/packages/devtools/readme/README.html) is required to install TEffectR.
5. TEffectR uses these R packages so you have to install all of them. You may visit the following websites to install them easily: 
    - [dplyr](https://dplyr.tidyverse.org/)
    - [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
    - [biomartr](https://cran.r-project.org/web/packages/biomartr/readme/README.html)


### How to install this R package ?
```

library(devtools)

devtools::install_github("karakulahg/ToolX")

```

### How does it work?

1. Load the library:
```

library(biomartr)
library(GenomicRanges)
library(dplyr)
library(valr)
library(ToolX)


```

2. 
```

rmsk<-biomartr::read_rm("../test/hg19.fa.out") # repeat annotation
input.file <- read.csv("../test/Annotated-Sorted_09_018U_004MSoton_MCF7-ZEB1_unind_Pol2_hs_i91_peaks.narrowPeak.tsv", header=TRUE, stringsAsFactors=FALSE, sep = "\t")

```
3. 
```

pathList <- list("Promoter" = "/home/nazmiye/Desktop/BIP/test/hg19.Promoter",
     "Exon" = "../test/hg19.Exons",
     "Intron" = "../test/hg19.Introns",
     "5UTR" = "../test/hg19.5Prime",
     "3UTR" = "../test/hg19.3Prime",
     "Downstream" = "../test/hg19.Downstream",
     "genomeSizePath" = "../test/hg19.chrom.sizes"
  )


```
4.
```

test <- EnrichPARs(inputPeakFile = input.file, pathList = pathList, numberOfShuffle = 2, repeatMaskerFile = rmsk)

FindMotifs(df = test, repeatMaskerFile = rmsk, outDir = "../test/", homerPath = "~/Downloads/Tools/homer/")

```

#### Session Info

```
> sessionInfo()
R version 3.6.0 (2019-04-26)
Platform: x86_64-redhat-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /usr/lib64/R/lib/libRblas.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
[1] compiler_3.6.0 tools_3.6.0   

```
