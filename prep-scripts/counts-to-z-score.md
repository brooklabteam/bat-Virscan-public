
# Converting peptide counts from PhIP-Seq to Z-scores

Here are directions for transforming peptide counts from the PhIP-Seq data into Z-scores for input into AVARDA.
We first received peptide counts in raw form (see subfolders within the "/raw-data/PhIP-Seq/") from gamma-poisson modeling of sequencing reads (*link to pipeline*). To convert these to Z-scores for input to AVARDA (*link to package*), we ran the script below.

*Emily, note that you need to include all of these .tar.gz zip files in this folder so that someone can exactly repeat your steps with only this repo.*

```{r install packages from source}
install.packages("VGAM")
install.packages("~/Desktop/UChi_Brook_Lab/bat_VirScan_main/Zscore_binning_pipeline/virScanR_0.1.0.9000.tar.gz", repos = NULL, type = "source")

install.packages("/Users/ecornelius/Desktop/UChi_Brook_Lab/bat_VirScan_main/Zscore_binning_pipeline/mmR_0.1.0.tar.gz",repos=NULL,type="source")
```


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r load libraries}
library(dplyr)
library(VGAM)
library(virScanR)
library(mmR)
library(data.table)
```

## Input data format:

"counts" table should have peptides as rows and samples as columns. "id" column represents peptide id's. "input" column is the combined counts for mock-IP wells. propTrim is the proportion of each bin to be removed from the extreme of the distribution. makeGroupSize is the size of each bin. cols_to_not_evaluate excludes columns from being used for z-score calculation.

*Emily, we need to include the raw counts csv files in here so that this is replicable.*

```{r set wd and load count data}
setwd = ("/Users/ecornelius/Desktop/UChi_Brook_Lab/bat_VirScan_main")
```

## Convert the tsv to csv prior to this step. 

You will need to do text to column with ',' separator

```{r load count data}
counts=read.csv("PhIP13_counts.csv")
```


```{r set parameters}
propTrim = 0.05

convert_params <-
  vs.set_params_convert(stat = "Z_score",
                        makeGroups = TRUE,
                        makeGroupSize = 300,
                        makeGroup_by_col = "input",
                        idCol = "id",
                        groupTogetherIfGreaterThanGroupSize = TRUE,
                        splitHighVarMeanGrps = TRUE,
                        cols_to_remove = NULL,
                        cols_to_not_evaluate = c("id","input"), 
                        propExtremesToRemove = propTrim,
                        removeTail = "both",
                        coresAcrossGroups = parallel::detectCores()-2,
                        returnAs = "data.table", 
                        returnParams = TRUE)
```


```{r convert to Z-score}
dat = 
  virScanR::vs.convert(
    data = counts, 
    paramsList = convert_params
  )

z <- dat$out

write.csv(z, file ="PhIP13_z_scores_binning.csv")
```

