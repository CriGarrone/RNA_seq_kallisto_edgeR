# RNA_seq_kallisto_edgeR

## Cleanup files
We will copy a renamed version of the abundance.tsv files from each experiment to a parent directory.
```{bash}
cd ~/kallisto_anya_rnaseq
mkdir consolidated_tsv
for folder in *_bam;do
  cp $folder/abundance.tsv consolidated_tsv/$folder.tsv
done

# mv is to move files, cp is to copy
```

## tximport to get count matrices / offsets 

 
helpful material: https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html


```{r}

BiocManager::install("tximport")
BiocManager::install("ensembldb")
BiocManager::install("AnnotationHub")

library(tximport)
library(ensembldb)
library(AnnotationHub)
library(AnnotationDbi)
library(tidyverse)

ah <- AnnotationHub()
ahDb <- query(ah, pattern = c("Mus Musculus", "EnsDb", 87))
ahDb

ahEdb <- ahDb[[1]] #retrieve the MusMusculus Db
ahEdb

txdb <- transcripts(ahEdb, return.type = "data.frame")
txdb

tx2gene <- select(txdb, tx_name, gene_id) 
#check order of tx_name and gene_id..
head(tx2gene)

files <- c("~/kallisto_anya_rnaseq/consolidated_tsv/mESC_1_bam.tsv","~/kallisto_anya_rnaseq/consolidated_tsv/mESC_2_bam.tsv", "~/kallisto_anya_rnaseq/consolidated_tsv/diff3_1_bam.tsv", "~/kallisto_anya_rnaseq/consolidated_tsv/diff3_2_bam.tsv", "~/kallisto_anya_rnaseq/consolidated_tsv/diff6_1_bam.tsv","~/kallisto_anya_rnaseq/consolidated_tsv/diff6_2_bam.tsv", "~/kallisto_anya_rnaseq/consolidated_tsv/diff12_1_bam.tsv", "~/kallisto_anya_rnaseq/consolidated_tsv/diff12_2_bam.tsv", "~/kallisto_anya_rnaseq/consolidated_tsv/diff24_1_bam.tsv", "~/kallisto_anya_rnaseq/consolidated_tsv/diff24_2_bam.tsv")

txi.kallisto.tsv <- tximport(files ,type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion = TRUE) #do i have to set the txOut = TRUE????
head(txi.kallisto.tsv$counts)

cts <- txi.kallisto.tsv$counts
normMat <- txi.kallisto.tsv$length
normMat <- normMat/exp(rowMeans(log(normMat)))

library(edgeR)
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
y <- DGEList(cts)
y$offset <- t(t(log(normMat)) + o)

# y should be ready for estimate dispersion functions 

```

## import countfiles into edgeR (no longer needed I guess) 
Several things to be done  
1) Create a lists of files to be imported
2) Import these files into R using edgeR::readDGE()
```{r}
# install edgeR if not yet installed
if (!requireNamespace("BiocManager", quietly = TRUE)) 
    install.packages("BiocManager")

BiocManager::install("edgeR")
library(edgeR)

# create a list of files
kallisto.files <- list.files("~/kallisto_anya_rnaseq/consolidated_tsv/", pattern = "*.tsv")
kallisto.files
kallisto.edgeR <- readDGE(kallisto.files, path = "~/kallisto_anya_rnaseq/consolidated_tsv/", columns=c(1,5))
View(kallisto.edgeR$samples)
View(kallisto.edgeR$counts)
```

## Filtering 

Link to tutorial = https://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day3/rnaSeq_DE.pdf

Here, we retain only those genes that are represented at least 1cpm reads in at least two samples (cpm=counts per million).

```{r}
View(kallisto.edgeR$counts)
countsPerMillion <- cpm(kallisto.edgeR$counts)
head(countsPerMillion)

countCheck <- countsPerMillion > 1
head(countCheck)

?rowSums
?which

keep <- which(rowSums(countCheck) >= 2)
dgList <- kallisto.edgeR[keep,]
summary(cpm(dgList))
View(dgList$counts)
View(kallisto.edgeR$counts)

```

## Normalisation (No longer needed??)

edgeR implements the trimmed mean of M-values (TMM) method.

```{r}
dgList <- calcNormFactors(dgList, method="TMM")
dgList$samples

eff.lib.size <- dgList$samples$lib.size*dgList$samples$norm.factors

NormalisedCounts <- cpm(dgList)
NormalisedCounts
pseudoNormalisedCounts <- log2(NormalisedCounts + 1)
boxplot(pseudoNormalisedCounts, col="gray", las=3)
plotMDS(pseudoNormalisedCounts)
```
