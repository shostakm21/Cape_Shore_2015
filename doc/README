# Cape_Shore_2015

```{r}
library(dada2)
library(ggplot2)
library(ape)
library(microbiome)
library(tidyverse)
library(DECIPHER)
library(reshape2)
library(dplyr)
library(ggpattern)
library(vegan)
library('cowplot')
```

```{r}
path <- ("/Users/maggieshostak/Desktop/Cape Shore Amplicon Data/Fastq_Files/")
list.files(path)
```

```{r}
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnFs
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
```

```{r}
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
list(sample.names)
```

```{r}
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
```

```{r}
#FWD <- "GTGYCAGCMGCCGCGGTAA"
#REV <- "CCGYCAATTYMTTTRAGTTT"
#trimLeft = c(FWD, REV)
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
```

```{r}
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(220,200), maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = TRUE, trimLeft = c(19,20))
#head(out)
```

```{r}
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
```

```{r}
#plotErrors(errF, nominalQ=TRUE)
#plotErrors(errR, nominalQ=TRUE)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
#dadaFs[[1]]

dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
#dadaRs[[1]]
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, justConcatenate=TRUE, verbose=TRUE)
head(mergers[[1]])
```

```{r}
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```

```{r}
taxa <- assignTaxonomy(seqtab.nochim, "/Users/maggieshostak/Desktop/Cape Shore Amplicon Data/Fastq_Files/silva_nr99_v138.1_train_set.fa.gz")
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
for (i in 1:dim(seqtab.nochim)[2]) {
asv_headers[i] <- paste(">ASV", i, sep="_")
}
asv_fasta <- c(rbind(asv_headers, asv_seqs))

asv_otu <- t(seqtab.nochim)

row.names(asv_otu) <- sub(">", "", asv_headers)

asv_tax <- taxa

row.names(asv_tax) <- sub(">", "", asv_headers)

otu_tax_table <- merge(asv_otu, asv_tax, by=0)

write(asv_fasta, "asv_fasta_final.fa")
write.table(asv_otu, "asv_otu_final.csv", sep=",", quote=F, col.names=NA)
write.table(asv_tax, "asv_tax_final.csv", sep=",", quote=F, col.names=NA)
write.table(otu_tax_table, "otu_tax_table_final.csv", sep=",", quote=F, col.names=NA)
```

```{r}
otu_counts <- read.csv("/Users/maggieshostak/Desktop/Cape\ Shore\ Amplicon\ Data/Capshore_Project_Dada2/asv_otu_final.csv") %>%
  pivot_longer(-ASV, names_to="sample_id", values_to = "count")
otu_counts

taxonomy <- read.csv("/Users/maggieshostak/Desktop/Cape\ Shore\ Amplicon\ Data/Capshore_Project_Dada2/asv_tax_final.csv")
taxonomy
```

```{r}
otu_rel_abund <- inner_join(taxonomy, otu_counts, by="ASV") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols=c("Kingdom", "Phylum", "Class", "Order", "Family", "ASV"),
         names_to="level",
         values_to="taxon")
otu_rel_abund

write.table(otu_rel_abund, "otu_rel_abund_final.csv", sep=",", quote=F, col.names=NA)

otu_rel_abund <- read.csv("/Users/maggieshostak/Desktop/Cape\ Shore\ Amplicon\ Data/Capshore_Project_Dada2/otu_rel_abund_final.csv")
```

```{r}
## Phylum
otu_rel_abund %>%
  filter(level=="Phylum") %>%
  group_by(sample_id, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(sample_id, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=sample_id, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=sample_id, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("phylum_stacked_barchart_capeshore.tiff", width=20, height=10)

## Class
otu_rel_abund %>%
  filter(level=="Class") %>%
  group_by(sample_id, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(sample_id, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=sample_id, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=sample_id, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("class_stacked_barchart_capeshore.tiff", width=25, height=10)

## Order
otu_rel_abund %>%
  filter(level=="Order") %>%
  group_by(sample_id, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(sample_id, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=sample_id, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=sample_id, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("order_stacked_barchart_capeshore.tiff", width=55, height=25, limitsize = FALSE)

## Family
otu_rel_abund %>%
  filter(level=="Family") %>%
  group_by(sample_id, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(sample_id, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=sample_id, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=sample_id, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("family_stacked_barchart_capeshore.tiff", width=50, height=30, limitsize = FALSE)
```
