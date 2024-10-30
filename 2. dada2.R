getwd()
setwd("/home/ke36dar/He/one_manu")

library(dada2)
packageVersion("dada2") #1.26
library(ShortRead)
packageVersion("ShortRead") #1.56.1
library(Biostrings)
packageVersion("Biostrings") #2.66.0
library(dplyr)
packageVersion("dplyr") #1.1.4
library(stringr)
packageVersion("stringr") #1.5.1

##################one_run
path <- "/home/ke36dar/He/one_manu"
list.files(path)

##Forward and reverse fastq
fnFs <- list.files(path, pattern="_R1.fastq.gz", full.names =TRUE)
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz", full.names =TRUE))
sample.names <- sapply(strsplit(basename(fnRs), "_"), '[', 1)
sample.names <- str_replace(sample.names,"341F-785R-16S-bac-", "")
sample.names <- str_replace(sample.names,"341F-785R-", "")
sample.names <- str_replace(sample.names,"-bac", "")
sample.names <- str_replace(sample.names,"PNK87-02", "-87")
sample.names <- str_replace(sample.names,"PNK84-02", "-84")
sample.names <- str_replace(sample.names,"-16S", "")
sample.names <- str_replace(sample.names,"A", "")
sample.names <- str_replace(sample.names,"PNK-", "")
sample.names <- str_replace(sample.names,"-0-2u", "")
sample.names <- str_replace(sample.names,"bac","")
sample.names <- str_replace(sample.names,"PNK","-")
sample.names <- str_replace(sample.names,"--","-")
sample.names <- str_replace(sample.names,"81-02", "81")
sample.names <- str_replace(sample.names,"78-02", "78")
sample.names <- str_replace(sample.names,"8-50", "8")
sample.names <- str_replace(sample.names,"7-50", "7")
sample.names <- str_replace(sample.names,".trim", "")
sample.names <- str_replace(sample.names,"-", ".PNK")
sample.names <- str_replace(sample.names,"S1.PNK", "S1.JSC")
sample.names <- str_replace(sample.names,"S2.PNK", "S2.JSC")
sample.names <- str_replace(sample.names,"JSCJSC", "JSC")

## inspection of read quality profiles
pdf("QualityProfilefnFs.pdf")
plotQualityProfile(fnFs[996:1008])
plotQualityProfile(fnFs[1009:1018])
#plotQualityProfile(fnFs[25:36])
#plotQualityProfile(fnFs[37:48])
#plotQualityProfile(fnFs[49])
#plotQualityProfile(fnFs[61:72])
#plotQualityProfile(fnFs[73:84])
dev.off()

pdf("QualityProfilefnRs.pdf")
plotQualityProfile(fnRs[996:1008])
plotQualityProfile(fnRs[1009:1018])
#plotQualityProfile(fnRs[25:36])
#plotQualityProfile(fnRs[37:48])
#plotQualityProfile(fnRs[49])
#plotQualityProfile(fnRs[61:72])
#plotQualityProfile(fnRs[73:84])
dev.off()


###filtering and trimming
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(265, 235),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)


###error rates

set.seed(100)

errF <- learnErrors(filtFs,nbases = 1e8, multithread=TRUE)
errR <- learnErrors(filtRs,nbases = 1e8, multithread=TRUE)


pdf("errFplot.pdf")
plotErrors(errF, nominalQ=TRUE)
dev.off()

pdf("errRplot.pdf")
plotErrors(errR, nominalQ=TRUE)
dev.off()

##Dereplication of identical reads
derepFs <- derepFastq(filtFs, verbose =TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

##sample inference

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]

###merging paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

head(mergers[[1]])

##constraction of a sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

saveRDS(seqtab, "seqtab.rds")

getN <- function(x)sum(getUniques(x))
track <- cbind(out, sapply(mergers, getN))
colnames(track)<-c("input", "filtered", "merged")
rownames(track) <- sample.names
head(track)
write.csv(track, "read_numbers.csv")
