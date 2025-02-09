library(dada2)
library(dplyr)

setwd("C:/Users/brean/Downloads/masters/Fodor/builtenv/37740013")
fnFs <- sort(list.files("./fastq", pattern = "_1.fastq", full.names = TRUE))

fnRs <- sort(list.files("./fastq", pattern = "_2.fastq", full.names = TRUE))


sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#### Inspect read quality profiles ####
plotQualityProfile(fnFs[1:4])

plotQualityProfile(fnRs[1:4])

#### Filter and trim ####
filtFs <- file.path("./filtered_fastq", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("./filtered_fastq", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE)

head(out)


a <- read.csv("./SraRunTable.txt")
o = out

write.csv(as.data.frame(out), "out.csv", row.names = FALSE)

#### Learn the Error Rates ####
errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)



#### Sample Inference ####

dadaFs <- dada(filtFs, err=errF, multithread=FALSE)

dadaRs <- dada(filtRs, err=errR, multithread=FALSE)

dadaFs[[1]]



#### Merge paired reads ####
mergers <- mergePairs(dadaFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])



#### Construct Sequence table ####
seqtab <- makeSequenceTable(mergers)
dim(seqtab)



# Inspect distribution of sequence lengths ####
table(nchar(getSequences(seqtab)))

#### Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)



#### Track reads through the pipeline ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


