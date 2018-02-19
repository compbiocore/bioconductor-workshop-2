
library(tximportData)
datadir <- system.file("extdata", package = "tximportData")
list.files(datadir)

samples <- read.table(file.path(datadir, "samples.txt"), header=TRUE)
samples

files <- file.path(datadir, "salmon", samples$run, "quant.sf")
names(files) <- paste0("sample",1:6)
all(file.exists(files))

library(tximport)
library(readr)
txi <- tximport(files, type="salmon", txOut=TRUE)
##txi <- tximport(files, type="salmon", tx2gene=tx2gene)
names(txi)
head(txi$counts)

tx2gene <- read.csv(file.path(datadir, "tx2gene.csv"))

suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
## the different units you could summarize by:
listColumns(EnsDb.Hsapiens.v86)
tx3gene <- transcripts(EnsDb.Hsapiens.v86, columns="symbol", 
                  return.type="DataFrame")
(tx3gene <- tx3gene[, 2:1])

library(EnsDb.Hsapiens.v86)

keytypes(EnsDb.Hsapiens.v86)

columns(EnsDb.Hsapiens.v86)

k <- keys(EnsDb.Hsapiens.v86, "TXNAME")
tx4gene <- select(EnsDb.Hsapiens.v86, keys=k, 
                  keytype = "TXNAME", columns = "SYMBOL")
tx4gene <- tx4gene[, 3:2]

# Check that this method is equivalent to the previous one:
all.equal(tx4gene[, 1], tx3gene[, 1])
all.equal(tx4gene[, 2], tx3gene[, 2])

txi.sum <- summarizeToGene(txi, tx2gene)

suppressPackageStartupMessages({
    library(airway)
    library(SummarizedExperiment)
})

data(airway)

airway

dim(airway)

assayNames(airway)

head(assay(airway), 3)

colSums(assay(airway))

rowRanges(airway)

str(metadata(rowRanges(airway)))

colData(airway)

suppressPackageStartupMessages(library(DESeq2))

dds <- DESeqDataSet(airway, design = ~ cell + dex)

nrow(dds)

dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

options(repr.plot.height=3, repr.plot.width=5)
lambda <- 10^seq(from = -1, to = 2, length = 1000)
cts <- matrix(rpois(1000*100, lambda), ncol = 100)
library(vsn)
meanSdPlot(cts, ranks = FALSE)

options(repr.plot.height=3, repr.plot.width=5)
log.cts.one <- log2(cts + 1)
meanSdPlot(log.cts.one, ranks = FALSE)

options(repr.plot.height=3, repr.plot.width=5)
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)
meanSdPlot(assay(rld), ranks=FALSE)

vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)
meanSdPlot(assay(vsd), ranks=FALSE)

dds <- estimateSizeFactors(dds)

ddsESF <- estimateSizeFactors(dds)
df1 <- data.frame(log2(counts(ddsESF, normalized=TRUE)[, 1:2] + 1))
df1$transformation <- "log2(x + 1)"
df2 <- data.frame(assay(rld)[, 1:2])
df2$transformation <- "rld"
df3 <- data.frame(assay(vsd)[, 1:2])
df3$transformation <- "vsd"
df <- rbind(df1, df2, df3)
colnames(df)[1:2] <- c("x", "y")
head(df)

options(repr.plot.width=6, repr.plot.height=2.5)
library(ggplot2)
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  

par(mfrow=c(1,3))
boxplot(log2(assay(ddsESF)+1), las=2, main="log2(x+1)")
boxplot(assay(rld), las=2, main="rld")
boxplot(assay(vsd), las=2, main="vsd")

sampleDists <- dist(t(assay(rld)))
sampleDists

library(pheatmap)
library(RColorBrewer)

options(repr.plot.height=3, repr.plot.width=5)
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$dex, rld$cell, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

library(PoiClaClu)
poisd <- PoissonDistance(t(counts(ddsESF)))

options(repr.plot.height=3, repr.plot.width=5)
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( rld$dex, rld$cell, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

plotPCA(rld, intgroup = c("dex", "cell"))

pcaData <- plotPCA(rld, intgroup = c( "dex", "cell"), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = dex, shape = cell)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()

mds <- cbind(as.data.frame(colData(rld)), cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = dex, shape = cell)) +
  geom_point(size = 3) + coord_fixed()

mdsPois <- cbind(as.data.frame(colData(ddsESF)),
   cmdscale(samplePoisDistMatrix))
ggplot(mdsPois, aes(x = `1`, y = `2`, color = dex, shape = cell)) +
  geom_point(size = 3) + coord_fixed()

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
keytypes(txdb)
k <- keys(txdb, "TXNAME")
columns(txdb)
tx4gene <- select(txdb, keys=k, keytype="TXNAME", column="GENEID")

genes(txdb)
exonsBy(txdb)

suppressPackageStartupMessages(library(AnnotationHub))
ah <- AnnotationHub()
query(ah, "TxDb")
txdb <- ah[["AH52245"]] #TxDb.Athaliana.BioMart.plantsmart22.sqlite

keytypes(txdb)

columns(txdb)

k <- keys(txdb, "TXNAME")
tx5gene <- select(txdb, keys=k, keytype="TXNAME", column="GENEID")
head(tx5gene)

library(GenomicFeatures)
gffFile <- system.file("extdata", "GFF3_files", "a.gff3", 
                       package="GenomicFeatures")
txdb <- makeTxDbFromGFF(file=gffFile,
            dataSource="partial gtf file for Tomatoes for testing",
            organism="Solanum lycopersicum")
keytypes(txdb)

columns(txdb)

k <- keys(txdb, "TXNAME")
tx6gene <- select(txdb, keys=k, keytype="TXNAME", column="GENEID")
head(tx5gene)

grep("makeTxDb", ls("package:GenomicFeatures"), val=TRUE)
