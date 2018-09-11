
suppressPackageStartupMessages({
    library(DiffBind)
    library(BiocFileCache)
})
bfc <- BiocFileCache(cache="~/chipseq-data")
DBfiles <- list.files(system.file("extra", package="DiffBind"), 
                     recursive = TRUE, full.names = TRUE)
DBfiles <- DBfiles[!DBfiles %in% bfcinfo(bfc)$rname]
for (i in seq_along(DBfiles))
  bfcadd(bfc, rname=DBfiles[i], rtype="local", action="copy")

tamfile <- bfcquery(bfc, "tamoxifen.csv")$fpath
read.csv(tamfile)

bfcquery(bfc, "peaks")$rpath

setwd(system.file("extra", package="DiffBind")) #necessary because `tamfile` contains relative paths
ta <- dba(sampleSheet=tamfile)
ta

names(ta)
class(ta$peaks)
head(ta$peaks[[1]])

ta$samples

options(repr.plot.width=5, repr.plot.height=5)
# the next line does not actually work, because the BAM files are not included in the package
# ta <- dba.count(ta, minOverlap=3)
# instead we load the counts:
data(tamoxifen_counts)
plot(tamoxifen)

tamoxifen$config$AnalysisMethod

ta2 <- dba.contrast(tamoxifen, categories=DBA_CONDITION)
ta2 <- dba.analyze(ta2)
ta2

tadb <- dba.report(ta2)
tadb
counts <- dba.report(ta2, bCounts=TRUE)

library(Homo.sapiens)
gn <- genes(Homo.sapiens, columns="SYMBOL")
summary(counts %over% gn)

table(countOverlaps(counts, gn))

count2 <- counts[countOverlaps(counts, gn) == 2]
gn2 <- gn[gn %over% count2]
width(gn2) / 1e3  #width in kb
gn2 <- gn2[order(ranges(gn2))]
gn2

genome(counts) <- "hg19"
library(rtracklayer)
session <- browserSession("UCSC")
genome(session) <- "hg19"
track(session, "counts") <- counts

rangemaxFC <- counts[which.max(abs(counts$Fold))]
browserView(session, range=rangemaxFC * 0.75)  #0.75 zoom-factor

gn2[3:4]

(plotrange <- reduce(gn2[3:4], ignore.strand=TRUE))

browserView(session, range=plotrange * 0.75)

library(ChIPseeker)
peakAnno <-  annotatePeak(counts,
           TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)
peakAnno

options(repr.plot.width=6, repr.plot.height=5)
plotAnnoPie(peakAnno)

suppressPackageStartupMessages(library(rtracklayer))
bl <- import("https://www.encodeproject.org/files/ENCFF419RSJ/@@download/ENCFF419RSJ.bed.gz", genome="hg19")
