
group <- factor( c(1, 1, 2, 2) )
model.matrix(~ group)

group <- c(1, 1, 2, 2)
model.matrix(~ group)

group <- factor(c(1,1,2,2,3,3))
model.matrix(~ group)

group <- factor(c(1,1,2,2,3,3))
group <- relevel(x=group, ref=3)
model.matrix(~ group)

diet <- factor(c(1,1,1,1,2,2,2,2))
sex <- factor(c("f","f","m","m","f","f","m","m"))
model.matrix(~ diet + sex)

model.matrix(~ diet + sex + diet:sex)

options(repr.plot.height=5, repr.plot.width=7)
plot(x=0:40, y=dnbinom(0:40, size=10, prob=0.5), 
     type="b", lwd=2, ylim=c(0, 0.15),
     xlab="Counts (k)", ylab="Probability density")
lines(x=0:40, y=dnbinom(0:40, size=20, prob=0.5), 
      type="b", lwd=2, lty=2, pch=2)
lines(x=0:40, y=dnbinom(0:40, size=10, prob=0.3),
      type="b", lwd=2, lty=3, pch=3)
lines(x=0:40, y=dpois(0:40, lambda=9), col="red")
lines(x=0:40, y=dpois(0:40, lambda=20), col="red")
legend("topright", lwd=c(2,2,2,1), lty=c(1:3,1), pch=c(1:3,-1), col=c(rep("black", 3), "red"),
       legend=c("n=10, p=0.5", "n=20, p=0.5", "n=10, p=0.3", "Poisson"))

options(repr.plot.height=3, repr.plot.width=3)
set.seed(1)
onesample=rpois(100, lambda=2)
xdens = seq(min(onesample), max(onesample), by=1)
ydens = length(onesample) * dpois(xdens, lambda=2)
par(mar=c(2, 4, 0.1, 0.1))
options(repr.plot.width=4, repr.plot.height=3)
res=hist(onesample, main="", prob=FALSE, col="lightgrey", xlab="", ylim=c(0, max(ydens)),
     breaks=seq(-0.5, round(max(onesample))+0.5, by=0.5))
lines(xdens, ydens, lw=2)

set.seed(1)
samplingdistr=replicate(1000, mean(rpois(100, lambda=2)))
par(mar=c(4, 4, 0.1, 0.1))
res=hist(samplingdistr, xlab="Means of the 100-counts", main="")
dy=density(samplingdistr, adjust=2)
## Note this next step is just an empirical way to put the density line approximately the same scale as the histogram:
dy$y = dy$y / max(dy$y) * max(res$counts)
lines(dy)

options(repr.plot.width=6, repr.plot.height=3)
par(mfrow=c(1,2))
hist(onesample, main="One sample of 100 counts\n from Poisson(lambda=2)", xlab="counts")
abline(v=2, col="red", lw=2)
res=hist(samplingdistr, xlab="Means of the 100-counts", main="1000 samples\n of n=100")
abline(v=2, col="red", lw=2)

options(repr.plot.width=6, repr.plot.height=3)
par(mfrow=c(1,3))
hist(replicate(1000, mean(rpois(30, lambda=2))), xlim=c(1.3, 3),
     xlab="Means of 30-counts", main="1000 samples\n of n=30")
abline(v=2, col="red", lw=2)
hist(replicate(1000, mean(rpois(100, lambda=2))), xlim=c(1.3, 3),
     xlab="Means of 100-counts", main="1000 samples\n of n=100")
abline(v=2, col="red", lw=2)
hist(replicate(1000, mean(rpois(500, lambda=2))), xlim=c(1.3, 3),
     xlab="Means of 500-counts", main="1000 samples\n of n=500")
abline(v=2, col="red", lw=2)

options(repr.plot.width=6, repr.plot.height=3)
par(mfrow=c(1,3))
set.seed(5)
onesample=rlnorm(n=100, sdlog=1.5)
hist(onesample, main="100 observations\n from lognormal pop'n", xlab="counts")
sampling30=replicate(1000, mean(rlnorm(n=30, sdlog=1.5)))
hist(sampling30, xlab="Means of the 100 obs.", main="1000 samples\n of n=30")
sampling1000=replicate(1000, mean(rlnorm(n=1000, sdlog=1.5)))
hist(sampling1000, xlab="Means of the 1000 obs.", main="1000 samples\n of n=1000")

suppressPackageStartupMessages({
    library(DESeq2)
    library(airway)
})
data(airway)
dds <- DESeqDataSet(airway, design = ~ cell + dex)
dds <- DESeq(dds)

res <- results(dds)
res

res <- results(dds, contrast=c("dex","trt","untrt"))

mcols(res, use.names = TRUE)

summary(res)

res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)

resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)

results(dds, contrast = c("cell", "N061011", "N61311"))

sum(res$pvalue < 0.05, na.rm=TRUE)

sum(!is.na(res$pvalue))

sum(res$padj < 0.1, na.rm=TRUE)

## [1] 4816

resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])

head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])

res <- lfcShrink(dds, contrast=c("dex","trt","untrt"), res=res)
plotMA(res, ylim = c(-5, 5))

res.noshr <- results(dds)
plotMA(res.noshr, ylim = c(-5, 5))

plotMA(res, ylim = c(-5,5))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

volcanoPlot <- function(res, lfc=2, pval=0.01){
    tab = data.frame(logFC = res$log2FoldChange, negLogPval = -log10(res$pvalue))
    plot(tab, pch = 16, cex = 0.6, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue))
    signGenes = (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))
    points(tab[signGenes, ], pch = 16, cex = 0.8, col = "red") 
    abline(h = -log10(pval), col = "green3", lty = 2) 
    abline(v = c(-lfc, lfc), col = "blue", lty = 2) 
    mtext(paste("pval =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1) 
    mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)
}
volcanoPlot(res)

topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("dex"))

library("ggbeeswarm")
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("dex","cell"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = dex, y = count, color = cell)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)

ggplot(geneCounts, aes(x = dex, y = count, color = cell, group = cell)) +
  scale_y_log10() + geom_point(size = 3) + geom_line()

rld <- rlog(dds, blind=FALSE, fitType="mean")
suppressPackageStartupMessages(library(genefilter))
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 20)

library(pheatmap)
mat  <- assay(rld)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[, c("cell","dex")])
pheatmap(mat, annotation_col = anno)

resGR <- results(dds, lfcThreshold = 1, format = "GRanges")

resGR

library(org.Hs.eg.db)
resGR$symbol <- mapIds(org.Hs.eg.db, names(resGR), "SYMBOL", "ENSEMBL")

suppressPackageStartupMessages(library(Gviz))

window <- resGR[topGene] + 1e6
strand(window) <- "*"
resGRsub <- resGR[resGR %over% window]
naOrDup <- is.na(resGRsub$symbol) | duplicated(resGRsub$symbol)
resGRsub$group <- ifelse(naOrDup, names(resGRsub), resGRsub$symbol)

status <- factor(ifelse(resGRsub$padj < 0.1 & !is.na(resGRsub$padj),
                     "sig", "notsig"))

options(repr.plot.height=5)
options(ucscChromosomeNames = FALSE)
g <- GenomeAxisTrack()
a <- AnnotationTrack(resGRsub, name = "gene ranges", feature = status)
d <- DataTrack(resGRsub, data = "log2FoldChange", baseline = 0,
               type = "h", name = "log2 fold change", strand = "+")
plotTracks(list(g, d, a), groupAnnotation = "group",
           notsig = "grey", sig = "hotpink")

suppressPackageStartupMessages(library(sva))

dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ dex, colData(dds))
mod0 <- model.matrix(~   1, colData(dds))
svseq <- svaseq(dat, mod, mod0, n.sv = 2)

svseq$sv

options(repr.plot.height=5)
par(mfrow = c(2, 1), mar = c(3,5,3,1))
for (i in 1:2) {
  stripchart(svseq$sv[, i] ~ dds$cell, vertical = TRUE, main = paste0("SV", i))
  abline(h = 0)
 }

ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
design(ddssva) <- ~ SV1 + SV2 + dex

DESeq(ddssva)
