
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(monocle))

load(file="~/workshop_materials/scrna_workshop_data/pbmc_3k_sce_filt.rda")

sce.filt <- computeSumFactors(sce.filt)
summary(sizeFactors(sce.filt))

plot(sizeFactors(sce.filt), sce.filt$total_counts/1e6, log="xy",ylab="Library size (millions)", xlab="Size factor")

tail(exprs(sce.filt))

sce.filt <- normalize(sce.filt,exprs_values="counts")

tail(exprs(sce.filt))

var.fit.nospike <- trendVar(sce.filt, parametric=TRUE, use.spikes=FALSE, span=0.2)

var.out.nospike <- decomposeVar(sce.filt, var.fit.nospike)
head(var.out.nospike)

ggplot(var.out.nospike, aes(x=tech, y=bio, color=p.value))+geom_point()

plot(var.out.nospike$mean, var.out.nospike$total, pch=16, cex=0.6, 
    xlab="Mean log-expression", ylab="Variance of log-expression")
curve(var.fit.nospike$trend(x), col="dodgerblue", lwd=2, add=TRUE)
#points(var.out.nospike$mean[cur.spike], var.out.nospike$total[cur.spike], col="red", pch=16)

names(colData(sce.filt))

design <- model.matrix(~cycle_phases,colData(sce.filt))
head(design)

system.time(var.fit.nospike.dsgn <- trendVar(sce.filt, parametric=TRUE, span=0.4, design=design))

system.time(var.out.nospike.dsgn <- decomposeVar(sce.filt, var.fit.nospike.dsgn))
head(var.out.nospike.dsgn)

#g1 <-ggplot(var.out.nospike, aes(x=tech, y=bio, color=FDR))+geom_point()

#var.out.nospike.sig_bio <- rep("yes", nrow(var.out.nospike))
#var.out.nospike.sig_bio[var.out.nospike$FDR > 0.005] <-"NO"
#g2 <- ggplot(var.out.nospike, aes(x=tech, y=bio, color=var.out.nospike.sig_bio))+geom_point()
#multiplot(g1,g2,cols = 1)

plot(var.out.nospike$bio, var.out.nospike.dsgn$bio, pch=16, cex=0.6, 
    xlab="biological variance(no regression)", ylab="biological variance (cell-cycle)")
abline(0,1)

fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))

chosen.genes <- order(var.out.nospike$bio, decreasing=TRUE)[1:20]
plotExpression(sce.filt, rownames(var.out.nospike)[chosen.genes], alpha=0.05, jitter="jitter") + fontsize

chosen.genes <- order(var.out.nospike.dsgn$bio, decreasing=TRUE)[1:20]
plotExpression(sce.filt, rownames(var.out.nospike.dsgn)[chosen.genes], alpha=0.05, jitter="jitter") + fontsize

system.time(sce2 <- denoisePCA(sce.filt, technical=var.fit.nospike.dsgn$trend) )
dim(reducedDim(sce2, "PCA")) 

sce2_lr <- denoisePCA(sce2, technical=var.fit.nospike$trend, value="lowrank") 
assayNames(sce2_lr)

IRdisplay::display_html('<iframe src="https://distill.pub/2016/misread-tsne/" width=1000, height=500> </iframe>')

plotReducedDim(sce2, use_dimred="PCA", ncomponents=3, colour_by="FTL")

    out5 <- plotTSNE(sce2, use_dimred="PCA", perplexity=5, colour_by="FTL", 
                     rand_seed=100) + fontsize + ggtitle("Perplexity = 50")

out5

out10 <- plotTSNE(sce2, use_dimred="PCA", perplexity=10, colour_by="total_features",
    rand_seed=100) + fontsize + ggtitle("Perplexity = 10")
out10

out20 <- plotTSNE(sce2, use_dimred="PCA", perplexity=20, colour_by="LST1",
    rand_seed=100) + fontsize + ggtitle("Perplexity = 20")
out20

plotDiffusionMap(sce2, use_dimred="PCA", sigma=25, colour_by="cycle_phases") + fontsize

pbmc <- CreateSeuratObject(raw.data = exprs(sce.filt),project = "10X_PBMC")

pbmc@meta.data <- cbind(pbmc@meta.data,as.data.frame(colData(sce.filt)))

pbmc@calc.params$NormalizeData <-list(assay.type="RNA",normalization.method="LogNormalize",scale.factor=NULL)

head(pbmc@meta.data)

pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR,
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = pbmc@var.genes)

pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "pct_counts_MT","cycle_phases"))

pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

PrintPCA(object = pbmc, pcs.print = 1:5, genes.print = 5, use.full = FALSE)

VizPCA(object = pbmc, pcs.use = 1:2)

PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2)

PCHeatmap(object = pbmc, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

PCHeatmap(object = pbmc, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)

PCElbowPlot(object = pbmc)

pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10, 
                     resolution = 0.6, print.output = 0, save.SNN = TRUE)

pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE)

# note that you can set do.label=T to help label individual clusters
TSNEPlot(object = pbmc, dark.theme=T)

#save(pbmc,file="~/workshop_materials/scrna_workshop_data/pbmc_from_sce.rds")

current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)

new.cluster.ids <- c("CD4 T cells", "CD14+ Monocytes",
                     "B cells", "CD8 T cells", "FCGR3A+ Monocytes",
                     "NK cells", "Dendritic cells", "Megakaryocytes")
pbmc@ident <- plyr::mapvalues(x = pbmc@ident, from = current.cluster.ids, to = new.cluster.ids)
TSNEPlot(object = pbmc, do.label = TRUE, pt.size = 0.5, dark.theme=T)
