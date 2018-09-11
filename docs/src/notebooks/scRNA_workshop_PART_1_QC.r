
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(monocle))

library(IRdisplay)

suppressPackageStartupMessages(library(scran))

suppressPackageStartupMessages(library(scater))

pbmc.data <- Read10X(data.dir = "/home/ubuntu/workshop_materials/scrna_workshop_data/filtered_gene_bc_matrices/hg19/")

pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 0, min.genes = 0, project = "10X_PBMC")

gene_map <- read.table("/home/ubuntu/workshop_materials/scrna_workshop_data/filtered_gene_bc_matrices/hg19/genes.tsv",
                       header=F)
names(gene_map)<-c("Ens_ID","name")
head(gene_map, n=2)

as.matrix(head(pbmc@raw.data, n=2))

sce <- SingleCellExperiment(list(counts=as.matrix(pbmc@raw.data)))

str(sce)

exprs(sce) <- log2(calculateCPM(sce, use.size.factors = FALSE) + 1)  #SCATER

#assays(sce)

head(rowSums(exprs(sce)), n=2)

head(colnames(sce), n=2)

isSpike(sce, "MT") <- grepl("^MT-", rownames(sce)) #SCATER

sce <- calculateQCMetrics(sce,feature_controls = list(MT = isSpike(sce, "MT") )) #SCATER

cbind(names(colData(sce)))

cbind(names(rowData(sce)))

ens_pos <- match(rownames(counts(sce)), gene_map$name)
rownames(sce@assays[['counts']]) <- gene_map$Ens_ID[ens_pos]

cbind(rownames(counts(sce))[1:6],as.character(gene_map$name[1:6]))

hg.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))

load(file="/home/ubuntu/workshop_materials/scrna_workshop_data/pbmc_3k_cyclone_assigned.rda")

head(assigned$scores)

table(assigned$phases)

newcoldat <- colData(sce)

newcoldat$cycle_phases <- assigned$phases

newcoldat<- cbind(newcoldat,assigned$scores)

colData(sce) <- newcoldat

load(file="/home/ubuntu/workshop_materials/scrna_workshop_data/pbmc_3k_sce_raw.rda")

keep_feature <- rowSums(counts(sce) > 0 ) > 0
table(keep_feature)

sce.filt <- sce[keep_feature,]

dim(sce.filt)

g1 = ggplot(as.data.frame(colData(sce.filt))) 
g1 = g1 + geom_histogram(aes(x=pct_counts_MT),bins=200,fill = "grey")
#g1 = g1 + geom_density(aes(x=pct_counts_MT)) 
g1

mito.keep <- !(isOutlier(sce.filt$pct_counts_MT, nmads=3, type="higher"))
table(mito.keep)

filter_by_MT <- colData(sce.filt)$pct_counts_MT < 10 
table(filter_by_MT)

threshold = median(log(colData(sce.filt)$total_counts)) - 3*mad(log(colData(sce.filt)$total_counts))
threshold

libsize.keep <- !isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
table(libsize.keep)

g1 = ggplot(as.data.frame(colData(sce))) 
g1 = g1 + geom_histogram(aes(x=total_counts, y= ..density..),bins=200,fill = "grey")
g1 = g1 + geom_density(aes(x=total_counts)) 
g1 + geom_vline(xintercept = exp(threshold), color="red")

g1 = ggplot(as.data.frame(colData(sce))) 
g1 = g1 + geom_histogram(aes(x=log(total_counts), y= ..density..),bins=200,fill = "grey")
g1 = g1 + geom_density(aes(x=log(total_counts)))
g1 + geom_vline(xintercept = threshold, color="red")

threshold = median(log(colData(sce)$total_features)) - 3*mad(log(colData(sce)$total_features))
threshold
                   
feature.keep <- !isOutlier(sce$total_features, nmads=3, type="lower", log=TRUE)
table(feature.keep)

g1 = ggplot(as.data.frame(colData(sce))) 
g1 = g1 + geom_histogram(aes(x=log(total_features), y= ..density..),bins=200,fill = "grey")
g1 = g1 + geom_density(aes(x=log(total_features)))
g1 + geom_vline(xintercept = threshold, color="red")

sce.filt$use <- (
    # sufficient features (genes)
    feature.keep &
    # sufficient molecules counted
    libsize.keep &
    # remove cells with unusual number of reads in MT genes
    mito.keep
)

table(sce.filt$use)

filter_genes <- apply(counts(sce.filt[ ,sce.filt$use]), 1, 
    function(x)
        {
        length(x[x > 1]) >= 2
        }
    )

rowData(sce.filt)$use <- filter_genes
table(rowData(sce.filt)$use)

table(colData(sce.filt)$use)

plotQC(sce.filt, type = "highest-expression", exprs_values = "counts")

plotExprsFreqVsMean(sce.filt)

plotMDS(sce.filt,exprs_value="logcounts", colour_by = "cycle_phases")

plotExpression(sce.filt, features=1:20)

colnames(colData(sce.filt))

plotColData(sce.filt, aes(x = log10_total_counts, y = pct_counts_MT, color=cycle_phases))

plotDiffusionMap(sce.filt)
#, color_by="cycle_phases")
#, aes(x = log10_total_counts, y = total_features, colour = log10_mean_counts))

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
pbmc <- AddMetaData(object = pbmc, metadata = colData(sce), col.name = colnames(colData(sce)))

pbmc@meta.data <- pbmc@meta.data[,1:3]
colnames(pbmc@meta.data)

save(sce.filt, file="~/workshop_materials/scrna_workshop_data/pbmc_3k_sce_filt.rda")

rm(pbmc)
rm(pbmc.data)
