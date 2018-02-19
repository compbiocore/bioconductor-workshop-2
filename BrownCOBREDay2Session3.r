
deTable <-
     matrix(c(28, 142, 501, 12000),
            nrow = 2,
            dimnames = list(c("DE", "Not.DE"),
                            c("In.gene.set", "Not.in.gene.set")))
deTable

fisher.test(deTable, alternative = "greater")

suppressPackageStartupMessages(library(EnrichmentBrowser))

library(ALL)
data(ALL)

ind.bs <- grep("^B", ALL$BT)
ind.mut <- which(ALL$mol.biol %in% c("BCR/ABL", "NEG"))
sset <- intersect(ind.bs, ind.mut)
all.eset <- ALL[, sset]

dim(all.eset)

exprs(all.eset)[1:4,1:4]

all.eset <- probe.2.gene.eset(all.eset) 
head(names(all.eset))

library(airway)
data(airway)

air.eset <- airway[grep("^ENSG", names(airway)), ]
dim(air.eset)

assay(air.eset)[1:4,1:4]

all.eset$GROUP <- ifelse(all.eset$mol.biol == "BCR/ABL", 1, 0)
table(all.eset$GROUP)

air.eset$GROUP <- ifelse(colData(airway)$dex == "trt", 1, 0)
table(air.eset$GROUP)

air.eset$BLOCK <- airway$cell
table(air.eset$BLOCK)

all.eset <- de.ana(all.eset)
rowData(all.eset, use.names=TRUE)

air.eset <- de.ana(air.eset, de.method="edgeR")

rowData(air.eset, use.names=TRUE)

kegg.gs <- get.kegg.genesets("hsa")

go.gs <- get.go.genesets(org="hsa", onto="BP", mode="GO.db")

data.dir <- system.file("extdata", package="EnrichmentBrowser")
gmt.file <- file.path(data.dir, "hsa_kegg_gs.gmt")
hsa.gs <- parse.genesets.from.GMT(gmt.file)
length(hsa.gs)

hsa.gs[1:2]

ora.all <- sbea(method="ora", eset=all.eset, gs=hsa.gs, perm=0, alpha=0.2)
gs.ranking(ora.all)

ea.browse(ora.all)

air.eset <- map.ids(air.eset, org="hsa", from="ENSEMBL", to="ENTREZID")

ora.air <- sbea(method="ora", eset=air.eset, gs=hsa.gs, perm=0)
gs.ranking(ora.air)

gsea.all <- sbea(method="gsea", eset=all.eset, gs=hsa.gs, perm=1000)  

gs.ranking(gsea.all)

gsea.air <- sbea(method="gsea", eset=air.eset, gs=hsa.gs, perm=100)  

roast.air <- sbea(method="roast", eset=air.eset, gs=hsa.gs)
gs.ranking(roast.air)  

sbea.methods()

pwys <- file.path(data.dir, "hsa_kegg_pwys.zip")
hsa.grn <- compile.grn.from.kegg(pwys)
head(hsa.grn)

spia.all <- nbea(method="spia", eset=all.eset, gs=hsa.gs, grn=hsa.grn, alpha=0.2)
gs.ranking(spia.all)

ggea.all <- nbea(method="ggea", eset=all.eset, gs=hsa.gs, grn=hsa.grn)
gs.ranking(ggea.all)

nbea.methods()

res.list <- list(ora.all, gsea.all)
comb.res <- comb.ea.results(res.list)
gs.ranking(comb.res)

suppressPackageStartupMessages(library(regioneR))

cpgHMM <- toGRanges("http://www.haowulab.org/software/makeCGI/model-based-cpg-islands-hg19.txt")
cpgHMM <- filterChromosomes(cpgHMM, chr.type="canonical")
cpgHMM <- sort(cpgHMM)
cpgHMM

promoters <- toGRanges("http://gattaca.imppc.org/regioner/data/UCSC.promoters.hg19.bed")
promoters <- filterChromosomes(promoters, chr.type="canonical")
promoters <- sort(promoters)
promoters

cpg <- cpgHMM[seqnames(cpgHMM) %in% c("chr21", "chr22")]
prom <- promoters[seqnames(promoters) %in% c("chr21", "chr22")]

# pt <- overlapPermTest(cpg, prom, genome="hg19", ntimes=100, per.chromosome=TRUE, count.once=TRUE)
download.file("https://www.dropbox.com/s/dlqq0m99wpu5xab/BrownCOBREDay2Session3_pt.rds?raw=1", destfile="BrownCOBREDay2Session3_pt.rds")
pt <- readRDS("BrownCOBREDay2Session3_pt.rds")
pt

summary(pt[[1]]$permuted)
