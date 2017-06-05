
## cd scratch
## startnode10
## cp -R /depot/nihomics/data/bootcamp2017/unit1/ .
## cd unit1
## jupyter notebook --browser=Firefox

source("http://bioconductor.org/biocLite.R")
biocLite()

biocLite("GEOquery")
biocLite("affy")
biocLite("limma")
biocLite("hgu133plus2.db")
biocLite("Homo.sapiens")

install.packages("hexbin")
install.packages("RColorBrewer")
install.packages("corrplot")

library(GEOquery)

my.gse <- "GSE14987"

if(!file.exists("geo_downloads")) dir.create("geo_downloads")
if(!file.exists("results"))  dir.create("results")

list.files(r=TRUE, all=TRUE)

my.geo.gse <- getGEO(GEO=my.gse, filename=NULL, destdir="geo_downloads", GSElimits=NULL, GSEMatrix=TRUE, 
                     AnnotGPL=FALSE, getGPL=FALSE)


list.files(recursive=TRUE)

ls()

class(my.geo.gse)

length(my.geo.gse)

names(my.geo.gse)

my.geo.gse <- my.geo.gse[[1]]

class(my.geo.gse)

show(my.geo.gse)

head(exprs(my.geo.gse))

head(pData(my.geo.gse))

colnames(pData(my.geo.gse))

pData(my.geo.gse)$data_processing[1]

my.file <- file.path("geo_downloads", my.gse)

my.condition <- file.exists(my.file)

my.condition

if(!my.condition)
getGEOSuppFiles(my.gse, makeDirectory=TRUE, baseDir="geo_downloads")

list.files(recursive=TRUE)

read.delim(file.path("geo_downloads", my.gse, "filelist.txt"), as.is=TRUE)

untar(file.path("geo_downloads", my.gse, paste0(my.gse, "_RAW.tar")), exdir=file.path("geo_downloads", my.gse, "CEL"))
list.files(file.path("geo_downloads", my.gse, "CEL"))

cel.path <- file.path("geo_downloads", my.gse, "CEL/")
cel.path

my.cels <- list.files(cel.path, pattern=".CEL")
my.cels

my.pdata <- as.data.frame(pData(my.geo.gse), stringsAsFactors=FALSE)
head(my.pdata)

colnames(my.pdata)

head(my.pdata[, c("title", "geo_accession", "description")])

my.pdata <- my.pdata[, c("title", "geo_accession", "description")]
head(my.pdata)

rownames(my.pdata) %in% my.cels

head(rownames(my.pdata))

head(my.cels)

temp.rownames <- paste(rownames(my.pdata), ".CEL.gz", sep="")
temp.rownames

temp.rownames %in% my.cels

rownames(my.pdata) <- temp.rownames
rownames(my.pdata) %in% my.cels
rm(temp.rownames)

my.pdata

write.table(my.pdata, file=paste0(cel.path, my.gse, "_SelectPhenoData.txt"), sep="\t", quote=FALSE)
list.files(recursive=TRUE)

library(affy)

my.affy <- ReadAffy(celfile.path=cel.path, phenoData=paste(cel.path, paste0(my.gse, "_SelectPhenoData.txt"), sep=""))

show(my.affy)

head(exprs(my.affy))

dim(exprs(my.affy))

pData(my.affy)

table(pData(my.affy)$description)

pData(my.affy)$sample.levels <- c(rep("con", 3), rep("erbb2", 3), rep("egf", 3))
table(pData(my.affy)$sample.levels)
table(pData(my.affy)$description, pData(my.affy)$sample.levels)
pData(my.affy)$sample.labels <- c(paste("con", 1:3, sep="."), paste("erbb2", 1:3, sep="."), paste("egf", 1:3, sep="."))
table(pData(my.affy)$sample.labels)

cbind(pData(my.affy)$title, pData(my.affy)$sample.labels)

my.rma <- rma(my.affy, normalize=FALSE, background=FALSE)

show(my.rma)

head(exprs(my.rma))

head(exprs(my.geo.gse))

pData(my.rma)

pData(my.rma)$sample.levels <- as.factor(pData(my.rma)$sample.levels)
pData(my.rma)$sample.levels <- relevel(pData(my.rma)$sample.levels, ref="con")
levels(pData(my.rma)$sample.levels)

plotDensity(exprs(my.rma))

library(limma)

library(RColorBrewer)
display.brewer.all(colorblindFriendly=TRUE)

level.pal <- c("blue", "red","orange")
level.cols <- level.pal[pData(my.rma)$sample.levels]
level.cols

plotDensities(exprs(my.rma), legend=FALSE, col=level.cols, main="Arrays Not Normalized")
legend("topright", legend=levels(pData(my.rma)$sample.levels), fill=level.pal)

boxplot(exprs(my.rma), las=2, names=pData(my.rma)$sample.labels, outline=FALSE, col=level.cols, 
        main="Arrays Not Normalized")

pdf(file="results/DensityNoNorm.pdf", w=6, h=6)
plotDensities(exprs(my.rma), legend=FALSE, col=level.cols, main="Arrays Not Normalized")
legend("topright", legend=levels(pData(my.rma)$sample.levels), fill=level.pal)
dev.off()

pdf(file="results/DensityNoNorm.pdf", w=6, h=6)
boxplot(exprs(my.rma), las=2, names=pData(my.rma)$sample.labels, outline=FALSE, col=level.cols, 
        main="Arrays Not Normalized")
dev.off()

my.rma <- rma(my.affy, normalize=TRUE, background=TRUE)
pData(my.rma)$sample.levels <- as.factor(pData(my.rma)$sample.levels)
pData(my.rma)$sample.levels <- relevel(pData(my.rma)$sample.levels, ref="con")

plotDensities(exprs(my.rma), legend=FALSE, col=level.cols, main="Arrays Normalized")
legend("topright", legend=levels(pData(my.rma)$sample.levels), fill=level.pal)

boxplot(exprs(my.rma), las=2, names=pData(my.rma)$sample.labels, outline=FALSE, col=level.cols, 
        main="Arrays Normalized")

pdf(file="results/DensityNorm.pdf", w=6, h=6)
plotDensities(exprs(my.rma), legend=FALSE, col=level.cols, main="Arrays Normalized")
legend("topright", legend=levels(pData(my.rma)$sample.levels), fill=level.pal)
dev.off()

pdf(file="results/DensityNorm.pdf", w=6, h=6)
boxplot(exprs(my.rma), las=2, names=pData(my.rma)$sample.labels, outline=F, col=level.cols, main="Arrays Normalized")
dev.off()

rma.comp <- as.data.frame(exprs(my.rma))
head(rma.comp)

colnames(rma.comp) <- pData(my.rma)$sample.labels
head(rma.comp)

rma.comp$PROBEID <- rownames(rma.comp)
head(rma.comp)

dim(rma.comp)

if(colnames(rma.comp)[1] != "PROBEID") rma.comp <- rma.comp[, c(10, 1:9)]
head(rma.comp)

write.table(rma.comp, file=paste0("results/", my.gse, "_RMA_Norm_Complete.txt"), sep="\t", quote=FALSE, 
            row.names=FALSE)

my.calls <- mas5calls(my.affy)
head(exprs(my.calls))

table(exprs(my.calls)[, 1])

table(exprs(my.calls)[3, ])

table(exprs(my.calls)[3, ], pData(my.calls)$sample.levels)

present <- apply(exprs(my.calls), 1, function(x)(sum(x == "P")))

table(present)

prop.table(table(present)) * 100

#the number in this line of code should reflect the number of biological replicates; get rid of the crap
table(present >= 3)

plotDensities(exprs(my.rma)[present >= 3, ], col=level.cols, legend=F, main="Present >= 3")

plotDensities(exprs(my.rma)[present < 3, ], col=level.cols, legend=F, main="Present < 3")

dim(my.rma)

length(present)

if(nrow(my.rma) == length(present)) my.rma <- my.rma[present >= 3, ]

dim(my.rma)

temp.grep <- grep("_x_", rownames(exprs(my.rma)), fixed=TRUE)
if(length(temp.grep) > 0) my.rma <- my.rma[-temp.grep]
dim(my.rma)

temp.grep <- grep("AFFX", rownames(exprs(my.rma)), fixed=TRUE)
if(length(temp.grep) > 0) my.rma <- my.rma[-temp.grep]
dim(my.rma)

library(hexbin)

my.level <- c("con")
hexplom(exprs(my.rma[, pData(my.rma)$sample.levels %in% my.level]))

rma.cor <- cor(exprs(my.rma), method="pearson")
rma.cor

gene.mean <- apply(exprs(my.rma), 1, mean)
gene.sd <- apply(exprs(my.rma), 1, sd)
zscores <- sweep(exprs(my.rma), 1, gene.mean, "-")
zscores <- sweep(zscores, 1, gene.sd, "/")

zscore.cor <- cor(zscores, method="pearson")
zscore.cor

library(corrplot)

rownames(rma.cor) <- colnames(rma.cor) <- pData(my.rma)$sample.labels
rownames(zscore.cor) <- colnames(zscore.cor) <- pData(my.rma)$sample.labels

corrplot(rma.cor, method="ellipse", order="original")

corrplot(rma.cor, method="ellipse", order="hclust", hclust.method="average")

corrplot(zscore.cor, method="ellipse", order="hclust", hclust.method="average")

my.dist <- as.dist(1 - zscore.cor)

my.hclust <- hclust(my.dist, method="average")

plot(my.hclust, cex=0.75, main="Comparison of Biological Replicates")

layout(matrix(c(1, 2), ncol=1), heights=c(0.3, 0.7))
par(mar=c(0, 6.5, 4, 6.6))
plot((as.dendrogram(my.hclust)), type="rectangle", horiz=FALSE, leaflab="none",yaxt="n",xaxt="n", xaxs="i", yaxs="i")
par(mar=c(4, 4, 0, 4))
corrplot(zscore.cor, method="ellipse", order="hclust", hclust.method="average")


top.genes <- 1:250
gene.sd.sort <- sort(gene.sd, decreasing=TRUE)
zscores.top <- zscores[names(gene.sd.sort)[top.genes], ]
zscores.top.cor <- cor(zscores.top, method="pearson")
rownames(zscores.top.cor) <- colnames(zscores.top.cor) <- pData(my.rma)$sample.labels
my.dist <- as.dist(1 - zscores.top.cor)
my.hclust <- hclust(my.dist, method="average")
layout(matrix(c(1, 2), ncol=1), heights=c(0.3, 0.7))
par(mar=c(0, 6.5, 4, 6.6))
plot((as.dendrogram(my.hclust)), type="rectangle", horiz=FALSE, leaflab="none",yaxt="n",xaxt="n", xaxs="i", yaxs="i")
par(mar=c(4, 4, 0, 4))
corrplot(zscores.top.cor, method="ellipse", order="hclust", hclust.method="average")

sample.no <- 1:9
plotMDS(exprs(my.rma)[, sample.no], labels=pData(my.rma)$sample.labels[sample.no], top=250, 
        gene.selection="pairwise", main="MDS Plot to Compare Replicates")

pdf(file="results/MDS_plot.pdf", w=6, h=6)
plotMDS(exprs(my.rma), labels=pData(my.rma)$sample.labels, top=500, gene.selection="common", 
        main="MDS Plot to Compare Replicates")
dev.off()

my.design <- model.matrix(~0 + sample.levels, pData(my.rma))
my.design[1:9, ]

rownames(my.design) <- pData(my.rma)$sample.labels
colnames(my.design) <- levels(pData(my.rma)$sample.levels)
my.design[1:9, ]

my.fit <- lmFit(my.rma, my.design)
my.fit

my.contrasts <- makeContrasts(erbb2vcon = erbb2-con, egfvcon = egf-con,levels = my.design)
my.contrasts

contrast.fits <- sapply(colnames(my.contrasts), function(x)(contrasts.fit(my.fit, contrasts=my.contrasts[, x])))
length(contrast.fits)

names(contrast.fits)

contrast.fits["erbb2vcon"]

contrast.ebs <- lapply(contrast.fits, function(x)(eBayes(x, proportion=0.1, trend=FALSE, robust=FALSE)))
contrast.ebs[["erbb2vcon"]]

contrast.tts <- lapply(contrast.ebs, function(x)(topTable(x, adjust="BH", number=length(x$coefficients), 
                                                          sort.by="none")))
head(contrast.tts[["erbb2vcon"]])

contrast.tests <- lapply(contrast.ebs, function(x)(decideTests(x, method="separate", adjust.method="BH", 
                                                               p.value=0.05, lfc=0)))
table(contrast.tests[["erbb2vcon"]])

tests.mat <- do.call(cbind, contrast.tests)
colnames(tests.mat) <- names(contrast.tests)
head(tests.mat)

ma.cols <- brewer.pal(9, "RdBu")[c(2, 9)]
my.contrast <- "erbb2vcon"

##Note: may need to use the form limma::plotMA to get the plotMA function from limma.  This is usually not necessary
plotMA(contrast.ebs[[my.contrast]], status=contrast.tests[[my.contrast]], hl.pch=20, hl.col=ma.cols, 
       main=paste(my.gse, my.contrast, sep=" "))

ma.cols <- c(brewer.pal(11, "RdBu")[10], "grey50", brewer.pal(11, "RdBu")[2])
my.cols <- ma.cols[as.factor(contrast.tests[[my.contrast]][, 1])]
table(my.cols)

plot(contrast.tts[[my.contrast]]$logFC, -log10(contrast.tts[[my.contrast]]$adj.P.Val), col=my.cols, 
     pch=19, main=paste(my.gse, my.contrast, sep=" "), xlab="log2 Fold Change", ylab="-log10 Adjusted P-Value")

up.probes <- rownames(contrast.tests[[my.contrast]])[contrast.tests[[my.contrast]][, 1] == 1]
plotMDS(exprs(my.rma)[up.probes, ], pch=rep(21:24, 6), bg=level.cols, labels=NULL, top=length(up.probes), 
        main=paste(my.gse, "MDS Plot to Compare Up-Regulated Genes"), cex=1.5, cex.main=0.9, sub=my.contrast)
legend1 <- c("array1", "array2", "array3")
legend2 <- c("con","erbb2","egf")
legend("top", horiz=T, bty="n", legend=legend1, pch=c(21:24), cex=0.8)
legend("bottom", horiz=T, bty="n", legend=legend2, fill=level.pal, cex=0.8)

down.probes <- rownames(contrast.tests[[my.contrast]])[contrast.tests[[my.contrast]][, 1] == -1]
plotMDS(exprs(my.rma)[down.probes, ], pch=rep(21:24, 6), bg=level.cols, labels=NULL, top=length(down.probes), 
        main=paste(my.gse, "MDS Plot to Compare Down-Regulated Genes"), cex=1.5, cex.main=0.9, sub=my.contrast)
legend1 <- c("array1", "array2", "array3")
legend2 <- c("con","erbb2","egf")
legend("top", horiz=T, bty="n", legend=legend1, pch=c(21:24), cex=0.8)
legend("bottom", horiz=T, bty="n", legend=legend2, fill=level.pal, cex=0.8)

deg.probes <- union(up.probes, down.probes)
plotMDS(exprs(my.rma)[deg.probes, ], pch=rep(21:24, 6), bg=level.cols, labels=NULL, top=length(deg.probes), 
        main=paste(my.gse, "MDS Plot to Compare DEGs"), cex=1.5, cex.main=0.9, sub=my.contrast)
legend1 <- c("array1", "array2", "array3")
legend2 <- c("con","erbb2","egf")
legend("top", horiz=T, bty="n", legend=legend1, pch=c(21:24), cex=0.8)
legend("bottom", horiz=T, bty="n", legend=legend2, fill=level.pal, cex=0.8)

not.deg.probes <- setdiff(rownames(exprs(my.rma)), deg.probes)
plotMDS(exprs(my.rma)[not.deg.probes, ], pch=rep(21:24, 6), bg=level.cols, labels=NULL, top=1000, 
        main=paste(my.gse, "MDS Plot to Compare Unchanged Genes"), cex=1.5, cex.main=0.9, sub=my.contrast)
legend1 <- c("array1", "array2", "array3")
legend2 <- c("con","erbb2","egf")
legend("top", horiz=T, bty="n", legend=legend1, pch=c(21:24), cex=0.8)
legend("bottom", horiz=T, bty="n", legend=legend2, fill=level.pal, cex=0.8)

test.sum <- apply(tests.mat, 1, function(x)(sum(abs(x))))
table(test.sum)

deg.probes <- rownames(tests.mat)[test.sum > 0]
head(deg.probes)

my.probe <- deg.probes[100]
boxplot(exprs(my.rma)[my.probe, ] ~ pData(my.rma)$sample.levels, col=level.pal, main=my.probe, ylab="RMA Value", 
        xlab="Treatment")

library (hgu133plus2.db)

my.probe <- deg.probes[100]
gene.symbol <- select(hgu133plus2.db, key=my.probe, keytype="PROBEID", columns="SYMBOL")
gene.symbol <- gene.symbol$SYMBOL
boxplot(exprs(my.rma)[my.probe, ] ~ pData(my.rma)$sample.levels, col=level.pal, main=gene.symbol, ylab="RMA Value", 
        xlab="Treatment")
stripchart(exprs(my.rma)[my.probe, ] ~ pData(my.rma)$sample.levels, vertical=TRUE, method="jitter", jitter=0.1, 
          add=TRUE, pch=1)

library (hgu133plus2.db)

hgu133plus2.db

columns(hgu133plus2.db)

keytypes(hgu133plus2.db)

single.gene <- select(hgu133plus2.db, keys=rownames(contrast.tts[["erbb2vcon"]])[1], keytype="PROBEID", 
                    columns=c("ACCNUM", "SYMBOL", "ALIAS"))
dim(single.gene)
head(single.gene)

gene.data <- select(hgu133plus2.db, keys=rownames(contrast.tts[["erbb2vcon"]]), keytype="PROBEID", 
                    columns=c("ENTREZID", "GENENAME", "SYMBOL"))
head(gene.data)

dim(gene.data)

dim(contrast.tts[["erbb2vcon"]])

summary(duplicated(gene.data$PROBEID))

summary(duplicated(gene.data$SYMBOL))

summary(duplicated(gene.data$ENTREZID))

table(duplicated(gene.data$SYMBOL), duplicated(gene.data$PROBEID), useNA="always", dnn=list("gene", "probe"))

table(duplicated(gene.data$ENTREZID), duplicated(gene.data$PROBEID), useNA="always", dnn=list("entrez", "probe"))

table(duplicated(gene.data$ENTREZID), duplicated(gene.data$SYMBOL), useNA="always", dnn=list("entrez", "symbol"))

deg.data <- gene.data[gene.data$PROBEID %in% deg.probes, ]
rownames(deg.data) <- NULL
head(deg.data, 10)

my.symbol <- deg.data$SYMBOL[7]
#my.symbol <- "PLK1"

my.symbol

my.probe <- deg.data$PROBEID[deg.data$SYMBOL %in% my.symbol][1] ##insures single probeset
boxplot(exprs(my.rma)[my.probe, ] ~ pData(my.rma)$sample.levels, col=level.pal, main=my.symbol, ylab="RMA Value", 
        xlab="Treatment")
stripchart(exprs(my.rma)[my.probe, ] ~ pData(my.rma)$sample.levels, vertical=TRUE, method="jitter", jitter=0.1, 
          add=TRUE, pch=1)

my.results <- lapply(names(contrast.tts), function(x)(cbind(contrast.tts[[x]], contrast.tests[[x]])))
names(my.results) <- names(contrast.tts)
for(i in 1:length(my.results)){
  my.results[[i]]$PROBEID <- rownames(my.results[[i]])
  colnames(my.results[[i]])[7] <- "test"
}

my.results <- lapply(my.results, function(x)(merge(gene.data, x, by="PROBEID")))
head(my.results[[1]])
for(i in 1:length(my.results)){
  write.table(my.results[[i]], file=paste0("results/",my.gse,"_", names(my.results)[i], "_Limma_Out.txt"), 
              row.names=FALSE, sep="\t", quote=FALSE)
}

my.coeff <- as.data.frame(my.fit$coefficients)
head(my.coeff)

my.coeff$PROBEID <- rownames(my.coeff)
head(my.coeff)

my.coeff <- merge(gene.data, my.coeff, by="PROBEID")
head(my.coeff)

table(duplicated(my.coeff$PROBEID))

write.table(my.coeff, file=paste0("results/", my.gse, "_Limma_Coeff.txt"), sep="\t", quote=F)

rma.filt <- as.data.frame(exprs(my.rma))
rma.filt$PROBEID <- rownames(rma.filt)
rma.filt <- merge(gene.data, rma.filt, by="PROBEID")
head(rma.filt)

write.table(rma.filt, file=paste0("results/", my.gse, "_RMA_Norm_Filtered.txt"), sep="\t", quote=FALSE)

write.table(tests.mat, file=paste0("results/",my.gse,"_GeneCalls.txt"), sep="\t", quote=F)

sessionInfo()


