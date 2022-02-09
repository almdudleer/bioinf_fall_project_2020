# Title     : Differential expression of miRNA from human pituitary adenoma
# Objective : Find activated pathways
# Created by: Alexey.Serdyukov
# Created on: 17.12.2021
library(ggplot2)
library(readr)
library(dplyr)
library(data.table)
library(tidyr)
library(Biobase)
library(pcaExplorer)
library(pheatmap)
library(DESeq2)
library(limma)
library(ggrepel)
library(glue)

# Reading
counts <- read.table("build/final_counts/mirna_counts.tsv", header = T, sep = "\t")
cond <- read.table("data/clinical_data.txt", header = T, sep = "\t", row.names = 1)

tags <- sapply(colnames(counts), function(x) sub(".miRNAcount", "", basename(x)))
colnames(counts) <- tags

# Apply same ordering to conditions rows and counts columns
counts <- counts[, c("miRNA", rownames(cond))]

head(counts)
head(cond)
nrow(counts)

# Preprocessing
cond.convert.cols <- function(cond) {
  res <- copy(cond)
  res["Sex"][res["Sex"] == "F"] <- FALSE
  res["Sex"][res["Sex"] == "M"] <- TRUE
  res[res == "no" | res == "np" | res == "not norm"] <- FALSE
  res[res == "yes" | res == "norm"] <- TRUE
  res[res == "x"] <- NA
  res[] <- lapply(res, function(x) type.convert(as.character(x), as.is = TRUE))
  colnames(res) <- tolower(colnames(res))
  res <- res %>%
    dplyr::rename(
      male = sex,
      recurrence = recurrence..tumor.growth.,
      recurrence.2 = recurrence2..tumor.growth.,
      ifr.levels.norm = ifr.levels,
      tumor.size = size.of.tumor,
      tumor.volume = volume.of.tumor
    )
  res
}

## Rename columns and change column types
cond.prep <- cond.convert.cols(cond)
head(cond.prep)

## Calculate read counts per sample
col_sums <- colSums(counts[-which(names(counts) == "miRNA")])
col_sums

## Plot read counts per sample
counts[-which(names(counts) == "miRNA")] %>%
  gather %>%
  ggplot(aes(x = key, y = value)) +
  stat_summary(fun = sum, geom = "bar", fill = "darkblue") +
  labs(
    title = "Counts per sample",
    x = "Sample",
    y = "Counts"
  )

## Check that there is no NAs
table(is.na(counts))


# Create an expression set
es <- ExpressionSet(as.matrix(counts[, -1]))
fData(es) <- counts["miRNA"]
rownames(es) <- fData(es)$miRNA
pData(es) <- cond.prep
head(fData(es))
head(pData(es))
head(exprs(es))

## See some expression examples
exprs(es)[which(fData(es)$miRNA == "hsa-miR-99a-5p"),]
exprs(es)[grep("miR-34b", fData(es)$miRNA),]

# Filtering by mean expression > 1.0
counts.flt <- counts
counts.flt$mean.expr <- rowMeans(counts[, names(counts) != "miRNA"])
counts.flt <- counts.flt[counts.flt$mean.expr > 1.0,]
counts.flt <- counts.flt[order(counts.flt$mean.expr, decreasing = TRUE),]

nrow(counts.flt)
head(counts.flt)
counts.flt[grep("miR-34b", counts.flt$miRNA),]

# 480/2656 miRNAs left


# Quantile normalization
es.qnorm <- es

exprs(es.qnorm) <- normalizeBetweenArrays(log2(exprs(es.qnorm) + 1), method = "quantile")

fData(es.qnorm)$mean <- apply(exprs(es.qnorm), 1, mean)
es.qnorm <- es.qnorm[order(fData(es.qnorm)$mean, decreasing = TRUE),]
head(exprs(es.qnorm))

es.qnorm.top <- es.qnorm[1:480,]
# write.gct(es.qnorm.top, file="es_qnorm_top480.gct")

# PCA
pca <- prcomp(t(exprs(es.qnorm.top)))
pca.df <- cbind(cond.prep, pca$x)
for (target in colnames(cond.prep)) {
  print(
    ggplot(as.data.frame(pca.df), aes_string("PC1", "PC2", color = target)) +
      geom_point() +
      geom_label_repel(aes(label = rownames(pca.df)), force = 5, force_pull = 3) +
      ggtitle(glue("PCA: {target}"))
  )
}


# Clustering
bwr <- colorRampPalette(c("blue", "white", "red"))(50)

## Heatmap without clustering
pheatmap(exprs(es.qnorm.top), color = bwr, show_rownames = F, cluster_cols = F, cluster_rows = F)

## Clustered heatmap
pheatmap(exprs(es.qnorm.top), color = bwr, show_rownames = F)

## Kmeans
phmap <- pheatmap(exprs(es.qnorm.top), kmeans_k = 13)
pheatmap(exprs(es.qnorm.top)[order(phmap$kmeans$cluster),], color = bwr, cluster_rows = F, show_rownames = F)


# Scaling and clustering
df <- exprs(es)
df.qnorm <- normalizeBetweenArrays(log2(df + 1), method = "quantile")
means <- rowMeans(df.qnorm)
df.qnorm <- as.data.frame(df.qnorm)
df.qnorm$means <- means
df.qnorm <- df.qnorm[order(df.qnorm$mean, decreasing = TRUE),]
df.qnorm.top <- df.qnorm[0:480,]
df.qnorm.top$means <- NULL

df.qnorm.top.z <- t(scale(t(df.qnorm.top)))[-1,]

## Scaled heatmap without clustering
pheatmap(df.qnorm.top.z, color = bwr, show_rownames = F, cluster_cols = F, cluster_rows = F)

## Clustered scaled heatmap
pheatmap(df.qnorm.top.z, color = bwr, show_rownames = F)

## Scaled kmeans
phmap.z <- pheatmap(df.qnorm.top.z, kmeans_k = 13)
pheatmap(df.qnorm.top.z[order(phmap.z$kmeans$cluster),], color = bwr, cluster_rows = F, show_rownames = F)

# After PCA and clustering, N58, N74, N18, N47, N20, N48 are candidate outliers


# Removing outliers
outliers <- "N58"
counts.clean <- counts.flt[, !(names(counts.flt) %in% outliers)] %>% as.data.table()
cond.clean <- cond.prep[!(rownames(cond.prep) %in% outliers),]

es.clean <- ExpressionSet(as.matrix(counts.clean[, !c("mean.expr", "miRNA")]))
fData(es.clean) <- counts.clean[, "miRNA"]
rownames(es.clean) <- fData(es.clean)$miRNA
pData(es.clean) <- cond.clean

# DESeq2
deSeq <- function(es, target, contrast) {
  es.ds <- es[, !is.na(pData(es)[target])]
  es.ds <- es.ds[, order(pData(es.ds)[, target])]
  ds.countData <- exprs(es.ds)
  ds.colData <- pData(es.ds)
  ds.colData[sapply(ds.colData, is.logical)] <-
    lapply(ds.colData[sapply(ds.colData, is.logical)], as.factor)
  ds.colData[sapply(ds.colData, is.numeric)] <-
    lapply(ds.colData[sapply(ds.colData, is.numeric)], as.factor)

  dds <- DESeqDataSetFromMatrix(
    countData = ds.countData,
    colData = ds.colData,
    design = as.formula(glue("~ {target}"))
  )
  dds <- DESeq(dds)

  vst <- varianceStabilizingTransformation(dds)
  print(
    plotPCA(vst, intgroup = target) +
      ggtitle(glue("PCA: {target}")) +
      geom_label_repel(aes(label = colnames(vst)), force = 5, force_pull = 5)
  )

  de <- results(dds, contrast = c(target, contrast))
  de <- data.table(ID = rownames(de), as.data.table(de))
  de <- de[order(de$log2FoldChange),]

  # Heatmap
  phmap.de <- as.data.frame(df.qnorm.top.z) %>% select(rownames(pData(es.ds)))
  pheatmap(
    phmap.de[intersect(de$ID, rownames(df.qnorm.top.z)),],
    color = bwr, show_rownames = F, cluster_cols = F, cluster_rows = F,
    main = glue("DE heatmap: {target}")
  )

  list(de = de, dds = dds, vst = vst, target = target)
}

dir.create("build/de/", showWarnings = F)
dir.create("build/de/significant", showWarnings = F)
for (target in colnames(pData(es.clean))[-1]) {
  es.ds <- es.clean[, !is.na(pData(es.clean)[target])]
  contrast <- unlist(lapply(unique(pData(es.ds)[target]), as.character), use.names = FALSE)
  contrast <- contrast[order(contrast, decreasing = TRUE)]
  print(glue("Performing DE for target: {target}, with contrast: {glue_collapse(contrast,  sep = ', ')}"))
  de <- deSeq(es.ds, target, contrast)
  de.df <- de$de
  upreg <- de.df[de.df$padj < 0.05 & de.df$log2FoldChange > 0, c("ID", "padj", "log2FoldChange")]
  downreg <- de.df[de.df$padj < 0.05 & de.df$log2FoldChange < 0, c("ID", "padj", "log2FoldChange")]
  fwrite(de.df, file = glue("build/de/{target}_{contrast[1]}_vs_{contrast[2]}.tsv"), sep = "\t")
  if (nrow(upreg) > 0) {
    fwrite(upreg, file = glue("build/de/significant/{target}_{contrast[1]}_vs_{contrast[2]}.up.tsv"), sep = "\t")
  }
  if (nrow(downreg) > 0) {
    fwrite(downreg, file = glue("build/de/significant/{target}_{contrast[1]}_vs_{contrast[2]}.down.tsv"), sep = "\t")
  }
}

# TODO: volcano plots
# EnhancedVolcano(res,
#     lab = rownames(res),
#     x = 'log2FoldChange',
#     y = 'pvalue')
