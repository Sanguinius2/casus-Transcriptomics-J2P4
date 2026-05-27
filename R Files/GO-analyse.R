# zet dit naar eigen directorie
setwd("C:/Users/sande/OneDrive - NHL Stenden/school/jaar 2/Periode 4/casus/human refseq")
getwd()
# Laad de required libraries
library(tidyverse)
library(goseq)
library(GO.db)
library(biomaRt)

# Lees DESeq2 resultaten in
results <- read.csv(
  "ResultatenDOCENTGEGEVEN.csv",
  row.names = 1,
  sep = " "
)

head(results)
view(results)

sigGenes <- as.integer(
  !is.na(results$padj) &
    results$padj < 0.05 &
    results$log2FoldChange > 1
)

# Gebruik gen IDs als namen
names(sigGenes) <- rownames(results)

table(sigGenes)
BiocManager::install(c(
  "goseq",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "org.Hs.eg.db",
  "biomaRt"
))

library(biomaRt)

mart <- useMart(
  "ensembl",
  dataset = "hsapiens_gene_ensembl"
)

geneLengths <- getBM(
  attributes = c(
    "hgnc_symbol",
    "start_position",
    "end_position"
  ),
  filters = "hgnc_symbol",
  values = rownames(results),
  mart = mart
)

# Bereken genlengte
geneLengths$length <- geneLengths$end_position -
  geneLengths$start_position

head(geneLengths)
#2
lengthVector <- geneLengths$length

names(lengthVector) <- geneLengths$hgnc_symbol

# Match volgorde met eerdere results bestand
lengthVector <- lengthVector[rownames(results)]
keep <- !is.na(lengthVector)

sigGenes_filtered <- sigGenes[keep]
lengthVector_filtered <- lengthVector[keep]
#test of het goed werkt
sum(is.na(lengthVector_filtered))
#is goed? ga door
summary(lengthVector)
pwf <- nullp(
  DEgenes = sigGenes_filtered,
  genome = "hg38",
  id = "geneSymbol",
  bias.data = lengthVector_filtered
)
#NEXT STEP
goResults <- goseq(
  pwf,
  genome = "hg38",
  id = "geneSymbol",
  test.cats = c("GO:BP")
)

head(goResults)
#ggplot time
library(ggplot2)
library(dplyr)

topGO <- goResults %>%
  arrange(over_represented_pvalue) %>%
  slice_head(n = 10) %>%
  mutate(
    hitsPerc = numDEInCat * 100 / numInCat
  )

ggplot(
  topGO,
  aes(
    x = hitsPerc,
    y = reorder(term, hitsPerc),
    colour = over_represented_pvalue,
    size = numDEInCat
  )
) +
  geom_point() +
  theme_bw() +
  labs(
    title = "Top 10 GO Terms",
    x = "Hits (%)",
    y = "GO term"
  )
dev.copy(png, 'GO-analyse-Resultaten.png')
dev.off()


#Pathway analyse
resultaten <- read.csv(
  "ResultatenDOCENTGEGEVEN.csv",
  row.names = 1,
  sep = " "
)
resultaten[1] <- NULL
resultaten[2:5] <- NULL
#NEW data codex from ct
geneList <- resultaten$log2FoldChange

names(geneList) <- rownames(resultaten)

head(geneList)
# install als nog niet gedaan
BiocManager::install("org.Hs.eg.db")
#Convert gene symbols to Entrez IDs
library(org.Hs.eg.db)
library(AnnotationDbi)

geneIDs <- mapIds(
  org.Hs.eg.db,
  keys = rownames(resultaten),
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)
#gene vector
geneData <- resultaten$log2FoldChange

names(geneData) <- geneIDs

# remove NA IDs
geneData <- geneData[!is.na(names(geneData))]
#run pathway on hsa
library(pathview)

pathview(
  gene.data = geneData,
  pathway.id = "hsa05323",
  species = "hsa",
  gene.idtype = "entrez",
  limit = list(gene = 5)
)
#another one for top 1 GO
pathview(
  gene.data = geneData,
  pathway.id = "hsa04662",
  species = "hsa",
  gene.idtype = "entrez"
)
# RA pathway
pathview(
  gene.data = geneData,
  pathway.id = "hsa05323",
  species = "hsa",
  gene.idtype = "entrez"
)
#library voor zoeken
library(KEGGREST)

keggList("pathway", "hsa")
