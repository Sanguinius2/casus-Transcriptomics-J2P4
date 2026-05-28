# Set this to your own directory
setwd("C:/Users/sande/OneDrive - NHL Stenden/school/jaar 2/Periode 4/casus/human refseq")
getwd()
#install the necessary packages for this part
BiocManager::install(c(
  "goseq",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "org.Hs.eg.db",
  "biomaRt"
))

# Load de required libraries
library(tidyverse)
library(goseq)
library(GO.db)
library(biomaRt)
library(ggplot2)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(pathview)
library(KEGGREST)

# Read the DESeq2 results
results <- read.csv(
  "ResultatenDOCENTGEGEVEN.csv",
  row.names = 1,
  sep = " "
)

head(results)
view(results)
#keep the significant genes apart from the other ones, wich we don't use
sigGenes <- as.integer(
  !is.na(results$padj) &
    results$padj < 0.05 &
    results$log2FoldChange > 1
)

# Use gen IDs as names
names(sigGenes) <- rownames(results)

table(sigGenes)
#get the gene lengths
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

# calculate gene length
geneLengths$length <- geneLengths$end_position -
  geneLengths$start_position

head(geneLengths)
#match the gene length to your own results
lengthVector <- geneLengths$length

names(lengthVector) <- geneLengths$hgnc_symbol

# Match the order with earlier results file
lengthVector <- lengthVector[rownames(results)]
#removes genes without a length(NA)
keep <- !is.na(lengthVector)
sigGenes_filtered <- sigGenes[keep]
lengthVector_filtered <- lengthVector[keep]
#test if everything went correctly, answer should be 0
sum(is.na(lengthVector_filtered))
#call upon the gene lengths from a database
summary(lengthVector)
pwf <- nullp(
  DEgenes = sigGenes_filtered,
  genome = "hg38",
  id = "geneSymbol",
  bias.data = lengthVector_filtered
)
#GO-analysis
goResults <- goseq(
  pwf,
  genome = "hg38",
  id = "geneSymbol",
  test.cats = c("GO:BP")
)
#extra check
head(goResults)
#makes a plot from the GO-analysis 
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


#Pathway analysis
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

#Convert gene symbols to Entrez IDs

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
#KEGG pathway finder, not mandatory
keggList("pathway", "hsa")
