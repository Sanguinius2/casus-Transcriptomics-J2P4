#eerst naar juiste directorie,dus waar je al je bestanden hebt staan
setwd("C:/Users/sande/OneDrive - NHL Stenden/school/jaar 2/Periode 4/casus/verwerkte data/.BAM files")
#controle of het inderdaad juiste mapje is
getwd()
#install alle packages die nodig zijn
install.packages('BiocManager')
BiocManager::install('Rsubread')
BiocManager::install('Rsamtools')
BiocManager::install("DESeq2")
BiocManager::install("KEGGREST")
BiocManager::install("EnhancedVolcano")
BiocManager::install("pathview")
BiocManager::install(c(
  "goseq",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "org.Hs.eg.db",
  "biomaRt"
))
# Laad de benodigde libraries voor dit R script
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
library(DESeq2)
library(KEGGREST)
library(EnhancedVolcano)
library(pathview)
library(BiocManager)
library(Rsubread)
library(Rsamtools)
#laad hett humaane genoom, dat je hebt gedownload. zie README voor meer info
buildindex(
  basename = 'ref_human_genome',
  reference = 'GCF_000001405.40_GRCh38.p14_genomic.fna',
  memory = 4000,
  indexSplit = TRUE)
#align sample met humaan genoom en zet het in .BAM files
align.human1 <- align(index = "ref_human_genome", readfile1 = "SRR4785819_1_subset40k.FASTQ", readfile2 = "SRR4785819_2_subset40k.FASTQ", output_file = "Ref_human_genome1.BAM")
align.human2 <- align(index = "ref_human_genome", readfile1 = "SRR4785820_1_subset40k.FASTQ", readfile2 = "SRR4785820_2_subset40k.FASTQ", output_file = "Ref_human_genome2.BAM")
align.human3 <- align(index = "ref_human_genome", readfile1 = "SRR4785828_1_subset40k.FASTQ", readfile2 = "SRR4785828_2_subset40k.FASTQ", output_file = "Ref_human_genome3.BAM")
align.human4 <- align(index = "ref_human_genome", readfile1 = "SRR4785831_1_subset40k.FASTQ", readfile2 = "SRR4785831_2_subset40k.FASTQ", output_file = "Ref_human_genome4.BAM")
align.human_RA1 <- align(index = "ref_human_genome", readfile1 = "SRR4785979_1_subset40k.FASTQ", readfile2 = "SRR4785979_2_subset40k.FASTQ", output_file = "Ref_human_genome_RA1.BAM")
align.human_RA2 <- align(index = "ref_human_genome", readfile1 = "SRR4785980_1_subset40k.FASTQ", readfile2 = "SRR4785980_2_subset40k.FASTQ", output_file = "Ref_human_genome_RA2.BAM")
align.human_RA3 <- align(index = "ref_human_genome", readfile1 = "SRR4785986_1_subset40k.FASTQ", readfile2 = "SRR4785986_2_subset40k.FASTQ", output_file = "Ref_human_genome_RA3.BAM")
align.human_RA4 <- align(index = "ref_human_genome", readfile1 = "SRR4785988_1_subset40k.FASTQ", readfile2 = "SRR4785988_2_subset40k.FASTQ", output_file = "Ref_human_genome_RA4.BAM")

#zet de .BAM files in een variable en maak hier een count matrix van
allsampleshumangenome <- c("Ref_human_genome1.BAM", "Ref_human_genome2.BAM", "Ref_human_genome3.BAM", "Ref_human_genome4.BAM", "Ref_human_genome_RA1.BAM", "Ref_human_genome_RA2.BAM", "Ref_human_genome_RA3.BAM", "Ref_human_genome_RA4.BAM")
count_matrix <- featureCounts(
  files = allsampleshumangenome,
  annot.ext = "genomic.gtf",
  isPairedEnd = TRUE,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "gene", 
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE
)

counts <- count_matrix$counts
View(counts)
head(counts)
colnames(counts) <- c("Ref_human_genome1.BAM", "Ref_human_genome2.BAM", "Ref_human_genome3.BAM", "Ref_human_genome4.BAM", "Ref_human_genome_RA1.BAM", "Ref_human_genome_RA2.BAM", "Ref_human_genome_RA3.BAM", "Ref_human_genome_RA4.BAM")
head(counts)
write.csv(counts, "Ref_Human_genome.csv")
#REFORIGINAL
counts <- read.csv("Ref_Human_genome.csv", row.names = 1)
# Vanaf hier wordt een andere dataset gebruikt met makkelijkere data om mee te werken
# Deze matrix heeft dezzelfde info, maar de stappen gaan snellen, dus makkelijker uit te leggen wat overal gebeurd
# Bij het volgen van dit script kun je gewoon de eigen dataset gebbruiken en de code regel hieronder negeren
counts <- read.table(
  "count_matrix_RA.txt",
  header = TRUE
  )

#namen voor de tabel en zet deze namen in de tabel met bijhoorende samples
treatment <- c("Ref_norm", "Ref_norm", "Ref_norm", "Ref_norm", "Ref_reuma", "Ref_reuma", "Ref_reuma", "Ref_reuma")
treatment_table <- data.frame(treatment)
rownames(treatment_table) <- c("SRR4785819","SRR4785820",	"SRR4785828",	"SRR4785831",	"SRR4785979",	"SRR4785980",	"SRR4785986",	"SRR4785988")
#maak de DESeq Data set
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = treatment_table,
                              design = ~ treatment)
# start the analysis met DESeq
dds <- DESeq(dds)
resultaten <- results(dds)
#maak tabel van resultaten
write.table(resultaten, file = 'Resultaten.csv', row.names = TRUE, col.names = TRUE)
#kijk hoeveel genen zijn upregulated(>1) of downregulated(<-1)
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange > 1, na.rm = TRUE)
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange < -1, na.rm = TRUE)
#kijk wat de meest significante genes zijn
hoogste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = TRUE), ]
laagste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = FALSE), ]
laagste_p_waarde <- resultaten[order(resultaten$padj, decreasing = FALSE), ]
hoogste_fold_change
laagste_fold_change
laagste_p_waarde
#maak een volcano plot om de meest significante genes te laten zien
EnhancedVolcano(resultaten,
                lab = rownames(resultaten),
                x = 'log2FoldChange',
                y = 'padj')
#save de volcano plot, zelf naam van png aanpassen kan
dev.copy(png, 'Volcanoplot.png', 
         width = 8,
         height = 10,
         units = 'in',
         res = 500)
dev.off()

# Read the DESeq2 results
results <- read.csv(
  "Resultaten.csv",
  row.names = 1,
  sep = " "
)

head(results)
view(results)
#scheid de significante genes van de niet significante, vanag hier gebruiken we alleen de significante genes
sigGenes <- as.integer(
  !is.na(results$padj) &
    results$padj < 0.05 &
    results$log2FoldChange > 1
)

# gebruik gen IDs als naam
names(sigGenes) <- rownames(results)

table(sigGenes)
#get gene lengtes van de dataset
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
#match de gene lengte aan de eigen resultaten
lengthVector <- geneLengths$length

names(lengthVector) <- geneLengths$hgnc_symbol

# Match deze volgorde met de eerdere resultaten file
lengthVector <- lengthVector[rownames(results)]
#removes genes without a length(NA)
keep <- !is.na(lengthVector)
sigGenes_filtered <- sigGenes[keep]
lengthVector_filtered <- lengthVector[keep]
#kijk of elke stap goed is gegaan, antwoord op deze code moet 0 zijn
sum(is.na(lengthVector_filtered))
#store de gene lengte van de database in een variabele, hier is hg38 gebruikt, kan ook andere zijn
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
#extra check, staat er iets in de net gemaakt GoResults
head(goResults)
#maak de plot van de GO-analysis 
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


#Pathway analysis, zet eerder gemaakte resultaten opnieuw in de resultaten variable
resultaten <- read.csv(
  "Resultaten.csv",
  row.names = 1,
  sep = " "
)
resultaten[1] <- NULL
resultaten[2:5] <- NULL
#zet de foldchange in een variable
geneList <- resultaten$log2FoldChange
#geef dezelfde rownames als de resultaten table
names(geneList) <- rownames(resultaten)

head(geneList)

#zet gene symbols naar Entrez IDs om

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

# verwijder NA IDs
geneData <- geneData[!is.na(names(geneData))]

#KEGG pathway finder, niet nodig, maar handig om te gebruiken
keggList("pathway", "hsa")

#run pathway on hsa, het kan zijn dat je zelf andere significcante genen krijft
#pas het pathway.id dan aan en kijk of de species nog steeds klopt
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