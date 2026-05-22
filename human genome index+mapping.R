setwd("C:/Users/sande/OneDrive - NHL Stenden/school/jaar 2/Periode 4/casus/human refseq")
getwd()
library(Rsubread)
library(BiocManager)
library(Rsamtools)

buildindex(
  basename = 'ref_human_genome',
  reference = 'GCF_000001405.40_GRCh38.p14_genomic.fna',
  memory = 4000,
  indexSplit = TRUE)

align.human1 <- align(index = "ref_human_genome", readfile1 = "SRR4785819_1_subset40k.FASTQ", readfile2 = "SRR4785819_2_subset40k.FASTQ", output_file = "Ref_human_genome1.BAM")
align.human2 <- align(index = "ref_human_genome", readfile1 = "SRR4785820_1_subset40k.FASTQ", readfile2 = "SRR4785820_2_subset40k.FASTQ", output_file = "Ref_human_genome2.BAM")
align.human3 <- align(index = "ref_human_genome", readfile1 = "SRR4785828_1_subset40k.FASTQ", readfile2 = "SRR4785828_2_subset40k.FASTQ", output_file = "Ref_human_genome3.BAM")
align.human4 <- align(index = "ref_human_genome", readfile1 = "SRR4785831_1_subset40k.FASTQ", readfile2 = "SRR4785831_2_subset40k.FASTQ", output_file = "Ref_human_genome4.BAM")
align.human_RA1 <- align(index = "ref_human_genome", readfile1 = "SRR4785979_1_subset40k.FASTQ", readfile2 = "SRR4785979_2_subset40k.FASTQ", output_file = "Ref_human_genome_RA1.BAM")
align.human_RA2 <- align(index = "ref_human_genome", readfile1 = "SRR4785980_1_subset40k.FASTQ", readfile2 = "SRR4785980_2_subset40k.FASTQ", output_file = "Ref_human_genome_RA2.BAM")
align.human_RA3 <- align(index = "ref_human_genome", readfile1 = "SRR4785986_1_subset40k.FASTQ", readfile2 = "SRR4785986_2_subset40k.FASTQ", output_file = "Ref_human_genome_RA3.BAM")
align.human_RA4 <- align(index = "ref_human_genome", readfile1 = "SRR4785988_1_subset40k.FASTQ", readfile2 = "SRR4785988_2_subset40k.FASTQ", output_file = "Ref_human_genome_RA4.BAM")


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
#REFVANDOCENT
counts <- read.table(
  "count_matrix_RA.txt",
  header = TRUE
  )
library(DESeq2)
library(KEGGREST)
library(EnhancedVolcano)
library(pathview)


treatment <- c("Ref_norm", "Ref_norm", "Ref_norm", "Ref_norm", "Ref_reuma", "Ref_reuma", "Ref_reuma", "Ref_reuma")
treatment_table <- data.frame(treatment)
rownames(treatment_table) <- c("SRR4785819","SRR4785820",	"SRR4785828",	"SRR4785831",	"SRR4785979",	"SRR4785980",	"SRR4785986",	"SRR4785988")

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = treatment_table,
                              design = ~ treatment)

dds <- DESeq(dds)
resultaten <- results(dds)
write.table(resultaten, file = 'ResultatenDOCENTGEGEVEN.csv', row.names = TRUE, col.names = TRUE)
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange > 1, na.rm = TRUE)
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange < -1, na.rm = TRUE)
hoogste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = TRUE), ]
laagste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = FALSE), ]
laagste_p_waarde <- resultaten[order(resultaten$padj, decreasing = FALSE), ]
hoogste_fold_change
laagste_fold_change
laagste_p_waarde

EnhancedVolcano(resultaten,
                lab = rownames(resultaten),
                x = 'log2FoldChange',
                y = 'padj')
dev.copy(png, 'VolcanoplotDOCENTGEGEVEN.png', 
         width = 8,
         height = 10,
         units = 'in',
         res = 500)
dev.off()

