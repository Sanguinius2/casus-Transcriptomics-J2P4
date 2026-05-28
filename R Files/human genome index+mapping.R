#eerst naar juiste directorie,dus waar je al je bestanden hebt staan
setwd("C:/Users/sande/OneDrive - NHL Stenden/school/jaar 2/Periode 4/casus/human refseq")
#controle of het inderdaad juiste mapje is
getwd()
#install all the needed packages
install.packages('BiocManager')
BiocManager::install('Rsubread')
BiocManager::install('Rsamtools')
BiocManager::install("DESeq2")
BiocManager::install("KEGGREST")
BiocManager::install("EnhancedVolcano")
BiocManager::install("pathview")
#Load the needed packages for this part
library(BiocManager)
library(Rsubread)
library(Rsamtools)
#load the human genome downloaded from online sources
buildindex(
  basename = 'ref_human_genome',
  reference = 'GCF_000001405.40_GRCh38.p14_genomic.fna',
  memory = 4000,
  indexSplit = TRUE)
#align sample with human genome and store in .BAM files
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
#REFVANDOCENT, this is another database with better values, you can keep using the previous one
counts <- read.table(
  "count_matrix_RA.txt",
  header = TRUE
  )
#load all the needed packages for this part
library(DESeq2)
library(KEGGREST)
library(EnhancedVolcano)
library(pathview)

#names for table+ putting names in the table with the data
treatment <- c("Ref_norm", "Ref_norm", "Ref_norm", "Ref_norm", "Ref_reuma", "Ref_reuma", "Ref_reuma", "Ref_reuma")
treatment_table <- data.frame(treatment)
rownames(treatment_table) <- c("SRR4785819","SRR4785820",	"SRR4785828",	"SRR4785831",	"SRR4785979",	"SRR4785980",	"SRR4785986",	"SRR4785988")
#make the DESeq Data set
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = treatment_table,
                              design = ~ treatment)
# start the analysis
dds <- DESeq(dds)
resultaten <- results(dds)
#make a table of the results
write.table(resultaten, file = 'ResultatenDOCENTGEGEVEN.csv', row.names = TRUE, col.names = TRUE)
#check how many genes really changed due to Rheuma
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange > 1, na.rm = TRUE)
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange < -1, na.rm = TRUE)
#check what the most significant genes are
hoogste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = TRUE), ]
laagste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = FALSE), ]
laagste_p_waarde <- resultaten[order(resultaten$padj, decreasing = FALSE), ]
hoogste_fold_change
laagste_fold_change
laagste_p_waarde
#makes a volcano plot for visualizing the significant genes
EnhancedVolcano(resultaten,
                lab = rownames(resultaten),
                x = 'log2FoldChange',
                y = 'padj')
#saves the volcano plot on your device
dev.copy(png, 'VolcanoplotDOCENTGEGEVEN.png', 
         width = 8,
         height = 10,
         units = 'in',
         res = 500)
dev.off()

