############################################################
# RNA-seq Analysis Pipeline for NPM1 Knockdown vs Control
############################################################

# Set working directory
cd /home/xxm_xxm/CJX_workspace
cd RNA-SEQ
ls

############################################################
# Install Conda Environment
############################################################

# Install Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# Check conda version
conda --version

# Create RNA-seq environment
conda create -n rnaseq python=3.10 -y

# Activate environment
conda activate rnaseq

# Install commonly used RNA-seq tools
conda install -c bioconda -c conda-forge \
fastqc multiqc star subread hisat2 salmon samtools -y

# Check installation
fastqc --version
STAR --version
featureCounts -v
samtools --version

############################################################
# Exit / Re-enter Environment
############################################################

conda deactivate
conda activate rnaseq

############################################################
# Create Output Directories
############################################################

mkdir merged qc align count

############################################################
# Merge Sequencing Lanes
############################################################

cat KD-1_L1_1.fq.gz KD-1_L3_1.fq.gz > merged/KD-1_R1.fastq.gz
cat KD-1_L1_2.fq.gz KD-1_L3_2.fq.gz > merged/KD-1_R2.fastq.gz

cat KD-2_L1_1.fq.gz KD-2_L3_1.fq.gz > merged/KD-2_R1.fastq.gz
cat KD-2_L1_2.fq.gz KD-2_L3_2.fq.gz > merged/KD-2_R2.fastq.gz

cat KD-3_L1_1.fq.gz KD-3_L3_1.fq.gz > merged/KD-3_R1.fastq.gz
cat KD-3_L1_2.fq.gz KD-3_L3_2.fq.gz > merged/KD-3_R2.fastq.gz

cat NC-1_L1_1.fq.gz NC-1_L3_1.fq.gz > merged/NC-1_R1.fastq.gz
cat NC-1_L1_2.fq.gz NC-1_L3_2.fq.gz > merged/NC-1_R2.fastq.gz

cat NC-2_L1_1.fq.gz NC-2_L3_1.fq.gz > merged/NC-2_R1.fastq.gz
cat NC-2_L1_2.fq.gz NC-2_L3_2.fq.gz > merged/NC-2_R2.fastq.gz

cat NC-3_L1_1.fq.gz NC-3_L3_1.fq.gz > merged/NC-3_R1.fastq.gz
cat NC-3_L1_2.fq.gz NC-3_L3_2.fq.gz > merged/NC-3_R2.fastq.gz

# Check merged files
ls -lh merged

############################################################
# Quality Control
############################################################

fastqc merged/*.fastq.gz -o qc -t 8
multiqc qc -o qc

# Download the generated HTML report
multiqc_report.html

############################################################
# Build Reference Genome
############################################################

mkdir -p reference/star_index
cd reference

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/GRCh38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz

gunzip GRCh38.primary_assembly.genome.fa.gz
gunzip gencode.v48.annotation.gtf.gz

############################################################
# Build STAR Index
############################################################

STAR \
--runThreadN 6 \
--runMode genomeGenerate \
--genomeDir ~/CJX_workspace/RNA-SEQ/reference/star_index \
--genomeFastaFiles ~/CJX_workspace/RNA-SEQ/reference/GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile ~/CJX_workspace/RNA-SEQ/reference/gencode.v48.annotation.gtf \
--sjdbOverhang 149

# Check index files
ls -lh ~/CJX_workspace/RNA-SEQ/reference/star_index

############################################################
# Read Alignment
############################################################

cd ~/CJX_workspace/RNA-SEQ

STAR \
--runThreadN 6 \
--genomeDir reference/star_index \
--readFilesIn merged/NC-1_R1.fastq.gz merged/NC-1_R2.fastq.gz \
--readFilesCommand zcat \
--outFileNamePrefix align/NC-1_ \
--outSAMtype BAM SortedByCoordinate

# Check mapping summary
cat align/NC-1_Log.final.out

# Align remaining samples
for sample in NC-2 NC-3 KD-1 KD-2 KD-3
do
STAR \
--runThreadN 6 \
--genomeDir reference/star_index \
--readFilesIn merged/${sample}_R1.fastq.gz merged/${sample}_R2.fastq.gz \
--readFilesCommand zcat \
--outFileNamePrefix align/${sample}_ \
--outSAMtype BAM SortedByCoordinate
done

############################################################
# Summarize Alignment Metrics
############################################################

for f in align/*_Log.final.out
do
echo "========== $f =========="
grep -E "Number of input reads|Average input read length|Uniquely mapped reads %|Average mapped length|% of reads mapped to multiple loci|% of reads mapped to too many loci|% of reads unmapped: too many mismatches|% of reads unmapped: too short|% of reads unmapped: other" $f
echo
done

# KD-2 showed abnormal unique mapping rate (44.38%), considered as an outlier candidate

############################################################
# Gene Counting
############################################################

featureCounts \
-T 4 \
-p \
--countReadPairs \
-B \
-C \
-a reference/gencode.v48.annotation.gtf \
-o count/gene_counts.txt \
align/*_Aligned.sortedByCoord.out.bam

# Check output files
ls count
head count/gene_counts.txt


############################################################
# Differential Expression Analysis in R
############################################################

R

library(DESeq2)
library(tidyverse)

# Read count matrix
counts <- read.delim("count/gene_counts.txt", comment.char="#")

# Extract count matrix
count_mat <- counts[,7:ncol(counts)]
rownames(count_mat) <- counts$Geneid

# Rename sample columns
colnames(count_mat) <- c("KD1","KD2","KD3","NC1","NC2","NC3")

# Sample metadata
coldata <- data.frame(
  row.names = colnames(count_mat),
  group = c("KD","KD","KD","NC","NC","NC")
)

# Construct DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = round(count_mat),
  colData = coldata,
  design = ~ group
)

# Filter low-expression genes
dds <- dds[rowSums(counts(dds)) >= 10,]

# Run differential expression analysis
dds <- DESeq(dds)

res <- results(dds, contrast=c("group","KD","NC"))
res <- as.data.frame(res)

# Save differential expression results for all samples
write.csv(res, "DEG_results_allSamples.csv")

############################################################
# PCA Plot for All Samples
############################################################

vsd <- vst(dds)

png("PCA_plot.png", width=1800, height=1400, res=200)
plotPCA(vsd, intgroup="group")
dev.off()

############################################################
# Remove Outlier Sample and Re-run Analysis
############################################################

# KD2 was considered an outlier candidate based on mapping/counting metrics and PCA
keep <- c("KD1","KD3","NC1","NC2","NC3")

count_mat2 <- count_mat[, keep]

coldata2 <- data.frame(
  row.names = keep,
  group = c("KD","KD","NC","NC","NC")
)

dds2 <- DESeqDataSetFromMatrix(
  countData = round(count_mat2),
  colData = coldata2,
  design = ~ group
)

dds2 <- dds2[rowSums(counts(dds2)) >= 10,]
dds2 <- DESeq(dds2)

res2 <- results(dds2, contrast=c("group","KD","NC"))
write.csv(as.data.frame(res2), "DEG_results_remove_KD2.csv")

############################################################
# PCA Plot After Removing KD2
############################################################

vsd2 <- vst(dds2)

png("PCA_remove_KD2.png", width=1800, height=1400, res=220)
plotPCA(vsd2, intgroup="group")
dev.off()

############################################################
# Identify Significant Differentially Expressed Genes
############################################################

res_df <- as.data.frame(res2)

# Stringent threshold
sig <- subset(res_df,
              padj < 0.05 &
              abs(log2FoldChange) > 1)

nrow(sig)

############################################################
# Annotate Ensembl IDs with Gene Symbols
############################################################

library(org.Hs.eg.db)
library(AnnotationDbi)

# Copy results table
res_annot <- res_df

# Add Ensembl ID without version suffix
res_annot$ENSEMBL <- sub("\\..*", "", rownames(res_df))

# Map Ensembl IDs to SYMBOL and GENENAME
anno <- select(org.Hs.eg.db,
               keys = unique(res_annot$ENSEMBL),
               keytype = "ENSEMBL",
               columns = c("SYMBOL","GENENAME"))

# Merge annotation back into result table
res_annot <- merge(res_annot, anno, by="ENSEMBL", all.x=TRUE)

# Check NPM1 result to confirm knockdown efficiency
subset(res_annot, SYMBOL == "NPM1")

# Save annotated DEG table
write.csv(res_annot, "DEG_results_remove_KD2_annotated.csv", row.names=FALSE)

############################################################
# Define a More Permissive DEG Set for Exploratory Analysis
############################################################

sig2 <- subset(res_annot,
               pvalue < 0.05 &
               abs(log2FoldChange) > 0.5)

nrow(sig2)

up2 <- subset(sig2, log2FoldChange > 0.5)
down2 <- subset(sig2, log2FoldChange < -0.5)

nrow(up2)
nrow(down2)

############################################################
# GO Enrichment Analysis
############################################################

library(clusterProfiler)
library(org.Hs.eg.db)

ego_up <- enrichGO(
  gene = up2$SYMBOL,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH"
)

ego_down <- enrichGO(
  gene = down2$SYMBOL,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH"
)

png("GO_up.png", width=1800, height=1400, res=220)
dotplot(ego_up, showCategory=15)
dev.off()

png("GO_down.png", width=1800, height=1400, res=220)
dotplot(ego_down, showCategory=15)
dev.off()

############################################################
# Volcano Plot
############################################################

# Re-define up/down sets if needed
up2 <- subset(sig2, log2FoldChange > 0.5)
down2 <- subset(sig2, log2FoldChange < -0.5)

library(ggplot2)
library(ggrepel)

vol_df <- res_annot

# Remove rows with missing values
vol_df <- subset(vol_df, !is.na(pvalue) & !is.na(log2FoldChange))

# Classify genes using the same threshold as sig2
vol_df$group <- "NS"
vol_df$group[vol_df$pvalue < 0.05 & vol_df$log2FoldChange > 0.5] <- "Up"
vol_df$group[vol_df$pvalue < 0.05 & vol_df$log2FoldChange < -0.5] <- "Down"

# Use raw P value for volcano plot visualization
vol_df$logP <- -log10(vol_df$pvalue)

# Highlight selected genes
label_genes <- c("NPM1")
lab_df <- subset(vol_df, SYMBOL %in% label_genes)

p_vol2 <- ggplot(vol_df, aes(x=log2FoldChange, y=logP, color=group)) +
  geom_point(size=1.8, alpha=0.85) +
  scale_color_manual(values=c("Down"="#4C78A8", "NS"="grey80", "Up"="#D64B4B")) +
  geom_vline(xintercept=c(-0.5, 0.5), linetype=2, linewidth=0.4, color="grey40") +
  geom_hline(yintercept=-log10(0.05), linetype=2, linewidth=0.4, color="grey40") +
  geom_text_repel(
    data=lab_df,
    aes(label=SYMBOL),
    size=4,
    box.padding=0.3,
    point.padding=0.2,
    max.overlaps=50,
    show.legend=FALSE
  ) +
  labs(
    title="NPM1 knockdown vs NC",
    x="log2 Fold Change",
    y=expression(-log[10]("P value"))
  ) +
  theme_classic(base_size=15) +
  theme(
    plot.title = element_text(hjust=0.5, face="bold"),
    legend.title = element_blank(),
    legend.position = "right"
  )

ggsave("Volcano_sig2_101genes.pdf", p_vol2, width=7.5, height=6.2, dpi=300)

############################################################
# GO Bar Plots
############################################################

library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)

go_up_df <- as.data.frame(ego_up)
go_down_df <- as.data.frame(ego_down)

top_up <- go_up_df %>%
  arrange(p.adjust) %>%
  slice(1:5) %>%
  mutate(direction = "Upregulated genes")

top_down <- go_down_df %>%
  arrange(p.adjust) %>%
  slice(1:5) %>%
  mutate(direction = "Downregulated genes")

top_up$minusLog10FDR <- -log10(top_up$p.adjust)
top_down$minusLog10FDR <- -log10(top_down$p.adjust)

top_up$Description <- str_wrap(top_up$Description, width = 28)
top_down$Description <- str_wrap(top_down$Description, width = 28)

top_up$Description <- factor(top_up$Description, levels = rev(top_up$Description))
top_down$Description <- factor(top_down$Description, levels = rev(top_down$Description))

p_go_up <- ggplot(top_up, aes(x = minusLog10FDR, y = Description)) +
  geom_col(fill = "#D64B4B", width = 0.72) +
  theme_classic(base_size = 14) +
  labs(
    title = paste0("Upregulated genes (n=", nrow(up2), ")"),
    x = expression(-log[10]("FDR")),
    y = NULL
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 11)
  )

p_go_down <- ggplot(top_down, aes(x = minusLog10FDR, y = Description)) +
  geom_col(fill = "#4C78A8", width = 0.72) +
  theme_classic(base_size = 14) +
  labs(
    title = paste0("Downregulated genes (n=", nrow(down2), ")"),
    x = expression(-log[10]("FDR")),
    y = NULL
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 11)
  )

ggsave("GO_bar_up.png", p_go_up, width = 7.2, height = 4.8, dpi = 300)
ggsave("GO_bar_down.png", p_go_down, width = 7.2, height = 4.8, dpi = 300)

p_go_combined <- p_go_up / p_go_down
ggsave("GO_bar_combined.pdf", p_go_combined, width = 8, height = 8, dpi = 300)

############################################################
# Hallmark GSEA
############################################################

library(clusterProfiler)
library(enrichplot)
library(msigdbr)
library(dplyr)
library(ggplot2)

# Build ranked gene list using gene symbols
gsea_df <- res_annot %>%
  dplyr::filter(!is.na(SYMBOL), SYMBOL != "", !is.na(log2FoldChange)) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::slice_max(order_by = abs(log2FoldChange), n = 1) %>%
  dplyr::ungroup()

gene_list <- gsea_df$log2FoldChange
names(gene_list) <- gsea_df$SYMBOL
gene_list <- sort(gene_list, decreasing = TRUE)

# Retrieve Hallmark gene sets
msig_h <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- msig_h %>% dplyr::select(gs_name, gene_symbol)

# Run GSEA
gsea_h <- GSEA(
  geneList = gene_list,
  TERM2GENE = hallmark_list,
  pvalueCutoff = 1,
  verbose = FALSE
)

gsea_res <- as.data.frame(gsea_h)
write.csv(gsea_res, "GSEA_Hallmark_results.csv", row.names = FALSE)

# Save overview dot plot
png("GSEA_Hallmark_dotplot.png", width=1800, height=1400, res=220)
dotplot(gsea_h, showCategory = 15, split = ".sign") + facet_grid(.~.sign)
dev.off()

# Inspect GSEA summary table
gsea_res[, c("ID","NES","p.adjust")] %>% arrange(p.adjust)

############################################################
# Hallmark GSEA Bar Plot
############################################################

library(ggplot2)
library(dplyr)
library(stringr)
library(enrichplot)
library(pheatmap)
library(tibble)

gsea_plot_df <- gsea_res %>%
  dplyr::filter(ID %in% c(
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_MYC_TARGETS_V1",
    "HALLMARK_MYC_TARGETS_V2",
    "HALLMARK_MITOTIC_SPINDLE",
    "HALLMARK_DNA_REPAIR",
    "HALLMARK_MTORC1_SIGNALING"
  )) %>%
  mutate(
    ID = gsub("HALLMARK_", "", ID),
    ID = gsub("_", " ", ID),
    ID = factor(ID, levels = rev(ID)),
    direction = ifelse(NES > 0, "Activated", "Suppressed")
  )

p_gsea_bar <- ggplot(gsea_plot_df, aes(x = NES, y = ID, fill = direction)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, linetype = 1, linewidth = 0.4, color = "black") +
  scale_fill_manual(values = c("Activated" = "#D64B4B", "Suppressed" = "#4C78A8")) +
  labs(
    title = "Hallmark pathways altered by NPM1 knockdown",
    x = "Normalized Enrichment Score (NES)",
    y = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    legend.title = element_blank(),
    legend.position = "right"
  )

ggsave("Nature_GSEA_barplot.pdf", p_gsea_bar, width = 5, height = 5.6, dpi = 300)

############################################################
# Hallmark GSEA Enrichment Plots
############################################################

pdf("GSEA_E2F_MYC_G2M_P53.pdf", width = 7, height = 8)

gseaplot2(
  gsea_h,
  geneSetID = "HALLMARK_E2F_TARGETS",
  title = "E2F targets",
  base_size = 10,
  pvalue_table = TRUE,
  ES_geom = "line",
  rel_heights = c(1.3,0.3,0.7)
)

gseaplot2(
  gsea_h,
  geneSetID = "HALLMARK_MYC_TARGETS_V1",
  title = "MYC targets",
  base_size = 10,
  pvalue_table = TRUE,
  ES_geom = "line",
  rel_heights = c(1.3,0.3,0.7)
)

gseaplot2(
  gsea_h,
  geneSetID = "HALLMARK_G2M_CHECKPOINT",
  title = "G2/M checkpoint",
  base_size = 10,
  pvalue_table = TRUE,
  ES_geom = "line",
  rel_heights = c(1.3,0.3,0.7)
)

gseaplot2(
  gsea_h,
  geneSetID = "HALLMARK_P53_PATHWAY",
  title = "P53 pathway",
  base_size = 10,
  pvalue_table = TRUE,
  ES_geom = "line",
  rel_heights = c(1.3,0.3,0.7)
)

dev.off()

############################################################
# Key Gene Heatmap
############################################################

key_genes <- c(
  "NPM1",
  "MYC","E2F1","CCNB1","CDK1","MKI67","PCNA",
  "CDKN1A","GADD45A","ATF3"
)

library(pheatmap)
library(dplyr)

# Extract selected genes from annotated DEG table
key_df <- subset(res_annot, SYMBOL %in% key_genes)
key_ids <- unique(key_df$ENSEMBL)

# Extract expression matrix from vst-transformed data
mat <- assay(vsd2)
mat_df <- data.frame(
  ENSEMBL = rownames(mat),
  mat,
  check.names = FALSE
)

mat_df$ENSEMBL_clean <- sub("\\..*", "", mat_df$ENSEMBL)

# Merge with gene symbols
heat_df <- merge(mat_df, key_df[, c("ENSEMBL","SYMBOL")],
                 by.x = "ENSEMBL_clean", by.y = "ENSEMBL")

# Remove duplicate symbols
heat_df <- heat_df[!duplicated(heat_df$SYMBOL), ]

# Extract sample columns
sample_cols <- rownames(coldata2)
mat2 <- as.matrix(heat_df[, sample_cols])
rownames(mat2) <- heat_df$SYMBOL

# Z-score transformation by gene
mat2 <- t(scale(t(mat2)))

# Column annotation
ann_col <- data.frame(group = coldata2$group)
rownames(ann_col) <- rownames(coldata2)

# Fix sample order
desired_order <- c("NC1","NC2","NC3","KD1","KD3")
desired_order <- desired_order[desired_order %in% colnames(mat2)]
mat2 <- mat2[, desired_order, drop = FALSE]
ann_col <- ann_col[desired_order, , drop = FALSE]

# Fix gene order
gene_order <- key_genes[key_genes %in% rownames(mat2)]
mat2 <- mat2[gene_order, , drop = FALSE]

png("Heatmap_key_genes_reference_style.png", width = 1200, height = 1800, res = 220)

pheatmap(
  mat2,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = ann_col,
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color = "white",
  color = colorRampPalette(c("#3B7DDD", "black", "#D94B4B"))(100),
  fontsize_row = 12,
  fontsize_col = 12,
  cellwidth = 40,
  cellheight = 28,
  annotation_colors = list(
    group = c("NC" = "#4DAF4A", "KD" = "#E41A1C")
  )
)

dev.off()

############################################################
# Exit R and Deactivate Environment
############################################################

q()
conda deactivate
