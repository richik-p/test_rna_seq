import os
import argparse

# generate a file in r script to run deseq2
parser = argparse.ArgumentParser(description='RNA-Seq analysis pipeline using STAR for alignment.')
parser.add_argument('--output_dir', required=True, help='Output directory for results')

args = parser.parse_args()
output_directory = args.output_dir

commands = []

'''
library(DESeq2)
library(pheatmap) 
library(RColorBrewer)
'''

commands.append('library("DESeq2")')
commands.append('library("ggplot2")')
commands.append('library("pheatmap")')


'''
# Read files
raw_counts <- read.table("2wk_filtered_counts.tsv", header = TRUE, row.names = 1)
# make sure counts are integers, needed for running DESeq2
countdata_2wk <- round(countdata_2wk)
'''

commands.append('raw_counts <- read.table("'+os.path.join(output_directory, "RSEM_counts", "raw_counts.txt")+'", header=TRUE, row.names=1, sep="\\t")')
commands.append('raw_counts <- round(raw_counts)')


# get number of replicates in each group from file names
condition1 = ""
condition2 = ""
files = os.listdir(os.path.join(output_directory, "RSEM_counts"))
files.sort()
for file in files:
    if file.endswith('genes.results'):
        if "normal" in file:
            condition1 = "normal"
        elif "tumor" in file:
            condition2 = "tumor"

# make metadata file as below
'''
Make metadata like this:
##                     Group Replicate    sampleid
## CK_6wk_1     No_treatment      Rep1    CK_6wk_1
## CK_6wk_2     No_treatment      Rep2    CK_6wk_2
## CK_6wk_3     No_treatment      Rep3    CK_6wk_3
## CKp25_6wk_1 p25_treatment      Rep1 CKp25_6wk_1
## CKp25_6wk_2 p25_treatment      Rep2 CKp25_6wk_2
## CKp25_6wk_3 p25_treatment      Rep3 CKp25_6wk_3
'''

commands.append('metadata <- data.frame(Group = c(rep(\"'+condition1+'\", 3), rep(\"'+condition2+'\", 3)))')
commands.append('rownames(metadata) <- colnames(raw_counts)')
commands.append('metadata$Replicate <- c("Rep1", "Rep2", "Rep3", "Rep1", "Rep2", "Rep3")')
commands.append('metadata$sampleid <- rownames(metadata)')
commands.append('metadata$Group <- factor(metadata$Group)')

'''
ddsMat_2wk <- DESeqDataSetFromMatrix(countData = countdata_2wk,colData = metadata_2wk, design = ~Group)
'''

commands.append('dds <- DESeqDataSetFromMatrix(countData = raw_counts, colData = metadata, design = ~ Group)')

'''
ddsMat_2wk <- DESeq(ddsMat_2wk)
'''

commands.append('dds <- DESeq(dds)')

'''
# Get results from testing with FDR adjust pvalues
res_6wk <- results(ddsMat_6wk, pAdjustMethod = "fdr", alpha = 0.01)
summary(res_6wk)
'''
commands.append('alpha <- 0.05')
commands.append('res <- results(dds, alpha = alpha)')
commands.append('summary(res)')

'''
# Subset for only significant genes (q < 0.01)
res_sig_2wk <- subset(res_2wk, padj < 0.01)

# Convert all samples to rlog
ddsMat_2wk_rlog <- rlog(ddsMat_2wk, blind = FALSE)

# Gather top 30 genes and make matrix
mat <- assay(ddsMat_2wk_rlog[row.names(res_sig_2wk)])[1:30, ]

# Choose which column variables you want to annotate the columns by.
annotation_col = data.frame(
  Group = factor(colData(ddsMat_2wk_rlog)$Group),
  Replicate = factor(colData(ddsMat_2wk_rlog)$Replicate),
  row.names = colData(ddsMat_2wk_rlog)$sampleid
)

# Specify colors you want to annotate the columns by.
ann_colors = list(
  Group = c("No_treatment" = "lightblue", "p25_treatment" = "magenta"),
  Replicate = c(Rep1 = "darkred", Rep2 = "forestgreen", Rep3 = "pink")
)

# Make Heatmap with pheatmap function.
# See more in documentation for customization
pheatmap(mat = mat, 
         color = colorRampPalette(brewer.pal(9, "YlOrBr"))(255), 
         scale = "row", 
         annotation_col = annotation_col, 
         annotation_colors = ann_colors, 
         fontsize = 8, main = "Top DE genes at 2 weeks (rlog transform)",
         show_colnames = F,filename = "Heatmap_2wk.pdf")

'''

commands.append('res_sig <- subset(res, padj < alpha)')
commands.append('rld <- rlog(dds, blind=FALSE)')
commands.append('mat <- assay(rld[row.names(res_sig)])[1:30, ]')
commands.append('annotation_col = data.frame(Group = factor(colData(rld)$Group), row.names = colData(rld)$sampleid)')
commands.append('ann_colors = list(Group = c("'+condition1+'" = "blue", "'+condition2+'" = "red"))')
commands.append('pheatmap(mat = mat, color = colorRampPalette(brewer.pal(9, "YlOrBr"))(255), scale = "row", annotation_col = annotation_col, annotation_colors = ann_colors, fontsize = 8, main = "Top DE genes at 2 weeks (rlog transform)", show_colnames = F, filename = '+os.path.join(output_directory, "deseq", "heatmap.pdf")+')')

'''
# Get results from testing with FDR adjust pvalues for 6wk results
output_res_6wk = as.data.frame(res_6wk)
resordered_6wk = data.frame(output_res_6wk[order(output_res_6wk$padj, na.last=NA),])
print(paste("Number of significant DEGs (6wk):", sum(resordered_6wk$padj < 0.01 & abs(resordered_6wk$log2FoldChange) >= 1, na.rm=TRUE)))
'''

commands.append('output_res <- as.data.frame(res)')
commands.append('resordered <- data.frame(output_res[order(output_res$padj, na.last=NA),])')
commands.append('print(paste("Number of significant DEGs:", sum(resordered$padj < alpha & abs(resordered$log2FoldChange) >= 1, na.rm=TRUE)))')

# generate a file in r script to run deseq2
with open("DESeq2_script.r", 'w') as f: # TODO change output file location later 
    f.write("\n".join(commands))