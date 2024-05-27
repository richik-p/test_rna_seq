# create anaconda env with deseq2 and ggplot2 installed
# conda create -n deseq2 -c bioconda r-deseq2 r-ggplot2 r-essentials r-base
# done 

import gen_bash
import os
import sys
import pandas as pd
import numpy as np
import subprocess

# generate a file in r script to run deseq2
output_directory = gen_bash.output_directory
# subprocess.run(["mkdir", "-p", os.path.join(output_directory, "DESeq2")])

commands = []

commands.append('library("DESeq2")')
commands.append('library("ggplot2")')
commands.append('library("pheatmap")')
commands.append('library("RColorBrewer")')

commands.append('raw_counts <- read.table(os.path.join(output_directory, "raw_counts.txt"), header=TRUE, row.names=1)')

# get number of replicates in each group from file names
condition1 = ""
condition2 = ""
for file in os.listdir(os.path.join(output_directory, "RSEM_counts")).sort():
    if file.endswith('genes.results'):
        if "normal" in file:
            condition1 = "normal"
        else:
            condition2 = "tumor"
            break

# add alpha for padj cutoff by user input
commands.append('alpha <- 0.05') # get alpha from user input

# create metadata files in r to commands list
commands.append('metadata <- data.frame(condition = c(rep ('+condition1+', 3), rep('+condition2+', 3)))')
commands.append('rownames(metadata) <- colnames(raw_counts)')

# create DESeqDataSet object
commands.append('dds <- DESeqDataSetFromMatrix(countData = raw_counts, colData = metadata, design = ~ condition)')

# convert to factors
commands.append('dds$condition <- factor(dds$condition), levels = c('+condition1+', '+condition2+'))')

# run DESeq2
commands.append('dds <- DESeq(dds)')

commands.append('res <- results(dds), name = '+condition2+'_vs_'+condition1+'), alpha = alpha)') # get alpha from user input

# subset out for only significant genes padj < 0.05
commands.append('resSig <- res[which(res$padj < alpha), ]')

# convert to rlog
commands.append('rld <- rlog(dds, blind=FALSE)') # check if blind should be true or false

# get top 50 genes 
number_of_top_genes = 50
commands.append('topGenes <- head(order(rowMeans(counts(dds, normalized=TRUE),decreasing=TRUE), number_of_top_genes)') # TODO variable for number of top genes

# plot heatmap
commands.append('pheatmap(assay(rld)[topGenes,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=dds$condition, annotation_colors=list(condition=c('+condition1+'="blue", '+condition2+'="red"))), filename=os.path.join(output_directory, "DESeq2", "heatmap.pdf"))') # TODO give option to change file type pdf vs png 

# plot MA plot
commands.append('plotMA(res, main="DESeq2", ylim=c(-2,2))')
commands.append('abline(h=c(-2,2), col="blue")')
commands.append('abline(h=0, col="red")')
commands.append('abline(v=0, col="red")')
commands.append('dev.copy(pdf, os.path.join(output_directory, "DESeq2", "MA_plot.pdf"))')
commands.append('dev.off()')

# plot PCA
commands.append('plotPCA(rld, intgroup="condition")')
commands.append('dev.copy(pdf, os.path.join(output_directory, "DESeq2", "PCA_plot.pdf"))')
commands.append('dev.off()')

# LFC shrinkage
commands.append('resLFC <- lfcShrink(dds, coef='+condition2+'_'+condition1+', res=res, type="apeglm")')
commands.append('resLFC <- as.data.frame(resLFC)')
commands.append('write.table(resLFC, file=os.path.join(output_directory, "DESeq2", "LFC_shrinkage.txt"), sep="\t")')


# plot Volcano
commands.append('plot(resLFC$log2FoldChange, -log10(resLFC$padj), pch=20, cex=0.3, col="black")')
commands.append('points(resLFC$log2FoldChange[which(resLFC$padj < alpha)], -log10(resLFC$padj[which(resLFC$padj < alpha)]), pch=20, cex=0.3, col="red")')
commands.append('dev.copy(pdf, os.path.join(output_directory, "DESeq2", "Volcano_plot.pdf"))')
commands.append('dev.off()')

# save results
commands.append('write.table(resLFC, file=os.path.join(output_directory, "DESeq2", "DESeq2_results.txt"), sep="\t")')

# subset out for only significant genes padj < 0.05 and log2FC > 1 and write to file
log2FC_cutoff = 1 # TODO give option to change log2FC cutoff
commands.append('resSig <- resLFC[which(resLFC$padj < alpha & abs(resLFC$log2FoldChange) > log2FC_cutoff), ]') # TODO give option to change log2FC cutoff and padj cutoff
commands.append('write.table(resSig, file=os.path.join(output_directory, "DESeq2", "DESeq2_significant_results.txt"), sep="\t")')

# generate a file in r script to run deseq2
with open(os.path.join(output_directory, "DESeq2", "DESeq2_script.R"), 'w') as f: # TODO change output file location later 
    f.write("\n".join(commands))

# run the script
# subprocess.run(["Rscript", os.path.join(output_directory, "DESeq2", "DESeq2_script.R")])