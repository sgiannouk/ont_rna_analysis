### Stavros Giannoukakos ###
### ONT EXPLORATORY ANALYSIS


args <- commandArgs(TRUE)
if (length(args) == 5) {
  # Input filtered matrix (output of step 8)
  matrix <- args[1]
  # CSV file used for running TALON
  input_groups <- args[2]
  # Output direcotry where all stats will be saves
  main_outdir <- args[3]
  # Minimum gene counts
  minGeneExpr <- args[4]
  # Top N genes to be used for the heatmap
  n_top <- args[5]
} else {
  cat("ERROR - The number of input arguments is not correct...\nEXITING!\n")
  quit()
}

# matrix <- "/Users/stavris/Desktop/Projects/silvia_ont_umc/talon_analysis_reimplementation/filt_talon_abundance.tsv"
# input_groups <- "/Users/stavris/Desktop/Projects/silvia_ont_umc/talon_analysis_reimplementation/talon_input.csv"
# main_outdir <- "/Users/stavris/Desktop/Projects/silvia_ont_umc/talon_analysis_reimplementation/diffExpr_analysis"
# minGeneExpr <- 10  # Minimum gene counts
# n_top <- 50


library("edgeR")
library("dplyr")
library("plotly")
library("DESeq2")
library("ggplot2")
library("DRIMSeq")
library("heatmaply")
library("RColorBrewer")

Sys.setenv("plotly_username"="sgiannouk")
Sys.setenv("plotly_api_key"="MV5szwYMXkYMTiSifC1h")



### EXPLORATORY ANALYSIS ###
print("PERFORMING ON GENE LEVEL EXPLORATORY ANALYSIS")

outdir <- file.path(main_outdir, "ExploratoryAnalysis")
dir.create(outdir, showWarnings = FALSE)
setwd(outdir)

# Input the filtered expression matrix
expr_file <- read.csv(matrix)

# Reading input csv file containing the groups
groups <- read.csv(input_groups, header=F)[ ,1:3]
groups <- groups[order(groups$V2), ] # Order data frame

sampletypevalues <- factor(unique(groups$V2)) # Obtaining the sample groups

matfile <- expr_file[ ,c(1,10:length(expr_file))]  # Obtaining only the annot_gene_id and the counts

# Summarising all reads per gene name and removing the duplicate  rows
matfile <- data.frame(dplyr::group_by(matfile, annot_gene_id) %>% dplyr::summarise_all(sum))

# Making the annot_transcript_id as row names and removing it as from 1st column
row.names(matfile) <- matfile$annot_gene_id; matfile$annot_gene_id <- NULL

# Applying a basic filtering step, where genes with less 
# than readCountMinThreshold will be excluded from the analysis.
# This step is recommeneded by ONT
matfile <- matfile[rowSums(matfile) > minGeneExpr, ]

# Designing the data's factors which indicate the experimental group for each sample
samplefactors <- data.frame(row.names=groups$V1, condition = factor(groups$V2, levels=sampletypevalues), batch = factor(groups$V3))

# Constructing the DESeqDataSet object which is the starting point of the analysis
if (length(unique(samplefactors$batch) == 1)) {
  dds = DESeqDataSetFromMatrix(countData = matfile, colData = samplefactors, design = ~ condition)
} else {
  dds = DESeqDataSetFromMatrix(countData = matfile, colData = samplefactors, design = ~ batch + condition) }
rm(samplefactors)

## PCA PLOT
# DESeq2 offers transformations for count data that stabilize the variance across the mean:
# the regularize logarithm (rlog) and the variance stabilizing transformation (VST).
# Constructing the DESeqDataSet object which is the starting point of the analysis
vsd <- vst(dds, blind=T)
data_pca <- plotPCA(vsd, intgroup = "condition", returnData=TRUE)
percentVar <- round(100 * attr(data_pca, "percentVar"))
pca <- ggplot(data_pca, aes(PC1, PC2, color=condition, label=name)) + geom_point(size=2) +
       scale_color_manual(name="Groups", values=c("#66CC99", "#877598")) +
       theme_bw() +
       theme(panel.border = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
       theme(legend.position = "bottom", legend.justification = "center") +
       ggtitle("Principal Component Analysis")+
       xlab(paste0("PC1: ",percentVar[1],"% variance")) +
       ylab(paste0("PC2: ",percentVar[2],"% variance"))
ggsave(file=paste(outdir, "/explAnalysis_PCAPlot.png",sep=""), width = 10, height = 6, units = "in", dpi = 1200)
ptly <- ggplotly(pca, originalData = T, dynamicTicks = T) #%>% layout(yaxis = list(tickformat = "%"))
htmlwidgets::saveWidget(ptly, paste(outdir, "/explAnalysis_PCAPlot.html",sep=""))
rm(dds, vsd, data_pca, percentVar, pca, ptly)

check.integer <- function(N){ !grepl("[^[:digit:]]", format(N,  digits = 20, scientific = FALSE)) }

data <- DGEList(counts=matfile, group=groups$V2)  # Summarise the input data


# Colouring the different conditions
if (length(sampletypevalues) == 2) {
  col_condition <- c("#7FC97F", "#BEAED4")[data$samples$group]
} else {
  col_condition <- c(brewer.pal(n = length(sampletypevalues), name = "Accent"))[data$samples$group]
}

# Examine the distributions of the raw counts by plotting the log2CPM of the counts
print("Checking the distribution of the read counts on the log2 scale...")
png(paste(outdir,"/explAnalysis_Log2DistPlot.png",sep=""), units='px', height=900, width=1600, res=90)
# Check distributions of samples using boxplots
par(mar=c(8.1, 4.1, 4.1, 2.1))
boxplot(cpm(data$counts, prior.count=2, log=TRUE),col=col_condition, xlab="", ylab="Log2 counts per million", las=2)
# Adding a blue horizontal line that corresponds to the median log2CPM
abline(h=median(cpm(data$counts, prior.count=2, log=TRUE)), col="slategrey", lwd=2)
title("Boxplots of log2CPMs (unnormalised)")
dev.off()

# An MDSplot is a visualisation of a principle components analysis, which determines the greatest sources of variation in the data.
# If the experiment is well controlled and has worked well, what we hope to see is that the greatest sources of variation are the
# treatments/groups we are interested in. It is also an incredibly useful tool for quality control and checking for outliers.
print("Creating the MultiDimensional Scaling plot...")
png(paste(outdir,"/explAnalysis_plotMDS.png",sep=""), units='px', height=900, width=1600, res=90)
plotMDS(data, col=col_condition)
title(main = "MultiDimensional Scaling plot\n(distances approximate the log2 fold changes between the samples)")
dev.off()
rm(data)

# Calculating counts per million 
myCPM <- cpm(matfile)
# Discarding lowly expressed genes
thresh <- myCPM > 0.5
keep <- rowSums(thresh) >= 2
counts.keep <- matfile[keep, ]
# Importing data to DESeq2 data matrix
y <- DGEList(counts.keep)
# Get log2 counts per million
logcounts <- cpm(y,prior.count = 0.1, log=TRUE)
var_genes <- apply(logcounts, 1, var)
# Selecting the most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))
highly_variable_lcpm <- logcounts[select_var, ]


if (check.integer(n_top)) {
  # Obtaining the n top genes
  top_genes <- highly_variable_lcpm[1:n_top, ]
} else {
  # Obtaining the percentage of top genes
  top_genes <- highly_variable_lcpm[1:(nrow(highly_variable_lcpm)*n_top), ]
}

colnames(groups) <- c("Samples","Groups","Batch")
groups$Batch <- NULL
# Heatmap title
heatmap_title = paste("Top",n_top,"most variable genes across samples",sep=" ")

# Creating the heatmaps. If samples are more than 100, then sample labels are being omitted
if (length(groups) <= 100) {
  # Heatmap of top selected normalised genes
  heatmaply(top_genes, file = paste(outdir,"/explAnalysis_heatmap.html",sep =""),
            limits = NULL, colors = brewer.pal(11,"Spectral"), scale = "row", main = heatmap_title,
            key.title=NULL, col_side_colors = data.frame(groups), hide_colorbar = FALSE,
            column_text_angle=60, fontsize_col = 9, fontsize_row = 8,  showticklabels=c(TRUE,TRUE))
  heat_map <- heatmaply(top_genes, limits = NULL, colors = brewer.pal(11,"Spectral"), scale = "row",
                        main = heatmap_title, key.title=NULL, col_side_colors = data.frame(groups),
                        hide_colorbar = FALSE, column_text_angle=60, fontsize_col = 9, fontsize_row = 8,
                        showticklabels=c(TRUE,TRUE))
} else {
  # Heatmap of top selected genes, No sample-names
  heatmaply(top_genes, file = paste(outdir,"/explAnalysis_heatmap.html",sep =heatmap_title),
            limits = NULL, colors = brewer.pal(11, "Spectral"), scale = "row", main = "Top",
            key.title=NULL, col_side_colors = data.frame(groups), hide_colorbar = FALSE,
            showticklabels = c(FALSE, TRUE), fontsize_row = 8)
  heat_map <- heatmaply(top_genes, limits = NULL, colors = brewer.pal(11, "Spectral"), scale = "row",
                        main = heatmap_title, key.title=NULL, col_side_colors = data.frame(groups),
                        hide_colorbar = FALSE, fontsize_row = 8, showticklabels = c(FALSE, TRUE))

}

plotly_IMAGE(heat_map, width = 1200, height = 800, format = "png", out_file = paste(outdir,"/explAnalysis_heatmap.png",sep =""))
rm(counts.keep, heat_map, highly_variable_lcpm, logcounts, myCPM, thresh, top_genes, y, col_condition, heatmap_title, keep, 
   n_top, select_var, var_genes, check.integer)


### Checking number of  transcripts per condition
# Obtaining the counts table
matfile <- expr_file[ ,c(2,1,10:length(expr_file))]  # Obtaining only the annot_gene_id and the counts
colnames(matfile)[1:2] <- c("feature_id","gene_id")

group_samples <- data.frame(sample_id = groups$Samples, condition = groups$Groups)

# Create a dmDSdata object that is the starting point of DRIMSeq
drimseq <- dmDSdata(counts = matfile, samples = group_samples)

# Filtering of lowly expressed genes
# The parameters are the suggested by ONT
drimseq <- dmFilter(drimseq, 
                    #min_samps_gene_expr = minSampsGeneExpr, 
                    min_gene_expr = minGeneExpr)

# Obtaining the counts after dmFiltering
filtered_counts <- counts(drimseq)

# Selecting the samples from the first group
selected_group1 <- cbind(filtered_counts[,c(1,2)], filtered_counts[ ,which(colnames(filtered_counts) %in% group_samples$sample_id[group_samples$condition==sampletypevalues[1]])])
# Removing unexpressed genes in this group (genes with 0 in all samples)
selected_group1 <- selected_group1[rowSums(selected_group1[ ,which(colnames(selected_group1) %in% group_samples$sample_id[group_samples$condition==sampletypevalues[1]])]) > 0, ] 
# Transform to count table
selected_group1 <- data.frame(table(selected_group1$gene_id))
# Collapsing and counting transcripts
selected_group1 <- data.frame(table(selected_group1$Freq))
selected_group1$group <- as.character(sampletypevalues[1])  # Renaming all genes to samplegroup

# Selecting the samples from the second group
selected_group2 <- cbind(filtered_counts[,c(1,2)], filtered_counts[ ,which(colnames(filtered_counts) %in% group_samples$sample_id[group_samples$condition==sampletypevalues[2]])])
# Removing unexpressed genes in this group (genes with 0 in all samples)
selected_group2 <- selected_group2[rowSums(selected_group2[ ,which(colnames(selected_group2) %in% group_samples$sample_id[group_samples$condition==sampletypevalues[2]])]) > 0, ]
# Transform to count table
selected_group2 <- data.frame(table(selected_group2$gene_id))
# Collapsing and counting transcripts
selected_group2 <- data.frame(table(selected_group2$Freq))
selected_group2$group <- as.character(sampletypevalues[2])  # Renaming all genes to samplegroup

# Merging the melted data frames
merged_melted_tables <- rbind(selected_group1, selected_group2)
rm(selected_group1, selected_group2)

ggplot(merged_melted_tables, aes(x = Var1, y = Freq, fill = group)) + 
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual("Groups", values = c("#66CC99", "#877598")) +
  theme_bw() +
  theme(legend.position="bottom") +
  theme(panel.border = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Number of transcript usage per condition") +
  xlab("Number of transcprits") +
  ylab("Frequency")
ggsave(file=paste(outdir, "/explAnalysis_TrNumPerCondition.png",sep=""), width = 10, height = 6, units = "in", dpi = 1200)
