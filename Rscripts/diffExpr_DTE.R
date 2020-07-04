### Stavros Giannoukakos ###
### ONT ANALYSIS FOR DIFFERENTIAL TRANSCRIPT EXPRESSION (DTE)


args <- commandArgs(TRUE)
if (length(args) == 7) {
  # Input filtered matrix (output of step 8)
  matrix <- args[1]
  # CSV file used for running TALON
  input_groups <- args[2]
  # Output direcotry where all stats will be saves
  main_outdir <- args[3]
  # Transcripts expressed in minimum this many samples
  minSampsFeatureExpr <- args[4]
  # Minimum transcript counts
  minFeatureExpr <- args[5]
  # Adjusted p-value threshold for differential expression analysis
  adjPValueThreshold <- args[6]
  # Minimum required log2 fold change for differential expression analysis
  lfcThreshold <- args[7]
} else {
  cat("ERROR - The number of input arguments is not correct...\nEXITING!\n")
  quit()
}

# matrix <- "/Users/stavris/Desktop/Projects/silvia_ont_umc/talon_analysis_reimplementation/filt_talon_abundance.tsv"
# input_groups <- "/Users/stavris/Desktop/Projects/silvia_ont_umc/talon_analysis_reimplementation/talon_input.csv"
# main_outdir <- "/Users/stavris/Desktop/Projects/silvia_ont_umc/talon_analysis_reimplementation/diffExpr_analysis"
# minSampsFeatureExpr <- 3  # Transcripts expressed in minimum this many samples
# minFeatureExpr <- 10  # Minimum transcript counts
# adjPValueThreshold <- 0.05
# lfcThreshold <- 1


library("dplyr")
library("edgeR")
library("DRIMSeq")
library("ggplot2")
library("reshape")



##### DIFFERENTIAL TRANSCRIPT EXPRESSION (DTE) ANALYSIS USING DRIMSEQ/EDGER #####
print("RUNNING DIFFERENTIAL TRANSCRIPT EXPRESSION (DTE) ANALYSIS USING DRIMSEQ/EDGER")

outdir <- file.path(main_outdir, "DiffTranscriptExpr")
dir.create(outdir, showWarnings = FALSE)
setwd(outdir)

# Input the filtered expression matrix
expr_file <- read.csv(matrix)
print(paste("Total number of unique trancripts:", length(unique(expr_file$annot_transcript_id)), sep=" " ))

# Obtain number of input samples
n <- length(expr_file[ ,10:length(expr_file)])
print(paste("Total number of input samples:", n, sep=" " ))
# Reading input csv file containing the groups
group_samples <- read.csv(input_groups, header=F)[ ,1:3]
group_samples <- group_samples[order(group_samples$V2), ] # Order data frame
# Preparing a dataframe containing information about the samples 
group_samples <- data.frame(sample_id = group_samples$V1, condition = group_samples$V2, batch=group_samples$V3)

sampletypevalues <- factor(unique(group_samples$condition)) # Obtaining the sample groups
print(paste("Total number of input groups: ", 
            length(sampletypevalues)," (", sampletypevalues[1], " and ", sampletypevalues[2],")", sep="" ))


# Obtaining the counts table
matfile <- expr_file[ ,c(2,1,10:length(expr_file))]  # Obtaining only the annot_transcript_id and the counts
colnames(matfile)[1:2] <- c("feature_id","gene_id")

# Create a dmDSdata object that is the starting point of DRIMSeq
drimseq <- dmDSdata(counts = matfile, samples = group_samples)
print("General input stats:")
print(drimseq)   # General stats

# Check what is the minimal number of replicates per condition
print("Number of replicates per condition")
print(table((DRIMSeq::samples(drimseq))$condition))

# Filtering of lowly expressed transcript
# The parameters are the suggested by ONT
drimseq <- dmFilter(drimseq, 
                    min_samps_feature_expr = minSampsFeatureExpr,
                    min_feature_expr = minFeatureExpr)


# Obtaining the counts after dmFiltering
filtered_counts <- counts(drimseq)

# Printing several stats
out = nrow(matfile)  - length(drimseq)
tot = nrow(matfile)
perc = round((out/tot) * 100, 1)
print(paste("In total ", out,
            " out of ", tot,  
            " (", perc, "%)"," transcripts did not pass the minimum thersholds", sep=""))
print(paste("We are continuing the analysis with ", length(unique(filtered_counts$feature_id))," transcripts",sep=""))

# Creating with design matrix
if (length(unique(group_samples$batch) == 1)) {
  design <- model.matrix( ~ condition, data = DRIMSeq::samples(drimseq))
} else {
  design <- model.matrix( ~ condition + batch, data = DRIMSeq::samples(drimseq)) }
rm(n, out, tot, perc, minSampsFeatureExpr, minFeatureExpr, matfile)



# Removing gene_id column
transcript_filt_counts <- filtered_counts; transcript_filt_counts$gene_id <- NULL

# Making the annot_transcript_id as row names and removing it as from 1st column
row.names(transcript_filt_counts) <- transcript_filt_counts$feature_id; transcript_filt_counts$feature_id <- NULL

# Storing the raw read counts table in a simple list-based data object called a DGEList.
edgeR_table <- DGEList(transcript_filt_counts, group=group_samples$condition)

# Normalisation for RNA composition by finding a set of scaling factors for the library sizes 
# that minimize the log-fold changes between the samples for most transcripts. The default method
# for computing these scale factors uses a trimmed mean of M-values (TMM) between each pair of samples.
edgeR_table <- calcNormFactors(edgeR_table)

# Estimating common dispersion and tagwise dispersions 
edgeR_table <- estimateDisp(edgeR_table, design)

# Plotting the per-transcript dispersion estimates
png(paste(outdir,"/edgeR_DispersionPlot.png",sep=""), units='px', height=900, width=1600, res=90)
plotBCV(edgeR_table)
title(main = "Dispersion Estimates")
dev.off()

# Performing quasi-likelihood F-tests
fit <- glmQLFit(edgeR_table, design)

# Plotting the quasi-likelihood dispersion  estimates
png(paste(outdir,"/edgeR_QLDispersionPlot.png",sep=""), units='px', height=900, width=1600, res=90)
plotQLDisp(fit)
title(main = "QL Dispersion Estimates")
dev.off()

qlf <- glmQLFTest(fit)

# Obtaining the results table and sorting by p-value
edger_res_overall <- topTags(qlf, n=Inf, sort.by="PValue")[[1]]
edger_res <- edger_res_overall[,c(1,4,5)];  colnames(edger_res) <- c("log2FoldChange",  "pval", "padj")

# Output basic stats
print(paste("Results indicating the transcript set with log2FoldChange value > |", lfcThreshold,"| and with adjasted p-value < ", adjPValueThreshold, sep=""))
is.de <- decideTests(qlf, p.value=adjPValueThreshold, lfc=lfcThreshold)
print(summary(is.de))


# Selecting the samples
selected_samples <- colnames(transcript_filt_counts[ ,(which(group_samples$condition==sampletypevalues[1] | group_samples$condition==sampletypevalues[2]))])
# Obtaining the list of features with absolute log2FoldChange greater than 1 and p adjusted value lower than the user-input value
results_selected <- edger_res[(abs(edger_res$log2FoldChange)>=lfcThreshold) & (edger_res$padj<=adjPValueThreshold), ]
# Combing the normalised data along with statistical analysis results ("log2FoldChange", "pval", "padj")
results_selected <- data.frame(merge(cpm(edgeR_table)[ ,selected_samples], results_selected, by=0))
names(results_selected)[1] <- "transcript_id"

# Manupulation of expr_file in order to 
# merge add the transcript names in the final table
expr_file <- expr_file[ ,c(2,4)]
expr_file <- expr_file[row.names(unique(expr_file[ ,c(1,2)])), ]
row.names(expr_file) <- expr_file$annot_transcript_id; expr_file$annot_transcript_id <- NULL

# Adding transcript_names in the final table
results_selected <- data.frame(merge(expr_file, results_selected, by.x=0, by.y="transcript_id"))
# Renaming the first two columns
colnames(results_selected)[1] <- "transcript_id"; colnames(results_selected)[2] <- "transcript_name"
# Sorting the table by adjasted p-value
results_selected <- results_selected[order(results_selected$pval), ]

# Exporting the normalised results table containing the selected features with log2FoldChange greater than 1 and adjusted p value lower than the user-input value
write.table(results_selected, file=paste(outdir,"/",sampletypevalues[1],"VS",sampletypevalues[2],"_edgeR_topTranscriptsBelow", gsub("[.]", "", adjPValueThreshold), "LFC", gsub("[.]", "", lfcThreshold), ".csv", sep=""), sep="\t", row.names = F, quote=FALSE)

# Creating the overall table with cpm values and stats from EdgeR
total_results <- data.frame(merge(cpm(edgeR_table)[ ,selected_samples], edger_res, by=0))
names(total_results)[1] <- "transcript_id"
# Adding transcript_names in the final table
total_results <- data.frame(merge(expr_file, total_results, by.x=0, by.y="transcript_id"))
colnames(total_results)[1] <- "transcript_id"; colnames(total_results)[2] <- "transcript_name"
# Sorting the table by adjasted p-value
total_results <- total_results[order(total_results$padj), ]

# Exporting the normalised results table containing the selected features with log2FoldChange greater than 1 and adjusted p value lower than the user-input value
write.table(total_results, file=paste(outdir,"/",sampletypevalues[1],"VS",sampletypevalues[2],"_edgeR_allTranscripts.csv", sep=""), sep="\t", row.names = F, quote=FALSE)


# Plotting log-fold change against log-counts per million, with DE transcripts highlighted
png(paste(outdir,"/edgeR_MDPlot.png",sep=""), units='px', height=900, width=1600, res=90)
plotMD(qlf)
abline(h=c(-lfcThreshold, lfcThreshold), col="blue")
dev.off()


# MA-plot - Drawing the expression levels over the exons to highlight differential exon usage
MAPlotData <- qlf$table
logUp <- which(MAPlotData$logFC >= lfcThreshold)
logDown <- which(MAPlotData$logFC <= -lfcThreshold)
withStat <- which(MAPlotData$PValue <= adjPValueThreshold)
colours <- c(noDifference="dimgray", upRegulated="indianred3", downRegulated="mediumseagreen")
gene <- rep("noDifference", nrow(MAPlotData))
gene[logUp[logUp %in% withStat]] <- "upRegulated"
gene[logDown[logDown %in% withStat]] <- "downRegulated"
ggplot(data.frame(MAPlotData), aes(y=logFC, x=logCPM)) + 
        geom_point(size=1.2) + 
        geom_hline(yintercept = -lfcThreshold, color="mediumseagreen") + 
        geom_hline(yintercept = lfcThreshold, color="indianred3") +
        theme_bw() +
        aes(colour=gene) + 
        scale_colour_manual(name="Transcripts", values=colours) +
        xlab("log(CountsPerMillion)") +
        ylab("log(FoldChange)") +
        ggtitle("MA plot - log(FC) vs. log(CPM) on transcript level data")
ggsave(file=paste(outdir,"/edgeR_MAplot.png",sep=""), width = 10, height = 6, units = "in", dpi = 1200)
rm(MAPlotData, logUp, logDown, withStat, colours, gene)


# Volcano plot
# Compute significance, with a maximum of 350 for the p-values set to 0 due to limitation of computation precision
edger_res_overall <- merge(expr_file, edger_res_overall, by = 0)
row.names(edger_res_overall) <- edger_res_overall$Row.names; edger_res_overall$Row.names <- NULL; colnames(edger_res_overall)[1] <- "transcript_name"
edger_res_overall$significant <- (-log10(edger_res_overall$FDR))
edger_res_overall[is.infinite(edger_res_overall$significant), "significant"] <- 350
# Select transcripts with a defined p-value (DESeq assigns NA to some transcripts)
genes.to.plot <- !is.na(edger_res_overall$PValue)
## Volcano plot of adjusted p-values
cols <- densCols(edger_res_overall$logFC, edger_res_overall$significant)
cols[edger_res_overall$PValue==0] <- "purple"
edger_res_overall$pch <- 19
edger_res_overall$pch[edger_res_overall$PValue == 0] <- 6
print("Generating the volcano plot...")
png(paste(outdir,"/edgeR_volcanoPlot.png",sep=""), units='px', height=900, width=1600, res=100)
plot(edger_res_overall$logFC,
     edger_res_overall$significant,
     col=cols,
     panel.first=grid(),
     main="Volcano plot",
     xlab="Effect size: log2(Fold-Change)",
     ylab="-log10(adjusted p-value)",
     pch=edger_res_overall$pch, cex=0.4)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(adjPValueThreshold), col="brown")
## Plot the names of a reasonable number of transcripts, by selecting those begin not only significant but also having a strong effect size
gn.selected <- (abs(edger_res_overall$logFC)>=lfcThreshold & edger_res_overall$FDR<=adjPValueThreshold)
text(edger_res_overall$logFC[gn.selected], edger_res_overall$significant[gn.selected],lab=edger_res_overall$transcript_name[gn.selected], cex=0.6)
dev.off()
rm(genes.to.plot, cols, gn.selected)
edger_res_overall$significant <- NULL; edger_res_overall$pch  <- NULL; edger_res_overall$transcript_name <- NULL
rm(is.de,  fit, qlf, results_selected, transcript_filt_counts, edger_res)


### TOP 30 DE TRANSCRIPTS
# Plotting the normalized count values for the top 30 differentially expressed transcripts (by padj values)
## Order results by padj values and get the first 30 transcripts
edger_res_overall <- edger_res_overall[order(edger_res_overall$PValue), ]
top30_sigGenes <- row.names(edger_res_overall)[1:30]
## Normalized counts for top 20 significant transcripts
top30_sigNorm <- data.frame(cpm(edgeR_table)) %>% filter(rownames(cpm(edgeR_table)) %in% top30_sigGenes)
top30_sigNorm <- data.frame(merge(expr_file, top30_sigNorm, by=0)); top30_sigNorm$Row.names <- NULL
colnames(top30_sigNorm)[1] <- "transcript_name"
# Gathering the columns to have normalized counts to a single column
melt_top30_sigNorm <- melt(top30_sigNorm)
melt_top30_sigNorm$group <- gsub("_[^_]+$", "", melt_top30_sigNorm$variable)
melt_top30_sigNorm$value[melt_top30_sigNorm$value == 0] <- 0.1

ggplot(melt_top30_sigNorm, aes(x = transcript_name, y = value, color = group)) +
  geom_point(size=2.5) +
  scale_y_log10(labels = function(x) format(x, scientific = F)) +
  xlab("Transcripts") +
  ylab("log10(CPM)") +
  ggtitle("Top 30 Significant DE Transcripts") +
  theme_bw() +
  scale_colour_manual(name="", values=c("#66CC99", "#877598")) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="bottom")
ggsave(file=paste(outdir,"/edgeR_top30MostSingificantTranscripts.png",sep=""), width = 10, height = 6, units = "in", dpi = 1200)
rm(top30_sigGenes, melt_top30_sigNorm, top30_sigNorm)

