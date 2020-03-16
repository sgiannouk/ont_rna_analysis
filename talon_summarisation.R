### Stavros Giannoukakos ###
### ONT ANALYSIS / TALON summary analysis 

args <- commandArgs(TRUE)
if (length(args) == 3) {
  # Input matrix as it is outputted from TALON
  matrix <- args[1]
  # Output direcotry where all stats will be saves
  outdir <- args[2]
  # CSV file used for running TALON
  input_groups <- args[3] 
} else {
  cat("ERROR - The number of input arguments is not correct...\nEXITING!\n")
  quit()
}

# outdir <- "/Users/stavris/Desktop/Projects/silvia_ont_umc/batch3_output_for_presentation/talon_analysis"
# matrix <- "/Users/stavris/Desktop/Projects/silvia_ont_umc/batch3_output_for_presentation/talon_analysis/talon_abundance_talon_abundance.tsv"
# input_groups <- "/Users/stavris/Desktop/Projects/silvia_ont_umc/batch3_output_for_presentation/talon_analysis/talon_input.csv"

library("dplyr")
library("plotly")
library("htmltools")
library("webshot")
library("formattable")
library("ggpubr")
theme_set(theme_pubr())

setwd(outdir)

# Input the edited expression matrix
expr_file <- read.delim(matrix)
initial_transcripts <- data.frame(initial_transcripts = nrow(expr_file))
print(paste("Total number of input transcripts:", initial_transcripts, sep=" " ))
# Obtain number of input samples
n <- length(expr_file[ ,12:length(expr_file)])
# Applying basic filterring step 
filter <- 20 #3*n
filtered_table <- expr_file[rowSums(expr_file[ ,12:length(expr_file)]) > filter, ]
filtered_transcripts <- data.frame(filtered_transcripts = nrow(filtered_table))
print(paste("Applying basic filtering step. Transcripts with less counts than", filter, "in all", n, "samples are being discarded.", sep=" " ))
print(paste("Total number of transcripts after filtering:", filtered_transcripts, sep=" " ))

# Obtaining the important annotation to export the stats 
filt_table_overview  <- filtered_table[,3:11]
# Output the overview stats
write.csv(filt_table_overview, file=paste(outdir,"talon_filtered_table.csv",sep="/"), quote=F, row.names=F)
rm(expr_file, matrix)

################## OVERALL ANALYSIS ################## 
# Extracting basic stats
unique_genes <- data.frame(unique_genes = length(unique(filt_table_overview$annot_gene_id)))
unique_transcripts <- data.frame(unique_transcripts = length(unique(filt_table_overview$annot_transcript_id)))
gene_novelty <- count(filt_table_overview[!duplicated(filt_table_overview$annot_gene_id), ], gene_novelty)
transcript_novelty <- count(filt_table_overview, transcript_novelty)

novel_transcripts <- data.frame(novel_transcripts = sum(transcript_novelty$n))

general_table <- t(data.frame(initial_transcripts, 
                            filtered_transcripts, 
                            unique_genes,
                            known_genes = gene_novelty[which(gene_novelty$gene_novelty=="Known"), 2][[1]],
                            novel_genes = gene_novelty[which(gene_novelty$gene_novelty=="Antisense"), 2][[1]] + 
                                          gene_novelty[which(gene_novelty$gene_novelty=="Intergenic"), 2][[1]],
                            unique_transcripts,
                            known_transcripts = transcript_novelty[which(transcript_novelty$transcript_novelty=="Known"), 2][[1]],
                            novel_transcripts = unique_transcripts[[1]]-
                                                transcript_novelty[which(transcript_novelty$transcript_novelty=="Known"), 2][[1]]))

general_table <- data.frame(Id = 1:8, Attributes = row.names(general_table), Remarks = general_table[,1])
row.names(general_table) <- NULL
rm(initial_transcripts, filtered_transcripts, unique_genes, unique_transcripts, novel_transcripts)

# Extracting and Plotting tables containing the above stats
write.csv(general_table, file=paste(outdir, "general_stats_table.csv",sep="/"),quote=F, row.names=F)
write.csv(gene_novelty, file=paste(outdir, "genes_stats_table.csv",sep="/"),quote=F, row.names=F)
write.csv(transcript_novelty, file=paste(outdir, "transcripts_stats_table.csv",sep="/"),quote=F, row.names=F)

# Plotting the general overview 
general_table_sub <- general_table[2:nrow(general_table), 2:length(general_table)]
# Divide each column by total number of transcripts (getting percentages)
general_table_sub$Remarks <- general_table_sub$Remarks/general_table_sub$Remarks[1]

legend_colors <- c("filtered_transcripts"="lavenderblush4", 
                   "known_genes"="#C4961A",
                   "known_transcripts"="#52854C",
                   "novel_genes"="#D16103",
                   "novel_transcripts"="#4E84C4", 
                   "unique_genes"="#C3D7A4",
                   "unique_transcripts"="lavenderblush4")

p1 <- ggplot(general_table_sub, aes(x = Attributes, y = Remarks, fill = Attributes)) + 
      geom_bar(stat = "identity", position = "dodge") +
      geom_text(aes(label=ifelse(round(Remarks, 3) < 1,round(Remarks, 3),"")), position = position_dodge(0.9), vjust = -0.3, size = 3.0, color = "dimgrey") +
      scale_fill_manual("", values = legend_colors) +
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "None") +
      ggtitle("Overall TALON stats")+
      ylab("Percentage of reads (%)") +
      xlab("") 

# Plotting gene/transcript novelty
p2 <- ggplot(transcript_novelty, aes(x = transcript_novelty, y = n, fill = transcript_novelty)) + 
      geom_bar(stat = "identity", position = "dodge") +
      geom_text(aes(label=n), position = position_dodge(0.9), vjust = -0.3, size = 3.0, color = "dimgrey") +
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "None") +
      ggtitle("Overall TALON transcripts novelty")+
      ylab("Reads") +
      xlab("") 

p3 <- ggplot(gene_novelty, aes(x = gene_novelty, y = n, fill = gene_novelty)) + 
      geom_bar(stat = "identity", position = "dodge") +
      geom_text(aes(label=n), position = position_dodge(0.9), vjust = -0.3, size = 3.0, color = "dimgrey") +
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "None") +
      ggtitle("Overall TALON gene novelty")+
      ylab("Reads") +
      xlab("") 

figure <- ggarrange(p1, ggarrange(p3, p2, ncol = 2, labels = c("B", "C")), nrow = 2, labels = "A") 
ggsave(file=paste(outdir, "general_stats_fig.png",sep="/"), width = 10, height = 6, units = "in", dpi = 1200)

# figure2 <- subplot(p1, subplot(p3, p2), nrows = 2 , heights = c(0.4, 0.5), margin = 0.1)
# ptly <- ggplotly(figure2, originalData = T, dynamicTicks = T) %>% layout(yaxis = list(tickformat = "%"))
# htmlwidgets::saveWidget(ptly, paste(outdir, "general_stats_fig.html", sep="/"))

# Outputting the 3 major stats tables 
gene_novelty <- data.frame(Id = 1:3, Gene_novelty = gene_novelty$gene_novelty, Counts = gene_novelty$n)
transcript_novelty <- data.frame(Id = 1:7, Transcript_novelty = transcript_novelty$transcript_novelty, Counts = transcript_novelty$n)

export_formattable <- function(f, file, width = "40%", height = NULL, background = "white", delay = 0.2) {
                               w <- as.htmlwidget(f, width = width, height = height)
                               path <- html_print(w, background = background, viewer = NULL)
                               url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
                               webshot(url,file = file, selector = ".formattable_widget", delay = delay, zoom = 20)}

gtable <- formattable(general_table, list(Remarks = color_bar("cadetblue")))
export_formattable(gtable, paste(outdir, "general_stats_table.png",sep="/"))

gn <- formattable(gene_novelty, list(Counts = color_bar("khaki")))
export_formattable(gn, paste(outdir, "genes_stats_table.png",sep="/"))

tn <- formattable(transcript_novelty, list(Counts = color_bar("lightsteelblue")))
export_formattable(tn, paste(outdir, "transcripts_stats_table.png",sep="/"))
rm(list=setdiff(ls(), c("filtered_table", "n", "outdir", "input_groups")))


################## PER GROUP ANALYSIS ################## 
# Reading input csv file
groups <- read.csv(input_groups, header=F)[ ,1:2]
groups <- groups[order(groups$V1), ]  # Order data frame

samplegroup <- unique(groups$V2) # Obtaining the sample groups


for(i in 1:length(samplegroup)) { 
      
      # Selecting the samples
      selected_samples <- groups$V1[which(groups$V2==samplegroup[i])]
      selected_samples <-as.character(groups$V1[which(groups$V2==samplegroup[i])])
      group_name <- samplegroup[i]  # Group 
      
      # Selecting the matching columns
      selected_filtered_table <- data.frame(filtered_table[1:(length(filtered_table)-n)], filtered_table[ ,selected_samples])
      initial_transcripts <- nrow(selected_filtered_table)
      # Filter 0 rows
      selected_filtered_table <- selected_filtered_table[rowSums(selected_filtered_table[ ,12:(11+length(selected_samples))]) > 0, ]
      filtered_transcripts <- nrow(selected_filtered_table)
      # Output the overview stats
      write.csv(selected_filtered_table, file=paste(outdir, paste("talon", group_name, "selected_filtered_table.csv", sep="_"),sep="/"), quote=F, row.names=F)
      # Extracting stats and plotting
      unique_genes <- data.frame(unique_genes = length(unique(selected_filtered_table$annot_gene_id)))
      unique_transcripts <- data.frame(unique_transcripts = length(unique(selected_filtered_table$annot_transcript_id)))
      gene_novelty <- count(selected_filtered_table[!duplicated(selected_filtered_table$annot_gene_id), ], gene_novelty)
      transcript_novelty <- count(selected_filtered_table, transcript_novelty)
      
      novel_transcripts <- data.frame(novel_transcripts = sum(transcript_novelty$n))
      
      general_table <- t(data.frame(initial_transcripts, 
                                    filtered_transcripts, 
                                    unique_genes,
                                    known_genes = gene_novelty[which(gene_novelty$gene_novelty=="Known"), 2][[1]],
                                    novel_genes = gene_novelty[which(gene_novelty$gene_novelty=="Antisense"), 2][[1]] + 
                                      gene_novelty[which(gene_novelty$gene_novelty=="Intergenic"), 2][[1]],
                                    unique_transcripts,
                                    known_transcripts = transcript_novelty[which(transcript_novelty$transcript_novelty=="Known"), 2][[1]],
                                    novel_transcripts = unique_transcripts[[1]]-
                                      transcript_novelty[which(transcript_novelty$transcript_novelty=="Known"), 2][[1]]))
      
      general_table <- data.frame(Id = 1:8, Attributes = row.names(general_table), Remarks = general_table[,1])
      row.names(general_table) <- NULL
      rm(initial_transcripts, filtered_transcripts, unique_genes, unique_transcripts, novel_transcripts)
      
      # Extracting and Plotting tables containing the above stats
      write.csv(general_table, file=paste(outdir, paste(group_name, "general_stats_table.csv", sep="_"), sep="/"),quote=F, row.names=F)
      write.csv(gene_novelty, file=paste(outdir, paste(group_name, "genes_stats_table.csv", sep="_"), sep="/"),quote=F, row.names=F)
      write.csv(transcript_novelty, file=paste(outdir, paste(group_name, "transcripts_stats_table.csv", sep="_"), sep="/"),quote=F, row.names=F)
      
      # Plotting the general overview 
      general_table_sub <- general_table[2:nrow(general_table), 2:length(general_table)]
      # Divide each column by total number of transcripts (getting percentages)
      general_table_sub$Remarks <- general_table_sub$Remarks/general_table_sub$Remarks[1]
      
      legend_colors <- c("filtered_transcripts"="lavenderblush4", 
                         "known_genes"="#C4961A",
                         "known_transcripts"="#52854C",
                         "novel_genes"="#D16103",
                         "novel_transcripts"="#4E84C4", 
                         "unique_genes"="#C3D7A4",
                         "unique_transcripts"="lavenderblush4")
      
      p1 <- ggplot(general_table_sub, aes(x = Attributes, y = Remarks, fill = Attributes)) + 
        geom_bar(stat = "identity", position = "dodge") +
        geom_text(aes(label=ifelse(round(Remarks, 3) < 1,round(Remarks, 3),"")), position = position_dodge(0.9), vjust = -0.3, size = 3.0, color = "dimgrey") +
        scale_fill_manual("", values = legend_colors) +
        theme_bw() +
        theme(panel.border = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "None") +
        ggtitle(paste(group_name, "Overall TALON stats", sep=" "))+
        ylab("Percentage of reads (%)") +
        xlab("") 
      
      # Plotting gene/transcript novelty
      p2 <- ggplot(transcript_novelty, aes(x = transcript_novelty, y = n, fill = transcript_novelty)) + 
        geom_bar(stat = "identity", position = "dodge") +
        geom_text(aes(label=n), position = position_dodge(0.9), vjust = -0.3, size = 3.0, color = "dimgrey") +
        theme_bw() +
        theme(panel.border = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "None") +
        ggtitle(paste(group_name, "Overall TALON transcripts novelty", sep=" "))+
        ylab("Reads") +
        xlab("") 
      
      p3 <- ggplot(gene_novelty, aes(x = gene_novelty, y = n, fill = gene_novelty)) + 
        geom_bar(stat = "identity", position = "dodge") +
        geom_text(aes(label=n), position = position_dodge(0.9), vjust = -0.3, size = 3.0, color = "dimgrey") +
        theme_bw() +
        theme(panel.border = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "None") +
        ggtitle(paste(group_name, "Overall TALON gene novelty", sep=" "))+
        ylab("Reads") +
        xlab("") 
      
      figure <- ggarrange(p1, ggarrange(p3, p2, ncol = 2, labels = c("B", "C")), nrow = 2, labels = "A") 
      ggsave(file=paste(outdir, paste(group_name, "general_stats_fig.png", sep="_"),sep="/"), width = 10, height = 6, units = "in", dpi = 1200)
      
      # Outputting the 3 major stats tables 
      gene_novelty <- data.frame(Id = 1:3, Gene_novelty = gene_novelty$gene_novelty, Counts = gene_novelty$n)
      transcript_novelty <- data.frame(Id = 1:7, Transcript_novelty = transcript_novelty$transcript_novelty, Counts = transcript_novelty$n)
      
      export_formattable <- function(f, file, width = "40%", height = NULL, background = "white", delay = 0.2) {
        w <- as.htmlwidget(f, width = width, height = height)
        path <- html_print(w, background = background, viewer = NULL)
        url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
        webshot(url,file = file, selector = ".formattable_widget", delay = delay, zoom = 20)}
      
      gtable <- formattable(general_table, list(Remarks = color_bar("cadetblue")))
      export_formattable(gtable, paste(outdir, paste(group_name, "general_stats_table.png", sep="_"), sep="/"))
      
      gn <- formattable(gene_novelty, list(Counts = color_bar("khaki")))
      export_formattable(gn, paste(outdir, paste(group_name, "genes_stats_table.png", sep="_"), sep="/"))
      
      tn <- formattable(transcript_novelty, list(Counts = color_bar("lightsteelblue")))
      export_formattable(tn, paste(outdir, paste(group_name, "transcripts_stats_table.png", sep="_"), sep="/"))
}



