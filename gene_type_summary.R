### Stavros Giannoukakos ###
### ONT ANALYSIS gene type analysis 

args <- commandArgs(TRUE)
if (length(args) == 2) {
  current_dir = args[1]
  outdir <- args[2]
} else {
  cat("ERROR - The number of input arguments is not correct...\nEXITING!\n")
  quit()
}

library("plotly")
library("ggplot2")
library("reshape2")
library("RColorBrewer")
setwd(outdir)

# Input the edited expression matrix
expr_file <- read.csv(paste(current_dir, "perGene_expression_matrix.csv", sep="/"), header=TRUE, row.names=1)
# Convert sample columns into numeric
expr_file[ ,2:ncol(expr_file)] <- sapply(expr_file[ ,2:ncol(expr_file)], as.numeric)

# Divide each column by its sum (getting percentages)
expr_file[ ,2:ncol(expr_file)] <- as.data.frame(lapply(expr_file[ ,2:ncol(expr_file)], function(x) x/sum(x)))


# Prepare the matrix for plotting
mtx <- melt(expr_file, id.vars="gene_type")

colors <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA","#E6AB02", "#A6761D", 
            "#1B9E77", "#D7C1B1", "#689030", "#652926", "#66A61E", "#56B4E9", 
            "#5E738F", "#D1A33D", "#3F4921", "#7FDCC0", "#C84248")
num <- length(unique(expr_file$gene_type))

# ggplot 
p <- ggplot(mtx, aes(x = variable , y = value, fill = gene_type)) +
     geom_bar(stat = "identity", position = "stack") +
     theme_bw() +
     scale_fill_manual("", values = colorRampPalette(colors)(num)) +
     scale_y_continuous(labels = scales::percent) +
     theme(legend.position = "bottom") +
     ylab("Percentage of reads (%)") +
     xlab("") +
     ggtitle("Gene type summarisation")
ggsave(file="gene_type_summarisation.png", width = 10, height = 6, units = "in", dpi = 1200)

ptly <- ggplotly(p, originalData = T, dynamicTicks = T) %>% 
        layout(yaxis = list(tickformat = "%"))

htmlwidgets::saveWidget(ptly, "gene_type_summarisation.html")
