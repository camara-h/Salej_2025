# List of packages you want to ensure are installed
packages_needed <- c("readxl", "tidyverse", "ggpubr")

# Find which packages are not installed yet
packages_to_install <- packages_needed[!(packages_needed %in% installed.packages()[, "Package"])]

# Install the missing packages
if (length(packages_to_install) > 0) {
  install.packages(packages_to_install)
}

knitr::opts_chunk$set(echo = TRUE)
library(readxl)
library(tidyverse)
library(ggpubr)

# Formatting the metadata
meta <- read_excel("MassSpectrometryAnalysis TSH-BAT.xlsx", sheet = 2, col_names = T)
colnames(meta) <- meta[19, ]
meta <- meta[20:35, ]
head(meta)

df <- rbind(meta[, 1:6], meta[, 7:12])
df <- cbind(df, meta[, 13:ncol(meta)])
names(df)[1] <- "Sample"
df$subject <- gsub(".* (.*)", "\\1", df$Sample)
df$tissue <- gsub("(.*) (.*)", "\\1", df$Sample)
df <- df %>% select(-"NA")

names(df) <- gsub(" ", "_", names(df))
names(df) <- gsub("/", ".per.", names(df))
names(df) <- gsub("\\[|\\]", "", names(df))

# Formatting the metadata
meta <- read_excel("Metadata_for_correlation.xlsx", col_names = T)

meta <- meta %>% separate_wider_delim(cols = Sample, delim = " ", names = c("Tissue", "Donor_ID"))
meta <- meta %>% rename(
  "Serum_Free_T3" = "Free T3",
  "Serum_Free_T4" = "Free T4",
  "Serum_TSH" = "TSH",
  "tissue" = "Tissue",
  "subject" = "Donor_ID"
)

names(meta)
head(meta)
names(meta)

names(meta) <- gsub("[ /]", "_", names(meta))
plot.meta <- meta %>% select(-`log[TSH]`)
plot.meta$Serum_Free_T3 <- as.numeric(plot.meta$Serum_Free_T3)


correct_mg_string <- function(x) {
  # First, replace '.per.' with '/'
  x <- gsub(".per.", "/", x)

  # First, replace '_mg' with '/mg'
  x <- gsub("_mg", "/mg", x)

  # Next, replace '_μg' with '/μg'
  x <- gsub("_μg", "/μg", x)

  # Next, replace 'T3_' with 'T3/'
  x <- gsub("T3_", "T3/", x)

  # Next, replace 'T4_' with 'T4/'
  x <- gsub("T4_", "T4/", x)

  # Finally, replace any remaining '_' with a space
  x <- gsub("_", " ", x)

  return(x)
}

# Pseudocode:
# This part of the script needs to get as input: 1) a gene or vector of genes; 2) The CPM expression matrix + Metadata information;
# The output will be a set of correlation graphs for each gene to the variable of interest. The correlation needs to be colored by depot and shaped by subject;
# In the graphs we should have the R2 and p-value. The graphs can be plotted using ggscatter (from ggpubr).
#
# Step 1: Load CPM matrix
# Step 2: Merge CPM matrix and Metadata by sample (Rows = Samples, Cols = Genes + Metadata)
# Step 3: Get the list of genes (This sould be the only step with user interaction)
# Step 4: Plot the correlation from selected genes with all metadata variables


# Assuming the CPM matrix is stored as a CSV file
cpm_matrix <- read.csv("./logCPM_for_plotting.csv", row.names = 1)

cpm_matrix$gene[cpm_matrix$gene == ""] <- rownames(cpm_matrix)[cpm_matrix$gene == ""]
gene_names <- cpm_matrix %>% select(gene)
t_cpm_matrix <- t(cpm_matrix[, -1])

# Load metadata, assuming it is also stored in a CSV file
f.df <- df %>%
  filter(!is.na(Sample)) %>%
  column_to_rownames("Sample")
rownames(f.df) <- gsub("Deep ", "Deep", rownames(f.df))
rownames(f.df) <- gsub("SubQ ", "SC", rownames(f.df))
data_merged <- merge(t_cpm_matrix, f.df, by = "row.names", all = TRUE)
data_merged <- data_merged %>% column_to_rownames("Row.names")

# Create the rownames of meta
meta <- meta %>%
  mutate(Rows = paste0(tissue, subject)) %>%
  mutate(Rows = gsub("SubQ", "SC", Rows)) %>%
  column_to_rownames("Rows")
data_merged <- merge(t_cpm_matrix, meta, by = "row.names", all = TRUE)



### START HERE###

# PLOTTING THE DATA - Start inputting your genes of interest here

PlotGeneCorrelation <- function(genes_of_interest = "UCP1", FDR_threshold = 0.15) {
  # Extract statistics from BMI_temp_covar
  gene_stats <- read.csv("./BMI_temp_covars/gene_stats.csv")
  # Extract JBC information and add to the plot

  names(meta)
  names(gene_stats)

  # Format the names
  names(gene_stats) <- gsub("T3_per_mg_tissue", "T3_mg_tissue", names(gene_stats))
  names(gene_stats) <- gsub("T4_per_mg_tissue", "T4_mg_tissue", names(gene_stats))
  names(gene_stats) <- gsub("T3_per_mcg_protein", "T3_μg_protein", names(gene_stats))
  names(gene_stats) <- gsub("T4_per_mcg_protein", "T4_μg_protein", names(gene_stats))
  names(gene_stats) <- gsub("TSH", "Serum_TSH", names(gene_stats))
  names(gene_stats) <- gsub("FreeT3", "Serum_Free_T3", names(gene_stats))
  names(gene_stats) <- gsub("FreeT4", "Serum_Free_T4", names(gene_stats))

  # Remove deep and sc
  cols_to_remove <- grep("Deep_and_SC", names(gene_stats), value = T)
  gene_stats <- gene_stats %>% select(-all_of(cols_to_remove))

  names(meta)

  # DEFINE YOUR GENES OF INTEREST

  # This part is to potentially deal with duplicated genes
  ens_genes_of_interest <- NULL
  sym_genes_of_interest <- NULL

  for (gene in genes_of_interest) {
    # Check if the gene is in the list of genes of interest
    if (!any(gene_names$gene %in% gene)) {
      print(paste0("Gene ", gene, " was not found. Please check spelling.Maybe you meant ", grep(paste0("^", substr(gene, 1, 3)), gene_names$gene, value = T)))
      next # Skip to the next iteration of the loop
    }
    ens_genes_of_interest <- c(ens_genes_of_interest, rownames(gene_names)[gene_names$gene %in% gene])
    sym_genes_of_interest <- c(sym_genes_of_interest, gene_names$gene[gene_names$gene %in% gene])
  }


  # Loop over each gene
  for (n in 1:length(ens_genes_of_interest)) {
    gene <- sym_genes_of_interest[n]
    ens_gene <- ens_genes_of_interest[n]

    for (variable in names(meta)[c(3, 4, 7:9)]) {
      # Prepare data for plotting
      plot_data <- data_merged %>%
        select(contains(ens_gene), tissue, subject, all_of(variable)) %>%
        na.omit() %>% # Remove rows with NA
        filter(.data[[variable]] > 0) # Remove 0 values
      plot_data[[variable]] <- as.numeric(plot_data[[variable]])

      # Create plot title
      plot_title <- correct_mg_string(paste(gene, " X ", variable))
      max_y <- max(plot_data[ens_gene])

      # Collect axis information (This plot is only used to get the frame of the plot)
      p <- ggscatter(plot_data,
        x = variable, y = ens_gene,
        add = "reg.line", # Add regression line
        conf.int = TRUE, # Add confidence interval
        color = "tissue", palette = c("#743f1e", "#ffd045"), # Color by groups "cyl"
        label = "subject",
        mean.point = TRUE, # Change point shape by groups "cyl"
        title = plot_title, xlab = variable, ylab = paste0(gene, " (log2CPM)")
      ) +
        stat_cor(
          method = "pearson", aes(color = tissue),
          label.x.npc = "left", label.y = c(max_y * 2, max_y - 1)
        ) +
        scale_x_continuous(limits = c(min(plot_data[variable]), max(plot_data[variable]))) +
        labs(subtitle = ens_gene)

      # Define y-axis range
      p_built <- ggplot_build(p)
      y_range <- plyr::ldply(p_built$layout$panel_params, function(x) x$y.range)
      max_y <- y_range[1, 2]
      min_y <- y_range[1, 1]
      y_range <- abs(min_y - max_y)


      # Define X axis range
      min_x <- min(plot_data[variable])
      max_x <- max(plot_data[variable])

      # Rename SC
      plot_data <- plot_data %>%
        rename("Tissue" = "tissue") %>%
        mutate(Tissue = case_when(
          Tissue == "SubQ" ~ "SC",
          TRUE ~ Tissue
        ))
      ### Extract the stats for annotation=====================
      # Separate the columns by statistics
      gene_stats_list <- list(
        "FDR" = grep("\\.FDR", names(gene_stats), value = T),
        "slope" = grep("\\.slope", names(gene_stats), value = T),
        "pval" = grep("\\.p", names(gene_stats), value = T)
      )

      # Extract the columns
      deep_stat_col <- list()
      sc_stat_col <- list()
      for (columns in names(gene_stats_list)) {
        cols <- grep(variable, gene_stats_list[[columns]], value = T)
        deep_stat_col[[columns]] <- grep("Deep", cols, value = T)
        sc_stat_col[[columns]] <- grep("SC", cols, value = T)
      }

      # Pull the value of the statistics
      deep_stat <- list()
      sc_stat <- list()
      for (stat in names(deep_stat_col)) {
        deep_stat[[stat]] <- gene_stats %>%
          filter(X == ens_gene) %>%
          select(all_of(deep_stat_col[[stat]])) %>%
          pull(!!sym(deep_stat_col[[stat]]))
        sc_stat[[stat]] <- gene_stats %>%
          filter(X == ens_gene) %>%
          select(all_of(sc_stat_col[[stat]])) %>%
          pull(!!sym(sc_stat_col[[stat]]))
      }

      # Define the text format to use in the graph
      # FDR_threshold = 0.15

      # DEEP
      if (deep_stat[["FDR"]] < FDR_threshold) {
        deep.face <- "bold.italic"
      } else {
        deep.face <- "plain"
      }

      # SC
      if (sc_stat[["FDR"]] < FDR_threshold) {
        sc.face <- "bold.italic"
      } else {
        sc.face <- "plain"
      }

      # Write the annotations
      deep_annot <- paste0(
        sprintf("Slope: %.2f", deep_stat[["slope"]]),
        sprintf("; p-val: %.3f", deep_stat[["pval"]]),
        sprintf("; FDR: %.2f", deep_stat[["FDR"]])
      )

      sc_annot <- paste0(
        sprintf("Slope: %.2f", sc_stat[["slope"]]),
        sprintf("; p-val: %.3f", sc_stat[["pval"]]),
        sprintf("; FDR: %.2f", sc_stat[["FDR"]])
      )

      # Plot=========================================================
      p <- ggscatter(plot_data,
        x = variable, y = ens_gene,
        add = "reg.line", # Add regression line
        conf.int = TRUE, # Add confidence interval
        color = "Tissue", palette = c("#E66100", "#5D3A9B"), # Color by tissue
        # label = "subject", #Add subject ID labels
        mean.point = FALSE, # Add a point with the average expression of the group
        title = plot_title, xlab = correct_mg_string(variable), ylab = paste0(gene, " (log2CPM)")
      ) +
        # stat_cor(method = "pearson", aes(color = tissue),
        # label.x.npc = "left", label.y = c(max_y,max_y-(y_range*0.05)) #Adjust the position of the regression stats
        # )    +
        scale_x_continuous(limits = c(min_x, max_x)) +
        # labs(subtitle = ens_gene)+ #This is for when we want to show the ENS ID
        # guides(label = "none") +
        ## The code below should be useless based on the theme_set above
        theme_classic() +
        theme(
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 18),
          plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 16),
          strip.text = element_text(size = 16)
        ) +
        # ADD THE ANNOTATIONS
        geom_text(aes(x = min_x, y = max_y, label = deep_annot),
          size = 4.5, color = "#E66100", fontface = deep.face, hjust = 0
        ) +
        geom_text(aes(x = min_x, y = max_y, label = sc_annot),
          size = 4.5, color = "#5D3A9B", fontface = sc.face, hjust = 0, vjust = 2
        )


      ggsave(paste0("./figures/Paper_figures/Gene_correlations_with_JBC_stats/", variable, "/", gene, "_", ens_gene, ".png"), device = "png", plot = p, width = 6, height = 6, dpi = 300, create.dir = T)

      print(p)
    }
  }
}
