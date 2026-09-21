
# Packages -----

# Data Manipulation
library(tidyverse)
library(readxl)
library(writexl)


# Statistical Analysis
library(vegan)
library(phyloseq)
library(fantaxtic)
library(rstatix)


# Global variables -------

# Seed
set.seed(1234)

color_year <- c(
  "y2013" = "#97C1A9",
  "y2014" = "#FA9BCF",
  "y2015" = "#CDD0F8"
)
color_fraction <- c(
  "0.2µM" = "#B8D8F7",
  "3µM" = "#FFCCB8",
  "20µM" = "#F9E1AB"
)
color_group <- c(
  "Before" = "darksalmon",
  "Peak" = "darkorchid4",
  "After" = "cadetblue3"
)
level_fraction <- c("0.2µM", "3µM", "20µM")

# Output
fig_path <-   here::here("..", "Analysis", "Final Figures")
res_path <- here::here("..", "Analysis", "Final Figures", "Tables")
phyloseq_path <- here::here("..", "Data", "Phyloseq")

# Alexandrium ID
alex_asv <- "49d06a3ae3c9ca45d30c80ff48fc7caf"

# I. Import --------
phyloseq_protist_filtered_20µM <- readRDS(paste(phyloseq_path, "phyloseq_protist_filtered_20µM.R", sep =
                                                  "/"))


phyloseq_alexandrium <- subset_taxa(phyloseq_protist_filtered_20µM, genus == "Alexandrium")
phyloseq_alexandrium

standf = function(x) ((x / sum(x)))
phyloseq_alexandrium_relative <- transform_sample_counts(phyloseq_alexandrium, standf)


relative_abundance <- rowSums(phyloseq_alexandrium@otu_table) / sum(rowSums(phyloseq_alexandrium@otu_table))
taxonomy <- phyloseq_alexandrium@tax_table %>% 
  as.data.frame() %>% 
  mutate(taxonomy = paste(family, genus, specie, sep = " ")) %>% 
  pull(taxonomy)
asv_id <- taxa_names(phyloseq_alexandrium)


write_xlsx(data.frame(asv_id, taxonomy, relative_abundance),
           paste(res_path, "Supplementary_Table2.xlsx",sep = "/"))
