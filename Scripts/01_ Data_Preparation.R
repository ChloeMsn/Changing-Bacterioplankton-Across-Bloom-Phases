# Data Preparation

# Packages -----

# Data Manipulation
library(tidyverse)
library(readxl)
library(writexl)

# Statistical Analysis
library(phyloseq)


# Global variables -------

# Seed
set.seed(1234)

# Data
env_path <- here::here("..", "Data", "Metadata", "Data_Daoulex_modified.xlsx")

asv_bacteria_path <-   here::here("..", "Data", "16S", "02_report", "final_asv_table.tsv")
meta_bacteria_path <- here::here("..", "Data", "Metadata", "sample_file_16S.xlsx")

asv_protist_path <- here::here("..", "Data", "18S", "02_report", "final_asv_table.tsv")
meta_protist_path <- here::here("..", "Data" , "Metadata", "sample_file_18S.xlsx")


# Output
data_path <- here::here("..", "Data")

# Alexandrium ID
alex_asv <- "49d06a3ae3c9ca45d30c80ff48fc7caf"

# I. Import ----

## 1. ASV -------
asv_table_bacteria <- read.table(
  asv_bacteria_path,
  header = T,
  row.names = 1,
  sep = "\t"
)

asv_table_protist <- read.table(asv_protist_path, header = T, row.names = 1)


## 2. Metadata ------

meta_bacteria <- as.data.frame(read_excel(meta_bacteria_path, sheet = "metadata"))
rownames(meta_bacteria) <- meta_bacteria$sampleid

meta_protist <- as.data.frame(read_excel(meta_protist_path, sheet = "metadata"))
rownames(meta_protist) <- meta_protist$sampleid

## 3. Taxonomy -----

splitTaxonomy <- function(asv_table, organism_type) {
  if (organism_type == "bacteria") {
    tax_table <- as.data.frame(stringr::str_split_fixed(asv_table[, "taxonomy"], ";", n =
                                                          7),
                               row.names = rownames(asv_table))
    colnames(tax_table) <- c("kingdom",
                             "phylum",
                             "class",
                             "order",
                             "family",
                             "genus",
                             "specie")
  }
  else{
    tax_table <- as.data.frame(stringr::str_split_fixed(asv_table$taxonomy, ";", n =
                                                          9),
                               row.names = rownames(asv_table))
    colnames(tax_table) <- c(
      "kingdom",
      "supergroup",
      "division",
      "subdivision",
      "class",
      "order",
      "family",
      "genus",
      "specie"
    )
  }
  
  return(tax_table)
}

tax_table_bacteria <- splitTaxonomy(asv_table_bacteria, "bacteria")
tax_table_protist <- splitTaxonomy(asv_table_protist, "protist")

## 4. Environment ------

env <- read_xlsx(env_path)


# II. Format ------

## 1. Metadata -------

meta_bacteria <- meta_bacteria %>%
  filter(sampleid != "DA13_17") %>% # Samples not analysed
  mutate(
    `Sampling Date` = as.Date(`Sampling Date`, format = "%Y-%m-%d"),
    Alexandrium = unlist(asv_table_protist[alex_asv, sampleid_18S])
  )

meta_protist <- meta_protist %>%
  filter(sampleid != "DA13_25") %>%
  mutate(
    `Sampling Date` = as.Date(`Sampling Date`, format = "%Y-%m-%d"),
    Alexandrium = unlist(asv_table_protist[alex_asv, sampleid_18S])
  )

## 2. Taxonomy -------

filterTaxonomy <- function(tax_table, asv_table, organism_type) {
  if (organism_type == "bacteria") {
    asv_filtered <- tax_table %>%
      dplyr::filter(kingdom == "Bacteria") %>%
      dplyr::filter(order != "Chloroplast" &
                      family != "Mitochondria") %>%
      rownames()
  }
  else{
    asv_filtered <- tax_table %>%
      dplyr::filter(kingdom == "Eukaryota") %>%
      dplyr::filter(division != "Streptophyta") %>%
      dplyr::filter(!subdivision %in% c("Metazoa", "Fungi")) %>%
      dplyr::filter(!class %in% c("Florideophyceae", "Ulvophyceae", "Phaeophyceae")) %>%
      rownames()
  }
  
  asv_table_filtered <- asv_table[asv_filtered, !colnames(asv_table) %in% colnames(tax_table)]
  
  return(list(asv_table_filtered, tax_table[asv_filtered, ]))
}


asv_table_bacteria_filtered <- filterTaxonomy(tax_table_bacteria, asv_table_bacteria, "bacteria")[[1]]
tax_table_bacteria_filtered <- filterTaxonomy(tax_table_bacteria, asv_table_bacteria, "bacteria")[[2]]
# ncol(asv_table_bacteria_filtered) - 1
# nrow(asv_table_bacteria_filtered)
# nrow(tax_table_bacteria_filtered)


asv_table_protist_filtered <- filterTaxonomy(tax_table_protist, asv_table_protist, "protist")[[1]]
tax_table_protist_filtered <- filterTaxonomy(tax_table_protist, asv_table_protist, "protist")[[2]]
# ncol(asv_table_protist_filtered) - 1
# nrow(asv_table_protist_filtered)
# nrow(tax_table_protist_filtered)

## 3. Environment -------

env$`Sampling Date` <- as.Date(env$`Sampling Date`, format = "%d/%m/%Y")
meta_bacteria_env <- left_join(meta_bacteria, env, by = "Sampling Date")
meta_protist_env <- left_join(meta_protist, env, by = "Sampling Date")
rownames(meta_bacteria_env) <- meta_bacteria_env$sampleid
rownames(meta_protist_env) <- meta_protist_env$sampleid


# III. Export to Phyloseq  ------


## 1. All samples  --------

phyloseq_bacteria_raw <- phyloseq(
  otu_table(asv_table_bacteria %>%
              select(-taxonomy), taxa_are_rows = T),
  tax_table(as.matrix(tax_table_bacteria)),
  sample_data(meta_bacteria_env)
)

phyloseq_protist_raw <- phyloseq(
  otu_table(asv_table_protist %>%
              select(-taxonomy), taxa_are_rows = T),
  tax_table(as.matrix(tax_table_protist)),
  sample_data(meta_protist_env)
)

phyloseq_bacteria_filtered <- phyloseq(
  otu_table(
    asv_table_bacteria_filtered %>%
      select(-taxonomy),
    taxa_are_rows = T
  ),
  tax_table(as.matrix(tax_table_bacteria_filtered)),
  sample_data(meta_bacteria_env)
)

phyloseq_protist_filtered <- phyloseq(
  otu_table(
    asv_table_protist_filtered %>%
      select(-taxonomy),
    taxa_are_rows = T
  ),
  tax_table(as.matrix(tax_table_protist_filtered)),
  sample_data(meta_protist_env)
)

saveRDS(
  phyloseq_bacteria_raw,
  paste(data_path, "Phyloseq", "phyloseq_bacteria_raw.R", sep = "/")
)
saveRDS(
  phyloseq_bacteria_filtered,
  paste(data_path, "Phyloseq", "phyloseq_bacteria_filtered.R", sep =
          "/")
)

saveRDS(
  phyloseq_protist_raw,
  paste(data_path, "Phyloseq", "phyloseq_protist_raw.R", sep = "/")
)
saveRDS(
  phyloseq_protist_filtered,
  paste(data_path, "Phyloseq", "phyloseq_protist_filtered.R", sep =
          "/")
)

## 2. Per Size Fraction -------

subsetFraction <- function(ps, fraction) {
  # print(fraction)
  
  # Select Samples 
  ps_fraction <- subset_samples(ps, Fraction == "20µM" & !is.na(A..minutum.Cells.L))
  print(ps_fraction)
  
  # Remove ASV with sum = 0
  ps_fraction_asv <- prune_taxa(taxa_sums(ps_fraction) > 0, ps_fraction)
  print(ps_fraction_asv)
  
  # Remove Samples with no reads
  ps_fraction_asv_samples <- prune_samples(sample_sums(ps_fraction_asv) > 0, ps_fraction_asv)
  print(ps_fraction_asv_samples)
  
  # Warning if samples are removed
  if (nsamples(ps_fraction) != nsamples(ps_fraction_asv_samples)) {
    print("WARNING: SAMPLES WITH NO READS")
  }
  return(ps_fraction_asv_samples)
}

### 20µM ------

phyloseq_bacteria_raw_20µM <- subsetFraction(phyloseq_bacteria_raw, "20µM")
phyloseq_bacteria_filtered_20µM <- subsetFraction(phyloseq_bacteria_filtered, "20µM")
phyloseq_protist_raw_20µM <- subsetFraction(phyloseq_protist_raw, "20µM")
phyloseq_protist_filtered_20µM <- subsetFraction(phyloseq_protist_filtered, "20µM")


saveRDS(
  phyloseq_bacteria_raw_20µM,
  paste(data_path, "Phyloseq", "phyloseq_bacteria_raw_20µM.R", sep =
          "/")
)
saveRDS(
  phyloseq_bacteria_filtered_20µM,
  paste(
    data_path,
    "Phyloseq",
    "phyloseq_bacteria_filtered_20µM.R",
    sep = "/"
  )
)

saveRDS(
  phyloseq_protist_raw_20µM,
  paste(data_path, "Phyloseq", "phyloseq_protist_raw_20µM.R", sep =
          "/")
)
saveRDS(
  phyloseq_protist_filtered_20µM,
  paste(
    data_path,
    "Phyloseq",
    "phyloseq_protist_filtered_20µM.R",
    sep = "/"
  )
)

### 3µM ------
phyloseq_bacteria_raw_3µM <- subsetFraction(phyloseq_bacteria_raw, "3µM")
phyloseq_bacteria_filtered_3µM <- subsetFraction(phyloseq_bacteria_filtered, "3µM")
phyloseq_protist_raw_3µM <- subsetFraction(phyloseq_protist_raw, "3µM")
phyloseq_protist_filtered_3µM <- subsetFraction(phyloseq_protist_filtered, "3µM")


saveRDS(
  phyloseq_bacteria_raw_3µM,
  paste(data_path, "Phyloseq", "phyloseq_bacteria_raw_3µM.R", sep =
          "/")
)
saveRDS(
  phyloseq_bacteria_filtered_3µM,
  paste(
    data_path,
    "Phyloseq",
    "phyloseq_bacteria_filtered_3µM.R",
    sep = "/"
  )
)

saveRDS(
  phyloseq_protist_raw_3µM,
  paste(data_path, "Phyloseq", "phyloseq_protist_raw_3µM.R", sep =
          "/")
)
saveRDS(
  phyloseq_protist_filtered_3µM,
  paste(
    data_path,
    "Phyloseq",
    "phyloseq_protist_filtered_3µM.R",
    sep = "/"
  )
)


### 02µM ------
phyloseq_bacteria_raw_02µM <- subsetFraction(phyloseq_bacteria_raw, "02µM")
phyloseq_bacteria_filtered_02µM <- subsetFraction(phyloseq_bacteria_filtered, "02µM")
phyloseq_protist_raw_02µM <- subsetFraction(phyloseq_protist_raw, "02µM")
phyloseq_protist_filtered_02µM <- subsetFraction(phyloseq_protist_filtered, "02µM")


saveRDS(
  phyloseq_bacteria_raw_02µM,
  paste(data_path, "Phyloseq", "phyloseq_bacteria_raw_02µM.R", sep =
          "/")
)
saveRDS(
  phyloseq_bacteria_filtered_02µM,
  paste(
    data_path,
    "Phyloseq",
    "phyloseq_bacteria_filtered_02µM.R",
    sep = "/"
  )
)

saveRDS(
  phyloseq_protist_raw_02µM,
  paste(data_path, "Phyloseq", "phyloseq_protist_raw_02µM.R", sep =
          "/")
)
saveRDS(
  phyloseq_protist_filtered_02µM,
  paste(
    data_path,
    "Phyloseq",
    "phyloseq_protist_filtered_02µM.R",
    sep = "/"
  )
)
