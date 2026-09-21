# Taxonomic Composition

# Packages -----

# Data Manipulation
library(tidyverse)
library(readxl)
library(writexl)

# Statistical Analysis
library(phyloseq)
library(metagMisc)
library(fantaxtic)

library(ggalluvial)
# Global variables -------

# Seed
set.seed(1234)

# Data
phyloseq_path <- here::here("..", "Data", "Phyloseq")

# Output
fig_path <-   here::here("..", "Analysis", "Final Figures")
res_path <- here::here("..", "Analysis", "Final Figures", "Tables")
data_path <- here::here("..", "Data")

# Alexandrium ID
alex_asv <- "49d06a3ae3c9ca45d30c80ff48fc7caf"

color_group <- c(
  "Before" = "darksalmon",
  "Peak" = "darkorchid4",
  "After" = "cadetblue3"
)

relabel_group <- c(
  "Before" = "Before Peak", "Peak" = "Peak", "After" = "After Peak"
)

level_group <- c("Before", "Peak", "After")

# Theme -------

my_theme <- theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
                  axis.text.y = element_text(size = 10),
                  strip.text = element_text(size = 10),
                  legend.text = element_text(size = 10),
                  legend.title = element_text(size = 10))

# I. Import ---------

phyloseq_bacteria_filtered_20µM <- readRDS(paste(phyloseq_path, "phyloseq_bacteria_filtered_20µM.R", sep =
                                                   "/"))
phyloseq_protist_filtered_20µM <- readRDS(paste(phyloseq_path, "phyloseq_protist_filtered_20µM.R", sep =
                                                  "/"))



sample_data(phyloseq_bacteria_filtered_20µM) <- sample_data(
  sample_data(phyloseq_bacteria_filtered_20µM) %>%
    data.frame() %>%
    arrange(Sampling.Date) %>%
    mutate(
      Group =
        c(
          "Before",
          "Before",
          "Peak",
          "Peak",
          "Peak",
          "After",
          "After",
          "Before",
          "Before",
          "Peak",
          "Peak",
          "After",
          "Peak",
          "Peak",
          "After",
          "After",
          "After",
          "After",
          "After",
          "After",
          "After",
          "After",
          "After",
          "After",
          "After",
          "After",
          "Before",
          "Before",
          "Before",
          "Peak",
          "After",
          "After",
          "Before",
          "Peak",
          "After"
        ),
      Group = factor(Group, levels = c("Before", "Peak", "After"))
    )
)

sample_data(phyloseq_protist_filtered_20µM) <- sample_data(
  sample_data(phyloseq_protist_filtered_20µM) %>%
    data.frame() %>%
    arrange(Sampling.Date) %>%
    mutate(
      Group =
        c(
          "Before",
          "Before",
          "Peak",
          "Peak",
          "Peak",
          "After",
          "After",
          "Before",
          "Before",
          "Peak",
          "Peak",
          "After",
          "Peak",
          "Peak",
          "After",
          "After",
          "After",
          "After",
          "After",
          "After",
          "After",
          "After",
          "After",
          "After",
          "After",
          "After",
          "Before",
          "Before",
          "Before",
          "Peak",
          "After",
          "After",
          "Before",
          "Peak",
          "After"
        ),
      Group = factor(Group, levels = c("Before", "Peak", "After"))
    )
)



# II. Harmonize Taxonomy --------

# Empty Tax Rank filled with "Unknown"
harmoniseTaxonomy <- function(ps) {
  tax_table(ps)[tax_table(ps) == ""] <- "Unassigned"
  return(ps)
}
phyloseq_bacteria_filtered_harmo_taxo_20µM <- harmoniseTaxonomy(phyloseq_bacteria_filtered_20µM)
phyloseq_protist_filtered_harmo_taxo_20µM <- harmoniseTaxonomy(phyloseq_protist_filtered_20µM)


# III. Plot ---------

my_tax_rank <- "genus"
n_top_taxa <- 10


## 1. Protist ---------

phyloseq_protist_filtered_harmo_taxo_merged_20µM <- merge_samples(phyloseq_protist_filtered_harmo_taxo_20µM,
                                                                  "Group")

top_taxa_protist_merged <- fantaxtic::top_taxa(
  phyloseq_protist_filtered_harmo_taxo_merged_20µM,
  tax_level = my_tax_rank,
  n_taxa = n_top_taxa
)

df_barplot_protist_merged_20µM <- plot_nested_bar(
  ps_obj = top_taxa_protist_merged$ps_obj,
  top_level = "class",
  nested_level = "genus",
  include_rank = T
)$data

color_pal_protist <- c(
  "Other" = "lightgrey",
  "Alexandrium" =  "darkblue",
  "Gonyaulax" = "deepskyblue2",
  "Heterocapsa" = "dodgerblue",
  "Unassigned" = "cornflowerblue",
  "Chaetoceros" = "darkseagreen1",
  "Leptocylindrus" = "darkseagreen3",
  "Thalassiosira" = "darkseagreen4",
  "Teleaulax" = "coral",
  "Tintinnopsis_05"  = "pink",
  "Tintinnopsis_07" = "lightpink3"
)

barplot_protist_merged_20µM <- df_barplot_protist_merged_20µM %>%
  group_by(Sample) %>%
  mutate(rel_ab = Abundance / sum(Abundance)) %>%
  ggplot(aes(x = factor(Sample, levels = level_group), y = rel_ab, fill = genus)) +
  geom_col() +
  scale_fill_manual(values = color_pal_protist) +
  scale_x_discrete(labels = as_labeller(relabel_group)) +
  theme_test() +
  my_theme +
  labs(x = "", y = "Relative Sequence Abundance", fill = "Genus")

legend_barplot_protist_merged_20µM <- df_barplot_protist_merged_20µM %>%
  ggplot(aes(y = factor(genus, level = rev(
    names(color_pal_protist)
  )), fill = genus)) +
  geom_bar() +
  scale_fill_manual(values = color_pal_protist) +
  theme_void() +
  theme(axis.text.y = element_text(size = 12), legend.position = "none")

sankey_protist_20µM <- df_barplot_protist_merged_20µM %>% 
  group_by(Sample) %>% 
  mutate(Relative_Abundance = (Abundance) / sum(Abundance)) %>% 
  ggplot(
    aes(x = as.factor(Sample),
        stratum = genus,
        alluvium = genus,
        y = Relative_Abundance,
        fill = genus,
        color = genus)) +
  geom_flow(alpha = 1) +
  geom_stratum(alpha = 1, color = NA) +
  scale_fill_manual(values = color_pal_protist) +
  scale_color_manual(values = color_pal_protist) +
  scale_x_discrete(limits = level_group,
                   labels = as_labeller(relabel_group)) +
  theme_test() +
  my_theme +
  theme(
    legend.position = "none"
  ) +
  labs(x = "", y = "Relative Sequence Abundance")
sankey_protist_20µM

## 2. Bacteria --------

phyloseq_bacteria_filtered_harmo_taxo_merged_20µM <- merge_samples(phyloseq_bacteria_filtered_harmo_taxo_20µM,
                                                                   "Group")


top_taxa_nested_bacteria_merged <- fantaxtic::nested_top_taxa(
  phyloseq_bacteria_filtered_harmo_taxo_merged_20µM,
  top_tax_level = "class",
  nested_tax_level = "family",
  n_top_taxa = 4,
  n_nested_taxa = 3
)

pal <- c("Bacteroidia" = "#D55E00", "Alphaproteobacteria" = "#0072B2", "Gammaproteobacteria" = "lavender",
         "Desulfobulbia" = "#CC79A7", "Other" = "grey")

barplot_bacteria_nested_merged_20µM <- plot_nested_bar(
  ps_obj = top_taxa_nested_bacteria_merged$ps_obj,
  top_level = "class",
  nested_level = "family",
  include_rank = T,
  relative_abundances = T, pal = pal
) +
  scale_x_discrete(limits = level_group,
                   labels = as_labeller(relabel_group)) +
  theme_test() +
  my_theme +
  theme(
    legend.position = "none"
  ) +
  labs(x = "", y = "Relative Sequence Abundance", fill = "Genus")

color_pal_bacteria_nested <- barplot_bacteria_nested_merged_20µM$data %>%
  select(subgroup_colour) %>%
  unique() %>%
  pull(subgroup_colour)
names <- barplot_bacteria_nested_merged_20µM$data %>%
  select(family) %>%
  unique() %>%
  pull(family)

names(color_pal_bacteria_nested) <- names


levels <- c(
  "Crocinitomicaceae",
  "Flavobacteriaceae",
  "NS9 marine group",
  "Other Bacteroidia",
  "Paracoccaceae",
  "SAR116 clade",
  "Sphingomonadaceae",
  "Other Alphaproteobacteria",
  "Alteromonadaceae",
  "Halieaceae",
  "Vibrionaceae",
  "Other Gammaproteobacteria",
  "Desulfobulbaceae",
  "Desulfocapsaceae",
  "Other"
)

legend_barplot_bacteria_nested_merged_20µM <- barplot_bacteria_nested_merged_20µM$data %>%
  ggplot(aes(y = factor(family, levels = rev(levels)), fill = family)) +
  geom_bar() +
  scale_fill_manual(values = color_pal_bacteria_nested) +
  theme_void() +
  theme(axis.text.y = element_text(size = 12), legend.position = "none")


sankey_bacteria_20µM <- barplot_bacteria_nested_merged_20µM$data %>% 
  group_by(Sample) %>% 
  mutate(Relative_Abundance = (Abundance) / sum(Abundance),
         family = factor(family, levels = levels)) %>% 
  ggplot(
    aes(x = as.factor(Sample),
        stratum = family,
        alluvium = family,
        y = Relative_Abundance,
        fill = family,
        color = family)) +
  geom_flow(alpha = 1) +
  geom_stratum(color = NA) +
  scale_fill_manual(values = color_pal_bacteria_nested) +
  scale_color_manual(values = color_pal_bacteria_nested) +
  scale_x_discrete(limits = level_group,
                   labels = as_labeller(relabel_group)) +
  theme_test() +
  my_theme +
  theme(
    legend.position = "none"
  ) +
  labs(x = "", y = "Relative Sequence Abundance")
sankey_bacteria_20µM

# IV. Export --------

## 1. Plot ---------
png(paste(fig_path, "Figure_3B.png", sep = "/"),
    res = 200,
    height = 900,
    width = 900)
print(sankey_bacteria_20µM)
dev.off()

svg(paste(fig_path, "Figure_3B.svg", sep = "/"))
print(sankey_bacteria_20µM)
dev.off()


png(paste(fig_path, "Figure_S1B.png", sep = "/"),
    res = 200,
    height = 900,
    width = 900)
print(sankey_protist_20µM)
dev.off()

svg(paste(fig_path, "Figure_S1B.svg", sep = "/"))
print(sankey_protist_20µM)
dev.off()

## 2. Legend -------

svg(paste(fig_path, "Fig_3B_legend.svg", sep = "/"))
print(legend_barplot_bacteria_nested_merged_20µM)
dev.off()

svg(paste(fig_path, "Fig_S1B_legend.svg", sep = "/"))
print(legend_barplot_protist_merged_20µM)
dev.off()



## 3. Data ---------


write_xlsx(
  barplot_protist_merged_20µM$data %>%
    select(Sample, rel_ab, group_subgroup) %>% 
    pivot_wider(names_from = Sample,
                values_from = rel_ab),
  paste(res_path, "Data_Figure_S1B.xlsx", sep =
          "/")
)


write_xlsx(
  barplot_bacteria_nested_merged_20µM$data %>%
    group_by(Sample) %>%
    mutate(rel_ab = Abundance / sum(Abundance)) %>% 
    select(Sample, rel_ab, group_subgroup) %>% 
    pivot_wider(names_from = Sample,
                values_from = rel_ab),
  paste(res_path, "Data_Figure_3B.xlsx", sep =
          "/")
)
