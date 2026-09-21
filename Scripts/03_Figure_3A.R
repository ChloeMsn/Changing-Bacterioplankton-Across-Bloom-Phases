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
library(rstatix)

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

relabel_group <- c("Before" = "Before Peak",
                   "Peak" = "Peak",
                   "After" = "After Peak")

# Theme -------

my_theme <- theme(
  axis.text.x = element_text(
    size = 10,
    angle = 90,
    hjust = 1
  ),
  axis.text.y = element_text(size = 10),
  # strip.text = element_text(size = 10),
  # legend.text = element_text(size = 10),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 10),
  strip.background = element_rect(
    fill = "white",
    color = "black",
    linewidth = 1
  ),
  strip.text = element_text(face = "bold", size = 12),
  panel.border = element_rect(
    color = "black",
    fill = NA,
    linewidth = 1
  )
)


# I. Import ---------

phyloseq_bacteria_filtered_20µM <- readRDS(paste(phyloseq_path, "phyloseq_bacteria_filtered_20µM.R", sep =
                                                   "/"))
phyloseq_protist_filtered_20µM <- readRDS(paste(phyloseq_path, "phyloseq_protist_raw_20µM.R", sep =
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

# II. Alpha Diversity -------

df_alpha_bacteria_20µM <- plot_richness(phyloseq_bacteria_filtered_20µM,
                                        measures = c("Observed", "Shannon", "Simpson"))$data
df_alpha_protist_20µM <- plot_richness(phyloseq_protist_filtered_20µM,
                                       measures = c("Observed", "Shannon", "Simpson"))$data

write_xlsx(
  df_alpha_bacteria_20µM %>%
    pivot_wider(names_from = variable, values_from = value) %>%
    select(sampleid, Sampling.Date, Group, Observed, Shannon, Simpson),
  paste(res_path, "Data_Figure_3A.xlsx", sep =
          "/")
)

write_xlsx(
  df_alpha_protist_20µM %>%
    pivot_wider(names_from = variable, values_from = value) %>%
    select(sampleid, Sampling.Date, Group, Observed, Shannon, Simpson),
  paste(res_path, "Data_Figure_S1A.xlsx", sep =
          "/")
)


# III. Plot -------

## 1. Bacteria --------

# Figure 3A
signif_bacteria_20µM <- df_alpha_bacteria_20µM %>%
  mutate(
    Index = stringr::str_replace_all(variable, "Observed", "Species Richness"),
    Group = Group
  ) %>%
  group_by(Index) %>%
  dunn_test(value ~ Group) %>%
  add_xy_position(x = "Group", scales = "free_y")

png(
  paste(fig_path, "Figure_3A.png", sep = "/"),
  res = 300,
  height = 1500,
  width = 2200
)
df_alpha_bacteria_20µM %>%
  mutate(Index = stringr::str_replace_all(variable, "Observed", "Species Richness")) %>%
  ggplot() +
  geom_boxplot(aes(
    x = as.factor(Group),
    y = value,
    fill = as.factor(Group),
    color = as.factor(Group)
  ),
  alpha = 0.6) +
  scale_x_discrete(labels = as_labeller(relabel_group)) +
  scale_fill_manual(labels =
                      as_labeller(relabel_group), values = color_group) +
  scale_color_manual(values = color_group, labels = as_labeller(relabel_group)) +
  facet_wrap(~ Index, scales = "free_y") +
  ggpubr::stat_pvalue_manual(signif_bacteria_20µM, hide.ns = T) +
  theme_test() +
  my_theme +
  labs(x = "",
       y = "Value",
       fill = "Bloom Stage",
       color = "Bloom Stage")
dev.off()

## 2. Protist -------

signif_protist_20µM <- df_alpha_protist_20µM %>%
  mutate(
    Index = stringr::str_replace_all(variable, "Observed", "Species Richness"),
    Group = Group
  ) %>%
  group_by(Index) %>%
  dunn_test(value ~ Group) %>%
  add_xy_position(x = "Group", scales = "free_y")

# Supplementary Figure 1
png(
  paste(fig_path, "Figure_S1A.png", sep = "/"),
  res = 300,
  height = 1500,
  width = 2200
)
df_alpha_protist_20µM %>%
  mutate(Index = stringr::str_replace_all(variable, "Observed", "Species Richness")) %>%
  ggplot() +
  geom_boxplot(aes(
    x = as.factor(Group),
    y = value,
    fill = as.factor(Group),
    color = as.factor(Group)
  ),
  alpha = 0.6) +
  scale_fill_manual(values = color_group, labels = as_labeller(relabel_group)) +
  scale_color_manual(values = color_group, labels = as_labeller(relabel_group)) +
  scale_x_discrete(labels = as_labeller(relabel_group)) +
  facet_wrap(~ Index, scales = "free_y") +
  ggpubr::stat_pvalue_manual(signif_protist_20µM, hide.ns = T) +
  theme_test() +
  my_theme +
  labs(
    x = "",
    y = "Value",
    fill = "Bloom Stage",
    color = "Bloom Stage",
    y = ""
  )
dev.off()
