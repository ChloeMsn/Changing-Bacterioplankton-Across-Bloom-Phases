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
phyloseq_path <- here::here("..", "Data", "Phyloseq")

# Output
fig_path <-   here::here("..", "Analysis", "Final Figures")
res_path <- here::here("..", "Analysis", "Final Figures")

color_cluster <- c(
  "1" = "darksalmon",
  "2" = "darkorchid4",
  "3" = "cadetblue3"
)

color_group <- c(
  "Before" = "darksalmon",
  "Peak" = "darkorchid4",
  "After" = "cadetblue3"
)

year_label <- c("y2013" = "2013", "y2014" = "2014", "y2015" = "2015")

# I. Import ---------
phyloseq_bacteria_filtered_20µM <- readRDS(paste(phyloseq_path, "phyloseq_bacteria_filtered_20µM.R", sep =
                                                   "/"))

# Attribute a "bloom stage" variable for each sample 
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
        )
    )
)

# II. Plot --------

png(paste(fig_path, "Figure_2.png", sep = "/"),
    res = 200,
    width = 1500,
    height = 1000)
sample_data(phyloseq_bacteria_filtered_20µM) %>%
  ggplot(aes(
    x = `Progression.throught.the.bloom`,
    y = log(`A..minutum.Cells.L`),
    color = factor(Group, levels = c("Before", "Peak", "After")),
    group = "Year"
  )) +
  geom_line(linewidth = 0.5, color = "black", linetype = "dotted") +
  geom_point(size = 3) +
  facet_grid(~ Year,
             labeller = as_labeller(year_label)) +
  theme_classic() +
  scale_color_manual(values = color_group, labels = as_labeller(c(
    "Before" = "Before Peak",
    "Peak" = "Peak",
    "After" = "After Peak"
  ))) +
  labs(color = "Bloom Stage", x = "Progression through sampling (%)", y = "log(A.minutum Cells/L)") +
  theme(legend.text = element_text(size = 10),
        strip.background = element_rect(fill="white", color="black", linewidth =1), 
        strip.text = element_text(face="bold", size=12),                       
        panel.border = element_rect(color="black", fill=NA, linewidth=1))
dev.off()