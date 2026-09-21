# Indicspecies

# Packages -----

# Data Manipulation
library(tidyverse)
library(readxl)
library(writexl)


# Statistical Analysis
library(vegan)
library(phyloseq)
library(rstatix)
library(indicspecies)

# Seed
set.seed(1234)


# Data
phyloseq_path <- here::here("..", "Data", "Phyloseq")
fig_path <-   here::here("..", "Analysis", "Final Figures")
res_path <- here::here("..", "Analysis", "Final Figures", "Tables")

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



# Alexandrium ID
alex_asv <- "49d06a3ae3c9ca45d30c80ff48fc7caf"


label_group <- c("Before" = "Before Peak",
                 "Peak" = "Peak",
                 "After" = "After Peak")


# I. Import --------
phyloseq_bacteria_filtered_20µM <- readRDS(paste(phyloseq_path, "phyloseq_bacteria_filtered_20µM.R", sep =
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

phyloseq_bacteria_filtered_CSS_20µM <- metagMisc::phyloseq_transform_css(phyloseq_bacteria_filtered_20µM)




# saveRDS(phyloseq_bacteria_filtered_CSS_20µM, paste(phyloseq_path, "phyloseq_bacteria_filtered_CSS_20µM.R", sep =
#                                                      "/"))



# II. Indicspecies --------

indicspAnalysis <- function(df_dbrda, asv_table_CSS, cluster_column) {
  samples <- df_dbrda$sampleid
  clusters <- df_dbrda[, cluster_column]
  indicsp <- indicspecies::multipatt(t(asv_table_CSS[, samples]), clusters, duleg = T)
  
  return(indicsp$sign %>%
           filter(p.value <= 0.05))
}



indicsp_Group <- indicspAnalysis(
  data.frame(
    sample_data(phyloseq_bacteria_filtered_CSS_20µM),
    check.names = F
  ),
  otu_table(phyloseq_bacteria_filtered_CSS_20µM),
  "Group"
)

indicsp_Group %>% 
  nrow()

indicsp_Group %>% 
  group_by(index) %>% 
  summarise(count = n())

df_indicsp_tax <- cbind(indicsp_Group,
                        tax_table(phyloseq_bacteria_filtered_CSS_20µM)[rownames(indicsp_Group),]) %>% 
  rownames_to_column("ASV_ID")

write_xlsx(df_indicsp_tax, paste(res_path, "Data_Figure_5.xlsx", sep ="/"))
write_xlsx(df_indicsp_tax %>% 
             group_by(index, kingdom, phylum, class, order, family) %>% 
             summarise(count = n()),
           paste(res_path, "Data_Figure_5_count.xlsx", sep ="/"))

asv_indicsp <- c(
  df_indicsp_tax %>% filter(index == 1) %>% arrange(class, family, genus) %>%  pull(ASV_ID),
  df_indicsp_tax %>% filter(index == 2) %>%  arrange(class, family, genus) %>%  pull(ASV_ID),
  df_indicsp_tax %>% filter(index == 3) %>%  arrange(class, family, genus) %>%  pull(ASV_ID)
)

asv_indicsp_before <- df_indicsp_tax %>% filter(index == 1) %>%  pull(ASV_ID)
asv_indicsp_peak <- df_indicsp_tax %>% filter(index == 2) %>%  pull(ASV_ID)
asv_indicsp_after <- df_indicsp_tax %>% filter(index == 3) %>%  pull(ASV_ID)


# III.  Heatmap ----------

order_samples <- data.frame(sample_data(phyloseq_bacteria_filtered_CSS_20µM)) %>%
  arrange(Group, sampleid) %>%
  pull(sampleid)

order_asv <- rev(asv_indicsp)

label_asv <- data.frame(tax_table(phyloseq_bacteria_filtered_CSS_20µM)) %>%
  rownames_to_column("ASV_ID") %>%
  filter(ASV_ID %in% order_asv) %>%
  arrange(ASV_ID) %>%
  mutate(
    genus = ifelse(genus == "", "Unk. Genus", genus),
    family = ifelse(family == "", "Unk. Family", family),
    label = paste0(class, " ", family, " <i>", genus, "</i>")
  ) %>%
  pull(label)

names(label_asv) <-  data.frame(tax_table(phyloseq_bacteria_filtered_CSS_20µM)) %>%
  rownames_to_column("ASV_ID") %>%
  filter(ASV_ID %in% order_asv) %>%
  arrange(ASV_ID) %>%
  pull(ASV_ID)


facet_label_x <- c("Before" = "<span style='color:darksalmon'>Before Peak</span>" ,
                   "Peak" = "<span style='color:darkorchid4'>Peak</span>" ,
                   "After" = "<span style='color:cadetblue3'>After Peak</span>")
facet_label_y <- c(                 "1" = "<span style='color:darksalmon'>Before Peak Indicators</span>" ,
                                    "2" = "<span style='color:darkorchid4'>Peak Indicators</span>" ,
                                    "3" = "<span style='color:cadetblue3'>After Peak Indicators</span>")


svg(
  paste(fig_path, "Figure_5.svg", sep = "/"),
  height = 18,
  width = 12
)
data.frame(otu_table(phyloseq_bacteria_filtered_CSS_20µM)) %>%
  rownames_to_column("ASV_ID") %>%
  pivot_longer(cols = -ASV_ID,
               names_to = "sampleid",
               values_to = "Abundance") %>%
  filter(ASV_ID %in% asv_indicsp) %>%
  filter(sampleid %in% order_samples) %>%
  left_join(data.frame(sample_data(phyloseq_bacteria_filtered_CSS_20µM)) , by = "sampleid") %>%
  left_join(data.frame(tax_table(phyloseq_bacteria_filtered_CSS_20µM)) %>%
              rownames_to_column("ASV_ID"), by = "ASV_ID") %>%
  mutate(
    status = case_when(
      ASV_ID %in% asv_indicsp_before ~ "1",
      ASV_ID %in% asv_indicsp_peak ~ "2",
      ASV_ID %in% asv_indicsp_after ~ "3"
    ),
    Group = factor(Group, levels = c("Before", "Peak", "After"))
  ) %>%
  ggplot(aes(
    x = as.factor(Sampling.Date),
    y = factor(ASV_ID, levels = order_asv),
    # y = Label,
    fill = Abundance
  )) +
  geom_tile() +
  facet_grid(status ~ Group,
             labeller = labeller(
               Group  = as_labeller(facet_label_x),
               status = as_labeller(facet_label_y)), scales = "free", space = "free_y") +
  # scale_fill_viridis_c() +
  scale_fill_gradientn(
    colors = viridis::viridis(5),
    values = scales::rescale(c(0, 5, 7, 9, 12)),
    breaks = c(0, 5, 7, 9, 12),
    limits = c(0, 12)
  ) +
  scale_y_discrete(labels = as_labeller(label_asv)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, size = 8),
    axis.text.y = ggtext::element_markdown(size = 8),
    strip.text.x = ggtext::element_markdown(size = 12),
    strip.text.y = ggtext::element_markdown(size = 12),
    strip.background.y = element_blank()
  ) +
  labs(x = "", y = "", fill = "Reads Abundance")
dev.off()
