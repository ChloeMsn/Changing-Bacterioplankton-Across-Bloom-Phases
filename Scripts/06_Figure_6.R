# Co-occurrence Networks


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


# Network
library(igraph)

# Plot
library(UpSetR)

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
phyloseq_bacteria_filtered_20µM <- readRDS(paste(phyloseq_path, "phyloseq_bacteria_filtered_20µM.R", sep =
                                                   "/"))
phyloseq_protist_filtered_20µM <- readRDS(paste(phyloseq_path, "phyloseq_protist_filtered_20µM.R", sep =
                                                  "/"))
df_indicsp <- read_excel(paste(res_path, "Data_Figure_5.xlsx", sep = "/"))


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


# II. Co-occurrence NW ---------

# CLR transformation 
my_clr <- function(x) {
  as.vector(compositions::clr(x))
}

pseudo_count <- 1

networkCluster <- function(cluster,
                           meta,
                           asv_table_protist_fraction,
                           asv_table_bacteria_fraction,
                           threshold) {
  # Select Samples
  samples_cluster_bacteria <- meta %>%
    filter(Group == cluster) %>%
    pull(sampleid)
  samples_cluster_protist <- meta %>%
    filter(Group == cluster) %>%
    pull(sampleid_18S)
  
  
  asv_table_bacteria_clr <- dplyr::mutate_all((asv_table_bacteria_fraction[which(rowSums(asv_table_bacteria_fraction[, samples_cluster_bacteria]) > 0), samples_cluster_bacteria]
                                               + pseudo_count), my_clr)
  
  asv_table_protist_clr <- dplyr::mutate_all((asv_table_protist_fraction[which(rowSums(asv_table_protist_fraction[, samples_cluster_protist]) > 0), samples_cluster_protist]
                                              + pseudo_count), my_clr)
  
  # Matrix with Bacteria + Alexandrium Reads
  bacteria_matrix <- as.data.frame(t(asv_table_bacteria_clr[, samples_cluster_bacteria]))
  bacteria_matrix$alex <- as.vector(asv_table_protist_clr[alex_asv, samples_cluster_protist])
  
  # Spearman Correlation
  corr_matrix <- Hmisc::rcorr(as.matrix(bacteria_matrix), type = "spearman")
  corr_matrix_r <- corr_matrix$r
  corr_matrix_p <- corr_matrix$P
  
  # Coefficients
  df_r <- corr_matrix_r %>%
    as.data.frame() %>%
    rownames_to_column("from") %>%
    pivot_longer(-from, names_to = "to", values_to = "R")
  
  # Pvalues
  df_p <- corr_matrix_p %>%
    as.data.frame() %>%
    rownames_to_column("from") %>%
    pivot_longer(-from, names_to = "to", values_to = "p")
  
  # Format as a df and keep only significant correlations
  df_corr <- left_join(df_r, df_p, by = c("from", "to")) %>%
    filter(from != to, # Remove diagonal
           from > to)  %>% # Keep only upper matrix
    mutate(padjust =  p.adjust(p, method = "BH")) %>%
    filter(padjust <= 0.05, abs(R) >= threshold)
  
  print("alex" %in% df_corr$from)
  print("alex" %in% df_corr$to)
  
  return(df_corr)
}

threshold <- 0.6
df_corr_1 <- networkCluster(
  "Before",
  data.frame(sample_data(phyloseq_bacteria_filtered_20µM), check.names = F),
  as.data.frame(otu_table(phyloseq_protist_filtered_20µM), check.names = F),
  as.data.frame(otu_table(phyloseq_bacteria_filtered_20µM), check.names = F),
  threshold
)
df_corr_2 <- networkCluster(
  "Peak",
  data.frame(sample_data(phyloseq_bacteria_filtered_20µM), check.names = F),
  as.data.frame(otu_table(phyloseq_protist_filtered_20µM), check.names = F),
  as.data.frame(otu_table(phyloseq_bacteria_filtered_20µM), check.names = F),
  threshold
)
df_corr_3 <- networkCluster(
  "After",
  data.frame(sample_data(phyloseq_bacteria_filtered_20µM), check.names = F),
  as.data.frame(otu_table(phyloseq_protist_filtered_20µM), check.names = F),
  as.data.frame(otu_table(phyloseq_bacteria_filtered_20µM), check.names = F),
  threshold
)


asv_associated_1 <- df_corr_1 %>% filter(from == "alex" | to == "alex") %>% 
  mutate(asv = str_remove_all(paste0(from,to), "alex|alex")) %>% 
  pull(asv)
tax_table(phyloseq_bacteria_filtered_20µM)[asv_associated_1,]

asv_associated_3 <- df_corr_3 %>% filter(from == "alex" | to == "alex") %>% 
  mutate(asv = str_remove_all(paste0(from,to), "alex|alex")) %>% 
  pull(asv)
tax_table(phyloseq_bacteria_filtered_20µM)[asv_associated_3,]



ecount(graph_from_data_frame(df_corr_1 %>% 
                               select(from, to)))

ecount(graph_from_data_frame(df_corr_2 %>% 
                               select(from, to)))

ecount(graph_from_data_frame(df_corr_3 %>% 
                               select(from, to)))

# II. Edge Betweenness ----------

# Compute edge beetweenness
betweennessEdges <- function(df_corr,
                             asv_table_bacteria,
                             meta,
                             cluster) {
  # Create Global Graph
  plot <- graph_from_edgelist(el = as.matrix(df_corr[, 1:2]), directed = F)
  # Width of edge = Spearman Coeff
  E(plot)$width <- abs(df_corr$R)
  # Edge name
  E(plot)$name <- paste(df_corr$from, df_corr$to, sep = "_")
  
  # print("alex" %in% V(plot)$name)
  # print(grep("alex", E(plot)$name))
  
  # Create Subgraphs and compute betweenness
  samples_cluster <- meta %>%
    filter(Group == cluster) %>%
    pull(sampleid)
  
  e_betweenness <- sampleid <- edge_id <- c()
  for (sample in samples_cluster) {
    # Induce Subgraph
    asv_present_sample <- rownames(asv_table_bacteria)[asv_table_bacteria[, sample] > 0]
    asv_present <- intersect(V(plot)$name, c(asv_present_sample, "alex"))
    sub <- induced_subgraph(plot, asv_present)
    
    # Remove isolated nodes
    isolated_nodes <- which(degree(sub) == 0)
    sub_bis <- delete_vertices(sub, isolated_nodes)
    # print("alex" %in% V(sub_bis)$name)
    # print(grep("alex",  E(sub_bis)$name))
    
    e_betweenness <- append(e_betweenness,
                            edge_betweenness(sub, directed = F, cutoff = -1))
    
    sampleid <- append(sampleid, rep(sample, times = ecount(sub)))
    edge_id <- append(edge_id, E(sub)$name)
    
  }
  
  df_eb <- data.frame(sampleid, edge_id, e_betweenness) %>%
    mutate(
      from = str_split_fixed(edge_id, n = 2, pattern = "_")[, 1],
      to = str_split_fixed(edge_id, n = 2, pattern = "_")[, 2]
    ) %>%
    left_join(meta, by = "sampleid")
  
  return(df_eb)
  
}

## 1. Compute --------

df_eb_1 <- betweennessEdges(
  df_corr_1,
  as.data.frame(otu_table(phyloseq_bacteria_filtered_20µM), check.names = F),
  data.frame(sample_data(phyloseq_bacteria_filtered_20µM), check.names = F),
  "Before"
)
df_eb_2 <- betweennessEdges(
  df_corr_2,
  as.data.frame(otu_table(phyloseq_bacteria_filtered_20µM), check.names = F),
  data.frame(sample_data(phyloseq_bacteria_filtered_20µM), check.names = F),
  "Peak"
)
df_eb_3 <- betweennessEdges(
  df_corr_3,
  as.data.frame(otu_table(phyloseq_bacteria_filtered_20µM), check.names = F),
  data.frame(sample_data(phyloseq_bacteria_filtered_20µM), check.names = F),
  "After"
)

tax <- as.data.frame(tax_table(phyloseq_bacteria_filtered_20µM)) %>%
  mutate(ASV_ID = rownames(tax_table(phyloseq_bacteria_filtered_20µM)))

df_eb_1 %>% 
  filter(from =="alex" | to =="alex") %>% 
  group_by(edge_id, Year) %>% 
  summarise(count = n()) %>% 
  arrange(edge_id, Year)

df_eb_3 %>% 
  filter(from =="alex" | to =="alex") %>% 
  group_by(edge_id, Year) %>% 
  summarise(count = n()) %>% 
  arrange(edge_id, Year)


## 2. Plot ---------
df_eb_1 %>%
  filter(from == "alex" | to == "alex") %>%
  mutate(ASV_ID = stringr::str_remove_all(edge_id, "_alex|alex_")) %>%
  left_join(tax, by = "ASV_ID") %>% 
  pull(class) %>% 
  unique()
df_eb_3 %>%
  filter(from == "alex" | to == "alex") %>%
  mutate(ASV_ID = stringr::str_remove_all(edge_id, "_alex|alex_")) %>%
  left_join(tax, by = "ASV_ID") %>% 
  pull(class) %>% 
  unique()

color <- c("Bacteroidia" = "#D55E00",
           "Gammaproteobacteria" = "#B497D6",
           "Bacteriovoracia" = "#87AE73",
           "Desulfobacteria" =  "#FF6EC7")

order_ASV_p1 <- df_eb_1 %>%
  filter(from == "alex" | to == "alex") %>%
  mutate(ASV_ID = stringr::str_remove_all(edge_id, "_alex|alex_")) %>%
  left_join(tax, by = "ASV_ID") %>%
  mutate(genus = ifelse(genus == "", "Un.", genus),
         label = paste0("paste('", family, " ', italic('", genus, "'))")) %>%
  arrange(class) %>%
  pull(label) %>%
  unique()

p1 <- df_eb_1 %>%
  filter(from == "alex" | to == "alex") %>%
  mutate(ASV_ID = stringr::str_remove_all(edge_id, "_alex|alex_")) %>%
  left_join(tax, by = "ASV_ID") %>%
  mutate(genus = ifelse(genus == "", "Un.", genus),
         label = paste0("paste('", family, " ', italic('", genus, "'))")) %>%
  group_by(ASV_ID) %>%
  mutate(
    count = n(),
    .groups = "drop",
    mean = mean(e_betweenness),
    sd = sd(e_betweenness)
  ) %>%
  ggplot(aes(x = e_betweenness, y = factor(label, levels = rev(order_ASV_p1)))) +
  geom_point(alpha = 0.7,
             size = 2,
             position = position_jitter(width = 0.01),
             color = "grey") +
  geom_point(aes(x = mean, color = class), size = 4, alpha = 0.5) +
  geom_errorbar(aes(
    xmin = mean - sd,
    xmax = mean + sd,
    color = class
  ), width = 0.3, linewidth = 0.7) +
  scale_color_manual(values = color) +
  scale_y_discrete(labels = function(x) parse(text = x )) +
  theme_test() +
  theme(
    axis.text = ggtext::element_markdown(size = 16),
    axis.title = element_text(size = 16),
    legend.position = "bottom",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.justification = "left",
    legend.box.just = "left"
  )+
  labs(x = "Edge Betweenness", y = "", color = "Class")
p1

order_ASV_p2 <- df_eb_3 %>%
  filter(from == "alex" | to == "alex") %>%
  mutate(ASV_ID = stringr::str_remove_all(edge_id, "_alex|alex_")) %>%
  left_join(tax, by = "ASV_ID") %>%
  mutate(genus = ifelse(genus == "", "Un.", genus),
         label =  paste0("paste('", family, " ', italic('", genus, "'))")) %>%
  arrange(class) %>%
  pull(label) %>%
  unique()

p2 <- df_eb_3 %>%
  filter(from == "alex" | to == "alex") %>%
  mutate(ASV_ID = stringr::str_remove_all(edge_id, "_alex|alex_")) %>%
  left_join(tax, by = "ASV_ID") %>%
  mutate(label =  paste0("paste('", family, " ', italic('", genus, "'))")) %>%
  group_by(ASV_ID) %>%
  mutate(
    count = n(),
    .groups = "drop",
    mean = mean(e_betweenness),
    sd = sd(e_betweenness)
  ) %>%
  ggplot(aes(x = e_betweenness, y = factor(label, levels = rev(order_ASV_p2)))) +
  geom_point(alpha = 0.7,
             size = 2,
             position = position_jitter(width = 0.01),
             color = "grey") +
  geom_point(aes(x = mean, color = class), size = 4, alpha = 0.5) +
  geom_errorbar(aes(
    xmin = mean - sd,
    xmax = mean + sd,
    color = class
  ), width = 0.2, linewidth = 0.7) +
  scale_color_manual(values = color) +
  scale_y_discrete(labels = function(x) parse(text = x )) +
  theme_test() +
  theme(
    axis.text = ggtext::element_markdown(size = 16),
    axis.title = element_text(size = 16),
    legend.position = "bottom",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.justification = "left",
    legend.box.just = "left"
  )+
  labs(x = "Edge Betweenness", y = "", color = "Class")
p2


data_p1 <- df_eb_1 %>%
  filter(from == "alex" | to == "alex") %>%
  mutate(ASV_ID = stringr::str_remove_all(edge_id, "_alex|alex_")) %>%
  group_by(ASV_ID) %>%
  summarise(mean = mean(e_betweenness),
            sd = sd(e_betweenness)) %>% 
  left_join(tax, by = "ASV_ID")

data_p3 <- df_eb_3 %>%
  filter(from == "alex" | to == "alex") %>%
  mutate(ASV_ID = stringr::str_remove_all(edge_id, "_alex|alex_")) %>%
  group_by(ASV_ID) %>%
  summarise(mean = mean(e_betweenness),
            sd = sd(e_betweenness)) %>% 
  left_join(tax, by = "ASV_ID")


# III. UpSet ---------

## 1. Before -----------


# ASV of each Year
ps_2013 <- subset_samples(phyloseq_bacteria_filtered_20µM,
                          Year == "y2013" & Group == "Before")
asv_2013 <- taxa_names(prune_taxa(taxa_sums(ps_2013) > 0 , ps_2013))
ps_2014 <- subset_samples(phyloseq_bacteria_filtered_20µM,
                          Year == "y2014" & Group == "Before")
asv_2014 <- taxa_names(prune_taxa(taxa_sums(ps_2014) > 0 , ps_2014))
ps_2015 <- subset_samples(phyloseq_bacteria_filtered_20µM,
                          Year == "y2015" & Group == "Before")
asv_2015 <- taxa_names(prune_taxa(taxa_sums(ps_2015) > 0 , ps_2015))

core <- unique(intersect(asv_2013, intersect(asv_2014, asv_2015)))


# ASV associated with Alexandrium
asv_associated <- unique(c(
  df_corr_1 %>%
    filter(from == "alex" | to == "alex") %>%
    mutate(
      edge_id = paste(from, to, sep = "_"),
      ASV_ID = stringr::str_remove_all(edge_id, "_alex|alex_")
    ) %>%
    pull(ASV_ID)
))

# ASV indicators of the start
asv_indicsp_1 <- df_indicsp %>%
  filter(index == 1) %>%
  pull(ASV_ID)

tax[intersect(asv_associated, core),]
df_corr_1 %>% 
  filter(from == "alex" | to == "alex")

# upset Plot
list_asv <- list(
  "Before Peak" = asv_indicsp_1,
  "2013" = asv_2013,
  "2014" = asv_2014,
  "2015" = asv_2015,
  "Associated" = asv_associated
)

mat_upset <- fromList(list_asv)


p3 <- UpSetR::upset(
  mat_upset,
  nset = 5,
  keep.order = T,
  show.numbers = "yes",
  set_size.show = T,
  set_size.scale_max = 1700,
  text.scale = 1.75,
  sets = c("2015", "2014", "2013", "Before Peak", "Associated")
)

p3

## 2. After ---------



# ASV of each year
ps_2013 <- subset_samples(phyloseq_bacteria_filtered_20µM,
                          Year == "y2013" & Group == "After")
asv_2013 <- taxa_names(prune_taxa(taxa_sums(ps_2013) > 0 , ps_2013))
ps_2014 <- subset_samples(phyloseq_bacteria_filtered_20µM,
                          Year == "y2014" & Group == "After")
asv_2014 <- taxa_names(prune_taxa(taxa_sums(ps_2014) > 0 , ps_2014))
ps_2015 <- subset_samples(phyloseq_bacteria_filtered_20µM,
                          Year == "y2015" & Group == "After")
asv_2015 <- taxa_names(prune_taxa(taxa_sums(ps_2015) > 0 , ps_2015))
core <- unique(intersect(asv_2013, intersect(asv_2014, asv_2015)))


# ASV associated directly with Alexandrium
asv_associated <- df_corr_3 %>%
  filter(from == "alex" | to == "alex") %>%
  mutate(
    edge_id = paste(from, to, sep = "_"),
    ASV_ID = stringr::str_remove_all(edge_id, "_alex|alex_")
  ) %>%
  pull(ASV_ID)

# ASV indicator
asv_indicsp_3 <- df_indicsp %>%
  filter(index == 3) %>%
  pull(ASV_ID)

tax[intersect(asv_associated, core),]
df_corr_3 %>% 
  filter(from == "alex" | to == "alex")


# Upset Plot

list_asv <- list(
  "After Peak" = asv_indicsp_3,
  "2013" = asv_2013,
  "2014" = asv_2014,
  "2015" = asv_2015,
  "Associated" = asv_associated
)
mat_upset <- fromList(list_asv)


p4 <- UpSetR::upset(
  mat_upset,
  nset = 5,
  keep.order = T,
  show.numbers = "yes",
  set_size.show = T,
  set_size.scale_max = 2500,
  text.scale = 1.75,
  sets = c("2015", "2014", "2013", "After Peak", "Associated")
)


# Figure X.  -------

# Convert in usable plot with ggarrange
p3_cowplot <- cowplot::plot_grid(NULL,
                                 p3$Main_bar,
                                 p3$Sizes ,
                                 p3$Matrix,
                                 nrow = 2,
                                 align = 'hv',
                                 rel_widths =  c(1.5, 1.25, 1.5,1),
                                 rel_heights = c(1.5,1.5, 0.9,0.9))
p4_cowplot <- cowplot::plot_grid(NULL,
                                 p4$Main_bar,
                                 p4$Sizes,
                                 p4$Matrix,
                                 nrow = 2,
                                 align = 'hv',
                                 rel_widths =  c(1.5, 1.25, 1.5,1),
                                 rel_heights = c(1.5,1.5, 0.9,0.9))


p_ab <- ggpubr::ggarrange(p1, p2,
                  align = "hv",
                  labels = c("A.", "B."),
                  legend = "none")
p_cd <- ggpubr::ggarrange(p3_cowplot, p4_cowplot,
                  align = "hv",
                  labels = c("C.", "D."))


svg(paste(fig_path, "Figure_6.svg", sep ="/"),
    height = 9,
    width = 16)
ggpubr::ggarrange(p_ab,
                  p_cd,
                  align = "hv",
                  ncol = 1,
                  heights = c(0.75,1.25))
dev.off()


