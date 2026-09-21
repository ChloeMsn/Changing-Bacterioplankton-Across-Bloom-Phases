# Beta Diversity

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
library(rdacca.hp)

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
relabel_group <- c("Before" = "Before Peak",
                   "Peak" = "Peak",
                   "After" = "After Peak")
label_year <- c("y2013" = "2013",
                "y2014" = "2014",
                "y2015" = "2015")
level_group <- c("Before", "Peak", "After")


# Data
phyloseq_path <- here::here("..", "Data", "Phyloseq")
fig_path <-   here::here("..", "Analysis", "Final Figures")
res_path <- here::here("..", "Analysis", "Final Figures", "Tables")


# Theme ---------
my_theme <- theme(
  axis.text.x = element_text(size = 10, hjust = 1),
  axis.text.y = element_text(size = 10),
  strip.text = element_text(size = 10),
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 10)
)

# I. Import -------

phyloseq_bacteria_filtered <- readRDS(paste(phyloseq_path, "phyloseq_bacteria_filtered.R", sep =
                                              "/"))
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




# II. CSS Transformation --------

phyloseq_bacteria_filtered_CSS_20µM <- metagMisc::phyloseq_transform_css(phyloseq_bacteria_filtered_20µM)

# III. BC Distance -------

dist_bacteria_20µM <- distance(phyloseq_bacteria_filtered_CSS_20µM, "bray")

# IV. Env Matrix ----------
# Format Environmental Matrix
getEnv <- function(meta_table) {
  env <- meta_table %>%
    select(
      `Tidal.coefficient`,
      `Ammonium..µM.`,
      `Nitrogen.oxide...µM.`,
      `Phosphate..µM.`,
      `Silicate...µM.`,
      `Mignonne.outflow..m.3.s..1.`,
      `Temperature...C.`,
      Salinity,
      `Irradiance..W.m..2.`,
      `A..minutum.Cells.L`,
      # `Progression.throught.the.bloom`,
      Alexandrium
    ) %>%
    scale() %>%
    as.data.frame()
  return(env)
}

env_dbrda_bacteria_20µM <- getEnv(data.frame(
  sample_data(phyloseq_bacteria_filtered_CSS_20µM),
  check.names = F
))


# V. dbRDA --------

dbRDA <- function(dist, env, meta_table) {
  # dbRDA
  dbRDA <- vegan::dbrda(dist ~ ., data = env, na.action = na.omit)
  
  # extract scores
  df_dbrda <- as.data.frame(vegan::scores(dbRDA)$sites)
  df_env <- as.data.frame(vegan::scores(dbRDA)$biplot)
  # print(rownames(df_env))
  
  # select significant variable
  signif_env <- anova.cca(dbRDA, by = "terms")
  signif_env <- rownames(signif_env)[which(signif_env$`Pr(>F)` <= 0.05)]
  # print(signif_env)
  
  df_env <- df_env[signif_env, ]
  df_env$variable <- stringr::str_remove_all(signif_env, "\\`")
  
  df_dbrda <- left_join(df_dbrda %>%
                          rownames_to_column("sampleid"), meta_table, by =
                          "sampleid")
  
  
  return(list(df_dbrda, df_env, dbRDA))
}

getVarExplaineEigenval <- function(ord) {
  var_expl <- round(eigenvals(ord) / sum(abs(eigenvals(ord))), 4) * 100
  return(var_expl)
}

## 1. Models ---------

dbRDA_bacteria_20µM <- dbRDA(
  dist_bacteria_20µM,
  env_dbrda_bacteria_20µM %>%
    select(-Alexandrium),
  data.frame(
    sample_data(phyloseq_bacteria_filtered_CSS_20µM),
    check.names = F
  )
)


## 2. Test Significance of the model ---------

res_dbrda_bacteria_20µM <- anova.cca(dbRDA_bacteria_20µM[[3]])
res_dbrda_bacteria_20µM_terms <- anova.cca(dbRDA_bacteria_20µM[[3]], by =
                                             "terms")
res_dbrda_bacteria_20µM_axis <- anova.cca(dbRDA_bacteria_20µM[[3]], by =
                                            "axis")
contribution_dbRDA_20µM <- rdacca.hp(
  dv = dist_bacteria_20µM,
  iv = env_dbrda_bacteria_20µM %>%
    select(
      -Alexandrium,
      -`Tidal.coefficient`,
      -`Mignonne.outflow..m.3.s..1.`,
      -`Irradiance..W.m..2.`
    ),
  method = "dbRDA"
)

## 3. Plot ------

xlabel_20µM <- paste0("dbRDA1 (",
                      getVarExplaineEigenval(dbRDA_bacteria_20µM[[3]])[1],
                      "%)")
ylabel_20µM <- paste0("dbRDA1 (",
                      getVarExplaineEigenval(dbRDA_bacteria_20µM[[3]])[2],
                      "%)")

### a) Progression -------

dbRDA_env_20µM_progression <- ggplot() +
  geom_point(
    data = dbRDA_bacteria_20µM[[1]],
    aes(
      x = dbRDA1,
      y = dbRDA2,
      color = `Progression.throught.the.bloom`,
      shape = Year
    ),
    size = 4
  ) +
  geom_segment(
    data = dbRDA_bacteria_20µM[[2]],
    aes(
      x = 0,
      xend = dbRDA1,
      y = 0,
      yend = dbRDA2
    ),
    arrow = arrow(length = unit(0.2, "cm")),
    linewidth = 1
  ) +
  ggrepel::geom_text_repel(data = dbRDA_bacteria_20µM[[2]],
                           aes(x = dbRDA1, y = dbRDA2, label = variable),
                           force = 1.5) +
  scale_color_viridis_c(option = "viridis", direction = -1) +
  scale_shape_manual(values = c(15, 16, 17), labels = as_labeller(label_year)) +
  theme_test() +
  my_theme +
  theme(legend.position = "bottom") +
  labs(x = xlabel_20µM, y = ylabel_20µM, color = "Progression through sampling")

### b) Bloom Stage -------

dbRDA_env_20µM_Group <- ggplot() +
  geom_point(data = dbRDA_bacteria_20µM[[1]],
             aes(
               x = dbRDA1,
               y = dbRDA2,
               color = (Group),
               shape = Year
             ),
             size = 4) +
  geom_segment(
    data = dbRDA_bacteria_20µM[[2]],
    aes(
      x = 0,
      xend = dbRDA1,
      y = 0,
      yend = dbRDA2
    ),
    arrow = arrow(length = unit(0.2, "cm")),
    linewidth = 1
  ) +
  ggrepel::geom_text_repel(data = dbRDA_bacteria_20µM[[2]],
                           aes(x = dbRDA1, y = dbRDA2, label = variable),
                           force = 1.5,
                           size = 4) +
  scale_color_manual(values = color_group, labels =
                       as_labeller(relabel_group)) +
  scale_shape_manual(values = c(15, 16, 17), labels = as_labeller(label_year)) +
  theme_test() +
  my_theme +
  theme(legend.position = "bottom") +
  labs(x = xlabel_20µM, y = ylabel_20µM, color = "Bloom Stage")


## 4. Export -------

### a) Data --------

write.table(
  res_dbrda_bacteria_20µM,
  paste(res_path, "Data_Figure_4_ANOVA.tsv", sep = "/")
)
write.table(
  res_dbrda_bacteria_20µM_terms,
  paste(res_path, "Data_Figure_4_ANOVA_terms.tsv", sep =
          "/")
)
write.table(
  res_dbrda_bacteria_20µM_axis,
  paste(res_path, "Data_Figure_4_ANOVA_axis.tsv", sep =
          "/")
)
write.table(
  contribution_dbRDA_20µM$Hier.part,
  paste(res_path, "Data_Figure_4_ANOVA_contribution.tsv", sep = "/")
)


### b) Plot --------

svg(paste(fig_path, "Figure_4.svg", sep = "/"))
dbRDA_env_20µM_Group
dev.off()


# VI. PERMANOVA ---------


permanova_bacteria_20µM <- adonis2(dist_bacteria_20µM ~ Group * Year,
                                   data = dbRDA_bacteria_20µM[[1]],
                                   by = "terms")

# anova(betadisper(
#   dist_bacteria_20µM,
#   sample_data(phyloseq_bacteria_filtered_CSS_20µM)$Group
# ))
# anova(betadisper(
#   dist_bacteria_20µM,
#   sample_data(phyloseq_bacteria_filtered_CSS_20µM)$Year
# ))

write.table(
  data.frame(permanova_bacteria_20µM),
  paste(res_path, "Data_Figure_4_PERMANOVA.tsv", sep = "/")
)

