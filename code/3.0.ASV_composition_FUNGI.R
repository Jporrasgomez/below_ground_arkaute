






rm(list = ls(all.names = TRUE))  
pacman::p_unload(pacman::p_loaded(), character.only = TRUE) 
pacman::p_load(dplyr,reshape2,tidyverse, tidyr, lubridate, ggplot2, ggpubr,
               ggrepel, gridExtra, patchwork, vegan, eulerr, ggthemes) 

#rm(list = ls(all.names = TRUE))
#pacman::p_load(dplyr, tidyverse, DT, viridis, ggrepel, codyn, vegan, eulerr, ggplot2, ggthemes, ggpubr, ggforce )#
##source("ASV_num/1.first_script.R"); rm(list = setdiff(ls(), c("fungi_ASV_sampling", "biomass_imp", "biomass_noimp")))


source("code/palettes_labels.R")
palette <- palette_CB
labels <- labels1


fungi_ASV_raw <- read.csv("data/fungi_ASV_raw.csv") %>%
  select(treatment, plot, sampling, ASV_num, abundance) 


#fungi_ASV_filter <- fungi_ASV_raw %>% 
#  filter(abundance > 0.001)


####### DIFFERENCES AT SAMPLING x TREATMENT LEVEL #############

fungi_ASV_sampling <- fungi_ASV_raw %>% 
  group_by(sampling, treatment, ASV_num) %>% 
  summarize(
    abundance_mean_sampling = mean(abundance)
  ) %>% 
  mutate(id = paste0(as.character(treatment), "/" , as.character(sampling)))

treats <- unique(fungi_ASV_sampling$treatment)

list1 <- list()
gglist1 <- list()

for(i in seq_along(treats)){
  
  
  list1[[i]] <- subset(fungi_ASV_sampling, treatment == treats[i])
  
  sp_wide <- list1[[i]] %>%
    pivot_wider(id_cols = sampling,
                names_from = ASV_num,
                values_from = abundance_mean_sampling,
                values_fill = list(abundance_mean_sampling = 0),
    ) %>% 
    column_to_rownames(var = "sampling") %>% 
    arrange(as.numeric(rownames(.)))
  
  
  #Relative abundance
  
  # Perform NMDS using Bray-Curtis distance
  nmds_res <- metaMDS(sp_wide, distance = "bray", k = 2, trymax = 250, maxit = 999) 
  
  # Extract NMDS sample scores
  nmds_samples <- as.data.frame(scores(nmds_res, display = "sites"))
  
  # Extract NMDS species scores (optional)
  nmds_species <- as.data.frame(scores(nmds_res, display = "species"))
  
  gglist1[[i]] <- ggplot() +
    #geom_text_repel(data = nmds_species %>% 
    #                  rownames_to_column(var = "sp"),
    #                aes(x = NMDS1, y = NMDS2, label = sp),
    #                color = "grey",
    #                max.overlaps = 30) +
    geom_point(data = nmds_samples %>% 
                 rownames_to_column(var = "sampling"),
               aes(x = NMDS1, y = NMDS2),
               size = 1.5) +
    geom_text_repel(data = nmds_samples %>% 
                      rownames_to_column(var = "sampling"),
                    aes(x = NMDS1, y = NMDS2, label = sampling),
                    max.overlaps = 9) +
    geom_path(data = nmds_samples %>% 
                rownames_to_column(var = "sampling"),
              aes(x = NMDS1, y = NMDS2)) +
    geom_hline(aes(yintercept = 0), color = "gray52", linetype = "dashed") +
    geom_vline(aes(xintercept = 0), color = "gray52", linetype = "dashed") +
    labs(title = paste("NMDS using Bray-Curtis:", treats[i], sep = " "),
         subtitle = paste0("Stress = ", round(nmds_res$stress, 3)),
         x = "NMDS1",
         y = "NMDS2")
  
  
}


ggarrange(
  gglist1[[1]],
  gglist1[[2]],
  gglist1[[3]],
  gglist1[[4]], 
  ncol = 2, nrow = 2)


#https://rpubs.com/CPEL/NMDS
#We will do so using the metaMDS function. When running an NMDS you will have to identify your distance matrix
#(island.spp_distmat), your distance metric which should match the distance metric from your distance matrix 
#(here we use "bray"), your selected number of dimensions ("k"), your max number of iterations (usually **"maxit = 999"**),
#and the maximum number of random starts (usually "trymax = 250"). You may need to do a couple of runs with different 
#trymax values, especially if you are working with community data with a lot pf 0’s. Finally,
#wascores is a method of calculating species scores, default is TRUE. Check the metaMDS help file for
#other options to further customize your ordination if necessary.

#As a rule of thumb literature has identified the following cut-off values for stress-level:
#  
#  Higher than 0.2 is poor (risks for false interpretation).
#0.1 - 0.2 is fair (some distances can be misleading for interpretation).
#0.05 - 0.1 is good (can be confident in inferences from plot).
#Less than 0.05 is excellent (this can be rare).



########### 1. VISUALIZATION: NMDS ##########################

i = 1
# 1: no abundance_mean_sampling transformation
# 2: Hellinger transformation of abundance_mean_sampling



sp_wide_sampling <- fungi_ASV_sampling %>%
  pivot_wider(id_cols = c(sampling, treatment, id),
              names_from = ASV_num,
              values_from = abundance_mean_sampling,
              values_fill = list(abundance_mean_sampling = 0)) %>% 
  ungroup() %>% 
  as_tibble()

# create a distance matrix 
abundance_mean_sampling_data_sampling <- sp_wide_sampling %>% 
  select(-treatment, -sampling, -id)

# Optional step in which I transform abundance_mean_sampling data sampling in Hellinger abundance_mean_samplings
abundance_mean_sampling_data_sampling_hellinger <- vegan::decostand(abundance_mean_sampling_data_sampling, method = "hellinger")

abundance_mean_sampling_list <- list(abundance_mean_sampling_data_sampling, abundance_mean_sampling_data_sampling_hellinger)
# Compute Bray-Curtis distance matrix

distance_matrix_sampling_bc <- vegan::vegdist(as.data.frame(abundance_mean_sampling_list[[i]]), method = "bray")

# Run NMDS (2 dimensions, 100 tries)
nmds_bc_sampling <- metaMDS(distance_matrix_sampling_bc, k = 2, trymax = 250, maxit = 999)

# Extract NMDS coordinates
nmds_df_sampling <- data.frame(
  NMDS1 = nmds_bc_sampling$points[, 1],
  NMDS2 = nmds_bc_sampling$points[, 2],
  treatment = sp_wide_sampling$treatment,
  sampling = sp_wide_sampling$sampling
)

# Arrange by sampling order
nmds_df_sampling <- nmds_df_sampling %>% arrange(sampling)




# Plot NMDS results using ggplot

ggnmds_alltreatments <- 
  ggplot(nmds_df_sampling, aes(x = NMDS1, y = NMDS2, color = treatment)) +
  
  
  stat_ellipse(geom = "polygon", aes(fill = treatment),
               alpha = 0.12, show.legend = FALSE, level = 0.68) + 
  
  geom_point(size = 2, show.legend = T) +
  
  geom_text_repel(aes(label = sampling),
                  max.overlaps = 8,
                  size = 4.5,
                  show.legend = F) +
  
  geom_hline(yintercept = 0, color = "gray52", linetype = "dashed") +
  
  geom_path(aes(group = treatment), linewidth = 0.5, alpha = 0.2) +
  
  geom_vline(xintercept = 0, color = "gray52", linetype = "dashed") +
  
  
  scale_color_manual(values = palette_CB, labels = labels, guide = "legend") +
  
  scale_fill_manual(values = palette_CB, guide = "none" ) +
  
  scale_shape_manual(values = point_shapes, guide = "none") +
  
  labs(
    #title = "NMDS Bray-Curtis: mean abundance_mean_sampling of species at sampling level",
    subtitle = paste0("Stress = ", round(nmds_bc_sampling$stress, 3)),
    x = "NMDS1", y = "NMDS2", color = " ") +
  theme1

print(ggnmds_alltreatments)



#ggsave("results/Plots/protofinal/species_composition_sampling_LABELS.png", plot = ggnmds_alltreatments, dpi = 300)





########### 1.2. STATISTICAL ANALYSIS: PERMANOVA ##########################

adonis_sampling <- adonis2(
  distance_matrix_sampling_bc ~ treatment,  # puedes agregar más variables si quieres
  data = sp_wide_sampling,                 # debe tener las variables explicativas
  permutations = 999,                      # número de permutaciones
  method = "bray"
)

# Mostrar resultados
print(adonis_sampling)
# Hay un efecto significativo del tratamiento sobre la composición de especies (p = 0.001). 
# El tratamiento explica aproximadamente el 33.7% de la variación en la composición.

bd <- betadisper(distance_matrix_sampling_bc, sp_wide_sampling$treatment)
anova(bd)
# El resultado de ANOVA para las dispersiónes dentro de grupos (tratamientos) es significativo (p = 0.0004).
# Esto significa que la variabilidad o dispersión dentro de al menos un grupo es diferente respecto a otros grupos.
permutest(bd) 
plot(bd)
boxplot(bd)

permutest(bd, pairwise = TRUE)
TukeyHSD(bd)

library(pairwiseAdonis)

# Ejecutamos las comparaciones por pares
pw_adonis <- pairwise.adonis(
  x           = distance_matrix_sampling_bc,                 # tu matriz de distancias Hellinger–Bray
  factors     = sp_wide_sampling$treatment,  # factor con los cuatro tratamientos
  perm        = 999,                         # número de permutaciones
  p.adjust.m  = "BH"                         # corrección de p por Benjamini–Hochberg
)

print(pw_adonis)







# Chat gpt:

##| Unlike Principal Coordinate Analysis (PCoA) or Principal Component Analysis (PCA),
##| where each axis has an explained variance, Non-Metric Multidimensional Scaling (NMDS)
##| does not provide a direct percentage of variance explained per axis. 
##| However, you can estimate how much each NMDS axis contributes to the representation of 
##| the distances by using the correlation between the NMDS axes and the original distance matrix.
##| You can calculate the squared correlation (R²) between the original distance matrix and each NMDS axis.
##|  This gives you an approximation of how well each axis represents the distances.

cor1 <- cor(distance_matrix_sampling_bc, dist(nmds_bc_sampling$points[,1]), method = "pearson") 
cor2 <- cor(distance_matrix_sampling_bc, dist(nmds_bc_sampling$points[,2]), method = "pearson") 

# Compute percentage explained by each axis
explained_NMDS1 <- cor1^2 / (cor1^2 + cor2^2) * 100
explained_NMDS2 <- cor2^2 / (cor1^2 + cor2^2) * 100

# Print results
cat("NMDS1 explains:", round(explained_NMDS1, 2), "%\n")
cat("NMDS2 explains:", round(explained_NMDS2, 2), "%\n")



















































list_sorensen <- list()
list_sorensen_df <- list()
samps <- sort(unique(fungi_ASV_sampling$sampling))

for(i in 1:length(samps)){
  
  sp_wide <- fungi_ASV_sampling %>% 
    filter(sampling == samps[i]) %>% 
    pivot_wider(id_cols = c(sampling, treatment),
                names_from = ASV_num,
                values_from = abundance_mean_sampling,
                values_fill = list(abundance_mean_sampling = 0)) %>% 
    ungroup() %>% 
    as_tibble()
  
  
  abundance_mean_sampling_data <- sp_wide %>%
    select(-treatment, -sampling)
  
  sorensen <- vegdist(abundance_mean_sampling_data, method = "bray", binary = TRUE)
  sorensen <- as.matrix(sorensen)  
  rownames(sorensen) <- sp_wide$id
  colnames(sorensen) <- sp_wide$id
  sorensen[upper.tri(sorensen)] <- NA
  
  sorensen_df <- sorensen %>% 
    as.data.frame() %>%
    rownames_to_column(var = "row_name") %>%
    pivot_longer(-row_name, names_to = "col_name", values_to = "value") %>% 
    filter(!is.na(value))
  
  list_sorensen[[i]] <- sorensen
  list_sorensen_df[[i]] <- sorensen_df
  
  
}


print(list_sorensen)


sorensen_df <- bind_rows(list_sorensen_df) %>% 
  filter(!value == "0") %>%
  separate(col_name, into = c("treatment_x", "sampling_x"), sep = "/") %>%
  separate(row_name, into = c("treatment_y", "sampling_y"), sep = "/") %>% 
  select(-sampling_x) %>% 
  rename(sampling = sampling_y) %>% 
  mutate(sampling = factor(as.numeric(sampling), levels = sort(unique(as.numeric(sampling))))) %>% 
  arrange(sampling) %>% 
  mutate(comparison = paste0(treatment_x, "-", treatment_y)) %>% 
  select(-treatment_x, -treatment_y) %>% 
  mutate(comparison = ifelse(comparison == "p-c", "c-p", comparison),
         comparison = ifelse(comparison == "w-c", "c-w", comparison),
         comparison = ifelse(comparison == "wp-c", "c-wp", comparison))

{gg_sorensen <- 
    sorensen_df %>% 
    filter(comparison %in% c("c-p", "c-w", "c-wp")) %>% 
    ggplot(aes(x = sampling, y = value, color = comparison, group = comparison)) + 
    facet_wrap(~ comparison, labeller = labeller(
      RR_descriptor = as_labeller(labels_RR)), ncol = 1) + 
    geom_point() + 
    geom_line() +
    geom_smooth(
      se = TRUE, aes(color = comparison, fill = comparison),
      method = "lm", span = 0.6, alpha = 0.2 ) + 
    scale_color_manual(values = c("c-p" = "#0077FF", "c-w" = "#E0352F", "c-wp" = "#A238A2")) +
    scale_fill_manual(values = c("c-p" = "#0077FF", "c-w" = "#E0352F", "c-wp" = "#A238A2")) +
    labs( y = "Beta-diversity Sorensen") + 
    scale_x_discrete(breaks = function(x) x[seq(1, length(x), by = 2)]) + 
    labs(x = "Sampling") + 
    theme1+
    theme(legend.position = "bottom")
  print(gg_sorensen)
  #ggsave("results/Plots/protofinal/sorensen_c.png", plot = gg_sorensen, dpi = 300)
  }


{gg_sorensen_wp <- 
    sorensen_df %>% 
    filter(comparison %in% c("p-wp", "w-wp")) %>% 
    ggplot(aes(x = sampling, y = value, color = comparison, group = comparison)) + 
    facet_wrap(~ comparison) + 
    geom_point() + 
    geom_line() +
    geom_smooth(
      se = TRUE, aes(color = comparison, fill = comparison),
      method = "lm", span = 0.6, alpha = 0.2 ) + 
    scale_color_manual(values = c("w-wp" = "#D08A00", "p-wp" = "#3A3A3A")) +
    scale_fill_manual(values = c("w-wp" = "#D08A00", "p-wp" = "#3A3A3A")) +
    labs( y = "Beta-diversity Sorensen") + 
    scale_x_discrete(breaks = function(x) x[seq(1, length(x), by = 2)]) + 
    labs(x = "Sampling") + 
    theme(legend.position = "bottom")
  print(gg_sorensen_wp)
  #ggsave("results/Plots/protofinal/sorensen_wp.png", plot = gg_sorensen_wp, dpi = 300)
}



