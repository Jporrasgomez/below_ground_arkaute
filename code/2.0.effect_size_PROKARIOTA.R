



rm(list = ls(all.names = TRUE))  
pacman::p_unload(pacman::p_loaded(), character.only = TRUE) 
pacman::p_load(dplyr,reshape2,tidyverse, tidyr, lubridate, ggplot2, ggpubr, gridExtra, readxl, patchwork) 

source("code/palettes_labels.R")

palette <- palette_CB

prokariota_data <- read.csv("data/data_prokariota.csv") %>%  select(-X)

prokariota_data %>% 
  pivot_longer(
    cols = c(-marker, -date, -date_label_noyear, -year, -sampling, -plot, -treatment), 
    names_to = "variable", 
    values_to = "values"
  ) %>% 
  ggplot(aes(y = values, x = "")) + 
  facet_wrap(~ variable, scale = "free") + 
  geom_boxplot()

# There are outliers. However, each dot represent a replicate. For dynamic analysis
# we cannot lose points because we only have 4. 




# Geary test modification proposed in Lajeunesse 2015. We need to test
# how suitable the data is for log response ratio analysis and their delta or sigma
# modifications (also proposed in the same paper).If the test is below 3, data is not suitable

geary_test_treatment <- prokariota_data %>%
  pivot_longer(
    cols = c(-marker, -date, -date_label_noyear, -year, -sampling, -plot, -treatment),
    values_to = "value", 
    names_to = "variable"
  ) %>% 
  filter(sampling != "1") %>% 
  group_by(treatment, variable) %>% 
  summarize(
    n = n(),
    mean_variable = mean(value, na.rm = T), 
    sd_variable = sd(value, na.rm = T)
  ) %>% 
  mutate(geary_test_value = (mean_variable/sd_variable) *((4 * n^1.5) / (1 + 4 * n))) %>% 
  mutate(geary_test_outcome = ifelse(geary_test_value >= 3, paste0("TRUE"), paste0("FALSE")))


geary_test_sampling <- prokariota_data %>%
  pivot_longer(
    cols = c(-marker, -date, -date_label_noyear, -year, -sampling, -plot, -treatment),
    values_to = "value", 
    names_to = "variable"
  ) %>% 
  group_by(treatment, sampling, variable) %>% 
  summarize(
    n = n(),
    mean_variable = mean(value, na.rm = T), 
    sd_variable = sd(value, na.rm = T)
  ) %>% 
  mutate(geary_test_value = (mean_variable/sd_variable) *((4 * n^1.5) / (1 + 4 * n))) %>% 
  mutate(geary_test_outcome = ifelse(geary_test_value >= 3, paste0("TRUE"), paste0("FALSE")))

false_sampling <- geary_test_sampling %>% 
  filter(geary_test_outcome == "FALSE")
unique(false_sampling$variable)

nrow(false_sampling)/nrow(geary_test_sampling) *100

false_cases_sampling <- false_sampling %>%
  ungroup() %>% 
  select(variable, treatment, sampling, geary_test_value, geary_test_outcome) %>% 
  group_by(variable, treatment) %>% 
  summarize(
    n = n(),
    samplings = list(as.numeric(sampling)),
    .groups = "drop"
  )



  

############ EFFECT SIZE #############



{
  variables_prokariota <- 
  (prokariota_data %>% 
         select(-plot , -marker, -treatment, -date, -date_label_noyear, -year, -sampling, -goods_coverage) %>% 
         colnames())
prokariota_no0 <- prokariota_data %>% 
  filter( sampling != "0") 
prokariota <- prokariota_data


limits_main_variables <- variables_prokariota
labels_main_variables <- c(
  "nASV"                                          = "Number ASV", 
  "observed_features"                             = "Number ASV novogene", 
  "chao1"                                         = "Estimated richness (chao1)",
  "dominance"                                     = "Dominance", 
  "pielou_e"                                      = "Evenness", 
  "shannon"                                       = "Shannon", 
  "simpson"                                       = "Simpson", 
  "aerobic_ammonia_oxidation"                     = "Aer. Ammonia Ox.",
  "aerobic_chemoheterotrophy"                     = "Aer. Chemoheterotrophy", 
  "aliphatic_non_methane_hydrocarbon_degradation" = "ANMHD", 
  "animal_parasites_or_symbionts"                 = "Animal parasites/symbionts", 
  "aromatic_compound_degradation"                 = "Aromatic comp. deg.", 
  "aromatic_hydrocarbon_degradation"              = "Aromatic HC deg.",
  "cellulolysis"                                  = "Cellulolysis",
  "chemoheterotrophy"                             = "Chemoheterotrophy",
  "chitinolysis"                                  = "Chitinolysis",
  "denitrification"                               = "Denitrification",
  "fermentation"                                  = "Fermentation",
  "human_pathogens_all"                           = "Human Pathogens",
  "hydrocarbon_degradation"                       = "HC Degradation",
  "iron_respiration"                              = "Iron Respiration",
  "nitrate_denitrification"                       = "Nitrate Denitrif.",
  "nitrate_reduction"                             = "Nitrate Reduction",
  "nitrate_respiration"                           = "Nitrate Respir.",
  "nitrification"                                 = "Nitrification",
  "nitrite_denitrification"                       = "Nitrite Denitrif.",
  "nitrite_respiration"                           = "Nitrite Respir.",
  "nitrogen_respiration"                          = "Nitrogen Respir.",
  "nitrous_oxide_denitrification"                 = "N2O Denitrif.",
  "Others"                                        = "Others",
  "predatory_or_exoparasitic"                     = "Predatory/Exoparasitic",
  "respiration_of_sulfur_compounds"               = "Sulfur Comp. Respir.",
  "sulfur_respiration"                            = "Sulfur Respir.",
  "ureolysis"                                     = "Ureolysis",
  "xylanolysis"                                   = "Xylanolysis"
)




source("code/functions/new_aggregated.R")
source("code/functions/new_dynamics.R")
source("code/functions/gg_aggregated_function_2.R")
source("code/functions/gg_dynamics_function2.R")



list_prokariota <- list(prokariota_no0, prokariota)
comparissons <- c("p_vs_c", "w_vs_c", "wp_vs_c")

list_agg <- list()
list_dyn <- list()
for(i in seq_along(variables_prokariota)){ 
    

    LRR_agg(list_prokariota[[1]], variables_prokariota[i])
    
    list_agg[[i]] <- effsize_data %>% 
      #mutate(
      #  eff_value = round(eff_value, 2),
      #  lower_limit = round(lower_limit, 2),
      #  upper_limit = round(upper_limit, 2)
      #) %>% 
      select(eff_descriptor, variable, eff_value, lower_limit, upper_limit, null_effect, analysis)
    
  
    
    LRR_dynamics(list_prokariota[[2]], variables_prokariota[i])
    
    list_dyn[[i]] <- effsize_dynamics_data
    
}

agg <- do.call(rbind, list_agg)
dyn <- do.call(rbind, list_dyn) %>% 
  mutate(
    date_label_noyear = factor(
      date_label_noyear,
      levels = unique(date_label_noyear[order(date)]),
      ordered = TRUE
    )
  )


   
list_results <- list()
results_prokariota <- list()


    lvls1 <- limits_main_variables[1:6]
    labs1 <- unname(labels_main_variables[1:6])
    
    lvls2 <- limits_main_variables[7:12]
    labs2 <- unname(labels_main_variables[7:12])
    
    lvls3 <- limits_main_variables[13:18]
    labs3 <- unname(labels_main_variables[13:18])
    
    lvls4 <- limits_main_variables[19:24]
    labs4 <- unname(labels_main_variables[19:24])
    
    lvls5 <- limits_main_variables[25:30]
    labs5 <- unname(labels_main_variables[25:30])
    
    lvls6 <- limits_main_variables[31:35]
    labs6 <- unname(labels_main_variables[31:35])
    
    
list_levels <- list(lvls1, lvls2, lvls3, lvls4, lvls5, lvls6)
list_labs   <- list(labs1, labs2, labs3, labs4, labs5, labs6)
    
list_results_c <- list()
list_results_wp <- list()
  

  for (i in seq_along(1:6)){
    
    gg_eff_agg_c2 <- agg %>% 
      filter(eff_descriptor %in% comparissons,
             variable %in% list_levels[[i]]) %>% 
      mutate(
        eff_descriptor = factor(eff_descriptor, levels = comparissons),
        variable       = factor(variable, levels = list_levels[[i]],
                                labels = list_labs[[i]])
      ) %>% 
      ggagg2(
        palette   = palette_RR_CB,
        labels    = labels_RR2,
        colorline = "grey50",
        limitvar  = lvls,
        labelvar  = labs, 
        breaks_axix_y = 3
      )
    
    
    gg_eff_dynamics_c2<- dyn %>% 
      filter(eff_descriptor %in% comparissons) %>% 
      filter(variable %in% list_levels[[i]]) %>%  
      mutate(
        variable = factor(variable, 
                          levels = list_levels[[i]], 
                          labels = list_labs[[i]])) %>% 
      ggdyn2(palette_RR_CB,
             labels_RR2, 
             "grey50",
             position = position_dodge(width = 0.5),
             asterisk = 8, 
             caps = position_dodge(width = 0.5)$width)
    
    
    gg_control <-
      (gg_eff_agg_c2 + 
         gg_eff_dynamics_c2 + theme (legend.position = "none") + 
         plot_layout(guides = "collect",
                     widths = c(1, 8))) +
      plot_annotation(theme = theme(legend.position = "bottom"))
    
    list_results_c[[i]] <- gg_control
    
    
    
    
    ########### COMBINED / PERTURBATION    ###
    
    pos_dod_wp_agg <- position_dodge2(width = 0.1, preserve = "single")
    pos_dod_wp_dyn <- position_dodge2(width = 4, preserve = "single")
    
    
    
    gg_eff_agg_wp2 <- agg %>% 
      filter(eff_descriptor == "wp_vs_p",
             variable %in% list_levels[[i]]) %>% 
      mutate(
        variable = factor(variable, levels = list_levels[[i]],
                          labels = list_labs[[i]])
      ) %>% 
      ggagg2(
        palette   = palette_RR_wp,
        labels    = labels_RR_wp,
        colorline = p_CB,
        limitvar  = lvls,
        labelvar  = labs, 
        breaks_axix_y = 2
      )
    
    
    
    
    gg_eff_dynamics_wp2<- dyn %>% 
      filter(eff_descriptor %in% c("wp_vs_p")) %>% 
      filter(variable %in% list_levels[[i]]) %>%  
      mutate(
        variable = factor(variable, 
                          levels = list_levels[[i]], 
                          labels = list_labs[[i]])) %>% 
      
      ggdyn2(palette_RR_wp,
             labels_RR_wp2, 
             p_CB,
             position = position_dodge(width = 0.5),
             asterisk = 8, 
             caps = position_dodge(width = 0.5)$width)
    
    
    gg_wp <-
      (gg_eff_agg_wp2 + 
         gg_eff_dynamics_wp2 + theme (legend.position = "none") + 
         plot_layout(guides = "collect",
                     widths = c(1, 8))) +
      plot_annotation(theme = theme(legend.position = "bottom"))
    
    
    list_results_wp[[i]] <- gg_wp
    
  
  }
}


list_results_c[[1]]
list_results_c[[2]]
list_results_c[[3]]
list_results_c[[4]]
list_results_c[[5]]
list_results_c[[6]]


list_results_wp[[1]]
list_results_wp[[2]]
list_results_wp[[3]]
list_results_wp[[4]]
list_results_wp[[5]]
list_results_wp[[6]]




ggsave("results/plots/prokariota_control_1.png", plot = list_results_c[[1]], dpi = 300)
ggsave("results/plots/prokariota_control_2.png", plot = list_results_c[[2]], dpi = 300)
ggsave("results/plots/prokariota_control_3.png", plot = list_results_c[[3]], dpi = 300)
ggsave("results/plots/prokariota_control_4.png", plot = list_results_c[[4]], dpi = 300)
ggsave("results/plots/prokariota_control_5.png", plot = list_results_c[[5]], dpi = 300)
ggsave("results/plots/prokariota_control_6.png", plot = list_results_c[[6]], dpi = 300)


ggsave("results/plots/prokariota_wp_1.png", plot = list_results_wp[[1]], dpi = 300)
ggsave("results/plots/prokariota_wp_2.png", plot = list_results_wp[[2]], dpi = 300)
ggsave("results/plots/prokariota_wp_3.png", plot = list_results_wp[[3]], dpi = 300)
ggsave("results/plots/prokariota_wp_4.png", plot = list_results_wp[[4]], dpi = 300)
ggsave("results/plots/prokariota_wp_5.png", plot = list_results_wp[[5]], dpi = 300)
ggsave("results/plots/prokariota_wp_6.png", plot = list_results_wp[[6]], dpi = 300)







