



rm(list = ls(all.names = TRUE))  
pacman::p_unload(pacman::p_loaded(), character.only = TRUE) 
pacman::p_load(dplyr,reshape2,tidyverse, tidyr, lubridate, ggplot2, ggpubr, gridExtra, readxl, patchwork) 

source("code/palettes_labels.R")

palette <- palette_CB

temp_vwc_data <- read.csv("data/processed_data/temp_vwc_data.csv")
prokariota_data <- read.csv("data/data_prokariota.csv") %>%  select(-X)
fungi_data <- read.csv("data/data_fungi.csv") %>%  select(-X)
SR_data <- read.csv("data/soil_respiration_micro.csv") %>% 
  select(-sampling, -date)


variables_prokariota <- 
  (prokariota_data %>% 
     select(-plot , -marker, -treatment, -date, -date_label_noyear, -year, -sampling, -goods_coverage) %>% 
     colnames())

variables_fungi <- 
  (fungi_data %>% 
     select(-plot , -marker, -treatment, -date, -date_label_noyear, -year, -sampling, -goods_coverage) %>% 
     colnames())



prokariota_SR <- prokariota_data %>% 
  right_join(SR_data) %>% 
  left_join(temp_vwc_data)


fungi_SR <- fungi_data %>% 
  right_join(SR_data) %>% 
  left_join(temp_vwc_data)




  fungi_SR %>% 
  ggplot(aes(x = SRL_rate, y = mean_t_ground)) +  
  geom_point(size = 2) + 
  labs(x = "SRL_rate", y = "mean_t_ground") + 
  geom_smooth(method = "lm", se = FALSE) +
  
  stat_regline_equation(
    mapping     = aes(label = after_stat(eq.label)),
    formula     = y ~ x,
    label.x.npc = 0.2,
    label.y.npc = 0.95,
    size        = 5,
    show.legend = FALSE
  ) +
  
  stat_regline_equation(
    mapping     = aes(label = after_stat(rr.label)),
    formula     = y ~ x,
    label.x.npc = 0.2,
    label.y.npc = 0.80,
    size        = 5,
    show.legend = FALSE
  ) +
  
  stat_cor(
    mapping     = aes(label = after_stat(p.label)),
    method      = "pearson",      
    label.x.npc = 0.2,           
    label.y.npc = 0.60,           
    size        = 5,
    show.legend = FALSE
  ) +
  
  theme_otc






library(purrr)


list_prokariota <- list()

for(i in seq_along(variables_prokariota)){
  
  
  var_actual <- variables_prokariota[i]
  
  list_prokariota[[i]]  <- 
    prokariota_SR %>% 
  
    ggplot(aes(x = SRL_rate, y = .data[[var_actual]])) +  
    geom_point(size = 2) + 
    labs(x = "SRL_rate", y = var_actual) + 
    geom_smooth(method = "lm", se = FALSE) +
    
    stat_regline_equation(
      mapping     = aes(label = after_stat(eq.label)),
      formula     = y ~ x,
      label.x.npc = 0.2,
      label.y.npc = 0.95,
      size        = 5,
      show.legend = FALSE
    ) +
    
    stat_regline_equation(
      mapping     = aes(label = after_stat(rr.label)),
      formula     = y ~ x,
      label.x.npc = 0.2,
      label.y.npc = 0.80,
      size        = 5,
      show.legend = FALSE
    ) +
    
    stat_cor(
      mapping     = aes(label = after_stat(p.label)),
      method      = "pearson",      
      label.x.npc = 0.2,           
      label.y.npc = 0.60,           
      size        = 5,
      show.legend = FALSE
    ) +
    
    theme_otc
}



# Sustituye 'list_result' por el nombre real de tu lista de gráficos
walk(list_prokariota, print)











list_fungi <- list()

for(i in seq_along(variables_fungi)){
  
  
  var_actual <- variables_fungi[i]
  
  list_fungi[[i]]  <- 
    fungi_SR %>% 
    
    ggplot(aes(x = SRL_rate, y = .data[[var_actual]])) +  
    geom_point(size = 2) + 
    labs(x = "SRL_rate", y = var_actual) + 
    geom_smooth(method = "lm", se = FALSE) +
    
    stat_regline_equation(
      mapping     = aes(label = after_stat(eq.label)),
      formula     = y ~ x,
      label.x.npc = 0.2,
      label.y.npc = 0.95,
      size        = 5,
      show.legend = FALSE
    ) +
    
    stat_regline_equation(
      mapping     = aes(label = after_stat(rr.label)),
      formula     = y ~ x,
      label.x.npc = 0.2,
      label.y.npc = 0.80,
      size        = 5,
      show.legend = FALSE
    ) +
    
    stat_cor(
      mapping     = aes(label = after_stat(p.label)),
      method      = "pearson",      
      label.x.npc = 0.2,           
      label.y.npc = 0.60,           
      size        = 5,
      show.legend = FALSE
    ) +
    
    theme_otc
}


# Sustituye 'list_result' por el nombre real de tu lista de gráficos
walk(list_fungi, print)





















list_prokariota_temp <- list()

for(i in seq_along(variables_prokariota)){
  
  
  var_actual <- variables_prokariota[i]
  
  list_prokariota_temp[[i]]  <- 
    prokariota_SR %>% 
    
    ggplot(aes(x = mean_t_ground, y = .data[[var_actual]])) +  
    geom_point(size = 2) + 
    labs(x = "Mean ground temperature", y = var_actual) + 
    geom_smooth(method = "lm", se = FALSE) +
    
    stat_regline_equation(
      mapping     = aes(label = after_stat(eq.label)),
      formula     = y ~ x,
      label.x.npc = 0.2,
      label.y.npc = 0.95,
      size        = 5,
      show.legend = FALSE
    ) +
    
    stat_regline_equation(
      mapping     = aes(label = after_stat(rr.label)),
      formula     = y ~ x,
      label.x.npc = 0.2,
      label.y.npc = 0.80,
      size        = 5,
      show.legend = FALSE
    ) +
    
    stat_cor(
      mapping     = aes(label = after_stat(p.label)),
      method      = "pearson",      
      label.x.npc = 0.2,           
      label.y.npc = 0.60,           
      size        = 5,
      show.legend = FALSE
    ) +
    
    theme_otc
}

# Sustituye 'list_result' por el nombre real de tu lista de gráficos
walk(list_prokariota_temp, print)







list_fungi_temp <- list()

for(i in seq_along(variables_fungi)){
  
  
  var_actual <- variables_fungi[i]
  
  list_fungi_temp[[i]]  <- 
    fungi_SR %>% 
    
    ggplot(aes(x = mean_t_ground, y = .data[[var_actual]])) +  
    geom_point(size = 2) + 
    labs(x = "Mean ground temperature", y = var_actual) + 
    geom_smooth(method = "lm", se = FALSE) +
    
    stat_regline_equation(
      mapping     = aes(label = after_stat(eq.label)),
      formula     = y ~ x,
      label.x.npc = 0.2,
      label.y.npc = 0.95,
      size        = 5,
      show.legend = FALSE
    ) +
    
    stat_regline_equation(
      mapping     = aes(label = after_stat(rr.label)),
      formula     = y ~ x,
      label.x.npc = 0.2,
      label.y.npc = 0.80,
      size        = 5,
      show.legend = FALSE
    ) +
    
    stat_cor(
      mapping     = aes(label = after_stat(p.label)),
      method      = "pearson",      
      label.x.npc = 0.2,           
      label.y.npc = 0.60,           
      size        = 5,
      show.legend = FALSE
    ) +
    
    theme_otc
}


# Sustituye 'list_result' por el nombre real de tu lista de gráficos
walk(list_fungi_temp, print)


