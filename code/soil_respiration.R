

rm(list = ls(all.names = TRUE))  
pacman::p_unload(pacman::p_loaded(), character.only = TRUE) 
pacman::p_load(dplyr,reshape2,tidyverse, tidyr, readr, lubridate, ggplot2, ggpubr, gridExtra, readxl, patchwork) 

source("code/palettes_labels.R")

plots <- read.csv("data/plots.csv") %>% 
  select(plot, treatment_code) %>% 
  rename(treatment = treatment_code)



sampling_dates <- read.csv("data/sampling_dates.csv") %>%
  mutate(
    date              = ymd(date), 
    day               = day(date), 
    month             = month(date, label = TRUE),
    year              = year(date),
    label_micro       = as.factor(label_micro),
    date_label        = factor(format(date, "%d-%b-%y"), 
                               levels = format(sort(unique(date)), "%d-%b-%y"), 
                               ordered = TRUE),
    
    date_label_noyear = factor(format(date, "%b %d"), 
                               levels = format(sort(unique(date)), "%b %d"), 
                               ordered = TRUE),
  ) %>% 
  
  select(date, day, month, year, label_micro, date_label, date_label_noyear)
  #filter(!is.na(label_micro))

sampling_dates$sampling <- as.numeric(row.names(sampling_dates)) - 1




variable_labels <- c("Mtype", "date", "time", "plot", "rec_no", "co2_ref", "AP", "flow", "h20_ref", 
                     "RHT", "o2_ref", "error_code", "aux_voltage", "PAR", "temp_soil", "temp_air", 
                     "soil_moisture", "process", "DC", "DT", "SRL_rate", "SRQ_rate")


vector_files <- list.files(path = "data/soil_respiration", pattern = ".TXT", full.names = TRUE)
sampling_labels <- c("0", "1", "10", "11", "13", "15", "16", "18", "19", 2:9)
  
list_result <- list()
list_plots <- list()

for(i in seq_along(vector_files)){
  
  lines <- read_lines(vector_files[i]) 
  
  clean_lines <- lines[!stringr::str_detect(lines, "Zero|Start|End")]
  
  list_result[[i]] <- read_csv(I(clean_lines), col_names = variable_labels) %>% 
    rename(date_EGM = date) %>% 
    mutate(sampling = as.numeric(paste0(sampling_labels[i]))) %>% 
    left_join(plots) %>% 
    left_join(sampling_dates[sampling_dates$sampling == sampling_labels[i], ])
  
  
  list_plots[[i]] <- 
   list_result[[i]] %>%
    ggplot() +
    facet_wrap(~plot, scale = "free")+
    geom_point(aes(x = time, y = SRL_rate), color = "black", alpha = 0.5)+
    geom_point(aes(x = time, y = SRQ_rate), color = "blue", alpha = 0.5)+
    labs(title = paste0(sampling_labels[i])) +
    theme(legend.position = "bottom")
 
  
}



list_plots[[1]]
list_plots[[2]]
list_plots[[3]] # Review this
list_plots[[4]]
list_plots[[5]]
list_plots[[6]]
list_plots[[7]]
list_plots[[8]] # Review this
list_plots[[9]]
list_plots[[10]]
list_plots[[11]]
list_plots[[12]]
list_plots[[13]]
list_plots[[14]] # Review this
list_plots[[15]]
list_plots[[16]]
list_plots[[17]]


soil_respiration_raw <- do.call(rbind, list_result)

soil_respiration_point <- soil_respiration_raw %>% 
  group_by(sampling, date, date_EGM, date_label_noyear, treatment, plot) %>% 
  filter(time == max(time)) %>% 
  reframe(
    time = time,
    label_micro = label_micro,
    SRL_rate = SRL_rate,
    SRQ_rate = SRQ_rate
  )

  

soil_respiration_point %>% 
  #filter(!is.na(label_micro)) %>% 
  ggplot(aes(x = treatment, y = SRL_rate, fill = treatment)) +
  geom_boxplot() +
  scale_fill_manual(values = palette_CB) + 
  theme_otc

soil_respiration_point %>% 
  #filter(!is.na(label_micro)) %>% 
  group_by(sampling, treatment) %>% 
  mutate(
    SRL_rate_mean = mean(SRL_rate),
    SRL_rate_sd = sd(SRL_rate)
  ) %>% 
  ggplot() +
  geom_point(aes(x = date, y = SRL_rate, color = treatment)) + 
  geom_point(aes(x = date, y = SRL_rate_mean, color = treatment), size = 3) + 
  geom_line(aes(x = date, y = SRL_rate_mean, color = treatment)) + 
  scale_color_manual(values = palette_CB) 
   

soil_respiration_micro <- soil_respiration_point %>% 
  filter(!is.na(label_micro))

soil_respiration_micro %>%  write.csv("data/soil_respiration_micro.csv", row.names = FALSE)


source("code/functions/eff_size_LRR_function.R")
source("code/functions/eff_size_dynamics_LRR_function.R")
source("code/functions/gg_aggregated_function_2.R")
source("code/functions/gg_dynamics_function2.R")



{

SR_no0 <- soil_respiration_point %>% filter(sampling != "0")
SR_dyn <- soil_respiration_point



LRR_agg(SR_no0, "SRL_rate")

agg <- effsize_data %>% 
  mutate(
    eff_value = round(eff_value, 2),
    lower_limit = round(lower_limit, 2),
    upper_limit = round(upper_limit, 2)
  ) %>% 
  select(eff_descriptor, variable, eff_value, lower_limit, upper_limit, null_effect)



LRR_dynamics(SR_dyn, "SRL_rate")

dyn <- effsize_dynamics_data %>% 
  mutate(
    date_label_noyear = factor(
      date_label_noyear,
      levels = unique(date_label_noyear[order(date)]),
      ordered = TRUE
    )
  )



comparissons <- c("p_vs_c", "w_vs_c", "wp_vs_c")

lvls <-"SRL_rate"
labs <- c("SRL_rate" = "Soil respiration rate")


gg_eff_agg_c2 <- agg %>% 
  filter(eff_descriptor %in% comparissons,
         variable == "SRL_rate") %>% 
  mutate(
    eff_descriptor = factor(eff_descriptor, levels = comparissons),
    variable       = factor(variable, levels = lvls,
                            labels = labs)
  ) %>% 
  ggagg2(
    palette   = palette_RR_CB,
    labels    = labels_RR2,
    colorline = "grey50",
    limitvar  = lvls,
    labelvar  = labs, 
    breaks_axix_y = 4
  )


gg_eff_dynamics_c2<- dyn %>% 
  filter(eff_descriptor %in% comparissons) %>% 
  filter(variable == "SRL_rate") %>%  
  mutate(
    variable       = factor(variable, levels = lvls,
                            labels = labs)
    ) %>% 
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
                 widths = c(1, 10))) +
  plot_annotation(theme = theme(legend.position = "bottom"))



########### COMBINED / PERTURBATION    ###

# <- position_dodge2(width = 0.1, preserve = "single")
#pos_dod_wp_dyn <- position_dodge2(width = 4, preserve = "single")



gg_eff_agg_wp2 <- agg %>% 
  filter(eff_descriptor == "wp_vs_p",
         variable == "SRL_rate") %>% 
  mutate(
    variable       = factor(variable, levels = lvls,
                            labels = labs)
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
  filter(variable == "SRL_rate") %>%  
  mutate(
    variable = factor(variable, levels = lvls,
                                       labels = labs)
    ) %>% 
  
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
                 widths = c(1, 10))) +
  plot_annotation(theme = theme(legend.position = "bottom"))


}



gg_control
gg_wp


ggsave("results/plots/soil_respiration_c.png", plot = gg_control, dpi = 300)
ggsave("results/plots/soil_respiration_wp.png", plot = gg_wp, dpi = 300)







