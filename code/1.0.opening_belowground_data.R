

# HACER UNA BASE DE DATOS ÚNICA PARA FUNGI CON TODAS LAS VARIABLES Y OTRA PARA PROKARIOTA



rm(list = ls(all.names = TRUE))  
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)
pacman::p_load(dplyr,reshape2,tidyverse, tidyr, lubridate, ggplot2, ggpubr, gridExtra, readxl, readr) 


# Sampling dates and plots


labels_micro <- read.csv("data/labels_micro.csv") %>%
  select(plot, dna_epp_FINAL_LABEL_BACTERIA_BATCH, dna_epp_FINAL_LABEL_FUNGHI_BATCH) %>% 
  rename(sample_prokariota = dna_epp_FINAL_LABEL_BACTERIA_BATCH, 
         sample_fungi = dna_epp_FINAL_LABEL_FUNGHI_BATCH)

info_plots <- read.csv("data/plots.csv") %>% 
  select(treatment_code, plot) %>% 
  rename(treatment = treatment_code)




sampling_dates <- read.csv("data/sampling_dates.csv") %>%
  filter(!is.na(label_micro)) %>% 
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
                               ordered = TRUE)
  ) %>% 
  
  select(date, day, month, year, label_micro, date_label, date_label_noyear)





info_micro <- merge(labels_micro, info_plots) %>% 
  mutate(
    label_micro = substr(sample_prokariota, 1, 1)) %>% 
  merge(sampling_dates, by = "label_micro") %>% 
  mutate(label_micro = as.factor(label_micro)) %>% 
  mutate(sampling = as.integer(label_micro))






info_micro %>%  write.csv("data/metadata.csv")


info_fungi <- info_micro %>% 
  select(-sample_prokariota) %>%
  rename(sample = sample_fungi)

info_prokariota <- info_micro %>% 
  select(-sample_fungi) %>%
  rename(sample = sample_prokariota)




##### 1.  Opening and arranging data ABUNDANCE #####


  prokariota_raw1 <- read.delim(
  "data/16S_1st_batch_relative_abundance/featureTable.sample.total.relative.xls",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
) 
colnames(prokariota_raw1)[1] <- "ASV_num"

#abundancias <- prokariota_raw1[, sapply(prokariota_raw1, is.numeric)]
#sum(abundancias > 0)

prokariota_raw2 <- read.delim(
  "data/16S_2nd_batch_relative_abundance/featureTable.sample.total.relative.xls",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
) 

colnames(prokariota_raw2)[1] <- "ASV_num"

fungi_raw1 <- read.delim(
  "data/ITS_1st_batch_relative_abundance/featureTable.sample.total.relative.xls",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

colnames(fungi_raw1)[1] <- "ASV_num"

fungi_raw2 <- read.delim(
  "data/ITS_2nd_batch_relative_abundance/featureTable.sample.total.relative.xls",
  header = TRUE,
  sep = "\t",
  check.names = FALSE)

  colnames(fungi_raw2)[1] <- "ASV_num"

  
list_data <- list(prokariota_raw1, prokariota_raw2, fungi_raw1, fungi_raw2)
list_merge <- list(info_prokariota, info_prokariota, info_fungi, info_fungi)
markers <- c("16S", "16S", "ITS", "ITS")
list_result <- list()

for(i in seq_along(list_data)){
  
  list_result[[i]] <- list_data[[i]] %>% 
    pivot_longer(
      cols = -c(Taxonomy, ASV_num),   # Keep Taxonomy and X.OTU_num unchanged
      names_to = "sample",              # Name for the new "variable" column
      values_to = "abundance"               # Name for the new "value" column
    )  %>%
    separate(
      Taxonomy,
      into = c("k", "p", "c", "o", "f", "g", "s"),
      sep = ";",
      fill = "right",
      remove = FALSE
    ) %>%
    mutate(
      across(
        c(k, p, c, o, f, g),
        ~ sub("^[a-z]__", "", .x)
      )
    ) %>% 
    rename(
      kingdom = k, 
      phyllum = p, 
      class = c, 
      order = o, 
      family = f, 
      genus = g
    ) %>% 
    select(-Taxonomy) %>% 
    filter(abundance >0) %>% 
    mutate( marker = markers[i]) %>% 
    # At "pivot_longering" a matrix where rows are ASV abundances,
    # we now find a lot of ASV where abundance == 0 for those
    # ASV that were present at some samples, but not in others. 
    
    merge(list_merge[[i]]) 
  
  print("done")
}

data_micro_abundance <- do.call(rbind, list_result)


prokariota_nASV <- data_micro_abundance %>% 
  filter(marker == "16S") %>% 
  group_by(sampling, plot, treatment, date) %>% 
  summarize(nASV = n(),
            abundance = sum(abundance),
            .groups = "drop")

length(prokariota_nASV$abundance) == sum(prokariota_nASV$abundance)



fungi_nASV <- data_micro_abundance %>% 
  filter(marker == "ITS") %>% 
  group_by(sampling, plot, treatment, date) %>% 
  summarize(nASV = n(),
            abundance = sum(abundance),
            .groups = "drop") 

length(fungi_nASV$abundance) == sum(fungi_nASV$abundance)



##### 2. Opening dominance data #######
# Actually this data contains the number of ASV calculated in the previous data (observed features)


dominance_prokariota_1 <- read_tsv("data/16S_1stbatch_dominance.txt")
dominance_prokariota_2 <- read_tsv("data/16S_2ndbatch_dominance.txt")

prokariota_dominance <- bind_rows(dominance_prokariota_1, dominance_prokariota_2) %>%
  as.data.frame() %>% 
  rename(sample = Sample_Name) %>% 
  left_join(info_prokariota) %>% 
  mutate(marker = paste0("16S")) 



dominance_fungi_1 <- read_tsv("data/ITS_1stbatch_dominance.txt")
dominance_fungi_2 <- read_tsv("data/ITS_2ndbatch_dominance.txt")

fungi_dominance <- bind_rows(dominance_fungi_1, dominance_fungi_2) %>%
  as.data.frame() %>% 
  rename(sample = Sample_Name) %>% 
  left_join(info_fungi) %>% 
  mutate(marker = paste0("ITS")) 

dominance_micro <- rbind(prokariota_dominance, fungi_dominance) %>% 
  select(-sample, -label_micro)

#dominance_micro %>%  write.csv("data/data_micro_dominance.csv")





##### 3. Opening and arranging data FUNCTIONAL TRAITS #####


prokariota_ft_raw1 <- read.delim(
  "data/16S_1st_batch_relative_funtraits/faprotax.sample.anno.relative.xls",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
) 


#abundancias <- prokariota_ft_raw1[, sapply(prokariota_ft_raw1, is.numeric)]
#sum(abundancias > 0)

prokariota_ft_raw2 <- read.delim(
  "data/16S_2nd_batch_relative_funtraits/faprotax.sample.anno.relative.xls",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
) 

# There is also the option of "guild" instead of mode, that provide information. But we do not want
# that right?

fungi_ft_raw1 <- read.delim(
  "data/ITS_1st_batch_relative_funtraits/funguild.sample.mode.relative.xls",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

fungi_ft_raw2 <- read.delim(
  "data/ITS_2nd_batch_relative_funtraits/funguild.sample.mode.relative.xls",
  sep = "\t",
  check.names = FALSE)



list_data <- list(prokariota_ft_raw1, prokariota_ft_raw2, fungi_ft_raw1, fungi_ft_raw2)
list_merge <- list(info_prokariota, info_prokariota, info_fungi, info_fungi)
markers <- c("16S", "16S", "ITS", "ITS")
list_result <- list()

for(i in seq_along(list_data)){
  
  list_result[[i]] <- list_data[[i]] %>% 
    pivot_longer(
      cols = -c(description, feature),   # Keep Taxonomy and X.OTU_num unchanged
      names_to = "sample",              # Name for the new "variable" column
      values_to = "abundance"               # Name for the new "value" column
    )  %>%
    filter(abundance >0) %>% 
    mutate( marker = markers[i]) %>% 

    
    merge(list_merge[[i]]) 
  
  print("done")
}


data_micro_funtraits <- do.call(rbind, list_result)

prokariota_funtraits0 <- data_micro_funtraits %>% 
  filter(marker == "16S")

fungi_funtraits0 <- data_micro_funtraits %>% 
  filter(marker == "ITS")


sort(unique(prokariota_funtraits0$feature))
sort(unique(fungi_funtraits0$feature))


prokariota_funtraits_checking <- prokariota_funtraits0 %>% 
  select(-description) %>% 
  pivot_wider(
    names_from = "feature", 
    values_from = "abundance"
  ) 
summary(prokariota_funtraits_checking)


prokariota_funtraits <- prokariota_funtraits_checking %>% 
  select(where(~sum(is.na(.)) <= 12)) # I delete all funtraits with more than 12 NAs (10% of data)

cat(
  ncol(prokariota_funtraits_checking) - ncol(prokariota_funtraits), "out of", 
  length(unique(prokariota_funtraits0$feature)), " prokariota functional traits
have been removed beacuse they cointained more than 12 NA's (~10% of the data)")

# Checking if the NAs present are concentrated in some sampling x treatment 
prokariota_funtraits %>% 
  pivot_longer(
    cols = c(-sample, -marker, -label_micro, -plot,
               -treatment, -date, -day, -month, -year, -sampling,
             -date_label, - date_label_noyear),
    names_to = "funtrait", 
    values_to = "value"
  ) %>% 
  filter(is.na(value)) %>% 
  group_by(treatment, sampling, funtrait) %>% 
  summarize(
    n_na = n()
  ) %>% 
  filter(n_na > 1) %>% 
  print()

#Only in sampling 6, treatment p there are 2 NA in the funtrait "chitinolysis"


# List of functional traits excluded 
sort(prokariota_funtraits_checking %>% 
       select(where(~sum(is.na(.)) >= 12)) %>% 
       colnames())

# List of functional traits included
vector_prokariota_funtraits <- 
  sort(prokariota_funtraits %>% 
       select(-plot ,-sample, -marker, -label_micro, -treatment, -date, -day,
              -month, -year, -sampling,-date_label, - date_label_noyear) %>% 
       colnames())
print(vector_prokariota_funtraits)



fungi_funtraits_checking <- fungi_funtraits0 %>% 
  select(-description) %>% 
  pivot_wider(
    names_from = "feature", 
    values_from = "abundance"
  ) 
summary(fungi_funtraits_checking)

fungi_funtraits <- fungi_funtraits_checking %>% 
  select(where(~sum(is.na(.)) <= 12))
summary(prokariota_funtraits)


cat(
  ncol(fungi_funtraits_checking) - ncol(fungi_funtraits), "out of", 
  length(unique(fungi_funtraits0$feature)), "fungi functional traits 
  have been removed beacuse they cointained more than 12 NA's (~10% of the data)")

# Checking if the NAs present are concentrated in some sampling x treatment 
fungi_funtraits %>% 
  pivot_longer(
    cols = c(-sample, -marker, -label_micro, -plot,
             -treatment, -date, -day, -month, -year, -sampling, -date_label, - date_label_noyear),
    names_to = "funtrait", 
    values_to = "value"
  ) %>% 
  filter(is.na(value)) %>% 
  group_by(treatment, sampling, funtrait) %>% 
  summarize(
    n_na = n()
  ) %>% 
  filter(n_na > 1) %>% 
  print()
# Only in sampling 3 for wp treatment there are 2 NA for Saprotoph-Symbiotroph


# List of fungi functional traits included
vector_fungi_funtraits <- 
  sort(fungi_funtraits %>% 
         select(-plot ,-sample, -marker, -label_micro, -treatment, -date, -day, -month,
                -year, -sampling, -date_label, - date_label_noyear) %>% 
         colnames())








############ MERGING TABLES ##########


prokariota <- prokariota_nASV %>% 
  left_join(prokariota_dominance) %>% 
  left_join(prokariota_funtraits) %>% 
  select(
         marker, date, date_label_noyear, year, sampling, plot, treatment,
         nASV, observed_features, chao1, dominance, goods_coverage,
         pielou_e, shannon, simpson, all_of(vector_prokariota_funtraits))
str(prokariota)
  
prokariota %>%  write.csv("data/data_prokariota.csv")


fungi <- fungi_nASV %>% 
  left_join(fungi_dominance) %>% 
  left_join(fungi_funtraits) %>% 
  select(
    marker, date, date, date_label_noyear, year, sampling, plot, treatment,
    nASV, observed_features, chao1, dominance, goods_coverage,
    pielou_e, shannon, simpson, all_of(vector_fungi_funtraits))
str(fungi)

fungi %>%  write.csv("data/data_fungi.csv")






