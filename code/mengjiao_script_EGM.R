





library(googledrive)
library("nlme")
library("car")
library("tidyverse")
library("ggpubr")

##################### download data from google drive ###################
# create local folder
out_dir <- "data/EGM"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- shared_drive_find(pattern = "RECODYN") %>% drive_ls %>% # find recodyn drive
  filter(name == "Data") %>% drive_ls %>% # go into data
  filter(name == "EGM-5") %>% drive_ls %>% # go into EGM-5 folder
  filter(name == "raw_data") %>% drive_ls %>% 
  filter(grepl("^2.*\\.TXT$", name)) %>%   # get file start with 2
  arrange(name)

print(files)
file_names <- files$name

# download files
if (nrow(files) > 0) {
  for (i in seq_len(nrow(files))) {
    drive_download(
      file = files$id[[i]],
      path = file.path(out_dir, file_names[i]),
      overwrite = TRUE
    )
  }
} else {
  message("No matching files found.")
}

######################## data clean #####################
# combine data to one file
data_path <- "data/EGM"

combined_data <- list.files(data_path, pattern = "\\.TXT$", full.names = TRUE) %>%
  lapply(function(file) {
    data <- read.csv(file, header = TRUE)
    data$SourceFile <- file 
    return(data)}) %>%
  bind_rows %>%
  as_tibble

# correct date

fixes <- list(
  march = list(
    orig = combined_data %>%
      filter(Date == "28/03/24", Time != "", !is.na(Time)) %>%
      summarise(dt = dmy_hms(paste("28/03/2024", first(Time)))) %>%
      pull(dt),
    corr = dmy_hms("07/04/2024 10:15:00")
  ),
  
  june_2024 = list(
    orig = combined_data %>%
      filter(Date == "20/06/24", Time != "", !is.na(Time)) %>%
      summarise(dt = dmy_hms(paste("20/06/2024", first(Time)))) %>%
      pull(dt),
    corr = dmy_hms("26/06/2024 09:30:00")
  ),
  
  june_2025 = list(
    orig = combined_data %>%
      filter(Date == "29/06/25", Time != "", !is.na(Time)) %>%
      summarise(dt = dmy_hms(paste("29/06/2025", first(Time)))) %>%
      pull(dt),
    corr = dmy_hms("17/07/2025 09:30:00")
  ),
  
  july_2025 = list(
    orig = combined_data %>%
      filter(Date == "29/07/25", Time != "", !is.na(Time)) %>%
      summarise(dt = dmy_hms(paste("29/07/2025", first(Time)))) %>%
      pull(dt),
    corr = dmy_hms("17/07/2025 09:30:00")   
  ),
  
  august = list(
    orig = combined_data %>%
      filter(Date == "02/08/24", Time != "", !is.na(Time)) %>%
      summarise(dt = dmy_hms(paste("02/08/2024", first(Time)))) %>%
      pull(dt),
    corr = dmy_hms("08/08/2024 09:30:00")
  )
)

fixes <- lapply(fixes, function(i) {
  i$difference <- i$corr - i$orig
  i
})



fixed_data <- combined_data %>%
  mutate(
    date_time = dmy_hms(paste(Date, Time), quiet = TRUE),
    date_only = dmy(Date, quiet = TRUE),
    
    shift_sec = case_when(
      date_only <  dmy("30/03/2024")  ~ as.numeric(fixes$march$difference, units = "secs"),
      date_only == dmy("20/06/2024")  ~ as.numeric(fixes$june_2024$difference, units = "secs"),
      date_only == dmy("02/08/2024")  ~ as.numeric(fixes$august$difference, units = "secs"),
      date_only == dmy("29/06/2025")  ~ as.numeric(fixes$june_2025$difference, units = "secs"),
      date_only == dmy("29/07/2025")  ~ as.numeric(fixes$july_2025$difference, units = "secs"),
      TRUE ~ 0
    ),
    
    fixed_date_time = date_time + dseconds(shift_sec),
    
    Date_original = Date,
    Time_original = Time,
    
    Date = as.Date(fixed_date_time),
    Time = hour(fixed_date_time) +
      minute(fixed_date_time) / 60 +
      second(fixed_date_time) / 3600
  ) %>%
  mutate(
    sampling_number = case_when(
      Date == as.Date("2024-04-07") ~ "202402",
      Date == as.Date("2024-04-26") ~ "202403",
      Date == as.Date("2024-06-26") ~ "202406",
      Date == as.Date("2024-08-08") ~ "202408",
      Date == as.Date("2024-09-18") ~ "202410",
      Date == as.Date("2025-04-02") ~ "202502", 
      Date == as.Date("2025-04-25") ~ "202503",
      Date == as.Date("2025-07-17") ~ "202507",
      Date == as.Date("2025-08-09") ~ "202508", 
      Date == as.Date("2025-09-19") ~ "202510",
      Date == as.Date("2025-09-20") ~ "202510",
      TRUE ~ NA_character_
    )
  )



# select all the rows between start and end 
selected_data <- fixed_data %>%
  mutate(
    tag_clean = trimws(Tag.M5.),
    Interval = cumsum(tag_clean == "Start")
  ) %>%
  group_by(Interval) %>%
  mutate(
    row_in_interval = row_number(),
    end_row = match("End", tag_clean)
  ) %>%
  ungroup() %>%
  # some start doesn't have end row
  filter(
    Interval > 0,              # remove mess rows before start
    !is.na(end_row),           # keep rows with end
    row_in_interval <= end_row 
  ) %>%
  # fixes
  mutate(
    Plot_No = ifelse(Interval == 425 & Plot_No == 182, 184, Plot_No),
    Plot_No = ifelse(Interval == 424 & Plot_No == 181, 182, Plot_No),
    Plot_No = ifelse(Interval == 885 & Plot_No == 222, 212, Plot_No)
  ) %>%
  # remove plots with Plot_No < 10
  filter(Plot_No > 10) %>%
  # separate Plot_No into cage and subplot
  mutate(
    subplot = as.factor(Plot_No %% 10),
    cage = Plot_No %/% 10
  ) %>%
  # remove Date with NA: "start" rows
  drop_na(sampling_number)

# make sure each plot has 4 different measurement
check_subplot <- selected_data %>%
  mutate(
    subplot = as.factor(Plot_No %% 10),                
    cage = Plot_No %/% 10
  ) %>%
  drop_na(Plot_No) %>%
  distinct(sampling_number, cage, subplot) %>%
  group_by(sampling_number, cage) %>%
  summarise(
    n_subplot = n(),
    subplot_set = paste(sort(subplot), collapse = ","),
    .groups = "drop"
  ) %>%
  filter(n_subplot != 4 | subplot_set != "1,2,3,4")

check_subplot

subplot_summary <- selected_data %>%
  distinct(sampling_number, cage, subplot) %>%
  group_by(sampling_number, cage) %>%
  summarise(
    subplot_set = paste(sort(subplot), collapse = ","),
    n_subplot = n(),
    .groups = "drop"
  ) %>%
  filter(n_subplot>4)

subplot_summary

cage_summary <- selected_data %>%
  distinct(sampling_number, cage) %>%
  group_by(sampling_number) %>%
  summarise(
    cage_set = paste(sort(cage), collapse = ","),
    n_cage = n(),
    .groups = "drop"
  )

cage_summary

# remove all the replicated data in each plot
data_unique <- selected_data %>%
  group_by(sampling_number, Plot_No) %>%
  slice_max(Interval) %>%
  ungroup()

unique(data_unique$end_row) # all 83 rows

# Tsoil and Msoil are 0 at some time-points
# remove these 0 rows and then take the average values
mean_soil <- data_unique %>%
  group_by(Interval) %>%
  mutate(
    Tsoil_avg = ifelse(Tag.M5. == "R5", mean(Tsoil[Tsoil != 0], na.rm = TRUE), NA),
    Tair_avg = ifelse(Tag.M5. == "R5", mean(Tair[Tair != 0], na.rm = TRUE), NA),
    Msoil_avg = ifelse(Tag.M5. == "R5", mean(Msoil[Msoil != 0], na.rm = TRUE), NA)
  ) %>%
  ungroup

# renew data in R5 row with the averaged data
R5_new <- mean_soil %>%
  mutate(
    Tsoil = ifelse(Tag.M5. == "R5", Tsoil_avg, Tsoil),
    Tair = ifelse(Tag.M5. == "R5", Tair_avg, Tair),
    Msoil = ifelse(Tag.M5. == "R5", Msoil_avg, Msoil)
  ) %>%
  dplyr::select(-c(Tsoil_avg, Tair_avg, Msoil_avg)) 

# select all rows whose TagM5 is "R5"
select_R5 <- R5_new %>%
  filter(Tag.M5. == "R5") %>%
  # mutate(Time = as.numeric(hms::as_hms(Time)) / 3600) %>%
  dplyr::select(Date, Time, fixed_date_time, sampling_number, Plot_No, cage, subplot,
                Tsoil, Tair, Msoil, P4, Interval) %>%
  ungroup %>%
  rename(SRL = P4, DateTime = fixed_date_time)

check_subplot <- select_R5 %>%
  distinct(sampling_number, cage, subplot) %>%
  group_by(sampling_number, cage) %>%
  summarise(
    n_subplot = n(),
    subplot_set = paste(sort(subplot), collapse = ","),
    .groups = "drop"
  ) %>%
  filter(n_subplot != 4 | subplot_set != "1,2,3,4")

check_subplot

# ggplot(select_R5) + 
#   aes(x = Date, y = Time, color = sampling_number) +
#   geom_point() +
#   theme_bw() 

# mean value of each cage and collar
collar_R5 <- select_R5 %>%
  group_by(sampling_number, cage, subplot) %>%
  summarise(
    Tsoil = mean(Tsoil, na.rm = TRUE),
    Tair = mean(Tair, na.rm = TRUE),
    Msoil = mean(Msoil, na.rm = TRUE),
    SRL = mean(SRL, na.rm = TRUE))

# mean value of each cage
mean_R5 <- select_R5 %>%
  group_by(sampling_number, cage) %>%
  summarise(
    Tsoil = mean(Tsoil, na.rm = TRUE),
    Tair = mean(Tair, na.rm = TRUE),
    Msoil = mean(Msoil, na.rm = TRUE),
    SRL = mean(SRL, na.rm = TRUE))



ggsave("Figures/CWM herbivory_sampling number.png", width = 6, height = 4.5)

saveRDS(mean_R5, "data/mean_R5.RData")
mean_R5 <- readRDS("data/mean_R5.RData")