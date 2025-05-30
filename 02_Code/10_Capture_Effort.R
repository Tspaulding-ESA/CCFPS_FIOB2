#####################################
#### Determine Daily Labor Hours ####
#####################################
library(tidyverse)


files <- list.files(file.path("01_Data","Input","CCF_Effort_Data"), 
                    full.names = TRUE)

effort_list <- list()
for(i in 1:length(files)){
  effort_list[[i]] <- read.csv(files[i])
}

effort <- bind_rows(effort_list)

effort_2 <- effort %>%
  mutate(
    activity_hours = ifelse(is.na(activity_hours),effort_hrs,activity_hours),
    date = lubridate::parse_date_time(date,
      orders = c("%Y%m%d","%m%d%Y")),
    survey = case_when(
      year(date) %in% c(2016:2018) ~ "PRES",
      year(date) %in% c(2019:2020) ~ "PFRS",
      year(date) %in% c(2021:2023) ~ "EPFRRS"
    ),
    gear = case_when(
      method == "electrofishing" ~ "E-fishing",
      method == "kodiak" ~ "Kodiak",
      method == "hoop" ~ "Hoop",
      method == "fyke" ~ "Fyke",
      method == "hook-and-line" ~ "Hook-Line",
      method %in% c("seine","Beach Seine") ~ "Seine",
      method == "lampara" ~ "Lampara",
      method == "processing" ~ "Processing",
      method == "transport" ~ "Transport"
    )
    ) %>%
  rowwise() %>%
  mutate(crew = str_replace_all(crew,", ,",",")) %>%
  mutate(crew = str_replace_all(crew,", ",",")) %>%
  mutate(crew = strsplit(crew, ",")) %>%
  unnest(crew) %>%
  mutate(crew = str_replace(crew, "\\s*","")) %>%
  mutate(crew = str_replace(crew, "\\t*","")) %>%
  mutate(crew = trimws(crew, which = "both", whitespace = "\\s*")) %>%
  filter(!is.na(crew)|survey %in% c("PRES","EPFRRS")|!is.null(crew)|crew != "") %>%
  group_by(date, survey, boat, gear, crew_count, start_time, end_time, hrs, activity_hours,
           labor_hours, effort_hrs)%>%
  summarise(crew = list(unique(crew))) %>%
  ungroup() %>%
  rowwise() %>%
    mutate(crew_length = length(crew)) %>%
  mutate(crew_count = ifelse(is.na(crew_count),crew_length, crew_count))


# Build LU for average activity hours
activity_hours <- effort_2 %>%
  group_by(gear) %>%
  summarise(avg_act_hrs = mean(activity_hours, na.rm = TRUE)) %>%
  mutate(avg_act_hrs = ifelse(is.na(avg_act_hrs),
                              mean(avg_act_hrs, na.rm = TRUE),
                              avg_act_hrs))

activity_effort <- effort_2 %>%
  left_join(activity_hours) %>%
  mutate(activity_hours = ifelse(is.na(activity_hours),avg_act_hrs,activity_hours)) %>%
  mutate(labor_hours = activity_hours * crew_count) %>%
  group_by(date, gear) %>%
  summarize(labor_hours = sum(labor_hours, na.rm = TRUE))

saveRDS(activity_effort, file.path("01_Data","Input","effort.rds"))         
