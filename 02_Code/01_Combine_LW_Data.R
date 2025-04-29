################## Retrieve Removed Fish Length-Weights ####################
library(tidyverse)
library(readxl)

# Grab the PRES data and reformat
pres_2016 <- read_xlsx(file.path("01_Data","Input","PRES_catch_data_all_years.xlsx"),
          sheet = "2016")

names(pres_2016) <- c("date","species","fork_length_mm","weight_lbs","mort","comments")

pres_2016 <- pres_2016 %>%
  filter(species == "STB") %>%
  mutate(weight_kg = round(weight_lbs*0.45359237,2),
         mort = ifelse(mort == "Y",1,0)) %>%
  select(date, species, fork_length_mm,weight_kg)

pres_2017 <-  read_xlsx(file.path("01_Data","Input","PRES_catch_data_all_years.xlsx"),
                        sheet = "2017")

names(pres_2017) <- c("date","month", "common_name","species_grp",
                      "fork_length_mm","total_length_mm","weight_lbs")
  
pres_2017 <- pres_2017 %>%
  filter(common_name == "Striped Bass") %>%
  mutate(weight_kg = round(weight_lbs*0.45359237,2),
        species = "STB") %>%
  select(date, species, fork_length_mm, weight_kg)


pres_2018 <-  read_xlsx(file.path("01_Data","Input","PRES_catch_data_all_years.xlsx"),
                        sheet = "2018")

names(pres_2018) <- c("sample_id", "catch_id", "date","month", "boat_id",
                      "location", "position", "common_name","species",
                      "fork_length_mm","total_length_mm","mort","weight_lbs",
                      "comments")

pres_2018 <- pres_2018 %>%
  filter(species == "STB") %>%
  mutate(weight_kg = round(weight_lbs*0.45359237,2)) %>%
  select(date, species, fork_length_mm, weight_kg)

pres <- bind_rows(pres_2016, pres_2017, pres_2018) %>%
  mutate(survey = "PRES")

# Grab the PFRS data and reformat
pfrs <-  read_xlsx(file.path("01_Data","Input","pfrs_processing_data.xlsx"),
                        sheet = 1)

names(pfrs) <- c("date", "crew", "time","gear", "species",
                 "fork_length_mm","weight_kg", "mort", "sample_id",
                 "transport_id","start_trans", "end_trans","comments")

pfrs <- pfrs %>%
  filter(species == "SB") %>%
  mutate(species = "STB",
         survey = "PFRS") %>%
  select(survey, date, species, fork_length_mm, weight_kg)

# Grab the EPFRRS data and reformat
epfrrs <- read_csv(file.path("01_Data","Input","EPFRRS.csv"))

names(epfrrs) <- c("date", "species", "fork_length_mm")

epfrrs <- epfrrs %>%
  filter(species == "striped-bass") %>%
  mutate(species = "STB",
         survey = "EPFRRS",
         date = mdy(date))

# Combine all

catch_comb <- bind_rows(pres, pfrs, epfrrs)

ggplot(data = catch_comb, aes(x = fork_length_mm, y = weight_kg))+
  geom_point(aes(color = survey),
             shape = 21)
