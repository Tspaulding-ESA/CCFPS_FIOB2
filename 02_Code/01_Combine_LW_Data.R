################## Retrieve Removed Fish Length-Weights ####################
library(tidyverse)
library(readxl)

data_path = file.path("01_Data", "Input")

pred_weight <- function(mod, fl){
  coefs = unname(coef(mod))
  exp(coefs[1] + log(fl)*coefs[2])
}

# Grab the PRES data and reformat
pres_2016 <- read_xlsx(file.path(data_path, "PRES_catch_data_all_years.xlsx"),
                       sheet = "2016")

names(pres_2016) <- c("date","species","fork_length_mm","weight_lbs","mort","comments")

pres_2016 <- pres_2016 %>%
  filter(species == "STB") %>%
  mutate(weight_lbs = ifelse(fork_length_mm == 89 & weight_lbs == 0.2, 0.02, weight_lbs),
         weight_lbs = ifelse(fork_length_mm == 388 & weight_lbs == 0.15, 1.5, weight_lbs)) %>%
  select(date, species, fork_length_mm, weight_lbs)

pres_2017 <-  read_xlsx(file.path(data_path, "PRES_catch_data_all_years.xlsx"),
                        sheet = "2017")

names(pres_2017) <- c("date","month", "common_name","species_grp",
                      "fork_length_mm","total_length_mm","weight_lbs")

pres_2017 <- pres_2017 %>%
  filter(common_name == "Striped Bass") %>%
  mutate(species = "STB") %>%
  select(date, species, fork_length_mm, weight_lbs)

pres_2018 <-  read_xlsx(file.path(data_path, "PRES_catch_data_all_years.xlsx"),
                        sheet = "2018")

names(pres_2018) <- c("sample_id", "catch_id", "date","month", "boat_id",
                      "location", "position", "common_name","species",
                      "fork_length_mm","total_length_mm","mort","weight_lbs",
                      "comments")

pres_2018 <- pres_2018 %>%
  filter(species == "STB") %>%
  select(date, species, fork_length_mm, weight_lbs)

pres <- bind_rows(pres_2016, pres_2017, pres_2018) %>%
  mutate(survey = "PRES",
         gear = "E-fishing",
         weight_kg = round(weight_lbs * 0.45359237, 2)) %>% 
  select(-weight_lbs)

# Grab the PFRS data and reformat
pfrs <- read_xlsx(file.path(data_path, "pfrs_processing_data.xlsx"),
                  sheet = 1)

names(pfrs) <- c("date", "crew", "time","gear", "species",
                 "fork_length_mm","weight_kg", "mort", "sample_id",
                 "transport_id","start_trans", "end_trans","comments")

pfrs <- pfrs %>%
  filter(species == "SB") %>%
  mutate(species = "STB",
         survey = "PFRS",
         gear = ifelse(gear == "lampara", "Lampara", gear),
         weight_kg = ifelse(fork_length_mm == 122 & weight_kg == 0.35, 0.035, weight_kg),
         weight_kg = ifelse(fork_length_mm == 525 & weight_kg == 0.26, 2.6, weight_kg),
         weight_kg = ifelse(fork_length_mm == 225 & weight_kg == 0.9, 0.09, weight_kg),
         weight_kg = ifelse(fork_length_mm == 215 & weight_kg == 0.75, 0.075, weight_kg),
         fork_length_mm = ifelse(fork_length_mm == 34 & weight_kg == 0.48, 340, fork_length_mm),
         fork_length_mm = ifelse(fork_length_mm == 305 & weight_kg == 0.06, 205, fork_length_mm),
         fork_length_mm = ifelse(fork_length_mm == 368 & weight_kg == 0.12, 268, fork_length_mm),
         fork_length_mm = ifelse(fork_length_mm == 149 & weight_kg == 0.18, 249, fork_length_mm),
         fork_length_mm = ifelse(fork_length_mm == 138 & weight_kg == 0.13, 238, fork_length_mm),
         fork_length_mm = ifelse(fork_length_mm == 255 & weight_kg == 0.05, 155, fork_length_mm)) %>%
  select(survey, date, gear, species, fork_length_mm, weight_kg)

# Grab the EPFRRS data and reformat
epfrrs_gear = lapply(list.files(data_path, "_Metadata", full.names = TRUE), 
                     \(x){select(read.csv(x), processing_record_id, gear = SAMPLE_METHOD)}) %>%
  bind_rows() %>%
  filter(!is.na(processing_record_id) & processing_record_id != "") %>%
  distinct()

epfrrs = read.csv(file.path(data_path, "FishMeasurements.csv")) %>%
  filter(COMMON_NAME == "striped-bass") %>%
  left_join(epfrrs_gear) %>% 
  mutate(species = "STB",
         survey = "EPFRRS",
         date = ymd(SAMPLING_DATETIME),
         gear = case_when(is.na(gear) & hl_record_id != "" ~ "Hook-Line",
                          gear == "boat-electrofishing" ~ "E-fishing",
                          gear == "seining" ~ "Seine",
                          gear == "hoop" ~ "Hoop",
                          gear == "kodiak-trawl" ~ "Kodiak",
                          .default = gear)) %>% 
  select(survey, date, gear, species, fork_length_mm = FORK_LENGTH_MM)

# Combine all

catch_comb <- bind_rows(pres, pfrs, epfrrs) %>% 
  filter(!(is.na(fork_length_mm) & is.na(weight_kg))) 

ggplot(data = catch_comb, aes(x = log(fork_length_mm), y = log(weight_kg))) +
  geom_point(aes(color = survey), shape = 21) +
  facet_wrap(~gear)

catch_comb_lw = filter(catch_comb, !is.na(fork_length_mm) & !is.na(weight_kg))

# robust regression to reduce influence of outliers
rr_mod = MASS::rlm(log(weight_kg) ~ log(fork_length_mm), data = catch_comb_lw)

catch_comb_lw = mutate(catch_comb_lw, weight_kg_pred = pred_weight(rr_mod, fork_length_mm))

ggplot(data = catch_comb_lw, aes(x = fork_length_mm)) +
  geom_point(aes(y = weight_kg), size = 2, alpha = 0.2) +
  geom_line(aes(y = weight_kg_pred))

# 0.25 weight cutoff comes from Kimmerer et al. 2005
catch_comb_lw_out = catch_comb_lw %>% 
  filter(rr_mod$w < 0.25) %>% 
  mutate(diff = weight_kg_pred - weight_kg,
         prop_diff = diff/weight_kg)

p = ggplot(data = catch_comb_lw_out, 
           aes(x = log(fork_length_mm), y = log(weight_kg), col = prop_diff,
               text = paste("fl:", fork_length_mm, 
                            "\nwt:", weight_kg,
                            "\nwt_pred:", weight_kg_pred))) +
  geom_point() +
  scale_color_continuous(type = "viridis")

plotly::ggplotly(p)

# compare lengths that were used in length-weight (which does not include outliers)
# to fork lengths with no weight measurement
catch_comb_comp = catch_comb %>% 
  mutate(in_lw = !is.na(fork_length_mm) & !is.na(weight_kg))

# nothing jumping out at me that seems like it needs closer inspection
ggplot(catch_comb_comp, aes(x = gear, y = fork_length_mm, fill = in_lw)) +
  geom_boxplot(alpha = 0.4)




# # below this point I was exploring the idea that the shorter fork lengths were driving 
# # the poor for longer fork lengths (based on log-log model)
# # not worth pursuing further but retaining here for FYI
# fl_thresh = 482
# 
# rr_mod_short = MASS::rlm(log(weight_kg) ~ log(fork_length_mm), 
#                          data = filter(catch_comb_lw, fork_length_mm < fl_thresh))
# 
# rr_mod_long = MASS::rlm(log(weight_kg) ~ log(fork_length_mm), 
#                         data = filter(catch_comb_lw, fork_length_mm >= fl_thresh))
# 
# catch_comb_lw = catch_comb_lw %>% 
#   mutate(weight_kg_pred_short = pred_weight(rr_mod_short, fork_length_mm),
#          weight_kg_pred_long = pred_weight(rr_mod_long, fork_length_mm))
# 
# # different relationship fits best for short and long fish
# # maybe obvious (statistical) reason for this?
# ggplot(data = catch_comb_lw, aes(x = fork_length_mm)) +
#   geom_point(aes(y = weight_kg), size = 2, alpha = 0.2) +
#   geom_line(aes(y = weight_kg_pred_short), col = "red") +
#   geom_line(aes(y = weight_kg_pred_long), col = "blue")
