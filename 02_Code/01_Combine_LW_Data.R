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
  # mutate(weight_lbs = ifelse(fork_length_mm == 89 & weight_lbs == 0.2, 0.02, weight_lbs),
  #        weight_lbs = ifelse(fork_length_mm == 388 & weight_lbs == 0.15, 1.5, weight_lbs)) %>%
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
         # weight_kg = ifelse(fork_length_mm == 122 & weight_kg == 0.35, 0.035, weight_kg),
         # weight_kg = ifelse(fork_length_mm == 525 & weight_kg == 0.26, 2.6, weight_kg),
         # weight_kg = ifelse(fork_length_mm == 225 & weight_kg == 0.9, 0.09, weight_kg),
         # weight_kg = ifelse(fork_length_mm == 215 & weight_kg == 0.75, 0.075, weight_kg),
         # fork_length_mm = ifelse(fork_length_mm == 34 & weight_kg == 0.48, 340, fork_length_mm),
         # fork_length_mm = ifelse(fork_length_mm == 305 & weight_kg == 0.06, 205, fork_length_mm),
         # fork_length_mm = ifelse(fork_length_mm == 368 & weight_kg == 0.12, 268, fork_length_mm),
         # fork_length_mm = ifelse(fork_length_mm == 149 & weight_kg == 0.18, 249, fork_length_mm),
         # fork_length_mm = ifelse(fork_length_mm == 138 & weight_kg == 0.13, 238, fork_length_mm),
         # fork_length_mm = ifelse(fork_length_mm == 255 & weight_kg == 0.05, 155, fork_length_mm)
         ) %>%
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
# different relationship fits best for short and long fish
fl_thresh = 400

rr_mod_short = MASS::rlm(log(weight_kg) ~ log(fork_length_mm),
                         data = filter(catch_comb_lw, fork_length_mm < fl_thresh))

rr_mod_long = MASS::rlm(log(weight_kg) ~ log(fork_length_mm),
                        data = filter(catch_comb_lw, fork_length_mm >= fl_thresh))

catch_comb_lw_resid <- catch_comb_lw %>%
  filter(fork_length_mm < fl_thresh) %>%
  mutate(residual = rr_mod_short$w) %>%
  bind_rows(filter(catch_comb_lw, fork_length_mm >= fl_thresh) %>%
              mutate(residual = rr_mod_long$w)) %>%
  mutate(outlier = residual < 0.25)

# 0.25 weight cutoff comes from Kimmerer et al. 2005
catch_comb_lw_resid = catch_comb_lw_resid %>%
  mutate(weight_kg_pred = case_when(
    fork_length_mm < fl_thresh ~ pred_weight(rr_mod_short, fork_length_mm),
    fork_length_mm >= fl_thresh ~ pred_weight(rr_mod_long, fork_length_mm)
  )) %>%
  mutate(diff = weight_kg_pred - weight_kg,
         prop_diff = diff/weight_kg) %>%
  arrange(outlier)

max_outlier <- max(catch_comb_lw_resid$fork_length_mm[catch_comb_lw_resid$outlier == TRUE])

(full_plot <- ggplot(data = catch_comb_lw_resid, aes(x = fork_length_mm)) +
  geom_point(aes(y = weight_kg, color = outlier, size = outlier), alpha = 0.5) +
  geom_line(aes(y = weight_kg_pred), color = "blue")+
  scale_size_manual(breaks = c(FALSE,TRUE),
                    values = c(1.5,2))+
  scale_color_manual(breaks = c(FALSE,TRUE),
                       values = c("grey50","tomato"))+
  scale_x_continuous(breaks = seq(0,1300,100))+
  scale_y_continuous(breaks = seq(0,25,2))+
  theme_classic()+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = -45, hjust = 0)
  )+
    labs(x = "Fork Length (mm)", y = "Weight (kg)"))

(inset_plot <- ggplot(data = catch_comb_lw_resid %>% filter(fork_length_mm <= max_outlier), 
                     aes(x = fork_length_mm)) +
  geom_point(aes(y = weight_kg, color = outlier), size = 2, alpha = 0.5) +
  geom_line(aes(y = weight_kg_pred), color = "blue")+
  scale_color_manual(breaks = c(FALSE,TRUE),
                     values = c("grey50","tomato"))+
  scale_x_continuous(breaks = seq(0,1300,100))+
  scale_y_continuous(breaks = seq(0,3,1))+
  theme_classic()+
  theme(
    legend.position = "none"
  )+
    labs(x = "Fork Length (mm)", y = "Weight (kg)"))

(lw_plot <- full_plot+
  annotation_custom(
    ggplotGrob(inset_plot),
    xmin = 0, xmax = 750,
    ymin = 9, ymax = 25
  )+
  annotate("rect", xmin = 25, ymin = -1, xmax = max_outlier, ymax = 4,
           color = "grey30", fill = NA, linewidth = 0.75, alpha = 0.5)+
  annotate("rect", xmin = 0, ymin = 9, xmax = 775, ymax = 25,
           color = "grey30", fill = NA, linewidth = 0.75, alpha = 0.5)+
  annotate("segment", x = 25, y = 4, xend = 0, yend = 9,
           color = "grey50", linewidth = 0.5, alpha = 0.5)+
  annotate("segment", x = max_outlier, y = 4, xend = 775, yend = 9,
           color = "grey50", linewidth = 0.5, alpha = 0.5))

ggsave(file.path("01_Data","Output","Figures","LW_outliers.png"), plot = lw_plot,
       height = 4, width = 4, units = "in", dpi = 300)

catch_comb <- catch_comb %>%
  left_join(select(catch_comb_lw_resid, date:weight_kg, outlier)) %>%
  replace_na(list("outlier" = FALSE))

write.csv(catch_comb, file.path("01_Data","Input","SBlw.csv"))