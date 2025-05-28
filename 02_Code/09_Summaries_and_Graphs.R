################################################################################
#
# 5. Additional Output
#
################################################################################
library(tidyverse)

TCHNconsmatVA <- readRDS(file.path("01_Data","Output","Tables","TCHNconsmatVA.rds")) 

TCHNconsmatVA %>%
  ungroup() %>%
  group_by(cap_wy, outlier) %>%
  summarise(total_plt1_mean = sum(plt1_mean),
            total_plt1_lwr = sum(plt1_lwr),
            total_plt1_upr = sum(plt1_upr),
            total_p1_mean = sum(p1_mean),
            total_p1_lwr = sum(p1_lwr),
            total_p1_upr = sum(p1_upr)) %>%
  ggplot()+
  geom_ribbon(aes(x = cap_wy, ymin = total_plt1_lwr, ymax = total_plt1_upr),
              fill = "blue", alpha = 0.1)+
  geom_line(aes(x = cap_wy, y = total_plt1_mean),
            color = "blue")+
  facet_grid(rows = vars(outlier), scales = "free_y")+
  labs(x = "Capture Year",
       y = "Saved Chinook Biomass (g)")+
  theme_classic()

TCHNconsmatVA %>%
  ungroup() %>%
  group_by(age, outlier) %>%
  summarise(total_plt1_mean = sum(plt1_mean),
            total_plt1_lwr = sum(plt1_lwr),
            total_plt1_upr = sum(plt1_upr),
            total_p1_mean = sum(p1_mean),
            total_p1_lwr = sum(p1_lwr),
            total_p1_upr = sum(p1_upr)) %>%
  ggplot()+
  geom_ribbon(aes(x = age, ymin = total_plt1_lwr, ymax = total_plt1_upr),
              fill = "blue", alpha = 0.1)+
  geom_line(aes(x = age, y = total_plt1_mean),
            color = "blue")+
  facet_grid(rows = vars(outlier), scales = "free_y")+
  labs(x = "Estimated Age",
       y = "Saved Chinook Biomass (g)")+
  theme_classic()

TCHNconsmatVA %>%
  ungroup() %>%
  group_by(gear) %>%
  summarise(total_plt1_mean = sum(plt1_mean),
            total_plt1_lwr = sum(plt1_lwr),
            total_plt1_upr = sum(plt1_upr),
            total_p1_mean = sum(p1_mean),
            total_p1_lwr = sum(p1_lwr),
            total_p1_upr = sum(p1_upr),
            n = n()) %>%
  ggplot()+
  geom_pointrange(aes(x = gear, y= total_plt1_mean, ymin = total_plt1_lwr, ymax = total_plt1_upr),
              )+
  geom_point(aes(x = gear, y = total_plt1_mean),
            color = "blue")+
  geom_text(aes( x = gear,y = -100000, label = paste("n =",n)))+
  labs(x = "Sample Gear",
       y = "Saved Chinook Biomass (g)")+
  theme_classic()

TCHNconsmatVA %>%
  ungroup() %>%
  group_by(month = month(cap_date)) %>%
  ggplot()+
  geom_histogram(aes(x = fork_length_mm), binwidth = 10)+
  facet_grid(rows = vars(month))+
  labs(x = "Fork Length (mm)",
       y = "Saved Chinook Biomass (g)")+
  theme_classic()


