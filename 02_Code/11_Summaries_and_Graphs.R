################################################################################
#
# 5. Additional Output
#
################################################################################
library(tidyverse)

TCHNconsmatVA <- readRDS(file.path("01_Data","Output","Tables","TCHNconsmatVA.rds")) 
effort <- readRDS(file.path("01_Data","Input","effort.rds"))

TCHNconsmatVA <- TCHNconsmatVA %>% 
  ungroup() %>%
  filter(outlier == FALSE) %>%
  group_by(cap_wy, survey, gear, exit_timing, cap_date) %>%
  summarise(total_plt1_mean = sum(plt1_mean),
            total_plt1_lwr = sum(plt1_lwr),
            total_plt1_upr = sum(plt1_upr),
            total_p1_mean = sum(p1_mean),
            total_p1_lwr = sum(p1_lwr),
            total_p1_upr = sum(p1_upr)) %>%
  left_join(effort, by = c("cap_date" = "date", "gear"))%>%
  left_join(effort %>%
              filter(gear %in% c("Transport","Processing")) %>%
              rename("method" = gear,
                     "activity_hours" = labor_hours) %>%
              select(-weight), 
                     by = c("cap_date" = "date")) %>%
  pivot_wider(names_from = "method", values_from = "activity_hours") %>%
  mutate(effort_plt1_mean = total_plt1_mean/(labor_hours + (Processing*weight) + (Transport*weight)),
         effort_plt1_lwr = total_plt1_lwr/(labor_hours + (Processing*weight) + (Transport*weight)),
         effort_plt1_upr = total_plt1_upr/(labor_hours  + (Processing*weight) + (Transport*weight)),
         effort_p1_mean  = total_p1_mean/(labor_hours  + (Processing*weight) + (Transport*weight)),
         effort_p1_lwr = total_p1_lwr/(labor_hours  + (Processing*weight) + (Transport*weight)),
         effort_p1_upr = total_p1_upr/(labor_hours  + (Processing*weight) + (Transport*weight)))

TCHNconsmatVA %>%
  ungroup() %>%
  filter(exit_timing == "median") %>%
  group_by(cap_wy) %>%
  summarise(effort_plt1_mean = sum(effort_plt1_mean),
            effort_plt1_lwr = sum(effort_plt1_lwr),
            effort_plt1_upr = sum(effort_plt1_upr),
            effort_p1_mean = sum(effort_p1_mean),
            effort_p1_lwr = sum(effort_p1_lwr),
            effort_p1_upr = sum(effort_p1_upr)) %>%
  ggplot()+
  geom_ribbon(aes(x = cap_wy, ymin = effort_plt1_lwr, ymax = effort_plt1_upr),
              fill = "blue", alpha = 0.1)+
  geom_line(aes(x = cap_wy, y = effort_plt1_mean),
            color = "blue")+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     labels = scales::label_comma())+
  labs(x = "Capture Year",
       y = "Saved Chinook Biomass (g) / Hour of Labor")+
  theme_classic()

TCHNconsmatVA %>%
  ungroup() %>%
  filter(exit_timing == "median") %>%
  group_by(gear) %>%
  summarise(effort_plt1_mean = sum(effort_plt1_mean),
            effort_plt1_lwr = sum(effort_plt1_lwr),
            effort_plt1_upr = sum(effort_plt1_upr),
            effort_p1_mean = sum(effort_p1_mean),
            effort_p1_lwr = sum(effort_p1_lwr),
            effort_p1_upr = sum(effort_p1_upr)) %>%
  ggplot()+
  geom_pointrange(aes(x = gear, y = effort_plt1_mean,
                      ymin = effort_plt1_lwr, ymax = effort_plt1_upr))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     labels = scales::label_comma())+
  labs(x = "Gear Type",
       y = "Saved Chinook Biomass (g) / Hour of Labor")+
  theme_classic()
