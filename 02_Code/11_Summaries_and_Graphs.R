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
  group_by(outlier, cap_wy, survey, gear, exit_timing, cap_date) %>%
  summarise(total_plt1_mean = sum(plt1_mean),
            total_plt1_lwr = sum(plt1_lwr),
            total_plt1_upr = sum(plt1_upr),
            total_p1_mean = sum(p1_mean),
            total_p1_lwr = sum(p1_lwr),
            total_p1_upr = sum(p1_upr)) %>%
  left_join(effort, by = c("cap_date" = "date", "gear"))%>%
  mutate(effort_plt1_mean = total_plt1_mean/(gear_hours + proc_hours + trans_hours),
         effort_plt1_lwr = total_plt1_lwr/(gear_hours + proc_hours + trans_hours),
         effort_plt1_upr = total_plt1_upr/(gear_hours  + proc_hours + trans_hours),
         effort_p1_mean  = total_p1_mean/(gear_hours  + proc_hours + trans_hours),
         effort_p1_lwr = total_p1_lwr/(gear_hours  + proc_hours + trans_hours),
         effort_p1_upr = total_p1_upr/(gear_hours  + proc_hours + trans_hours))

TCHNconsmatVA %>%
  ungroup() %>%
  filter(exit_timing == "median") %>%
  group_by(outlier, cap_wy) %>%
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
  facet_grid(rows = vars(outlier)#, scales = "free_y"
             )+
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

TCHNconsmatVA %>%
  ungroup() %>%
  filter(exit_timing == "median") %>%
  group_by(month = month(cap_date),
           cap_wy) %>%
  summarise(effort_plt1_mean = sum(effort_plt1_mean),
            effort_plt1_lwr = sum(effort_plt1_lwr),
            effort_plt1_upr = sum(effort_plt1_upr),
            effort_p1_mean = sum(effort_p1_mean),
            effort_p1_lwr = sum(effort_p1_lwr),
            effort_p1_upr = sum(effort_p1_upr)) %>%
  mutate(month = ifelse(month > 9, month - 9, month + 3)) %>%
  ggplot()+
  geom_ribbon(aes(x = month, ymin = effort_plt1_lwr, ymax = effort_plt1_upr),
              fill = "blue", alpha = 0.1)+
  geom_line(aes(x = month, y = effort_plt1_mean),
            color = "blue")+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5),
                     labels = scales::label_comma())+
  scale_x_continuous(breaks = c(1:9),
                     labels = month.abb[c(10:12,1:6)])+
  facet_grid(rows = vars(cap_wy), scales = "free_y")+
  labs(x = "Capture Month",
       y = "Saved Chinook Biomass (g) / Hour of Labor")+
  theme_classic()
