################################################################################
#
# 5. Additional Output
#
################################################################################

TCHNconsmatVA <- readRDS(file.path("01_Data","Output","Tables","TCHNconsmatVA.rds")) 

TCHNconsmatVA %>%
  ungroup() %>%
  group_by(cap_wy) %>%
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
  theme_classic()

TCHNconsmatVA %>%
  ungroup() %>%
  group_by(age) %>%
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
  labs(x = "Estimated Age",
       y = "Saved Chinook Biomass (g)")+
  theme_classic()


