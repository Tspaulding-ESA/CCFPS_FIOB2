##Timing Variability

library(tidyverse)
library(survival)

surv_mod <- readRDS(file.path("01_Data","Input","age_emmigration_model.rds"))

df = data.frame(age = c(1:10),
                AgeBin = c("1-2","1-2","3-5","3-5","3-5","6+","6+","6+","6+","6+"))


probs <- data.frame(7*round(predict(surv_mod, newdata = df, 
                            type = "quantile", p = c(0.25,0.5,0.75)),0))
names(probs) <- c("early","median","late")

em_probs <- bind_cols(df,probs) %>%
  pivot_longer(cols = early:late, names_to = "exit_timing", 
               values_to = "exit_days")

saveRDS(em_probs,file.path("01_Data","Input","emigration_probs.rds"))

# Age ID
age_label <- rep(c("1-2","1-2","3-5","3-5","3-5","6+","6+","6+","6+","6+"), 1000)

# Predict linear predictor
lp <- rep(predict(surv_mod, newdata = df, type = "lp"), times = 1000)

# Extract scale parameter (sigma)
sigma <- rep(surv_mod$scale,10000)

# Create a sequence of time values to evaluate survival function
time_seq <- rep(seq(1, 1000, by = 1),times = 10)  # adjust upper limit as needed

# Survival function for Weibull model: 
# S(t) = 1 - F(t), and F(t) = pnorm((log(t) - lp) / sigma) if Gaussian
# But for Weibull in survreg(), use:
# S(t) = exp(-exp((log(t) - lp) / sigma))

surv_probs <- 1-(exp(-exp((log(time_seq) - lp) / sigma)))

# Construct a data frame for plotting
plot_data <- data.frame(age = age_label, time = time_seq, survival = surv_probs)


ggplot(plot_data, aes(x = time, y = survival, color = age)) +
  geom_line(size = 1.2) +
  labs(
    title = "Predicted Exit Curve (Weibull)",
    x = "Time (days post-capture)",
    y = "Probability of Exit"
  ) +
  scale_y_continuous(labels = scales::percent)+
  facet_grid(rows = vars(age))+
  theme_classic()+
  theme(
    panel.grid.major.y = element_line(color = "grey80")
  )

ggsave(file.path("01_Data","Output","Figures","Predicted_Exit_Curve.png"), device = "png",
       width = 6, height = 4, units = "in", dpi = 300)
