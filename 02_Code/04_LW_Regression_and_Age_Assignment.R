library(tidyverse)

# predator l-w regression
LWdata <- read.csv(file.path("01_Data","Input","SBlw.csv"))

# fit the regression 	
LWdata <- LWdata %>%
  mutate(weight_g = weight_kg*1000)
LW_log <- lm(log(weight_g) ~ log(fork_length_mm), data = LWdata)			

# plot 	
xx <- min(LWdata$fork_length_mm, na.rm = T):max(LWdata$fork_length_mm, na.rm = T)

pred_log <- data.frame(predict(LW_log, newdata = data.frame("fork_length_mm" = xx), 
                 interval = "prediction")) %>%
  bind_cols(xx)
names(pred_log) <- c("fit","lwr","upr","fork_length_mm")

# predator sizes and age classes

# for the few removed predators that don't have lengths, we assume the 
#   mean size of fish removed on that day (exp-mean-log)

LWdataA <- LWdata
tofillin <- which(is.na(LWdata$fork_length_mm) == T)

for(i in 1:length(tofillin)){
  datefished <-LWdata$date[tofillin[i]]
  otherlengths <- na.omit(LWdata$fork_length_mm[
    which(LWdata$date == datefished)])
  assumedlength <- round(exp(mean(log(otherlengths))))
  LWdataA$fork_length_mm[tofillin[i]] <- assumedlength
}

# for any predators that don't have weights, we use the nls regession

# Fit all fish to new weights using NLS
LW_nls <- nls(weight_g ~ a * fork_length_mm^b, 
              data = LWdataA, 
              start = list(a = 3.28e-5,
                           b = 2.85))

# Then, bootstrap and simulate the data to create a prediction interval
pred_wt_bt <- nlraa::boot_nls(LW_nls, fitted) ## This takes about 7s
pred_wt_bt_pred <- nlraa::summary_simulate(t(pred_wt_bt$t))

pred_wt <- data.frame(method = "nls-bootstrap", 
                             "fork_length_mm" = LWdataA$fork_length_mm, 
                             fit = pred_wt_bt_pred[,1], 
                             lwr = pred_wt_bt_pred[,3],
                             upr = pred_wt_bt_pred[,4])

LWdataA <- LWdataA %>%
  bind_cols(select(pred_wt,-fork_length_mm)) %>%
  mutate(weight_g = round(ifelse(is.na(weight_g), fit, weight_g),0))

# Plot
ggplot()+
  geom_point(data = LWdata, aes(x = fork_length_mm, y = weight_g),
             shape = 21, fill = NA, color = "grey20")+
  geom_line(data = pred_log, aes(x = fork_length_mm, y =exp(fit)),
            color = "blue")+
  geom_ribbon(data = pred_log, aes(x = fork_length_mm, ymin = exp(lwr), ymax = exp(upr)),
              fill = "blue",alpha = 0.05)+
  geom_line(data = pred_wt, aes(x = fork_length_mm, y = fit),
            color = "tomato")+
  geom_ribbon(data = pred_wt, aes(x = fork_length_mm, ymin = lwr, ymax = upr),
              fill = "tomato",alpha = 0.1, inherit.aes = FALSE)+
  theme_classic()


# determine the age class of each individual sampled

#source(file.path("2. Code","Age_Analysis.R"))
stb.key <- read_rds(file.path("01_Data","Input","SB_ALKey.rds"))
FSA::alkPlot(stb.key)

## Assign Ages to all Striped Bass using age-length key ====================
stb.len <- LWdataA %>%
  mutate(fork_length_cm = fork_length_mm/10)

LWdataA <- FSA::alkIndivAge(key = stb.key, 
                                     formula = ~fork_length_cm, 
                                     data = stb.len, 
                                     type = "SR", seed = 42069)

saveRDS(LWdataA, file.path("01_Data","Output","LWdataA.rds"))
saveRDS(LW_nls, file.path("01_Data","Output","LW_nls.rds"))
