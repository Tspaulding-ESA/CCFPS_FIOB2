###################### Estimate Age based on length #######################
library(tidyverse)
library(FSA)

## Retrieve Length at capture from release data ===========================
release_data <- read.csv(file.path("01_Data","Input","ReleaseData.csv"))

# Develop Age-length key based on a methods from 
# fishR Vignette (D.Ogle 2013), based on methods from 
# Iserman and Knight (2015)

## Striped Bass Aging =====================================================
# Building the age-length key table ####
set.seed(42069)

AgeData <- read.csv(file.path("01_Data","Input","AgedFish.csv"), header = T, stringsAsFactors = F)

stb.age <- AgeData[AgeData$Species == "Striped Bass"# & AgeData$Length_cm >20 & AgeData$Age < 8
                   ,] %>%
  select(Length_cm, Age) %>%
  arrange(Age)

x <- stb.age$Age
y <- stb.age$Length_cm

alk.lm <- lm(y ~ I(log(x)))

predicted.intervals <- data.frame("Age" = c(1,unique(stb.age$Age)),
                                  predict(alk.lm,data.frame(x=c(1,unique(stb.age$Age))),
                                          interval='confidence', 
                                          level=0.95))
names(predicted.intervals) <- c("Age","fit","lwr_95","upr_95")
predicted.intervals <- predicted.intervals %>%
  left_join(data.frame("Age" = c(1,unique(stb.age$Age)),
                                  predict(alk.lm,data.frame(x=c(1,unique(stb.age$Age))),
                                          interval='confidence', 
                                          level=0.995)))
names(predicted.intervals) <- c("Age","fit","lwr_95","upr_95","lwr_995","upr_995")

ggplot()+
  geom_ribbon(data = predicted.intervals,
              aes(x = Age, ymin = lwr_95, ymax = upr_95),
              fill = "lightblue")+
  geom_point(data = stb.age, aes(x = Age, y = Length_cm),
             shape = 21, fill = NA, color = "black", size = 2)+
  geom_line(data = predicted.intervals,
            aes(x = Age, y = fit),
            color = "blue")+
  geom_line(data = predicted.intervals,
            aes(x = Age, y = lwr_995),
            color = "grey20",
            linetype = "dashed")+
  geom_line(data = predicted.intervals,
            aes(x = Age, y = upr_995),
            color = "grey20",
            linetype = "dashed")+
  scale_x_continuous(name = "Age",
                     breaks = c(1:8))+
  scale_y_continuous(name = "Length (cm)")+
  theme_classic()

ggsave(file.path("01_Data","Output","Figures","Age-LengthKey.png"),
       width = 6, height = 4, units = "in", dpi = 300)

x_new <- rep(seq(1,8, by = 1), 100)

alk.lm$fitted.values <- predict(alk.lm, data.frame(x = x_new))

y_new <- simulate(alk.lm)[,1]

 plot(x_new, y_new, col = 'red')
 points(stb.age$Age, stb.age$Length_cm, pch = 16, )
 lines(c(1,unique(stb.age$Age)), predicted.intervals[,2], col = 'red', lwd = 3)
 lines(c(1,unique(stb.age$Age)), predicted.intervals[,3], col = 'black', lwd = 1, lty = 2)
 lines(c(1,unique(stb.age$Age)), predicted.intervals[,4], col = 'black', lwd = 1, lty = 2)

alk_new <- tibble(Age = x_new, Length_cm = y_new) %>%
  mutate()

# Find the minimum size for use in lencat function
FSA::Summarize(~Length_cm, data = alk_new, digits = 1)

# use lencat function to define length categories for each fish starting at min 24, and increasing by 2 
stb.age1 <- FSA::lencat(~Length_cm, data = alk_new, startcat=0, w=5)

stb.raw <- with(stb.age1, table(LCat,Age))
stb.key <- prop.table(stb.raw, margin = 1)
round(stb.key,2)
saveRDS(stb.key, file.path("01_Data","Input","SB_ALKey.rds"))

# Cleanup #################################################################
rm(alk_new, alk.lm, AgeData, predicted.intervals, release_data, stb.age, 
   stb.age1, stb.raw, stb.key, x, x_new, y, y_new)