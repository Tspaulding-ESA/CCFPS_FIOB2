# for consumption: digitized data for cmax, 
# digitized data for temperature dependence

CmaxW <- read.csv(file.path("01_Data","Input","CmaxWt.csv"))
CtempA1 <- read.csv(file.path("01_Data","Input","CtempA1.csv"))
CtempA2 <- read.csv(file.path("01_Data","Input","CtempA2.csv"))
CtempA3 <- read.csv(file.path("01_Data","Input","CtempA3.csv"))

cparams <- data.frame(
  ck1 = c(0.262, 0.255, rep(0.323, 8)),
  ck4 = c(0.85, 0.9, rep(0.85, 8)),
  ctl = c(30, 32, rep(30, 8)),
  cto = c(19, 18, rep(15, 8)),
  cq = c(6.6, 6.6, rep(7.4, 8)),
  ctm = c(28, 29, rep(28, 8)))

pvals <- c(rep(0.69, 2), rep(0.73, 1), rep(0.68, 1), rep(0.62, 6))
refsize <- c(100, 350, rep(1000, 8))

source(file.path("02_Code","Functions.R"))

#consumption functions with variation

# refit the Cmax allometric relationship
ConAllo_log <- lm(log(Cmax) ~ log(Wt), data = CmaxW)
con_preds_log <- predict(ConAllo_log, newdata = list(Wt = 5:2700), 
                         interval = "prediction") %>%
  bind_cols(Wt = 5:2700)

ConAllo_nls <- nls(Cmax ~ a * Wt^b, 
                   data = na.omit(CmaxW), 
                   start = list(a = 0.302,
                                b = -0.252)) # using starting values from Juniper's Code

con_preds_nls <- xgxr::predict.nls(ConAllo_nls,
                                   newdata = data.frame("Wt" = 5:2700), 
                                   se.fit = TRUE, interval = "confidence", 
                                   level = 0.995)$fit %>%
  bind_cols(Wt = 5:2700)
names(con_preds_nls) <- c("fit","lwr","upr","Wt")

ggplot()+
  geom_ribbon(data = con_preds_log,
              aes(x = Wt, ymin = exp(lwr), ymax = exp(upr)),
              fill = "blue", alpha = 0.05)+
  geom_line(data = con_preds_log,
            aes(x = Wt, y = exp(fit)),
            color = "blue")+
  geom_ribbon(data = con_preds_nls,
              aes(x = Wt, ymin = lwr, ymax = upr),
              fill = "tomato", alpha = 0.05)+
  geom_line(data = con_preds_nls,
            aes(x = Wt, y = fit),
            color = "tomato")+
  geom_point(data = CmaxW,
             aes(x = Wt, y = Cmax),
             shape = 21)+
  theme_classic()+
  labs(x = "Wet Weight (g)", y = "Cmax (g/g/d)")

ggsave(file.path("01_Data","Output", "Figures", "Cmax.png"), device = "png", width = 6,
       height = 4, units = "in", dpi = 300)

# refit the age-specific temperature dependence relationships
#   because gam predict only gives the error of the fit, we expand the 
#   error to include residuals (to generate prediction interval)

# age 1 (100 g)
g1 <- mgcv::gam(C ~ s(Temp), data = CtempA1)
summary(g1)
xp <- seq(6, 30, 0.1)
gpy <- predict(g1, newdata = list(Temp = xp), type = 'response', 
               se.fit = T) %>%
  bind_cols("Temp" = xp)

gpy <- gpy %>%
  mutate(lwr = fit - 1.96 * (se.fit / sqrt(summary(g1)[10]$r.sq)),
         upr = fit + 1.96 * (se.fit / sqrt(summary(g1)[10]$r.sq)))


yp <- Cft(cmax = (0.3021 * 100^(-0.2523)), Te = xp, ck1 = 0.262, 
          ck4 = 0.85, ctl = 30, cto = 19, cq = 6.6, ctm = 28) %>%
  bind_cols("Temp" = xp)

ggplot()+
  geom_ribbon(data = gpy, aes(x = Temp, ymin = lwr, ymax = upr),
              fill = "grey90")+
  geom_point(data = CtempA1, aes(x = Temp,  y = C),
             shape = 21, size = 2)+
  geom_line(data = gpy, aes(x = Temp, y = fit))+
  geom_line(data = yp, aes(x = Temp, y = `...1`),
            color = "blue", linetype = "dashed")+
  theme_classic()

ggsave(file.path("01_Data","Output", "Figures", "CT_age_1.png"), device = "png", width = 6,
       height = 4, units = "in", dpi = 300)

# age 2 (350 g)

g2 <- mgcv::gam(C ~ s(Temp), data = CtempA2)
summary(g2)
xp <- seq(6, 30, 0.1)
gpy <- predict(g2, newdata = list(Temp = xp), type = 'response', 
               se.fit = T) %>%
  bind_cols("Temp" = xp)

gpy <- gpy %>%
  mutate(lwr = fit - 1.96 * (se.fit / sqrt(summary(g2)[10]$r.sq)),
         upr = fit + 1.96 * (se.fit / sqrt(summary(g2)[10]$r.sq)))

yp <- Cft(cmax = (0.3021 * 350^(-0.2523)), Te = xp, ck1 = 0.255, 
          ck4 = 0.9, ctl = 32, cto = 18, cq = 6.6, ctm = 29) %>%
  bind_cols("Temp" = xp)

ggplot()+
  geom_ribbon(data = gpy, aes(x = Temp, ymin = lwr, ymax = upr),
              fill = "grey90")+
  geom_point(data = CtempA1, aes(x = Temp,  y = C),
             shape = 21, size = 2)+
  geom_line(data = gpy, aes(x = Temp, y = fit))+
  geom_line(data = yp, aes(x = Temp, y = `...1`),
            color = "blue", linetype = "dashed")+
  theme_classic()

ggsave(file.path("01_Data","Output", "Figures", "CT_age_2.png"), device = "png", width = 6,
       height = 4, units = "in", dpi = 300)

# age 3 (1000 g)
g3 <- mgcv::gam(C ~ s(Temp), data = CtempA3)
summary(g3)
xp <- seq(6, 30, 0.1)
gpy <- predict(g3, newdata = list(Temp = xp), type = 'response', 
               se.fit = T) %>%
  bind_cols("Temp" = xp)

gpy <- gpy %>%
  mutate(lwr = fit - 1.96 * (se.fit / sqrt(summary(g3)[10]$r.sq)),
         upr = fit + 1.96 * (se.fit / sqrt(summary(g3)[10]$r.sq)))

yp <- Cft(cmax = (0.3021 * 1000^(-0.2523)), Te = xp, ck1 = 0.323, 
          ck4 = 0.85, ctl = 30, cto = 15, cq = 7.4, ctm = 28) %>%
  bind_cols("Temp" = xp)

ggplot()+
  geom_ribbon(data = gpy, aes(x = Temp, ymin = lwr, ymax = upr),
              fill = "grey90")+
  geom_point(data = CtempA1, aes(x = Temp,  y = C),
             shape = 21, size = 2)+
  geom_line(data = gpy, aes(x = Temp, y = fit))+
  geom_line(data = yp, aes(x = Temp, y = `...1`),
            color = "blue", linetype = "dashed")+
  theme_classic()

ggsave(file.path("01_Data","Output", "Figures", "CT_age_3.png"), device = "png", width = 6,
       height = 4, units = "in", dpi = 300)

saveRDS(g1, file.path("01_Data","Input","g1.rds"))
saveRDS(g2, file.path("01_Data","Input","g2.rds"))
saveRDS(g3, file.path("01_Data","Input","g3.rds"))
gs <- list(g1, g2, g3, g3, g3, g3, g3, g3, g3, g3)
# Build Lookup Tables for faster processing
# Predicted Cmax based on Age and Temp
pred_cmax_lu <- list()
for (i in 1:length(gs)){
  pred_cmax_lu[[i]] <- predict(gs[i][[1]], newdata = data.frame("Temp" = TE),
                               type = "response",se.fit = T) %>%
    bind_cols("age" = i,
              "TE" = TE)
}
pred_cmax_lu <- bind_rows(pred_cmax_lu) %>%
  mutate(TE = round(TE,2)) %>%
  group_by(age, TE) %>%
  summarise(fit = mean(fit),
            se.fit = mean(se.fit))

names(pred_cmax_lu) <- c("age", "TE","cmax_mean", "cmax_se.fit")
saveRDS(pred_cmax_lu, file.path("01_Data","Input","pred_cmax_lu.rds"))

# Residuals of the Predicted Cmax relationship
r_sq_lu <- list()
for (i in 1:length(gs)){
  r_sq_lu[[i]] <- sqrt(summary(gs[i][[1]])[10]$r.sq)
}
saveRDS(r_sq_lu, file.path("01_Data","Input","r_sq_lu.rds"))

# Reference Cmax
ref_cmax_lu <- list()
for (i in 1:length(gs)){
  ref_cmax_lu[[i]] <- xgxr::predict.nls(ConAllo,
                                        newdata = data.frame("Wt" = refsize[i]), 
                                        se.fit = TRUE, 
                                        interval = "confidence", 
                                        level = 0.995)$fit %>%
    mutate("age" = i,
           "ref_cmax_sds" = (fit-lwr)/1.96)
}
ref_cmax_lu <- bind_rows(ref_cmax_lu) %>%
  distinct() %>%
  select(age, "ref_cmax_mean" = fit, ref_cmax_sds)
saveRDS(ref_cmax_lu, file.path("01_Data","Input","ref_cmax_lu.rds"))

# Allometric Consumption Relationship
WW = round(seq(0,max(LWdataA$upr)+100,0.01),2)

con_allo_lu <- xgxr::predict.nls(ConAllo,
                                 newdata = data.frame("Wt" = WW),
                                 se.fit = TRUE, interval = "confidence",
                                 level = 0.95)$fit

con_allo_lu <- con_allo_lu %>%
  bind_cols("WW" = WW) %>%
  filter(!is.na(fit)) %>%
  mutate(cmaxi_sds = (fit-lwr)/1.96) %>%
  select(WW, "cmaxi_mean" = fit, cmaxi_sds)

saveRDS(con_allo_lu,file.path("01_Data","Input","ConAlloLU.rds"))


