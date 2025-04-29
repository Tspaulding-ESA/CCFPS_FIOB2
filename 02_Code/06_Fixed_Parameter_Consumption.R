################################################################################
#
# 3. Fixed Parameter Calculations
#
################################################################################

LWdataA <- readRDS(file.path("01_Data","Output","LWdataA.rds"))

# set up parameters based on age classes

cparams <- data.frame(
  ck1 = c(0.262, 0.255, rep(0.323, 8)),
  ck4 = c(0.85, 0.9, rep(0.85, 8)),
  ctl = c(30, 32, rep(30, 8)),
  cto = c(19, 18, rep(15, 8)),
  cq = c(6.6, 6.6, rep(7.4, 8)),
  ctm = c(28, 29, rep(28, 8)))

pvals <- c(rep(0.69, 2), rep(0.73, 1), rep(0.68, 1), rep(0.62, 6))

# Basic Consumption, p < 1

# for each predator, calculate the amount of food they would have eaten
#   over each day, accounting for temperature, but ignoring growth

TCHNconsmat_plt1 <- matrix(0, nrow = nrow(LWdataA), ncol = nrow(dates))

for(i in 1:nrow(LWdataA)){
  
  cmax <- 0.3021 * LWdataA$weight_g[i] ^ -0.2523
  p <- pvals[LWdataA$age[i]]
  
  ck1 <- cparams[LWdataA$age[i], "ck1"] 
  ck4 <- cparams[LWdataA$age[i], "ck4"] 
  ctl <- cparams[LWdataA$age[i], "ctl"] 
  cto <- cparams[LWdataA$age[i], "cto"] 
  cq <- cparams[LWdataA$age[i], "cq"] 
  ctm <- cparams[LWdataA$age[i], "ctm"] 
  Te <- TeT$Temp
  totalcons <- cmax * p *
    (ck1 * exp((1 / (cto - cq)) * 
                 log(0.98 * ((1 - ck1) / (0.02 * ck1))) *
                 (Te - cq))) / 
    (1 + ck1 * (exp((1 / (cto - cq)) * 
                      log(0.98 * ((1 - ck1) / (0.02 * ck1))) * 
                      (Te - cq)) - 1)) * 
    (ck4 * exp((1 / (ctl - ctm)) * 
                 log(0.98 * ((1 - ck4) / (0.02 * ck4))) * 
                 (ctl - Te))) / 
    (1 + ck4 * (exp((1 / (ctl - ctm)) * 
                      log(0.98 * ((1 - ck4) / (0.02 * ck4))) * 
                      (ctl - Te)) - 1)) * WW
  chncons <- totalcons * dietf["CHN"]
  
  TCHNconsmat_plt1[i, ] <- chncons
  
}


TotalBasicDirectOutput_plt1 <- data.frame(Date, gramsCHNconsumed = 
                                            round(apply(TCHNconsmat_plt1, 2, sum), 3))
write.csv(TotalBasicDirectOutput_plt1, 
          "Output/TotalBasicDirectCHNConsumption_plt1.csv", 
          row.names = F)

TotalFullBasicDirectOutput_plt1 <- data.frame(round(TCHNconsmat_plt1, 4))
colnames(TotalFullBasicDirectOutput_plt1) <- Date
write.csv(TotalFullBasicDirectOutput_plt1, 
          "Output/TotalFullBasicDirectCHNConsumption_plt1.csv", 
          row.names = F)


CHNconsmat_plt1 <- TCHNconsmat_plt1

for(i in 1:length(ACs)){
  SpecDate <- as.Date(fieldsizesA$Date[i], format = "%m/%d/%Y")
  ZDs <- which(TeT$Date <= SpecDate)
  CHNconsmat_plt1[i, ZDs] <- 0   
}



BasicDirectOutput_plt1 <- data.frame(Date, gramsCHNconsumed = 
                                       round(apply(CHNconsmat_plt1, 2, sum), 3))
write.csv(BasicDirectOutput_plt1, 
          "Output/BasicDirectCHNConsumption_plt1.csv", 
          row.names = F)

FullBasicDirectOutput_plt1 <- data.frame(round(CHNconsmat_plt1, 4))
colnames(FullBasicDirectOutput_plt1) <- Date
write.csv(FullBasicDirectOutput_plt1, 
          "Output/FullBasicDirectCHNConsumption_plt1.csv", 
          row.names = F)

# Basic Consumption, p = 1.0

# for each predator, calculate the amount of food they would have eaten
#   over each day, accounting for temperature, but ignoring growth
#   - assuming that p = 1 (predators would eat at max rate)

TCHNconsmat_p1 <- matrix(0, nrow = length(ACs), ncol = length(Date))

for(i in 1:length(ACs)){
  
  LL <- fieldsizesA$LengthMeasured[i]
  WW <- exp(predict(LWreg, newdata = list(L = LL)))
  
  cmax <- 0.3021 * WW ^ -0.2523
  p <- 1
  
  ck1 <- cparams[ACs[i], "ck1"] 
  ck4 <- cparams[ACs[i], "ck4"] 
  ctl <- cparams[ACs[i], "ctl"] 
  cto <- cparams[ACs[i], "cto"] 
  cq <- cparams[ACs[i], "cq"] 
  ctm <- cparams[ACs[i], "ctm"] 
  Te <- TeT$Temp
  totalcons <- cmax * p *
    (ck1 * exp((1 / (cto - cq)) * 
                 log(0.98 * ((1 - ck1) / (0.02 * ck1))) *
                 (Te - cq))) / 
    (1 + ck1 * (exp((1 / (cto - cq)) * 
                      log(0.98 * ((1 - ck1) / (0.02 * ck1))) * 
                      (Te - cq)) - 1)) * 
    (ck4 * exp((1 / (ctl - ctm)) * 
                 log(0.98 * ((1 - ck4) / (0.02 * ck4))) * 
                 (ctl - Te))) / 
    (1 + ck4 * (exp((1 / (ctl - ctm)) * 
                      log(0.98 * ((1 - ck4) / (0.02 * ck4))) * 
                      (ctl - Te)) - 1)) * WW
  chncons <- totalcons * dietf["CHN"]
  
  TCHNconsmat_p1[i, ] <- chncons
  
}


TotalBasicDirectOutput_p1 <- data.frame(Date, gramsCHNconsumed = 
                                          round(apply(TCHNconsmat_p1, 2, sum), 3))
write.csv(TotalBasicDirectOutput_p1, 
          "Output/TotalBasicDirectCHNConsumption_p1.csv", 
          row.names = F)

TotalFullBasicDirectOutput_p1 <- data.frame(round(TCHNconsmat_p1, 4))
colnames(TotalFullBasicDirectOutput_p1) <- Date
write.csv(TotalFullBasicDirectOutput_p1, 
          "Output/TotalFullBasicDirectCHNConsumption_p1.csv", 
          row.names = F)


CHNconsmat_p1 <- TCHNconsmat_p1

for(i in 1:length(ACs)){
  SpecDate <- as.Date(fieldsizesA$Date[i], format = "%m/%d/%Y")
  ZDs <- which(TeT$Date <= SpecDate)
  CHNconsmat_p1[i, ZDs] <- 0   
}



BasicDirectOutput_p1 <- data.frame(Date, gramsCHNconsumed = 
                                     round(apply(CHNconsmat_p1, 2, sum), 3))
write.csv(BasicDirectOutput_p1, 
          "Output/BasicDirectCHNConsumption_p1.csv", 
          row.names = F)

FullBasicDirectOutput_p1 <- data.frame(round(CHNconsmat_p1, 4))
colnames(FullBasicDirectOutput_p1) <- Date
write.csv(FullBasicDirectOutput_p1, 
          "Output/FullBasicDirectCHNConsumption_p1.csv", 
          row.names = F)
