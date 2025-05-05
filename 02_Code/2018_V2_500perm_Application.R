################################################################################
#
# Estimation of Potential Predation By Bass Removed from Clifton Court Forebay
#
# Code written by Juniper Simonis of DAPPER Stats
#
# 10/2017
#
# Adapted by Taylor Spaulding 
# 
# 05/2025
#
################################################################################
 
################################################################################
#
# TABLE OF CONTENTS
#
#  1. Preliminaries
#
#  2. Prep data and functions
#
#  3. Fixed Parameter Calculations
#
#  4. Variable Parameter Calculations
#
#  5. Additional Output
#
################################################################################

################################################################################
#
# 1. Preliminaries
#
################################################################################

  # source function scripts and load packages

    # source the script with the custom-written functions

      source(file.path("02_Code","Functions.R"))

    # mgcv is used for fitting and predicting with generalized additive models 

      # if(("mgcv" %in% rownames(installed.packages())) == F)
      #   install.packages("mgcv")
      library(mgcv)

  # load data files

    # for predator allometry: 
    #  field measurements, length at age 
    #  and the digitized data from the length-weight regression	

      fieldsizes <- read.csv(file.path("01_Data","Input","FieldSizes.csv"))
      #source(02_Age_Analysis)
      readRDS(file.path("01_Data","Input","SB_ALKey.rds"))
      LWdata <- read.csv(file.path("01_Data","Input","SBlw.csv"))

    # for predator diet: historic qpcr data, historic % data, historic FO data

      qpcr <- read.csv(file.path("01_Data","Input","QPCR.csv"))
      histPdiet <- read.csv(file.path("01_Data","Input","HistoricDiet.csv"))
      histFOdiet <- read.csv(file.path("01_Data","Input","HistoricFODiet.csv"))

    # temperature data
      
      TempT <- read.csv(file.path("01_Data","Input","Temps.csv")) %>%
        mutate(date = mdy(Date))
      
    # for consumption: digitized data for cmax, 
    # digitized data for temperature dependence

      CmaxW <- read.csv(file.path("01_Data","Input","CmaxWt.csv"))
      CtempA1 <- read.csv(file.path("01_Data","Input","CtempA1.csv"))
      CtempA2 <- read.csv(file.path("01_Data","Input","CtempA2.csv"))
      CtempA3 <- read.csv(file.path("01_Data","Input","CtempA3.csv"))
	
################################################################################
#
# 2. Prep Data and functions
#
################################################################################


  # Date window

      DW <- data.frame(Year = c(2018),  
              FirstDate = as.Date(c("2018-1-8")), 
              LastDate = as.Date(c("2018-5-5")))

  # generate mean temps for each day

      Year <- c(rep(2018, 118))
      DateX <- c(as.Date(as.Date("1-8-18", format = "%m-%d-%y"):
                        as.Date("5-5-18", format = "%m-%d-%y"), 
                        origin = "1970-01-01"))
      Date <- format(c(as.Date(as.Date("1-8-18", format = "%m-%d-%y"):
                        as.Date("5-5-18", format = "%m-%d-%y"), 
                        origin = "1970-01-01")), format = "%m/%d/%Y")
      Day <- c(1:118)

      TeTv <- rep(NA, length(Day))
      for(i in 1:length(Day)){

        TeTv[i] <- round(mean(TempT$Temp[which(TempT$Date == DateX[i])]),2)
      }

      TeT <- data.frame(Year, Date = DateX, Day, Temp = TeTv)	

  # predator l-w regression

    # fit the regression 	

      L <- LWdata$L
      W <- LWdata$W
      LWreg <- lm(log(W) ~ log(L))			

      # plot 	
		
        windows(8,6)
        par(mar = c(4.5,7,1,1))
        plot((W) ~ (L), xlim = c(100, 1100), ylim = c(0, 20500), xlab = "", 
             ylab = "", cex.axis = 1.5, las = 1, 
             col = rgb(0, 0, 0, alpha = 0.6), 
             lwd = 2, cex = 1.75)
        xx <- 110:1080
        preds <- predict(LWreg, newdata = list(L = xx), 
                          interval = "prediction")
        points((xx), exp(preds[, "fit"]), type = 'l', lwd = 6, 
                col = rgb(0.2, 0.2, 0.2, alpha = 0.8))
        points((xx), exp(preds[, "lwr"]), type = 'l', lwd = 3, lty = 3, 
                col = rgb(0.2, 0.2, 0.2, alpha = 0.8))
        points((xx), exp(preds[, "upr"]), type = 'l', lwd = 3, lty = 3, 
                col = rgb(0.2, 0.2, 0.2, alpha = 0.8))
        mtext(side = 2, "Weight (g)", line = 5, cex = 2.5)
        mtext(side = 1, "Length (mm)", line = 3, cex = 2)
        yy <- (3.28e-5) * (xx^2.85)
        points(yy ~ xx, type = 'l', col = rgb(0, 0.4, 0.7), lty = 2, lwd = 4)
        savePlot("Output/LengthWeight.tiff", type = "tiff")

  # predator diet frequency of occurrence

    # integrate the qpcr data through time

      qpcrT <- apply(qpcr[ , 4:15], 2, sum)

    # generate mean and sd of probabilities for each prey item on logit 
    #   scale and then back transformed

      FOm <- qpcrT
      FOse <- qpcrT
      FObtm <- qpcrT
      FObtlw <- qpcrT
      FObtup <- qpcrT
      StPs <- qpcrT	
	
      for(j in 2:length(qpcrT)){
        if(qpcrT[j] > 0){
          m1 <- NULL
          QPCRm <- glm(c(rep(0, qpcrT[1] - qpcrT[j]), 
                          rep(1, qpcrT[j])) ~ 1, 
                       family = 'binomial')
          FOm[j] <- as.numeric(coef(QPCRm))
          FOse[j] <- summary(QPCRm)[12]$coefficients[2]
        }else{
          FOm[j] <- -Inf
        }
        StPs[j] <- qpcrT[j] / qpcrT[1]
        FObtlw[j] <- ilogit(FOm[j] - 1.96 * FOse[j])
        FObtup[j] <- ilogit(FOm[j] + 1.96 * FOse[j])
        FObtm[j] <- ilogit(FOm[j])
      }

    # diet fractions for the fixed analyses

      pNd <- histPdiet$percNmean
      pWd <- histPdiet$percWTmean
      pFOd <- c(FObtm[2:12], histFOdiet$FOmean)
      RII <- pNd + pWd + pFOd
      dietf <- RII/sum(RII)
      names(dietf) <- histPdiet$Prey
      dietf

    # diet fractions for the variable analyses

      # number of permutations

        nperm <- 100000

      # set seed for reproducability

        set.seed(123)

      dietfVA <- matrix(NA, nrow = nperm, ncol = length(dietf))
      colnames(dietfVA) <- histPdiet$Prey
      for(i in 1:nperm){

        pNdi <- rnorm(14, histPdiet$percNmean, histPdiet$percNsd)
        pNdi[which(pNdi < 0)] <- 0
        pNdi[which(pNdi > 1)] <- 1
        pWdi <- rnorm(14, histPdiet$percWTmean, histPdiet$percWTsd)
        pWdi[which(pWdi < 0)] <- 0
        pWdi[which(pWdi > 1)] <- 1
        pFOdi <- c(ilogit(rnorm(11, as.numeric(FOm[2:12]), 
                   as.numeric(FOse[2:12]))), 
                   rnorm(3, histFOdiet$FOmean, histFOdiet$FOsd)) 
        pFOdi[which(pFOdi < 0)] <- 0
        pFOdi[which(pFOdi > 1)] <- 1
        RIIi <- pNdi + pWdi + pFOdi
        dietfVA[i, ] <- RIIi/sum(RIIi)
      }

      summary(dietfVA[, "CHN"])
      sd(dietfVA[, "CHN"])

    # for use, we just pull out the CHN fractions

      dietfVCHN <- dietfVA[, "CHN"]
      rm("dietfVA")

  # predator sizes and age classes

    # for the few removed predators that don't have sizes, we assume the 
    #   mean size of fish removed on that day (exp-mean-log)

      fieldsizesA <- fieldsizes
      tofillin <- which(is.na(fieldsizes$LengthMeasured) == T)

      for(i in 1:length(tofillin)){

        datefished <- fieldsizes$Date[tofillin[i]]
        otherlengths <- na.omit(fieldsizes$LengthMeasured[
                                   which(fieldsizes$Date == datefished)])
        assumedlength <- round(exp(mean(log(otherlengths))))
        fieldsizesA$LengthMeasured[tofillin[i]] <- assumedlength
      }

    # determine the age class of each individual sampled

      ACs <- rep(NA, nrow(fieldsizesA))

      for(i in 1:length(ACs)){
          ACs[i] <- which(laa$MinLength <= fieldsizesA$LengthMeasured[i] &
                          laa$MaxLength >= fieldsizesA$LengthMeasured[i])
      } 

  # consumption functions with variation

    # refit the Cmax allometric relationship

      cm <- CmaxW$Cmax
      wt <- CmaxW$Wt
      ConAllo <- lm(log(cm) ~ log(wt))
      preds <- predict(ConAllo, newdata = list(wt = 5:2700), 
                       interval = "prediction")

      windows(10,6)
      par(mar = c(4.5, 6, 1, 1))
      plot((cm) ~ (wt), xlim = c(0, (2800)), ylim=c(0, (0.3)), xlab = "", 
           ylab = "", cex.axis = 1.5, las = 1, col = rgb(0, 0, 0, alpha = 0.6), 
           lwd = 2, cex = 1.75)
      points((5:2700), exp(preds[, "fit"]), type = 'l', lwd = 6, 
           col = rgb(0.2, 0.2, 0.2, alpha = 1))
      points((5:2700), exp(preds[, "lwr"]), type = 'l', lwd = 3, lty = 3, 
           col=rgb(0.2, 0.2, 0.2, alpha = 1))
      points((5:2700), exp(preds[, "upr"]), type = 'l', lwd = 3, lty = 3, 
           col=rgb(0.2, 0.2, 0.2, alpha = 1))
      mtext(side = 2, "Cmax (g/g/d)", line = 4, cex = 2.5)
      mtext(side = 1, "Wet weight (g)", line = 3, cex = 2)
      yy <- 0.302 * (5:2700) ^ (-0.252)
      points((5:2700), yy, type = 'l', col=rgb(0.2, 0.6, 0.8, alpha = 1), 
           lty = 2, lwd = 4)
      savePlot("Output/Cmax.tiff", type="tiff")

    # refit the age-specific temperature dependence relationships
    #   because gam predict only gives the error of the fit, we expand the 
    #   error to include residuals (to generate prediction interval)

      # age 1 (100 g)

        cc <- CtempA1$C
        tt <- CtempA1$Temp
        g1 <- gam(cc ~ s(tt))
        summary(g1)
        xp <- seq(6, 30, 0.1)
        gpy <- predict(g1, newdata = list(tt = xp), type = 'response', 
                        se.fit = T)
        windows(8,6)
        par(mar = c(4.5, 6, 1, 1))
        plot(cc ~ tt, ylim = c(0, 0.2), xlim = c(5, 30), xlab = "", ylab = "", 
                cex.axis = 1.5, las = 1, col = rgb(0, 0, 0, alpha = 0.6), 
                lwd = 2, cex = 1.75)
        points(xp, gpy$fit, type='l', lwd=6, col=rgb(0.2,0.2,.2, alpha=1))
        points(xp, gpy$fit + 1.96 * (gpy$se.fit / sqrt(summary(g1)[10]$r.sq)), 
                type = 'l', lwd = 3, lty = 3, 
                col = rgb(0.2, 0.2, 0.2, alpha = 1))
        points(xp, gpy$fit - 1.96 * (gpy$se.fit / sqrt(summary(g1)[10]$r.sq)), 
                type = 'l', lwd = 3, lty = 3, 
                col = rgb(0.2, 0.2, 0.2, alpha = 1))
        mtext(side = 2, "Cmax (g/g/d)", line = 4, cex = 2.0)
        mtext(side = 1, "Temperature (C)", line = 3, cex = 1.75)
        yp <- Cft(cmax = (0.3021 * 100^(-0.2523)), Te = xp, ck1 = 0.262, 
                        ck4 = 0.85, ctl = 30, cto = 19, cq = 6.6, ctm = 28)
        points(yp ~ xp, type = 'l', col = rgb(0.2, 0.6, 0.8, alpha = 1), 
                lty = 2, lwd = 4)
        savePlot("Output/CTa1.tiff", type = "tiff")

      # age 2 (350 g)

        cc <- CtempA2$C
        tt <- CtempA2$Temp
        g2 <- gam(cc ~ s(tt))
        summary(g2)
        xp <- seq(6, 30, 0.1)
        gpy <- predict(g2, newdata = list(tt = xp), type = 'response', 
                        se.fit = T)
        windows(8,6)
        par(mar = c(4.5, 6, 1, 1))
        plot(cc ~ tt, ylim = c(0, 0.2), xlim = c(5, 30), xlab = "", ylab = "", 
                cex.axis = 1.5, las = 1, col = rgb(0, 0, 0, alpha = 0.6), 
                lwd = 2, cex = 1.75)
        points(xp, gpy$fit, type='l', lwd=6, col=rgb(0.2,0.2,.2, alpha=1))
        points(xp, gpy$fit + 1.96 * (gpy$se.fit / sqrt(summary(g2)[10]$r.sq)), 
                type = 'l', lwd = 3, lty = 3, 
                col = rgb(0.2, 0.2, 0.2, alpha = 1))
        points(xp, gpy$fit - 1.96 * (gpy$se.fit / sqrt(summary(g2)[10]$r.sq)), 
                type = 'l', lwd = 3, lty = 3, 
                col = rgb(0.2, 0.2, 0.2, alpha = 1))
        mtext(side = 2, "Cmax (g/g/d)", line = 4, cex = 2.0)
        mtext(side = 1, "Temperature (C)", line = 3, cex = 1.75)
        yp <- Cft(cmax = (0.3021 * 350^(-0.2523)), Te = xp, ck1 = 0.255, 
                        ck4 = 0.9, ctl = 32, cto = 18, cq = 6.6, ctm = 29)
        points(yp ~ xp, type = 'l', col = rgb(0.2, 0.6, 0.8, alpha = 1), 
                lty = 2, lwd = 4)
        savePlot("Output/CTa2.tiff", type = "tiff")

      # age 3 (1000 g)

        cc <- CtempA3$C
        tt <- CtempA3$Temp
        g3 <- gam(cc ~ s(tt))
        summary(g3)
        xp <- seq(6, 30, 0.1)
        gpy <- predict(g3, newdata = list(tt = xp), type = 'response', 
                        se.fit = T)
        windows(8,6)
        par(mar = c(4.5, 6, 1, 1))
        plot(cc ~ tt, ylim = c(0, 0.2), xlim = c(5, 30), xlab = "", ylab = "", 
                cex.axis = 1.5, las = 1, 
                col = rgb(0, 0, 0, alpha = 0.6), lwd = 2, cex = 1.75)
        points(xp, gpy$fit, type='l', lwd=6, col=rgb(0.2,0.2,.2, alpha=1))
        points(xp, gpy$fit + 1.96 * (gpy$se.fit / sqrt(summary(g3)[10]$r.sq)), 
                type = 'l', lwd = 3, lty = 3, 
                col = rgb(0.2, 0.2, 0.2, alpha = 1))
        points(xp, gpy$fit - 1.96 * (gpy$se.fit / sqrt(summary(g3)[10]$r.sq)), 
                type = 'l', lwd = 3, lty = 3, 
                col = rgb(0.2, 0.2, 0.2, alpha = 1))
        mtext(side = 2, "Cmax (g/g/d)", line = 4, cex = 2.0)
        mtext(side = 1, "Temperature (C)", line = 3, cex = 1.75)
        yp <- Cft(cmax = (0.3021 * 1000^(-0.2523)), Te = xp, ck1 = 0.323, 
                ck4 = 0.85, ctl = 30, cto = 15, cq = 7.4, ctm = 28)
        points(yp ~ xp, type = 'l', col = rgb(0.2, 0.6, 0.8, alpha = 1), 
                lty = 2, lwd = 4)
        savePlot("Output/CTa3.tiff", type = "tiff")


################################################################################
#
# 3. Fixed Parameter Calculations
#
################################################################################

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

      TCHNconsmat_plt1 <- matrix(0, nrow = length(ACs), ncol = length(Date))

      for(i in 1:length(ACs)){

        LL <- fieldsizesA$LengthMeasured[i]
        WW <- exp(predict(LWreg, newdata = list(L = LL)))

        cmax <- 0.3021 * WW ^ -0.2523
        p <- pvals[ACs[i]]

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



################################################################################
#
# 4. Variable Parameter Calculations
#
################################################################################

  # Total Basic with Direct Accumulation, p < 1

    # for each predator, calculate the amount of food they would have eaten
    #   over each day, accounting for temperature, but ignoring
    #   predator growth
    # repeat NPERM amount of times

    # because gam predict only gives the error of the fit, 
    #    we expand the error to include residuals/generate prediction interval

      set.seed(123)

      gs <- list(g1, g2, g3, g3, g3, g3, g3, g3, g3, g3)
      refsize <- c(100, 350, rep(1000, 8))

      NPERM <- 750

      TCHNconsmatVA_plt1 <- array(0, dim = c(length(ACs), length(Date), NPERM))

      for(i in 1:length(ACs)){

        LL <- fieldsizesA$LengthMeasured[i]
        PWs <- predict(LWreg, newdata = list(L = LL), se.fit = T,
                       interval = "prediction")
        WW <- rnorm(NPERM, exp(PWs$fit[1, "fit"]), 
                    (exp(PWs$fit[1, "fit"]) - exp(PWs$fit[1, "lwr"])) / 1.96)

        TE <- TeT$Temp

        predcmax <- predict(gs[ACs[i]][[1]], newdata = list(tt = TE), 
                       type = 'response', se.fit = T)
        predcmaxmeans <- (predcmax$fit)
        predcmaxsds <- (predcmax$se.fit/sqrt(summary(gs[ACs[i]][[1]])[10]$r.sq))

        predcmaxref <- predict(ConAllo, newdata = list(wt = refsize[ACs[i]]),
                           se.fit = T, interval = "prediction")
        predcmaxrefmeans <- exp(predcmaxref$fit[1, "fit"])
        predcmaxrefsds <- (exp(predcmaxref$fit[1, "fit"]) - 
                           exp(predcmaxref$fit[1, "lwr"])) / 1.96


        for(j in 1:length(Date)){

          predcmaxi <- predict(ConAllo, newdata = list(wt = WW),
                               se.fit = T, interval = "prediction")
          predcmaximeans <- exp(predcmaxi$fit[, "fit"])
          predcmaxisds <- (exp(predcmaxi$fit[, "fit"]) - 
                           exp(predcmaxi$fit[, "lwr"])) / 1.96

          TCHNconsmatVA_plt1[i, j, ] <- rnorm(NPERM, predcmaxmeans[j],
                                   predcmaxsds[j]) / 
                                  rnorm(NPERM, predcmaxrefmeans, 
                                   predcmaxrefsds) * 
                                  rnorm(NPERM, predcmaximeans, 
                                   predcmaxisds) * 
                                  WW * pvals[ACs[i]] *
                                  sample(dietfVCHN, NPERM, replace = T)
        }

      }

      for(i in 1:length(ACs)){
        for(j in 1:length(Date)){
          TCHNconsmatVA_plt1[i, j, which(TCHNconsmatVA_plt1[i, j, ] < 0)] <- 0
        }
      }


      CHNconsmatVA_plt1 <- TCHNconsmatVA_plt1 

      for(i in 1:length(ACs)){
        SpecDate <- as.Date(fieldsizesA$Date[i], format = "%m/%d/%Y")
        ZDs <- which(TeT$Date <= SpecDate)
        CHNconsmatVA_plt1[i, ZDs, ] <- 0   
      }


  # Total Basic with Direct Accumulation, p = 1

    # for each predator, calculate the amount of food they would have eaten
    #   over each day, accounting for temperature, but ignoring
    #   predator growth
    # repeat NPERM amount of times

    # because gam predict only gives the error of the fit, 
    #    we expand the error to include residuals/generate prediction interval

      set.seed(456)

      TCHNconsmatVA_p1 <- array(0, dim = c(length(ACs), length(Date), NPERM))

      for(i in 1:length(ACs)){

        LL <- fieldsizesA$LengthMeasured[i]
        PWs <- predict(LWreg, newdata = list(L = LL), se.fit = T,
                       interval = "prediction")
        WW <- rnorm(NPERM, exp(PWs$fit[1, "fit"]), 
                    (exp(PWs$fit[1, "fit"]) - exp(PWs$fit[1, "lwr"])) / 1.96)

        TE <- TeT$Temp

        predcmax <- predict(gs[ACs[i]][[1]], newdata = list(tt = TE), 
                       type = 'response', se.fit = T)
        predcmaxmeans <- (predcmax$fit)
        predcmaxsds <- (predcmax$se.fit/sqrt(summary(gs[ACs[i]][[1]])[10]$r.sq))

        predcmaxref <- predict(ConAllo, newdata = list(wt = refsize[ACs[i]]),
                           se.fit = T, interval = "prediction")
        predcmaxrefmeans <- exp(predcmaxref$fit[1, "fit"])
        predcmaxrefsds <- (exp(predcmaxref$fit[1, "fit"]) - 
                           exp(predcmaxref$fit[1, "lwr"])) / 1.96


        for(j in 1:length(Date)){

          predcmaxi <- predict(ConAllo, newdata = list(wt = WW),
                               se.fit = T, interval = "prediction")
          predcmaximeans <- exp(predcmaxi$fit[, "fit"])
          predcmaxisds <- (exp(predcmaxi$fit[, "fit"]) - 
                           exp(predcmaxi$fit[, "lwr"])) / 1.96

          TCHNconsmatVA_p1[i, j, ] <- rnorm(NPERM, predcmaxmeans[j],
                                   predcmaxsds[j]) / 
                                  rnorm(NPERM, predcmaxrefmeans, 
                                   predcmaxrefsds) * 
                                  rnorm(NPERM, predcmaximeans, 
                                   predcmaxisds) * 
                                  WW * 1 *
                                  sample(dietfVCHN, NPERM, replace = T)
        }

      }

      for(i in 1:length(ACs)){
        for(j in 1:length(Date)){
          TCHNconsmatVA_p1[i, j, which(TCHNconsmatVA_p1[i, j, ] < 0)] <- 0
        }
      }


      CHNconsmatVA_p1 <- TCHNconsmatVA_p1 

      for(i in 1:length(ACs)){
        SpecDate <- as.Date(fieldsizesA$Date[i], format = "%m/%d/%Y")
        ZDs <- which(TeT$Date <= SpecDate)
        CHNconsmatVA_p1[i, ZDs, ] <- 0   
      }


################################################################################
#
# 5. Additional Output
#
################################################################################



     CHNconsavoidVA_p1_l95 <- rep(NA, length(Date))
     CHNconsavoidVA_p1_u95 <- rep(NA, length(Date))
     CHNconsavoidVA_p1_m <- rep(NA, length(Date))
     for(i in 1:length(Date)){
       tt <- apply(CHNconsmatVA_p1[, i,], 2, sum)
       CHNconsavoidVA_p1_l95[i] <- quantile(tt, 0.025)
       CHNconsavoidVA_p1_u95[i] <- quantile(tt, 0.975)
       CHNconsavoidVA_p1_m[i] <- mean(tt)
     }

     
    CHNconsrem_p1 <- data.frame(Date, 
                       round(data.frame(L95 = CHNconsavoidVA_p1_l95, 
                                        Mean = CHNconsavoidVA_p1_m, 
                                        U95 = CHNconsavoidVA_p1_u95), 3))

    write.csv(CHNconsrem_p1, "Output/CHNConsRemovedp1Var.csv", row.names = F)


     CHNconsTVA_p1_l95 <- rep(NA, length(Date))
     CHNconsTVA_p1_u95 <- rep(NA, length(Date))
     CHNconsTVA_p1_m <- rep(NA, length(Date))
     for(i in 1:length(Date)){
       tt <- apply(TCHNconsmatVA_p1[, i,], 2, sum)
       CHNconsTVA_p1_l95[i] <- quantile(tt, 0.025)
       CHNconsTVA_p1_u95[i] <- quantile(tt, 0.975)
       CHNconsTVA_p1_m[i] <- mean(tt)
     }

     
    CHNconstot_p1 <- data.frame(Date, 
                       round(data.frame(L95 = CHNconsTVA_p1_l95, 
                                        Mean = CHNconsTVA_p1_m, 
                                        U95 = CHNconsTVA_p1_u95), 3))

    write.csv(CHNconstot_p1, "Output/CHNConsTotp1Var.csv", row.names = F)



     CHNconsavoidVA_plt1_l95 <- rep(NA, length(Date))
     CHNconsavoidVA_plt1_u95 <- rep(NA, length(Date))
     CHNconsavoidVA_plt1_m <- rep(NA, length(Date))
     for(i in 1:length(Date)){
       tt <- apply(CHNconsmatVA_plt1[, i,], 2, sum)
       CHNconsavoidVA_plt1_l95[i] <- quantile(tt, 0.025)
       CHNconsavoidVA_plt1_u95[i] <- quantile(tt, 0.975)
       CHNconsavoidVA_plt1_m[i] <- mean(tt)
     }

     
    CHNconsrem_plt1 <- data.frame(Date, 
                       round(data.frame(L95 = CHNconsavoidVA_plt1_l95, 
                                        Mean = CHNconsavoidVA_plt1_m, 
                                        U95 = CHNconsavoidVA_plt1_u95), 3))

    write.csv(CHNconsrem_plt1, "Output/CHNConsRemovedplt1Var.csv", 
              row.names = F)


     CHNconsTVA_plt1_l95 <- rep(NA, length(Date))
     CHNconsTVA_plt1_u95 <- rep(NA, length(Date))
     CHNconsTVA_plt1_m <- rep(NA, length(Date))
     for(i in 1:length(Date)){
       tt <- apply(TCHNconsmatVA_plt1[, i,], 2, sum)
       CHNconsTVA_plt1_l95[i] <- quantile(tt, 0.025)
       CHNconsTVA_plt1_u95[i] <- quantile(tt, 0.975)
       CHNconsTVA_plt1_m[i] <- mean(tt)
     }

     
    CHNconstot_plt1 <- data.frame(Date, 
                       round(data.frame(L95 = CHNconsTVA_plt1_l95, 
                                        Mean = CHNconsTVA_plt1_m, 
                                        U95 = CHNconsTVA_plt1_u95), 3))

    write.csv(CHNconstot_plt1, "Output/CHNConsTotplt1Var.csv", row.names = F)




    TotalFullVA_p1_m <- matrix(NA, length(ACs), length(Date))

    for(i in 1:length(Date)){

      TotalFullVA_p1_m[ , i] <- apply(TCHNconsmatVA_p1[, i, ], 1, mean)
    }

    write.csv(TotalFullVA_p1_m, "Output/FullCHNConsTotp1Var_m.csv", 
              row.names = F)


    RemFullVA_p1_m <- matrix(NA, length(ACs), length(Date))

    for(i in 1:length(Date)){

      RemFullVA_p1_m[ , i] <- apply(CHNconsmatVA_p1[, i, ], 1, mean)
    }

    write.csv(RemFullVA_p1_m, "Output/FullCHNConsRemovedp1Var_m.csv", 
               row.names = F)



    TotalFullVA_plt1_m <- matrix(NA, length(ACs), length(Date))

    for(i in 1:length(Date)){

      TotalFullVA_plt1_m[ , i] <- apply(TCHNconsmatVA_plt1[, i, ], 1, mean)
    }

    write.csv(TotalFullVA_plt1_m, "Output/FullCHNConsTotplt1Var_m.csv", 
              row.names = F)


    RemFullVA_plt1_m <- matrix(NA, length(ACs), length(Date))

    for(i in 1:length(Date)){

      RemFullVA_plt1_m[ , i] <- apply(CHNconsmatVA_plt1[, i, ], 1, mean)
    }

    write.csv(RemFullVA_plt1_m, "Output/FullCHNConsRemovedplt1Var_m.csv", 
                row.names = F)




    TotalFullVA_p1_u95 <- matrix(NA, length(ACs), length(Date))

    for(i in 1:length(Date)){

      TotalFullVA_p1_u95[ , i] <- apply(TCHNconsmatVA_p1[, i, ], 1, 
                                          quantile, 0.975)
    }

    write.csv(TotalFullVA_p1_u95, "Output/FullCHNConsTotp1Var_u95.csv", 
                                        row.names = F)


    RemFullVA_p1_u95 <- matrix(NA, length(ACs), length(Date))

    for(i in 1:length(Date)){

      RemFullVA_p1_u95[ , i] <- apply(CHNconsmatVA_p1[, i, ], 1,
                                         quantile, 0.975)
    }

    write.csv(RemFullVA_p1_u95, "Output/FullCHNConsRemovedp1Var_u95.csv", 
                                          row.names = F)



    TotalFullVA_plt1_u95 <- matrix(NA, length(ACs), length(Date))

    for(i in 1:length(Date)){

      TotalFullVA_plt1_u95[ , i] <- apply(TCHNconsmatVA_plt1[, i, ], 1, 
                                        quantile, 0.975)
    }

    write.csv(TotalFullVA_plt1_u95, "Output/FullCHNConsTotplt1Var_u95.csv", 
                                         row.names = F)


    RemFullVA_plt1_u95 <- matrix(NA, length(ACs), length(Date))

    for(i in 1:length(Date)){

      RemFullVA_plt1_u95[ , i] <- apply(CHNconsmatVA_plt1[, i, ], 1, 
                                       quantile, 0.975)
    }

    write.csv(RemFullVA_plt1_u95, "Output/FullCHNConsRemovedplt1Var_u95.csv", 
                                           row.names = F)



    TCHNconsmatVA_p1_total <- rep(NA, NPERM)
    CHNconsmatVA_p1_total <- rep(NA, NPERM)
    TCHNconsmatVA_plt1_total <- rep(NA, NPERM)
    CHNconsmatVA_plt1_total <- rep(NA, NPERM)


    for(i in 1:NPERM){

      TCHNconsmatVA_p1_total[i] <- sum(TCHNconsmatVA_p1[ , , i])
      CHNconsmatVA_p1_total[i] <- sum(CHNconsmatVA_p1[ , , i])
      TCHNconsmatVA_plt1_total[i] <- sum(TCHNconsmatVA_plt1[ , , i])
      CHNconsmatVA_plt1_total[i] <- sum(CHNconsmatVA_plt1[ , , i])
    }

    totaltable <- cbind(rep(c("1", "lessthan1"), each = 2), 
                        rep(c("Total", "Removed"), 2),
                        round(rbind(
                              c(quantile(TCHNconsmatVA_p1_total, 0.025), 
                                mean(TCHNconsmatVA_p1_total),
                                quantile(TCHNconsmatVA_p1_total, 0.975)),
                              c(quantile(CHNconsmatVA_p1_total, 0.025), 
                                mean(CHNconsmatVA_p1_total),
                                quantile(CHNconsmatVA_p1_total, 0.975)),
                              c(quantile(TCHNconsmatVA_plt1_total, 0.025), 
                                mean(TCHNconsmatVA_plt1_total),
                                quantile(TCHNconsmatVA_plt1_total, 0.975)),
                              c(quantile(CHNconsmatVA_plt1_total, 0.025), 
                                mean(CHNconsmatVA_plt1_total),
                                quantile(CHNconsmatVA_plt1_total, 0.975))),
                              1)
                          )

    colnames(totaltable) <- c("p", "Type", "Lower95CI", "Mean", "Upper95CI")
    write.csv(totaltable, "Output/totalconsumptionestimates.csv", 
               row.names = F)


    CHNconsmatVA_p1_total_size <- matrix(NA, NPERM, 10)
    TCHNconsmatVA_p1_total_size <- matrix(NA, NPERM, 10)
    CHNconsmatVA_plt1_total_size <- matrix(NA, NPERM, 10)
    TCHNconsmatVA_plt1_total_size <- matrix(NA, NPERM, 10)

    for(i in 1:NPERM){

      for(j in 1:10){
 
        spots <- which(ACs == j)
        CHNconsmatVA_p1_total_size[i, j] <- 
                                 sum(CHNconsmatVA_p1[spots , , i])
        TCHNconsmatVA_p1_total_size[i, j] <- 
                                 sum(TCHNconsmatVA_p1[spots , , i])
        CHNconsmatVA_plt1_total_size[i, j] <-
                                sum(CHNconsmatVA_plt1[spots , , i])
        TCHNconsmatVA_plt1_total_size[i, j] <-
                                sum(TCHNconsmatVA_plt1[spots , , i])
      }
    }




    windows(14, 10)
    par(mar = c(5, 8, 1, 1))
    plot(1, 1, type = 'n', bty = 'L', xlab = "", ylab = "", 
         xaxt = "n", yaxt = "n", ylim = c(0, 110000),
         xlim = c(0.5, 10.5))
    for(i in 1:10){

      points(runif(NPERM, i - 0.2, i + 0.2), 
             CHNconsmatVA_p1_total_size[ , i],
             cex = 1.5, col = rgb(0.7, 0.7, 0.7, 0.5))
      points(c(i, i), quantile(CHNconsmatVA_p1_total_size[ , i], 
            c(0.25, 0.975)), lwd = 2, type = 'l')
      points(i, mean(CHNconsmatVA_p1_total_size[ , i]), cex = 1.5, 
             pch = 16, col = "white")
      points(i, mean(CHNconsmatVA_p1_total_size[ , i]), cex = 1.5, 
             pch = 1, col = 1, lwd = 3)

      spots <- which(ACs == i)
      yy <- sum(CHNconsmat_p1[spots, ])

      points(c(i - 0.25, i + 0.25), rep(yy, 2), type = 'l', 
             col = rgb(0, 0.1, 0.8, 0.5), lwd = 3)
    }

    axis(1, at = 1:10, cex.axis = 1.5)
    mtext(side = 1, line = 3.25, "Age Class", cex = 2)
    axis(2, las = 1, cex.axis = 1.5, at = seq(0, 100000, 20000),
         labels = c("0", "20000", "40000", "60000", "80000", "100000"))
    mtext(side = 2, line = 5.5, "Total Chinook Consumed (g)", cex = 2)
    savePlot("Output/AgeBasedRemovedConsp1.tiff", type = "tiff")



    plot(1, 1, type = 'n', bty = 'L', xlab = "", ylab = "", 
         xaxt = "n", yaxt = "n", ylim = c(0, 60000),
         xlim = c(0.5, 10.5))
    for(i in 1:10){

      points(runif(NPERM, i - 0.2, i + 0.2), 
             CHNconsmatVA_plt1_total_size[ , i],
             cex = 1.5, col = rgb(0.7, 0.7, 0.7, 0.5))
      points(c(i, i), quantile(CHNconsmatVA_plt1_total_size[ , i], 
            c(0.25, 0.975)), lwd = 2, type = 'l')
      points(i, mean(CHNconsmatVA_plt1_total_size[ , i]), cex = 1.5, 
             pch = 16, col = "white")
      points(i, mean(CHNconsmatVA_plt1_total_size[ , i]), cex = 1.5, 
             pch = 1, col = 1, lwd = 3)

      spots <- which(ACs == i)
      yy <- sum(CHNconsmat_plt1[spots, ])

      points(c(i - 0.25, i + 0.25), rep(yy, 2), type = 'l', 
             col = rgb(0, 0.1, 0.8, 0.5), lwd = 3)
    }

    axis(1, at = 1:10, cex.axis = 1.5)
    mtext(side = 1, line = 3.25, "Age Class", cex = 2)
    axis(2, las = 1, cex.axis = 1.5, at = seq(0, 60000, 10000),
         labels = c("0", "10000", "20000", "30000", 
                    "40000", "50000", "60000"))
    mtext(side = 2, line = 5.5, "Total Chinook Consumed (g)", cex = 2)
    savePlot("Output/AgeBasedRemovedConsplt1.tiff", type = "tiff")



    plot(1, 1, type = 'n', bty = 'L', xlab = "", ylab = "", 
         xaxt = "n", yaxt = "n", ylim = c(0, 300000),
         xlim = c(0.5, 10.5))
    for(i in 1:10){

      points(runif(NPERM, i - 0.2, i + 0.2), 
             TCHNconsmatVA_p1_total_size[ , i],
             cex = 1.5, col = rgb(0.7, 0.7, 0.7, 0.5))
      points(c(i, i), quantile(TCHNconsmatVA_p1_total_size[ , i], 
            c(0.25, 0.975)), lwd = 2, type = 'l')
      points(i, mean(TCHNconsmatVA_p1_total_size[ , i]), cex = 1.5, 
             pch = 16, col = "white")
      points(i, mean(TCHNconsmatVA_p1_total_size[ , i]), cex = 1.5, 
             pch = 1, col = 1, lwd = 3)

      spots <- which(ACs == i)
      yy <- sum(TCHNconsmat_p1[spots, ])

      points(c(i - 0.25, i + 0.25), rep(yy, 2), type = 'l', 
             col = rgb(0, 0.1, 0.8, 0.5), lwd = 3)
    }

    axis(1, at = 1:10, cex.axis = 1.5)
    mtext(side = 1, line = 3.25, "Age Class", cex = 2)
    axis(2, las = 1, cex.axis = 1.5, at = seq(0, 300000, 50000),
         labels = c("0", "50000", "100000", "150000", "200000", 
                      "250000", "300000"))
    mtext(side = 2, line = 5.5, "Total Chinook Consumed (g)", cex = 2)
    savePlot("Output/AgeBasedTotalConsp1.tiff", type = "tiff")



    plot(1, 1, type = 'n', bty = 'L', xlab = "", ylab = "", 
         xaxt = "n", yaxt = "n", ylim = c(0, 200000),
         xlim = c(0.5, 10.5))
    for(i in 1:10){

      points(runif(NPERM, i - 0.2, i + 0.2), 
             TCHNconsmatVA_plt1_total_size[ , i],
             cex = 1.5, col = rgb(0.7, 0.7, 0.7, 0.5))
      points(c(i, i), quantile(TCHNconsmatVA_plt1_total_size[ , i], 
            c(0.25, 0.975)), lwd = 2, type = 'l')
      points(i, mean(TCHNconsmatVA_plt1_total_size[ , i]), cex = 1.5, 
             pch = 16, col = "white")
      points(i, mean(TCHNconsmatVA_plt1_total_size[ , i]), cex = 1.5, 
             pch = 1, col = 1, lwd = 3)

      spots <- which(ACs == i)
      yy <- sum(TCHNconsmat_plt1[spots, ])

      points(c(i - 0.25, i + 0.25), rep(yy, 2), type = 'l', 
             col = rgb(0, 0.1, 0.8, 0.5), lwd = 3)
    }

    axis(1, at = 1:10, cex.axis = 1.5)
    mtext(side = 1, line = 3.25, "Age Class", cex = 2)
    axis(2, las = 1, cex.axis = 1.5, at = seq(0, 200000, 50000),
         labels = c("0",  "50000", "100000", "150000", "200000"))
    mtext(side = 2, line = 5.5, "Total Chinook Consumed (g)", cex = 2)
    savePlot("Output/AgeBasedTotalConsplt1.tiff", type = "tiff")




