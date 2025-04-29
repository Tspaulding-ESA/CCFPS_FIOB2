################################################################################
#
# Fish Bioenergetics Model Functions
#
# Code written by Juniper Simonis of DAPPER Stats
#
# 10/2017
#
################################################################################

################################################################################
#
# TABLE OF CONTENTS
#
# 1. General functions
# 2. Bioenergetics functions
#   a. consumption
#   b. respiration
#   c. Egestion
#   d. Excretion
#   e. Specific dynamic action
# 3. Modeling functions
#   a. daily growth
#   b. time series of growth
#   c. estimation of p value
#   d. daily growth with consumption
#   e. time series of growth with consumption
#   f. weights
# 4. Sensitivity analysis functions
#   a. individual
#   b. cohort
#   c. no variance
#
################################################################################

################################################################################
#
# 1. General functions
#
################################################################################


  # logit 

    logit <- function(x, ...){
        log(x / (1 - x))
    }

  # inverse logit

    ilogit <- function(x, ...){
        exp(x) / (1 + exp(x))
    }

################################################################################
#
# 2. Bioenergetics functions
#
################################################################################

  # a. Consumption
  
    # calculate consumption quickly based on parameter values and cmax

      Cft <- function(cmax, Te, ck1, ck4, ctl, cto, cq, ctm, ...){

          output <- cmax * 
                    (ck1 * exp((1 / (cto - cq)) * 
                               log(0.98 * ((1 - ck1) / (0.02 * ck1))) *
                               (Te-cq))) / 
                    (1 + ck1 * (exp((1 / (cto - cq)) * 
                                    log(0.98 * ((1 - ck1) / (0.02 * ck1))) * 
                                    (Te - cq)) - 1)) * 
                    (ck4 * exp((1 / (ctl - ctm)) * 
                               log(0.98 * ((1 - ck4) / (0.02 * ck4))) * 
                               (ctl - Te))) / 
                    (1 + ck4 * (exp((1 / (ctl - ctm)) * 
                                    log(0.98 * ((1 - ck4) / (0.02 * ck4))) * 
                                (ctl - Te)) - 1))
       
          return(output)

      }

    # consumption function for the bioenergetics sensitivity analysis

      Consumption <- function(W, p, temp, Allo, CfT, SAToggle, 
                               StandardSize, ...){

          pCfT <- predict(CfT, newdata = list(tt = temp), 
                            type = 'response', se.fit = T)
          predCM <- predict(Allo, newdata = list(wt = W), 
                            interval = "prediction")
          predCMss <- predict(Allo, newdata = list(wt = StandardSize), 
                            interval = "prediction")
          prC1 <- rnorm(length(W), pCfT$fit, 
                            SAToggle["Consumption"] * pCfT$se.fit)
          prCs <- exp(rnorm(length(W), (predCMss[, "fit"]), 
                              SAToggle["Consumption"] * (((predCMss[, "fit"]) - 
                                        (predCMss[, "lwr"])) / 1.96)))
          prCM1 <- exp(rnorm(length(W), (predCM[, "fit"]), 
                              SAToggle["Consumption"] * (((predCM[, "fit"]) - 
                                        (predCM[, "lwr"])) / 1.96)))
          Cons <- (prC1 / prCs) * prCM1 * p

          return(Cons)
      }	

  # b. Respiration

    # respiration function for the bioenergetics sensitivity analysis

      Respiration <- function(W, temp, SAToggle, RespResM, ...){

          RA <- rnorm(length(W), 0.0028, 0.0003 * SAToggle["Respiration"])
          RB <- rnorm(length(W), -0.218, 0.016 * SAToggle["Respiration"])
          RQ <- rnorm(length(W), 0.076, 0.0028 * SAToggle["Respiration"])
          ACTIVITY <- rnorm(length(W), 1.649, 0.105 * SAToggle["Respiration"])
          pR1 <- (RA * W^RB) * (exp(RQ * temp)) * ACTIVITY

          # nab parameters from regression equation to fit quickly

            b1 <- coefficients(RespResM)[1]
            b2 <- coefficients(RespResM)[2]
            b3 <- coefficients(RespResM)[3]
            b4 <- coefficients(RespResM)[4]
            b5 <- coefficients(RespResM)[5]
            b6 <- coefficients(RespResM)[6]
            b7 <- coefficients(RespResM)[7]
            b8 <- coefficients(RespResM)[8]

          ResidSD <- b1 + b2 * temp + b3 * temp^2 + b4 * W + b5 * W^2 +
                     b6 * W^3 + b7 * W^4 + b8 * W^5
          Resp <- rnorm(length(W), pR1, ResidSD * SAToggle["Respiration"])

          return(Resp)
      }

  # c. Egestion

    # Egestion function for the bioenergetics sensitivity analysis

      Egestion <- function(Con, SAToggle, ...){
		
        mF <- 0.106
        sdF <- 0.0049 * as.numeric(SAToggle["Egestion"])
        Egest <- rnorm(length(Con), mF, sdF) * Con

        return(Egest)
      }

  # d. Excretion

    # Excretion function for the bioenergetics sensitivity analysis

      Excretion <- function(Con, Ege, SAToggle, ...){
		
        mU <- 0.088
        sdU <- 0.0099 * as.numeric(SAToggle["Excretion"])
        Excre <- rnorm(length(Con), mU, sdU) * (Con - Ege)

        return(Excre)
      }

  # e. Specific dynamic action

    # SDA function for the bioenergetics sensitivity analysis

      SpecificDynamicAction <- function(Con, Ege, SAToggle, ...){
			
        mS <- 0.159
        sdS <- 0.047 * as.numeric(SAToggle["SDA"])
        SDA <- rnorm(length(Con), mS, sdS) * (Con - Ege)

        return(SDA)
      }



################################################################################
#
# 3. Modeling functions
#
################################################################################

  # a. daily growth

    DailyGrowth <- function(W0, p, Te, ConAllo, ConTempF, StandSize, 
                            respresSDm, SAToggle, PreyEDtable, histPdiet, 
                            histFOdiet, FOmv, FOsev, ...){

        OxyCon <- 13560
        Cday <- Consumption(W = W0, temp = Te[1], p, Allo = ConAllo, 
                                CfT = ConTempF, StandardSize = StandSize, 
                                SAToggle)
        Rday <- Respiration(W = W0, temp = Te[1], SAToggle, 
                                RespResM = respresSDm)
        Fday <- Egestion(Con = Cday, SAToggle)
        Uday <- Excretion(Con = Cday, Ege = Fday, SAToggle)
        SDAday <- SpecificDynamicAction(Con = Cday, Ege = Fday, SAToggle)
        PredED <- rnorm(length(W0), 6488, 483 * as.numeric(SAToggle["PredED"]))

        PreyED <- matrix(NA, nrow = length(W0), ncol = nrow(PreyEDtable))
        for(i in 1:ncol(PreyED)){
          PreyED[,i] <- rnorm(length(W0), PreyEDtable[i,2], 
                        PreyEDtable[i,3] * as.numeric(SAToggle["PreyED"]))
        }

        DietProp <- matrix(0, nrow = length(W0), ncol = nrow(PreyEDtable))
        pNd <- matrix(0, nrow = length(W0), ncol = nrow(PreyEDtable))
        pWd <- matrix(0, nrow = length(W0), ncol = nrow(PreyEDtable))
        pFOd <- matrix(0, nrow = length(W0), ncol = nrow(PreyEDtable))

        for(i in 1:11){
          pNd[,i] <- rnorm(length(W0), histPdiet$percNmean[i], 
                        histPdiet$percNsd[i] * as.numeric(SAToggle["DietP"]))
          pWd[,i] <- rnorm(length(W0), histPdiet$percNmean[i], 
                        histPdiet$percNsd[i] * as.numeric(SAToggle["DietP"]))
          pFOd[,i] <- ilogit(rnorm(length(W0), FOmv[i], 
                        FOsev[i] * as.numeric(SAToggle["DietFO"])))
        }
        for(i in 12:14){
          pNd[,i] <- rnorm(length(W0), histPdiet$percNmean[i], 
                        histPdiet$percNsd[i]*as.numeric(SAToggle["DietP"]))
          pWd[,i] <- rnorm(length(W0), histPdiet$percNmean[i], 
                        histPdiet$percNsd[i]*as.numeric(SAToggle["DietP"]))
          pFOd[,i] <- rnorm(length(W0), histFOdiet$FOmean[i-11], 
                        histFOdiet$FOmean[i-11]*as.numeric(SAToggle["DietFO"]))
        }

        pNd[pNd < 0] <- 0
        pWd[pWd < 0] <- 0
        pFOd[pFOd < 0] <- 0
        RII <- pNd + pWd + pFOd
        DietProp <- RII / apply(RII, 1, sum)
        EDxDP <- PreyED * DietProp

      # Consumption in specific joules

        Cday_specificjoules <- Cday * apply(EDxDP, 1, sum)

      # Growth in grams per day, pre spawning

        Growth <- ((Cday_specificjoules - ((Fday + Uday + SDAday) * 
                   (apply(EDxDP, 1, sum)) + (Rday * OxyCon))) / PredED) * W0

      return(Growth)

    }


  # b. time series of growth

    TSGrowth <- function(WI, Tes,pv,Days, WaterYr, respresSDmodel, ConAllo, 
                         ConTempF, StandSizeTable, AgeClass, SAToggle, 
                         PreyEDtable, histPdiet, histFOdiet, 
                         FOmLUT, FOseLUT, ...){

        AC <- AgeClass[1]
        StandSize <- StandSizeTable[which(StandSizeTable[, "AC"] == AC), "SS"]
        ConTempF <- ConTempFunList[[AC]]
        Wmat <- data.frame(matrix(NA, ncol = Days, nrow = length(WI)))
        Wmat[, 1] <- WI

        for(i in 1:(Days-1)){
          W1 <- Wmat[,i]
          spd <- which(FOmLUT[, 1] == WaterYr & FOmLUT[, 2] == i)
          FOmv <- as.numeric(FOmLUT[spd, 3:13])
          FOsev <- as.numeric(FOseLUT[spd, 3:13])
          GI <- DailyGrowth(W0 = W1, Te = Tes[i], p = pv, 
                              respresSDm = respresSDmodel, ConAllo, ConTempF, 
                              StandSize, SAToggle, PreyEDtable, histPdiet, 
                              histFOdiet, FOmv, FOsev)
          W2 <- W1+GI
          Wmat[,i+1] <- W2
        }	

        return(Wmat)
    }


  # c. estimation of p value

    EstPVal <- function(W0, WF, Days, Te, Wyear, pvalues, respresSDmodel, 
                        ConAllo, ConTempF, StandSizeTable, AgeClass, 
                        SAToggle, PreyEDtable, histPdiet, histFOdiet, 
                        FOmLUT, FOseLUT, ...){

        pval <- pvalues
        TSpredFW <- matrix(NA, nrow = length(W0), ncol = length(pval))

        for(j in 1:length(pval)){
                xx <- TSGrowth(WI = W0, Tes = Te, Days, WaterYr = Wyear, 
                        pv = pval[j], respresSDmodel, ConAllo, ConTempF, 
                        StandSizeTable, AgeClass, SAToggle, PreyEDtable, 
                        histPdiet, histFOdiet, FOmLUT, FOseLUT)
                TSpredFW[, j] <- xx[, Days]
        }

        pvpreds <- rep(NA, length(WF))

        for(i in 1:length(WF)){
                x1 <- abs(WF[i] - TSpredFW[i,])
                pvpreds[i] <- pval[which(x1 == min(x1))]
        }

        return(pvpreds)
    }



  # d. daily growth with consumption

    DailyGrowthWCR <- function(W0, p, Te, ConAllo, ConTempF, StandSize, 
                                respresSDm, SAToggle, PreyEDtable, 
                                histPdiet, histFOdiet, FOmv, FOsev,...){

        OxyCon <- 13560
        Cday <- Consumption(W = W0, temp = Te[1], p, Allo = ConAllo, 
                            CfT = ConTempF, StandardSize = StandSize, SAToggle)
        Rday <- Respiration(W = W0, temp = Te[1], SAToggle, 
                            RespResM = respresSDm)
        Fday <- Egestion(Con = Cday, SAToggle)
        Uday <- Excretion(Con = Cday, Ege = Fday, SAToggle)
        SDAday <- SpecificDynamicAction(Con = Cday, Ege = Fday, SAToggle)
        PredED <- rnorm(length(W0), 6488, 483 * as.numeric(SAToggle["PredED"]))
        PreyED <- matrix(NA, nrow = length(W0), ncol = nrow(PreyEDtable))
        for(i in 1:ncol(PreyED)){
          PreyED[, i] <- rnorm(length(W0), PreyEDtable[i, 2], 
                           PreyEDtable[i, 3] * as.numeric(SAToggle["PreyED"]))
        }
        DietProp <- matrix(0, nrow = length(W0), ncol = nrow(PreyEDtable))
        pNd <- matrix(0, nrow = length(W0), ncol = nrow(PreyEDtable))
        pWd <- matrix(0, nrow = length(W0), ncol = nrow(PreyEDtable))
        pFOd <- matrix(0, nrow = length(W0), ncol = nrow(PreyEDtable))
        for(i in 1:11){
          pNd[, i] <- rnorm(length(W0), histPdiet$percNmean[i], 
                         histPdiet$percNsd[i] * as.numeric(SAToggle["DietP"]))
          pWd[, i] <- rnorm(length(W0), histPdiet$percNmean[i], 
                         histPdiet$percNsd[i] * as.numeric(SAToggle["DietP"]))
          pFOd[, i] <- ilogit(rnorm(length(W0), FOmv[i], 
                        FOsev[i] * as.numeric(SAToggle["DietFO"])))
        }
        for(i in 12:14){
          pNd[, i] <- rnorm(length(W0), histPdiet$percNmean[i], 
                         histPdiet$percNsd[i] * as.numeric(SAToggle["DietP"]))
          pWd[, i] <- rnorm(length(W0), histPdiet$percNmean[i], 
                         histPdiet$percNsd[i] * as.numeric(SAToggle["DietP"]))
          pFOd[, i] <- rnorm(length(W0), histFOdiet$FOmean[i - 11], 
                         histFOdiet$FOmean[i - 11] * 
                           as.numeric(SAToggle["DietFO"]))
        }
        pNd[pNd < 0] <- 0
        pWd[pWd < 0] <- 0
        pFOd[pFOd < 0] <- 0
        RII <- pNd + pWd + pFOd
        DietProp <- RII / apply(RII, 1, sum)
        EDxDP <- PreyED * DietProp

        # Consumption in specific joules

          Cday_specificjoules <- Cday * apply(EDxDP, 1, sum)

        # Growth in grams per day, pre spawning

          Growth <- ((Cday_specificjoules - ((Fday + Uday + SDAday) * 
                        (apply(EDxDP, 1, sum)) + (Rday * OxyCon))) / PredED) *
                         W0

        # Prey Consumed in gs

          PC <- DietProp * Cday * W0

        return(data.frame(Growth,PC))
    }


  # e. time series of growth with consumption

    TSGrowthWCR <- function(WI, Tes, pv, Days, WaterYr, respresSDmodel, 
                            ConAllo, ConTempF, StandSizeTable, AgeClass, 
                            SAToggle, PreyEDtable, histPdiet, histFOdiet, 
                            FOmLUT, FOseLUT, ...){

        StandSize <- StandSizeTable[which(StandSizeTable[,"AC"] == AgeClass),
                                    "SS"]
        ConTempF <- ConTempFunList[[AgeClass]]
        Wmat <- data.frame(matrix(NA, ncol = Days, nrow = length(WI)))
        PCmat <- data.frame(matrix(0, ncol = nrow(histPdiet), 
                                   nrow = length(WI)))
        Wmat[,1] <- WI
        for(i in 1:(Days - 1)){
          W1 <- Wmat[,i]
          spd <- which(FOmLUT[,1] == WaterYr & FOmLUT[,2] == i)
          FOmv <- as.numeric(FOmLUT[spd,3:13])
          FOsev <- as.numeric(FOseLUT[spd,3:13])
          oo <- DailyGrowthWCR(W0 = W1, Te = Tes[i], p = pv, 
                               respresSDm = respresSDmodel, ConAllo, 
                               ConTempF, StandSize, SAToggle, PreyEDtable, 
                               histPdiet, histFOdiet, FOmv, FOsev)
          W2 <- W1 + oo$Growth
          Wmat[, i+1] <- W2
          PCmat <- PCmat + oo[, 2:ncol(oo)]
        }	
		
      return(PCmat)
    }

    TSGrowthWCRts <- function(WI, Tes, pv, Days, WaterYr, respresSDmodel, 
                              ConAllo, ConTempF, StandSizeTable, AgeClass, 
                              SAToggle, PreyEDtable, histPdiet, histFOdiet, 
                              FOmLUT, FOseLUT,...){

      StandSize <- StandSizeTable[which(StandSizeTable[, "AC"] == AgeClass), 
                                    "SS"]
      ConTempF <- ConTempFunList[[AgeClass]]
      Wmat <- data.frame(matrix(NA, ncol = Days, nrow = length(WI)))
      PCmat <- data.frame(matrix(0, ncol = nrow(histPdiet), nrow = Days))
      Wmat[, 1] <- WI

      for(i in 1:(Days - 1)){
        W1 <- Wmat[, i]
        spd <- which(FOmLUT[, 1] == WaterYr & FOmLUT[, 2] == i)
        FOmv <- as.numeric(FOmLUT[spd, 3:13])
        FOsev <- as.numeric(FOseLUT[spd, 3:13])
        oo <- DailyGrowthWCR(W0 = W1, Te = Tes[i], p = pv, 
                              respresSDm = respresSDmodel, ConAllo, ConTempF, 
                              StandSize, SAToggle, PreyEDtable, histPdiet, 
                              histFOdiet, FOmv, FOsev)
        W2 <- W1 + oo$Growth
        Wmat[, i + 1] <- W2
        PCmat[(i + 1), ] <- oo[, 2:ncol(oo)]
      }	

      return(PCmat)
    }

  # f. weights through time

    CohortWts <- function(SampleLengths, SampleWYs, SampleDates, AgeC, 
                          DateWindows, LAA, LengthGrowthRates, LWreg, ConAllo, 
                          nr1, Temps, ConTempF, respresSDmodel, PreyEDtable,
                          PVAL, histPdiet, histFOdiet, FOmLUT, 
                          FOseLUT, SAToggle, ...){

      # determine the first and last date of the date window for the WY, 
      #    and then translate the date of sampling to the day of sampling 
      #    and the days after sampling

        FD <- as.character(DateWindows$FirstDate[which(DateWindows$Year 
                              %in% SampleWYs[1])])
        LD <- as.character(DateWindows$LastDate[which(DateWindows$Year 
                              %in% SampleWYs[1])])
        DayOfSampling <- as.numeric(SampleDates - as.Date(FD) + 1)
        DaysAfterSampling <- as.numeric(as.Date(LD) - SampleDates)

      # Grab the relevant temperatures

        TEMPS <- Temps$Temp[which(Temps$Year %in% SampleWYs[1])]

      # how many days

        DAYS <- DayOfSampling + DaysAfterSampling

      # determine the age class of the individual

        AgeClass <- rep(AgeC, length(SampleLengths))

      # estimate the beginning and ending lengths (with or without variation)

        GR <- LengthGrowthRates$LengthGrowthMean[AgeClass]	
        SD <- LengthGrowthRates$LengthGrowthSD[AgeClass] * 
                as.numeric(SAToggle["LengthGrowth"])
        Mini <- SampleLengths - (DayOfSampling - 1) * GR
        SDini <- (DayOfSampling-1)*SD
        Mfinal <- SampleLengths + DaysAfterSampling * GR
        SDfinal <- DaysAfterSampling * SD
        lengs <- Mini + (0:(DAYS - 1)) * GR
        wts <- exp(predict(LWreg, newdata = list(L = lengs)))

      return(wts)
    }



#
# 4. Sensitivity analysis functions
#

  # a. individual

    IndivSA <- function(SampleLength, SampleWY, SampleDate, DateWindows, LAA,
                        LengthGrowthRates, LWreg, ConAllo, nr1, Temps, ConTempF, 
                        respresSDmodel, PreyEDtable, PVAL, histPdiet, 
                        histFOdiet, FOmLUT, FOseLUT, ...){

      # determine the first and last date of the date window for the WY, 
      #  and then translate the date of sampling to the day of sampling 
      #  and the days after sampling

        FD <- as.character(DateWindows$FirstDate[which(DateWindows$Year 
                              %in% SampleWY)])
        LD <- as.character(DateWindows$LastDate[which(DateWindows$Year 
                              %in% SampleWY)])
        DayOfSampling <- as.numeric(SampleDate - as.Date(FD) + 1)
        DaysAfterSampling <- as.numeric(as.Date(LD) - SampleDate)

      # Grab the relevant temperatures

        TEMPS <- Temps$Temp[which(Temps$Year %in% SampleWY)]

      # how many days

        DAYS <- DayOfSampling + DaysAfterSampling

      # determine the age class of the individual

        AgeClass <- which(LAA$MinLength <= SampleLength & 
                          LAA$MaxLength > SampleLength)	

      # estimate the beginning and ending lengths

        GR <- LengthGrowthRates$LengthGrowthMean[AgeClass]	
        SD <- LengthGrowthRates$LengthGrowthSD[AgeClass] *
                        as.numeric(SAToggle["LengthGrowth"])

        Mini <- SampleLength - (DayOfSampling - 1) * GR
        SDini <- (DayOfSampling - 1) * SD
        Mfinal <- SampleLength + DaysAfterSampling * GR
        SDfinal <- DaysAfterSampling * SD

      # use the beginning and ending lengths to estimate beginning and 
      #   ending weights generating nr1 lengths, each generating a weight 
      #   drawn from a distribution with a mean and standard deviation

        predIL <- rnorm(nr1, Mini, SDini)
        predIWa <- predict(LWreg, newdata = list(L = predIL), 
                           interval = "prediction")
        predIWmean <- exp(predIWa[, "fit"])
        predIWsd <- ((exp(predIWa[, "fit"]) - exp(predIWa[, "lwr"])) / 1.96) * 
                      as.numeric(SAToggle["LengthWeight"])
        predIW <- rnorm(nr1, predIWmean, predIWsd)

        predFL <- rnorm(nr1, Mfinal, SDfinal)
        predFWa <- predict(LWreg, newdata = list(L = predFL), 
                           interval = "prediction")
        predFWmean <- exp(predFWa[, "fit"])
        predFWsd <- ((exp(predFWa[, "fit"]) - exp(predFWa[, "lwr"])) / 1.96) *
                      as.numeric(SAToggle["LengthWeight"])
        predFW <- rnorm(nr1, predFWmean, predFWsd)

      # estimate p value for each weight

        pp <- EstPVal(W0 = predIW, WF = predFW, Days = DAYS, Te = TEMPS, 
                      Wyear = SampleWY, pvalues = PVAL, respresSDmodel, 
                      ConAllo, ConTempF, StandSizeTable, AgeClass, SAToggle, 
                      PreyEDtable, histPdiet, histFOdiet, FOmLUT, FOseLUT)

      # estimate consumption at those p values through the time series

        consest <- TSGrowthWCR(WI = predIW, Tes = TEMPS, pv = pp, 
                              Days = DAYS, WaterYr = SampleWY, 
                              respresSDmodel, ConAllo, ConTempF, 
                              StandSizeTable, AgeClass, SAToggle, 
                              PreyEDtable, histPdiet, histFOdiet, 
                              FOmLUT, FOseLUT)

        colnames(consest)<-as.character(PreyEDtable$PreyClass)

      return(consest)	
	
  }

  # b. cohort

    CohortSA <- function(SampleLengths, SampleWYs, SampleDates, AgeC, 
                         DateWindows, LAA, LengthGrowthRates, LWreg, 
                         ConAllo, nr1, Temps, ConTempF, respresSDmodel, 
                         PreyEDtable, PVAL, histPdiet, histFOdiet, 
                         FOmLUT, FOseLUT,SAToggle, ...){

      # determine the first and last date of the date window for the WY, 
      #  and then translate the date of sampling to the day of sampling 
      #  and the days after sampling

        FD <- as.character(DateWindows$FirstDate[which(DateWindows$Year %in% 
                                SampleWYs[1])])
        LD <- as.character(DateWindows$LastDate[which(DateWindows$Year %in% 
                                SampleWYs[1])])
        DayOfSampling <- as.numeric(SampleDates - as.Date(FD) + 1)
        DaysAfterSampling <- as.numeric(as.Date(LD) - SampleDates)

      # Grab the relevant temperatures

        TEMPS <- Temps$Temp[which(Temps$Year %in% SampleWYs[1])]

      # how many days

        DAYS <- DayOfSampling + DaysAfterSampling

      # determine the age class of the individual

        AgeClass <- rep(AgeC, length(SampleLengths))

      # estimate the beginning and ending lengths (with or without variation)

        GR <- LengthGrowthRates$LengthGrowthMean[AgeClass]	
        SD <- LengthGrowthRates$LengthGrowthSD[AgeClass] * 
                        as.numeric(SAToggle["LengthGrowth"])

        Mini <- SampleLengths - (DayOfSampling - 1) * GR
        SDini <- (DayOfSampling - 1) * SD
        Mfinal <- SampleLengths + DaysAfterSampling * GR
        SDfinal <- DaysAfterSampling * SD

      # draw nr1 from the individuals, this sets it up so that the same 
      #  individual is used for the corresponding initial and final lengths

        indis <- 1:length(SampleLengths)
        IndsToUse <- sample(indis, nr1, replace = T)

      # use the beginning and ending lengths to estimate beginning 
      #  and ending weights generating nr1 lengths, each generating a 
      #  weight drawn from a distribution with a mean and standard 
      #  deviation (which was already set to 0 if needed)

        predIL <- rnorm(nr1, Mini[IndsToUse], SDini[IndsToUse])
        predIWa <- predict(LWreg, newdata = list(L = predIL), 
                           interval = "prediction")
        predIWmean <- exp(predIWa[, "fit"])
        predIWsd <- ((exp(predIWa[, "fit"]) - exp(predIWa[, "lwr"])) / 1.96) * 
                        as.numeric(SAToggle["LengthWeight"])
        predIW <- rnorm(nr1, predIWmean, predIWsd)

        predFL <- rnorm(nr1, Mfinal[IndsToUse], SDfinal[IndsToUse])
        predFWa <- predict(LWreg, newdata = list(L = predFL), 
                           interval = "prediction")
        predFWmean <- exp(predFWa[,"fit"])
        predFWsd <- ((exp(predFWa[,"fit"]) - exp(predFWa[, "lwr"])) / 1.96) *
                        as.numeric(SAToggle["LengthWeight"])
        predFW <- rnorm(nr1, predFWmean, predFWsd)

      # estimate p value for each weight

        pp <- EstPVal(W0 = predIW, WF = predFW, Days = DAYS[1], Te = TEMPS, 
                      Wyear = SampleWYs[1], pvalues = PVAL, respresSDmodel, 
                      ConAllo, ConTempF, StandSizeTable, AgeClass = AgeC, 
                      SAToggle, PreyEDtable, histPdiet, histFOdiet, 
                      FOmLUT, FOseLUT)

      # estimate consumption at those p values through the time series
 
        consest <- TSGrowthWCR(WI = predIW, Tes = TEMPS, pv = pp,Days = DAYS[1], 
                               WaterYr = SampleWYs[1], respresSDmodel, ConAllo, 
                               ConTempF, StandSizeTable, AgeClass = AgeC, 
                               SAToggle, PreyEDtable, histPdiet, histFOdiet, 
                               FOmLUT, FOseLUT)
        colnames(consest) <- as.character(PreyEDtable$PreyClass)

      return(consest)	
	
    }

  # c. no variance
	
    CohortConsTS <- function(SampleLengths, SampleWYs, SampleDates, AgeC, 
                             DateWindows, LAA, LengthGrowthRates, LWreg, 
                             ConAllo, nr1, Temps, ConTempF, respresSDmodel, 
                             PreyEDtable, PVAL, histPdiet, histFOdiet, 
                             FOmLUT, FOseLUT, SAToggle, ...){

      # determine the first and last date of the date window for the WY, 
      #   and then translate the date of sampling to the day of sampling 
      #   and the days after sampling

        FD <- as.character(DateWindows$FirstDate[which(DateWindows$Year %in% 
                           SampleWYs[1])])
        LD <- as.character(DateWindows$LastDate[which(DateWindows$Year %in% 
                           SampleWYs[1])])
        DayOfSampling <- as.numeric(SampleDates - as.Date(FD) + 1)
        DaysAfterSampling <- as.numeric(as.Date(LD) - SampleDates)

      # Grab the relevant temperatures

        TEMPS <- Temps$Temp[which(Temps$Year %in% SampleWYs[1])]

      # how many days

        DAYS <- DayOfSampling+DaysAfterSampling

      # determine the age class of the individual

        AgeClass <- rep(AgeC, length(SampleLengths))

      # estimate beginning and ending lengths (with or without variation)

        GR <- LengthGrowthRates$LengthGrowthMean[AgeClass]	
        SD <- LengthGrowthRates$LengthGrowthSD[AgeClass] * 
                      as.numeric(SAToggle["LengthGrowth"])

        Mini <- SampleLengths - (DayOfSampling - 1) * GR
        SDini <- (DayOfSampling - 1) * SD
        Mfinal <- SampleLengths + DaysAfterSampling * GR
        SDfinal <- DaysAfterSampling * SD

      # draw nr1 from the individuals, this sets it up so that the same 
      #   individual is used for the corresponding initial and final lengths

        indis <- 1:length(SampleLengths)
        IndsToUse <- sample(indis, nr1, replace = T)

      # use the beginning and ending lengths to estimate beginning and 
      #   ending weights generating nr1 lengths, each generating a 
      #   weight drawn from a distribution with a mean and standard 
      #   deviation (which was already set to 0 if needed)

        predIL <- rnorm(nr1, Mini[IndsToUse], SDini[IndsToUse])
        predIWa <- predict(LWreg, newdata = list(L = predIL), 
                           interval = "prediction")
        predIWmean <- exp(predIWa[, "fit"])
        predIWsd <- ((exp(predIWa[, "fit"]) - exp(predIWa[, "lwr"])) / 1.96) *
                           as.numeric(SAToggle["LengthWeight"])
        predIW <- rnorm(nr1, predIWmean, predIWsd)

        predFL <- rnorm(nr1, Mfinal[IndsToUse], SDfinal[IndsToUse])
        predFWa <- predict(LWreg, newdata=list(L = predFL), 
                           interval = "prediction")
        predFWmean <- exp(predFWa[, "fit"])
        predFWsd <- ((exp(predFWa[, "fit"]) - exp(predFWa[, "lwr"])) / 1.96) *
                           as.numeric(SAToggle["LengthWeight"])
        predFW <- rnorm(nr1, predFWmean, predFWsd)

      # estimate p value for each weight

        pp <- EstPVal(W0 = predIW, WF = predFW, Days = DAYS[1], Te = TEMPS, 
                      Wyear = SampleWYs[1], pvalues = PVAL, respresSDmodel, 
                      ConAllo, ConTempF, StandSizeTable, AgeClass = AgeC, 
                      SAToggle, PreyEDtable, histPdiet, histFOdiet, 
                      FOmLUT, FOseLUT)

      # estimate consumption at those p values through the time series

        consest <- TSGrowthWCRts(WI = predIW, Tes = TEMPS, pv = pp, 
                                 Days = DAYS[1], WaterYr = SampleWYs[1], 
                                 respresSDmodel, ConAllo, ConTempF, 
                                 StandSizeTable, AgeClass = AgeC, SAToggle, 
                                 PreyEDtable, histPdiet, histFOdiet, 
                                 FOmLUT, FOseLUT)
        colnames(consest) <- as.character(PreyEDtable$PreyClass)
	
      print(paste("AC ", AgeClass, " WY ", SampleWYs[1], " p = ", pp, sep = ""))
      return(consest)	

    }
