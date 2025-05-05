# for predator diet: historic qpcr data, historic % data, historic FO data

qpcr <- read.csv(file.path("01_Data","Input","QPCR.csv"))
histPdiet <- read.csv(file.path("01_Data","Input","HistoricDiet.csv"))
histFOdiet <- read.csv(file.path("01_Data","Input","HistoricFODiet.csv"))

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
  FObtlw[j] <- boot::inv.logit(FOm[j] - 1.96 * FOse[j])
  FObtup[j] <- boot::inv.logit(FOm[j] + 1.96 * FOse[j])
  FObtm[j] <- boot::inv.logit(FOm[j])
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
  pFOdi <- c(boot::inv.logit(rnorm(11, as.numeric(FOm[2:12]), 
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

saveRDS(dietf, file.path("01_Data","Input","diet_fraction.rds"))
saveRDS(dietfVCHN, file.path("01_Data","Input","diet_fraction_variable_analysis.rds"))
