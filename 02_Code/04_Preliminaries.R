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

# Setup Dates/Seasons
dates <- data.frame("date" = seq.Date(ymd("2016-10-01"),ymd("2023-09-30"), 
                                      by = "1 day")) %>%
  mutate(month = month(date), 
         year = year(date),
         doy = yday(date),
         wy_doy = case_when(
           month > 9 & year %% 4 == 0 ~ doy - 274,
           month > 9 & year %% 4 != 0 ~ doy - 273,
           month <= 9 ~ doy + 92
         ),
         wy_year = case_when(
           month > 9 ~ year,
           month <= 9 ~ year -1
         )
         )
