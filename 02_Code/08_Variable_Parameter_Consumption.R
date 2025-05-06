################################################################################
#
# 4. Variable Parameter Calculations
#
################################################################################

# Total Basic with Direct Accumulation, p < 1

# for each predator, calculate the amount of food they would have eaten
#   over each day, accounting for temperature, but ignoring
#   predator growth
#   repeat NPERM amount of times

# because gam predict only gives the error of the fit, 
#    we expand the error to include residuals/generate prediction interval
library(tidyverse)
library(lubridate)
library(foreach)
library(parallel)
library(doSNOW)
library(progress)
set.seed(123)

# Bring in the Data
LWdataA <- readRDS(file.path("01_Data","Output","LWdataA.rds"))
dietfVCHN <- readRDS(file.path("01_Data","Input","diet_fraction_variable_analysis.rds"))
TeT <- readRDS(file.path("01_Data","Input","temps.rds")) %>%
  group_by(WaterYear, Date) %>%
  summarize(Temp = mean(Value, na.rm = TRUE))

#Bring in the consumption functions and parameter lookup tables
g1 <- readRDS(file.path("01_Data","Input","g1.rds"))
g2 <- readRDS(file.path("01_Data","Input","g2.rds"))
g3 <- readRDS(file.path("01_Data","Input","g3.rds"))
gs <- list(g1, g2, g3, g3, g3, g3, g3, g3, g3, g3)

pred_cmax_lu <- readRDS(file.path("01_Data","Input","pred_cmax_lu.rds"))
r_sq_lu <- readRDS(file.path("01_Data","Input","r_sq_lu.rds"))
ref_cmax_lu <- readRDS(file.path("01_Data","Input","ref_cmax_lu.rds"))
con_allo_lu <- readRDS(file.path("01_Data","Input","ConAlloLU.rds"))
pvals <- c(rep(0.69, 2), rep(0.73, 1), rep(0.68, 1), rep(0.62, 6))

NPERM <- 100

LWdataA <- LWdataA %>%
  mutate(cap_wy = esaRmisc::water_year(date)) %>%
  rename("cap_date" = date)

# Setup the parallel processing
# Cores
n_cores <- detectCores()
cl <- makeCluster(n_cores - 2)
doSNOW::registerDoSNOW(cl)

TCHNconsmatVA <-list()
for(y in unique(LWdataA$cap_wy)){
  LWdataA_y <- LWdataA[LWdataA[["cap_wy"]] == y,]
  TeT_y <- TeT[TeT[["WaterYear"]] %in% c(y,y+1),]
  
  print(y)
  
  # Progress Bar
  iter = nrow(LWdataA_y)
  pb <- progress_bar$new(
    format = "[:bar] :elapsed | eta: :eta",
    total = iter, 
    width = 60)
  
  # allowing progress bar to be used in foreach -----------------------------
  progress <- function(n){
    pb$tick()
  } 
  
  opts <- list(progress = progress)
  
  fish_list <- list()
  foreach(f = 1:nrow(LWdataA_y), .options.snow = opts) %dopar% {
    
    #Subset the temperature data for 70 weeks post-capture
    TE <- TeT_y[TeT_y$Date <= (ymd(LWdataA_y[f,]$cap_date) + 490),]$Temp
    
    LWdata_f <- LWdataA_y[f,] |>
      dplyr::mutate(WW = list(round(rnorm(NPERM,  fit, 
                                   (fit - lwr) / 1.96),2))) |>
      dplyr::bind_cols(TeT_y) |>
      tidyr::unnest(WW) |>
      dplyr::mutate(WW = round(WW,3)) |>
      dplyr::rename("TE" = Temp) |>
      dplyr::left_join(pred_cmax_lu) |>
      dplyr::mutate(cmax_sds = cmax_se.fit/r_sq_lu[age][[1]]) |>
      dplyr::left_join(ref_cmax_lu) |>
      dplyr::left_join(con_allo_lu) |>
      dplyr::mutate(TCHN_cons_plt1 = rnorm(NPERM, cmax_mean, cmax_sds) / 
                      rnorm(NPERM, ref_cmax_mean,ref_cmax_sds) * 
                      rnorm(NPERM, cmaxi_mean, cmaxi_sds) * 
                      WW * pvals[age] *
                      sample(dietfVCHN, NPERM, replace = T),
                    TCHN_cons_p1 = rnorm(NPERM, cmax_mean, cmax_sds) / 
                      rnorm(NPERM, ref_cmax_mean,ref_cmax_sds) * 
                      rnorm(NPERM, cmaxi_mean, cmaxi_sds) * 
                      WW * 1 * #pval == 1
                      sample(dietfVCHN, NPERM, replace = T))
    
    LWdata_f <- LWdata_f |>
      #If the predicted value is below 0, fix to 0
      dplyr::mutate(TCHN_cons_plt1 = ifelse(TCHN_cons_plt1 < 0, 0, TCHN_cons_plt1),
                    TCHN_cons_p1 = ifelse(TCHN_cons_p1 < 0, 0, TCHN_cons_p1),
                    # if the date is <= date of capture, set consumption to 0 as well
                    TCHN_cons_plt1 = ifelse(Date <= cap_date, NA, TCHN_cons_plt1),
                    TCHN_cons_p1 = ifelse(Date <= cap_date, NA, TCHN_cons_p1)) %>%
      group_by(X, species, survey, gear, cap_wy, cap_date, fork_length_mm, weight_g) %>%
      summarize(plt1_mean = mean(TCHN_cons_plt1, na.rm = TRUE),
                plt1_lwr = quantile(TCHN_cons_plt1, 0.05, na.rm = TRUE),
                plt1_upr = quantile(TCHN_cons_plt1, 0.05, na.rm = TRUE),
                p1_mean = mean(TCHN_cons_p1, na.rm = TRUE),
                p1_lwr = quantile(TCHN_cons_p1, 0.05, na.rm = TRUE),
                p1_upr = quantile(TCHN_cons_p1, 0.05, na.rm = TRUE))
    
    
    fish_list[[f]] <- LWdata_f
  }
  TCHNconsmatVA[[y]] <- fish_list 
}
stopCluster(cl = cl)

saveRDS(TCHNconsmatVA,file.path("01_Data","Output","Tables","TCHNconsmatVA.rds"))
