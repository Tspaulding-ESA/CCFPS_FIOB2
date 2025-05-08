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

n_fish = 1000
# Bring in the Data
LWdataA <- readRDS(file.path("01_Data","Output","LWdataA.rds"))

## Subset for testing purposes
LWdataA = LWdataA[sample(seq(nrow(LWdataA)), n_fish),]

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

NPERM <- 10

LWdataA <- LWdataA %>%
  mutate(cap_wy = cfs.misc::water_year(date)) %>%
  rename("cap_date" = date) |>
  mutate(cap_date_ymd = lubridate::ymd(cap_date)) |> # Calculate once on vector
  mutate(fit_sd =  (fit - lwr) / 1.96) # Calculate vectorized for later

# Setup the parallel processing
# Cores
n_cores <- detectCores()
## cl <- makeCluster(n_cores - 4)
## doSNOW::registerDoSNOW(cl)

  # Do these ops once rather than 1x per fish
  TeT_tmp <- TeT |>
    dplyr::rename("TE" = Temp) |>
    dplyr::mutate(TE = round(TE,2)) |>
    dplyr::left_join(pred_cmax_lu) |>
    dplyr::left_join(ref_cmax_lu)

TCHNconsmatVA <-list()
Rprof("Add cols")

for(y in unique(LWdataA$cap_wy)){
  LWdataA_y <- LWdataA[LWdataA[["cap_wy"]] == y,]
  TeT_y <- TeT_tmp[TeT_tmp[["WaterYear"]] %in% c(y,y+1),]
  
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
  
  ## fish_list <- foreach(f = 1:nrow(LWdataA_y), .options.snow = opts) %dopar% {
  fish_list <- parallel::mclapply(seq(nrow(LWdataA_y)), function(f) {
    
    #Subset the temperature data for 70 weeks post-capture
    TE <- TeT_y[TeT_y$Date <= (LWdataA_y[f,]$cap_date_ymd + 490),]$TE 

    
    LWdata_f <- LWdataA_y[f,] |>
      dplyr::mutate(WW = list(round(rnorm(NPERM,  fit, fit_sd),2))) |>
      dplyr::left_join(TeT_y) |>
      tidyr::unnest(WW) |>
      dplyr::mutate(WW = round(WW,3)) |>
      dplyr::mutate(cmax_sds = cmax_se.fit/r_sq_lu[age][[1]]) |>
      dplyr::left_join(con_allo_lu)
    
    # Calculating the variation in consumption requires calculation by row
    tmp_list <- lapply(seq(nrow(LWdata_f)), function(i) {
      tmp <- LWdata_f[i,]
      TCHN_cons_plt1 = rnorm(NPERM, tmp$cmax_mean, tmp$cmax_sds) / 
        rnorm(NPERM, tmp$ref_cmax_mean,tmp$ref_cmax_sds) * 
        rnorm(NPERM, tmp$cmaxi_mean, tmp$cmaxi_sds) * 
        tmp$WW * pvals[tmp$age] * #pva; based on age
        sample(dietfVCHN, NPERM, replace = TRUE)
      TCHN_cons_p1 = rnorm(NPERM, tmp$cmax_mean, tmp$cmax_sds) / 
        rnorm(NPERM, tmp$ref_cmax_mean,tmp$ref_cmax_sds) * 
        rnorm(NPERM, tmp$cmaxi_mean, tmp$cmaxi_sds) * 
        tmp$WW * 1 * #pval == 1
        sample(dietfVCHN, NPERM, replace = TRUE)
      dplyr::bind_cols(tmp, "TCHN_cons_p1" = TCHN_cons_p1, "TCHN_cons_plt1" = TCHN_cons_plt1)
      
    })
    # rejoin everything back together
    LWdata_f <- dplyr::bind_rows(tmp_list)
    #rm(tmp_list)
    
    # Continue
    LWdata_f <- LWdata_f |>
      #If the predicted value is below 0, fix to 0
      dplyr::mutate(TCHN_cons_plt1 = ifelse(TCHN_cons_plt1 < 0, 0, TCHN_cons_plt1),
                    TCHN_cons_p1 = ifelse(TCHN_cons_p1 < 0, 0, TCHN_cons_p1))|>
      # if the date is <= date of capture, set consumption to 0 as well
      dplyr::mutate(TCHN_cons_plt1 = ifelse(Date <= cap_date_ymd, NA, TCHN_cons_plt1),
                    TCHN_cons_p1 = ifelse(Date <= cap_date_ymd, NA, TCHN_cons_p1)) |>
      dplyr::filter(!is.na(TCHN_cons_p1)) |>
      dplyr::ungroup() |>
      #Summarise by finding the mean, upper, and lower consumption values for each date for each fish
      dplyr::group_by(X, species, survey, gear, cap_wy, cap_date, fork_length_mm, age, Date) |>
      dplyr::summarize(plt1_mean = mean(TCHN_cons_plt1, na.rm = TRUE),
                       plt1_lwr = quantile(TCHN_cons_plt1, 0.025, na.rm = TRUE),
                       plt1_upr = quantile(TCHN_cons_plt1, 0.975, na.rm = TRUE),
                       p1_mean = mean(TCHN_cons_p1, na.rm = TRUE),
                       p1_lwr = quantile(TCHN_cons_p1, 0.025, na.rm = TRUE),
                       p1_upr = quantile(TCHN_cons_p1, 0.975, na.rm = TRUE))
    
    LWdata_f
  }, mc.cores = detectCores() - 4L )
  
  TCHNconsmatVA[[as.character(y)]] <- bind_rows(fish_list)
}

Rprof(NULL)
## stopCluster(cl = cl)

TCHNconsmatVA <- bind_rows(TCHNconsmatVA)

saveRDS(TCHNconsmatVA,file.path("01_Data","Output","Tables","TCHNconsmatVA.rds"))
