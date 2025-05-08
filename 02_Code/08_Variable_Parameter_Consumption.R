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

NPERM <- 10

# Bring in the Data
LWdataA <- readRDS(file.path("01_Data","Output","LWdataA.rds"))

dietfVCHN <- readRDS(file.path("01_Data","Input","diet_fraction_variable_analysis.rds"))
TeT <- readRDS(file.path("01_Data","Input","temps.rds")) |>
  group_by(WaterYear, Date) |>
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

LWdataA <- LWdataA |>
  ## mutate(cap_wy = cfs.misc::water_year(date)) |>
  mutate(cap_wy = esaRmisc::water_year(date)) |>
  rename("cap_date" = date) |>
  mutate(cap_date_ymd = lubridate::ymd(cap_date)) |> # Calculate once on vector
  mutate(fit_sd =  (fit - lwr) / 1.96) # Calculate vectorized for later

# Having in list adds unnecessary overhead
r_sq_lu <- data.frame(age = c(1:10),
                      r_sq = unlist(r_sq_lu))

# Setup the parallel processing
# Cores
n_cores <- detectCores()
cl <- makeSOCKcluster(n_cores - 4)
doSNOW::registerDoSNOW(cl)
snow::clusterExport(cl = cl, list = c("con_allo_lu", "NPERM","pvals", 
                                      "dietfVCHN"))

# Do these ops once rather than 1x per fish
TeT_tmp <- TeT |>
  dplyr::rename("TE" = Temp) |>
  dplyr::mutate(TE = round(TE,2)) |>
  dplyr::left_join(pred_cmax_lu) |>
  dplyr::left_join(ref_cmax_lu) |>
  dplyr::left_join(r_sq_lu)

# lapply() is more efficient than growing a list
TCHNconsmatVA = lapply(unique(LWdataA$cap_wy), function(y) {
  LWdataA_y <- LWdataA[LWdataA[["cap_wy"]] == y,]
  TeT_y <- TeT_tmp[TeT_tmp[["WaterYear"]] %in% c(y,y+1),]
  
  print(y)

  # Progress bar adds overhead 
  if(TRUE){
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
  }
  
  fish_list <- foreach(f = 1:nrow(LWdataA_y), .options.snow = opts) %dopar% {
  #fish_list <- parallel::mclapply(seq(nrow(LWdataA_y)), function(f) {
    
    #Subset the temperature data for 70 weeks post-capture
    TE <- TeT_y[TeT_y$Date <= (LWdataA_y[f,]$cap_date_ymd + 490),]$TE 

    # Lots of time here - rearranging gives small gain
    LWdata_f <- LWdataA_y[f,] |>
      dplyr::left_join(TeT_y) |>
      dplyr::mutate(cmax_sds = cmax_se.fit/r_sq) |>
      dplyr::mutate(WW = list(round(rnorm(NPERM,  fit, fit_sd),2))) |> #Round to 2 to match
      tidyr::unnest(WW) |>
      dplyr::left_join(con_allo_lu)
    
    n <- nrow(LWdata_f)
    
    TCHN_cons_plt1 <- rnorm(NPERM * n, rep(LWdata_f$cmax_mean, NPERM), rep(LWdata_f$cmax_sds, NPERM)) / 
      rnorm(NPERM * n, LWdata_f$ref_cmax_mean[1], LWdata_f$ref_cmax_sds[1]) * # constant - only need one
      rnorm(NPERM * n, rep(LWdata_f$cmaxi_mean, NPERM), rep(LWdata_f$cmaxi_sds, NPERM)) * 
      rep(LWdata_f$WW, NPERM) * pvals[rep(LWdata_f$age, NPERM)] * #pva; based on age
      sample(dietfVCHN, NPERM * NPERM, replace = TRUE)
    TCHN_cons_p1 <- rnorm(NPERM * n, rep(LWdata_f$cmax_mean, NPERM), rep(LWdata_f$cmax_sds, NPERM)) / 
      rnorm(NPERM * n, LWdata_f$ref_cmax_mean[1],LWdata_f$ref_cmax_sds[1]) * 
      rnorm(NPERM * n, rep(LWdata_f$cmaxi_mean, NPERM), rep(LWdata_f$cmaxi_sds, NPERM)) * 
      rep(LWdata_f$WW, NPERM) * #1 * #pval == 1
      sample(dietfVCHN, NPERM * NPERM, replace = TRUE)
    
    LWdata_f <- dplyr::bind_cols(LWdata_f[rep(seq(nrow(LWdata_f)), NPERM),],
                                "TCHN_cons_p1" = TCHN_cons_p1,
                                "TCHN_cons_plt1" = TCHN_cons_plt1)

    ## ifelse() is inefficient
    #If the predicted value is below 0, fix to 0
    LWdata_f$TCHN_cons_plt1[LWdata_f$TCHN_cons_plt1 < 0] <- 0
    LWdata_f$TCHN_cons_p1[LWdata_f$TCHN_cons_p1 < 0] <- 0

    # if the date is <= date of capture, set consumption to 0 as well
    # Since we remove the NA values later in the calc, we can do this as a one-op
    LWdata_f <- LWdata_f |> dplyr::filter(!(Date <= cap_date_ymd) & !is.na(TCHN_cons_p1))

    # Lots of time here
    # Continue
    LWdata_f <- LWdata_f |>
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
  }
#  }, mc.cores = detectCores() - 4L )
  
  bind_rows(fish_list)
})

stopCluster(cl = cl)
 
TCHNconsmatVA <- bind_rows(TCHNconsmatVA)

saveRDS(TCHNconsmatVA,file.path("01_Data","Output","Tables","TCHNconsmatVA.rds"))
