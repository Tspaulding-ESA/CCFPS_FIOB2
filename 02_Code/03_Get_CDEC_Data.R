library(cder)
library(dplyr)
library(lubridate)
library(esaRmisc)

clc_temp_raw = cdec_query("CLC", sensors = 146, durations = "H",
                          start.date = "2014-09-30", end.date = "2023-10-01")

clc_temp = clc_temp_raw %>% 
  mutate(WaterYear = water_year(DateTime),
         Date = date(DateTime),
         DateTime_UTC = with_tz(DateTime, "UTC"),
         Time = substring(DateTime, 12, 19),
         Time_UTC = substring(DateTime_UTC, 12, 19),
         across(c(Time, Time_UTC), ~ ifelse(.x == "", "00:00:00", .x))) %>% 
  filter(WaterYear > 2014 & WaterYear < 2024) %>% 
  select(WaterYear, Date, DateTime, DateTime_UTC, Time, Time_UTC, Value, SensorUnits)
