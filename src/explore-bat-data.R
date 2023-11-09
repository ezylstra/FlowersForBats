# Exploring bat data
# Erin Zylstra
# ezylstra@arizona.edu
# 2023-11-08

library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)
library(ggplot2)
library(geodata)
library(terra)
library(tidyterra)

rm(list = ls())

# Load BMGRE data
gr <- read.csv("data/agfd/agfd-bmgre-obs.csv", na.strings = c(NA, ""),
                    col.names = c("site", "date", "call_count", "pa"))
gr_locs <- read.csv("data/agfd/agfd-bmgre-locations.csv",
                    col.names = c("site_name", "datum", "zone", "east", "north"))

# Load backyard data
byb <- read.csv("data/agfd/agfd-backyard-daily.csv", na.strings = c(NA, ""),
                col.names = c("loc", "loc_yr", "date", "fluid_drop", 
                              "wind", "bat_obs", "comments"))
byb_expids <- read.csv("data/agfd/agfd-backyard-expert-ids.csv",
                       col.names = c("loc", "loc_yr", "Lepto", "Choero"))
byb_locs <- read.csv("data/agfd/agfd-backyard-locations.csv",
                     na.strings = c("000000", "0000000", 
                                    "Needs", "No", "Info", "Update"),
                     col.names = c("loc", "east", "north", "eastnorth", "zone"))

# BMGRE data
  count(gr, site)
  gr_locs$site <- c("ThanksgivingDay", "ThompsonTank", "Mohawk2", "Tank499", 
                    "AFAF", "BenderSpring")
  gr <- gr %>%
    mutate(obsdate = mdy(date),
           yr = year(obsdate),
           yday = yday(obsdate)) %>%
    select(-date)

  # Summarize data by year
  count(gr, pa, call_count)
  gr_yr <- gr %>%
    group_by(yr) %>%
    summarize(n_nights = length(pa),
              night_first = min(yday),
              night_last = max(yday),
              n_dets = sum(pa),
              prop_dets = round(n_dets / length(obsdate), 2),
              n_dets_100192 = sum(pa[yday %in% 100:192]),
              det_first = ifelse(n_dets > 0, min(yday[pa == 1]), NA),
              det_last = ifelse(n_dets > 0, max(yday[pa == 1]), NA),
              mn_ct = round(mean(call_count[pa == 1]), 2),
              max_ct = max(call_count)) %>%
    data.frame()
  gr_yr
    # Data/Detections in 2013-2015 are really sparse, probably better to skip.
    # Interesting that there are detections almost year round. False detections?
    # Notes in the report that LLNBs have potential to be present on BMGRE from
    # April 10th (yday 100) to July 11th (192)
  
  # Summarize data by year and site
  gr_siteyr <- gr %>%
    filter(yr > 2015) %>%
    group_by(site, yr) %>%
    summarize(n_nights = length(pa),
              night_first = min(yday),
              night_last = max(yday),
              n_dets = sum(pa),
              prop_dets = round(n_dets / length(obsdate), 2),
              n_dets_100192 = sum(pa[yday %in% 100:192]),
              det_first = ifelse(n_dets > 0, min(yday[pa == 1]), NA),
              det_last = ifelse(n_dets > 0, max(yday[pa == 1]), NA),
              mn_ct = round(mean(call_count[pa == 1]), 2),
              max_ct = max(call_count),
              .groups = "keep") %>%
    data.frame()
  gr_siteyr
  
  # Visualize this...
  # Horizontal bar for each site, year, with vertical bands (white = not
  # operating; gray = no detections; light color = call = 1; dark = call > 1)
  
  # Map with locations of 6 sites
  
# Backyard data
  # Clean up data and remove any observations if there's no location info
  byb <- byb %>%
    mutate(obsdate = dmy(date),
           yr = year(obsdate),
           yday = yday(obsdate)) %>%
    select(-c(date, loc_yr)) %>%
    left_join(select(byb_locs, loc, east, north), by = "loc") %>%
    filter(!is.na(east))
  
  count(byb, fluid_drop, bat_obs)
    # These data are very messy....
    # not sure what bat_obs = -1 indicates (Often associated with fluid drop)
  
  # Attach expert IDs (for year, not date)
  byb_expids <- byb_expids %>%
    mutate(yr = paste0("20", str_sub(loc_yr, -2, -1))) %>%
    filter(str_sub(yr, -1, -1) != "-") %>%
    mutate(yr = as.numeric(yr)) %>%
    distinct(loc, yr, Lepto)
  # Couple loc/yrs with both TRUE AND FALSE. Will remove TRUE
  byb_expids[which(duplicated(byb_expids[,1:2])),]
  filter(byb_expids, loc == 525, yr == 2018)
  filter(byb_expids, loc == 1005, yr == 2017)
  byb_expids <- byb_expids %>%
    filter(!(loc == 525 & yr == 2018 & Lepto == TRUE)) %>%
    filter(!(loc == 1005 & yr == 2017 & Lepto == TRUE))

  byb <- byb %>%
    left_join(byb_expids, by = c("loc", "yr"))
  count(distinct(byb, loc, yr, Lepto), Lepto)
    # 15 FALSE, 139 TRUE, 1251 NAs.


