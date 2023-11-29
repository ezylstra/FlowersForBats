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
  
  gr_locsv <- vect(gr_locs, geom = c("east", "north"), 
                   crs = "+proj=utm +zone=12 +datum=WGS84 +units=m")

  # Super simple plots...
  # Download shapefiles with administrate boundaries if not yet done
  # us <- gadm("USA", path = "data/", level = 1)
  # mx <- gadm("MEX", path = "data/", level = 1)
  us <- readRDS("data/gadm/gadm41_USA_1_pk.rds")
  mx <- readRDS("data/gadm/gadm41_MEX_1_pk.rds")
  
  states <- c("Arizona", "California", "New Mexico", "Texas", "Nevada", "Utah", 
              "Colorado", "Kansas", 'Oklahoma')
  states_mx <- c("Sonora", "Chihuahua", "Coahuila", "Baja California")
  us.states <- us[us$NAME_1 %in% states,]
  mx.states <- mx[mx$NAME_1 %in% states_mx,]
  
  dem_files <- list.files("data/dem", full.names = TRUE)
  dem_list <- NULL
  for (i in 1:length(dem_files)) {
    dem_list[[i]] <- rast(dem_files[i])
  }
  dem_coll <- sprc(dem_list)
  dem <- merge(dem_coll)

  # Put everything in the same crs
  gr_locsll <- project(gr_locsv, crs(dem))
  us.states <- project(us.states, crs(dem))
  mx.states <- project(mx.states, crs(dem))
  
  # Plot
  ggplot() +
    geom_spatvector(data = us.states, fill = NA) +    
    # geom_spatraster(data = dem) +
    ylim(31, 34) +
    xlim(-114.5, -109) +
    geom_spatvector(data = gr_locsll, col = "red") +
    theme_bw()
  
  
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
  
  # Extract just data for sites where Lepto == TRUE (in one or more years)
  # But excluding years if Lepto == FALSE
  lsites <- sort(unique(byb_expids$loc[which(byb_expids$Lepto == TRUE)]))
  bybl <- byb %>%
    filter(loc %in% lsites) %>%
    filter(is.na(Lepto) | Lepto == TRUE)
  
  # Plot locations
  byb_locsv <- vect(byb_locs, geom = c("east", "north"), 
                    crs = "+proj=utm +zone=12 +datum=WGS84 +units=m")
  byb_locsll <- project(byb_locsv, crs(dem))
  
  bybl_locs <- filter(byb_locs, loc %in% bybl$loc)
  bybl_locsv <- vect(bybl_locs, geom = c("east", "north"), 
                     crs = "+proj=utm +zone=12 +datum=WGS84 +units=m")
  bybl_locsll <- project(bybl_locsv, crs(dem))
  
  az <- us[us$NAME_1 == "Arizona",]
  
  feeders <- ggplot() +
    geom_spatraster(data = dem) +
    scale_fill_whitebox_c(palette = "soft") +
    geom_spatvector(data = az, fill = NA) +
    ylim(31.3, 32.8) +
    xlim(-111.7, -109) +
    geom_spatvector(data = byb_locsll, col = "black") +
    geom_spatvector(data = bybl_locsll, col = "red") +
    theme_bw() +
    labs(fill = "Elevation (m)", x = "", y = "") +
    theme(legend.position = "bottom",
          axis.title = element_blank(),
          legend.box.spacing = unit(0.1, "pt"))
  feeders
  # ggsave("output/maps/feeders.png",
  #        feeders, device = "png", width = 6.5, height = 5.5,
  #        units = "in", dpi = 300)
  # # Lots in Tucson, but some also in Sierra Vista area. Wondering if the 
  # distribution of sites where Leptos were confirmed is due to ease of access...

  # Simplify fluid variable
  bybl <- bybl %>%
    mutate(fluid = ifelse(fluid_drop %in% c("Drained", "Much lower", "MuchLower"), 2, 
                          ifelse(fluid_drop %in% c("A little lower", "LittleLower"),
                                 1, 0)))
  count(bybl, fluid_drop, fluid)
  
  # Day of the year when changes in fluid levels
  bybl_fl <- bybl %>%
    group_by(fluid) %>%
    summarize(n = length(fluid),
              n_bat = sum(bat_obs == "Yes" & !is.na(bat_obs)),
              n_nobat = sum(bat_obs %in% c("No", "0")),
              doy_lcl = round(quantile(yday, prob = 0.05)),
              doy_mn = round(mean(yday)),
              doy_ucl = round(quantile(yday, prob = 0.95))) %>%
    data.frame()
  bybl_fl

  # 95% of observations with significant fluid drops occurred between:
  parse_date_time(x = bybl_fl$doy_lcl[bybl_fl$fluid == 2], orders = "j")
  parse_date_time(x = bybl_fl$doy_ucl[bybl_fl$fluid == 2], orders = "j")
  # August 14th to October 21st
  
  par(mfrow = c(3, 1))
  hist(bybl$yday[bybl$fluid == 0], breaks = 50, xlim = c(1, 365), main = "fluid 0")
  hist(bybl$yday[bybl$fluid == 1], breaks = 50, xlim = c(1, 365), main = "fluid 1") 
  hist(bybl$yday[bybl$fluid == 2], breaks = 50, xlim = c(1, 365), main = "fluid 2")
  
  # Could see how things change if we look at all sites (not worrying about
  # expert species ID) or at only those sites below 32 or 31.7 deg lat...
