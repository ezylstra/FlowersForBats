# Flowering phenology, based on NPN data
# Erin Zylstra
# ezylstra@arizona.edu
# 2023-11-08

library(dplyr)
library(lubridate)
library(tidyr)

rm(list = ls())

# Load data (created in explore-npn-observations.R)
dat <- read.csv("data/FlowersforBats_CleanedObs.csv")
# Each row is a "complete" observation of an individual plant/patch on a given
# date, including phenophases (1/0/NA) and intensities (ordinal value/NA). 
# There is a maximum of one observation of each individual per day.
dat <- dat %>% mutate(obsdate = ymd(obsdate))

# Want to subset data so we're just looking at plants that are in LLNB range.
# Based on Fig 4 in the LLNB SSA, the area of interest extends to about 32.7 deg
# latitude in AZ and NM, from -108 deg longitude and west. There was only one  
# individual plant/patch observed in Sonora and that was near Guaymas, so we'll 
# exclude that for now.

# Look at distribution of latitudes in AZ
hist(unique(dat$lat[dat$state == "AZ"]), breaks = 50)

# Restrict geographic area of interest
dat <- dat %>% 
  filter(state == "AZ", lat < 32.7) %>%
  select(-state)

# Summarize information for each individual plant/patch (n = 649)
indiv <- dat %>%
  group_by(site_id, site_name, ind_id, plant_nickname, spp, species_id, patch,
           lat, lon, elev) %>%
  summarize(nyrs = length(unique(yr)),
            yr_first = min(yr),
            yr_last = max(yr),
            nobs = length(yr),
            .groups = "keep") %>%
  data.frame()

count(indiv, spp)
count(indiv, spp, patch)
count(indiv, nyrs)
  # Most individuals only observed in 1 or 2 years (n = 448, 100), but up to 11

# Summarize patch, location, elevation data for each plant species
plants <- indiv %>%
  group_by(spp, species_id) %>%
  summarize(ninds = length(ind_id),
            npatch = sum(!is.na(patch)),
            lon_min = min(lon),
            lon_max = max(lon),
            lat_min = min(lat),
            lat_max = max(lat),
            elev_min = min(elev),
            elev_max = max(elev),
            elev_mn = round(mean(elev)),
            .groups = "keep") %>%
  mutate(prop_patch = round(npatch / ninds, 2), .after = npatch) %>%
  data.frame()
select(plants, -species_id)
# Low % of individuals described as patches (2% for saguaro, 5-33% for agaves)
# A. palmeri and A. chrysantha (goldenflower) with higher mean elevations than 
  # other species (1500, 1361 vs 819-1053 m)

# To look at start/end for each phenophase, add (within yr) observation numbers
dat <- dat %>%
  arrange(indobs) %>%
  mutate(indyr = paste0(ind_id, "_", yr),
         obsnum_yr = sequence(rle(as.character(indyr))$lengths))
dat$days_since <- 
  c(NA, as.numeric(dat$obsdate[2:nrow(dat)] - dat$obsdate[1:(nrow(dat) - 1)]))
dat$days_since[dat$obsnum_yr == 1] <- NA

# Summarizing observations of each individual each year?
ind_yr <- dat %>%
  group_by(ind_id, yr, spp) %>%
  summarize(obs_first = min(obsdate),
            obs_last = max(obsdate),
            days_since_mn = round(mean(days_since, na.rm = TRUE)),
            nobs = max(obsnum_yr),
            nflowers = sum(flowers, na.rm = TRUE),
            nopen = sum(flowers_open, na.rm = TRUE),
            nfruit = sum(fruit, na.rm = TRUE),
            nripe = sum(fruit_ripe, na.rm = TRUE),
            flowers_first = ifelse(sum(flowers, na.rm = TRUE) > 0,
                                   yday(min(obsdate[which(flowers == 1)])), NA),
            open_first = ifelse(sum(flowers_open, na.rm = TRUE) > 0,
                                yday(min(obsdate[which(flowers_open == 1)])), NA),
            .groups = "keep") %>%
  data.frame()

# Number of observations per year, across species
hist(ind_yr$nobs, breaks = 50)

# Number of observations per year by species
ind_yr %>%
  group_by(spp) %>%
  summarize(n_indiv = length(unique(ind_id)),
            n_indivyrs = length(spp),
            n_oneobsperyr = sum(nobs == 1),
            n_multobsperyr = sum(nobs > 1),
            n_5obsperyr = sum(nobs > 4)) %>%
  mutate(propmult = round(n_multobsperyr / n_indiv, 2),
         prop5 = round(n_5obsperyr / n_indiv, 2)) %>%
  data.frame()

# First day of year with flowers or buds, open flowers
par(mfrow = c(2, 1))
hist(ind_yr$flowers_first, breaks = 50, xlim = c(0, 365))
hist(ind_yr$open_first, breaks = 50, xlim = c(0, 365))

# First day of year with flowers, by species (with 20 or more individuals)
par(mfrow = c(3, 2), mar = c(2,4,4,1))
hist(ind_yr$flowers_first[ind_yr$spp == "saguaro"], breaks = 25, xlim = c(0, 365),
     main = "Saguaro", xlab = "")
hist(ind_yr$flowers_first[ind_yr$spp == "American century plant"], breaks = 25, xlim = c(0, 365),
     main = "American century plant", xlab = "")
hist(ind_yr$flowers_first[ind_yr$spp == "Palmer's century plant"], breaks = 25, xlim = c(0, 365),
     main = "Palmer's century plant", xlab = "")
hist(ind_yr$flowers_first[ind_yr$spp == "goldenflower century plant"], breaks = 25, xlim = c(0, 365),
     main = "goldenflower century plant", xlab = "")
hist(ind_yr$flowers_first[ind_yr$spp == "Parry's agave"], breaks = 25, xlim = c(0, 365),
     main = "Parry's agave", xlab = "")

# First day of year with open flowers, by species (with 20 or more individuals)
par(mfrow = c(3, 2), mar = c(2,4,4,1))
hist(ind_yr$open_first[ind_yr$spp == "saguaro"], breaks = 25, xlim = c(0, 365),
     main = "Saguaro", xlab = "")
hist(ind_yr$open_first[ind_yr$spp == "American century plant"], breaks = 25, xlim = c(0, 365),
     main = "American century plant", xlab = "")
hist(ind_yr$open_first[ind_yr$spp == "Palmer's century plant"], breaks = 25, xlim = c(0, 365),
     main = "Palmer's century plant", xlab = "")
hist(ind_yr$open_first[ind_yr$spp == "goldenflower century plant"], breaks = 25, xlim = c(0, 365),
     main = "goldenflower century plant", xlab = "")
hist(ind_yr$open_first[ind_yr$spp == "Parry's agave"], breaks = 25, xlim = c(0, 365),
     main = "Parry's agave", xlab = "")

# Plotting all observations dates when 50-74% of flowers are open, by species
par(mfrow = c(3, 2), mar = c(2,4,4,1))
hist(dat$doy[which(dat$spp == "saguaro" & dat$i_flowers_open == "4: 50-74%")], 
     breaks = 25, xlim = c(0, 365), main = "Saguaro", xlab = "")
hist(dat$doy[which(dat$spp == "American century plant" & dat$i_flowers_open == "4: 50-74%")], 
     breaks = 25, xlim = c(0, 365), main = "American century plant", xlab = "")
hist(dat$doy[which(dat$spp == "Palmer's century plant" & dat$i_flowers_open == "4: 50-74%")], 
     breaks = 25, xlim = c(0, 365), main = "Palmer's century plant", xlab = "")
hist(dat$doy[which(dat$spp == "goldenflower century plant" & dat$i_flowers_open == "4: 50-74%")], 
     breaks = 25, xlim = c(0, 365), main = "goldenflower century plant", xlab = "")
hist(dat$doy[which(dat$spp == "Parry's agave" & dat$i_flowers_open == "4: 50-74%")], 
     breaks = 25, xlim = c(0, 365), main = "Parry's agave", xlab = "")
  # Saguaro and Palmer's data look decent. Goldenflower marginal. Almost no info
  # for the other spp (American and Parry's)

# Need to look at things proportionally, since the number of observations in each
# phenophase likely reflects effort as much as differences in prevalence....

# Combine P. palmeri and P. parryi? Parryi does seem to bloom a bit earlier,
# but maybe okay?


