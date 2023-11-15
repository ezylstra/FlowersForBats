# Flowering phenology, based on NPN data
# Erin Zylstra
# ezylstra@arizona.edu
# 2023-11-15

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

# Removing organpipe and 2 agave species from dataset (A. americana and deserti)
# (see notes and maps in explore-npn-observations.R), and give them a short name
dat <- dat %>%
  filter(!spp %in% c("American century plant", "desert agave", "organpipe cactus")) %>%
  rename(spp_common = spp) %>%
  mutate(spp = ifelse(spp_common == "saguaro", "C. gigantea",
                      ifelse(spp_common == "Palmer's century plant", "A. palmeri",
                             ifelse(spp_common == "Parry's agave", "A. parryi",
                                    "A. chrysanthus"))))

# Summarize information for each individual plant/patch (n = 649)
# Including some summaries of intensity data to potentially indicate whether or
# not patch was categorized correctly.
indiv <- dat %>%
  group_by(site_id, site_name, ind_id, plant_nickname, spp, species_id, patch,
           lat, lon, elev) %>%
  summarize(nyrs = length(unique(yr)),
            yr_first = min(yr),
            yr_last = max(yr),
            nobs = length(yr),
            flowers_up10 = sum(!is.na(i_flowers) & i_flowers %in% c("1: less than 3", "2: 3 to 10")),
            flowers_up100 = sum(!is.na(i_flowers) & i_flowers == "3: 11 to 100"),
            flowers_up1000 = sum(!is.na(i_flowers) & i_flowers == "4: 101 to 1,000"),
            flowers_up10000 = sum(!is.na(i_flowers) & i_flowers == "5: 1,001 to 10,000"),
            flowers_10000 = sum(!is.na(i_flowers) & i_flowers == "6: More than 10,000"),
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
  # A. parryi (1500, 1361 vs 1053 m)

# Summarize data for each species AND patch type
plantsp <- indiv %>%
  group_by(spp, patch) %>%
  summarize(ninds = length(ind_id),
            flowers_up10 = sum(flowers_up10 > 0),
            flowers_up100 = sum(flowers_up100 > 0), 
            flowers_up1000 = sum(flowers_up1000 > 0),
            flowers_up10000 = sum(flowers_up10000 > 0),
            flowers_10000 = sum(flowers_10000 > 0),
            .groups = "keep") %>%
  data.frame()
# People don't seem to be indicating that an individual is a patch very often...
filter(plantsp, spp == "A. palmeri")
  # only 15 of 292 with patch = 1. Lots of instances with patch = NA with 
  # number of flowers/infloresences > 100 or > 1,000.
filter(plantsp, spp == "A. parryi")
  # only 9 of 30 with patch = 1. Lots of instances with patch = NA with 
  # number of flowers/infloresences > 100 or > 1,000 and even > 10,0000
filter(plantsp, spp == "A. chrysanthus")
  # only 8 of 24 with patch = 1. Both patch = 1 and NA have instances with
  # number of flowers/infloresences > 10,0000
filter(plantsp, spp == "C. gigantea")
  # 5 of 249 with patch = 1 (which is expected/good). 
  # For both patch = 1 or NA, max number of flowers = 101 to 1000.

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
hist(ind_yr$nobs, breaks = 50, xlab = "Number of observations",
     main = "Number of observations of an individual per year")

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
hist(ind_yr$flowers_first, breaks = 50, xlim = c(0, 365), xlab = "Day of year",
     main = "First day of year with flowers")
hist(ind_yr$open_first, breaks = 50, xlim = c(0, 365), xlab = "Day of year",
     main = "First day of year with open flowers")

# First day of year with flowers, by species
par(mfrow = c(4, 1), mar = c(2, 4, 2, 1))
hist(ind_yr$flowers_first[ind_yr$spp == "C. gigantea"], breaks = 25, 
     xlim = c(0, 365), main = "C. gigantea", xlab = "")
hist(ind_yr$flowers_first[ind_yr$spp == "A. palmeri"], breaks = 25, 
     xlim = c(0, 365), main = "A. palmeri", xlab = "")
hist(ind_yr$flowers_first[ind_yr$spp == "A. parryi"], breaks = 25, 
     xlim = c(0, 365), main = "A. parryi", xlab = "")
hist(ind_yr$flowers_first[ind_yr$spp == "A. chrysanthus"], breaks = 25, 
     xlim = c(0, 365), main = "A. chrysanthus", xlab = "")

# First day of year with open flowers, by species
par(mfrow = c(4, 1), mar = c(2, 4, 2, 1))
hist(ind_yr$open_first[ind_yr$spp == "C. gigantea"], breaks = 25, 
     xlim = c(0, 365), main = "C. gigantea", xlab = "")
hist(ind_yr$open_first[ind_yr$spp == "A. palmeri"], breaks = 25, 
     xlim = c(0, 365), main = "A. palmeri", xlab = "")
hist(ind_yr$open_first[ind_yr$spp == "A. parryi"], breaks = 25, 
     xlim = c(0, 365), main = "A. parryi", xlab = "")
hist(ind_yr$open_first[ind_yr$spp == "A. chrysanthus"], breaks = 25, 
     xlim = c(0, 365), main = "A. chrysanthus", xlab = "")

# Plotting all observation dates when 50-74% of flowers are open, by species
par(mfrow = c(4, 1), mar = c(2, 4, 2, 1))
hist(dat$doy[which(dat$spp == "C. gigantea" & dat$i_flowers_open == "4: 50-74%")], 
     breaks = 25, xlim = c(0, 365), main = "C. gigantea", xlab = "")
hist(dat$doy[which(dat$spp == "A. palmeri" & dat$i_flowers_open == "4: 50-74%")], 
     breaks = 25, xlim = c(0, 365), main = "A. palmeri", xlab = "")
hist(dat$doy[which(dat$spp == "A. palmeri" & dat$i_flowers_open == "4: 50-74%")], 
     breaks = 25, xlim = c(0, 365), main = "A. parryi", xlab = "")
hist(dat$doy[which(dat$spp == "A. chrysanthus" & dat$i_flowers_open == "4: 50-74%")], 
     breaks = 25, xlim = c(0, 365), main = "A. chrysanthus", xlab = "")

# Need to look at things proportionally, since the number of observations in 
# each phenophase likely reflects effort as much as differences in prevalence...

# Create an index for week and a variable to group all agave together
dat <- dat %>%
  mutate(wk = isoweek(obsdate),
         spp2 = ifelse(spp == "C. gigantea", "saguaro", "agave"))

# Calculate proportions by week and species and year
prop_sppyr <- dat %>%
  group_by(spp, yr, wk) %>%
  summarize(nobs = length(spp),
            nobs_flower = sum(!is.na(flowers)),
            prop_flower = round(sum(flowers, na.rm = TRUE) / nobs_flower, 2),
            nobs_open = sum(!is.na(flowers_open)),
            prop_open = round(sum(flowers_open, na.rm = TRUE) / nobs_open, 2),
            .groups = "keep") %>%
  # Only keep spp/yr/week when there are XX observations?
  data.frame()

# Calculate proportions by week and genus and year
prop_genusyr <- dat %>%
  group_by(spp2, yr, wk) %>%
  summarize(nobs = length(spp),
            nobs_flower = sum(!is.na(flowers)),
            prop_flower = round(sum(flowers, na.rm = TRUE) / nobs_flower, 2),
            nobs_open = sum(!is.na(flowers_open)),
            prop_open = round(sum(flowers_open, na.rm = TRUE) / nobs_open, 2),
            .groups = "keep") %>%
  # Only keep genus/yr/week when there are XX observations?
  data.frame()
            

