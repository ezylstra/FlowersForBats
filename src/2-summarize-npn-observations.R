# Flowering/fruiting phenology, based on NPN data
# Erin Zylstra
# ezylstra@arizona.edu
# 2023-11-18

library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)
library(cowplot)
library(mgcv)
library(geodata)
library(terra)
library(tidyterra)

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
# hist(unique(dat$lat[dat$state == "AZ"]), breaks = 50)

# Restrict geographic area of interest
dat <- dat %>% 
  filter(state == "AZ", lat < 32.7) %>%
  select(-state)

# Removing organ pipe and 2 agave species from dataset (A. americana and deserti)
# (see notes and maps in explore-npn-observations.R), and give them a short name
dat <- dat %>%
  filter(!spp %in% c("American century plant", "desert agave", "organpipe cactus")) %>%
  rename(spp_common = spp) %>%
  mutate(spp = ifelse(spp_common == "saguaro", "C. gigantea",
                      ifelse(spp_common == "Palmer's century plant", "A. palmeri",
                             ifelse(spp_common == "Parry's agave", "A. parryi",
                                    "A. chrysanthus"))))

# Summarize information for each individual plant/patch
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
  # Most individuals only observed in 1 or 2 years, but up to 11

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
# write.table(select(plants, -c(species_id, lon_min, lon_max, lat_min, lat_max)),
#             "clipboard", sep = "\t", row.names = FALSE)
# Low % of individuals described as patches (2% for saguaro, 5-33% for agaves)
# A. palmeri and A. chrysantha with higher mean elevations than 
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

# Create an index for week and a variable to group all agave together
dat <- dat %>%
  mutate(wk = isoweek(obsdate),
         spp2 = ifelse(spp == "C. gigantea", "saguaro", "agave"))

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
mean(ind_yr$nobs); median(ind_yr$nobs); range(ind_yr$nobs)

# Number of observations per year by species
ind_yr %>%
  group_by(spp) %>%
  summarize(n_indiv = length(unique(ind_id)),
            n_indivyrs = length(spp),
            nobs_mn = round(mean(nobs), 1),
            nobs_md = median(nobs),
            nobs_min = min(nobs),
            nobs_max = max(nobs),
            ave_days_since_mn = mean(days_since_mn, na.rm = TRUE),
            n_oneobsperyr = sum(nobs == 1),
            n_multobsperyr = sum(nobs > 1),
            n_5obsperyr = sum(nobs > 4)) %>%
  mutate(propmult = round(n_multobsperyr / n_indivyrs, 2),
         prop5 = round(n_5obsperyr / n_indivyrs, 2)) %>%
  data.frame()

# First day of year with flowers, by species
# par(mfrow = c(4, 1), mar = c(2, 4, 2, 1))
# hist(ind_yr$flowers_first[ind_yr$spp == "C. gigantea"], breaks = 25, 
#      xlim = c(0, 365), main = "C. gigantea", xlab = "")
# hist(ind_yr$flowers_first[ind_yr$spp == "A. palmeri"], breaks = 25, 
#      xlim = c(0, 365), main = "A. palmeri", xlab = "")
# hist(ind_yr$flowers_first[ind_yr$spp == "A. parryi"], breaks = 25, 
#      xlim = c(0, 365), main = "A. parryi", xlab = "")
# hist(ind_yr$flowers_first[ind_yr$spp == "A. chrysanthus"], breaks = 25, 
#      xlim = c(0, 365), main = "A. chrysanthus", xlab = "")

# First day of year with open flowers, by species
# par(mfrow = c(4, 1), mar = c(2, 4, 2, 1))
# hist(ind_yr$open_first[ind_yr$spp == "C. gigantea"], breaks = 25, 
#      xlim = c(0, 365), main = "C. gigantea", xlab = "")
# hist(ind_yr$open_first[ind_yr$spp == "A. palmeri"], breaks = 25, 
#      xlim = c(0, 365), main = "A. palmeri", xlab = "")
# hist(ind_yr$open_first[ind_yr$spp == "A. parryi"], breaks = 25, 
#      xlim = c(0, 365), main = "A. parryi", xlab = "")
# hist(ind_yr$open_first[ind_yr$spp == "A. chrysanthus"], breaks = 25, 
#      xlim = c(0, 365), main = "A. chrysanthus", xlab = "")

# Plotting all observation dates when 50-74% of flowers are open, by species
# par(mfrow = c(4, 1), mar = c(2, 4, 2, 1))
# hist(dat$doy[which(dat$spp == "C. gigantea" & dat$i_flowers_open == "4: 50-74%")], 
#      breaks = 25, xlim = c(0, 365), main = "C. gigantea", xlab = "")
# hist(dat$doy[which(dat$spp == "A. palmeri" & dat$i_flowers_open == "4: 50-74%")], 
#      breaks = 25, xlim = c(0, 365), main = "A. palmeri", xlab = "")
# hist(dat$doy[which(dat$spp == "A. palmeri" & dat$i_flowers_open == "4: 50-74%")], 
#      breaks = 25, xlim = c(0, 365), main = "A. parryi", xlab = "")
# hist(dat$doy[which(dat$spp == "A. chrysanthus" & dat$i_flowers_open == "4: 50-74%")], 
#      breaks = 25, xlim = c(0, 365), main = "A. chrysanthus", xlab = "")

# Summarize information and plot locations of different plant species
# Extract information for each plant/patch
indiv <- dat %>%
  group_by(site_id, site_name, ind_id, plant_nickname, spp, species_id, patch, 
           lat, lon, elev) %>%
  summarize(nyrs = length(unique(yr)),
            yr_first = min(yr),
            yr_last = max(yr),
            nobs = length(yr),
            .groups = "keep") %>%
  data.frame()

# Summarize patch, location, elevation data
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
            nyrs_mn = mean(nyrs),
            .groups = "keep") %>%
  mutate(prop_patch = round(npatch / ninds, 2), .after = npatch) %>%
  data.frame()
plants

# Plots
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
  az <- us[us$NAME_1 == "Arizona",]
  
  # Add DEMs to look at elevation patterns
  dem_files <- list.files("data/dem", full.names = TRUE)
  dem_list <- NULL
  for (i in 1:length(dem_files)) {
    dem_list[[i]] <- rast(dem_files[i])
  }
  
  dem_coll <- sprc(dem_list)
  dem <- merge(dem_coll)
  rm(dem_files, dem_list, dem_coll)
  dem_crop <- crop(dem, c(-113, -109, 31, 32.7))

# Plot altogether
allspp <- ggplot() +
  geom_spatraster(data = dem_crop) +
  scale_fill_whitebox_c(palette = "arid", guide = "colorbar") +
  ylim(31, 32.7) +
  xlim(-113, -109) +
  geom_spatvector(data = az, fill = NA) +
  geom_point(data = indiv,
             aes(x = lon, y = lat, group = spp, color = spp)) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.box.spacing = unit(0.5, "pt"))
# ggsave("output/maps/allspp_saz.png",
#        allspp, device = "png", width = 6.5, height = 4.5,
#        units = "in", dpi = 300)

# Zoom into Tucson area
allspp_tucson <- ggplot() +
  geom_spatraster(data = dem_crop) +
  scale_fill_whitebox_c(palette = "arid") +
  ylim(31.8, 32.6) +
  xlim(-111.5, -110.4) +
  geom_spatvector(data = az, fill = NA) +
  geom_point(data = indiv,
             aes(x = lon, y = lat, group = spp, color = spp)) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.box.spacing = unit(0.5, "pt"))
# ggsave("output/maps/allspp_tucson.png",
#        allspp_tucson, device = "png", width = 6.5, height = 5.5,
#        units = "in", dpi = 300)

# Faceted, with DEMs
spp_facet <- ggplot() +
  geom_spatraster(data = dem_crop) +
  scale_fill_whitebox_c(palette = "soft") +
  geom_spatvector(data = az, fill = NA) +
  ylim(31.3, 32.7) +
  xlim(-111.7, -109) +
  geom_point(data = indiv, 
             aes(x = lon, y = lat)) +
  facet_wrap(~spp) +
  theme_bw() +
  labs(fill = "Elevation (m)", x = "", y = "") +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        legend.box.spacing = unit(0.1, "pt"))
spp_facet
# ggsave("output/maps/allspp_facet.png",
#        spp_facet, device = "png", width = 6.5, height = 5.5,
#        units = "in", dpi = 300)

# Where are most of the A. palmeri individuals?
palms <- filter(indiv, spp == "A. palmeri") %>%
  mutate(tucson = ifelse(lat > 32 & lat < 32.4 & lon > (-111.5) & lon < (-110.8), 1, 0),
         rincons = ifelse(lat > 32.12 & lat < 32.4 & lon > (-110.65) & lon < (-110.5), 1, 0))
count(palms, tucson, rincons)
# 5 in Tucson Valley (in town)
# 129 in the Rincons (along trails)
# 158 elsewhere
palms$rincons <- ifelse(palms$rincons == 1, 2, 0)
palms$geog <- factor(palms$tucson + palms$rincons)
ggplot() +
  geom_spatraster(data = dem_crop) +
  scale_fill_whitebox_c(palette = "soft") +
  geom_spatvector(data = az, fill = NA) +
  ylim(31.3, 32.7) +
  xlim(-111.7, -109) +
  geom_point(data = palms, aes(x = lon, y = lat, group = geog, color = geog)) +
  theme_bw() +
  labs(fill = "Elevation (m)", x = "", y = "") +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        legend.box.spacing = unit(0.1, "pt"))
