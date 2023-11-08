# Initial exploration of NPN data
# Erin Zylstra
# ezylstra@arizona.edu
# 2023-11-07

library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)
library(geodata)
library(terra)
library(tidyterra)

rm(list = ls())

# Load data
orig <- read.csv("data/FlowersforBats_RawObs_20230906.csv",
                 na.strings = c("", NA, -9999, "-9999"))

dim(orig) # 153,494 rows
str(orig)

# Extract information about species
spp <- orig %>%
  group_by(Species_ID, Genus, Species, Common_Name) %>%
  summarize(nobs = length(Species_ID), .groups = "keep") %>%
  data.frame()

# Remove extraneous columns and reformat/rename columns
dat <- orig %>%
  select(-c(Update_Datetime, Genus, Species, Kingdom, Abundance_Value)) %>%
  mutate(obsdate = mdy(Observation_Date),
         yr = as.numeric(year(obsdate))) %>%
  select(-Observation_Date) %>%
  rename(obs_id = Observation_ID,
         person_id = ObservedBy_Person_ID,
         group = Partner_Group,
         site_id = Site_ID,
         site_name = Site_Name,
         lat = Latitude,
         lon = Longitude,
         elev = Elevation_in_Meters,
         state = State,
         species_id = Species_ID,
         spp = Common_Name,
         ind_id = Individual_ID,
         plant_nickname = Plant_Nickname,
         patch = Patch,
         pheno_id = Phenophase_ID,
         pheno_name = Phenophase_Description,
         doy = Day_of_Year,
         pheno_status = Phenophase_Status,
         intensity_id = Intensity_Category_ID,
         intensity = Intensity_Value,
         comments = Observation_Comments) 

# Removing records with pheno_status == -1 (indicates pheno_status unknown)
  dat <- dat %>% filter(pheno_status != -1)
  
# Removing records in 2009-2011 (since many of them have different phenophases)
# that might result in the data being recorded slightly differently
  dat <- dat %>% filter(yr > 2011)

# Given so few records of cardon, will remove them
  dat <- dat %>% filter(spp != "Mexican giant cactus")

# Observations by state
  count(dat, state) 
  # Vast majority in AZ, 900+ in CA, 100+ in NM, 30 in TX, 20 in Sonora, 375 NA
  nostate <- filter(dat, is.na(state))
  nostate %>%
    select(group, site_id, site_name, lat, lon, 
           ind_id, plant_nickname, patch) %>%
    distinct()
  # All state = NA have group = Leslie Canyon NWR and similar lat/lon
  # Can change state to AZ
  dat$state[is.na(dat$state)] <- "AZ"
  rm(nostate)

# Phenophases
  count(dat, pheno_id, pheno_name)
  # Lots of records for all 5 phenophases

# Intensities
  count(dat, pheno_name, intensity_id, intensity)
  # intensity_id isn't useful (can have diff IDs for same value and phenophase)
  count(dat, intensity)
  # 10 records with "More than 1,000".
  filter(dat, intensity == "More than 1,000")
  count(filter(dat, spp == "saguaro"), pheno_name, intensity)
  # These are for flowers/fruits on saguaros where the "1,001 to 10,000" option 
  # doesn't seem to be available.

  # Going to create ordered levels that are more logical
  intensities <- data.frame(intensity = c("Less than 3",
                                          "3 to 10",
                                          "11 to 100",
                                          "101 to 1,000",
                                          "More than 1,000",
                                          "1,001 to 10,000",
                                          "More than 10,000",
                                          "Less than 5%",
                                          "5-24%",
                                          "25-49%",
                                          "50-74%",
                                          "75-94%",
                                          "95% or more")) %>%
    mutate(intensity_level = paste0(c(1:5, 5:6, 1:6), ": ", intensity))
  dat <- left_join(dat, intensities, by = "intensity") %>%
    select(-intensity_id)
  # Checks:
  # count(dat, pheno_name, intensity_level, intensity)
  # count(dat, intensity_level, intensity)

# Removing obs_id and comment columns and then remove any duplicate entries
  dat <- dat %>%
    select(-c(comments, obs_id)) %>%
    distinct()
  
# Groups
  count(dat, group) # 37 groups, some with thousands of observations
  # Some with > 10,000 rows (eg, McDowell, Tohono Chul, UA Campus Arb)
  # Flowers for Bats: 1832 obs
  # Borderlands Restoration: 1424
  # Chiricahua NM, Coronado NM, Fort Bowie NHS: each with 3300+

# Sites
  sites <- dat %>%
    group_by(site_id, site_name, lat, lon, elev, state) %>%
    summarize(nobs = length(site_id),
              nspp = length(unique(spp)),
              nind = length(unique(ind_id)),
              .groups = "keep") %>%
    data.frame()
  summary(sites)
  count(sites, nspp) # Most sites have 1 spp; range = 1-5
  count(sites, nind) # Most sites have 1 individual plant/patch; range = 1-14
  filter(sites, nind > 4) # 25 sites with 5 or more individuals
  filter(sites, grepl("Ranch", site_name)) 
    # Sands Ranch and Brown's Ranch, Jane Rau have a lot of data

# Individuals/Patches
  inds <- dat %>%
    group_by(site_id, site_name, lat, lon, 
             spp, ind_id, plant_nickname, patch) %>%
    summarize(nobs = length(ind_id),
              .groups = "keep") %>%
    data.frame()
  # Check that there's no individual associated with multiple sites and no
  # individual with patch = 1 and NA 
  nrow(inds); length(unique(inds$ind_id)) # ok
  
  # Distribution of species among individuals
  count(inds, spp) # 311 saguaros; 303 A. palmeri; 96 A. americana; 34 A. parryi
  
  # Distribution of "patch" for each spp
  count(inds, spp, patch)
  # Most saguaros have patch = NA (as expected); 5 have patch = 1
  # Most agaves have patch = NA (which isn't expected)
  
  filter(inds, spp == "saguaro" & patch == 1)
  count(filter(dat, spp == "saguaro" & patch == 1), 
        pheno_name, pheno_status, intensity_level)
  filter(dat, spp == "saguaro" & patch == 1 & intensity == "101 to 1,000")
  # I think they might be making observations at the patch level for multiple 
  # saguaros in a few instances (eg, ind_id = 16042 @ the UA Krutch Garden)
  
  # With respect to the patch entry for agaves, not sure that people are 
  # entering things correctly.
  
# Finding and removing observations of the same individual on the same date 
  obs <- dat %>%
    group_by(site_id, site_name, group, person_id, spp, ind_id, patch, 
             obsdate, yr) %>%
    summarize(.groups = "keep",
              phenos = length(pheno_id)) %>%
    data.frame()
  head(obs)
  count(obs, phenos)
  # 28,216 of 31,025 ind/date/person combos have 5 entries, which makes sense 
  # since there are 5 phenophases. BUT there are > 2,000 combos with 1-4 entries 
  # and 18 combos with 6-8

  obsp <- dat %>%
    group_by(site_id, site_name, group, person_id, spp, ind_id, patch, 
             obsdate, yr, pheno_id, pheno_name) %>%
    summarize(.groups = "keep",
              nobs = length(pheno_id)) %>%
    data.frame()
  head(obsp)
  count(obsp, nobs) # 32 combos with nobs = 2
    filter(obsp, nobs == 2)
    filter(dat, ind_id == 172139 & obsdate == "2019-08-09")
    # One observation with pheno_status = 1 with intensity values and one with 0s
    filter(dat, ind_id == 117053 & obsdate == "2021-10-27")
    # Different people making observations of the same individual on the same date
  
  # Get rid of rows in dat that are observations of the same ind plant/patch 
  # by the same person on the same day. Remove those with lower pheno_status and 
  # intensity value
  obsp <- obsp %>%
    arrange(site_id, person_id, spp, ind_id, obsdate, pheno_id) %>%
    mutate(obsnum = row_number())
  dat <- dat %>%
    arrange(site_id, person_id, spp, ind_id, obsdate, pheno_id, 
            desc(pheno_status), desc(intensity_level)) %>%
    left_join(select(obsp, site_id, person_id, spp, ind_id, obsdate, pheno_id, obsnum),
              by = c("site_id", "person_id", "spp", "ind_id", "obsdate", "pheno_id")) %>%
    mutate(dups = sequence(rle(as.character(obsnum))$lengths))
  dupi <- which(dat$dups == 2)
  dat[sort(c(dupi, dupi - 1)), c("site_name", "spp", "ind_id", "pheno_name", 
                                 "pheno_status", "intensity_level", 
                                 "obsdate", "obsnum", "dups")]
  # Quick look at these duplicates suggests it's fine to remove all the rows 
  # with dups = 2 (we're keeping the entries with more advanced phenology or
  # more information)
  dat <- dat %>%
    filter(dups == 1) %>%
    select(-dups)
  
  # Remove any identical observations of the same individual on the same date
  # by different people
  dat <- dat %>%
    arrange(site_id, group, person_id, spp, ind_id, obsdate, pheno_id, 
            desc(pheno_status), desc(intensity_level)) %>%
    distinct(ind_id, obsdate, pheno_id, pheno_status, intensity_level, 
             .keep_all = TRUE)
  
  # Now look at observations (that aren't identical) of the same individual on 
  # the same date by different people
  obsp <- dat %>%
    group_by(site_id, site_name, spp, ind_id, patch, obsdate, yr, 
             pheno_id, pheno_name) %>%
    summarize(.groups = "keep",
              nobs = length(pheno_id)) %>%
    data.frame() %>%
    arrange(site_id, spp, ind_id, obsdate, pheno_id) %>%
    mutate(obsnum = row_number())
  dat <- dat %>%
    select(-obsnum) %>%
    arrange(site_id, spp, ind_id, obsdate, pheno_id, 
            desc(pheno_status), desc(intensity_level)) %>%
    left_join(select(obsp, site_id, spp, ind_id, obsdate, pheno_id, obsnum),
              by = c("site_id", "spp", "ind_id", "obsdate", "pheno_id")) %>%
    mutate(dups = sequence(rle(as.character(obsnum))$lengths)) 
  count(dat, dups) 
  # Looks like we have 535 dups = 2; 18 dups = 3; and 2 dups = 4
  # Not sure what to do with these, other than to pick one observation...
  # For now, we'll keep the entries with more advanced phenology
  dat <- dat %>%
    filter(dups == 1) %>%
    select(-c(dups, obsnum))
  
  # Now there are 1-5 observed phenophases per individual & date

# Now working with reduced dataset --------------------------------------------#
  
# Look for phenophase status inconsistencies
  # Create a dataframe with one row for every observation (ind & date)
  indobs <- dat %>%
    group_by(site_id, site_name, spp, ind_id, patch, obsdate) %>%
    summarize(nobs = length(site_id), .groups = "keep") %>%
    mutate(indobs = paste0(ind_id, "_", obsdate)) %>%
    data.frame()
  
  indobs_full <- expand.grid(indobs = indobs$indobs,
                             pheno_name = sort(unique(dat$pheno_name)),
                             KEEP.OUT.ATTRS = FALSE)
  dat$indobs <- paste0(dat$ind_id, "_", dat$obsdate)
  
  indobs_full <- indobs_full %>%
    left_join(select(dat, c(indobs, pheno_name, pheno_status, intensity_level)),
              by = c("indobs", "pheno_name")) %>%
    left_join(select(indobs, -nobs), by = "indobs") %>%
    arrange(site_id, ind_id, spp, obsdate, pheno_name)
  
  indobs_ph <- indobs_full %>%
    select(-intensity_level) %>%
    pivot_wider(names_from = pheno_name,
                values_from = pheno_status) %>%
    data.frame() %>%
    rename(flowers = Flowers.or.flower.buds,
           fruit = Fruits,
           flowers_open = Open.flowers,
           fruit_drop = Recent.fruit.or.seed.drop,
           fruit_ripe = Ripe.fruits)

  count(indobs_ph, flowers, flowers_open)
  # Every combination of 0/1/NA exists....
  # If flowers = 0, flowers_open can only be 0 or NA
  # If flowers = NA, there's no useful information EXCEPT if flowers_open has 
  # an intensity value, then the flowers should probably be 1.
  
  count(indobs_ph, fruit, fruit_ripe)
  # Every combination of 0/1/NA exists....
  # If fruit = 0, fruit_ripe can only be 0 or NA
  # If fruit = NA, there's no useful information EXCEPT if fruit_rip has 
  # an intensity value, then the fruit should probably be 1.

  indobs_int <- indobs_full %>%
    select(-pheno_status) %>%
    pivot_wider(names_from = pheno_name,
                values_from = intensity_level) %>%
    data.frame() %>%
    rename(i_flowers = Flowers.or.flower.buds,
           i_fruit = Fruits,
           i_flowers_open = Open.flowers,
           i_fruit_drop = Recent.fruit.or.seed.drop,
           i_fruit_ripe = Ripe.fruits)
  
  indobs_w <- left_join(indobs_ph,
                        select(indobs_int, -c(site_id, site_name, spp, ind_id,
                                              patch, obsdate)),
                        by = "indobs")
  
  count(indobs_w, flowers, flowers_open, i_flowers, i_flowers_open)
  # CHANGES NEEDED:
    # flowers = 1 if flowers_open = 1 and i_flowers_open != NA
    # flowers_open = 0 if flowers = 0 and i_flowers_open = NA
    # flowers_open = NA if flowers = NA
  # UNINFORMATIVE: flowers = NA and flowers_open = 0 (n = 249)
  # UNINFORMATIVE: flowers = NA and flowers_open = NA (n = 13)

  count(indobs_w, fruit, fruit_ripe, i_fruit, i_fruit_ripe)
  # CHANGES NEEDED: 
    # fruit = 1 if fruit_ripe == 1 and i_fruit_ripe != NA
    # fruit_ripe = 0 if fruit = 0 and i_fruit_ripe = NA
    # fruit_ripe = NA if fruit = NA
  # UNINFORMATIVE: fruit = NA and fruit_ripe = 0 (n = 46)
  # UNINFORMATIVE: fruit = NA and fruit_ripe = NA (n = 528)

  # Fixing inconsistencies in the wide individual-observation dataframe:
  for (i in 1:nrow(indobs_w)) {
    if (!is.na(indobs_w$flowers_open[i]) & indobs_w$flowers_open[i] == 1 & 
        !is.na(indobs_w$i_flowers_open[i])) {
      indobs_w$flowers[i] <- 1
    }
    if (!is.na(indobs_w$flowers[i]) & indobs_w$flowers[i] == 0 & 
        !is.na(indobs_w$flowers_open[i]) & indobs_w$flowers_open[i] == 1 & 
        is.na(indobs_w$i_flowers_open[i])) {
      indobs_w$flowers_open[i] <- 0
    }
    if (is.na(indobs_w$flowers[i])) {
      indobs_w$flowers_open[i] <- NA
    }
    if (!is.na(indobs_w$fruit_ripe[i]) & indobs_w$fruit_ripe[i] == 1 & 
        !is.na(indobs_w$i_fruit_ripe[i])) {
      indobs_w$fruit[i] <- 1
    }
    if (!is.na(indobs_w$fruit[i]) & indobs_w$fruit[i] == 0 & 
        !is.na(indobs_w$fruit_ripe[i]) & indobs_w$fruit_ripe[i] == 1 & 
        is.na(indobs_w$i_fruit_ripe[i])) {
      indobs_w$fruit_ripe[i] <- 0
    }
    if (is.na(indobs_w$fruit[i])) {
      indobs_w$fruit_ripe[i] <- NA
    }
  }
  # Checks:
  # count(indobs_w, flowers, flowers_open, i_flowers, i_flowers_open)
  # count(indobs_w, fruit, fruit_ripe, i_fruit, i_fruit_ripe)

  # Add information from the original dataframe (dat) to indobs dataframe
  indobs_w <- indobs_w %>%
    left_join(select(dat, indobs, person_id, group, lat, lon, elev, state,
                     species_id, plant_nickname, doy, yr), 
              by = "indobs", multiple = "any") %>%
    arrange(indobs)
  
  # Give it an easier name
  dfw <- indobs_w
  rm(indobs_w)
  
# Summarize information and plot locations of different plant species
  # Extract information for each plant/patch
  indiv <- dfw %>%
    group_by(site_id, site_name, ind_id, plant_nickname, spp, species_id, patch, 
             lat, lon, elev, state) %>%
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
              .groups = "keep") %>%
    mutate(prop_patch = round(npatch / ninds, 2), .after = npatch) %>%
    data.frame()
  plants
  
  # Do the same, but only for species in southern Arizona (below 33 deg lat)
  indiv_saz <- indiv %>%
    filter(state == "AZ", lat < 32.9)
  plants_saz <- indiv_saz %>%
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
  plants_saz

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

  # All plants
  # ggplot() + 
  #   geom_spatvector(data = us.states) +
  #   geom_spatvector(data = mx.states) +
  #   geom_point(data = indiv, aes(x = lon, y = lat, group = spp, color = spp)) + 
  #   theme(legend.position = "bottom",
  #         legend.title = element_blank())
  
  # Just southern Arizona
  az <- us[us$NAME_1 == "Arizona",]

  # All on one plot
  # ggplot() + 
  #   geom_spatvector(data = az) + 
  #   ylim(31, 32.7) +
  #   xlim(-113, -109) +
  #   geom_point(data = indiv_saz, 
  #              aes(x = lon, y = lat, group = spp, color = spp)) +
  #   theme_bw() +
  #   theme(legend.position = "bottom",
  #         legend.title = element_blank())
  
  # Faceted
  # ggplot() + 
  #   geom_spatvector(data = az) + 
  #   ylim(31, 32.7) +
  #   xlim(-113, -109) +
  #   geom_point(data = indiv_saz, aes(x = lon, y = lat)) +
  #   facet_wrap(~spp) +
  #   theme_bw()

# Add DEMs to look at elevation patterns
  # dem_files <- list.files("data/dem", full.names = TRUE) 
  # dem_list <- NULL
  # for (i in 1:length(dem_files)) {
  #   dem_list[[i]] <- rast(dem_files[i])
  # }
  # 
  # dem_coll <- sprc(dem_list)
  # dem <- merge(dem_coll)
  # dem <- project(dem, crs(us.states)) # takes a minute
  # rm(dem_files, dem_list, dem_coll)
  
  # Zoom into Tucson area
  # ggplot() + 
  #   geom_spatraster(data = dem) +
  #   ylim(31.8, 32.6) +
  #   xlim(-111.5, -110.4) +
  #   geom_spatvector(data = az, fill = NA) + 
  #   geom_point(data = indiv_saz, 
  #              aes(x = lon, y = lat, group = spp, color = spp)) +
  #   theme_bw() +
  #   theme(legend.position = "right", 
  #         legend.title = element_blank())
 
# Save cleaned up observation dataframe for use in other scripts
# write.csv(dfw, "data/FlowersForBats_CleanedObs.csv", row.names = FALSE)
