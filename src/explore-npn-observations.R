# Initial exploration of NPN data
# Erin Zylstra
# ezylstra@arizona.edu
# 2023-11-02

library(dplyr)
library(lubridate)

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

# Given so few records of cardon, will remove them
  dat <- dat %>% filter(spp != "Mexican giant cactus")

# Observations by state
  count(dat, state) 
  # Vast majority in AZ, 1000+ in CA, 100+ in NM, 30 in TX, 20 in Sonora, 420 NA
  nostate <- filter(dat, is.na(state))
  nostate %>%
    select(group, site_id, site_name, lat, lon, 
           ind_id, plant_nickname, patch) %>%
    distinct()
  # All state = NA have group = Leslie Canyon NWR and similar lat/lon
  # Can change state to AZ
  dat$state[is.na(dat$state)] <- "AZ"

# Phenophases
  count(dat, pheno_id, pheno_name)
  # Look just at those phases that have (1 location) in their name
  dat %>%
    filter(grepl("(1 location)", pheno_name)) %>%
    count(yr)
  count(dat, yr)
  # All observations in 2009 and 2010 are the (1 location) phenophases
  # Some observations in 2011 are the (1 location) phenophases
  # Combine these phases? (will see if it's worth it given the number of sites
  # observed in early years)

# Intensities
  count(dat, pheno_name, intensity_id, intensity)
  
# Groups
  count(dat, group) # 38 groups, some with thousands of observations
  # Some with > 10,000 rows (eg, McDowell, Tohono Chul, UA Campus Arb)
  # Flowers for Bats: 1846 obs
  # Borderlands Restoration: 1506
  # Chiricahua NM, Coronado NM, Fort Bowie NHS: each with 3500+

# Patch
  count(dat, patch)
  # 8806 values are 1, rest/majority are all NA
  count(dat, spp, patch)

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
  count(sites, nind) # Most sites have 1 ind; range = 1-14
  filter(sites, nind > 4)
  filter(sites, grepl("Ranch", site_name))

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
  count(inds, spp) # 330 saguaros; 303 A. palmeri; 97 A. americana; 34 A. parryi
  
  # Distribution of "patch" for each spp
  count(inds, spp, patch)
  # Most saguaros have patch = NA (as expected); 5 have patch = 1
  # Most agaves have patch = NA (which isn't expected)
  
  filter(inds, spp == "saguaro" & patch == 1)
  count(filter(dat, spp == "saguaro" & patch == 1), 
        pheno_name, pheno_status, intensity)
  filter(dat, spp == "saguaro" & patch ==1 & intensity == "101 to 1,000")
  # I think they might be making observations at the patch level for multiple 
  # saguaros in a few instances (eg, ind_id = 16042 @ the UA Krutch Garden)
  
  # With respect to the patch entry for agaves, not sure that people are 
  # entering things correctly.
  
# Phenophase status (excluding those old "1 location" phases)
  obs <- dat %>%
    filter(!grepl("(1 location)", pheno_name)) %>%
    group_by(site_id, site_name, group, person_id, spp, ind_id, patch, 
             obsdate, yr) %>%
    summarize(.groups = "keep",
              phenos = length(pheno_id)) %>%
    data.frame()
  head(obs)
  count(obs, phenos)
  # 28,770 of 31,209 ind/date combos have 5 entries, which makes sense since 
  # there are 5 phenophases (flowers, open flowers, fruits, ripe fruits, drop)
  # BUT there are > 2,000 combos with 1-4 entries and 200 combos with 6-15
  
  obsp <- dat %>%
    filter(!grepl("(1 location)", pheno_name)) %>%
    group_by(site_id, site_name, group, person_id, spp, ind_id, patch, 
             obsdate, yr, pheno_id, pheno_name) %>%
    summarize(.groups = "keep",
              phenos = length(pheno_id)) %>%
    data.frame()
  
  count(obsp, phenos) # 951 combos with phenos = 2, 21 with phenos = 3
  filter(obsp, phenos == 3)
  filter(dat, ind_id == 74154 & obsdate == "2014-11-17")
  # Are these just duplicate entries/observations?
  
  # Removing what are essentially duplicate entries/observations...
  dat <- dat %>%
    distinct(ind_id, obsdate, pheno_id, pheno_status, intensity, 
             .keep_all = TRUE)
  
  obsp <- dat %>%
    filter(!grepl("(1 location)", pheno_name)) %>%
    group_by(site_id, site_name, group, person_id, spp, ind_id, patch, 
             obsdate, yr, pheno_id, pheno_name) %>%
    summarize(.groups = "keep",
              phenos = length(pheno_id)) %>%
    data.frame()
  
  count(obsp, phenos) # Now just 40 combos with phenos = 2
  filter(obsp, phenos > 1)
  filter(dat, ind_id == 172139 & obsdate == "2019-08-09") 
    # In this case, 2 entries for each ind/date/phase, pheno_status = 0 and 1
  filter(dat, ind_id == 200993 & obsdate == "2020-09-13") 
    # Same here, except a few have pheno_status = -1
  
  # TODO: make sure that we've removed all duplicates and gotten rid of
  # questionable entries (eg, pheno_status = -1?)

