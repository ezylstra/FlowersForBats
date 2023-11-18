# Flowering phenology, based on NPN data
# Erin Zylstra
# ezylstra@arizona.edu
# 2023-11-16

library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)
library(cowplot)
library(mgcv)

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

# Need to look at things proportionally, since the number of observations in 
# each phenophase likely reflects effort as much as differences in prevalence...

# Create an index for week and a variable to group all agave together
dat <- dat %>%
  mutate(wk = isoweek(obsdate),
         spp2 = ifelse(spp == "C. gigantea", "saguaro", "agave"))

# ---------------------------------------------------------------------------- #
# Calculate proportions of plants/patches in each phenophase by week and year
# and use GAMs to describe seasonal patterns
# ---------------------------------------------------------------------------- #

# First, will want to remove multiple observations of the same individual in
# the same week. Sort so the more advanced phenophase gets kept (if more 
# than one value in a week)

prop_sppyr <- dat %>%
  arrange(spp, ind_id, yr, wk, desc(flowers), desc(flowers_open)) %>%
  distinct(ind_id, yr, wk, .keep_all = TRUE) %>%
  group_by(spp, yr, wk) %>%
  summarize(nobs = length(spp),
            nobs_flower = sum(!is.na(flowers)),
            prop_flower = sum(flowers, na.rm = TRUE) / nobs_flower,
            nobs_open = sum(!is.na(flowers_open)),
            prop_open = sum(flowers_open, na.rm = TRUE) / nobs_open,
            .groups = "keep") %>%
  mutate(date_generic = parse_date_time(paste(2024, wk, 1, sep="/"), "Y/W/w"),
         date_generic = as.Date(date_generic),
         doy = yday(date_generic),
         fyear = factor(yr, levels = as.character(2012:2023))) %>% 
  data.frame()

# Create a dataframe with proportions of all agaves (across species) in each 
# phenophase
prop_genusyr <- dat %>%
  arrange(spp2, ind_id, yr, wk, desc(flowers), desc(flowers_open)) %>%
  distinct(ind_id, yr, wk, .keep_all = TRUE) %>%
  group_by(spp2, yr, wk) %>%
  summarize(nobs = length(spp),
            nobs_flower = sum(!is.na(flowers)),
            prop_flower = sum(flowers, na.rm = TRUE) / nobs_flower,
            nobs_open = sum(!is.na(flowers_open)),
            prop_open = sum(flowers_open, na.rm = TRUE) / nobs_open,
            .groups = "keep") %>%
  mutate(date_generic = parse_date_time(paste(2024, wk, 1, sep="/"), "Y/W/w"),
         date_generic = as.Date(date_generic),
         doy = yday(date_generic),
         fyear = factor(yr, levels = as.character(2012:2023))) %>% 
  data.frame()

# Saguaros ------------------------------------------------------------------- #
sagdat <- filter(prop_sppyr, spp == "C. gigantea")

# Flowers or buds
  # Simple model with same seasonal smooth each year
  sagfl_m1 <- gam(prop_flower ~ s(doy), weights = nobs_flower,
                  data = sagdat, method = "REML", family = binomial)
  # Model with an independent seasonal smooth each year
  sagfl_m2 <- gam(prop_flower ~ s(doy, by = fyear), weights = nobs_flower,
                  data = sagdat, method = "REML", family = binomial)
  # Compare models
  summary(sagfl_m1)
  summary(sagfl_m2)
  AIC(sagfl_m1, sagfl_m2)
    # Model with different annual curves is better (by a lot)

  # Create dataframe for predictions from models with annual curves
  sag_preds <- data.frame(fyear = NA, wks = NA)
  for (yr in sort(unique(sagdat$yr))) {
    min_doy <- min(sagdat$doy[sagdat$yr == yr])
    max_doy <- max(sagdat$doy[sagdat$yr == yr])
    max_doy <- ifelse(max_doy > 350, 365, max_doy)
    doys <- seq(min_doy, max_doy)
    if (yr == min(sagdat$yr)) {
      sag_preds <- data.frame(fyear = yr, doy = doys)
    } else {
      sag_preds <- rbind(sag_preds, 
                         data.frame(fyear = yr, doy = doys))
    }
  }
  sag_preds$fyear <- factor(sag_preds$fyear, levels = as.character(2012:2023))
  
  # Make predictions based on more complex model
  sagfl_preds <- cbind(sag_preds,
                       as.data.frame(predict(sagfl_m2, 
                                             newdata = sag_preds,
                                             type = "response", 
                                             se.fit = TRUE))) %>%
    mutate(lcl = fit - 1.96 * se.fit,
           ucl = fit + 1.96 * se.fit)
  # Plot predictions
  sagfl_plotyr <- ggplot(data = sagfl_preds, 
         aes(x = doy, y = fit, group = fyear, color = fyear)) +
    scale_color_discrete(name = "Year") +
    geom_line() +
    labs(y = "Proportion of plants with flowers or buds",
         x = "Day of year") +
    annotate("text", x = 365, y = 0.98, label = "C. gigantea",
             hjust = 1, vjust = 1, fontface = 2)
  
  # Make predictions based on simpler model
  sagfl_preds1 <- data.frame(doy = 1:365)
  sagfl_preds1 <- cbind(sagfl_preds1,
                        as.data.frame(predict(sagfl_m1, 
                                             newdata = sagfl_preds1,
                                             type = "response", 
                                             se.fit = TRUE))) %>%
    mutate(lcl = fit - 1.96 * se.fit,
           ucl = fit + 1.96 * se.fit)
  # Plot predictions (with raw data)
  sagfl_plot <- ggplot(sagfl_preds1, aes(x = doy)) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), color = "gray", alpha = 0.3) +
    geom_line(aes(y = fit)) +
    labs(y = "Proportion of plants with flowers or buds",
         x = "Day of year") +
    geom_point(data = sagdat, 
               aes(x = doy, y = prop_flower, group = fyear, color = fyear),
               alpha = 0.4) +
    geom_line(data = sagdat, 
              aes(x = doy, y = prop_flower, group = fyear, color = fyear),
              alpha = 0.4) +
    scale_color_discrete(name = "Year") +
    annotate("text", x = 365, y = 0.98, label = "C. gigantea",
             hjust = 1, vjust = 1, fontface = 2)
    
# Open flowers
  # Simple model with same seasonal smooth each year
  sagof_m1 <- gam(prop_open ~ s(doy), weights = nobs_open,
                  data = sagdat, method = "REML", family = binomial)
  # Model with an independent seasonal smooth each year
  sagof_m2 <- gam(prop_open ~ s(doy, by = fyear), weights = nobs_open,
                  data = sagdat, method = "REML", family = binomial)
  # Compare models
  summary(sagof_m1)
  summary(sagof_m2)
  AIC(sagof_m1, sagof_m2)
    # Model with different annual curves is better (by a lot)
  
  # Make predictions based on more complex model
  sagof_preds <- cbind(sag_preds,
                       as.data.frame(predict(sagof_m2, 
                                             newdata = sag_preds,
                                             type = "response", 
                                             se.fit = TRUE))) %>%
    mutate(lcl = fit - 1.96 * se.fit,
           ucl = fit + 1.96 * se.fit)
  # Plot predictions
  sagof_plotyr <- ggplot(data = sagof_preds, 
                         aes(x = doy, y = fit, group = fyear, color = fyear)) +
    scale_color_discrete(name = "Year") +
    geom_line() +
    labs(y = "Proportion of flowers that are open",
         x = "Day of year") +
    annotate("text", x = 365, y = 0.98, label = "C. gigantea",
             hjust = 1, vjust = 1, fontface = 2)
  
  # Make predictions based on all years together
  sagof_preds1 <- data.frame(doy = 1:365)
  sagof_preds1 <- cbind(sagof_preds1,
                        as.data.frame(predict(sagof_m1, 
                                              newdata = sagof_preds1,
                                              type = "response", 
                                              se.fit = TRUE))) %>%
    mutate(lcl = fit - 1.96 * se.fit,
           ucl = fit + 1.96 * se.fit)
  # Plot predictions (with raw data)
  sagof_plot <- ggplot(sagof_preds1, aes(x = doy)) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), color = "gray", alpha = 0.3) +
    geom_line(aes(y = fit)) +
    labs(y = "Proportion of flowers that are open",
         x = "Day of year") +
    geom_point(data = sagdat, 
               aes(x = doy, y = prop_open, group = fyear, color = fyear),
               alpha = 0.4) +
    geom_line(data = sagdat, 
              aes(x = doy, y = prop_open, group = fyear, color = fyear),
              alpha = 0.4) +
    scale_color_discrete(name = "Year") +
    annotate("text", x = 365, y = 0.98, label = "C. gigantea",
             hjust = 1, vjust = 1, fontface = 2)

  # Look at plots with seasonal curve for all years
  plot_grid(sagfl_plot, sagof_plot, ncol = 1)
  # Look at plots with independent curves each year
  plot_grid(sagfl_plotyr, sagof_plotyr, ncol = 1)

# A. palmeri ----------------------------------------------------------------- #
palmdat <- filter(prop_sppyr, spp == "A. palmeri")
  
# Flowers or buds
  # Simple model with same seasonal smooth each year
  palmfl_m1 <- gam(prop_flower ~ s(doy), weights = nobs_flower,
                  data = palmdat, method = "REML", family = binomial)
  # Model with an independent seasonal smooth each year
  palmfl_m2 <- gam(prop_flower ~ s(doy, by = fyear), weights = nobs_flower,
                  data = palmdat, method = "REML", family = binomial)
  # Compare models
  summary(palmfl_m1)
  summary(palmfl_m2)
  AIC(palmfl_m1, palmfl_m2)
    # Model with different annual curves is better (by a lot)
  
  # Create dataframe for predictions from models with annual curves
  palm_preds <- data.frame(fyear = NA, wks = NA)
  for (yr in sort(unique(palmdat$yr))) {
    min_doy <- min(palmdat$doy[palmdat$yr == yr])
    max_doy <- max(palmdat$doy[palmdat$yr == yr])
    max_doy <- ifelse(max_doy > 350, 365, max_doy)
    doys <- seq(min_doy, max_doy)
    if (yr == min(palmdat$yr)) {
      palm_preds <- data.frame(fyear = yr, doy = doys)
    } else {
      palm_preds <- rbind(palm_preds, 
                         data.frame(fyear = yr, doy = doys))
    }
  }
  palm_preds$fyear <- factor(palm_preds$fyear, levels = as.character(2012:2023))
  
  # Make predictions based on more complex model
  palmfl_preds <- cbind(palm_preds,
                       as.data.frame(predict(palmfl_m2, 
                                             newdata = palm_preds,
                                             type = "response", 
                                             se.fit = TRUE))) %>%
    mutate(lcl = fit - 1.96 * se.fit,
           ucl = fit + 1.96 * se.fit)
  # Plot predictions
  palmfl_plotyr <- ggplot(data = palmfl_preds, 
                         aes(x = doy, y = fit, group = fyear, color = fyear)) +
    scale_color_discrete(name = "Year") +
    geom_line() +
    labs(y = "Proportion of patches with flowers or buds",
         x = "Day of year") +
    annotate("text", x = 365, y = 0.98, label = "A. palmeri",
             hjust = 1, vjust = 1, fontface = 2)
  
  # Make predictions based on simpler model
  palmfl_preds1 <- data.frame(doy = 1:365)
  palmfl_preds1 <- cbind(palmfl_preds1,
                        as.data.frame(predict(palmfl_m1, 
                                              newdata = palmfl_preds1,
                                              type = "response", 
                                              se.fit = TRUE))) %>%
    mutate(lcl = fit - 1.96 * se.fit,
           ucl = fit + 1.96 * se.fit)
  # Plot predictions (with raw data)
  palmfl_plot <- ggplot(palmfl_preds1, aes(x = doy)) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), color = "gray", alpha = 0.3) +
    geom_line(aes(y = fit)) +
    labs(y = "Proportion of patches with flowers or buds",
         x = "Day of year") +
    geom_point(data = palmdat, 
               aes(x = doy, y = prop_flower, group = fyear, color = fyear),
               alpha = 0.4) +
    geom_line(data = palmdat, 
              aes(x = doy, y = prop_flower, group = fyear, color = fyear),
              alpha = 0.4) +
    scale_color_discrete(name = "Year") +
    annotate("text", x = 365, y = 0.98, label = "A. palmeri",
             hjust = 1, vjust = 1, fontface = 2)
  
# Open flowers
  # Simple model with same seasonal smooth each year
  palmof_m1 <- gam(prop_open ~ s(doy), weights = nobs_open,
                  data = palmdat, method = "REML", family = binomial)
  # Model with an independent seasonal smooth each year
  palmof_m2 <- gam(prop_open ~ s(doy, by = fyear), weights = nobs_open,
                  data = palmdat, method = "REML", family = binomial)
  # Compare models
  summary(palmof_m1)
  summary(palmof_m2)
  AIC(palmof_m1, palmof_m2)
    # Model with different annual curves is better (by a lot)
  
  # Make predictions based on more complex model
  palmof_preds <- cbind(palm_preds,
                       as.data.frame(predict(palmof_m2, 
                                             newdata = palm_preds,
                                             type = "response", 
                                             se.fit = TRUE))) %>%
    mutate(lcl = fit - 1.96 * se.fit,
           ucl = fit + 1.96 * se.fit)
  # Plot predictions
  palmof_plotyr <- ggplot(data = palmof_preds, 
                         aes(x = doy, y = fit, group = fyear, color = fyear)) +
    scale_color_discrete(name = "Year") +
    geom_line() +
    labs(y = "Proportion of flowers that are open",
         x = "Day of year") +
    annotate("text", x = 365, y = 0.98, label = "A. palmeri",
             hjust = 1, vjust = 1, fontface = 2)
  
  # Make predictions based on all years together
  palmof_preds1 <- data.frame(doy = 1:365)
  palmof_preds1 <- cbind(palmof_preds1,
                        as.data.frame(predict(palmof_m1, 
                                              newdata = palmof_preds1,
                                              type = "response", 
                                              se.fit = TRUE))) %>%
    mutate(lcl = fit - 1.96 * se.fit,
           ucl = fit + 1.96 * se.fit)
  # Plot predictions (with raw data)
  palmof_plot <- ggplot(palmof_preds1, aes(x = doy)) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), color = "gray", alpha = 0.3) +
    geom_line(aes(y = fit)) +
    labs(y = "Proportion of flowers that are open",
         x = "Day of year") +
    geom_point(data = palmdat, 
               aes(x = doy, y = prop_open, group = fyear, color = fyear),
               alpha = 0.4) +
    geom_line(data = palmdat, 
              aes(x = doy, y = prop_open, group = fyear, color = fyear),
              alpha = 0.4) +
    scale_color_discrete(name = "Year") +
    annotate("text", x = 365, y = 0.98, label = "A. palmeri",
             hjust = 1, vjust = 1, fontface = 2)
  
  # Look at plots with seasonal curve for all years
  plot_grid(palmfl_plot, palmof_plot, ncol = 1)
  # Look at plots with independent curves each year
  plot_grid(palmfl_plotyr, palmof_plotyr, ncol = 1)

# A. parryi ------------------------------------------------------------------ #
parrdat <- filter(prop_sppyr, spp == "A. parryi")
  
# Flowers or buds
  # Simple model with same seasonal smooth each year
  parrfl_m1 <- gam(prop_flower ~ s(doy), weights = nobs_flower,
                   data = parrdat, method = "REML", family = binomial)
  # Model with an independent seasonal smooth each year
  parrfl_m2 <- gam(prop_flower ~ s(doy, by = fyear), weights = nobs_flower,
                   data = parrdat, method = "REML", family = binomial)
  # Compare models
  summary(parrfl_m1)
  summary(parrfl_m2)
  AIC(parrfl_m1, parrfl_m2)
    # Model with different annual curves is better
  
  # Create dataframe for predictions from models with annual curves
  parr_preds <- data.frame(fyear = NA, wks = NA)
  for (yr in sort(unique(parrdat$yr))) {
    min_doy <- min(parrdat$doy[parrdat$yr == yr])
    max_doy <- max(parrdat$doy[parrdat$yr == yr])
    max_doy <- ifelse(max_doy > 350, 365, max_doy)
    doys <- seq(min_doy, max_doy)
    if (yr == min(parrdat$yr)) {
      parr_preds <- data.frame(fyear = yr, doy = doys)
    } else {
      parr_preds <- rbind(parr_preds, 
                          data.frame(fyear = yr, doy = doys))
    }
  }
  parr_preds$fyear <- factor(parr_preds$fyear, levels = as.character(2012:2023))
  
  # Make predictions based on more complex model
  parrfl_preds <- cbind(parr_preds,
                        as.data.frame(predict(parrfl_m2, 
                                              newdata = parr_preds,
                                              type = "response", 
                                              se.fit = TRUE))) %>%
    mutate(lcl = fit - 1.96 * se.fit,
           ucl = fit + 1.96 * se.fit)
  # Plot predictions
  parrfl_plotyr <- ggplot(data = parrfl_preds, 
                          aes(x = doy, y = fit, group = fyear, color = fyear)) +
    scale_color_discrete(name = "Year") +
    geom_line() +
    labs(y = "Proportion of patches with flowers or buds",
         x = "Day of year") +
    annotate("text", x = 365, y = 0.98, label = "A. parryi",
             hjust = 1, vjust = 1, fontface = 2)
  
  # Make predictions based on simpler model
  parrfl_preds1 <- data.frame(doy = 1:365)
  parrfl_preds1 <- cbind(parrfl_preds1,
                         as.data.frame(predict(parrfl_m1, 
                                               newdata = parrfl_preds1,
                                               type = "response", 
                                               se.fit = TRUE))) %>%
    mutate(lcl = fit - 1.96 * se.fit,
           ucl = fit + 1.96 * se.fit)
  # Plot predictions (with raw data)
  parrfl_plot <- ggplot(parrfl_preds1, aes(x = doy)) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), color = "gray", alpha = 0.3) +
    geom_line(aes(y = fit)) +
    labs(y = "Proportion of patches with flowers or buds",
         x = "Day of year") +
    geom_point(data = parrdat, 
               aes(x = doy, y = prop_flower, group = fyear, color = fyear),
               alpha = 0.4) +
    geom_line(data = parrdat, 
              aes(x = doy, y = prop_flower, group = fyear, color = fyear),
              alpha = 0.4) +
    scale_color_discrete(name = "Year") +
    annotate("text", x = 365, y = 0.98, label = "A. parryi",
             hjust = 1, vjust = 1, fontface = 2)
  
# Open flowers
  # Simple model with same seasonal smooth each year
  parrof_m1 <- gam(prop_open ~ s(doy), weights = nobs_open,
                   data = parrdat, method = "REML", family = binomial)
  # Model with an independent seasonal smooth each year
  parrof_m2 <- gam(prop_open ~ s(doy, by = fyear), weights = nobs_open,
                   data = parrdat, method = "REML", family = binomial)
  # Compare models
  summary(parrof_m1)
  summary(parrof_m2)
  AIC(parrof_m1, parrof_m2)
    # Model with different annual curves is better
  
  # Make predictions based on more complex model
  parrof_preds <- cbind(parr_preds,
                        as.data.frame(predict(parrof_m2, 
                                              newdata = parr_preds,
                                              type = "response", 
                                              se.fit = TRUE))) %>%
    mutate(lcl = fit - 1.96 * se.fit,
           ucl = fit + 1.96 * se.fit)
  # Plot predictions
  parrof_plotyr <- ggplot(data = parrof_preds, 
                          aes(x = doy, y = fit, group = fyear, color = fyear)) +
    scale_color_discrete(name = "Year") +
    geom_line() +
    labs(y = "Proportion of flowers that are open",
         x = "Day of year") +
    annotate("text", x = 365, y = 0.98, label = "A. parryi",
             hjust = 1, vjust = 1, fontface = 2)
  
  # Make predictions based on all years together
  parrof_preds1 <- data.frame(doy = 1:365)
  parrof_preds1 <- cbind(parrof_preds1,
                         as.data.frame(predict(parrof_m1, 
                                               newdata = parrof_preds1,
                                               type = "response", 
                                               se.fit = TRUE))) %>%
    mutate(lcl = fit - 1.96 * se.fit,
           ucl = fit + 1.96 * se.fit)
  # Plot predictions (with raw data)
  parrof_plot <- ggplot(parrof_preds1, aes(x = doy)) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), color = "gray", alpha = 0.3) +
    geom_line(aes(y = fit)) +
    labs(y = "Proportion of flowers that are open",
         x = "Day of year") +
    geom_point(data = parrdat, 
               aes(x = doy, y = prop_open, group = fyear, color = fyear),
               alpha = 0.4) +
    geom_line(data = parrdat, 
              aes(x = doy, y = prop_open, group = fyear, color = fyear),
              alpha = 0.4) +
    scale_color_discrete(name = "Year") +
    annotate("text", x = 365, y = 0.98, label = "A. parryi",
             hjust = 1, vjust = 1, fontface = 2)
  
  # Look at plots with seasonal curve for all years
  plot_grid(parrfl_plot, parrof_plot, ncol = 1)
  # Look at plots with independent curves each year
  plot_grid(parrfl_plotyr, parrof_plotyr, ncol = 1)  

# A. chrysanthus ------------------------------------------------------------- #
chrydat <- filter(prop_sppyr, spp == "A. chrysanthus")
  
# Flowers or buds
  # Simple model with same seasonal smooth each year
  chryfl_m1 <- gam(prop_flower ~ s(doy), weights = nobs_flower,
                   data = chrydat, method = "REML", family = binomial)
  # Model with an independent seasonal smooth each year
  chryfl_m2 <- gam(prop_flower ~ s(doy, by = fyear), weights = nobs_flower,
                   data = chrydat, method = "REML", family = binomial)
  # Compare models
  summary(chryfl_m1)
  summary(chryfl_m2)
  AIC(chryfl_m1, chryfl_m2)
    # Model with different annual curves is better
  
  # Create dataframe for predictions from models with annual curves
  chry_preds <- data.frame(fyear = NA, wks = NA)
  for (yr in sort(unique(chrydat$yr))) {
    min_doy <- min(chrydat$doy[chrydat$yr == yr])
    max_doy <- max(chrydat$doy[chrydat$yr == yr])
    max_doy <- ifelse(max_doy > 350, 365, max_doy)
    doys <- seq(min_doy, max_doy)
    if (yr == min(chrydat$yr)) {
      chry_preds <- data.frame(fyear = yr, doy = doys)
    } else {
      chry_preds <- rbind(chry_preds, 
                          data.frame(fyear = yr, doy = doys))
    }
  }
  chry_preds$fyear <- factor(chry_preds$fyear, levels = as.character(2012:2023))
  
  # Make predictions based on more complex model
  chryfl_preds <- cbind(chry_preds,
                        as.data.frame(predict(chryfl_m2, 
                                              newdata = chry_preds,
                                              type = "response", 
                                              se.fit = TRUE))) %>%
    mutate(lcl = fit - 1.96 * se.fit,
           ucl = fit + 1.96 * se.fit)
  # Plot predictions
  chryfl_plotyr <- ggplot(data = chryfl_preds, 
                          aes(x = doy, y = fit, group = fyear, color = fyear)) +
    scale_color_discrete(name = "Year") +
    geom_line() +
    labs(y = "Proportion of patches with flowers or buds",
         x = "Day of year") +
    annotate("text", x = 365, y = 0.98, label = "A. chrysanthus",
             hjust = 1, vjust = 1, fontface = 2)
  
  # Make predictions based on simpler model
  chryfl_preds1 <- data.frame(doy = 1:365)
  chryfl_preds1 <- cbind(chryfl_preds1,
                         as.data.frame(predict(chryfl_m1, 
                                               newdata = chryfl_preds1,
                                               type = "response", 
                                               se.fit = TRUE))) %>%
    mutate(lcl = fit - 1.96 * se.fit,
           ucl = fit + 1.96 * se.fit)
  # Plot predictions (with raw data)
  chryfl_plot <- ggplot(chryfl_preds1, aes(x = doy)) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), color = "gray", alpha = 0.3) +
    geom_line(aes(y = fit)) +
    labs(y = "Proportion of patches with flowers or buds",
         x = "Day of year") +
    geom_point(data = chrydat, 
               aes(x = doy, y = prop_flower, group = fyear, color = fyear),
               alpha = 0.4) +
    geom_line(data = chrydat, 
              aes(x = doy, y = prop_flower, group = fyear, color = fyear),
              alpha = 0.4) +
    scale_color_discrete(name = "Year") +
    annotate("text", x = 365, y = 0.98, label = "A. chrysanthus",
             hjust = 1, vjust = 1, fontface = 2)
  
# Open flowers
  # Simple model with same seasonal smooth each year
  chryof_m1 <- gam(prop_open ~ s(doy), weights = nobs_open,
                   data = chrydat, method = "REML", family = binomial)
  # Model with an independent seasonal smooth each year
  chryof_m2 <- gam(prop_open ~ s(doy, by = fyear), weights = nobs_open,
                   data = chrydat, method = "REML", family = binomial)
  # Compare models
  summary(chryof_m1)
  summary(chryof_m2)
  AIC(chryof_m1, chryof_m2)
    # Model with different annual curves is better
  
  # Make predictions based on more complex model
  chryof_preds <- cbind(chry_preds,
                        as.data.frame(predict(chryof_m2, 
                                              newdata = chry_preds,
                                              type = "response", 
                                              se.fit = TRUE))) %>%
    mutate(lcl = fit - 1.96 * se.fit,
           ucl = fit + 1.96 * se.fit)
  # Plot predictions
  chryof_plotyr <- ggplot(data = chryof_preds, 
                          aes(x = doy, y = fit, group = fyear, color = fyear)) +
    scale_color_discrete(name = "Year") +
    geom_line() +
    labs(y = "Proportion of flowers that are open",
         x = "Day of year") +
    annotate("text", x = 365, y = 0.98, label = "A. chrysanthus",
             hjust = 1, vjust = 1, fontface = 2)
  
  # Make predictions based on all years together
  chryof_preds1 <- data.frame(doy = 1:365)
  chryof_preds1 <- cbind(chryof_preds1,
                         as.data.frame(predict(chryof_m1, 
                                               newdata = chryof_preds1,
                                               type = "response", 
                                               se.fit = TRUE))) %>%
    mutate(lcl = fit - 1.96 * se.fit,
           ucl = fit + 1.96 * se.fit)
  # Plot predictions (with raw data)
  chryof_plot <- ggplot(chryof_preds1, aes(x = doy)) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), color = "gray", alpha = 0.3) +
    geom_line(aes(y = fit)) +
    labs(y = "Proportion of flowers that are open",
         x = "Day of year") +
    geom_point(data = chrydat, 
               aes(x = doy, y = prop_open, group = fyear, color = fyear),
               alpha = 0.4) +
    geom_line(data = chrydat, 
              aes(x = doy, y = prop_open, group = fyear, color = fyear),
              alpha = 0.4) +
    scale_color_discrete(name = "Year") +
    annotate("text", x = 365, y = 0.98, label = "A. chrysanthus",
             hjust = 1, vjust = 1, fontface = 2)
  
  # Look at plots with seasonal curve for all years
  plot_grid(chryfl_plot, chryof_plot, ncol = 1)
  # Look at plots with independent curves each year
  plot_grid(chryfl_plotyr, chryof_plotyr, ncol = 1)  
  
# All agaves ----------------------------------------------------------------- #
agavedat <- filter(prop_genusyr, spp2 == "agave")
  
# Flowers or buds
  # Simple model with same seasonal smooth each year
  agavefl_m1 <- gam(prop_flower ~ s(doy), weights = nobs_flower,
                   data = agavedat, method = "REML", family = binomial)
  # Model with an independent seasonal smooth each year
  agavefl_m2 <- gam(prop_flower ~ s(doy, by = fyear), weights = nobs_flower,
                   data = agavedat, method = "REML", family = binomial)
  # Compare models
  summary(agavefl_m1)
  summary(agavefl_m2)
  AIC(agavefl_m1, agavefl_m2)
    # Model with different annual curves is better (by a lot)
  
  # Create dataframe for predictions from models with annual curves
  agave_preds <- data.frame(fyear = NA, wks = NA)
  for (yr in sort(unique(agavedat$yr))) {
    min_doy <- min(agavedat$doy[agavedat$yr == yr])
    max_doy <- max(agavedat$doy[agavedat$yr == yr])
    max_doy <- ifelse(max_doy > 350, 365, max_doy)
    doys <- seq(min_doy, max_doy)
    if (yr == min(agavedat$yr)) {
      agave_preds <- data.frame(fyear = yr, doy = doys)
    } else {
      agave_preds <- rbind(agave_preds, 
                          data.frame(fyear = yr, doy = doys))
    }
  }
  agave_preds$fyear <- factor(agave_preds$fyear, levels = as.character(2012:2023))
  
  # Make predictions based on more complex model
  agavefl_preds <- cbind(agave_preds,
                        as.data.frame(predict(agavefl_m2, 
                                              newdata = agave_preds,
                                              type = "response", 
                                              se.fit = TRUE))) %>%
    mutate(lcl = fit - 1.96 * se.fit,
           ucl = fit + 1.96 * se.fit)
  # Plot predictions
  agavefl_plotyr <- ggplot(data = agavefl_preds, 
                          aes(x = doy, y = fit, group = fyear, color = fyear)) +
    scale_color_discrete(name = "Year") +
    geom_line() +
    labs(y = "Proportion of patches with flowers or buds",
         x = "Day of year") +
    annotate("text", x = 365, y = 0.98, label = "Agave spp.",
             hjust = 1, vjust = 1, fontface = 2)
  
  # Make predictions based on simpler model
  agavefl_preds1 <- data.frame(doy = 1:365)
  agavefl_preds1 <- cbind(agavefl_preds1,
                         as.data.frame(predict(agavefl_m1, 
                                               newdata = agavefl_preds1,
                                               type = "response", 
                                               se.fit = TRUE))) %>%
    mutate(lcl = fit - 1.96 * se.fit,
           ucl = fit + 1.96 * se.fit)
  # Plot predictions (with raw data)
  agavefl_plot <- ggplot(agavefl_preds1, aes(x = doy)) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), color = "gray", alpha = 0.3) +
    geom_line(aes(y = fit)) +
    labs(y = "Proportion of patches with flowers or buds",
         x = "Day of year") +
    geom_point(data = agavedat, 
               aes(x = doy, y = prop_flower, group = fyear, color = fyear),
               alpha = 0.4) +
    geom_line(data = agavedat, 
              aes(x = doy, y = prop_flower, group = fyear, color = fyear),
              alpha = 0.4) +
    scale_color_discrete(name = "Year") +
    annotate("text", x = 365, y = 0.98, label = "Agave spp.",
             hjust = 1, vjust = 1, fontface = 2)
  
  # Open flowers
  # Simple model with same seasonal smooth each year
  agaveof_m1 <- gam(prop_open ~ s(doy), weights = nobs_open,
                   data = agavedat, method = "REML", family = binomial)
  # Model with an independent seasonal smooth each year
  agaveof_m2 <- gam(prop_open ~ s(doy, by = fyear), weights = nobs_open,
                   data = agavedat, method = "REML", family = binomial)
  # Compare models
  summary(agaveof_m1)
  summary(agaveof_m2)
  AIC(agaveof_m1, agaveof_m2)
  # Model with different annual curves is better
  
  # Make predictions based on more complex model
  agaveof_preds <- cbind(agave_preds,
                        as.data.frame(predict(agaveof_m2, 
                                              newdata = agave_preds,
                                              type = "response", 
                                              se.fit = TRUE))) %>%
    mutate(lcl = fit - 1.96 * se.fit,
           ucl = fit + 1.96 * se.fit)
  # Plot predictions
  agaveof_plotyr <- ggplot(data = agaveof_preds, 
                          aes(x = doy, y = fit, group = fyear, color = fyear)) +
    scale_color_discrete(name = "Year") +
    geom_line() +
    labs(y = "Proportion of flowers that are open",
         x = "Day of year") +
    annotate("text", x = 365, y = 0.98, label = "Agave spp.",
             hjust = 1, vjust = 1, fontface = 2)
  
  # Make predictions based on all years together
  agaveof_preds1 <- data.frame(doy = 1:365)
  agaveof_preds1 <- cbind(agaveof_preds1,
                         as.data.frame(predict(agaveof_m1, 
                                               newdata = agaveof_preds1,
                                               type = "response", 
                                               se.fit = TRUE))) %>%
    mutate(lcl = fit - 1.96 * se.fit,
           ucl = fit + 1.96 * se.fit)
  # Plot predictions (with raw data)
  agaveof_plot <- ggplot(agaveof_preds1, aes(x = doy)) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), color = "gray", alpha = 0.3) +
    geom_line(aes(y = fit)) +
    labs(y = "Proportion of flowers that are open",
         x = "Day of year") +
    geom_point(data = agavedat, 
               aes(x = doy, y = prop_open, group = fyear, color = fyear),
               alpha = 0.4) +
    geom_line(data = agavedat, 
              aes(x = doy, y = prop_open, group = fyear, color = fyear),
              alpha = 0.4) +
    scale_color_discrete(name = "Year") +
    annotate("text", x = 365, y = 0.98, label = "Agave spp.",
             hjust = 1, vjust = 1, fontface = 2)
  
  # Look at plots with seasonal curve for all years
  plot_grid(agavefl_plot, agaveof_plot, ncol = 1)
  # Look at plots with independent curves each year
  plot_grid(agavefl_plotyr, agaveof_plotyr, ncol = 1)   
  
# Can do some back-of-the-envelope (~inverse) calculations to get 
# estimates of "flowering period" and "open flower period" for each species:
# 1) select a threshold or minimum proportion of plants/patches with flowers 
#    (or open flowers) to designate the period (0.10?)
# 2) calculate the first and last doy when prediction proportion is closest to 
#    that threshold 
# 3) get an estimate of uncertainty for first or last date by finding doys 
#    giving: predicted prob = threshold +/- 1.96 * SE(for pred y @ threshold) 
  
# For saguaros, open flower period 
  # (sagof_m1 with sagof_preds1; sagof_m2 with sagof_preds)
  plot_grid(sagfl_plot, sagof_plot, ncol = 1)
  plot_grid(sagfl_plotyr, sagof_plotyr, ncol = 1)
  
  # Find doy with peak predicted value:
  peak_doy <- sagof_preds1$doy[which.max(sagof_preds1$fit)]
  
  threshold <- 0.1
  sagof_preds1s <- filter(sagof_preds1, doy < peak_doy)
  sagof_preds1e <- filter(sagof_preds1, doy > peak_doy)
  # Find first/last doy where predicted proportion closest to threshold
  doy_f <- sagof_preds1s$doy[which.min(abs(sagof_preds1s$fit - threshold))]
  doy_l <- sagof_preds1e$doy[which.min(abs(sagof_preds1e$fit - threshold))]
  # Find confidence limits around doy estimate for start of period
  ydoy_f <- which(sagof_preds1$doy == doy_f)
  yf_min <- threshold - 1.96 * sagof_preds1$se.fit[ydoy_f]
  yf_max <- threshold + 1.96 * sagof_preds1$se.fit[ydoy_f]
  doy_f_low <- sagof_preds1s$doy[which.min(abs(sagof_preds1s$fit - yf_min))] 
  doy_f_upp <- sagof_preds1s$doy[which.min(abs(sagof_preds1s$fit - yf_max))] 
  paste0("Start: ", doy_f, " (", doy_f_low, ", ", doy_f_upp, ")")
  # Find confidence limits around doy estimate for end of period
  ydoy_l <- which(sagof_preds1$doy == doy_l)
  yl_min <- threshold - 1.96 * sagof_preds1$se.fit[ydoy_l]
  yl_max <- threshold + 1.96 * sagof_preds1$se.fit[ydoy_l]
  doy_l_low <- sagof_preds1e$doy[which.min(abs(sagof_preds1e$fit - yl_max))] 
  doy_l_upp <- sagof_preds1e$doy[which.min(abs(sagof_preds1e$fit - yl_min))] 
  paste0("End: ", doy_l, " (", doy_l_low, ", ", doy_l_upp, ")")
  # Find confidence limits around doy estimate for peak
  yp <- sagof_preds1$lcl[sagof_preds1$doy == peak_doy]
  doy_p_low <- sagof_preds1s$doy[which.min(abs(sagof_preds1s$fit - yp))]
  doy_p_upp <- sagof_preds1e$doy[which.min(abs(sagof_preds1e$fit - yp))]
  paste0("Peak: ", peak_doy, " (", doy_p_low, ", ", doy_p_upp, ")")
  
# Evaluating trends (in start of open flower period, for example):

  # Can use the predicted date for start/end/peakcalculated above from each year 
  # in a simple linear regression (so we'd have sample sizes of 6 or 12, 
  # depending on the species)
  
  # In theory, we could use the first date for individual plants/patches, but the
  # problem is that some individuals weren't observed until after that period 
  # began. Could do some filtering to just use those individuals that were
  # observed on a weekly basis before flowering began or flowers opened, but 
  # I'm not sure how many individuals that would leave us with.
  # Use a better version of ind_yr dataframe?
  ind_yr <- dat %>%
    group_by(ind_id, yr, spp) %>%
    summarize(nobs = length(doy),
              obs_first = min(doy),
              obs_last = max(doy),
              days_since_mn = round(mean(days_since, na.rm = TRUE)),
              nflowers = sum(flowers, na.rm = TRUE),
              nopen = sum(flowers_open, na.rm = TRUE),
              nopen50 = sum(i_flowers_open == "4: 50-74%", na.rm = TRUE),
              nfruit = sum(fruit, na.rm = TRUE),
              nripe = sum(fruit_ripe, na.rm = TRUE),
              flowers_first = ifelse(sum(flowers, na.rm = TRUE) > 0,
                                     min(doy[which(flowers == 1)]), NA),
              open_first = ifelse(sum(flowers_open, na.rm = TRUE) > 0,
                                  min(doy[which(flowers_open == 1)]), NA),
              open50_first = ifelse(sum(i_flowers_open == "4: 50-74%", na.rm = TRUE) > 0,
                                    min(doy[which(i_flowers_open == "4: 50-74%")]), NA),
              flowers_last = ifelse(sum(flowers, na.rm = TRUE) > 0,
                                    max(doy[which(flowers == 1)]), NA),
              open_last = ifelse(sum(flowers_open, na.rm = TRUE) > 0,
                                 max(doy[which(flowers_open == 1)]), NA),
              open50_last = ifelse(sum(i_flowers_open == "4: 50-74%", na.rm = TRUE) > 0,
                                   max(doy[which(i_flowers_open == "4: 50-74%")]), NA),
              .groups = "keep") %>%
    data.frame()
  # Sample sizes (spp*yr)?
  count(filter(ind_yr, nobs > 4), spp)
  # If we require a minimum of 5 obs per year:
    # 34 and 37 for A. parryi, A. chrysanthus
    # 307 and 362 for A. palmeri and C. gigantea
  # Exclude those that weren't monitored until later in the year?
  hist(ind_yr$obs_first[ind_yr$nobs > 4], breaks = 25)
    # Maybe limit to those that have a first observation date before the 
    # estimated peak of the open flower period?
  
  # Test run: start of open flower phenophase in saguaros
  count(filter(ind_yr, spp == "C. gigantea"), yr) # 15-83 plants/yr
  summary_fl <- ind_yr %>%
    group_by(spp, yr) %>%
    summarize(nobs = length(ind_id),
              mn_flowers_first = round(mean(flowers_first, na.rm = TRUE)),
              mn_open_first = round(mean(open_first, na.rm = TRUE)),
              mn_open50_first = round(mean(open50_first, na.rm = TRUE)),
              mn_flowers_last = round(mean(flowers_last, na.rm = TRUE)),
              mn_open_last = round(mean(open_last, na.rm = TRUE)),
              mn_open50_last = round(mean(open50_last, na.rm = TRUE)),
              .groups = "keep") %>%
    data.frame()
    # Sample sizes for A. chrysanthus and parryi really not sufficient to do much
  
  sag_fl <- summary_fl %>%
    filter(spp == "C. gigantea") %>%
    pivot_longer(cols = mn_flowers_first:mn_open50_last, names_to = "event",
                 values_to = "doy") %>%
    mutate(firstlast = ifelse(grepl("first", event), "first", "last"),
           type = ifelse(grepl("50", event), "open50", 
                         ifelse(grepl("flowers", event), "flowers", "open"))) %>%
    data.frame()
  ggplot(data = sag_fl, aes(x = yr, y = doy, group = type, color = type)) +
    geom_point() + 
    geom_line(alpha = 0.2) +
    facet_grid(rows = vars(firstlast))
  # Flowering WAY later in 2013 and 2015 (especially initiation of)
  
  palm_fl <- summary_fl %>%
    filter(spp == "A. palmeri") %>%
    pivot_longer(cols = mn_flowers_first:mn_open50_last, names_to = "event",
                 values_to = "doy") %>%
    mutate(firstlast = ifelse(grepl("first", event), "first", "last"),
           type = ifelse(grepl("50", event), "open50", 
                         ifelse(grepl("flowers", event), "flowers", "open"))) %>%
    data.frame()
  ggplot(data = sag_fl, aes(x = yr, y = doy, group = type, color = type)) +
    geom_point() + 
    geom_line(alpha = 0.2) +
    facet_grid(rows = vars(firstlast))
  # Flowering WAY later in 2013 and 2015 (especially initiation of)
  
  library(lme4)
  m_sag <- lmer(open_first ~ I(yr - 2012) + (1 | ind_id),
                data = filter(ind_yr, spp == "C. gigantea"))
  summary(m_sag)
  m_sag2 <- lmer(open_first ~ I(yr - 2012) + (1 | ind_id),
                data = filter(ind_yr, spp == "C. gigantea", !yr %in% c(2013, 2015)))
  summary(m_sag2)  
  # Compare predictions
  cbind(predict(m_sag, newdata = data.frame(yr = 2012:2023, ind_id = 0), re.form = NA),
        predict(m_sag2, newdata = data.frame(yr = 2012:2023, ind_id = 0), re.form = NA))
  
  
  
  
  
  
# Do all the same calculations for fruit/fruit_ripe?  
# Will have to deal with fewer data (ie, more NAs)
count(dat, flowers); count(dat, flowers_open)
count(dat, fruit); count(dat, fruit_ripe)

# What about intensity data?
count(dat, spp, flowers_open, i_flowers_open)
  # Number of observations with 50-74% of flowers open across all years:
  # A. chrysanthus: 37
  # A. palmeri: 292
  # A. parryi: 5
  # C. gigantea: 100


