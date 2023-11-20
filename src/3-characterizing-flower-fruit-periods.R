# Characterizing flowering/fruiting periods
# Erin Zylstra
# ezylstra@arizona.edu
# 2023-11-19

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

#------------------------------------------------------------------------------#
# Calculate proportions of plants/patches in each phenophase by week and year,
# and use GAMs to describe seasonal patterns
#------------------------------------------------------------------------------#

# Note that we first want to remove multiple observations of the same individual 
# in the same week. Sort so the more advanced phenophase gets kept (if more than 
# one value in a week)
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

# Create a separate dataframe with proportions of all agaves (across species) in 
# each phenophase
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

# For each species (or for all agaves together), we're going to fit GAMs to the
# weekly proportions for a given phenophase. We'll then use the fitted model to
# predict the start, peak, and end day-of-the-year for a phenophase in one year 
# or across all years (after delineating a threshold/minimum proportion). Most 
# efficient to do this in a loop, saving estimates in a table and saving ggplot 
# objects in an output folder.

# Create empty table to hold estimates from GAM models for flowering phases
genusyr <- select(prop_genusyr, spp2, yr) %>%
  distinct() %>%
  filter(spp2 == "agave") %>%
  mutate(spp2 = "Agave spp.") %>%
  rename(taxa = spp2)
taxayr <- select(prop_sppyr, spp, yr) %>%
  distinct() %>%
  rename(taxa = spp) %>%
  rbind(., genusyr)
taxaallyr <- taxayr %>%
  group_by(taxa) %>%
  summarize(yr = paste0(min(yr), "-", max(yr))) %>%
  data.frame()
taxayr <- rbind(taxaallyr, taxayr) %>%
  arrange(taxa)
phase <- c("flowers", "flowers_open")
ests <- taxayr[rep(seq_len(nrow(taxayr)), each = length(phase)), ] %>%
  mutate(phase = rep(phase, nrow(taxayr)),
         start = NA,
         start_lcl = NA,
         start_ucl = NA,
         peak = NA,
         peak_lcl = NA,
         peak_ucl = NA,
         end = NA,
         end_lcl = NA,
         end_ucl = NA)
taxa <- unique(ests$taxa)
alphaline <- 0.3
alphapoly <- 0.4
figw <- 6.5
figh <- 3
# Threshold proportion to define start/end of phenophase
thresh <- 0.1

# for (taxon %in% taxa) {
taxon <- "C. gigantea"

  if (taxon != "Agave spp.") {
    df <- filter(prop_sppyr, spp == taxon)
  } else {
    df <- filter(prop_genusyr, spp2 == "agave")
  }
  
  units <- ifelse(taxon == "C. gigantea", "plants", "patches")
  sppcode <- ifelse(taxon == "C. gigantea", "SAGU", 
                    ifelse(taxon == "A. palmeri", "PALM",
                           ifelse(taxon == "A. parryi", "PARR",
                                  ifelse(taxon == "A. chrysanthus", "CHRY", "AGAV"))))
  
  # Flowers or buds phenophase ------------------------------------------------#
  # Model with same smooth each year
  fl1 <- gam(prop_flower ~ s(doy), weights = nobs_flower,
             data = df, method = "REML", family = binomial)
  # Model with an independent smooth each year
  fl <- gam(prop_flower ~ s(doy, by = fyear), weights = nobs_flower,
            data = df, method = "REML", family = binomial)

  # Make predictions from single smooth model for ggplot object
  fl1_preds <- data.frame(doy = 1:365)
  fl1_preds <- cbind(fl1_preds,
                     as.data.frame(predict(fl1, 
                                           newdata = fl1_preds,
                                           type = "response", 
                                           se.fit = TRUE))) %>%
    mutate(lcl = fit - 1.96 * se.fit,
           ucl = fit + 1.96 * se.fit)
  # Create and save ggplot object
  fl1_plot <- ggplot() +
    geom_point(data = df, aes(x = doy, y = prop_flower, group = fyear, color = fyear),
               alpha = alphaline, size = 1) +
    geom_line(data = df, aes(x = doy, y = prop_flower, group = fyear, color = fyear),
              alpha = alphaline) +
    # scale_color_discrete(name = "Year") +    
    geom_ribbon(data = fl1_preds, aes(x = doy, ymin = lcl, ymax = ucl), 
                color = "gray", alpha = alphapoly) +
    geom_line(data = fl1_preds, aes(x = doy, y = fit)) +
    labs(y = paste0("Proportion of ", units, " with flowers"), 
         x = "Day of year") +
    annotate("text", x = 365, y = 0.98, label = taxon,
             hjust = 1, vjust = 1, fontface = 2) +
    theme(text = element_text(size = 10),
          legend.text = element_text(size = 8))
  ggsave(paste0("output/gams/all-years/", sppcode, "-flowers.png"),
         fl1_plot, device = "png", width = figw, height = figh, 
         units = "in", dpi = 300)
  
  # Create a dataframe to hold predictions from more complex model
  yrpreds <- data.frame(fyear = NA, wks = NA)
  for (yr in sort(unique(df$yr))) {
    min_doy <- min(df$doy[df$yr == yr])
    max_doy <- max(df$doy[df$yr == yr])
    max_doy <- ifelse(max_doy > 350, 365, max_doy)
    doys <- seq(min_doy, max_doy)
    if (yr == min(df$yr)) {
      yrpreds <- data.frame(fyear = yr, doy = doys)
    } else {
      yrpreds <- rbind(yrpreds, data.frame(fyear = yr, doy = doys))
    }
  }
  yrpreds$fyear <- factor(yrpreds$fyear, levels = as.character(2012:2023))
  
  # Make predictions from model with annual smooths for ggplot object
  fl_preds <- cbind(yrpreds,
                    as.data.frame(predict(fl, 
                                          newdata = yrpreds,
                                          type = "response", 
                                          se.fit = TRUE))) %>%
    mutate(lcl = fit - 1.96 * se.fit,
           ucl = fit + 1.96 * se.fit)
  # Create and save ggplot object
  fl_plot <- ggplot(data = fl_preds, 
                    aes(x = doy, y = fit, group = fyear, color = fyear)) +
    # scale_color_discrete(name = "Year") +
    geom_line() +
    labs(y = paste0("Proportion of ", units, " with flowers"), 
         x = "Day of year") +
    annotate("text", x = 365, y = 0.98, label = taxon,
             hjust = 1, vjust = 1, fontface = 2) +
    theme(text = element_text(size = 10),
          legend.text = element_text(size = 8)) +
  ggsave(paste0("output/gams/indiv-years/", sppcode, "-flowers.png"),
         fl_plot, device = "png", width = figw, height = figh, 
         units = "in", dpi = 300)
    
  # Extract predictions (with associated uncertainty) for table
  ri <- which(ests$taxa == taxon & ests$phase == "flowers")[1]
  # Find doy with peak predicted value:
  peak_doy <- fl1_preds$doy[which.max(fl1_preds$fit)]
  ests[ri, "peak"] <- peak_doy
  # Split predictions into before/after peak
  fl1_preds1 <- filter(fl1_preds, doy < peak_doy)
  fl1_preds2 <- filter(fl1_preds, doy > peak_doy)
  # Find first/last doy where predicted proportion closest to threshold
  ests[ri, "start"] <- fl1_preds1$doy[which.min(abs(fl1_preds1$fit - thresh))]
  ests[ri, "end"] <- fl1_preds2$doy[which.min(abs(fl1_preds2$fit - thresh))]
  # Find confidence limits around doy estimate for start of period
  ydoy_f <- which(fl1_preds$doy == ests[ri, "start"])
  yf_min <- thresh - 1.96 * fl1_preds$se.fit[ydoy_f]
  yf_max <- thresh + 1.96 * fl1_preds$se.fit[ydoy_f]
  ests[ri, "start_lcl"] <- min(fl1_preds1$doy[fl1_preds1$fit > yf_min])
  ests[ri, "start_ucl"] <- min(fl1_preds1$doy[fl1_preds1$fit > yf_max])
  # Find confidence limits around doy estimate for end of period
  ydoy_l <- which(fl1_preds$doy == ests[ri, "end"])
  yl_min <- thresh - 1.96 * fl1_preds$se.fit[ydoy_l]
  yl_max <- thresh + 1.96 * fl1_preds$se.fit[ydoy_l]
  ests[ri, "end_lcl"] <- min(fl1_preds2$doy[fl1_preds2$fit < yl_max])
  ests[ri, "end_ucl"] <- min(fl1_preds2$doy[fl1_preds2$fit < yl_min])
  # Find confidence limits around doy estimate for peak
  yp <- fl1_preds$lcl[fl1_preds$doy == peak_doy]
  ests[ri, "peak_lcl"] <- fl1_preds1$doy[which.min(abs(fl1_preds1$fit - yp))]
  ests[ri, "peak_ucl"] <- fl1_preds2$doy[which.min(abs(fl1_preds2$fit - yp))]

  # Do the same for each year using fl_preds...
