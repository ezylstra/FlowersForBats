# Characterizing flowering/fruiting periods
# Erin Zylstra
# ezylstra@arizona.edu
# 2023-11-22

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
# Figure parameters
#------------------------------------------------------------------------------#

alphaline <- 0.3
alphapoly <- 0.4
figw <- 6.5
figh <- 3

#------------------------------------------------------------------------------#
# Calculate proportions of plants/patches in each flowering phenophase by week 
# and year, and use GAMs to describe seasonal patterns
#------------------------------------------------------------------------------#

# Note that we first want to remove multiple observations of the same individual 
# in the same week. Sort so the more advanced phenophase gets kept (if more than 
# one value in a week)
prop_sppyr_flowers <- dat %>%
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
prop_genusyr_flowers <- dat %>%
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
genusyr <- select(prop_genusyr_flowers, spp2, yr) %>%
  distinct() %>%
  filter(spp2 == "agave") %>%
  mutate(spp2 = "Agave spp.") %>%
  rename(taxa = spp2)
taxayr <- select(prop_sppyr_flowers, spp, yr) %>%
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
ests_flower <- taxayr[rep(seq_len(nrow(taxayr)), each = length(phase)), ] %>%
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
taxa <- unique(ests_flower$taxa)

# Eventually, we'll want to remove annual estimates for a taxa from table when 
# there were < 10 individuals observed that year. Figure out which years meet
# criteria
nobs <- dat %>%
  group_by(spp, yr) %>%
  summarize(nobs = length(unique(ind_id)),
            .groups = "keep") %>%
  rename(taxa = spp) %>%
  data.frame()
nobs2 <- dat %>%
  group_by(spp2, yr) %>%
  summarize(nobs = length(unique(ind_id)),
            .groups = "keep") %>%
  filter(spp2 == "agave") %>%
  rename(taxa = spp2) %>%
  mutate(taxa = "Agave spp.") %>%
  data.frame() 
nobs <- rbind(nobs, nobs2)

# Threshold proportion to define start/end of phenophase
thresholds <- c(0.1, 0.25)

for (thresh in thresholds) {
  
  # Threshold proportion on real scale
  threshl <- log(thresh / (1 - thresh))
  
  for (taxon in taxa) {
  
    if (taxon != "Agave spp.") {
      df <- filter(prop_sppyr_flowers, spp == taxon)
    } else {
      df <- filter(prop_genusyr_flowers, spp2 == "agave")
    }
    
    units <- ifelse(taxon == "C. gigantea", "plants", "patches")
    sppcode <- ifelse(taxon == "C. gigantea", "SAGU", 
                      ifelse(taxon == "A. palmeri", "PALM",
                             ifelse(taxon == "A. parryi", "PARR",
                                    ifelse(taxon == "A. chrysanthus", "CHRY", "AGAV"))))
    yrs <- sort(unique(df$yr))
  
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
                                     type = "link", 
                                     se.fit = TRUE))) %>%
      rename(fitl = fit) %>%
      mutate(lcll = fitl - 1.96 * se.fit,
             ucll = fitl + 1.96 * se.fit,
             fit = exp(fitl) / (1 + exp(fitl)),
             lcl = exp(lcll) / (1 + exp(lcll)),
             ucl = exp(ucll) / (1 + exp(ucll))) %>%
      select(-c(lcll, ucll))
  
    # Create and save ggplot object
    fl1_plot <- ggplot() +
      geom_point(data = df, aes(x = doy, y = prop_flower, group = fyear, color = fyear),
                 alpha = alphaline, size = 1) +
      geom_line(data = df, aes(x = doy, y = prop_flower, group = fyear, color = fyear),
                alpha = alphaline) +
      scale_color_discrete(name = "Year") +    
      geom_ribbon(data = fl1_preds, aes(x = doy, ymin = lcl, ymax = ucl), 
                  color = "gray", alpha = alphapoly, linetype = 0) +
      geom_line(data = fl1_preds, aes(x = doy, y = fit)) +
      labs(y = paste0("Proportion of ", units, " with flowers"), 
           x = "Day of year") +
      annotate("text", x = 365, y = 0.98, label = taxon,
               hjust = 1, vjust = 1, fontface = 2) +
      theme(text = element_text(size = 10),
            legend.text = element_text(size = 8))
    # ggsave(paste0("output/gams/all-years/", sppcode, "-flowers.png"),
    #        fl1_plot, device = "png", width = figw, height = figh, 
    #        units = "in", dpi = 300)
    
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
                                            type = "link", 
                                            se.fit = TRUE))) %>%
      rename(fitl = fit) %>%
      mutate(lcll = fitl - 1.96 * se.fit,
             ucll = fitl + 1.96 * se.fit,
             fit = exp(fitl) / (1 + exp(fitl)),
             lcl = exp(lcll) / (1 + exp(lcll)),
             ucl = exp(ucll) / (1 + exp(ucll))) %>%
      select(-c(lcll, ucll))
    # Create and save ggplot object
    fl_plot <- ggplot(data = fl_preds, 
                      aes(x = doy, y = fit, group = fyear, color = fyear)) +
      scale_color_discrete(name = "Year") +
      geom_line() +
      labs(y = paste0("Proportion of ", units, " with flowers"), 
           x = "Day of year") +
      annotate("text", x = 365, y = 0.98, label = taxon,
               hjust = 1, vjust = 1, fontface = 2) +
      theme(text = element_text(size = 10),
            legend.text = element_text(size = 8))
    # ggsave(paste0("output/gams/indiv-years/", sppcode, "-flowers.png"),
    #        fl_plot, device = "png", width = figw, height = figh, 
    #        units = "in", dpi = 300)
      
    # Extract predictions (with associated uncertainty) for table
    # Find appropriate row of table to insert estimates
    ri <- which(ests_flower$taxa == taxon & ests_flower$phase == "flowers")[1]
    # Find doy with peak predicted value:
    peak_doy <- fl1_preds$doy[which.max(fl1_preds$fit)]
    ests_flower[ri, "peak"] <- peak_doy
    # Split predictions into before/after peak
    fl1_preds1 <- filter(fl1_preds, doy < peak_doy)
    fl1_preds2 <- filter(fl1_preds, doy > peak_doy)
    # Find confidence limits around doy estimate for peak
    yp <- fl1_preds$lcl[fl1_preds$doy == peak_doy]
    ests_flower[ri, "peak_lcl"] <- fl1_preds1$doy[which.min(abs(fl1_preds1$fit - yp))]
    ests_flower[ri, "peak_ucl"] <- fl1_preds2$doy[which.min(abs(fl1_preds2$fit - yp))]
    # Find first doy where predicted proportion closest to threshold, with CI
    ests_flower[ri, "start"] <- fl1_preds1$doy[which.min(abs(fl1_preds1$fit - thresh))]
    ydoy_f <- which(fl1_preds$doy == ests_flower[ri, "start"])
    yf_minl <- threshl - 1.96 * fl1_preds$se.fit[ydoy_f]
    yf_maxl <- threshl + 1.96 * fl1_preds$se.fit[ydoy_f]
    yf_min <- exp(yf_minl) / (1 + exp(yf_minl))
    yf_max <- exp(yf_maxl) / (1 + exp(yf_maxl))
    ests_flower[ri, "start_lcl"] <- min(fl1_preds1$doy[fl1_preds1$fit > yf_min])
    ests_flower[ri, "start_ucl"] <- min(fl1_preds1$doy[fl1_preds1$fit > yf_max])
    # Find last doy where predicted proportion closest to threshold, with CI
    ests_flower[ri, "end"] <- fl1_preds2$doy[which.min(abs(fl1_preds2$fit - thresh))]
    ydoy_l <- which(fl1_preds$doy == ests_flower[ri, "end"])
    yl_minl <- threshl - 1.96 * fl1_preds$se.fit[ydoy_l]
    yl_maxl <- threshl + 1.96 * fl1_preds$se.fit[ydoy_l]
    yl_min <- exp(yl_minl) / (1 + exp(yl_minl))
    yl_max <- exp(yl_maxl) / (1 + exp(yl_maxl))
    ests_flower[ri, "end_lcl"] <- min(fl1_preds2$doy[fl1_preds2$fit < yl_max])
    ests_flower[ri, "end_ucl"] <- min(fl1_preds2$doy[fl1_preds2$fit < yl_min])
  
    # Do the same for each year using fl_preds, except don't calculate CIs 
    # because data in each year are often too sparse
    for (yr in yrs) {
      ri <- which(ests_flower$taxa == taxon & ests_flower$phase == "flowers" & ests_flower$yr == yr)
      fl_preds_yr <- filter(fl_preds, yr == fyear)
      # Find doy with peak predicted value:
      peak_doy <- fl_preds_yr$doy[which.max(fl_preds_yr$fit)]
      ests_flower[ri, "peak"] <- peak_doy
      # Split predictions into before/after peak
      fl_preds_yr1 <- filter(fl_preds_yr, doy < peak_doy)
      fl_preds_yr2 <- filter(fl_preds_yr, doy > peak_doy)
      # Find first doy where predicted proportion closest to threshold
      # (will exclude an early peak if it's a separate event [propr near 0 between])
      ests_flower[ri, "start"] <- max(fl_preds_yr1$doy[fl_preds_yr1$fit < thresh])
      # Find last doy where predicted proportion closest to threshold
      # (will exclude fall peak if it's separate event [props near 0 between])
      if (yr != 2023) {
        ests_flower[ri, "end"] <- min(fl_preds_yr2$doy[fl_preds_yr2$fit < thresh])
      }
    }
  
    # Open flowers phenophase ---------------------------------------------------#
    # Model with same smooth each year
    fo1 <- gam(prop_open ~ s(doy), weights = nobs_open,
               data = df, method = "REML", family = binomial)
    # Model with an independent smooth each year
    fo <- gam(prop_open ~ s(doy, by = fyear), weights = nobs_open,
              data = df, method = "REML", family = binomial)
    
    # Make predictions from single smooth model for ggplot object
    fo1_preds <- data.frame(doy = 1:365)
    fo1_preds <- cbind(fo1_preds,
                       as.data.frame(predict(fo1, 
                                             newdata = fo1_preds,
                                             type = "link", 
                                             se.fit = TRUE))) %>%
      rename(fitl = fit) %>%
      mutate(lcll = fitl - 1.96 * se.fit,
             ucll = fitl + 1.96 * se.fit,
             fit = exp(fitl) / (1 + exp(fitl)),
             lcl = exp(lcll) / (1 + exp(lcll)),
             ucl = exp(ucll) / (1 + exp(ucll))) %>%
      select(-c(lcll, ucll))
    # Create and save ggplot object
    fo1_plot <- ggplot() +
      geom_point(data = df, aes(x = doy, y = prop_open, group = fyear, color = fyear),
                 alpha = alphaline, size = 1) +
      geom_line(data = df, aes(x = doy, y = prop_open, group = fyear, color = fyear),
                alpha = alphaline) +
      scale_color_discrete(name = "Year") +    
      geom_ribbon(data = fo1_preds, aes(x = doy, ymin = lcl, ymax = ucl), 
                  color = "gray", alpha = alphapoly, linetype = 0) +
      geom_line(data = fo1_preds, aes(x = doy, y = fit)) +
      labs(y = paste0("Proportion of ", units, " with open flowers"), 
           x = "Day of year") +
      annotate("text", x = 365, y = 0.98, label = taxon,
               hjust = 1, vjust = 1, fontface = 2) +
      theme(text = element_text(size = 10),
            legend.text = element_text(size = 8))
    # ggsave(paste0("output/gams/all-years/", sppcode, "-flowers-open.png"),
    #        fo1_plot, device = "png", width = figw, height = figh, 
    #        units = "in", dpi = 300)
    
    # Make predictions from model with annual smooths for ggplot object
    fo_preds <- cbind(yrpreds,
                      as.data.frame(predict(fo, 
                                            newdata = yrpreds,
                                            type = "link", 
                                            se.fit = TRUE))) %>%
      rename(fitl = fit) %>%
      mutate(lcll = fitl - 1.96 * se.fit,
             ucll = fitl + 1.96 * se.fit,
             fit = exp(fitl) / (1 + exp(fitl)),
             lcl = exp(lcll) / (1 + exp(lcll)),
             ucl = exp(ucll) / (1 + exp(ucll))) %>%
      select(-c(lcll, ucll))
    # Create and save ggplot object
    fo_plot <- ggplot(data = fo_preds, 
                      aes(x = doy, y = fit, group = fyear, color = fyear)) +
      scale_color_discrete(name = "Year") +
      geom_line() +
      labs(y = paste0("Proportion of ", units, " with open flowers"), 
           x = "Day of year") +
      annotate("text", x = 365, y = 0.98, label = taxon,
               hjust = 1, vjust = 1, fontface = 2) +
      theme(text = element_text(size = 10),
            legend.text = element_text(size = 8))
    # ggsave(paste0("output/gams/indiv-years/", sppcode, "-flowers-open.png"),
    #        fo_plot, device = "png", width = figw, height = figh, 
    #        units = "in", dpi = 300)
    
    # Extract predictions (with associated uncertainty) for table
    # Find appropriate row of table to insert estimates
    ri <- which(ests_flower$taxa == taxon & ests_flower$phase == "flowers_open")[1]
    # Find doy with peak predicted value:
    peak_doy <- fo1_preds$doy[which.max(fo1_preds$fit)]
    ests_flower[ri, "peak"] <- peak_doy
    # Split predictions into before/after peak
    fo1_preds1 <- filter(fo1_preds, doy < peak_doy)
    fo1_preds2 <- filter(fo1_preds, doy > peak_doy)
    # Find confidence limits around doy estimate for peak
    yp <- fo1_preds$lcl[fo1_preds$doy == peak_doy]
    ests_flower[ri, "peak_lcl"] <- fo1_preds1$doy[which.min(abs(fo1_preds1$fit - yp))]
    ests_flower[ri, "peak_ucl"] <- fo1_preds2$doy[which.min(abs(fo1_preds2$fit - yp))]
    # Find first doy where predicted proportion closest to threshold, with CI
    ests_flower[ri, "start"] <- fo1_preds1$doy[which.min(abs(fo1_preds1$fit - thresh))]
    ydoy_f <- which(fo1_preds$doy == ests_flower[ri, "start"])
    yf_minl <- threshl - 1.96 * fo1_preds$se.fit[ydoy_f]
    yf_maxl <- threshl + 1.96 * fo1_preds$se.fit[ydoy_f]
    yf_min <- exp(yf_minl) / (1 + exp(yf_minl))
    yf_max <- exp(yf_maxl) / (1 + exp(yf_maxl))
    ests_flower[ri, "start_lcl"] <- min(fo1_preds1$doy[fo1_preds1$fit > yf_min])
    ests_flower[ri, "start_ucl"] <- min(fo1_preds1$doy[fo1_preds1$fit > yf_max])
    # Find last doy where predicted proportion closest to threshold, with CI
    ests_flower[ri, "end"] <- fo1_preds2$doy[which.min(abs(fo1_preds2$fit - thresh))]
    ydoy_l <- which(fo1_preds$doy == ests_flower[ri, "end"])
    yl_minl <- threshl - 1.96 * fo1_preds$se.fit[ydoy_l]
    yl_maxl <- threshl + 1.96 * fo1_preds$se.fit[ydoy_l]
    yl_min <- exp(yl_minl) / (1 + exp(yl_minl))
    yl_max <- exp(yl_maxl) / (1 + exp(yl_maxl))
    ests_flower[ri, "end_lcl"] <- min(fo1_preds2$doy[fo1_preds2$fit < yl_max])
    ests_flower[ri, "end_ucl"] <- min(fo1_preds2$doy[fo1_preds2$fit < yl_min])
    
    # Do the same for each year using fo_preds, except don't calculate CIs 
    # because data in each year are often too sparse
    for (yr in yrs) {
      ri <- which(ests_flower$taxa == taxon & ests_flower$phase == "flowers_open" & ests_flower$yr == yr)
      fo_preds_yr <- filter(fo_preds, yr == fyear)
      # Find doy with peak predicted value:
      peak_doy <- fo_preds_yr$doy[which.max(fo_preds_yr$fit)]
      ests_flower[ri, "peak"] <- peak_doy
      # Split predictions into before/after peak
      fo_preds_yr1 <- filter(fo_preds_yr, doy < peak_doy)
      fo_preds_yr2 <- filter(fo_preds_yr, doy > peak_doy)
      # Find first doy where predicted proportion closest to threshold
      # (will exclude an early peak if it's a separate event [propr near 0 between])
      ests_flower[ri, "start"] <- max(fo_preds_yr1$doy[fo_preds_yr1$fit < thresh])
      # Find last doy where predicted proportion closest to threshold
      # (will exclude fall peak if it's a separate event [props near 0 between])
      if (yr != 2023) {
        ests_flower[ri, "end"] <- min(fo_preds_yr2$doy[fo_preds_yr2$fit < thresh])
      }
    }
  }
  
  ests_flower_annual <- ests_flower %>%
    filter(!grepl("-", yr)) %>%
    mutate(yr = as.numeric(yr)) %>%
    left_join(nobs, by = c("taxa", "yr")) %>%
    filter(nobs > 9) %>%
    select(!contains("lcl")) %>%
    select(!contains("ucl"))
  ests_flower_annual
  # Save to file
  # file_annual <- paste0("output/gams/estimates-annual-thresh", 
  #                       thresh * 100, ".csv")
  # write.csv(ests_flower_annual,
  #           file_annual,
  #           row.names = FALSE)
  
  ests_flower_allyrs <- ests_flower %>%
    filter(grepl("-", yr))
  ests_flower_allyrs
  # Save to file
  # file_allyrs <- paste0("output/gams/estimates-allyrs-thresh", 
  #                       thresh * 100, ".csv")
  # write.csv(ests_flower_allyrs,
  #           file_allyrs,
  #           row.names = FALSE)
}

#------------------------------------------------------------------------------#
# Calculate proportions of plants/patches in each fruiting phenophase by week 
# and year, and use GAMs to describe seasonal patterns
#------------------------------------------------------------------------------#

# Note that we first want to remove multiple observations of the same individual 
# in the same week. Sort so the more advanced phenophase gets kept (if more than 
# one value in a week)
prop_sppyr_fruit <- dat %>%
  arrange(spp, ind_id, yr, wk, desc(fruit), desc(fruit_ripe)) %>%
  filter(!(is.na(fruit) & is.na(fruit_ripe))) %>%
  distinct(ind_id, yr, wk, .keep_all = TRUE) %>%
  group_by(spp, yr, wk) %>%
  summarize(nobs = length(spp),
            nobs_fruit = sum(!is.na(fruit)),
            prop_fruit = sum(fruit, na.rm = TRUE) / nobs_fruit,
            nobs_ripe = sum(!is.na(fruit_ripe)),
            prop_ripe = sum(fruit_ripe, na.rm = TRUE) / nobs_ripe,
            .groups = "keep") %>%
  mutate(date_generic = parse_date_time(paste(2024, wk, 1, sep="/"), "Y/W/w"),
         date_generic = as.Date(date_generic),
         doy = yday(date_generic),
         fyear = factor(yr, levels = as.character(2012:2023))) %>% 
  data.frame()

# Bats don't eat agave fruit, but going to include A. palmeri just to see what
# that looks like
prop_sppyr_fruit <- prop_sppyr_fruit %>%
  filter(spp %in% c("A. palmeri", "C. gigantea")) %>%
  filter(!(spp == "A. palmeri" & yr == 2017)) 

# For each species, we're going to fit GAMs to the
# weekly proportions for a given phenophase. We'll then use the fitted model to
# predict the start, peak, and end day-of-the-year for a phenophase in one year 
# or across all years (after delineating a threshold/minimum proportion). Most 
# efficient to do this in a loop, saving estimates in a table and saving ggplot 
# objects in an output folder.

# Create empty table to hold estimates from GAM models for fruiting phases
taxayr_fruit <- select(prop_sppyr_fruit, spp, yr) %>%
  distinct() %>%
  rename(taxa = spp)
taxaallyr_fruit <- taxayr_fruit %>%
  group_by(taxa) %>%
  summarize(yr = paste0(min(yr), "-", max(yr))) %>%
  data.frame()
taxayr_fruit <- rbind(taxaallyr_fruit, taxayr_fruit) %>%
  arrange(taxa)
phase <- c("fruit", "fruit_ripe")
ests_fruit <- taxayr_fruit[rep(seq_len(nrow(taxayr_fruit)), each = length(phase)), ] %>%
  mutate(phase = rep(phase, nrow(taxayr_fruit)),
         start = NA,
         start_lcl = NA,
         start_ucl = NA,
         peak = NA,
         peak_lcl = NA,
         peak_ucl = NA,
         end = NA,
         end_lcl = NA,
         end_ucl = NA)
taxa <- unique(ests_fruit$taxa)

for (thresh in thresholds) {
  
  # Threshold proportion on real scale
  threshl <- log(thresh / (1 - thresh))
  
  for (taxon in taxa) {
    
    df <- filter(prop_sppyr_fruit, spp == taxon)

    units <- ifelse(taxon == "C. gigantea", "plants", "patches")
    sppcode <- ifelse(taxon == "C. gigantea", "SAGU", "PALM")
    yrs <- sort(unique(df$yr))
    
    # Fruits phenophase -------------------------------------------------------#
    # Model with same smooth each year
    fr1 <- gam(prop_fruit ~ s(doy), weights = nobs_fruit,
               data = df, method = "REML", family = binomial)
    # Model with an independent smooth each year
    fr <- gam(prop_fruit ~ s(doy, by = fyear), weights = nobs_fruit,
              data = df, method = "REML", family = binomial)
    
    # Make predictions from single smooth model for ggplot object
    fr1_preds <- data.frame(doy = 1:365)
    fr1_preds <- cbind(fr1_preds,
                       as.data.frame(predict(fr1, 
                                             newdata = fr1_preds,
                                             type = "link", 
                                             se.fit = TRUE))) %>%
      rename(fitl = fit) %>%
      mutate(lcll = fitl - 1.96 * se.fit,
             ucll = fitl + 1.96 * se.fit,
             fit = exp(fitl) / (1 + exp(fitl)),
             lcl = exp(lcll) / (1 + exp(lcll)),
             ucl = exp(ucll) / (1 + exp(ucll))) %>%
      select(-c(lcll, ucll))
    
    # Create and save ggplot object
    fr1_plot <- ggplot() +
      geom_point(data = df, aes(x = doy, y = prop_fruit, group = fyear, color = fyear),
                 alpha = alphaline, size = 1) +
      geom_line(data = df, aes(x = doy, y = prop_fruit, group = fyear, color = fyear),
                alpha = alphaline) +
      scale_color_discrete(name = "Year") +    
      geom_ribbon(data = fr1_preds, aes(x = doy, ymin = lcl, ymax = ucl), 
                  color = "gray", alpha = alphapoly, linetype = 0) +
      geom_line(data = fr1_preds, aes(x = doy, y = fit)) +
      labs(y = paste0("Proportion of ", units, " with fruit"), 
           x = "Day of year") +
      annotate("text", x = 365, y = 0.98, label = taxon,
               hjust = 1, vjust = 1, fontface = 2) +
      theme(text = element_text(size = 10),
            legend.text = element_text(size = 8))
    # ggsave(paste0("output/gams/all-years/", sppcode, "-fruit.png"),
    #        fr1_plot, device = "png", width = figw, height = figh,
    #        units = "in", dpi = 300)
    
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
    fr_preds <- cbind(yrpreds,
                      as.data.frame(predict(fr, 
                                            newdata = yrpreds,
                                            type = "link", 
                                            se.fit = TRUE))) %>%
      rename(fitl = fit) %>%
      mutate(lcll = fitl - 1.96 * se.fit,
             ucll = fitl + 1.96 * se.fit,
             fit = exp(fitl) / (1 + exp(fitl)),
             lcl = exp(lcll) / (1 + exp(lcll)),
             ucl = exp(ucll) / (1 + exp(ucll))) %>%
      select(-c(lcll, ucll))
    # Create and save ggplot object
    fr_plot <- ggplot(data = fr_preds, 
                      aes(x = doy, y = fit, group = fyear, color = fyear)) +
      scale_color_discrete(name = "Year") +
      geom_line() +
      labs(y = paste0("Proportion of ", units, " with fruit"), 
           x = "Day of year") +
      annotate("text", x = 365, y = 0.98, label = taxon,
               hjust = 1, vjust = 1, fontface = 2) +
      theme(text = element_text(size = 10),
            legend.text = element_text(size = 8))
    # ggsave(paste0("output/gams/indiv-years/", sppcode, "-fruit.png"),
    #        fr_plot, device = "png", width = figw, height = figh,
    #        units = "in", dpi = 300)
    
    # Extract predictions (with associated uncertainty) for table
    # Find appropriate row of table to insert estimates
    ri <- which(ests_fruit$taxa == taxon & ests_fruit$phase == "fruit")[1]
    # Find doy with peak predicted value:
    peak_doy <- fr1_preds$doy[which.max(fr1_preds$fit)]
    ests_fruit[ri, "peak"] <- peak_doy
    # Split predictions into before/after peak
    fr1_preds1 <- filter(fr1_preds, doy < peak_doy)
    fr1_preds2 <- filter(fr1_preds, doy > peak_doy)
    # Find confidence limits around doy estimate for peak
    yp <- fr1_preds$lcl[fr1_preds$doy == peak_doy]
    ests_fruit[ri, "peak_lcl"] <- fr1_preds1$doy[which.min(abs(fr1_preds1$fit - yp))]
    ests_fruit[ri, "peak_ucl"] <- fr1_preds2$doy[which.min(abs(fr1_preds2$fit - yp))]
    # Find first doy where predicted proportion closest to threshold, with CI
    ests_fruit[ri, "start"] <- fr1_preds1$doy[which.min(abs(fr1_preds1$fit - thresh))]
    ydoy_f <- which(fr1_preds$doy == ests_fruit[ri, "start"])
    yf_minl <- threshl - 1.96 * fr1_preds$se.fit[ydoy_f]
    yf_maxl <- threshl + 1.96 * fr1_preds$se.fit[ydoy_f]
    yf_min <- exp(yf_minl) / (1 + exp(yf_minl))
    yf_max <- exp(yf_maxl) / (1 + exp(yf_maxl))
    ests_fruit[ri, "start_lcl"] <- min(fr1_preds1$doy[fr1_preds1$fit > yf_min])
    ests_fruit[ri, "start_ucl"] <- min(fr1_preds1$doy[fr1_preds1$fit > yf_max])
    # Find last doy where predicted proportion closest to threshold, with CI
    ests_fruit[ri, "end"] <- fr1_preds2$doy[which.min(abs(fr1_preds2$fit - thresh))]
    ydoy_l <- which(fr1_preds$doy == ests_fruit[ri, "end"])
    yl_minl <- threshl - 1.96 * fr1_preds$se.fit[ydoy_l]
    yl_maxl <- threshl + 1.96 * fr1_preds$se.fit[ydoy_l]
    yl_min <- exp(yl_minl) / (1 + exp(yl_minl))
    yl_max <- exp(yl_maxl) / (1 + exp(yl_maxl))
    ests_fruit[ri, "end_lcl"] <- min(fr1_preds2$doy[fr1_preds2$fit < yl_max])
    ests_fruit[ri, "end_ucl"] <- min(fr1_preds2$doy[fr1_preds2$fit < yl_min])
    
    # Do the same for each year using fl_preds, except don't calculate CIs 
    # because data in each year are often too sparse
    for (yr in yrs) {
      ri <- which(ests_fruit$taxa == taxon & ests_fruit$phase == "fruit" & ests_fruit$yr == yr)
      fr_preds_yr <- filter(fr_preds, yr == fyear)
      # Find doy with peak predicted value:
      peak_doy <- fr_preds_yr$doy[which.max(fr_preds_yr$fit)]
      ests_fruit[ri, "peak"] <- peak_doy
      # Split predictions into before/after peak
      fr_preds_yr1 <- filter(fr_preds_yr, doy < peak_doy)
      fr_preds_yr2 <- filter(fr_preds_yr, doy > peak_doy)
      # Find first doy where predicted proportion closest to threshold
      # (will exclude an early peak if it's a separate event [propr near 0 between])
      ests_fruit[ri, "start"] <- max(fr_preds_yr1$doy[fr_preds_yr1$fit < thresh])
      # Find last doy where predicted proportion closest to threshold
      # (will exclude fall peak if it's separate event [props near 0 between])
      if (yr != 2023) {
        ests_fruit[ri, "end"] <- min(fr_preds_yr2$doy[fr_preds_yr2$fit < thresh])
      }
    }
    
    # Ripe fruit phenophase ---------------------------------------------------#
    # Model with same smooth each year
    ri1 <- gam(prop_ripe ~ s(doy), weights = nobs_ripe,
               data = df, method = "REML", family = binomial)
    # Model with an independent smooth each year
    ri <- gam(prop_ripe ~ s(doy, by = fyear), weights = nobs_ripe,
              data = df, method = "REML", family = binomial)
    
    # Make predictions from single smooth model for ggplot object
    ri1_preds <- data.frame(doy = 1:365)
    ri1_preds <- cbind(ri1_preds,
                       as.data.frame(predict(ri1, 
                                             newdata = ri1_preds,
                                             type = "link", 
                                             se.fit = TRUE))) %>%
      rename(fitl = fit) %>%
      mutate(lcll = fitl - 1.96 * se.fit,
             ucll = fitl + 1.96 * se.fit,
             fit = exp(fitl) / (1 + exp(fitl)),
             lcl = exp(lcll) / (1 + exp(lcll)),
             ucl = exp(ucll) / (1 + exp(ucll))) %>%
      select(-c(lcll, ucll))
    # Create and save ggplot object
    ri1_plot <- ggplot() +
      geom_point(data = df, aes(x = doy, y = prop_ripe, group = fyear, color = fyear),
                 alpha = alphaline, size = 1) +
      geom_line(data = df, aes(x = doy, y = prop_ripe, group = fyear, color = fyear),
                alpha = alphaline) +
      scale_color_discrete(name = "Year") +    
      geom_ribbon(data = ri1_preds, aes(x = doy, ymin = lcl, ymax = ucl), 
                  color = "gray", alpha = alphapoly, linetype = 0) +
      geom_line(data = ri1_preds, aes(x = doy, y = fit)) +
      labs(y = paste0("Proportion of ", units, " with ripe fruit"), 
           x = "Day of year") +
      annotate("text", x = 365, y = 0.98, label = taxon,
               hjust = 1, vjust = 1, fontface = 2) +
      theme(text = element_text(size = 10),
            legend.text = element_text(size = 8))
    # ggsave(paste0("output/gams/all-years/", sppcode, "-fruit-ripe.png"),
    #        ri1_plot, device = "png", width = figw, height = figh,
    #        units = "in", dpi = 300)
    
    # Make predictions from model with annual smooths for ggplot object
    ri_preds <- cbind(yrpreds,
                      as.data.frame(predict(ri, 
                                            newdata = yrpreds,
                                            type = "link", 
                                            se.fit = TRUE))) %>%
      rename(fitl = fit) %>%
      mutate(lcll = fitl - 1.96 * se.fit,
             ucll = fitl + 1.96 * se.fit,
             fit = exp(fitl) / (1 + exp(fitl)),
             lcl = exp(lcll) / (1 + exp(lcll)),
             ucl = exp(ucll) / (1 + exp(ucll))) %>%
      select(-c(lcll, ucll))
    # Create and save ggplot object
    ri_plot <- ggplot(data = ri_preds, 
                      aes(x = doy, y = fit, group = fyear, color = fyear)) +
      scale_color_discrete(name = "Year") +
      geom_line() +
      labs(y = paste0("Proportion of ", units, " with ripe fruit"), 
           x = "Day of year") +
      annotate("text", x = 365, y = 0.98, label = taxon,
               hjust = 1, vjust = 1, fontface = 2) +
      theme(text = element_text(size = 10),
            legend.text = element_text(size = 8))
    # ggsave(paste0("output/gams/indiv-years/", sppcode, "-fruit-ripe.png"),
    #        ri_plot, device = "png", width = figw, height = figh,
    #        units = "in", dpi = 300)
    
    # Extract predictions (with associated uncertainty) for table
    # Find appropriate row of table to insert estimates
    rin <- which(ests_fruit$taxa == taxon & ests_fruit$phase == "fruit_ripe")[1]
    # Find doy with peak predicted value:
    peak_doy <- ri1_preds$doy[which.max(ri1_preds$fit)]
    ests_fruit[rin, "peak"] <- peak_doy
    # Split predictions into before/after peak
    ri1_preds1 <- filter(ri1_preds, doy < peak_doy)
    ri1_preds2 <- filter(ri1_preds, doy > peak_doy)
    # Find confidence limits around doy estimate for peak
    yp <- ri1_preds$lcl[ri1_preds$doy == peak_doy]
    ests_fruit[rin, "peak_lcl"] <- ri1_preds1$doy[which.min(abs(ri1_preds1$fit - yp))]
    ests_fruit[rin, "peak_ucl"] <- ri1_preds2$doy[which.min(abs(ri1_preds2$fit - yp))]
    # Find first doy where predicted proportion closest to threshold, with CI
    ests_fruit[rin, "start"] <- ri1_preds1$doy[which.min(abs(ri1_preds1$fit - thresh))]
    ydoy_f <- which(ri1_preds$doy == ests_fruit[rin, "start"])
    yf_minl <- threshl - 1.96 * ri1_preds$se.fit[ydoy_f]
    yf_maxl <- threshl + 1.96 * ri1_preds$se.fit[ydoy_f]
    yf_min <- exp(yf_minl) / (1 + exp(yf_minl))
    yf_max <- exp(yf_maxl) / (1 + exp(yf_maxl))
    ests_fruit[rin, "start_lcl"] <- min(ri1_preds1$doy[ri1_preds1$fit > yf_min])
    ests_fruit[rin, "start_ucl"] <- min(ri1_preds1$doy[ri1_preds1$fit > yf_max])
    # Find last doy where predicted proportion closest to threshold, with CI
    ests_fruit[rin, "end"] <- ri1_preds2$doy[which.min(abs(ri1_preds2$fit - thresh))]
    ydoy_l <- which(ri1_preds$doy == ests_fruit[rin, "end"])
    yl_minl <- threshl - 1.96 * ri1_preds$se.fit[ydoy_l]
    yl_maxl <- threshl + 1.96 * ri1_preds$se.fit[ydoy_l]
    yl_min <- exp(yl_minl) / (1 + exp(yl_minl))
    yl_max <- exp(yl_maxl) / (1 + exp(yl_maxl))
    ests_fruit[rin, "end_lcl"] <- min(ri1_preds2$doy[ri1_preds2$fit < yl_max])
    ests_fruit[rin, "end_ucl"] <- min(ri1_preds2$doy[ri1_preds2$fit < yl_min])
    
    # Do the same for each year using fo_preds, except don't calculate CIs 
    # because data in each year are often too sparse
    for (yr in yrs) {
      rin <- which(ests_fruit$taxa == taxon & ests_fruit$phase == "fruit_ripe" & ests_fruit$yr == yr)
      ri_preds_yr <- filter(ri_preds, yr == fyear)
      # Find doy with peak predicted value:
      peak_doy <- ri_preds_yr$doy[which.max(ri_preds_yr$fit)]
      ests_fruit[rin, "peak"] <- peak_doy
      # Split predictions into before/after peak
      ri_preds_yr1 <- filter(ri_preds_yr, doy < peak_doy)
      ri_preds_yr2 <- filter(ri_preds_yr, doy > peak_doy)
      # Find first doy where predicted proportion closest to threshold
      # (will exclude an early peak if it's a separate event [propr near 0 between])
      ests_fruit[rin, "start"] <- max(ri_preds_yr1$doy[ri_preds_yr1$fit < thresh])
      # Find last doy where predicted proportion closest to threshold
      # (will exclude fall peak if it's a separate event [props near 0 between])
      if (yr != 2023) {
        ests_fruit[rin, "end"] <- min(ri_preds_yr2$doy[ri_preds_yr2$fit < thresh])
      }
    }
  }
  
  ests_fruit_annual <- ests_fruit %>%
    filter(!grepl("-", yr)) %>%
    mutate(yr = as.numeric(yr)) %>%
    left_join(nobs, by = c("taxa", "yr")) %>%
    filter(nobs > 9) %>%
    select(!contains("lcl")) %>%
    select(!contains("ucl"))
  ests_fruit_annual
  # Save to file
  file_annual <- paste0("output/gams/estimates-fruit-annual-thresh", 
                        thresh * 100, ".csv")
  # write.csv(ests_fruit_annual,
  #           file_annual,
  #           row.names = FALSE)
  
  ests_fruit_allyrs <- ests_fruit %>%
    filter(grepl("-", yr))
  ests_fruit_allyrs
  # Save to file
  file_allyrs <- paste0("output/gams/estimates-fruit-allyrs-thresh", 
                        thresh * 100, ".csv")
  # write.csv(ests_fruit_allyrs,
  #           file_allyrs,
  #           row.names = FALSE)
}

# Fruiting data for A. palmeri are wacky -- fruit (and even ripe fruit) observed
# from late July through March...

#------------------------------------------------------------------------------#
# Create figures with estimates
#------------------------------------------------------------------------------#
flowers_allyrs_10 <- read.csv("output/gams/estimates-allyrs-thresh10.csv")
flowers_allyrs_25 <- read.csv("output/gams/estimates-allyrs-thresh25.csv")
fruit_allyrs_10 <- read.csv("output/gams/estimates-fruit-allyrs-thresh10.csv")
fruit_allyrs_25 <- read.csv("output/gams/estimates-fruit-allyrs-thresh25.csv")

# Plot onset, peak, end of flowering, fruiting for C. gigantea and A. palmeri
fl_10 <- flowers_allyrs_10 %>%
  dplyr::select(!contains("ucl")) %>%
  dplyr::select(!contains("lcl")) %>%
  filter(taxa %in% c("C. gigantea", "A. palmeri")) %>%
  mutate(taxa = factor(taxa, levels = c("A. palmeri", "C. gigantea")),
         taxa_num = as.numeric(taxa),
         threshold = "0.10") %>%
  data.frame()
fl_25 <- flowers_allyrs_25 %>%
  dplyr::select(!contains("ucl")) %>%
  dplyr::select(!contains("lcl")) %>%
  filter(taxa %in% c("C. gigantea", "A. palmeri")) %>%
  mutate(taxa = factor(taxa, levels = c("A. palmeri", "C. gigantea")),
         taxa_num = as.numeric(taxa),
         threshold = "0.25") %>%
  data.frame()
fr_10 <- fruit_allyrs_10 %>%
  dplyr::select(!contains("ucl")) %>%
  dplyr::select(!contains("lcl")) %>%
  filter(taxa == "C. gigantea") %>%
  mutate(taxa = factor(taxa),
         taxa_num = as.numeric(taxa),
         threshold = "0.10") %>%
  data.frame()
fr_25 <- fruit_allyrs_25 %>%
  dplyr::select(!contains("ucl")) %>%
  dplyr::select(!contains("lcl")) %>%
  filter(taxa == "C. gigantea") %>%
  mutate(taxa = factor(taxa),
         taxa_num = as.numeric(taxa),
         threshold = "0.25") %>%
  data.frame()

dfe <- rbind(fl_10, fl_25, fr_10, fr_25)

cap <- 0.2

startends <- ggplot(dfe) +
  geom_segment(aes(x = start, xend = end, y = taxa_num, yend = taxa_num, 
                   group = threshold, color = threshold)) +
  geom_segment(aes(x = start, xend = start, 
                   y = taxa_num - cap, yend = taxa_num + cap, 
                   group = threshold, color = threshold)) +
  geom_segment(aes(x = peak, xend = peak, 
                   y = taxa_num - cap, yend = taxa_num + cap, 
                   group = threshold, color = threshold), linewidth = 1.2) +
  geom_segment(aes(x = end, xend = end, 
                   y = taxa_num - cap, yend = taxa_num + cap, 
                   group = threshold, color = threshold)) +
  scale_color_manual(values = c("gray60", "dodgerblue")) +
  facet_grid(rows = "phase") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0.5, 2.5), breaks = 1:2, labels = levels(dfe$taxa)) +
  labs(x = "Day of year", y = "", color = "Proportion")
ggsave("output/gams/all-years/flowerfruit-start-ends.png", startends, device = "png", 
       width = 6.5, height = 4, units = "in", dpi = 300)
