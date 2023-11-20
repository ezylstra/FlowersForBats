# Assessing trends in flowering/fruiting periods
# Erin Zylstra
# ezylstra@arizona.edu
# 2023-11-19

library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)
library(cowplot)
library(lme4)
library(merTools)

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
  dplyr::select(-state)

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
# Using first/last dates for individual plant/patch data in trend analyses
#------------------------------------------------------------------------------#

# Summarize data for each individual in each year
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

# Number of observations per year, across species
hist(ind_yr$nobs, breaks = 50, xlab = "Number of observations",
     main = "Number of observations of an individual per year")

# Need to do some filtering so we're just using those individuals that were
# observed multiple times per year, starting before the particular phenophase
# began or peaked. 

# Number of observations per year by species
  ind_yr %>%
    group_by(spp) %>%
    summarize(n_indiv = length(unique(ind_id)),
              n_indivyrs = length(spp),
              n_oneobsperyr = sum(nobs == 1),
              n_multobsperyr = sum(nobs > 1),
              n_5obsperyr = sum(nobs > 4)) %>%
    mutate(propmult = round(n_multobsperyr / n_indivyrs, 2),
           prop5 = round(n_5obsperyr / n_indivyrs, 2)) %>%
    data.frame()
    # If we required a minimum of 5 obs per year:
    # 307 and 362 indiv*yr combinations for A. palmeri and C. gigantea
    # 34 and 37 indiv*yr combinations for A. parryi, A. chrysanthus
  
# Number of flowering observations for each species, each year  
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
  summary_fl
  # Even if we include all individuals, regardless of when and how many
  # observations in a year, there's too little to do much with A. parryi or 
  # A. chrysanthus
  # What we can do with flowering data:
    # C. gigantea: 2012-2023 (first dates), 2012-2022 (last dates)
    # A. palmeri: 2018-2023 (first dates), 2018-2022 (last dates)

# Exclude those that weren't monitored until later in the year?
  hist(ind_yr$obs_first[ind_yr$nobs > 4], breaks = 25)
  # Maybe limit to those that have a first observation date before the 
  # estimated peak of the open flower period?
  
# Was going to try and look at trends in the start/end of 50% open intensity
# values, but I'm not sure this makes sense and there's much less data to go on.
  
# Saguaro flowering -----------------------------------------------------------#

# Plotting mean first/last dates
  sag_fl <- summary_fl %>%
    filter(spp == "C. gigantea") %>%
    pivot_longer(cols = mn_flowers_first:mn_open50_last, names_to = "event",
                 values_to = "doy") %>%
    mutate(firstlast = ifelse(grepl("first", event), "Mean first date", "Mean last date"),
           type = ifelse(grepl("50", event), "open50", 
                         ifelse(grepl("flowers", event), "flowers", "open"))) %>%
    data.frame()
  ggplot(data = filter(sag_fl, !grepl("50", event)), 
         aes(x = yr, y = doy, group = type, color = type)) +
    geom_point() + 
    geom_line(alpha = 0.2) +
    labs(x = "Year", y = "Day of year") +
    scale_color_discrete(name = "C. gigantea") +
    facet_grid(rows = vars(firstlast))
    # Flowering much later in 2013 (and to a lesser degree 2015)

# Flower or bud phenophase: first date
  sag_fl_first <- lmer(flowers_first ~ I(yr - 2012) + (1 | ind_id),
                       data = filter(ind_yr, spp == "C. gigantea")) 
  summary(sag_fl_first)
  # Does removing 2013 change things?  Yes
  # (have to remove random effects b/c of singularity issues)
  sag_fl_first2 <- lm(flowers_first ~ I(yr - 2012),
                      data = filter(ind_yr, spp == "C. gigantea", !yr == 2013)) 
  summary(sag_fl_first2)
    # Significant trend towards earlier start to flower phenophase is largely 
    # driven by late start in 2013

# Open flower phenophase: first date
  sag_fo_first <- lmer(open_first ~ I(yr - 2012) + (1 | ind_id),
                       data = filter(ind_yr, spp == "C. gigantea")) 
  summary(sag_fo_first)
  # Does removing 2013 change things?  Yes
  sag_fo_first2 <- lmer(open_first ~ I(yr - 2012) + (1 | ind_id),
                        data = filter(ind_yr, spp == "C. gigantea", !yr == 2013)) 
  summary(sag_fo_first2)
    # Significant trend towards earlier open flower phenophase is largely driven 
    # by late start in 2013

# Flower or bud phenophase: last date
  sag_fl_last <- lmer(flowers_last ~ I(yr - 2012) + (1 | ind_id),
                       data = filter(ind_yr, spp == "C. gigantea", yr < 2023)) 
  summary(sag_fl_last)
  sag_fl_last2 <- lmer(flowers_last ~ I(yr - 2012) + (1 | ind_id),
                       data = filter(ind_yr, spp == "C. gigantea", !yr == 2013, yr < 2023)) 
  summary(sag_fl_last2)
    # Significant trend towards earlier end to flower phenophase largely driven by
    # late end in 2013
  
# Open flower phenophase: last date
  sag_fo_last <- lmer(open_last ~ I(yr - 2012) + (1 | ind_id),
                      data = filter(ind_yr, spp == "C. gigantea", yr < 2023)) 
  summary(sag_fo_last)
  # Does removing 2013 change things?  Yes
  sag_fo_last2 <- lmer(open_last ~ I(yr - 2012) + (1 | ind_id),
                       data = filter(ind_yr, spp == "C. gigantea", !yr == 2013, yr < 2023)) 
  summary(sag_fo_last2)
    # Significant trend towards earlier end to open flower phenophase largely 
    # driven by late end in 2013

# Plot predictions for flower phenophase
  newdat_first <- data.frame(yr = seq(2012, 2023, length = 100), ind_id = 0)
  sag_flf_preds <- predictInterval(merMod = sag_fl_first, 
                                   newdata = newdat_first,
                                   level = 0.95, n.sims = 1000,
                                   stat = "mean", type = "linear.prediction",
                                   include.resid.var = FALSE)
  sag_flf_preds2 <- predict(sag_fl_first2, 
                            newdata = newdat_first,
                            interval = "confidence", level = 0.95,
                            type = "response")
  newdat_last <- data.frame(yr = seq(2012, 2022, length = 80), ind_id = 0)
  sag_fll_preds <- predictInterval(merMod = sag_fl_last, 
                                   newdata = newdat_last,
                                   level = 0.95, n.sims = 1000,
                                   stat = "mean", type = "linear.prediction",
                                   include.resid.var = FALSE)
  sag_fll_preds2 <- predictInterval(merMod = sag_fl_last2, 
                                    newdata = newdat_last,
                                    level = 0.95, n.sims = 1000,
                                    stat = "mean", type = "linear.prediction",
                                    include.resid.var = FALSE)
  sag_flf_preds <- cbind(sag_flf_preds, yr = newdat_first$yr)
  sag_flf_preds2 <- as.data.frame(cbind(sag_flf_preds2, yr = newdat_first$yr))
  sag_fll_preds <- cbind(sag_fll_preds, yr = newdat_last$yr)
  sag_fll_preds2 <- cbind(sag_fll_preds2, yr = newdat_last$yr)
  
  sag_flf_plot <- ggplot() +
    geom_ribbon(data = sag_flf_preds, aes(x = yr, ymin = lwr, ymax = upr), 
                fill = "black", alpha = 0.3, linetype = 0) +
    geom_line(data = sag_flf_preds, aes(x = yr, y = fit)) +
    labs(x = "Year", y = "Mean first flowers day") +
    geom_ribbon(data = sag_flf_preds2, aes(x = yr, ymin = lwr, ymax = upr), 
                fill = "blue", alpha = 0.3, linetype = 0) +
    geom_line(data = sag_flf_preds2, aes(x = yr, y = fit), linetype = 2, color = "blue") +  
    theme(text = element_text(size = 10)) +
    annotate("text", x = 2023, y = 150, label = "C. gigantea - All years", 
             hjust = 1, vjust = 1, fontface = 1) +
    annotate("text", x = 2023, y = 146, label = "C. gigantea - Excluding 2013", 
             hjust = 1, vjust = 1, fontface = 1, color = "blue")
  
  sag_fll_plot <- ggplot() +
    geom_ribbon(data = sag_fll_preds, aes(x = yr, ymin = lwr, ymax = upr), 
                fill = "black", alpha = 0.3, linetype = 0) +
    geom_line(data = sag_fll_preds, aes(x = yr, y = fit)) +
    labs(x = "Year", y = "Mean last flowers day") +
    geom_ribbon(data = sag_fll_preds2, aes(x = yr, ymin = lwr, ymax = upr), 
                fill = "blue", alpha = 0.3, linetype = 0) +
    geom_line(data = sag_fll_preds2, aes(x = yr, y = fit), linetype = 2, color = "blue") +  
    theme(text = element_text(size = 10)) +
    xlim(c(2012, 2023)) +
    annotate("text", x = 2023, y = 200, label = "C. gigantea - All years", 
             hjust = 1, vjust = 1, fontface = 1) +
    annotate("text", x = 2023, y = 195, label = "C. gigantea - Excluding 2013", 
             hjust = 1, vjust = 1, fontface = 1, color = "blue")
  
  plot_grid(sag_flf_plot, sag_fll_plot, ncol = 1)
  
# Plot predictions for open flower phenophase
  sag_fof_preds <- predictInterval(merMod = sag_fo_first, 
                                   newdata = newdat_first,
                                   level = 0.95, n.sims = 1000,
                                   stat = "mean", type = "linear.prediction",
                                   include.resid.var = FALSE)
  sag_fof_preds2 <- predictInterval(merMod = sag_fo_first2, 
                                   newdata = newdat_first,
                                   level = 0.95, n.sims = 1000,
                                   stat = "mean", type = "linear.prediction",
                                   include.resid.var = FALSE)
  sag_fol_preds <- predictInterval(merMod = sag_fo_last, 
                                   newdata = newdat_last,
                                   level = 0.95, n.sims = 1000,
                                   stat = "mean", type = "linear.prediction",
                                   include.resid.var = FALSE)
  sag_fol_preds2 <- predictInterval(merMod = sag_fo_last2, 
                                    newdata = newdat_last,
                                    level = 0.95, n.sims = 1000,
                                    stat = "mean", type = "linear.prediction",
                                    include.resid.var = FALSE)
  sag_fof_preds <- cbind(sag_fof_preds, yr = newdat_first$yr)
  sag_fof_preds2 <- cbind(sag_fof_preds2, yr = newdat_first$yr)
  sag_fol_preds <- cbind(sag_fol_preds, yr = newdat_last$yr)
  sag_fol_preds2 <- cbind(sag_fol_preds2, yr = newdat_last$yr)
  
  sag_fof_plot <- ggplot() +
    geom_ribbon(data = sag_fof_preds, aes(x = yr, ymin = lwr, ymax = upr), 
                fill = "black", alpha = 0.3, linetype = 0) +
    geom_line(data = sag_fof_preds, aes(x = yr, y = fit)) +
    labs(x = "Year", y = "Mean first open flowers day") +
    geom_ribbon(data = sag_fof_preds2, aes(x = yr, ymin = lwr, ymax = upr), 
                fill = "blue", alpha = 0.3, linetype = 0) +
    geom_line(data = sag_fof_preds2, aes(x = yr, y = fit), linetype = 2, color = "blue") +  
    theme(text = element_text(size = 10)) +
    annotate("text", x = 2023, y = 168, label = "C. gigantea - All years", 
             hjust = 1, vjust = 1, fontface = 1) +
    annotate("text", x = 2023, y = 162, label = "C. gigantea - Excluding 2013", 
             hjust = 1, vjust = 1, fontface = 1, color = "blue")
  
  sag_fol_plot <- ggplot() +
    geom_ribbon(data = sag_fol_preds, aes(x = yr, ymin = lwr, ymax = upr), 
                color = "black", alpha = 0.3, linetype = 0) +
    geom_line(data = sag_fol_preds, aes(x = yr, y = fit)) +
    labs(x = "Year", y = "Mean last open flowers day") +
    geom_ribbon(data = sag_fol_preds2, aes(x = yr, ymin = lwr, ymax = upr), 
                fill = "blue", alpha = 0.3, linetype = 0) +
    geom_line(data = sag_fol_preds2, aes(x = yr, y = fit), linetype = 2, color = "blue") +  
    theme(text = element_text(size = 10)) +    
    xlim(c(2012, 2023)) +
    annotate("text", x = 2023, y = 191, label = "C. gigantea - All years", 
             hjust = 1, vjust = 1, fontface = 1) +
    annotate("text", x = 2023, y = 185, label = "C. gigantea - Excluding 2013", 
             hjust = 1, vjust = 1, fontface = 1, color = "blue")

  plot_grid(sag_fof_plot, sag_fol_plot, ncol = 1)

# A. palmeri flowering --------------------------------------------------------#

# Plotting mean first/last dates
palm_fl <- summary_fl %>%
  filter(spp == "A. palmeri", yr > 2017) %>%
  pivot_longer(cols = mn_flowers_first:mn_open50_last, names_to = "event",
               values_to = "doy") %>%
  mutate(firstlast = ifelse(grepl("first", event), "Mean first date", "Mean last date"),
         type = ifelse(grepl("50", event), "open50", 
                       ifelse(grepl("flowers", event), "flowers", "open"))) %>%
  data.frame()
ggplot(data = filter(palm_fl, !grepl("50", event)), 
       aes(x = yr, y = doy, group = type, color = type)) +
  geom_point() + 
  geom_line(alpha = 0.2) +
  labs(x = "Year", y = "Day of year") +
  scale_color_discrete(name = "A. palmeri") +
  facet_grid(rows = vars(firstlast))

# Flower or bud phenophase: first date
palm_fl_first <- lmer(flowers_first ~ I(yr - 2012) + (1 | ind_id),
                      data = filter(ind_yr, spp == "A. palmeri")) 
summary(palm_fl_first)
  # Trend towards earlier start to flower phenophase

# Open flower phenophase: first date
palm_fo_first <- lmer(open_first ~ I(yr - 2012) + (1 | ind_id),
                      data = filter(ind_yr, spp == "A. palmeri")) 
summary(palm_fo_first)
  # Trend towards earlier open flower phenophase

# Flower or bud phenophase: last date
palm_fl_last <- lmer(flowers_last ~ I(yr - 2012) + (1 | ind_id),
                     data = filter(ind_yr, spp == "A. palmeri")) 
summary(palm_fl_last)
  # Slight trend towards earlier end to flower phenophase

# Open flower phenophase: last date
palm_fo_last <- lmer(open_last ~ I(yr - 2012) + (1 | ind_id),
                     data = filter(ind_yr, spp == "A. palmeri")) 
summary(palm_fo_last)
  # No trend in end to open flower phenophase

# Plot predictions for flower phenophase
  newdat_first <- data.frame(yr = seq(2018, 2023, length = 100), ind_id = 0)
  palm_flf_preds <- predictInterval(merMod = palm_fl_first, 
                                   newdata = newdat_first,
                                   level = 0.95, n.sims = 1000,
                                   stat = "mean", type = "linear.prediction",
                                   include.resid.var = FALSE)
  newdat_last <- data.frame(yr = seq(2018, 2022, length = 80), ind_id = 0)
  palm_fll_preds <- predictInterval(merMod = palm_fl_last, 
                                   newdata = newdat_last,
                                   level = 0.95, n.sims = 1000,
                                   stat = "mean", type = "linear.prediction",
                                   include.resid.var = FALSE)
  palm_flf_preds <- cbind(palm_flf_preds, yr = newdat_first$yr) %>%
    mutate(Date = "First")
  palm_fll_preds <- cbind(palm_fll_preds, yr = newdat_last$yr) %>%
    mutate(Date = "Last")
  palm_fl_preds <- rbind(palm_flf_preds, palm_fll_preds)
  
  palm_fl_plot <- ggplot(palm_fl_preds, aes(x = yr)) +
    geom_ribbon(aes(x = yr, ymin = lwr, ymax = upr, group = Date, fill = Date), 
                alpha = 0.3, linetype = 0) +
    geom_line(aes(x = yr, y = fit, group = Date, color = Date)) +
    labs(x = "Year", y = "Mean day of year") +
    theme(text = element_text(size = 10)) +
    scale_color_manual(values = c("forest green", "blue")) +
    scale_fill_manual(values = c("forest green", "blue")) +
    theme_bw() +
    annotate("text", x = 2023, y = 252, label = "A. palmeri, flowers or buds", 
             hjust = 1, vjust = 1, fontface = 2)
  palm_fl_plot

# Plot predictions for open flower phenophase
  palm_fof_preds <- predictInterval(merMod = palm_fo_first, 
                                    newdata = newdat_first,
                                    level = 0.95, n.sims = 1000,
                                    stat = "mean", type = "linear.prediction",
                                    include.resid.var = FALSE)
  palm_fol_preds <- predictInterval(merMod = palm_fo_last, 
                                    newdata = newdat_last,
                                    level = 0.95, n.sims = 1000,
                                    stat = "mean", type = "linear.prediction",
                                    include.resid.var = FALSE)
  palm_fof_preds <- cbind(palm_fof_preds, yr = newdat_first$yr) %>%
    mutate(Date = "First")
  palm_fol_preds <- cbind(palm_fol_preds, yr = newdat_last$yr) %>%
    mutate(Date = "Last")
  palm_fo_preds <- rbind(palm_fof_preds, palm_fol_preds)

palm_fo_plot <- ggplot(palm_fo_preds, aes(x = yr)) +
  geom_ribbon(aes(x = yr, ymin = lwr, ymax = upr, group = Date, fill = Date), 
              alpha = 0.3, linetype = 0) +
  geom_line(aes(x = yr, y = fit, group = Date, color = Date)) +
  labs(x = "Year", y = "Mean day of year") +
  theme(text = element_text(size = 10)) +
  scale_color_manual(values = c("forest green", "blue")) +
  scale_fill_manual(values = c("forest green", "blue")) +
  theme_bw() +
  annotate("text", x = 2023, y = 252, label = "A. palmeri, open flowers", 
           hjust = 1, vjust = 1, fontface = 2)
palm_fo_plot

plot_grid(palm_fl_plot, palm_fo_plot, ncol = 1)


#------------------------------------------------------------------------------#
# Could use predicted start/end/peak dates from GAMs in trend analyses...
#------------------------------------------------------------------------------#
# Samples sizes would be very small (5-12)

