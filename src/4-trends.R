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
# Using predicted start/end/peak dates from GAMs in trend analyses
#------------------------------------------------------------------------------#
# Here, we're left with sample sizes of 12 (for saguaros) or fewer

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
    # If we require a minimum of 5 obs per year:
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
  # observations in a year, there's too little to do much with A. parryi 
  # What we can do with flowering data:
    # C. gigantea: 2012-2023 (first dates), 2012-2022 (last dates)
    # A. palmeri: 2018-2023 (first dates), 2018-2022 (last dates)
    # A. chrysanthus?? 2018-2023 (first date), 2018-2022 (last dates)

# Exclude those that weren't monitored until later in the year?
  hist(ind_yr$obs_first[ind_yr$nobs > 4], breaks = 25)
  # Maybe limit to those that have a first observation date before the 
  # estimated peak of the open flower period?
  
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
  ggplot(data = sag_fl, aes(x = yr, y = doy, group = type, color = type)) +
    geom_point() + 
    geom_line(alpha = 0.2) +
    labs(x = "Year", y = "Day of year") +
    scale_color_discrete(name = "C. gigantea") +
    facet_grid(rows = vars(firstlast))
    # Flowering much later in 2013 and 2015 (especially initiation of)

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
  sag_of_first <- lmer(open_first ~ I(yr - 2012) + (1 | ind_id),
                       data = filter(ind_yr, spp == "C. gigantea")) 
  summary(sag_of_first)
  # Does removing 2013 and 2015 change things?  Yes
  sag_of_first2 <- lmer(open_first ~ I(yr - 2012) + (1 | ind_id),
                       data = filter(ind_yr, spp == "C. gigantea", !yr %in% c(2013, 2015))) 
  summary(sag_of_first2)
  # Significant trend towards earlier open flower phenophase is largely driven 
  # by late starts in 2013 and 2015
  
# Open50 flower phenophase (intensity values for open flowers = 50-74%): first date
  sag_of50_first <- lmer(open50_first ~ I(yr - 2012) + (1 | ind_id),
                        data = filter(ind_yr, spp == "C. gigantea")) 
  summary(sag_of50_first)
  # Non-significant trend in first date with 50% of flowers open.
  
# Flower or bud phenophase: last date
  sag_fl_last <- lmer(flowers_last ~ I(yr - 2012) + (1 | ind_id),
                       data = filter(ind_yr, spp == "C. gigantea")) 
  summary(sag_fl_last)
  # Significant trend towards earlier end to flower phenophase
  
# Open flower phenophase: last date
  sag_of_last <- lmer(open_last ~ I(yr - 2012) + (1 | ind_id),
                      data = filter(ind_yr, spp == "C. gigantea")) 
  summary(sag_of_last)
  # Does removing 2013 change things?  Yes
  sag_of_last2 <- lmer(open_last ~ I(yr - 2012) + (1 | ind_id),
                       data = filter(ind_yr, spp == "C. gigantea", !yr == 2013)) 
  summary(sag_of_last2)
  # Significant trend towards earlier end to open flower phenophase driven
  # largely by late end in 2013
  
# Open50 flower phenophase (intensity values for open flowers = 50-74%): last date
  sag_of50_last <- lmer(open50_last ~ I(yr - 2012) + (1 | ind_id),
                        data = filter(ind_yr, spp == "C. gigantea")) 
  summary(sag_of50_last)
  # Non-significant trend in last date with 50% of flowers open.  

# Get predictions?
library(merTools)
PI <- predictInterval(merMod = sag_of_first, 
                      newdata = data.frame(yr = 2013:2023, ind_id = 0),
                      level = 0.95, n.sims = 1000,
                      stat = "mean", type = "linear.prediction",
                      include.resid.var = FALSE)




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
ggplot(data = palm_fl, aes(x = yr, y = doy, group = type, color = type)) +
  geom_point() + 
  geom_line(alpha = 0.2) +
  labs(x = "Year", y = "Day of year") +
  scale_color_discrete(name = "A. palmeri") +
  facet_grid(rows = vars(firstlast))


# A. chrysanthus flowering ----------------------------------------------------#

# Plotting mean first/last dates
chry_fl <- summary_fl %>%
  filter(spp == "A. chrysanthus", yr > 2017) %>%
  pivot_longer(cols = mn_flowers_first:mn_open50_last, names_to = "event",
               values_to = "doy") %>%
  mutate(firstlast = ifelse(grepl("first", event), "Mean first date", "Mean last date"),
         type = ifelse(grepl("50", event), "open50", 
                       ifelse(grepl("flowers", event), "flowers", "open"))) %>%
  data.frame()
ggplot(data = chry_fl, aes(x = yr, y = doy, group = type, color = type)) +
  geom_point() + 
  geom_line(alpha = 0.2) +
  labs(x = "Year", y = "Day of year") +
  scale_color_discrete(name = "A. chrysanthus") +
  facet_grid(rows = vars(firstlast))
  # Note that in several years, first and last date for open50 are the same...

