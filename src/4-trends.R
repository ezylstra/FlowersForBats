# Assessing trends in flowering/fruiting periods
# Erin Zylstra
# ezylstra@arizona.edu
# 2023-11-20

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

# Exclude those that weren't monitored early or late in the year (depending
# on the question)?
  hist(ind_yr$obs_first[ind_yr$nobs > 4], breaks = 25)
  hist(ind_yr$obs_last[ind_yr$nobs > 4], breaks = 25)
  hist(ind_yr$obs_first[ind_yr$nobs > 2], breaks = 25)
  hist(ind_yr$obs_last[ind_yr$nobs > 2], breaks = 25)

# Flowers period (based on GAMs)
  # C. gigantea average first/last date with > 10% individuals = 79/185 
  # C. gigantea average date with highest proportion of individuals = 129
  # A. palmeri average first/last date with > 10% individuals = 128/278
  # A. palmeri average date with highest proportion of individuals = 187
  
# Open flowers period (based on GAMs)
  # C. gigantea average first/last date with > 10% individuals = 104/174 
  # C. gigantea average date with highest proportion of individuals = 135
  # A. palmeri average first/last date with > 10% individuals = 169/269 
  # A. palmeri average date with highest proportion of individuals = 205

###############################################################################
  
# Saguaro flowering -----------------------------------------------------------#
# Use individuals that were observed at least 3 times in year
sagdat <- filter(ind_yr, spp == "C. gigantea", nobs > 2)

# Plotting mean first/last dates
summary_sag_flf <- sagdat %>%
  # Just individuals observed before open flower peak
  filter(obs_first < 135) %>%
  group_by(yr) %>%
  summarize(nobs_firstdates = length(ind_id),
            mn_flowers_first = round(mean(flowers_first, na.rm = TRUE)),
            mn_open_first = round(mean(open_first, na.rm = TRUE)),
            mn_open50_first = round(mean(open50_first, na.rm = TRUE))) %>%
  data.frame()
summary_sag_fll <- sagdat %>%
  # Just individuals observed after open flower peak
  filter(obs_last > 135) %>%
  group_by(yr) %>%
  summarize(nobs_lastdates = length(ind_id),
            mn_flowers_last = round(mean(flowers_last, na.rm = TRUE)),
            mn_open_last = round(mean(open_last, na.rm = TRUE)),
            mn_open50_last = round(mean(open50_last, na.rm = TRUE))) %>%
  data.frame()
summary_sag_fl <- left_join(summary_sag_flf, summary_sag_fll, by = "yr") %>%
  relocate(nobs_lastdates, .after = "nobs_firstdates") %>%
  pivot_longer(cols = mn_flowers_first:mn_open50_last, names_to = "event",
               values_to = "doy") %>%
  mutate(firstlast = ifelse(grepl("first", event), "Mean first date", "Mean last date"),
         type = ifelse(grepl("50", event), "open50", 
                       ifelse(grepl("flowers", event), "flowers", "open"))) %>%
  data.frame()
ggplot(data = summary_sag_fl, 
       aes(x = yr, y = doy, group = type, color = type)) +
  geom_point() + 
  geom_line(alpha = 0.2) +
  labs(x = "Year", y = "Day of year") +
  scale_color_discrete(name = "C. gigantea") +
  facet_grid(rows = vars(firstlast))

# Flower or bud phenophase: first date
  sag_fl_first <- lmer(flowers_first ~ I(yr - 2012) + (1 | ind_id), 
                       data = filter(sagdat, obs_first < 135)) 
  summary(sag_fl_first)
  confint(sag_fl_first)
  # Some evidence of trend earlier

# Open flower phenophase: first date
  sag_fo_first <- lmer(open_first ~ I(yr - 2012) + (1 | ind_id), 
                       data = filter(sagdat, obs_first < 135)) 
  summary(sag_fo_first)
  confint(sag_fo_first)
  # Some evidence of trend earlier

# Flower or bud phenophase: last date
  sag_fl_last <- lmer(flowers_last ~ I(yr - 2012) + (1 | ind_id),
                       data = filter(sagdat, obs_last > 135, yr < 2023)) 
  summary(sag_fl_last)
  confint(sag_fl_last)
  # Some evidence of trend earlier
  
# Open flower phenophase: last date
  sag_fo_last <- lmer(open_last ~ I(yr - 2012) + (1 | ind_id),
                      data = filter(sagdat, obs_last > 135, yr < 2023)) 
  summary(sag_fo_last)
  confint(sag_fo_last)
  # No evidence of trend

# Predictions: flowers or buds
  newdat_first <- data.frame(yr = seq(2012, 2023, length = 100), ind_id = 0)
  sag_flf_preds <- predictInterval(merMod = sag_fl_first, 
                                   newdata = newdat_first,
                                   level = 0.95, n.sims = 1000,
                                   stat = "mean", type = "linear.prediction",
                                   include.resid.var = FALSE)
  newdat_last <- data.frame(yr = seq(2012, 2022, length = 80), ind_id = 0)
  sag_fll_preds <- predictInterval(merMod = sag_fl_last, 
                                   newdata = newdat_last,
                                   level = 0.95, n.sims = 1000,
                                   stat = "mean", type = "linear.prediction",
                                   include.resid.var = FALSE)
  sag_flf_preds <- cbind(sag_flf_preds, yr = newdat_first$yr) %>%
    mutate(Date = "First")
  sag_fll_preds <- cbind(sag_fll_preds, yr = newdat_last$yr) %>%
    mutate(Date = "Last")
  sag_fl_preds <- rbind(sag_flf_preds, sag_fll_preds)

# Predictions: open flowers
  sag_fof_preds <- predictInterval(merMod = sag_fo_first, 
                                   newdata = newdat_first,
                                   level = 0.95, n.sims = 1000,
                                   stat = "mean", type = "linear.prediction",
                                   include.resid.var = FALSE)
  sag_fol_preds <- predictInterval(merMod = sag_fo_last, 
                                   newdata = newdat_last,
                                   level = 0.95, n.sims = 1000,
                                   stat = "mean", type = "linear.prediction",
                                   include.resid.var = FALSE)
  sag_fof_preds <- cbind(sag_fof_preds, yr = newdat_first$yr) %>%
    mutate(Date = "First")
  sag_fol_preds <- cbind(sag_fol_preds, yr = newdat_last$yr) %>%
    mutate(Date = "Last")
  sag_fo_preds <- rbind(sag_fof_preds, sag_fol_preds)

# Plot predictions
  sag_fl_plot <- ggplot(sag_fl_preds, aes(x = yr)) +
    geom_ribbon(aes(x = yr, ymin = lwr, ymax = upr, group = Date, fill = Date), 
                alpha = 0.3, linetype = 0) +
    geom_line(aes(x = yr, y = fit, group = Date, color = Date)) +
    labs(x = "Year", y = "Mean day of year") +
    scale_color_manual(values = c("forest green", "blue")) +
    scale_fill_manual(values = c("forest green", "blue")) +
    theme_bw() +
    annotate("text", x = 2023, y = 212, label = "C. gigantea, flowers or buds", 
             hjust = 1, vjust = 1, fontface = 2, size = 3) +
    theme(text = element_text(size = 10),
          legend.text = element_text(size = 8))
  sag_fo_plot <- ggplot(sag_fo_preds, aes(x = yr)) +
    geom_ribbon(aes(x = yr, ymin = lwr, ymax = upr, group = Date, fill = Date), 
                alpha = 0.3, linetype = 0) +
    geom_line(aes(x = yr, y = fit, group = Date, color = Date)) +
    labs(x = "Year", y = "Mean day of year") +
    scale_color_manual(values = c("forest green", "blue")) +
    scale_fill_manual(values = c("forest green", "blue")) +
    theme_bw() +
    annotate("text", x = 2023, y = 195, label = "C. gigantea, open flowers", 
             hjust = 1, vjust = 1, fontface = 2, size = 3) +
    theme(text = element_text(size = 10),
          legend.text = element_text(size = 8))
  sag_trends <- plot_grid(sag_fl_plot, sag_fo_plot, ncol = 1)
  ggsave("output/trends/SAGU-flowers.png",
         sag_trends, device = "png", width = 6.5, height = 5, 
         units = "in", dpi = 300)

# Flowering intensity: 50-74% of flowers open 
# (using first date since last dates are often the same or very similar)
  # Probably don't need to filter individuals based on number of observations.
  sagdat2 <- filter(ind_yr, spp == "C. gigantea")
  count(filter(sagdat2, !is.na(open50_first)), yr)
  # There are very few observations in some years. Eliminating 2012-2013, 
  # where we only have one individual
  sagdat2 <- filter(sagdat2, yr > 2013)
  sag_fo50_first <- lmer(open50_first ~ I(yr - 2014) + (1 | ind_id), 
                         data = sagdat2) 
  summary(sag_fo50_first)
  # No evidence of trend

  # Annual and overall estimates
  sag_fo50_yr <- lm(open50_first ~ factor(yr) - 1, data = sagdat2)
  summary(sag_fo50_yr)
  confint(sag_fo50_yr)
  sag_fo50_1 <- lm(open50_first ~ 1, data = sagdat2)
  summary(sag_fo50_1)
  confint(sag_fo50_1)
  # Overall mean DOY = 150 (SE = 3.9; 95% CI = 142-158)
  AIC(sag_fo50_yr); AIC(sag_fo50_1) # Year model a little better
  
  # Table:
  sag50 <- sagdat2 %>%
    filter(!is.na(open50_first)) %>%
    group_by(yr) %>%
    summarize(no_indiv = length(yr)) %>%
    mutate(yr = as.character(yr)) %>%
    rbind(data.frame(yr = "2014-2023", 
                     no_indiv = sum(!is.na(sagdat2$open50_first)))) %>%
    mutate(est = c(coef(sag_fo50_yr), coef(sag_fo50_1)),
           lcl = c(confint(sag_fo50_yr)[,1], confint(sag_fo50_1)[,1]),
           ucl = c(confint(sag_fo50_yr)[,2], confint(sag_fo50_1)[,2])) %>%
    mutate(spp = "C. gigantea", .before = yr) %>%
    data.frame()
  # write.table(sag50, "clipboard", sep = "\t", row.names = FALSE)
  
# A. palmeri flowering --------------------------------------------------------#
# Use individuals that were observed at least 3 times in year (and removing the
# one indvidual observed in 2017)
palmdat <- filter(ind_yr, spp == "A. palmeri", nobs > 2, yr > 2017)

# Plotting mean first/last dates
summary_palm_flf <- palmdat %>%
  # Just individuals observed before open flower peak
  filter(obs_first < 205) %>%
  group_by(yr) %>%
  summarize(nobs_firstdates = length(ind_id),
            mn_flowers_first = round(mean(flowers_first, na.rm = TRUE)),
            mn_open_first = round(mean(open_first, na.rm = TRUE)),
            mn_open50_first = round(mean(open50_first, na.rm = TRUE))) %>%
  data.frame()
summary_palm_fll <- palmdat %>%
  # Just individuals observed after open flower peak
  filter(obs_last > 205) %>%
  group_by(yr) %>%
  summarize(nobs_lastdates = length(ind_id),
            mn_flowers_last = round(mean(flowers_last, na.rm = TRUE)),
            mn_open_last = round(mean(open_last, na.rm = TRUE)),
            mn_open50_last = round(mean(open50_last, na.rm = TRUE))) %>%
  data.frame()
summary_palm_fl <- left_join(summary_palm_flf, summary_palm_fll, by = "yr") %>%
  relocate(nobs_lastdates, .after = "nobs_firstdates") %>%
  pivot_longer(cols = mn_flowers_first:mn_open50_last, names_to = "event",
               values_to = "doy") %>%
  mutate(firstlast = ifelse(grepl("first", event), "Mean first date", "Mean last date"),
         type = ifelse(grepl("50", event), "open50", 
                       ifelse(grepl("flowers", event), "flowers", "open"))) %>%
  data.frame()
ggplot(data = summary_palm_fl, 
       aes(x = yr, y = doy, group = type, color = type)) +
  geom_point() + 
  geom_line(alpha = 0.2) +
  labs(x = "Year", y = "Day of year") +
  scale_color_discrete(name = "A. palmeri") +
  facet_grid(rows = vars(firstlast))
  
# Flower or bud phenophase: first date
  palm_fl_first <- lmer(flowers_first ~ I(yr - 2012) + (1 | ind_id), 
                       data = filter(palmdat, obs_first < 205)) 
  summary(palm_fl_first)
  confint(palm_fl_first)
  # Trending earlier

# Open flower phenophase: first date
  palm_fo_first <- lmer(open_first ~ I(yr - 2012) + (1 | ind_id), 
                       data = filter(palmdat, obs_first < 205)) 
  summary(palm_fo_first)
  confint(palm_fo_first)
  # Some evidence of trend earlier (smaller effect)
  
# Flower or bud phenophase: last date
  palm_fl_last <- lmer(flowers_last ~ I(yr - 2012) + (1 | ind_id),
                      data = filter(palmdat, obs_last > 205, yr < 2023)) 
  summary(palm_fl_last)
  confint(palm_fl_last)
  # Evidence of trend for later dates
  
# Open flower phenophase: last date
  palm_fo_last <- lmer(open_last ~ I(yr - 2012) + (1 | ind_id),
                      data = filter(palmdat, obs_last > 205, yr < 2023)) 
  summary(palm_fo_last)
  confint(palm_fo_last)
  # Some evidence of trend for later dates
  
# Predictions: flowers or buds
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
  
  # Predictions: open flowers
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
  
# Plot predictions
  palm_fl_plot <- ggplot(palm_fl_preds, aes(x = yr)) +
    geom_ribbon(aes(x = yr, ymin = lwr, ymax = upr, group = Date, fill = Date), 
                alpha = 0.3, linetype = 0) +
    geom_line(aes(x = yr, y = fit, group = Date, color = Date)) +
    labs(x = "Year", y = "Mean day of year") +
    scale_color_manual(values = c("forest green", "blue")) +
    scale_fill_manual(values = c("forest green", "blue")) +
    theme_bw() +
    annotate("text", x = 2018, y = 250, label = "A. palmeri, flowers or buds", 
             hjust = 0, vjust = 1, fontface = 2, size = 3) +
    theme(text = element_text(size = 10),
          legend.text = element_text(size = 8))
  palm_fo_plot <- ggplot(palm_fo_preds, aes(x = yr)) +
    geom_ribbon(aes(x = yr, ymin = lwr, ymax = upr, group = Date, fill = Date), 
                alpha = 0.3, linetype = 0) +
    geom_line(aes(x = yr, y = fit, group = Date, color = Date)) +
    labs(x = "Year", y = "Mean day of year") +
    scale_color_manual(values = c("forest green", "blue")) +
    scale_fill_manual(values = c("forest green", "blue")) +
    theme_bw() +
    annotate("text", x = 2018, y = 245, label = "A. palmeri, open flowers", 
             hjust = 0, vjust = 1, fontface = 2, size = 3) +
    theme(text = element_text(size = 10),
          legend.text = element_text(size = 8))
  palm_trends <- plot_grid(palm_fl_plot, palm_fo_plot, ncol = 1)
  ggsave("output/trends/PALM-flowers.png",
          palm_trends, device = "png", width = 6.5, height = 5, 
         units = "in", dpi = 300)

# Flowering intensity: 50-74% of flowers open 
# (using first date)
  # Probably don't need to filter individuals based on number of observations.
  palmdat2 <- filter(ind_yr, spp == "A. palmeri", yr > 2017)
  palm_fo50_first <- lmer(open50_first ~ I(yr - 2012) + (1 | ind_id), 
                         data = palmdat2) 
  summary(palm_fo50_first)
  # Some evidence that peak has been occurring later in recent years
  # Positive slope (1.83), but with lots of uncertainty (SE = 0.64)
  
  # Annual and overall estimates
  palm_fo50_yr <- lm(open50_first ~ factor(yr), data = palmdat2)
  summary(palm_fo50_yr)
  # 2018-2019 lower than all the rest
  palm_fo50_yr <- lm(open50_first ~ factor(yr) - 1, data = palmdat2)
  summary(palm_fo50_yr)
  confint(palm_fo50_yr)
  palm_fo50_1 <- lm(open50_first ~ 1, data = palmdat2)
  summary(palm_fo50_1)
  confint(palm_fo50_1)
  # Overall mean DOY = 210 (SE = 1.3; 95% CI = 208-213)
  AIC(palm_fo50_yr); AIC(palm_fo50_1) # Model with year better.
  
  # Table:
  palm50 <- palmdat2 %>%
    filter(!is.na(open50_first)) %>%
    group_by(yr) %>%
    summarize(no_indiv = length(yr)) %>%
    mutate(yr = as.character(yr)) %>%
    rbind(data.frame(yr = "2018-2023", 
                     no_indiv = sum(!is.na(palmdat2$open50_first)))) %>%
    mutate(est = c(coef(palm_fo50_yr), coef(palm_fo50_1)),
           lcl = c(confint(palm_fo50_yr)[,1], confint(palm_fo50_1)[,1]),
           ucl = c(confint(palm_fo50_yr)[,2], confint(palm_fo50_1)[,2])) %>%
    mutate(spp = "A. palmeri", .before = yr) %>%
    data.frame()
  # write.table(palm50, "clipboard", sep = "\t", row.names = FALSE)

#------------------------------------------------------------------------------#
# Could use predicted start/end/peak dates from GAMs in trend analyses...
#------------------------------------------------------------------------------#
gam_ests <- read.csv("output/gams/estimates-annual.csv")
ests_sag <- filter(gam_ests, taxa == "C. gigantea")
ests_palm <- filter(gam_ests, taxa == "A. palmeri")

# Saguaros: flowers
m_sag_flf <- lm(start ~ I(yr - 2012), data = filter(ests_sag, phase == "flowers"))
summary(m_sag_flf); confint(m_sag_flf)
m_sag_flp <- lm(peak ~ I(yr - 2012), data = filter(ests_sag, phase == "flowers"))
summary(m_sag_flp); confint(m_sag_flp)
m_sag_fle <- lm(end ~ I(yr - 2012), data = filter(ests_sag, phase == "flowers"))
summary(m_sag_fle); confint(m_sag_fle)
# No evidence that start/peak change, but some evidence that end is earlier
# (beta = -7.6, P = 0.12)

# Saguaros: open flowers
m_sag_fof <- lm(start ~ I(yr - 2012), data = filter(ests_sag, phase == "flowers_open"))
summary(m_sag_fof); confint(m_sag_fof)
m_sag_fop <- lm(peak ~ I(yr - 2012), data = filter(ests_sag, phase == "flowers_open"))
summary(m_sag_fop); confint(m_sag_fop)
m_sag_foe <- lm(end ~ I(yr - 2012), data = filter(ests_sag, phase == "flowers_open"))
summary(m_sag_foe); confint(m_sag_foe)
# No evidence of any changes

# A. palmeri: flowers
m_palm_flf <- lm(start ~ I(yr - 2012), data = filter(ests_palm, phase == "flowers"))
summary(m_palm_flf); confint(m_palm_flf)
m_palm_flp <- lm(peak ~ I(yr - 2012), data = filter(ests_palm, phase == "flowers"))
summary(m_palm_flp); confint(m_palm_flp)
m_palm_fle <- lm(end ~ I(yr - 2012), data = filter(ests_palm, phase == "flowers"))
summary(m_palm_fle); confint(m_palm_fle)
# No evidence of changes

# A. palmeri: open flowers
m_palm_fof <- lm(start ~ I(yr - 2012), data = filter(ests_palm, phase == "flowers_open"))
summary(m_palm_fof); confint(m_palm_fof)
m_palm_fop <- lm(peak ~ I(yr - 2012), data = filter(ests_palm, phase == "flowers_open"))
summary(m_palm_fop); confint(m_palm_fop)
m_palm_foe <- lm(end ~ I(yr - 2012), data = filter(ests_palm, phase == "flowers_open"))
summary(m_palm_foe); confint(m_palm_foe)
# No evidence that start/peak change, but some evidence that open flowers phase
# is ending later (beta = 6.4, P = 0.11)
