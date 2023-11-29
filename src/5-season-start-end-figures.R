# Create figures with estimates
# Erin Zylstra
# ezylstra@arizona.edu
# 2023-11-29

library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)
library(cowplot)
library(mgcv)

rm(list = ls())

# Read in start, end, peak estimates from GAM models
flowers_allyrs_10 <- read.csv("output/gams/estimates-allyrs-thresh10.csv")
flowers_allyrs_25 <- read.csv("output/gams/estimates-allyrs-thresh25.csv")
fruit_allyrs_10 <- read.csv("output/gams/estimates-fruit-allyrs-thresh10.csv")
fruit_allyrs_25 <- read.csv("output/gams/estimates-fruit-allyrs-thresh25.csv")

# Get flowering/fruiting dates for 10% threshold
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
  dfe <- dfe %>%
    arrange(taxa, phase) %>%
    mutate(roost = ifelse(taxa == "C. gigantea", "maternity", "transient"),
           group = paste0(roost, " / ", taxa))
  
# Summaries of bat presence from various sources
bats <- read.csv("data/bat-dates.csv") %>%
  mutate(group = ifelse(roost == "maternity", "maternity / C. gigantea", 
                        "transient / A. palmeri"),
         yval = c(2, 3, 1, 3, 2, 1))

# Plotting parameters
  # Length of segment caps in figure
  cap <- 0.1

# Maternity roots in southwestern AZ
  ciga <- dfe %>%
    filter(taxa == "C. gigantea") %>%
    mutate(yval = rep(7:4, each = 2))
  
  matern_title <- expression(paste(bold("A. "), "Southwestern Arizona (maternity roosts)"))
  matern <- ggplot(ciga) +
    geom_segment(aes(x = start, xend = end, y = yval, yend = yval, 
                     group = threshold, color = threshold)) +
    geom_segment(aes(x = start, xend = start, 
                     y = yval - cap, yend = yval + cap, 
                     group = threshold, color = threshold)) +
    geom_segment(aes(x = peak, xend = peak, 
                     y = yval - cap, yend = yval + cap, 
                     group = threshold, color = threshold), linewidth = 1.2) +
    geom_segment(aes(x = end, xend = end, 
                     y = yval - cap, yend = yval + cap, 
                     group = threshold, color = threshold)) +
    scale_color_manual(values = c("gray60", "forestgreen")) + 
    geom_hline(yintercept = 3.5, color = "gray", linetype = 2) +
    annotate("rect", 
             xmin = bats$start[bats$roost == "maternity" & bats$yval == 1],
             xmax = bats$end[bats$roost == "maternity" & bats$yval == 1],
             ymin = 1 - cap, ymax = 1 + cap, fill = "black", color = NA) +
    annotate("rect", 
             xmin = bats$start[bats$roost == "maternity" & bats$yval == 2],
             xmax = bats$end[bats$roost == "maternity" & bats$yval == 2],
             ymin = 2 - cap, ymax = 2 + cap, fill = "black") +
    annotate("rect", 
             xmin = bats$start[bats$roost == "maternity" & bats$yval == 3],
             xmax = bats$end[bats$roost == "maternity" & bats$yval == 3],
             ymin = 3 - cap, ymax = 3 + cap, fill = "black") +
    annotate("text", x = 1, y = 7.25, hjust = 0, vjust = 1, label = "C. gigantea", 
             size = 3, fontface = "italic") +
    annotate("text", x = 1, y = 3.25, hjust = 0, vjust = 1, label = "Bat presence", 
             size = 3) +
    annotate("text", x = 110, y = 0.8, label = "*", size = 10) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          legend.position = c(0.98, 0.99),
          legend.justification = c(1, 1), 
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          legend.key.size = unit(0.5, "cm"),
          legend.spacing.y = unit(0.05, "cm"),
          plot.title = element_text(size = 10))  +
    scale_y_continuous(limits = c(0.75, 7.25), breaks = 1:7, 
                       labels = c("Walker", "BMGR", "SSA", "Ripe fruit", 
                                  "Fruit", "Open flowers", "Flowers")) +
    scale_x_continuous(limits = c(1, 365)) +
    labs(x = "Day of year", y = "", color = "Proportion") +
    ggtitle(matern_title)
  ggsave("output/maternity-roosts-overlap.png", matern, device = "png", 
         width = 6.5, height = 3.2, units = "in", dpi = 300)
  
# Transient roots in southeastern AZ  
  agpa <- dfe %>%
    filter(taxa == "A. palmeri") %>%
    mutate(yval = rep(5:4, each = 2))

  transient_title <- expression(paste(bold("B. "), "Southeastern Arizona (transient roosts)"))
  transient <- ggplot(agpa) +
    geom_segment(aes(x = start, xend = end, y = yval, yend = yval, 
                     group = threshold, color = threshold)) +
    geom_segment(aes(x = start, xend = start, 
                     y = yval - cap, yend = yval + cap, 
                     group = threshold, color = threshold)) +
    geom_segment(aes(x = peak, xend = peak, 
                     y = yval - cap, yend = yval + cap, 
                     group = threshold, color = threshold), linewidth = 1.2) +
    geom_segment(aes(x = end, xend = end, 
                     y = yval - cap, yend = yval + cap, 
                     group = threshold, color = threshold)) +
    scale_color_manual(values = c("gray60", "forestgreen")) + 
    geom_hline(yintercept = 3.5, color = "gray", linetype = 2) +
    annotate("rect", 
             xmin = bats$start[bats$roost == "transient" & bats$yval == 3],
             xmax = bats$end_ext[bats$roost == "transient" & bats$yval == 3],
             ymin = 3 - cap, ymax = 3 + cap, fill = "gray80") +
    annotate("rect", 
             xmin = bats$start[bats$roost == "transient" & bats$yval == 2],
             xmax = bats$end_ext[bats$roost == "transient" & bats$yval == 2],
             ymin = 2 - cap, ymax = 2 + cap, fill = "gray80") +
    annotate("rect", 
             xmin = bats$start[bats$roost == "transient" & bats$yval == 1],
             xmax = bats$end[bats$roost == "transient" & bats$yval == 1],
             ymin = 1 - cap, ymax = 1 + cap, fill = "black", color = NA) +
    annotate("rect", 
             xmin = bats$start[bats$roost == "transient" & bats$yval == 2],
             xmax = bats$end[bats$roost == "transient" & bats$yval == 2],
             ymin = 2 - cap, ymax = 2 + cap, fill = "black") +
    annotate("rect", 
             xmin = bats$start[bats$roost == "transient" & bats$yval == 3],
             xmax = bats$end[bats$roost == "transient" & bats$yval == 3],
             ymin = 3 - cap, ymax = 3 + cap, fill = "black") +
    annotate("text", x = 1, y = 5.25, hjust = 0, vjust = 1, label = "A. palmeri", 
             size = 3, fontface = "italic") +
    annotate("text", x = 1, y = 3.25, hjust = 0, vjust = 1, label = "Bat presence", 
             size = 3) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          legend.position = c(0.98, 0.99),
          legend.justification = c(1, 1), 
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          legend.key.size = unit(0.5, "cm"),
          legend.spacing.y = unit(0.05, "cm"),
          plot.title = element_text(size = 10)) +
    scale_y_continuous(limits = c(0.75, 5.25), breaks = 1:5, 
                       labels = c("Feeders", "Walker", "SSA", 
                                  "Open flowers", "Flowers")) +
    scale_x_continuous(limits = c(1, 365)) +
    labs(x = "Day of year", y = "", color = "Proportion") +
    ggtitle(transient_title)
  ggsave("output/transient-roost-overlap.png", transient, device = "png", 
         width = 6.5, height = 2.5, units = "in", dpi = 300)  
  
