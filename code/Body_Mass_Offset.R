###### Model body mass of F2 litters and determine offset to expected weight#####
#### BWalkenhorst 2026 ####

# clear everything and use garbage collector at start
rm(list = ls(all.names = TRUE)) 
gc() 

# Load required libraries
library(dplyr)
library(zoo)
library(lubridate)
library(ggplot2)
library(brms)
library(rstan)
library(tidyverse)
library(tidybayes)
library(ggokabeito) # colour palette
library(extrafont)# use font_import() on first use
library(bayestestR)

# Custom ggplot theme 
theme_clean <- function() {
  theme_minimal(base_family='Calibri') +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = rel(2), hjust = 0.5),
      axis.title = element_text(face = "bold", size = rel(1.75)),
      axis.text = element_text(face = "bold", size= rel(1.25)),
      strip.text = element_text(face = "bold", size = rel(1.5), color='white'),
      strip.background = element_rect(fill = "grey80", color = NA),
      legend.title = element_text(face = "bold", size = rel(1.25)),
      legend.text = element_text(face = 'italic', size = rel(1)))
}

# set seed to duplicate results
set.seed(42)

# half normal prior and weak normal for intercept
priors_halfnormal <- c(set_prior('normal(0,0.5)', class = 'b', lb = 0), set_prior("normal(0,1)", class = "Intercept"))

# load data
mass_data <- readRDS("../data/mass_data.rds")

#### MODEL ####
setwd("models/")

WEIGHT_AGE <- brms::brm(formula = scale(Weight) ~ scale(AGE_D) + SEX + scale(Rainfall30D) + (1|ID), 
                                     data = final_data, family = gaussian(link='identity'),
                                     chains = 4, iter = 10000, warmup = 2500, cores = 4, backend = "cmdstanr", 
                                     prior = priors_halfnormal,
                                     control = list(max_treedepth = 15, adapt_delta=0.999), init=0, 
                                     threads = threading(4),
                                     file ="WEIGHT_AGE_RAIN_F2")

#### Model details ####
WEIGHT_AGE <- readRDS("WEIGHT_AGE_RAIN_F2.rds")

summary(WEIGHT_AGE)


plot(WEIGHT_AGE)
pp_check(WEIGHT_AGE, ndraws=100)

# get the rope range  = -0.1 * SDy, 0.1 * SDy
# as its scaled, sd = 1 so -0.1 and 0.1 it is!
ropeRange <- c(-0.1* sd(scale(mass_data$Weight)), 0.1 * sd(scale(mass_data$Weight)))

describe_posterior(
  WEIGHT_AGE,
  effects = "all", #fixed vs all (for random effects)
  component = "all",
  rope_range = ropeRange,  
  test = c("p_direction", "p_significance"),
  centrality = "all",
  dispersion = TRUE
)


loo_R2(WEIGHT_AGE)

bayes_R2(WEIGHT_AGE)


performance::variance_decomposition(WEIGHT_AGE)

### predict body mass for F2 pups ####
setwd("../../data/")

# F2 data with average weight and rainfall
F2_data <- readRDS("../../data/PROP_data.rds")

# ensure correct format as copy
F2_data <- F2_data %>%
  mutate(TREATMENT = as.factor(TREATMENT), SEX = as.factor(SEX), ID = as.factor(ID))

# #reorder treatments:
F2_data$TREATMENT <- factor(F2_data$TREATMENT, levels = c("DC", "SC", "DT"))

copy_data <- F2_data
copy_data$AGE_D <- copy_data$REC_AGE_D
copy_data$Rainfall30D <- copy_data$MonthlyRainfall

# predict the body mass for recording days
weight_pred <- predict(WEIGHT_AGE, newdata = copy_data, seed=23, allow_new_levels=T)

saveRDS(weight_pred, "FLUT_RAIN_weight.rds")
rm(copy_flut)

#load predictions if not calculated
weight_pred <- readRDS("FLUT_RAIN_weight.rds")

############## determine BODY MASS OFFSET ####
F2_data$WEIGHT_DIFF <- F2_data$AvgWeight - F2_data$WEIGHT_PRED 

#use percentage
F2_data$WEIGHT_DIFF_PER <-as.numeric(((F2_data$AvgWeight - F2_data$WEIGHT_PRED)/F2_data$WEIGHT_PRED) *100)

# plot body mass offset
ggplot(F2_data, aes(x = REC_AGE_D, y = WEIGHT_DIFF_PER, color = TREATMENT)) +
  geom_point() +
 # geom_smooth(method='loess', )+
  geom_smooth(method='glm' )+
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_color_okabe_ito(order = c(2, 1, 7), name = "Offspring", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age",
       y = "Body mass offset (%)",
       color = "Grandmaternal treatment") +
  theme_clean()

saveRDS(F2_data, "PROP_data_bodymass.rds")




