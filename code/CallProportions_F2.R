#### Analysis script for F2 androgen study: CALL PROPORTIONS ###################
# BWalkenhorst, 2026 

#### SETUP ####
# clear everything and use garbage collector at start
rm(list = ls(all.names = TRUE)) #includes hidden objects.
gc() 

# load all necessary libraries
library(readxl) 
library(ggplot2) 
library(ggpubr) 
library(car)
library(tidyverse)
library(tidybayes) 
library(brms)
library(rstan)
library(bayestestR) #e.g. diagnostic_posterior
library(bayesplot)
library(ggokabeito) # colour palette
library(emmeans) # 
library(extrafont)# use font_import() on first use

set.seed(23)

# set working directory
setwd("../data")

PROP_data <- readRDS('PROP_data.rds')

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

get_age_vars <- function(){
  age_vars <- c((30- mean(PROP_data$REC_AGE))/sd(PROP_data$REC_AGE), 
                (75- mean(PROP_data$REC_AGE))/sd(PROP_data$REC_AGE),
                (120- mean(PROP_data$REC_AGE))/sd(PROP_data$REC_AGE))
  return(age_vars)
}

get_org_value_for_z <- function(z_value, column){
  return ((z_value*sd(column)) + mean(column))   
}

get_z_value_for_org <- function(org_value, column){
  return ((org_value - mean(column))/sd(column))   
}

# Data overview
# AGE
PROP_data %>%
  group_by(TREATMENT) %>%
  drop_na(REC_AGE)%>%
  reframe(
    Mean = mean(REC_AGE),
    SD = sd(REC_AGE),
    Min = min(REC_AGE),
    Max = max(REC_AGE)
  )

mean(PROP_data$REC_AGE, na.rm = TRUE)
sd(PROP_data$REC_AGE, na.rm = TRUE )

# body mass offset:
PROP_data %>%
  group_by(TREATMENT) %>%
  reframe(
    Mean = mean(WEIGHT_DIFF_PER),
    SD = sd(WEIGHT_DIFF_PER),
    Min = min(WEIGHT_DIFF_PER),
    Max = max(WEIGHT_DIFF_PER)
  )


# GS adults:
PROP_data %>%
  group_by(TREATMENT) %>%
  summarize(
    Mean = mean(GS_adults_REC, na.rm = T),
    SD = sd(GS_adults_REC, na.rm = T),
    Min = min(GS_adults_REC, na.rm = T),
    Max = max(GS_adults_REC, na.rm = T)
  )

mean(PROP_data$GS_adults_REC)
sd(PROP_data$GS_adults_REC)

# GS pups (<130 d):
PROP_data %>%
  group_by(TREATMENT) %>%
  summarize(
    Mean = mean(GS_pups_REC_litter, na.rm = T),
    SD = sd(GS_pups_REC_litter, na.rm = T),
    Min = min(GS_pups_REC_litter, na.rm = T),
    Max = max(GS_pups_REC_litter, na.rm = T)
  )

mean(PROP_data$GS_pups_REC_litter)
sd(PROP_data$GS_pups_REC_litter)

# GS (all):
PROP_data %>%
  group_by(TREATMENT) %>%
  summarize(
    Mean = mean(GROUPSIZE, na.rm = T),
    SD = sd(GROUPSIZE, na.rm = T),
    Min = min(GROUPSIZE, na.rm = T),
    Max = max(GROUPSIZE, na.rm = T)
  )

mean(PROP_data$GROUPSIZE)
sd(PROP_data$GROUPSIZE)

# COMPETITION:
PROP_data %>%
  group_by(TREATMENT) %>%
  summarize(
    Mean = mean(COMPETITION, na.rm = T),
    SD = sd(COMPETITION, na.rm = T),
    Min = min(COMPETITION, na.rm = T),
    Max = max(COMPETITION, na.rm = T)
  )

mean(PROP_data$COMPETITION)
sd(PROP_data$COMPETITION)

# COMPETITION normalised:
PROP_data %>%
  group_by(TREATMENT) %>%
  summarize(
    Mean = mean(COMP_NORM, na.rm = T),
    SD = sd(COMP_NORM, na.rm = T),
    Min = min(COMP_NORM, na.rm = T),
    Max = max(COMP_NORM, na.rm = T)
  )

mean(PROP_data$COMP_NORM)
sd(PROP_data$COMP_NORM)

# litter size
PROP_data %>%
  group_by(TREATMENT) %>%
  reframe(
    Mean = mean(LITTER_SIZE),
    SD = sd(LITTER_SIZE),
    Min = min(LITTER_SIZE),
    Max = max(LITTER_SIZE)
  )

mean(PROP_data$LITTER_SIZE)
sd(PROP_data$LITTER_SIZE)

# Monthly rainfall
PROP_data %>%
  group_by(TREATMENT) %>%
  reframe(
    Mean = mean(MonthlyRainfall),
    SD = sd(MonthlyRainfall),
    Min = min(MonthlyRainfall),
    Max = max(MonthlyRainfall)
  )

mean(PROP_data$MonthlyRainfall)
sd(PROP_data$MonthlyRainfall)

# overview plots 
# weight diff raw data plot
ggplot(PROP_data, aes(x = REC_AGE, y = WEIGHT_DIFF_PER, color = TREATMENT)) +
  geom_point() +
  # geom_smooth(method='loess')+
  geom_smooth(method='glm')+
  labs(x = "Age",
       y = "Body mass offset (%)",
       color = "Treatment") +
  scale_color_okabe_ito(order = c(2, 1, 7), name = "Grandmaternal\ntreatment", labels = c('DC', 'DT', 'SC')) +
  theme_clean()


# y distributions -RAW
plot(density(PROP_data$Sum_BEG/PROP_data$Total_calls), 
     xlab = "REP proportions",
     ylab = "Density")

plot(density(PROP_data$Sum_DIG/PROP_data$Total_calls), 
     xlab = "DIG proportions",
     ylab = "Density")

plot(density(PROP_data$Sum_CC/PROP_data$Total_calls), 
     xlab = "CC proportions",
     ylab = "Density")

sd(PROP_data$REC_AGE)
sd(PROP_data$WEIGHT_DIFF_PER)
sd(PROP_data$COMP_NORM)
sd(PROP_data$GROUPSIZE)
sd(PROP_data$MonthlyRainfall)


########################## BAYES MODELS ########################################
#  MULTIVARIATE MODEL ####

setwd('../code/models')

priors_refined <- c(
  set_prior("normal(0, 1.5)", class = "Intercept", resp = "SumBEG"),
  set_prior("normal(0, 1.5)", class = "Intercept", resp = "SumDIG"),
  set_prior("normal(0, 1.5)", class = "Intercept", resp = "SumCC"),
  
  set_prior("normal(0, 0.8)", class = "b", resp = "SumBEG"),
  set_prior("normal(0, 0.8)", class = "b", resp = "SumDIG"),
  set_prior("normal(0, 0.8)", class = "b", resp = "SumCC"),
  
  set_prior("student_t(3, 0, 1.5)", class = "sd", resp = "SumBEG"),
  set_prior("student_t(3, 0, 1.5)", class = "sd", resp = "SumDIG"),
  set_prior("student_t(3, 0, 1.5)", class = "sd", resp = "SumCC")
)

bf_REP <- bf(Sum_BEG | trials(Total_calls) ~ TREATMENT * AGE_z * MaternalRank + 
                AGE_z * SEX + WEIGHT_z + COMP_NORM_z + GS_z + RAIN_z + (1 | ID))

bf_DIG <- bf(Sum_DIG | trials(Total_calls) ~ TREATMENT * (AGE_z + I(AGE_z^2)) * MaternalRank + 
               + (AGE_z + I(AGE_z^2)) * SEX + WEIGHT_z + COMP_NORM_z + GS_z + RAIN_z + (1 | ID))

bf_CC <- bf(Sum_CC | trials(Total_calls) ~ TREATMENT * AGE_z * MaternalRank + 
              AGE_z * SEX + WEIGHT_z + COMP_NORM_z + GS_z + RAIN_z + (1 | ID))

multivar_formula <- mvbrmsformula(bf_REP, bf_DIG, bf_CC)

B_prop <- brms::brm(formula = multivar_formula,
                            data = PROP_data, family = beta_binomial(link='logit'),
                            chains = 4, iter = 6000, warmup = 1500, seed = 23, control = list(max_treedepth = 20, adapt_delta=0.99),
                            save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                            threads = threading(4),
                            prior = priors_refined,
                            file="B_prop")
B_prop<- add_criterion(B_prop, c("loo", "loo_R2"), moment_match = TRUE)

#### RESULTS: MULTIVAR model ####
B_prop <- readRDS("B_prop.rds")

summary(B_prop)

plot(B_prop)
pp_check(B_prop, ndraws = 100, resp='SumBEG')
pp_check(B_prop, ndraws = 100, resp='SumDIG')
pp_check(B_prop, ndraws = 100, resp='SumCC')

posterior <- describe_posterior(
  B_prop$fit,
  effects = "all", #fixed vs all (for random effects)
  component = "all",
  rope_range = c(-0.18, 0.18),  
  test = c("p_direction", "p_significance"),
  centrality = "all",
  #centrality = "median",
  dispersion = TRUE
)
(output <- posterior[1:70,])
saveRDS(output, file='posterior_desc.rds')

loo_R2(B_prop, moment_match=T)

bayes_R2(B_prop)

### REP: TAM ####
REP_ROPE <- c(-0.18, 0.18)
(mat_stat<- emtrends(B_prop, pairwise ~ TREATMENT:MaternalRank, var="AGE_z", resp='SumBEG'))

pd(mat_stat)

p_significance(mat_stat, threshold = REP_ROPE)

rm(mat_stat)

### REP at 30 days (intercept) ####
REP_ROPE <- c(-0.18, 0.18)
(REP_int <- emmeans(B_prop,
                           pairwise ~ TREATMENT * MaternalRank | AGE_z,
                           at = list(AGE_z = get_z_value_for_org(30, PROP_data$REC_AGE)),
                           resp='SumBEG' ))
                   
pd(REP_int)


p_significance(REP_int, threshold = REP_ROPE)

rm(REP_int)

# Plot emmeans: REP prop TAM ####
treat_mat_stat <- emmeans(B_prop,
                          ~TREATMENT * MaternalRank | AGE_z,
                          at = list(AGE_z = get_z_value_for_org(30, PROP_data$REC_AGE)),
                          resp='SumBEG')
emm_data <- as.data.frame(treat_mat_stat)
emm_data$REC_AGE <- emm_data$AGE_z * sd(PROP_data$REC_AGE) + mean(PROP_data$REC_AGE)
emm_data$REP_prop <- inv_logit_scaled(emm_data$emmean)
emm_data$HPD_low <- inv_logit_scaled(emm_data$lower.HPD)
emm_data$HPD_high <- inv_logit_scaled(emm_data$upper.HPD)
emm_data$TREATMENT <- factor(emm_data$TREATMENT, levels = c("DC", "SC", "DT"))

# 800* 500: EMMs_REP_prop_TAM
ggplot(emm_data, aes(x = as.factor(REC_AGE), y = REP_prop, color = TREATMENT)) +
  geom_point(position = position_dodge(width = 1), size = 3) +  # Show estimates as points
  geom_errorbar(aes(ymin = HPD_low, ymax = HPD_high), position = position_dodge(width = 1), width = 0.75) +
  scale_color_okabe_ito(order= c(2, 1, 7), name = "Grandmaternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_y_continuous(                   limits = c(0.0, 1.0),
                                        breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  labs(x = "Age (days)", y = "Repeat call proportion\n") +
  theme_clean()+
  facet_wrap(~MaternalRank)

rm(treat_mat_stat, emm_data)

### DIG slopes ####
DIG_ROPE <- c(-0.18, 0.18)
(treat_mat_stat <- emtrends(B_prop, pairwise ~ TREATMENT:MaternalRank, var = "AGE_z", max.degree = 2, resp = "SumDIG"))

pd(treat_mat_stat)


p_significance(treat_mat_stat, threshold = DIG_ROPE)

rm(treat_mat_stat)

### DIG: aged 75 days ####
DIG_ROPE <- c(-0.18, 0.18)
(treat_mat_stat <- emmeans(B_prop,
                          pairwise ~ TREATMENT * MaternalRank | AGE_z,
                           at = list(AGE_z = get_z_value_for_org(75, PROP_data$REC_AGE)),
                           resp='SumDIG' ))

pd(treat_mat_stat)

p_significance(treat_mat_stat, threshold = DIG_ROPE)

rm(treat_mat_stat)

# Plot emmeans: DIG prop TAM ####
treat_mat_stat <- emmeans(B_prop,
                           ~TREATMENT * MaternalRank | AGE_z,
                           at = list(AGE_z = get_z_value_for_org(75, PROP_data$REC_AGE)),
                           resp='SumDIG')
emm_data <- as.data.frame(treat_mat_stat)
emm_data$REC_AGE <- emm_data$AGE_z * sd(PROP_data$REC_AGE) + mean(PROP_data$REC_AGE)
emm_data$DIG_prop <- inv_logit_scaled(emm_data$emmean)
emm_data$HPD_low <- inv_logit_scaled(emm_data$lower.HPD)
emm_data$HPD_high <- inv_logit_scaled(emm_data$upper.HPD)
emm_data$TREATMENT <- factor(emm_data$TREATMENT, levels = c("DC", "SC", "DT"))

# 800* 500: EMMs_DIG_prop_TAM
ggplot(emm_data, aes(x = as.factor(REC_AGE), y = DIG_prop, color = TREATMENT)) +
  geom_point(position = position_dodge(width = 1), size = 3) +  # Show estimates as points
  geom_errorbar(aes(ymin = HPD_low, ymax = HPD_high), position = position_dodge(width = 1), width = 0.75) +
  scale_color_okabe_ito(order= c(2, 1, 7), name = "Grandmaternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_y_continuous(                   limits = c(0.0, 1.0),
                                        breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  labs(x = "Age (days)", y = "Digging call proportion\n") +
  theme_clean()+
  facet_wrap(~MaternalRank)

rm(treat_mat_stat, emm_data)

### CC TAM ####
CC_ROPE <- c(-0.18, 0.18)
(treat_mat <- emtrends(B_prop, pairwise ~ TREATMENT:MaternalRank, var="AGE_z", 
                       resp='SumCC'))

pd(treat_mat)

p_significance(treat_mat, threshold = CC_ROPE)

rm(treat_mat)

### CC at 120 days ####
CC_ROPE <- c(-0.18, 0.18)
(treat_mat_stat <- emmeans(B_prop,
                           pairwise ~ TREATMENT * MaternalRank | AGE_z,
                           at = list(AGE_z = get_z_value_for_org(120, PROP_data$REC_AGE)),
                           resp='SumCC'))
                      
pd(treat_mat_stat)

p_significance(treat_mat_stat, threshold = CC_ROPE)

rm(treat_mat_stat,  CC_ROPE, DIG_ROPE, REP_ROPE)

# Plot emmeans: CC prop TAM ####
treat_mat_stat <- emmeans(B_prop,
                          ~TREATMENT * MaternalRank | AGE_z,
                          at = list(AGE_z = get_z_value_for_org(120, PROP_data$REC_AGE)),
                          resp='SumCC')
emm_data <- as.data.frame(treat_mat_stat)
emm_data$REC_AGE <- emm_data$AGE_z * sd(PROP_data$REC_AGE) + mean(PROP_data$REC_AGE)
emm_data$CC_prop <- inv_logit_scaled(emm_data$emmean)
emm_data$HPD_low <- inv_logit_scaled(emm_data$lower.HPD)
emm_data$HPD_high <- inv_logit_scaled(emm_data$upper.HPD)
emm_data$TREATMENT <- factor(emm_data$TREATMENT, levels = c("DC", "SC", "DT"))

# 800* 500: EMMs_CC_prop_TAM
ggplot(emm_data, aes(x = as.factor(REC_AGE), y = CC_prop, color = TREATMENT)) +
  geom_point(position = position_dodge(width = 1), size = 3) +  # Show estimates as points
  geom_errorbar(aes(ymin = HPD_low, ymax = HPD_high), position = position_dodge(width = 1), width = 0.75) +
  scale_color_okabe_ito(order= c(2, 1, 7), name = "Grandmaternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_y_continuous(                   limits = c(0.0, 1.0),
                                        breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  labs(x = "Age (days)", y = "Close call proportion\n") +
  theme_clean()+
  facet_wrap(~MaternalRank)

rm(emm_data, treat_mat_stat)

#### Coefficient plots ####
posterior_desc <- readRDS('posterior_desc.rds')

REP_desc <- posterior_desc %>%
  filter(str_detect(Parameter, "SumBEG"))
DIG_desc <- posterior_desc %>%
  filter(str_detect(Parameter, "SumDIG"))
CC_desc <- posterior_desc %>%
  filter(str_detect(Parameter, "SumCC"))

REP_desc <- REP_desc[(2:18),]
DIG_desc <- DIG_desc[(2:25),]
CC_desc <- CC_desc[c(2:18),]

# define control vars
control_vars <- c('Monthly rainfall', 'Group size', 'Competition', 'Body mass offset',
                  'Age:Male', 'SUB mother:Male', 'Age^2:Male')

# REP Coef
# clean up labels:
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'b_SumBEG_', '')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'WEIGHT_z', 'Body mass offset')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'GS_z', 'Group size')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'AGE_z', 'Age')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'COMP_NORM_z', 'Competition')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'RAIN_z', 'Monthly rainfall')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'SEXM', 'Male')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'MaternalRankSUB', 'SUB mother')
REP_desc$TREATMENT <- ifelse(grepl("SC", REP_desc$Parameter), "SC",
                            ifelse(grepl("DT", REP_desc$Parameter), "DT", 
                                   ifelse(REP_desc$Parameter %in% control_vars, "control", "DC")))
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'TREATMENTDT', 'DT')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'TREATMENTSC', 'SC')

custom_order <- c('Monthly rainfall', 
                  'Group size', 
                  'Competition', 
                  'Body mass offset',
                  'SUB mother:Male', 
                  'Age:Male', 
                  'DT:Age:SUB mother','SC:Age:SUB mother','Age:SUB mother',
                  'DT:SUB mother','SC:SUB mother','SUB mother', 
                  'DT:Male','SC:Male',"Male",                   
                  'DT:Age','SC:Age',"Age", 
                  'DT','SC') 


# Update the order of TREATMENT factor levels
REP_desc$TREATMENT <- factor(REP_desc$TREATMENT, levels = c("DC", "SC", "DT", 'control'))
REP_desc$Parameter <- factor(REP_desc$Parameter, levels = custom_order)
# Coeff_REP 700 * 800
ggplot(REP_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 7, 8), name = "Effect group", labels = c('DC', 'SC', 'DT', 'Population-level\ncovariates'))+
  labs(x = "Parameter value", y = "Parameters") +
  theme_clean()

# DIG coeff
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'b_SumDIG_', '')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'WEIGHT_z', 'Body mass offset')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'GS_z', 'Group size')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'COMP_NORM_z', 'Competition')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'RAIN_z', 'Monthly rainfall')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'SEXM', 'Male')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'MaternalRankSUB', 'SUB mother')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'IAGE_zE2', 'Age^2')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'AGE_z', 'Age')
DIG_desc$TREATMENT <- ifelse(grepl("SC", DIG_desc$Parameter), "SC",
                            ifelse(grepl("DT", DIG_desc$Parameter), "DT", 
                                   ifelse(DIG_desc$Parameter %in% control_vars, "control", "DC")))
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'TREATMENTDT', 'DT')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'TREATMENTSC', 'SC')

custom_order <- c('Monthly rainfall', 
                  'Group size', 
                  'Competition', 
                  'Body mass offset',
                  'SUB mother:Male', 
                  'Age^2:Male', 
                  'Age:Male', 
                  'DT:Age^2:SUB mother','SC:Age^2:SUB mother','Age^2:SUB mother',
                  'DT:Age:SUB mother','SC:Age:SUB mother','Age:SUB mother',
                  'DT:SUB mother','SC:SUB mother','SUB mother', 
                  'DT:Male','SC:Male',"Male", 
                  'DT:Age^2','SC:Age^2',"Age^2", 
                  'DT:Age','SC:Age',"Age", 
                  'DT','SC') 


DIG_desc$Parameter <- factor(DIG_desc$Parameter, levels = custom_order)
DIG_desc$TREATMENT <- factor(DIG_desc$TREATMENT, levels = c("DC", "SC", "DT", "control"))

# Coeff_DIG 700*900
ggplot(DIG_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 7, 8), name = "Effect group", labels = c('DC', 'SC', 'DT', 'Population-level\ncovariates'))+
  labs(x = "Parameter value", y = "Parameters") +
  theme_clean()

# CC Coeff
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'b_SumCC_', '')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'WEIGHT_z', 'Body mass offset')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'GS_z', 'Group size')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'RAIN_z', 'Monthly rainfall')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'AGE_z', 'Age')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'COMP_NORM_z', 'Competition')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'SEXM', 'Male')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'MaternalRankSUB', 'SUB mother')
CC_desc$TREATMENT <- ifelse(grepl("SC", CC_desc$Parameter), "SC",
                                   ifelse(grepl("DT", CC_desc$Parameter), "DT", 
                                          ifelse(CC_desc$Parameter %in% control_vars, "control", "DC")))
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'TREATMENTDT', 'DT')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'TREATMENTSC', 'SC')

custom_order <- c('Monthly rainfall', 
                  'Group size', 
                  'Competition', 
                  'Body mass offset',
                  'SUB mother:Male', 
                  'Age:Male', 
                  'DT:Age:SUB mother','SC:Age:SUB mother','Age:SUB mother',
                  'DT:SUB mother','SC:SUB mother','SUB mother', 
                  'DT:Male','SC:Male',"Male",                   
                  'DT:Age','SC:Age',"Age", 
                  'DT','SC') 

CC_desc$Parameter <- factor(CC_desc$Parameter, levels = custom_order)
CC_desc$TREATMENT <- factor(CC_desc$TREATMENT, levels = c("DC", "SC", "DT", 'control'))

# Coeff_CC 700*800
ggplot(CC_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 7, 8), name = "Effect group", labels = c('DC', 'SC', 'DT', 'Population-level\ncovariates'))+
  labs(x = "Parameter value", y = "Parameters") +
  theme_clean()

rm(REP_desc, DIG_desc, CC_desc, posterior_desc, custom_order)

#### MODEL ONTOGENY PLOTS ####
B_prop <- readRDS('B_prop.rds')

# get all needed values
sd_age <- sd(PROP_data$REC_AGE)
mean_age <- mean(PROP_data$REC_AGE)

rec_age_c <- seq(30, 130, by=1)
age_z_vals <- seq(min(age_z_vals), max(age_z_vals), length.out = 20)

#predictions based on mean values
PROP_pred <- B_prop %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(PROP_data$TREATMENT),
                                    SEX = levels(PROP_data$SEX),
                                    MaternalRank = levels(PROP_data$MaternalRank),
                                    AGE_z = age_z_vals,
                                    WEIGHT_z = mean(PROP_data$WEIGHT_z),
                                    COMP_NORM_z = mean(PROP_data$COMP_NORM_z),
                                    GS_z = mean(PROP_data$GS_z),
                                    RAIN_z = mean(PROP_data$RAIN_z),
                                    Total_calls=1), 
              re_formula = NA,  robust = T)


#unscale AGE_z values:
PROP_pred$REC_AGE <- PROP_pred$AGE_z * sd_age + mean_age
# ensure right format
PROP_pred$Call_prop <- PROP_pred$.epred
PROP_pred$Call_type <- as.factor(PROP_pred$.category)
PROP_pred$SEX <- as.factor(PROP_pred$SEX)
PROP_pred$MaternalRank <- as.factor(PROP_pred$MaternalRank)
PROP_pred$TREATMENT <- factor(PROP_pred$TREATMENT, levels = c("DC", "SC", "DT"))

# just REP ###
REP_pred <- subset(PROP_pred, Call_type == 'SumBEG')
REP_pred$REP_prop <- REP_pred$Call_prop

# TAM 800x500
ggplot(REP_pred, aes(x = REC_AGE, y = REP_prop, color = TREATMENT, fill = TREATMENT)) +
  stat_lineribbon(.width = .95) +
  geom_point(data = PROP_data, size = 1.5, alpha = 1) +  # raw data
  scale_color_okabe_ito(order = c(2, 1, 7), name = "Grandmaternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 7), alpha = 0.2, name = "Grandmaternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Repeat call proportion\n",
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean() +
  facet_wrap(~MaternalRank)

# MAS 800x500
ggplot(REP_pred, aes(x = REC_AGE, y = REP_prop, color = MaternalRank, fill = MaternalRank)) +
  stat_lineribbon(.width = .95) +
  geom_point(data = PROP_data, size = 1.5, alpha = 1) +  # raw data
  scale_color_okabe_ito(order = c(2, 1), name = "Maternal\nstatus", labels = c('dominant', 'subordinate')) +
  scale_fill_okabe_ito(order = c(2, 1), alpha = 0.2, name = "Maternal\nstatus", labels = c('dominant', 'subordinate')) +
  labs(x = "Age (days)") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Repeat call proportion\n",
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean()+
  facet_wrap(~SEX)

rm(REP_pred)

# # just DIG
DIG_pred <- subset(PROP_pred, Call_type == 'SumDIG')
DIG_pred$DIG_prop <- DIG_pred$Call_prop

# TAM 800*500
ggplot(DIG_pred, aes(x = REC_AGE, y = DIG_prop, color = TREATMENT, fill = TREATMENT)) +
  stat_lineribbon(.width = .95) +
  geom_point(data = PROP_data, size = 1.5, alpha = 1) +  # raw data
  scale_color_okabe_ito(order = c(2, 1, 7), name = "Grandmaternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 7), alpha = 0.2, name = "Grandmaternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Digging call proportion\n",
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean() +
  facet_wrap(~MaternalRank)

# MAS 800*500
ggplot(DIG_pred, aes(x = REC_AGE, y = DIG_prop, color = MaternalRank, fill = MaternalRank)) +
  stat_lineribbon(.width = .95) +
  geom_point(data = PROP_data, size = 1.5, alpha = 1) +  # raw data
  scale_color_okabe_ito(order = c(2, 1), name = "Maternal\nstatus", labels = c('dominant', 'subordinate')) +
  scale_fill_okabe_ito(order = c(2, 1), alpha = 0.2, name = "Maternal\nstatus", labels = c('dominant', 'subordinate')) +
  labs(x = "Age (days)") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Digging call proportion\n",
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean()+
  facet_wrap(~SEX)

rm(DIG_pred)

# just CC
CC_pred <- subset(PROP_pred, Call_type == 'SumCC')
CC_pred$CC_prop <- CC_pred$Call_prop

# TAM 800*500
ggplot(CC_pred, aes(x = REC_AGE, y = CC_prop, color = TREATMENT, fill = TREATMENT)) +
  stat_lineribbon(.width = .95) +
  geom_point(data = PROP_data, size = 1.5, alpha = 1) +  # raw data
  scale_color_okabe_ito(order = c(2, 1, 7), name = "Grandmaternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 7), alpha = 0.2, name = "Grandmaternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Close call proportion\n",
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean() +
  facet_wrap(~MaternalRank)

# MAS
ggplot(CC_pred, aes(x = REC_AGE, y = CC_prop, color = MaternalRank, fill = MaternalRank)) +
  stat_lineribbon(.width = .95) +
  geom_point(data = PROP_data, size = 1.5, alpha = 1) +  # raw data
  scale_color_okabe_ito(order = c(2, 1), name = "Maternal\nstatus", labels = c('dominant', 'subordinate')) +
  scale_fill_okabe_ito(order = c(2, 1), alpha = 0.2, name = "Maternal\nstatus", labels = c('dominant', 'subordinate')) +
  labs(x = "Age (days)", y = "Close call proportion\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Close call proportion\n",
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean()+
  facet_wrap(~SEX)

rm(CC_pred)

#clean up
rm(list = ls(all.names = TRUE)) #includes hidden objects.
gc() 
