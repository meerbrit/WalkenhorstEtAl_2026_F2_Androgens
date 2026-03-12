### Analyse F2 milestone ages ####
### BWalkenhorst, 2026 #########################################################

#### SETUP ####
# clear everything and use garbage collector at start
rm(list = ls(all.names = TRUE)) #includes hidden objects.
gc() 

# load all necessary libraries
library(readxl) 
library(writexl)
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

# ###  load transition ages:   ONLY NEEDED ONCE
# GEN_data <- readRDS('../data/MILESTONES_POSTERIOR.rds')

# ensure F1 status is in correct format
# GEN_data$MAT_STAT <- factor(GEN_data$MAT_STAT, levels=c('DOM', 'SUB'))

# saveRDS(GEN_data, '../data/Milestones_F2.rds')

GEN_data <- readRDS('../data/Milestones_F2.rds')

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

#### Milestone ages overview ####
GEN_data %>%
  group_by(MAT_STAT, TREATMENT, SEX) %>%
  summarize(
    Mean = mean(SEMI_AGE, na.rm = T),
    SD = sd(SEMI_AGE, na.rm = T),
    Min = min(SEMI_AGE, na.rm = T),
    Max = max(SEMI_AGE, na.rm = T)
    )

GEN_data %>%
  group_by(MAT_STAT, TREATMENT, SEX) %>%
  summarize(
    Mean = mean(PEAK_DIG_AGE, na.rm = T),
    SD = sd(PEAK_DIG_AGE, na.rm = T),
    Min = min(PEAK_DIG_AGE, na.rm = T),
    Max = max(PEAK_DIG_AGE, na.rm = T)
  )

GEN_data %>%
  group_by(MAT_STAT, TREATMENT, SEX) %>%
  summarize(
    Mean = mean(FULL_AGE, na.rm = T),
    SD = sd(FULL_AGE, na.rm = T),
    Min = min(FULL_AGE, na.rm = T),
    Max = max(FULL_AGE, na.rm = T)
  )

################################################################################
################################################################################
### RESULTS ####
################################################################################
################################################################################

# Create datasets for each milestone - only needed once
setwd("../data/")
milestone_data <- readRDS('Milestones_F2.rds')
 
# SEMI_data <- milestone_data %>%
#   select(MAT_STAT, ID, TREATMENT, SEX, WEIGHT, COMP, GS, RAIN, MOTHER_ID, LITTER_CODE, .draw, SEMI_AGE) %>%
#   drop_na(SEMI_AGE) %>%
#   rename(Milestone_Age = SEMI_AGE) %>%
#   mutate(Milestone_Type = "SEMI")
# saveRDS(SEMI_data, 'SEMI_data_F2.rds')
# 
# PEAK_data <- milestone_data %>%
#   select(MAT_STAT, ID, TREATMENT, SEX, WEIGHT, COMP, GS, RAIN, MOTHER_ID, LITTER_CODE, .draw, PEAK_DIG_AGE) %>%
#   drop_na(PEAK_DIG_AGE) %>%
#   rename(Milestone_Age = PEAK_DIG_AGE) %>%
#   mutate(Milestone_Type = "PEAK")
# saveRDS(PEAK_data, 'PEAK_data_F2.rds')
# 
# FULL_data <- milestone_data %>%
#   select(MAT_STAT, ID, TREATMENT, SEX, WEIGHT, COMP, GS, RAIN, MOTHER_ID, LITTER_CODE, .draw, FULL_AGE) %>%
#   drop_na(FULL_AGE) %>%
#   rename(Milestone_Age = FULL_AGE) %>%
#   mutate(Milestone_Type = "FULL")
# saveRDS(FULL_data, 'FULL_data_F2.rds')

# Load data
SEMI_data <- readRDS('SEMI_data_F2.rds')
DIG_PEAK_data <- readRDS('PEAK_data_F2.rds')
FULL_data <- readRDS('FULL_data_F2.rds')

#### MODELS ####
setwd("../code/models/")

# REP -> DIG: SEMI ####

# plot(density(SEMI_data$Milestone_Age),
#      xlab = "Transition ages",
#      ylab = "Density")

# range(SEMI_data$Milestone_Age) 
# mean(SEMI_data$Milestone_Age) 
# sd(SEMI_data$Milestone_Age)# 

priors_SEMI_f2 <- c(
  set_prior("normal(75, 15)", class = "Intercept"),
  set_prior("normal(0, 15)", class = "b"), # +- 30 days
  set_prior("student_t(3, 0, 10)", class = "sd") # ID
)

SEMI_F2 <- brm(formula = Milestone_Age ~ TREATMENT * MAT_STAT * SEX + (1|ID),
              data = SEMI_data,
              family = student(link='identity'),
              chains = 4, iter = 5000, warmup = 1500, seed = 42234223,
              control = list(adapt_delta = 0.99, max_treedepth = 15),
              cores=4,
              backend = 'cmdstanr', init= 'random',
              prior = priors_SEMI_f2, #threads = threading(4),
              file="SEMI_F2"
)

# Results: SEMI ####
SEMI_F2 <- readRDS("SEMI_F2.rds")

summary(SEMI_F2)


plot(SEMI_F2)
pp_check(SEMI_F2, ndraws=100)

loo_R2(SEMI_F2, moment_match = T)

bayes_R2(SEMI_F2)

posterior_desc <- describe_posterior(
  SEMI_F2,
  effects = "all", #"fixed", 
  component = "all",
  rope_range = rope_range(SEMI_F2),  
  test = c("p_direction", "p_significance"),
  centrality = "all",
  dispersion = TRUE
)

# drop sigma
posterior_desc <- posterior_desc[-c(1, 13:112),]

# clean up labels:
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'b_', '')
posterior_desc$TREATMENT <- ifelse(grepl("TREATMENTSC", posterior_desc$Parameter), "SC",
                                   ifelse(grepl("TREATMENTDT", posterior_desc$Parameter), "DT", "DC"))
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTSC', 'SC')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTDT', 'DT')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'SEXM', 'Male')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'MAT_STATSUB', 'SUB mother')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'MAT_STATDOM', 'DOM mother')

custom_order <- c('DT:SUB mother:Male', 'SC:SUB mother:Male', "SUB mother:Male",
                  'DT:DOM mother:Male', 'SC:DOM mother:Male', "DOM mother:Male", 
                  'DT:Male', 'SC:Male', "Male", 
                  'DT:SUB mother','SC:SUB mother','SUB mother',
                  'DT:DOM mother','SC:DOM mother','DOM mother',
                  'DT','SC')

posterior_desc$Parameter <- factor(posterior_desc$Parameter, levels = custom_order)
posterior_desc$TREATMENT <- factor(posterior_desc$TREATMENT, levels = c("DC", "SC", "DT"))

# Coef_SEMI_F2 700*500
ggplot(posterior_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 7), name = "F0 treatment", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
  theme_clean()

### SEMI F2 EMMs####
(treat_mat<- emmeans(SEMI_F2, pairwise ~ TREATMENT:MAT_STAT))


pd(treat_mat)


p_significance(treat_mat, threshold = rope_range(SEMI_F2))

rm(treat_mat)

# Plot SEMI F2
treat_mat <- emmeans(SEMI_F2, ~ TREATMENT:MAT_STAT)
emm_data <- as.data.frame(treat_mat)
emm_data$TRANS_age <- emm_data$emmean
emm_data$HPD_low <- emm_data$lower.HPD
emm_data$HPD_high <- emm_data$upper.HPD
emm_data$TREATMENT <- factor(emm_data$TREATMENT, levels = c("DC", "SC", 'DT'))
emm_data$MAT_STAT <- factor(emm_data$MAT_STAT, levels = c("DOM", "SUB"))

milestones_org <- SEMI_data %>%
  group_by(ID, TREATMENT, SEX, MAT_STAT) %>%
  summarise(TRANS_age = median(Milestone_Age, na.rm = TRUE), .groups = "drop")

# 700* 500: SEMI_F2_TM
ggplot(emm_data, aes(x = TREATMENT, y = TRANS_age, color = TREATMENT, shape = MAT_STAT)) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +  # Show estimates as points
  geom_errorbar(aes(ymin = HPD_low, ymax = HPD_high), position = position_dodge(width = 0.7), width = 0.5) +
  scale_color_okabe_ito(order= c(2, 1, 7), name = "Grandmaternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_shape_discrete(name='Maternal status', labels = c('dominant', 'subordinate'))+
  scale_y_continuous(n.breaks=10)+
  scale_x_discrete(labels= c('DC', 'SC', 'DT'))+
  labs(x = "Grandmaternal\ntreatment", y = "Age (days)\n") +
  geom_point(data = milestones_org, aes(x = TREATMENT, y = TRANS_age, color=TREATMENT, shape = MAT_STAT), 
             position = position_dodge(width = 0.7), size = 3, alpha = 0.3)  +  
  theme_clean()

rm(treat_mat, emm_data)


################################################################################
## DIG -> CC ####

# plot(density(FULL_data$Milestone_Age),
#      xlab = "Transition ages",
#      ylab = "Density")

# range(FULL_data$Milestone_Age) 
# mean(FULL_data$Milestone_Age) 
# sd(FULL_data$Milestone_Age)

priors_FULL_f2 <- c(
  set_prior("normal(115, 15)", class = "Intercept"),
  set_prior("normal(0, 15)", class = "b"), # +- 30 days
  set_prior("student_t(3, 0, 10)", class = "sd") # ID
)

FULL_F2 <- brm(formula = Milestone_Age ~  TREATMENT * MAT_STAT * SEX + (1|ID),
               data = FULL_data,
               family = student(link='identity'),
               chains = 4, iter = 5000, warmup = 1500, seed = 42234223, 
               control = list(adapt_delta = 0.99, max_treedepth = 15),
               cores=4, 
               backend = 'cmdstanr', init= 'random',
               prior = priors_FULL_f2, #threads = threading(4),
               file="FULL_F2"
)

rm(priors_FULL, priors_FULL_f2)

# Results: FULL ####
FULL_F2 <- readRDS("FULL_F2.rds")

summary(FULL_F2)


plot(FULL_F2)
pp_check(FULL_F2, ndraws=100)

loo_R2(FULL_F2, moment_match = T)

bayes_R2(FULL_F2)


posterior_desc <- describe_posterior(
  FULL_F2,
  effects = "all", #"fixed", 
  component = "all",
  rope_range = rope_range(FULL_F2),  
  test = c("p_direction", "p_significance"),
  centrality = "all",
  dispersion = TRUE
)

# drop sigma
posterior_desc <- posterior_desc[-c(1, 13:112),]

# clean up labels:
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'b_', '')
posterior_desc$TREATMENT <- ifelse(grepl("TREATMENTSC", posterior_desc$Parameter), "SC",
                                   ifelse(grepl("TREATMENTDT", posterior_desc$Parameter), "DT", "DC"))
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTSC', 'SC')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTDT', 'DT')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'SEXM', 'Male')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'MAT_STATSUB', 'SUB mother')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'MAT_STATDOM', 'DOM mother')

custom_order <- c('DT:SUB mother:Male', 'SC:SUB mother:Male', "SUB mother:Male",
                  'DT:DOM mother:Male', 'SC:DOM mother:Male', "DOM mother:Male", 
                  'DT:Male', 'SC:Male', "Male", 
                  'DT:SUB mother','SC:SUB mother','SUB mother',
                  'DT:DOM mother','SC:DOM mother','DOM mother',
                  'DT','SC')

posterior_desc$Parameter <- factor(posterior_desc$Parameter, levels = custom_order)
posterior_desc$TREATMENT <- factor(posterior_desc$TREATMENT, levels = c("DC", "SC", "DT"))

# Coef_FULL_F2 700*500
ggplot(posterior_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 7), name = "F0 treatment", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
  theme_clean()

### FULL F2 EMMs####
(treat_mat<- emmeans(FULL_F2, pairwise ~ TREATMENT:MAT_STAT))

pd(treat_mat)

p_significance(treat_mat, threshold = rope_range(FULL_F2))

rm(treat_mat)

treat_mat <- emmeans(FULL_F2, ~ TREATMENT:MAT_STAT)
emm_data <- as.data.frame(treat_mat)
emm_data$TRANS_age <- emm_data$emmean
emm_data$HPD_low <- emm_data$lower.HPD
emm_data$HPD_high <- emm_data$upper.HPD
emm_data$TREATMENT <- factor(emm_data$TREATMENT, levels = c("DC", "SC", 'DT'))

milestones_org <- FULL_data %>%
  group_by(ID, TREATMENT, SEX, MAT_STAT) %>%
  summarise(TRANS_age = median(Milestone_Age, na.rm = TRUE), .groups = "drop")

# 700* 500: FULL_F2_TM
ggplot(emm_data, aes(x = TREATMENT, y = TRANS_age, color = TREATMENT, shape = MAT_STAT)) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +  # Show estimates as points
  geom_errorbar(aes(ymin = HPD_low, ymax = HPD_high), position = position_dodge(width = 0.7), width = 0.5) +
  scale_color_okabe_ito(order= c(2, 1, 7), name = "Grandmaternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_shape_discrete(name='Maternal status', labels = c('dominant', 'subordinate'))+
  scale_y_continuous(n.breaks=10)+
  scale_x_discrete(labels= c('DC', 'SC', 'DT'))+
  labs(x = "Grandmaternal\ntreatment", y = "Age (days)\n") +
  geom_point(data = milestones_org, aes(x = TREATMENT, y = TRANS_age, color=TREATMENT, shape = MAT_STAT), 
             position = position_dodge(width = 0.7), size = 3, alpha = 0.3)  +  
  theme_clean()

rm(treat_mat, emm_data, milestones_org)

### Peak DIG

# plot(density(DIG_PEAK_data$Milestone_Age),
#      xlab = "Transition ages",
#      ylab = "Density")

# range(DIG_PEAK_data$Milestone_Age) 
# mean(DIG_PEAK_data$Milestone_Age) 
# sd(DIG_PEAK_data$Milestone_Age)

priors_PEAK_f2 <- c(
  set_prior("normal(95, 10)", class = "Intercept"),
  set_prior("normal(0, 15)", class = "b"), # +- 30 days
  set_prior("student_t(3, 0, 10)", class = "sd") # ID
)

PEAK_F2 <- brm(formula = Milestone_Age ~ TREATMENT * MAT_STAT * SEX + (1|ID),
               data = DIG_PEAK_data,
               family = student(link='identity'),
               chains = 4, iter = 5000, warmup = 1500, seed = 42234223, 
               control = list(adapt_delta = 0.99, max_treedepth = 15),
               cores=4, 
               backend = 'cmdstanr', init= 'random',
               prior = priors_PEAK_f2, #threads = threading(4),
               file="PEAK_F2"
)

# Results: PEAK ####
PEAK_F2 <- readRDS("PEAK_F2.rds")

summary(PEAK_F2)

plot(PEAK_F2)
pp_check(PEAK_F2, ndraws=100)

loo_R2(PEAK_F2, moment_match = T)

bayes_R2(PEAK_F2)

posterior_desc <- describe_posterior(
  PEAK_F2,
  effects = "all", #"fixed", 
  component = "all",
  rope_range = rope_range(PEAK_F2),  
  test = c("p_direction", "p_significance"),
  centrality = "all",
  dispersion = TRUE
)

# drop sigma
posterior_desc <- posterior_desc[-c(1, 13:112),]

# clean up labels:
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'b_', '')
posterior_desc$TREATMENT <- ifelse(grepl("TREATMENTSC", posterior_desc$Parameter), "SC",
                                   ifelse(grepl("TREATMENTDT", posterior_desc$Parameter), "DT", "DC"))
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTSC', 'SC')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTDT', 'DT')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'SEXM', 'Male')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'MAT_STATSUB', 'SUB mother')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'MAT_STATDOM', 'DOM mother')

custom_order <- c('DT:SUB mother:Male', 'SC:SUB mother:Male', "SUB mother:Male",
                  'DT:DOM mother:Male', 'SC:DOM mother:Male', "DOM mother:Male", 
                  'DT:Male', 'SC:Male', "Male", 
                  'DT:SUB mother','SC:SUB mother','SUB mother',
                  'DT:DOM mother','SC:DOM mother','DOM mother',
                  'DT','SC')

posterior_desc$Parameter <- factor(posterior_desc$Parameter, levels = custom_order)
posterior_desc$TREATMENT <- factor(posterior_desc$TREATMENT, levels = c("DC", "SC", "DT"))

# Coef_PEAK_F2 700*500
ggplot(posterior_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 7), name = "F0 treatment", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
  theme_clean()

### PEAK F2 EMMs####
(treat_mat<- emmeans(PEAK_F2, pairwise ~ TREATMENT:MAT_STAT))


pd(treat_mat)


p_significance(treat_mat, threshold = rope_range(PEAK_F2))

rm(treat_mat)

treat_mat <- emmeans(PEAK_F2, ~ TREATMENT:MAT_STAT)
emm_data <- as.data.frame(treat_mat)
emm_data$TRANS_age <- emm_data$emmean
emm_data$HPD_low <- emm_data$lower.HPD
emm_data$HPD_high <- emm_data$upper.HPD
emm_data$TREATMENT <- factor(emm_data$TREATMENT, levels = c("DC", "SC", 'DT'))

milestones_org <- DIG_PEAK_data %>%
  group_by(ID, TREATMENT, SEX, MAT_STAT) %>%
  summarise(TRANS_age = median(Milestone_Age, na.rm = TRUE), .groups = "drop")

# 700* 500: PEAK_F2_TM
ggplot(emm_data, aes(x = TREATMENT, y = TRANS_age, color = TREATMENT, shape = MAT_STAT)) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +  # Show estimates as points
  geom_errorbar(aes(ymin = HPD_low, ymax = HPD_high), position = position_dodge(width = 0.7), width = 0.5) +
  scale_color_okabe_ito(order= c(2, 1, 7), name = "Grandmaternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_shape_discrete(name='Maternal status', labels = c('dominant', 'subordinate'))+
  scale_y_continuous(n.breaks=10)+
  scale_x_discrete(labels= c('DC', 'SC', 'DT'))+
  labs(x = "Grandmaternal\ntreatment", y = "Age (days)\n") +
  geom_point(data = milestones_org, aes(x = TREATMENT, y = TRANS_age, color=TREATMENT, shape = MAT_STAT), 
             position = position_dodge(width = 0.7), size = 3, alpha = 0.3)  +  
  theme_clean()

### cleanup ###
rm(list = ls(all.names = TRUE)) #includes hidden objects.
gc() 
