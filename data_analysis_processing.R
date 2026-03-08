#load packages
library(tidyverse)
library(ggplot2)
library(ggeffects)
library(glmmTMB)
library(DHARMa)
library(performance)
library(lme4)

#clean data

raw_data <- read.csv("Jay_data - DATA.csv") #import data

sjdf_clean <- raw_data %>% #remove stellar's jay
  filter(SPECIES != "STJA")

sjdf_clean <- sjdf_clean |> #change to 1s/0s
  mutate(MOB = if_else(MOB == "YES", 1, 0)) |>
  mutate(ALARM = if_else(ALARM == "YES", 1, 0)) |>
  mutate(FLEE = if_else(FLEE == "YES", 1, 0)) |>
  mutate(INTEREST = if_else(INTEREST == "YES", 1, 0))
  
sjdf_clean <- sjdf_clean |> #remove extraneous columns
  select(-ALARMS.GS, -TOTAL.GS, -TOTAL.VOCALIZATIONS, -GPS, 
          -DISTANCE..m., -HEIGHT..m., -NOTES, -CALLS.GS)

sjdf_clean <- sjdf_clean |> #rename latency
  rename(LATENCY.ALARM = LATENCY.ALARM..s., NUMBER.MOB = MOBS)

sjdf_clean <- sjdf_clean %>% #replace missing values with NAs
  mutate(across(everything(), as.character)) |> #make all columns characters
  mutate(across(everything(), ~na_if(.x, ""))) |>
  mutate(across(everything(), ~na_if(.x, "-"))) |>
  mutate(across(everything(), ~na_if(.x, "--"))) |>
  type.convert(as.is = TRUE)

sjdf_clean <- sjdf_clean |> #fix class of variables
  mutate(across(c(ALARM, MOB, FLEE, INTEREST, HYPOTENUSE, LATENCY.ALARM, 
                  GROUP.SIZE, NUMBER.MOB, NUMBER.ALARM), as.numeric)) |>
  mutate(across(c(SPECIES, SEASON, PLAYBACK, SITE, SUBSITE, TREATMENT), as.factor))
  


### data exploration






### data analysis

sjdf_clean$TREATMENT <- relevel(sjdf_clean$TREATMENT, ref = "HAWK")#make hawk
#the reference level for model estimates
sjdf_clean$SPECIES <- relevel(sjdf_clean$SPECIES, ref = "ISSJ")

sjdf_clean <- sjdf_clean |> #mutate the group size and hypotenuse variables to 
  #standardize
  mutate(GROUP.SIZE = as.numeric(scale(GROUP.SIZE))) |>
  mutate(HYPOTENUSE = as.numeric(scale(HYPOTENUSE)))


#alarm occurrence model
m.alarm <- glmer(ALARM ~  SPECIES + TREATMENT + GROUP.SIZE + HYPOTENUSE +
                 SEASON +
                 (1 | SUBSITE),
               data = sjdf_clean,
               family = binomial)  

check_collinearity(m.alarm) #check for colinearity between predictors

sim.m.alarm <- simulateResiduals(fittedModel = m.alarm) #check model fit
testDispersion(sim.m.alarm)
testZeroInflation(sim.m.alarm)
testOutliers(sim.m.alarm) 
  
summary(m.alarm) # get model output
exp(fixef(m.alarm)) #get odds-ratios

# Prediction grid
pred_data_alarm <- data.frame(
  SPECIES   = c("ISSJ", "CASJ"),   
  TREATMENT = "HAWK",              
  SEASON    = "SUMMER",          
  PLAYBACK  = "YES",               
  GROUP.SIZE   = 0,                 # scaled mean
  HYPOTENUSE  = 0,                  # scaled mean
  SUBSITE = NA                          # random effect ignored
)

# Predicted log-odds and standard errors
pred_link_alarm <- predict(m.alarm, newdata = pred_data_alarm, type = "link", 
                     se.fit = TRUE, re.form = NA)

# Convert log-odds ± 1.96*SE to probabilities
pred_data_alarm$pred_prob <- plogis(pred_link_alarm$fit)  # predicted probability
pred_data_alarm$lower <- plogis(pred_link_alarm$fit - 1.96 * pred_link_alarm$se.fit)  # lower 95% CI
pred_data_alarm$upper <- plogis(pred_link_alarm$fit + 1.96 * pred_link_alarm$se.fit)  # upper 95% CI




#mob occurrence model
m.mob <- glmer(MOB ~  SPECIES + TREATMENT + GROUP.SIZE + HYPOTENUSE +
                   SEASON +
                   (1 | SUBSITE),
                 data = sjdf_clean,
                 family = binomial)  

check_collinearity(m.mob) #check for colinearity between predictors

sim.m.mob <- simulateResiduals(fittedModel = m.mob) #check model fit
testDispersion(sim.m.mob)
testZeroInflation(sim.m.mob)
testOutliers(sim.m.mob) 

summary(m.mob) # get model output
exp(fixef(m.mob)) #get odds-ratios (convert estimates to interpret)


# Prediction grid
pred_data_mob <- data.frame(
  SPECIES   = c("ISSJ", "CASJ"),   
  TREATMENT = "HAWK",              
  SEASON    = "SUMMER",          
  PLAYBACK  = "YES",               
  GROUP.SIZE   = 0,                 # scaled mean
  HYPOTENUSE  = 0,                  # scaled mean
  SUBSITE = NA                          # random effect ignored
)

# Predicted log-odds and standard errors
pred_link_mob <- predict(m.mob, newdata = pred_data_mob, type = "link", 
                     se.fit = TRUE, re.form = NA)

# Convert log-odds ± 1.96*SE to probabilities
pred_data_mob$pred_prob <- plogis(pred_link_mob$fit)  # predicted probability
pred_data_mob$lower <- plogis(pred_link_mob$fit - 1.96 * pred_link_mob$se.fit)  # lower 95% CI
pred_data_mob$upper <- plogis(pred_link_mob$fit + 1.96 * pred_link_mob$se.fit)  # upper 95% CI




#flee occurrence model
m.flee <- glmer(FLEE ~  SPECIES + TREATMENT + GROUP.SIZE + HYPOTENUSE +
                   SEASON +
                   (1 | SUBSITE),
                 data = sjdf_clean,
                 family = binomial)  

check_collinearity(m.flee) #check for colinearity between predictors

sim.m.flee <- simulateResiduals(fittedModel = m.flee) #check model fit
testDispersion(sim.m.flee)
testZeroInflation(sim.m.flee)
testOutliers(sim.m.flee) 

summary(m.flee) # get model output
exp(fixef(m.flee)) #get odds-ratios


# Prediction grid
pred_data_flee <- data.frame(
  SPECIES   = c("ISSJ", "CASJ"),   
  TREATMENT = "HAWK",              
  SEASON    = "SUMMER",          
  PLAYBACK  = "YES",               
  GROUP.SIZE   = 0,                 # scaled mean
  HYPOTENUSE  = 0,                  # scaled mean
  SUBSITE = NA                          # random effect ignored
)

# Predicted log-odds and standard errors
pred_link_flee <- predict(m.flee, newdata = pred_data_flee, type = "link", 
                         se.fit = TRUE, re.form = NA)

# Convert log-odds ± 1.96*SE to probabilities
pred_data_flee$pred_prob <- plogis(pred_link_flee$fit)  # predicted probability
pred_data_flee$lower <- plogis(pred_link_flee$fit - 1.96 * pred_link_flee$se.fit)  # lower 95% CI
pred_data_flee$upper <- plogis(pred_link_flee$fit + 1.96 * pred_link_flee$se.fit)  # upper 95% CI




#interest occurrence model
m.interest <- glmer(INTEREST ~  SPECIES + TREATMENT + GROUP.SIZE + HYPOTENUSE +
                   SEASON +
                   (1 | SUBSITE),
                 data = sjdf_clean,
                 family = binomial)  

check_collinearity(m.interest) #check for colinearity between predictors

sim.m.interest <- simulateResiduals(fittedModel = m.interest) #check model fit
testDispersion(sim.m.interest)
testZeroInflation(sim.m.interest)
testOutliers(sim.m.interest) 

summary(m.interest) # get model output
exp(fixef(m.interest)) #get odds-ratios




pred_data_interest <- data.frame(
  SPECIES   = c("ISSJ", "CASJ"),   
  TREATMENT = "HAWK",              
  SEASON    = "SUMMER",          
  PLAYBACK  = "YES",               
  GROUP.SIZE   = 0,                 # scaled mean
  HYPOTENUSE  = 0,                  # scaled mean
  SUBSITE = NA                          # random effect ignored
)

# Predicted log-odds and standard errors
pred_link_interest <- predict(m.interest, newdata = pred_data_interest, type = "link", 
                          se.fit = TRUE, re.form = NA)

# Convert log-odds ± 1.96*SE to probabilities
pred_data_interest$pred_prob <- plogis(pred_link_interest$fit) # predicted probability
pred_data_interest$lower <- plogis(pred_link_interest$fit - 1.96 * pred_link_interest$se.fit) # lower 95% CI
pred_data_interest$upper <- plogis(pred_link_interest$fit + 1.96 * pred_link_interest$se.fit)  # upper 95% CI

#make graph of all four models
pred_data_alarm$MODEL <- "Model 1"
pred_data_mob$MODEL <- "Model 2"
pred_data_interest$MODEL <- "Model 3"
pred_data_flee$MODEL <- "Model 4"

all_pred <- bind_rows(pred_data_alarm, pred_data_mob, pred_data_interest, pred_data_flee)

labels <- data.frame(
  MODEL = c("Model 1", "Model 2", "Model 3", "Model 4"),
  PANEL_LABEL = c("a) *", "b) *", "c) *", "d)")
)

all_pred <- merge(all_pred, labels, by = "MODEL")
ggplot(all_pred, aes(x = SPECIES, y = pred_prob, color = SPECIES)) +
  geom_point(size = 2) +
  geom_errorbar(data = subset(all_pred, MODEL != "Model 4"),  # exclude model 4
                aes(ymin = lower, ymax = upper),
                width = 0.2) +
  geom_text(
    aes(label = PANEL_LABEL),
    x = -Inf,  # far left
    y = 1.0,  # just above top of y-axis (adjust as needed)
    hjust = -.5, # nudges left
    vjust = 1,
    inherit.aes = FALSE) +
  facet_wrap(~ MODEL) +
  ylim(0,1) +
  labs(
    x = "Species",
    y = "Predicted Probability of Behaviour"
  ) +
  theme_classic() +
  theme(strip.text = element_blank(), 
        panel.border = element_rect(color = "black", fill = NA, size = 1))





