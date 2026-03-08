#load packages
library(tidyverse)
library(dplyr)

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


  
write_csv(sjdf_clean, "clean_data.csv")




