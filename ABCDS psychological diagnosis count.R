library(tidyverse)
library(dplyr)
library(purrr)
library(glue)
library(readxl)
library(ggpp)
library(grid)

hh_df <- read.csv("Health_History_13May2024.csv")

hh_df <- hh_df[hh_df$event_code =="bl",]

depression <- table(hh_df$hh_depression)

print(depression)

bipolardis <- table(hh_df$hh_bd)

print(bipolardis)

anxiety <- table(hh_df$hh_ocd)

print(anxiety)

schiz <- table(hh_df$hh_schiz)

print(schiz)

dementia <- table(hh_df$hh_dementia)

print(dementia)

parkinsons <- table(hh_df$hh_park)

print(parkinsons)

seiz <- table(hh_df$hh_seizure)

print(seiz)
