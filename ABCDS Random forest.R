library(tidyverse)
library(dplyr)
library(purrr)
library(glue)
library(readxl)
library(ggpp)
library(grid)
library(party)
library(varImp)

dnsegdata_df <- read.csv("ABCDS_DnSeg_volumes.csv")
sclimbicdata_df <- read.csv("ABCDS_sclimbic_volumes.csv")
samseg_df <- read.csv("ABCDS_SAMSEG.csv")
dnseg_abcds_df <- read.csv("ABC-DS DnSeg DF for R 20231221.csv")
sclimbic_abcds_df <- read.csv("ABC-DS ScLimbic DF for R 20231221.csv")

sclimbic_abcds_df <- merge(sclimbic_abcds_df, sclimbicdata_df, by="SUBJECT")
sclimbic_abcds_df <- merge(sclimbic_abcds_df, samseg_df, by="SUBJECT")
sclimbic_abcds_df <- sclimbic_abcds_df[!duplicated(sclimbic_abcds_df$SUBJECT),]

sclimbic_abcds_df$totalBFavg <- (sclimbic_abcds_df$Right.Basal.Forebrain + 
  sclimbic_abcds_df$Left.Basal.Forebrain)/2

sclimbic_abcds_df$Sex <- factor(sclimbic_abcds_df$Gender)
sclimbic_abcds_df$avgchbfdivtiv <- sclimbic_abcds_df$totalBFavg/sclimbic_abcds_df$samseg_sbtiv

set.seed(123)
data.controls <- cforest_unbiased(ntree=2000, mtry=3)

sclimbic_data_forest <- cforest(totalBFavg ~ Amyloid..centiloids. + Age + Sex + samseg_sbtiv 
        + samseg_lesions, data=sclimbic_abcds_df, controls=data.controls)

sclimbictiv_data_forest <- cforest(avgchbfdivtiv ~ Amyloid..centiloids. + Age + Sex 
                                + samseg_lesions, data=sclimbic_abcds_df, controls=data.controls)

dnseg_abcds_df <- merge(dnseg_abcds_df, dnsegdata_df, by="SUBJECT")
dnseg_abcds_df <- merge(dnseg_abcds_df, samseg_df, by="SUBJECT")
dnseg_abcds_df <- dnseg_abcds_df[!duplicated(dnseg_abcds_df$SUBJECT),]

dnseg_abcds_df$totalch4avg <- (dnseg_abcds_df$Left_DnSeg + 
                                 dnseg_abcds_df$Left_DnSeg)/2

dnseg_abcds_df$Sex <- factor(dnseg_abcds_df$Gender)
dnseg_abcds_df$avgch4divtiv <- dnseg_abcds_df$totalch4avg/dnseg_abcds_df$samseg_sbtiv

data.controls <- cforest_unbiased(ntree=2000, mtry=3)

dnseg_data_forest <- cforest(totalch4avg ~ Amyloid..centiloids. + Age + Sex + samseg_sbtiv 
                                + samseg_lesions, data=dnseg_abcds_df, controls=data.controls)

dnsegtiv_data_forest <- cforest(avgch4divtiv ~ Amyloid..centiloids. + Age + Sex 
                             + samseg_lesions, data=dnseg_abcds_df, controls=data.controls)

library("writexl")
write_xlsx(dnseg_abcds_df, "dnseg_abcds.xlsx")
write_xlsx(sclimbic_abcds_df, "sclimbic_abcds.xlsx")
