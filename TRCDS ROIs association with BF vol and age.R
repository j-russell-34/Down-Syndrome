library(tidyverse)
library(dplyr)
library(purrr)
library(glue)
library(readxl)
library(ggpp)
library(grid)
library(tibble)
library(rstatix)

#import table/excels to df
feobv_suvr_df <- read.table("DSCHOL_FEOBV_SUVR.txt",
                            sep="\t", header=FALSE)
feobv_ba_df <- read.table("DSCHOL_FEOBV_SUVR_brodmann_2024-02-16.txt",
                            sep="\t", header=FALSE)
dnseg_vol <-read.csv("DSCHOL_DnSeg_volumes.csv")
sclimbic <- read.csv("DSCHOL_SCLIMBIC_2024-02-16.csv")
samseg_df <- read.csv("DSCHOL_samseg.csv")

samseg_df <- samseg_df %>% rename("VISIT" = "SESSTYPE")

#apply headers to columns
colnames(feobv_suvr_df) <- c("SUBJECT", "REF_REGION", "PVC", "ROI", "SUVR")
colnames(feobv_ba_df) <- c("SUBJECT", "REF_REGION", "PVC", "ROI", "SUVR")

#change from long to wide
feobv_suvr_df_wide <- reshape(data = feobv_suvr_df,
                              idvar = "SUBJECT",
                              v.names = c("SUVR"),
                              timevar = "ROI",
                              direction = "wide")

feobv_ba_df_wide <- reshape(data = feobv_ba_df,
                              idvar = "SUBJECT",
                              v.names = c("SUVR"),
                              timevar = "ROI",
                              direction = "wide")

#convert - to . in column headers
names(feobv_suvr_df_wide) <- make.names(names(feobv_suvr_df_wide))
names(feobv_ba_df_wide) <- make.names(names(feobv_ba_df_wide))

#remove ROIs that are not of interest
feobv_df <- feobv_suvr_df_wide %>%
  select(contains(c("Amygdala", "Hippocampus", "Insula", "Caudate", "Cerebellum.Cortex", 
                     "Putamen", "Thalamus", "Pons", "Vermis", "Pons", "SUBJECT", "PVC")))
colnames(feobv_ba_df_wide) <- gsub("SUVR.Seg000", "BA", colnames(feobv_ba_df_wide))
colnames(feobv_ba_df_wide) <- gsub("SUVR.Seg00", "BA", colnames(feobv_ba_df_wide))


feobv_ba_df_wide <- subset(feobv_ba_df_wide, select = -c(BA0))

#add column to feobv_df for timepoint of scan relative to trc_ds BL
feobv_df["VISIT"] <- c("Baseline", "Baseline", "Month16","Baseline",
                                 "Baseline", "Baseline", "Baseline", "Baseline",
                                 "Baseline", "Baseline", "Baseline", "Baseline")
feobv_ba_df_wide["VISIT"] <- c("Baseline", "Baseline", "Month16","Baseline",
                       "Baseline", "Baseline", "Baseline", "Baseline",
                       "Baseline", "Baseline", "Baseline", "Baseline")

#average Left and Right volumes ROIs

feobv_avg_df <- feobv_df %>% select(contains(c("SUBJECT","PVC","VISIT")))

roi_list <- list("Amygdala", "Hippocampus", "Insula", "Caudate", "Cerebellum.Cortex", 
             "Putamen", "Thalamus", "Pons", "Vermis", "Pons")

for (roi in roi_list){
  feobv_avg_df[roi] <- rowMeans(select(feobv_df, matches(roi)))
}

#sum BF volumes through both BF vol methods

sum <- c("left", "right")
dnseg_vol["Dn_Seg_Total"] <- rowSums(select(dnseg_vol, any_of(sum)))
dnseg_tot <- dnseg_vol %>% select(c("SUBJECT", "Dn_Seg_Total"))

sum_sclimbic <- c("Left.Basal.Forebrain", "Right.Basal.Forebrain")
sclimbic["bf_total"] <- rowSums(select(sclimbic, any_of(sum_sclimbic)))

sclimbic_tot <- sclimbic %>% select(c("SUBJECT", "SESSTYPE", "bf_total")) 
sclimbic_tot <- sclimbic_tot %>% rename("VISIT" = "SESSTYPE")

#create 4 data frames for each BF method and cortical and subcortical regions
feobv_sc_dnseg <- merge(feobv_avg_df,dnseg_tot, by="SUBJECT")
feobv_sc_sclimbic <- merge(feobv_avg_df, sclimbic_tot, by=c("SUBJECT", "VISIT"))

feobv_ba_sclimbic <- merge(feobv_ba_df_wide,sclimbic_tot, by=c("SUBJECT", "VISIT"))
feobv_ba_dnseg <- merge(feobv_ba_df_wide,dnseg_tot, by="SUBJECT")


#add total vol from samseg
feobv_sc_dnseg <- merge(feobv_sc_dnseg,samseg_df %>% select(c("SUBJECT", "VISIT", "samseg_sbtiv")), by=c("SUBJECT", "VISIT"))
feobv_sc_sclimbic <- merge(feobv_sc_sclimbic,samseg_df %>% select(c("SUBJECT", "VISIT", "samseg_sbtiv")), by=c("SUBJECT", "VISIT"))
feobv_ba_sclimbic <- merge(feobv_ba_sclimbic,samseg_df %>% select(c("SUBJECT", "VISIT", "samseg_sbtiv")), by=c("SUBJECT", "VISIT"))
feobv_ba_dnseg <- merge(feobv_ba_dnseg,samseg_df %>% select(c("SUBJECT", "VISIT", "samseg_sbtiv")), by=c("SUBJECT", "VISIT"))



#generate blank lists to add correlation analysis too
correlation_pvalue_dnseg <- list()
correlation_beta_dnseg <- list()
ROI_dnseg <- list()
correlation_r2_dnseg <- list()

sc_rois <- list("Amygdala", "Hippocampus", "Insula", "Caudate", "Cerebellum.Cortex", 
                "Putamen", "Thalamus", "Pons", "Vermis", "Pons")

#perform linear models to dnseg vol for sc rois
for (rois in sc_rois){
  model <- lm(feobv_sc_dnseg[,c(rois)] ~ Dn_Seg_Total + samseg_sbtiv, data=feobv_sc_dnseg)
  correlation_beta_dnseg <- append(correlation_beta_dnseg, list(summary(model)$coefficient[2]))
  correlation_pvalue_dnseg <- append(correlation_pvalue_dnseg, list(summary(model)$coefficient[11]))
  ROI_dnseg <- append(ROI_dnseg, list(rois))
  correlation_r2_dnseg <- append(correlation_r2_dnseg, list(summary(model)$r.squared))
}


#perform linear models to dnseg vol for ba rois
for (i in 5:45){
  model <- lm(feobv_ba_dnseg[,c(i)] ~ Dn_Seg_Total + samseg_sbtiv, data=feobv_ba_dnseg)
  correlation_beta_dnseg <- append(correlation_beta_dnseg, list(summary(model)$coefficient[2]))
  correlation_pvalue_dnseg <- append(correlation_pvalue_dnseg, list(summary(model)$coefficient[11]))
  ROI_dnseg <- append(ROI_dnseg, list(colnames(feobv_ba_dnseg)[i]))
  correlation_r2_dnseg <- append(correlation_r2_dnseg, list(summary(model)$r.squared))
}

#combine 3 lists to individual list and convert to single dataframe for correlations
df_list_dnseg <- list(ROI_dnseg, correlation_beta_dnseg, correlation_r2_dnseg, correlation_pvalue_dnseg)

dn_seg_lm_res <- as.data.frame(do.call(cbind, df_list_dnseg))

#apply headers to df
colnames(dn_seg_lm_res)[1] <- "ROI"
colnames(dn_seg_lm_res)[2] <- "Beta"
colnames(dn_seg_lm_res)[3] <-"R-squared"
colnames(dn_seg_lm_res)[4] <- "pvalue"



correlation_pvalue_sclimbic <- list()
correlation_beta_sclimbic <- list()
ROI_sclimbic <- list()
correlation_r2_sclimbic <- list()

#perform linear models to sclimbic vol for sc rois
for (rois in sc_rois){
  model <- lm(feobv_sc_sclimbic[,c(rois)] ~ bf_total + samseg_sbtiv, data=feobv_sc_sclimbic)
  correlation_beta_sclimbic <- append(correlation_beta_sclimbic, list(summary(model)$coefficient[2]))
  correlation_pvalue_sclimbic <- append(correlation_pvalue_sclimbic, list(summary(model)$coefficient[11]))
  ROI_sclimbic <- append(ROI_sclimbic, list(rois))
  correlation_r2_sclimbic <- append(correlation_r2_sclimbic, list(summary(model)$r.squared))
}


#perform linear models to sclimbic vol for ba rois
for (i in 5:45){
  model <- lm(feobv_ba_sclimbic[,c(i)] ~ bf_total + samseg_sbtiv, data=feobv_ba_sclimbic)
  correlation_beta_sclimbic <- append(correlation_beta_sclimbic, list(summary(model)$coefficient[2]))
  correlation_pvalue_sclimbic <- append(correlation_pvalue_sclimbic, list(summary(model)$coefficient[11]))
  ROI_sclimbic <- append(ROI_sclimbic, list(colnames(feobv_ba_sclimbic)[i]))
  correlation_r2_sclimbic <- append(correlation_r2_sclimbic, list(summary(model)$r.squared))
}

#combine 3 lists to individual list and convert to single dataframe for correlations
df_list_sclimbic <- list(ROI_sclimbic, correlation_beta_sclimbic, 
                         correlation_r2_sclimbic, correlation_pvalue_sclimbic)

sclimbic_lm_res <- as.data.frame(do.call(cbind, df_list_sclimbic))

#apply headers to df
colnames(sclimbic_lm_res)[1] <- "ROI"
colnames(sclimbic_lm_res)[2] <- "Beta"
colnames(sclimbic_lm_res)[3] <- "R-squared"
colnames(sclimbic_lm_res)[4] <- "pvalue"

#add age 
demo_df <- read.csv("DSChol_dems.csv")

feobv_sc_dnseg <- merge(feobv_sc_dnseg, demo_df, by="SUBJECT")
feobv_sc_sclimbic <- merge(feobv_sc_sclimbic,demo_df, by="SUBJECT")
feobv_ba_sclimbic <- merge(feobv_ba_sclimbic,demo_df, by="SUBJECT")
feobv_ba_dnseg <- merge(feobv_ba_dnseg,demo_df, by="SUBJECT")


#generate blank lists to add correlation analysis too
correlation_pvalue_age <- list()
correlation_beta_age <- list()
ROI_age <- list()
r2_age <- list()
intercept_age <-list()

sc_rois <- list("Amygdala", "Hippocampus", "Insula", "Caudate", "Cerebellum.Cortex", 
                "Putamen", "Thalamus", "Pons", "Vermis", "Pons")

#perform linear models to age for sc rois
for (rois in sc_rois){
  model <- lm(feobv_sc_dnseg[,c(rois)] ~ age, data=feobv_sc_dnseg)
  correlation_beta_age <- append(correlation_beta_age, list(summary(model)$coefficient[2]))
  correlation_pvalue_age <- append(correlation_pvalue_age, list(summary(model)$coefficient[8]))
  ROI_age <- append(ROI_age, list(rois))
  r2_age <- append(r2_age, list(summary(model)$r.squared))
  intercept_age <- append(intercept_age, list(summary(model)$coefficient[1]))
}


#perform linear models to age for ba rois
for (i in 5:45){
  model <- lm(feobv_ba_dnseg[,c(i)] ~ age, data=feobv_ba_dnseg)
  correlation_beta_age <- append(correlation_beta_age, list(summary(model)$coefficient[2]))
  correlation_pvalue_age <- append(correlation_pvalue_age, list(summary(model)$coefficient[8]))
  ROI_age <- append(ROI_age, list(colnames(feobv_ba_dnseg)[i]))
  r2_age <- append(r2_age, list(summary(model)$r.squared))
  intercept_age <- append(intercept_age, list(summary(model)$coefficient[1]))
}

#combine 3 lists to individual list and convert to single dataframe for correlations
df_list_age <- list(ROI_age, intercept_age, correlation_beta_age, r2_age, 
                    correlation_pvalue_age)

age_lm_res <- as.data.frame(do.call(cbind, df_list_age))

#apply headers to df
colnames(age_lm_res)[1] <- "ROI"
colnames(age_lm_res)[2] <- "Intercept"
colnames(age_lm_res)[3] <- "Beta"
colnames(age_lm_res)[4] <- "R-squared"
colnames(age_lm_res)[5] <- "pvalue"


#Associations with BF volume controling for age
#generate blank lists to add correlation analysis too
correlation_pvalue_dnseg_age <- list()
correlation_beta_dnseg_age <- list()
ROI_dnseg_age <- list()
r2_dnseg_age <-list()

#perform linear models to dnseg vol for sc rois
for (rois in sc_rois){
  model <- lm(feobv_sc_dnseg[,c(rois)] ~ Dn_Seg_Total + samseg_sbtiv + age, data=feobv_sc_dnseg)
  correlation_beta_dnseg_age <- append(correlation_beta_dnseg_age, list(summary(model)$coefficient[2]))
  correlation_pvalue_dnseg_age <- append(correlation_pvalue_dnseg_age, list(summary(model)$coefficient[14]))
  ROI_dnseg_age <- append(ROI_dnseg_age, list(rois))
  r2_dnseg_age <- append(r2_dnseg_age, list(summary(model)$r.squared))
}


#perform linear models to dnseg vol for ba rois
for (i in 5:45){
  model <- lm(feobv_ba_dnseg[,c(i)] ~ Dn_Seg_Total + samseg_sbtiv + age, data=feobv_ba_dnseg)
  correlation_beta_dnseg_age <- append(correlation_beta_dnseg_age, list(summary(model)$coefficient[2]))
  correlation_pvalue_dnseg_age <- append(correlation_pvalue_dnseg_age, list(summary(model)$coefficient[14]))
  ROI_dnseg_age <- append(ROI_dnseg_age, list(colnames(feobv_ba_dnseg)[i]))
  r2_dnseg_age <- append(r2_dnseg_age, list(summary(model)$r.squared))
}

#combine 3 lists to individual list and convert to single dataframe for correlations
df_list_dnseg_age <- list(ROI_dnseg_age, correlation_beta_dnseg_age, 
                          r2_dnseg_age, correlation_pvalue_dnseg_age)

dn_seg_age_lm_res <- as.data.frame(do.call(cbind, df_list_dnseg_age))

#apply headers to df
colnames(dn_seg_age_lm_res)[1] <- "ROI"
colnames(dn_seg_age_lm_res)[2] <- "Beta"
colnames(dn_seg_age_lm_res)[3] <- "R-squared"
colnames(dn_seg_age_lm_res)[4] <- "pvalue"



correlation_pvalue_sclimbic_age <- list()
correlation_beta_sclimbic_age <- list()
ROI_sclimbic_age <- list()
r2_sclimbic_age <-list()

#perform linear models to sclimbic vol for sc rois
for (rois in sc_rois){
  model <- lm(feobv_sc_sclimbic[,c(rois)] ~ bf_total + samseg_sbtiv + age, data=feobv_sc_sclimbic)
  correlation_beta_sclimbic_age <- append(correlation_beta_sclimbic_age, list(summary(model)$coefficient[2]))
  correlation_pvalue_sclimbic_age <- append(correlation_pvalue_sclimbic_age, list(summary(model)$coefficient[14]))
  ROI_sclimbic_age <- append(ROI_sclimbic_age, list(rois))
  r2_sclimbic_age <- append(r2_sclimbic_age, list(summary(model)$r.squared))
}


#perform linear models to sclimbic vol for ba rois
for (i in 5:45){
  model <- lm(feobv_ba_sclimbic[,c(i)] ~ bf_total + samseg_sbtiv +age, data=feobv_ba_sclimbic)
  correlation_beta_sclimbic_age <- append(correlation_beta_sclimbic_age, list(summary(model)$coefficient[2]))
  correlation_pvalue_sclimbic_age <- append(correlation_pvalue_sclimbic_age, list(summary(model)$coefficient[14]))
  ROI_sclimbic_age <- append(ROI_sclimbic_age, list(colnames(feobv_ba_sclimbic)[i]))
  r2_sclimbic_age <- append(r2_sclimbic_age, list(summary(model)$r.squared))
}

#combine 3 lists to individual list and convert to single dataframe for correlations
df_list_sclimbic_age <- list(ROI_sclimbic_age, correlation_beta_sclimbic_age,
                             r2_sclimbic_age, correlation_pvalue_sclimbic_age)

sclimbic_lm_res_age <- as.data.frame(do.call(cbind, df_list_sclimbic_age))

#apply headers to df
colnames(sclimbic_lm_res_age)[1] <- "ROI"
colnames(sclimbic_lm_res_age)[2] <- "Beta"
colnames(sclimbic_lm_res_age)[3] <- "R-squared"
colnames(sclimbic_lm_res_age)[4] <- "pvalue"

#corrected p-values
age_pvalues <- as.list(age_lm_res[["pvalue"]])
age_lm_res["adj_pvalue_fdr"] <- p.adjust(age_pvalues, method="fdr")

age_lm_res <- apply(age_lm_res,2,as.character)
write.csv(age_lm_res, "age_FEOBV_association.csv", row.names = FALSE)

sclimbic_pvalues <- as.list(sclimbic_lm_res[["pvalue"]])
sclimbic_lm_res["adj_pvalue_fdr"] <- p.adjust(sclimbic_pvalues, method="fdr")

sclimbic_lm_res <- apply(sclimbic_lm_res,2,as.character)
write.csv(sclimbic_lm_res, "sclimbic_FEOBV_association.csv", row.names = FALSE)

sclimbic_age_pvalues <- as.list(sclimbic_lm_res_age[["pvalue"]])
sclimbic_lm_res_age["adj_pvalue_fdr"] <- p.adjust(sclimbic_age_pvalues, method="fdr")

sclimbic_lm_res_age <- apply(sclimbic_lm_res_age,2,as.character)
write.csv(sclimbic_lm_res_age, "sclimbic_FEOBV_association_age_corr.csv", row.names = FALSE)

dnseg_pvalues <- as.list(dn_seg_lm_res[["pvalue"]])
dn_seg_lm_res["adj_pvalue_fdr"] <- p.adjust(dnseg_pvalues, method="fdr")

dn_seg_lm_res <- apply(dn_seg_lm_res,2,as.character)
write.csv(dn_seg_lm_res, "dnseg_FEOBV_association.csv", row.names = FALSE)

dnseg_age_pvalues <- as.list(dn_seg_age_lm_res[["pvalue"]])
dn_seg_age_lm_res["adj_pvalue_fdr"] <- p.adjust(dnseg_age_pvalues, method="fdr")

dn_seg_age_lm_res <- apply(dn_seg_age_lm_res,2,as.character)
write.csv(dn_seg_age_lm_res, "dnseg_FEOBV_association_age_corr.csv", row.names = FALSE)
