library(tidyverse)
library(dplyr)
library(purrr)
library(glue)
library(readxl)
library(ggpp)
library(grid)
library(tibble)
library(rstatix)

setwd("~/Study_data/Down Syndrome/TRCDS")

#import table/excels to df
feobv_suvr_df <- read.table("DSCHOL_FEOBV_SUVR.txt",
                            sep="\t", header=FALSE)
feobv_ba_df <- read.table("DSCHOL_FEOBV_SUVR_brodmann_2024-02-16.txt",
                          sep="\t", header=FALSE)
albin_df <- read.csv("albin2018.csv")
champ_suvr_df <- read.table("CHAMP_FEOBV_SUVR_2024-02-16.txt", sep='\t',
                            header=FALSE)
champ_ba_df <- read.table("CHAMP_FEOBV_SUVR_brodmann_2024-02-16.txt",
                          sep="\t", header=FALSE)

#apply headers to columns
colnames(feobv_suvr_df) <- c("SUBJECT", "REF_REGION", "PVC", "ROI", "SUVR")
colnames(feobv_ba_df) <- c("SUBJECT", "REF_REGION", "PVC", "ROI", "SUVR")
colnames(champ_suvr_df) <- c("SUBJECT", "REF_REGION", "PVC", "ROI", "SUVR")
colnames(champ_ba_df) <- c("SUBJECT", "REF_REGION", "PVC", "ROI", "SUVR")

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

champ_ba_df_wide <- reshape(data = champ_ba_df,
                              idvar = "SUBJECT",
                              v.names = c("SUVR"),
                              timevar = "ROI",
                              direction = "wide")
champ_suvr_df_wide <- reshape(data = champ_suvr_df,
                            idvar = "SUBJECT",
                            v.names = c("SUVR"),
                            timevar = "ROI",
                            direction = "wide")

#convert - to . in column headers
names(feobv_suvr_df_wide) <- make.names(names(feobv_suvr_df_wide))
names(feobv_ba_df_wide) <- make.names(names(feobv_ba_df_wide))
names(champ_ba_df_wide) <- make.names(names(champ_ba_df_wide))
names(champ_suvr_df_wide) <- make.names(names(champ_suvr_df_wide))

#remove ROIs that are not of interest
feobv_df <- feobv_suvr_df_wide %>%
  select(contains(c("Amygdala", "Hippocampus", "Insula", "Caudate", "Cerebellum.Cortex", 
                    "Putamen", "Thalamus", "Pons", "Vermis", "Pons", "SUBJECT", "PVC")))

champ_df <- champ_suvr_df_wide %>%
  select(contains(c("Amygdala", "Hippocampus", "Insula", "Caudate", "Cerebellum.Cortex", 
                    "Putamen", "Thalamus", "Pons", "Vermis", "Pons", "SUBJECT", "PVC")))

#add column to feobv_df for timepoint of scan relative to trc_ds BL
feobv_df["VISIT"] <- c("Baseline", "Baseline", "Month16","Baseline",
                       "Baseline", "Baseline", "Baseline", "Baseline",
                       "Baseline", "Baseline", "Baseline", "Baseline")
feobv_ba_df_wide["VISIT"] <- c("Baseline", "Baseline", "Month16","Baseline",
                               "Baseline", "Baseline", "Baseline", "Baseline",
                               "Baseline", "Baseline", "Baseline", "Baseline")

#average Left and Right volumes ROIs

feobv_avg_df <- feobv_df %>% select(contains(c("SUBJECT","PVC")))
champ_avg_df <- champ_df %>% select(contains(c("SUBJECT","PVC")))

roi_list <- list("Amygdala", "Hippocampus", "Insula", "Caudate", "Cerebellum.Cortex", 
                 "Putamen", "Thalamus", "Pons", "Vermis")

for (roi in roi_list){
  feobv_avg_df[roi] <- rowMeans(select(feobv_df, matches(roi)))
}

for (roi in roi_list){
  champ_avg_df[roi] <- rowMeans(select(champ_df, matches(roi)))
}


#subjectwise means of ROIS
subject_wise_roi_mean <- data.frame()
subject_wise_champ_roi_mean <- data.frame()

for (roi in roi_list){
  subject_wise_champ_roi_mean <- rbind(subject_wise_champ_roi_mean, 
                                 list(roi, mean(champ_avg_df[[roi]])))
}

for (roi in roi_list){
  subject_wise_roi_mean <- rbind(subject_wise_roi_mean, 
                                 list(roi, mean(feobv_avg_df[[roi]])))
}


#change names in BA df
colnames(feobv_ba_df_wide) <- gsub("SUVR.Seg000", "BA", colnames(feobv_ba_df_wide))
colnames(feobv_ba_df_wide) <- gsub("SUVR.Seg00", "BA", colnames(feobv_ba_df_wide))
colnames(champ_ba_df_wide) <- gsub("SUVR.Seg000", "BA", colnames(champ_ba_df_wide))
colnames(champ_ba_df_wide) <- gsub("SUVR.Seg00", "BA", colnames(champ_ba_df_wide))

feobv_ba_df_wide <- subset(feobv_ba_df_wide, select = -c(BA0))
champ_ba_df_wide <- subset(champ_ba_df_wide, select = -c(BA0))

#subjectwise means of BAs
for (ba in 4:44){
  subject_wise_roi_mean <- rbind(subject_wise_roi_mean, 
                                 list( colnames(feobv_ba_df_wide[ba]),
                                       mean(feobv_ba_df_wide[,c(ba)])))
}

for (ba in 4:44){
  subject_wise_champ_roi_mean <- rbind(subject_wise_champ_roi_mean, 
                                 list( colnames(champ_ba_df_wide[ba]),
                                       mean(champ_ba_df_wide[,c(ba)])))
}

colnames(subject_wise_champ_roi_mean) <- c("ROIs", "SUVR_CHAMP")
colnames(subject_wise_roi_mean) <- c("ROIs", "SUVR_DSCHOL")

#albin df change roi names to match DSCHOL
colnames(albin_df) <- gsub("ROI_BA", "BA", colnames(albin_df))
colnames(albin_df) <- gsub("_suvr", "", colnames(albin_df))

#albin df convert to tall
albin_df <- albin_df %>%
  pivot_longer(
    cols = "BA10":"Pons",
    names_to = "ROIs",
    values_to = "SUVR_albin"
  )

comparison_df <- merge(subject_wise_roi_mean, albin_df, all=TRUE, by="ROIs")
comparison_df <- merge(comparison_df, subject_wise_champ_roi_mean, all=TRUE, by ="ROIs")

albinvumcmodel <- lm (SUVR_DSCHOL ~ SUVR_albin,
                      data = comparison_df)
summary(albinvumcmodel)


#comment out annotation and remove annotation for paper figure
#grob <-grobTree(textGrob(expression(atop(beta == 1.11, p < "2x10"^-16)), x=0.8, y=0.2, hjust = 0,
 #                        gp=gpar(col="red", fontsize =14, fontface="bold")))

png(file="albin vs dschol.png", units="in", width=7, height=5, res=300)
ggplot(comparison_df,aes(y=SUVR_DSCHOL, x=SUVR_albin)) +
  geom_point(size=2, color="blue") + 
  geom_smooth(method='lm', color="black")+
  labs(y=expression(paste("VUMC DS Cohort ["^"18","F]-FEOBV SUVRs")),
       x= expression(paste("Albin (2018) published ["^"18","F]-FEOBV SUVRs")))+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=14), axis.title=element_text(size=16)) +
  expand_limits(x=c(0,10), y=c(0,10)) #+ annotation_custom(grob)
dev.off()

champdscholmodel <- lm (SUVR_DSCHOL ~ SUVR_CHAMP,
                      data = comparison_df)
summary(champdscholmodel)


# S.E.M. included as a column for CHAMP and DSCHOL in the comparison df
sem <- function(x) sd(x)/sqrt(length(x))

champ_sem_df <- data.frame()
feobv_sem_df <- data.frame()

for (ba in 4:44){
  champ_sem_df <- rbind(champ_sem_df, list( colnames(champ_ba_df_wide[ba]),
                              sem(champ_ba_df_wide[,c(ba)])))
}

for (ba in 4:44){
  feobv_sem_df <- rbind(feobv_sem_df, list( colnames(feobv_ba_df_wide[ba]),
                                            sem(feobv_ba_df_wide[,c(ba)])))
}

for (roi in roi_list){
 champ_sem_df <- rbind(champ_sem_df, list(roi, sem(champ_avg_df[[roi]])))
}

for (roi in roi_list){
  feobv_sem_df <- rbind(feobv_sem_df, list(roi, sem(feobv_avg_df[[roi]])))
}

colnames(feobv_sem_df) <- c("ROIs", "DSCHOL_SEM")
colnames(champ_sem_df) <- c("ROIs", "CHAMP_SEM")

comparison_df <- merge(comparison_df, feobv_sem_df, all=TRUE, by ="ROIs")
comparison_df <- merge(comparison_df, champ_sem_df, all=TRUE, by ="ROIs")

#Perform T-tests adding the t and p to the dataframe above
ttests <- data.frame()

for (ba in 4:44){
  ttestres <- t.test(feobv_ba_df_wide[,c(ba)], champ_ba_df_wide[,c(ba)],
                     alternative = "two.sided", paired= FALSE, 
                     var.equal = FALSE, conf.level =0.95)
  ttests <- rbind(ttests, list( colnames(feobv_ba_df_wide[ba]),
                                      ttestres[[1]], ttestres[3]))
}


feobv_df_clean <- feobv_df %>% select(contains(c("SUBJECT","PVC")))
champ_df_clean <- champ_df %>% select(contains(c("SUBJECT","PVC")))

for (roi in roi_list){
  feobv_df_clean[roi] <- select(feobv_avg_df, matches(roi))
}

for (roi in roi_list){
  champ_df_clean[roi] <- select(champ_avg_df, matches(roi))
}

#comment out annotation and remove annotation for paper figure
#grob <-grobTree(textGrob(expression(atop(beta == 1.15, p < "2x10"^-16)), x=0.8, y=0.2, hjust = 0,
#                       gp=gpar(col="red", fontsize =14, fontface="bold")))

png(file="champ vs dschol.png", units="in", width=7, height=5, res=300)
ggplot(comparison_df,aes(y=SUVR_DSCHOL, x=SUVR_CHAMP)) +
  geom_point(size=2, color="blue") + 
  geom_smooth(method='lm', color="black")+
  labs(y=expression(paste("VUMC DS Cohort ["^"18","F]-FEOBV SUVRs")), x= 
         expression(paste("VUMC Neurotypical Cohort ["^"18","F]-FEOBV SUVRs")))+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=14), axis.title=element_text(size=16)) +
  expand_limits(x=c(0,10), y=c(0,10)) #+ annotation_custom(grob)
dev.off()


for (roi in roi_list){
  ttestres <- t.test(feobv_df_clean[roi], champ_df_clean[roi],
                     alternative = "two.sided", paired= FALSE, 
                     var.equal = FALSE, conf.level =0.95)
  ttests <- rbind(ttests, list(roi, ttestres[[1]], ttestres[3]))
}

colnames(ttests) <- c("ROIs", "tvalue", "pvalue")

comparison_df <- merge(comparison_df, ttests, all=TRUE, by ="ROIs")

#correct p-values for multiple comparisons
pvalues <- as.list(comparison_df[["pvalue"]])

comparison_df["adj_pvalue_fdr"] <- p.adjust(pvalues, method="fdr")
comparison_df["adj_pvalue_bonferroni"]<- p.adjust(pvalues, method="bonferroni")

comparison_df <- apply(comparison_df,2,as.character)

write.csv(comparison_df, "albin_and_champ_comparison.csv", row.names = FALSE)

#sex comparison in DS cohort
ds_sex_df <- read.csv("DSChol_sex.csv")

feobv_sex_df <- merge(feobv_df_clean, ds_sex_df, all=TRUE, by ="SUBJECT")
feobv_ba_sex_df <- merge(feobv_ba_df_wide, ds_sex_df, all=TRUE, by ="SUBJECT")

feobv_males <- feobv_sex_df[feobv_sex_df$Sex == 'Male',]
feobv_females <- feobv_sex_df[feobv_sex_df$Sex == 'Female',]
feobv_males_ba <- feobv_ba_sex_df[feobv_ba_sex_df$Sex == 'Male',]
feobv_females_ba <- feobv_ba_sex_df[feobv_ba_sex_df$Sex == 'Female',]

feobv_male_avgs <- data.frame()
for (roi in roi_list){
  feobv_male_avgs <- rbind(feobv_male_avgs, list(roi, mean(feobv_males[[roi]])))
}


feobv_female_avgs <- data.frame()
for (roi in roi_list){
  feobv_female_avgs <- rbind(feobv_female_avgs, list(roi, mean(feobv_females[[roi]])))
}


for (ba in 4:44){
  feobv_male_avgs <- rbind(feobv_male_avgs, list(colnames(feobv_males_ba[ba]),
                              mean(feobv_males_ba[,c(ba)])))
}

colnames(feobv_male_avgs) <- c("ROIs", "SUVR_DS_Males")


for (ba in 4:44){
  feobv_female_avgs <- rbind(feobv_female_avgs, list(colnames(feobv_females_ba[ba]),
                                                       mean(feobv_females_ba[,c(ba)])))
}

colnames(feobv_female_avgs) <- c("ROIs", "SUVR_DS_Females")

comparison_sex_df <- merge(feobv_male_avgs, feobv_female_avgs, all=TRUE, by ="ROIs")

#s.e.m. calculations
feobv_sem_male_df <- data.frame()
feobv_sem_female_df <- data.frame()

for (ba in 4:44){
  feobv_sem_male_df <- rbind(feobv_sem_male_df, list( colnames(feobv_males_ba[ba]),
                                            sem(feobv_males_ba[,c(ba)])))
}

for (ba in 4:44){
  feobv_sem_female_df <- rbind(feobv_sem_female_df, list( colnames(feobv_females_ba[ba]),
                                            sem(feobv_females_ba[,c(ba)])))
}

for (roi in roi_list){
  feobv_sem_male_df <- rbind(feobv_sem_male_df, list(roi, sem(feobv_males[[roi]])))
}

for (roi in roi_list){
  feobv_sem_female_df <- rbind(feobv_sem_female_df, list(roi, sem(feobv_females[[roi]])))
}

colnames(feobv_sem_female_df) <- c("ROIs", "SUVR_DS_Females_SEM")
colnames(feobv_sem_male_df) <- c("ROIs", "SUVR_DS_Males_SEM")

comparison_sex_df <- merge(comparison_sex_df, feobv_sem_male_df, all=TRUE, by ="ROIs")
comparison_sex_df <- merge(comparison_sex_df, feobv_sem_female_df, all=TRUE, by ="ROIs")

ttestssex <- data.frame()

for (ba in 4:44){
  ttestressex <- t.test(feobv_males_ba[,c(ba)], feobv_females_ba[,c(ba)],
                     alternative = "two.sided", paired= FALSE, 
                     var.equal = FALSE, conf.level =0.95)
  ttestssex <- rbind(ttestssex, list( colnames(feobv_males_ba[ba]),
                                      ttestressex[[1]], ttestressex[3]))
}


for (roi in roi_list){
  ttestressex <- t.test(feobv_males[roi], feobv_females[roi],
                     alternative = "two.sided", paired= FALSE, 
                     var.equal = FALSE, conf.level =0.95)
  ttestssex <- rbind(ttestssex, list(roi, ttestressex[[1]], ttestressex[3]))
}


colnames(ttestssex) <- c("ROIs", "tvalue", "pvalue")

comparison_sex_df <- merge(comparison_sex_df, ttestssex, all=TRUE, by ="ROIs")

#correct p-values for multiple comparisons
pvalues_sex <- as.list(comparison_sex_df[["pvalue"]])

comparison_sex_df["adj_pvalue_fdr"] <- p.adjust(pvalues_sex, method="fdr")


#write sex comparison
comparison_sex_df <- apply(comparison_sex_df,2,as.character)

write.csv(comparison_sex_df, "DS_sex_comparison.csv", row.names = FALSE)
