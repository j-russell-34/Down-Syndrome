library(tidyverse)
library(dplyr)
library(purrr)
library(glue)
library(readxl)
library(ggpp)
library(grid)
library(tibble)


#cbf_df <- read.csv("DSCHOL_CBF.csv")
samseg_df <- read_excel("DSCHOL_samseg.xlsx")
dnseg_df <- read.csv("DSCHOL_DnSeg_volumes.csv")
cog_df <- read.csv("DSCHOL_Cognition.csv")
sclimbic <- read_excel("FS7sclimbic_v0.xlsx")
feobv_suvr_df <- read.table("DSCHOL_FEOBV_SUVR.txt",
                            sep="\t", header=FALSE)


#cbf_df <- subset(cbf_df, select = c("SUBJECT", "SESSTYPE","CH123_L_VOL", 
 #                                   "CH123_R_VOL", "CH4_L_VOL", "CH4_R_VOL"))
samseg_df <- subset(samseg_df, select = c("SUBJECT", "samseg_sbtiv"))
sclimbic <- subset(sclimbic, select = c("SUBJECT", "SESSTYPE",
                                        "Left-Basal-Forebrain", 
                                        "Right-Basal-Forebrain"))

colnames(feobv_suvr_df) <- c("SUBJECT", "REF_REGION", "PVC", "ROI", "SUVR")

#change from long to wide
feobv_suvr_df_wide <- reshape(data = feobv_suvr_df,
                              idvar = "SUBJECT",
                              v.names = c("SUVR"),
                              timevar = "ROI",
                              direction = "wide")

#add column to feobv_df for timepoint of scan relative to trc_ds BL
feobv_suvr_df_wide["VISIT"] <- c("Baseline", "Baseline", "M16","Baseline",
                                 "Baseline", "Baseline", "Baseline", "Baseline",
                                 "Baseline", "Baseline")

#remove ROIs that are not of interest
feobv_df <- feobv_suvr_df_wide %>%
  select(-contains(c("White-Matter", "Cerebellum", "Brain-Stem", "CSF", 
                     "choroid", "Air", "Skull", "Vermis", "Pons", "Head")))

#convert - to . in column headers
names(feobv_df) <- make.names(names(feobv_df))



#full_df <- merge(cbf_df, samseg_df, by="SUBJECT")
full_df <- merge(samseg_df, dnseg_df, by="SUBJECT")

#colnames(full_df)[2] <- "VISIT"
colnames(sclimbic)[2] <- "VISIT"

#full_df[3,2] = "M16"
sclimbic[3,2] = "M16"

full_df <- merge(full_df, sclimbic, by=c("SUBJECT"))

#bf_rois <- list("CH123_L_VOL", "CH123_R_VOL", "CH4_L_VOL", "CH4_R_VOL", 
 #               "left", "right", "Left_Basal_Forebrain", 
  #              "Right_Basal_Forebrain")


#for (x in bf_rois) {
#  full_df[x] = full_df[x]/full_df["samseg_sbtiv"]
#}

names(cog_df)[names(cog_df) == "record_id"] <- 
  "SUBJECT"


new_cog_df<- aggregate(.~SUBJECT, data = cog_df, FUN = max, na.rm=TRUE, na.action=NULL)

new_cog_df <- new_cog_df %>%
  select(-contains(c("SUBJECT")))

names(new_cog_df)[names(new_cog_df) == "id"] <- 
  "SUBJECT"

full_df <- merge(full_df, new_cog_df, by=c("SUBJECT"))

full_df<- transform(full_df, iq_composite = as.numeric(iq_composite),
                             total_fr = as.numeric(total_fr),
                             total_score = as.numeric(total_score))

names(full_df)[names(full_df) == "Left.Basal.Forebrain"] <- 
  "Left_Basal_Forebrain"
names(full_df)[names(full_df) == "Right.Basal.Forebrain"] <- 
  "Right_Basal_Forebrain"

model <- lm(total_fr ~ Left_Basal_Forebrain + iq_composite + samseg_sbtiv, data=full_df)
summary(model)

model <- lm(total_fr ~ Right_Basal_Forebrain + iq_composite + samseg_sbtiv, data=full_df)
summary(model)

model <- lm(total_fr ~ left + iq_composite + samseg_sbtiv, data=full_df)
summary(model)

model <- lm(total_fr ~ right + iq_composite, data=full_df)
summary(model)

model <- lm(total_score ~ Left_Basal_Forebrain + iq_composite + samseg_sbtiv, data=full_df)
summary(model)

model <- lm(total_score ~ Right_Basal_Forebrain + iq_composite + samseg_sbtiv, data=full_df)
summary(model)

model <- lm(total_score ~ left + iq_composite + samseg_sbtiv, data=full_df)
summary(model)

model <- lm(total_score ~ right + iq_composite + samseg_sbtiv, data=full_df)
summary(model)

#effect of L BF volume on DSMSE total when controlling for KBIT IQ

#generate custom annotation to be placed at absolute 0.1,0.9
#grob <-grobTree(textGrob(expression(atop(beta == 2.186*10^5, "p = 0.027")),
 #                        x=0.1, y=0.9, hjust = 0,
  #                       gp=gpar(col="red", fontsize =14, fontface="bold")))

#png(file="L BF vs DSMSE total TRC-DS.png", units="in", width=7, height=5, res=300)
#ggplot(full_df,aes(x= Left_Basal_Forebrain, y=total_score)) +
#  geom_point(size=2, color="blue") + 
#  geom_smooth(method='lm', color="black")+
#  labs(y="DSMSE Total Score", x= "Left BF Volume, Normalised to ICV")+ 
#  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
 #                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
  #                   axis.text=element_text(size=14), axis.title=element_text(size=16, face="bold"))+
#  annotation_custom(grob)
#dev.off()



feobv_full_df <- merge(full_df, feobv_df, by=c("SUBJECT"))

#function to calculate p-value of linear model
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


#remove ROIs that are not of interest
feobv_full_df <- feobv_full_df %>%
  select(-contains(c("White-Matter", "Cerebellum", "Brain-Stem", "CSF", 
                     "choroid", "Air", "Skull", "Vermis", "Pons", "Head")))

#generate blank lists to add correlation analysis too
correlation_pvalues <- list()
beta <- list()
ROI <- list()

#iterate through list of ROIs, lm with cog performance value
for(i in 76:159){
  model <- lm (total_fr ~ feobv_full_df[,c(i)] + iq_composite, 
               data = feobv_full_df)
  correlation_pvalues <- append(correlation_pvalues, summary(model)$coefficients[,4][2])
  beta <-append(beta, summary(model)$coefficients[,1][2])
  ROI <- append(ROI, list(colnames(feobv_full_df)[i]))
}

df_list_totalfr <- list(ROI, beta, correlation_pvalues)

totalfr_correlation_df <- as.data.frame(do.call(cbind, df_list_totalfr))

#apply headers to df
colnames(totalfr_correlation_df)[1] <- "ROI"
colnames(totalfr_correlation_df)[2] <- "beta"
colnames(totalfr_correlation_df)[3] <- "pvalue"

#generate blank lists to add correlation analysis too
correlation_pvalues <- list()
beta <- list()
ROI <- list()

#iterate through list of ROIs, lm with cog performance value
for(i in 76:159){
  model <- lm (total_score ~ feobv_full_df[,c(i)] + iq_composite, 
               data = feobv_full_df)
  beta <-append(beta, summary(model)$coefficients[,1][2])
  correlation_pvalues <- append(correlation_pvalues, summary(model)$coefficients[,4][2])
  ROI <- append(ROI, list(colnames(feobv_full_df)[i]))
}

df_list_dsmse <- list(ROI, beta, correlation_pvalues)

dsmse_correlation_df <- as.data.frame(do.call(cbind, df_list_dsmse))

#apply headers to df
colnames(dsmse_correlation_df)[1] <- "ROI"
colnames(dsmse_correlation_df)[2] <- "beta"
colnames(dsmse_correlation_df)[3] <- "pvalue"

sig_dsmse_correlation <- dsmse_correlation_df[which(dsmse_correlation_df$pvalue <= 0.1),]
sig_fr_correlation <- totalfr_correlation_df[which(totalfr_correlation_df$pvalue <= 0.1),]

full_df<- transform(full_df, personal_info_score = as.numeric(personal_info_score),
                    memory_object = as.numeric(memory_object),
                    memory_location = as.numeric(memory_location),
                    apraxia_score = as.numeric(apraxia_score),
                    language_score = as.numeric(language_score),
                    visuospatial_score= as.numeric(visuospatial_score),
                    knowledgeofexaminerscore = as.numeric(knowledgeofexaminerscore))

dsmebeta <- list()
dsmepval <- list()
subgroups <-list("Personal Information Score","Personal Information Score",
                 "Memory object", "Memory object", "Memory Location",
                 "Memory Location", "Apraxia Score", "Apraxia Score", 
                 "Language Score", "Language Score", "Visuospatial Score",
                 "Visuospatial Score", "Knowledge of the Examiner",
                 "Knowledge of the Examiner")

model <- lm(personal_info_score ~ Left_Basal_Forebrain + iq_composite + samseg_sbtiv, data=full_df)
summary(model)

dsmebeta <-append(dsmebeta, summary(model)$coefficients[,1][2])
dsmepval <- append(dsmepval, summary(model)$coefficients[,4][2])

model <- lm(personal_info_score ~ Right_Basal_Forebrain + iq_composite + samseg_sbtiv, data=full_df)
summary(model)

dsmebeta <-append(dsmebeta, summary(model)$coefficients[,1][2])
dsmepval <- append(dsmepval, summary(model)$coefficients[,4][2])

model <- lm(memory_object ~ Left_Basal_Forebrain + iq_composite + samseg_sbtiv, data=full_df)
summary(model)

dsmebeta <-append(dsmebeta, summary(model)$coefficients[,1][2])
dsmepval <- append(dsmepval, summary(model)$coefficients[,4][2])

model <- lm(memory_object ~ Right_Basal_Forebrain + iq_composite + samseg_sbtiv, data=full_df)
summary(model)

dsmebeta <-append(dsmebeta, summary(model)$coefficients[,1][2])
dsmepval <- append(dsmepval, summary(model)$coefficients[,4][2])

model <- lm(memory_location ~ Left_Basal_Forebrain + iq_composite+ samseg_sbtiv, data=full_df)
summary(model)

dsmebeta <-append(dsmebeta, summary(model)$coefficients[,1][2])
dsmepval <- append(dsmepval, summary(model)$coefficients[,4][2])

model <- lm(memory_location ~ Right_Basal_Forebrain + iq_composite + samseg_sbtiv, data=full_df)
summary(model)

dsmebeta <-append(dsmebeta, summary(model)$coefficients[,1][2])
dsmepval <- append(dsmepval, summary(model)$coefficients[,4][2])

model <- lm(apraxia_score ~ Left_Basal_Forebrain + iq_composite+ samseg_sbtiv, data=full_df)
summary(model)

dsmebeta <-append(dsmebeta, summary(model)$coefficients[,1][2])
dsmepval <- append(dsmepval, summary(model)$coefficients[,4][2])

model <- lm(apraxia_score ~ Right_Basal_Forebrain + iq_composite+ samseg_sbtiv, data=full_df)
summary(model)

dsmebeta <-append(dsmebeta, summary(model)$coefficients[,1][2])
dsmepval <- append(dsmepval, summary(model)$coefficients[,4][2])

model <- lm(language_score ~ Left_Basal_Forebrain + iq_composite+ samseg_sbtiv, data=full_df)
summary(model)

dsmebeta <-append(dsmebeta, summary(model)$coefficients[,1][2])
dsmepval <- append(dsmepval, summary(model)$coefficients[,4][2])

model <- lm(language_score ~ Right_Basal_Forebrain + iq_composite+ samseg_sbtiv, data=full_df)
summary(model)

dsmebeta <-append(dsmebeta, summary(model)$coefficients[,1][2])
dsmepval <- append(dsmepval, summary(model)$coefficients[,4][2])

model <- lm(visuospatial_score ~ Left_Basal_Forebrain + iq_composite+ samseg_sbtiv, data=full_df)
summary(model)

dsmebeta <-append(dsmebeta, summary(model)$coefficients[,1][2])
dsmepval <- append(dsmepval, summary(model)$coefficients[,4][2])

model <- lm(visuospatial_score ~ Right_Basal_Forebrain + iq_composite+ samseg_sbtiv, data=full_df)
summary(model)

dsmebeta <-append(dsmebeta, summary(model)$coefficients[,1][2])
dsmepval <- append(dsmepval, summary(model)$coefficients[,4][2])

model <- lm(knowledgeofexaminerscore ~ Left_Basal_Forebrain + iq_composite+ samseg_sbtiv, data=full_df)
summary(model)

dsmebeta <-append(dsmebeta, summary(model)$coefficients[,1][2])
dsmepval <- append(dsmepval, summary(model)$coefficients[,4][2])

model <- lm(knowledgeofexaminerscore ~ Right_Basal_Forebrain + iq_composite+ samseg_sbtiv, data=full_df)
summary(model)

dsmebeta <-append(dsmebeta, summary(model)$coefficients[,1][2])
dsmepval <- append(dsmepval, summary(model)$coefficients[,4][2])

df_list_dsmse_subs <- list(subgroups, dsmebeta, dsmepval)

dsmsesubs_correlation_df <- as.data.frame(do.call(cbind, df_list_dsmse_subs))

#apply headers to df
colnames(dsmsesubs_correlation_df)[1] <- "Sub-test"
colnames(dsmsesubs_correlation_df)[2] <- "beta"
colnames(dsmsesubs_correlation_df)[3] <- "pvalue"

#generate blank lists to add correlation analysis too
personalinfo_pvalues <- list()
personalinfobeta <- list()
personalinfoROI <- list()

feobv_full_df<- transform(feobv_full_df, personal_info_score = as.numeric(personal_info_score),
                    memory_object = as.numeric(memory_object),
                    memory_location = as.numeric(memory_location),
                    apraxia_score = as.numeric(apraxia_score),
                    language_score = as.numeric(language_score),
                    visuospatial_score= as.numeric(visuospatial_score),
                    knowledgeofexaminerscore = as.numeric(knowledgeofexaminerscore))

#iterate through list of ROIs, lm with cog performance value
for(i in 76:159){
  model <- lm (personal_info_score ~ feobv_full_df[,c(i)] + iq_composite, 
               data = feobv_full_df)
  personalinfo_pvalues <- append(personalinfo_pvalues, summary(model)$coefficients[,4][2])
  personalinfobeta <-append(personalinfobeta, summary(model)$coefficients[,1][2])
  personalinfoROI <- append(personalinfoROI, list(colnames(feobv_full_df)[i]))
}

df_list_personalinfo <- list(personalinfoROI, personalinfobeta, personalinfo_pvalues)

personalinfo_correlation_df <- as.data.frame(do.call(cbind, df_list_personalinfo))

#apply headers to df
colnames(personalinfo_correlation_df)[1] <- "ROI"
colnames(personalinfo_correlation_df)[2] <- "beta"
colnames(personalinfo_correlation_df)[3] <- "pvalue"

sig_personalinfo_correlation <- personalinfo_correlation_df[which(personalinfo_correlation_df$pvalue <= 0.1),]


memobj_pvalues <- list()
memobjbeta <- list()
memobjROI <- list()

for(i in 76:159){
  model <- lm (memory_object ~ feobv_full_df[,c(i)] + iq_composite, 
               data = feobv_full_df)
  memobj_pvalues <- append(memobj_pvalues, summary(model)$coefficients[,4][2])
  memobjbeta <-append(memobjbeta, summary(model)$coefficients[,1][2])
  memobjROI <- append(memobjROI, list(colnames(feobv_full_df)[i]))
}

df_list_memobj <- list(memobjROI, memobjbeta, memobj_pvalues)

memobj_correlation_df <- as.data.frame(do.call(cbind, df_list_memobj))

#apply headers to df
colnames(memobj_correlation_df)[1] <- "ROI"
colnames(memobj_correlation_df)[2] <- "beta"
colnames(memobj_correlation_df)[3] <- "pvalue"

sig_memoryobject_correlation <- memobj_correlation_df[which(memobj_correlation_df$pvalue <= 0.1),]



memloc_pvalues <- list()
memlocbeta <- list()
memlocROI <- list()


for(i in 76:159){
  model <- lm (memory_location ~ feobv_full_df[,c(i)] + iq_composite, 
               data = feobv_full_df)
  memloc_pvalues <- append(memloc_pvalues, summary(model)$coefficients[,4][2])
  memlocbeta <-append(memlocbeta, summary(model)$coefficients[,1][2])
  memlocROI <- append(memlocROI, list(colnames(feobv_full_df)[i]))
}

df_list_memloc <- list(memlocROI, memlocbeta, memloc_pvalues)

memloc_correlation_df <- as.data.frame(do.call(cbind, df_list_memloc))

#apply headers to df
colnames(memloc_correlation_df)[1] <- "ROI"
colnames(memloc_correlation_df)[2] <- "beta"
colnames(memloc_correlation_df)[3] <- "pvalue"

sig_memorylocation_correlation <- memloc_correlation_df[which(memloc_correlation_df$pvalue <= 0.1),]


apraxia_pvalues <- list()
apraxiabeta <- list()
apraxiaROI <- list()


for(i in 76:159){
  model <- lm (apraxia_score ~ feobv_full_df[,c(i)] + iq_composite, 
               data = feobv_full_df)
  apraxia_pvalues <- append(apraxia_pvalues, summary(model)$coefficients[,4][2])
  apraxiabeta <-append(apraxiabeta, summary(model)$coefficients[,1][2])
  apraxiaROI <- append(apraxiaROI, list(colnames(feobv_full_df)[i]))
}

df_list_apraxia <- list(apraxiaROI, apraxiabeta, apraxia_pvalues)

apraxia_correlation_df <- as.data.frame(do.call(cbind, df_list_apraxia))

#apply headers to df
colnames(apraxia_correlation_df)[1] <- "ROI"
colnames(apraxia_correlation_df)[2] <- "beta"
colnames(apraxia_correlation_df)[3] <- "pvalue"

sig_apraxia_correlation <- apraxia_correlation_df[which(apraxia_correlation_df$pvalue <= 0.1),]



language_pvalues <- list()
languagebeta <- list()
languageROI <- list()


for(i in 76:159){
  model <- lm (language_score ~ feobv_full_df[,c(i)] + iq_composite, 
               data = feobv_full_df)
  language_pvalues <- append(language_pvalues, summary(model)$coefficients[,4][2])
  languagebeta <-append(languagebeta, summary(model)$coefficients[,1][2])
  languageROI <- append(languageROI, list(colnames(feobv_full_df)[i]))
}

df_list_language <- list(languageROI, languagebeta, language_pvalues)

language_correlation_df <- as.data.frame(do.call(cbind, df_list_language))

#apply headers to df
colnames(language_correlation_df)[1] <- "ROI"
colnames(language_correlation_df)[2] <- "beta"
colnames(language_correlation_df)[3] <- "pvalue"

sig_language_correlation <- language_correlation_df[which(language_correlation_df$pvalue <= 0.1),]



vs_pvalues <- list()
vsbeta <- list()
vsROI <- list()


for(i in 76:159){
  model <- lm (visuospatial_score ~ feobv_full_df[,c(i)] + iq_composite, 
               data = feobv_full_df)
  vs_pvalues <- append(vs_pvalues, summary(model)$coefficients[,4][2])
  vsbeta <-append(vsbeta, summary(model)$coefficients[,1][2])
  vsROI <- append(vsROI, list(colnames(feobv_full_df)[i]))
}

df_list_vs <- list(vsROI, vsbeta, vs_pvalues)

vs_correlation_df <- as.data.frame(do.call(cbind, df_list_vs))

#apply headers to df
colnames(vs_correlation_df)[1] <- "ROI"
colnames(vs_correlation_df)[2] <- "beta"
colnames(vs_correlation_df)[3] <- "pvalue"

sig_visuospatial_correlation <- vs_correlation_df[which(vs_correlation_df$pvalue <= 0.1),]


examiner_pvalues <- list()
examinerbeta <- list()
examinerROI <- list()


for(i in 76:159){
  model <- lm (knowledgeofexaminerscore ~ feobv_full_df[,c(i)] + iq_composite, 
               data = feobv_full_df)
  examiner_pvalues <- append(examiner_pvalues, summary(model)$coefficients[,4][2])
  examinerbeta <-append(examinerbeta, summary(model)$coefficients[,1][2])
  examinerROI <- append(examinerROI, list(colnames(feobv_full_df)[i]))
}

df_list_examiner <- list(examinerROI, examinerbeta, examiner_pvalues)

examiner_correlation_df <- as.data.frame(do.call(cbind, df_list_examiner))

#apply headers to df
colnames(examiner_correlation_df)[1] <- "ROI"
colnames(examiner_correlation_df)[2] <- "beta"
colnames(examiner_correlation_df)[3] <- "pvalue"

sig_examinerknowledge_correlation <- examiner_correlation_df[which(examiner_correlation_df$pvalue <= 0.1),]

iq_pvalues <- list()
iqbeta <- list()
iqROI <- list()


for(i in 76:159){
  model <- lm (iq_composite ~ feobv_full_df[,c(i)], 
               data = feobv_full_df)
  iq_pvalues <- append(iq_pvalues, summary(model)$coefficients[,4][2])
  iqbeta <-append(iqbeta, summary(model)$coefficients[,1][2])
  iqROI <- append(iqROI, list(colnames(feobv_full_df)[i]))
}

df_list_iq <- list(iqROI, iqbeta, iq_pvalues)

iq_correlation_df <- as.data.frame(do.call(cbind, df_list_iq))

#apply headers to df
colnames(iq_correlation_df)[1] <- "ROI"
colnames(iq_correlation_df)[2] <- "beta"
colnames(iq_correlation_df)[3] <- "pvalue"

sig_iq_correlation <- iq_correlation_df[which(iq_correlation_df$pvalue <= 0.1),]





