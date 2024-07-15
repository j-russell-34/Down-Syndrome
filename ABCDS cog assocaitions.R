library(tidyverse)
library(dplyr)
library(purrr)
library(glue)
library(readxl)
library(ggpp)
library(grid)
library(tibble)
library(rstatix)

#import data
cuedrecall_df <- read.csv("Cued_Recall.csv")
dnseg_df <- read.csv("ABC-DS DnSeg DF for R 20231221.csv")
sclimbic_df <- read.csv("ABC-DS ScLimbic DF for R 20231221.csv")
catsanddogs_df <- read.csv("Cats_and_Dogs_task.csv")
dems_df <- read.csv("Participant_Demographics.csv")
dmsme_df <- read.csv("Down_Syndrome_Mental_Status_Exam.csv")
samseg_df <-read_csv("ABCDS_SAMSEG.csv")
sclimbicdata <- read_csv("ABCDS_sclimbic_volumes.csv")
dnsegdata <- read_csv("ABCDS_DnSeg_volumes.csv")

#select BL data
cuedrecall_df <- cuedrecall_df[which(cuedrecall_df$event_code == "bl"),]
catsanddogs_df <- catsanddogs_df[which(catsanddogs_df$event_code == "bl"),]
dems_df <- dems_df[which(dems_df$event_code == "bl"),]
dmsme_df <- dmsme_df[which(dmsme_df$event_code =="bl"),]

#select columns of interest, change intellectual disability
dems_df <- subset(dems_df, select = c("subject_label","pmiimp"))
colnames(dems_df) <- c("SUBJECT", "Intellectual_disability")

dems_df["Intellectual_disability"][dems_df["Intellectual_disability"] == 
                                     1] <- "Mild"
dems_df["Intellectual_disability"][dems_df["Intellectual_disability"] == 
                                     2] <- "Moderate"
dems_df["Intellectual_disability"][dems_df["Intellectual_disability"] == 
                                     3] <- "Severe"

cuedrecall_df["Total_Free_Recall"] <- cuedrecall_df["frt1"] + 
  cuedrecall_df["frt2"] +cuedrecall_df["frt3"]

cuedrecall_df <- subset(cuedrecall_df, select = c("subject_label", 
                                                  "Total_Free_Recall"))
colnames(cuedrecall_df)[which(names(cuedrecall_df) == "subject_label")] <-
  "SUBJECT"


dmsme_df <- subset(dmsme_df, select = c("subject_label", "dsme1", "dsto1"))
colnames(dmsme_df) <- c("SUBJECT", "DSMSE_Location_Memory", "DSMSE_Total")

catsanddogs_df["cats_dogs_switching"] <-catsanddogs_df["ttscd"] -
  catsanddogs_df["ttncd"]

catsanddogs_df <- subset(catsanddogs_df, select = c("subject_label",
                                                "cats_dogs_switching"))
colnames(catsanddogs_df)[which(names(catsanddogs_df) == "subject_label")] <-
  "SUBJECT"

sclimbic_df <- subset(sclimbic_df, select = -c(Gender, Tracer))
dnseg_df <- subset(dnseg_df, select = -c(Gender, Tracer))


dflist <- list(catsanddogs_df, cuedrecall_df,dems_df,dmsme_df)

sclimbic_full_df <- merge(sclimbic_df, sclimbicdata, by="SUBJECT")
dnseg_full_df <-merge(dnseg_df, dnsegdata, by="SUBJECT")

sclimbic_full_df <- merge(sclimbic_full_df, samseg_df, by="SUBJECT")
dnseg_full_df <-merge(dnseg_full_df, samseg_df, by="SUBJECT")

for (l in dflist){
  sclimbic_full_df <-merge(sclimbic_full_df,l, by="SUBJECT")
  dnseg_full_df <-merge(dnseg_full_df,l, by="SUBJECT")
}

names(sclimbic_full_df) <- gsub("\\-","_", names(sclimbic_full_df))

model <- lm(cats_dogs_switching ~ Left_Basal_Forebrain +
              Intellectual_disability + samseg_sbtiv, data = sclimbic_full_df)
summary(model)

model <- lm(cats_dogs_switching ~ Right_Basal_Forebrain +
              Intellectual_disability+ samseg_sbtiv, data = sclimbic_full_df)
summary(model)

model <- lm(Total_Free_Recall ~ Left_Basal_Forebrain +
              Intellectual_disability+ samseg_sbtiv, data = sclimbic_full_df)
summary(model)


#generate custom annotation to be placed at absolute 0.025,0.08
grob <-grobTree(textGrob(expression(atop(beta == 0.0386, "p = 0.00267")),
                         x=0.1, y=0.1, hjust = 0,
                         gp=gpar(col="red", fontsize =14, fontface="bold")))

png(file="Free recall vs left BF vol.png", units="in", width=7, height=5, res=300)
ggplot(sclimbic_full_df,aes(x=Left_Basal_Forebrain, 
                            y=Total_Free_Recall, col=Intellectual_disability)) +
  geom_point(size=2) + 
  geom_smooth(method='lm')+
  labs(y="Total Free Recall (Cued Recall Task)", x= expression("Left ChBF, mm"^3), 
       col="Premorbid\nIntellectual\nDisability")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=14), axis.title=element_text(size=16))+
  annotation_custom(grob)
dev.off()

model <- lm(Total_Free_Recall ~ Right_Basal_Forebrain +
              Intellectual_disability +samseg_sbtiv, data = sclimbic_full_df)
summary(model)

model <- lm(DSMSE_Total ~ Left_Basal_Forebrain +
              Intellectual_disability + samseg_sbtiv, data = sclimbic_full_df)
summary(model)

#generate custom annotation to be placed at absolute 0.1,0.1

grob <-grobTree(textGrob(expression(atop(beta == 0.0343, "p = 0.0122")),
                         x=0.1, y=0.1, hjust = 0,
                         gp=gpar(col="red", fontsize =14, fontface="bold")))

png(file="DSMSE Total vs left BF vol.png", units="in", width=7, height=5, res=300)
ggplot(sclimbic_full_df,aes(x=Left_Basal_Forebrain, 
                            y=DSMSE_Total, col=Intellectual_disability)) +
  geom_point(size=2) + 
  geom_smooth(method='lm')+
  labs(y="DSMSE Total Score", x= expression("Left ChBF, mm"^3), 
       col="Premorbid\nIntellectual\nDisability")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=14), axis.title=element_text(size=16))+
  annotation_custom(grob)
dev.off()

model <- lm(DSMSE_Total ~ Right_Basal_Forebrain +
              Intellectual_disability+ samseg_sbtiv, data = sclimbic_full_df)
summary(model)

#generate custom annotation to be placed at absolute 0.025,0.08
grob <-grobTree(textGrob(expression(atop(beta == 0.0501, "p = 0.0007")),
                         x=0.1, y=0.1, hjust = 0,
                         gp=gpar(col="red", fontsize =14, fontface="bold")))


png(file="DSMSE Total vs right BF vol.png", units="in", width=7, height=5, res=300)
ggplot(sclimbic_full_df,aes(x=Right_Basal_Forebrain, 
                            y=DSMSE_Total, col=Intellectual_disability)) +
  geom_point(size=2) + 
  geom_smooth(method='lm')+
  labs(y="DSMSE Total Score", x= expression("Right ChBF, mm"^3), 
       col="Premorbid\nIntellectual\nDisability")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=14), axis.title=element_text(size=16))+
  annotation_custom(grob)
dev.off()

model <- lm(cats_dogs_switching ~ Left_DnSeg +
              Intellectual_disability + samseg_sbtiv, data = dnseg_full_df)
summary(model)

model <- lm(cats_dogs_switching ~ Right_DnSeg +
              Intellectual_disability + samseg_sbtiv, data = dnseg_full_df)
summary(model)

model <- lm(Total_Free_Recall ~ Left_DnSeg +
              Intellectual_disability+ samseg_sbtiv, data = dnseg_full_df)
summary(model)

model <- lm(Total_Free_Recall ~ Right_DnSeg +
              Intellectual_disability + samseg_sbtiv, data = dnseg_full_df)
summary(model)

model <- lm(DSMSE_Total ~ Left_DnSeg +
              Intellectual_disability+ samseg_sbtiv, data = dnseg_full_df)
summary(model)

model <- lm(DSMSE_Total ~ Right_DnSeg +
              Intellectual_disability + samseg_sbtiv, data = dnseg_full_df)
summary(model)

res.aov <- sclimbic_full_df %>% anova_test(Right_Basal_Forebrain ~ 
                                             Age + Intellectual_disability + samseg_sbtiv)
get_anova_table(res.aov)

res.aov <- sclimbic_full_df %>% anova_test(Left_Basal_Forebrain ~ 
                                             Age + Intellectual_disability + samseg_sbtiv)
get_anova_table(res.aov)

res.aov <- dnseg_full_df %>% anova_test(Right_DnSeg ~ 
                                          Age + Intellectual_disability + samseg_sbtiv)
get_anova_table(res.aov)


res.aov <- dnseg_full_df %>% anova_test(Left_DnSeg ~ 
                                          Age + Intellectual_disability + samseg_sbtiv)
get_anova_table(res.aov)
