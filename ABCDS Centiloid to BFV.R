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
dems_df <- read.csv("DSChol_dems.csv")
centiloids <- read_excel("TRC-DS centiloids 20231215.xlsx")
sclimbic <- read_excel("FS7sclimbic_v0.xlsx")
cog_df <- read.csv("DSCHOL_Cognition.csv")

colnames(dems_df)[1] <- "SUBJECT" 

#combine dfs based on subject

#cbf_df <- subset(cbf_df, select = c("SUBJECT", "SESSTYPE","CH123_L_VOL", 
 #                                   "CH123_R_VOL", "CH4_L_VOL", "CH4_R_VOL"))
samseg_df <- subset(samseg_df, select = c("SUBJECT", "samseg_sbtiv"))
centiloids <- subset(centiloids, select = c("SUBJECT", "VISIT", "Centiloid"))
sclimbic <- subset(sclimbic, select = c("SUBJECT", "SESSTYPE",
                                        "Left-Basal-Forebrain", 
                                        "Right-Basal-Forebrain"))

names(sclimbic) <- gsub("\\-","_", names(sclimbic))

#full_df <- merge(cbf_df, samseg_df, by="SUBJECT")
full_df <- merge(samseg_df, dnseg_df, by="SUBJECT")
full_df <- merge(full_df, dems_df, by="SUBJECT")


colnames(sclimbic)[2] <- "VISIT"


sclimbic[3,2] = "M16"

full_df <- merge(full_df, centiloids, by=c("SUBJECT"))
full_df <- merge(full_df, sclimbic, by=c("SUBJECT", "VISIT"))

#linear model 

bf_rois <- list("CH123_L_VOL", "CH123_R_VOL", "CH4_L_VOL", "CH4_R_VOL", 
                "left", "right", "Left_Basal_Forebrain", 
                "Right_Basal_Forebrain")



names(full_df)[names(full_df) == "left"] <- "Dnseg_left"
names(full_df)[names(full_df) == "right"] <- "Dnseg_right"

#model <- lm(CH123_L_VOL ~ Centiloid, data=full_df)
#summary(model)



#ggplot(full_df, aes(x=Centiloid, y=CH123_L_VOL)) + geom_point()
#ggplot(full_df, aes(x=Centiloid, y=CH123_R_VOL)) + geom_point()
#ggplot(full_df, aes(x=Centiloid, y=CH4_L_VOL)) + geom_point()
#ggplot(full_df, aes(x=Centiloid, y=CH4_R_VOL)) + geom_point()

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

#Dnseg_left
model <- lm (Dnseg_left ~ Centiloid + samseg_sbtiv, data = full_df)
summary(model)

#R2 <- round(R2[[1]],digits= 4)
#pval <- round(pval[[1]],digits = 4)
#R2 <- sprintf("%.4f", R2)
#pval <- sprintf("%.4f", pval)
grob <-grobTree(textGrob(expression(atop(beta == 0.095, "p = 0.927")), x=0.85, y=0.1, hjust = 0,
                         gp=gpar(col="red", fontsize =14)))

png(file="Amyloid vs L nbM trc-ds.png", units="in", width=7, height=5, res=300)
ggplot(full_df,aes(x=Centiloid, y=Dnseg_left)) +
  geom_point(size=2, color="blue") + 
  geom_smooth(method='lm', color="black")+
  labs(y=expression("Left nbM volume, mm"^3), x= "Global Amyloid Accumulation (Centiloids)")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=14), axis.title=element_text(size=16))+
  annotation_custom(grob)
dev.off()

#DnSeg Right
model <- lm (Dnseg_right ~ Centiloid +samseg_sbtiv, data = full_df)
summary(model)


grob <-grobTree(textGrob(expression(atop(beta == 0.263, "p = 0.803")), x=0.85, y=0.1, hjust = 0,
                         gp=gpar(col="red", fontsize =14)))

png(file="Amyloid vs R nbM trc-ds.png", units="in", width=7, height=5, res=300)
ggplot(full_df,aes(x=Centiloid, y=Dnseg_right)) +
  geom_point(size=2, color="blue") + 
  geom_smooth(method='lm', color="black")+
  labs(y=expression("Right nbM volume, mm"^3), x= "Global Amyloid Accumulation (Centiloids)")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=14), axis.title=element_text(size=16))+
  annotation_custom(grob)
dev.off()

#ChBF left
model <- lm (Left_Basal_Forebrain ~ Centiloid +samseg_sbtiv, data = full_df)
summary(model)

grob <-grobTree(textGrob(expression(atop(beta == 1.056, "p = 0.223")), x=0.85, y=0.1, hjust = 0,
                         gp=gpar(col="red", fontsize =14)))


png(file="Amyloid vs L BF trc-ds.png", units="in", width=7, height=5, res=300)
ggplot(full_df,aes(x=Centiloid, y=Left_Basal_Forebrain)) +
  geom_point(size=2, color="blue") + 
  geom_smooth(method='lm', color="black")+
  labs(y=expression("Left BF volume, mm"^3), x= "Global Amyloid Accumulation (Centiloids)")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=14), axis.title=element_text(size=16))+
  annotation_custom(grob)
dev.off()

#ChBF Right
model <- lm (Right_Basal_Forebrain ~ Centiloid +samseg_sbtiv, data = full_df)
summary(model)

grob <-grobTree(textGrob(expression(atop(beta == 1.11, "p = 0.0707")), x=0.85, y=0.1, hjust = 0,
                         gp=gpar(col="red", fontsize =14)))


png(file="Amyloid vs R BF trc-ds.png", units="in", width=7, height=5, res=300)
ggplot(full_df,aes(x=Centiloid, y=Right_Basal_Forebrain)) +
  geom_point(size=2, color="blue") + 
  geom_smooth(method='lm', color="black")+
  labs(y=expression("Right BF volume, mm"^3), x= "Global Amyloid Accumulation (Centiloids)")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=14), axis.title=element_text(size=16))+
  annotation_custom(grob)
dev.off()

#amyloid age

model <- lm (Centiloid ~ dems_age, data = full_df)
R2 <- summary(model)$r.squared
pval <- lmp(model)

R2 <- round(R2[[1]],digits= 4)
pval <- round(pval[[1]],digits = 4)
R2 <- sprintf("%.4f", R2)
pval <- sprintf("%.4f", pval)
grob <-grobTree(textGrob(glue('R² = {R2}\np = {pval}'), x=0.1, y=0.1, hjust = 0,
                         gp=gpar(col="red", fontsize =14)))

png(file="Amyloid vs Age trc-ds.png", units="in", width=7, height=5, res=300)
ggplot(full_df,aes(x=dems_age, y=Centiloid)) +
  geom_point(size=2, color="blue") + 
  geom_smooth(method='lm', color="black")+
  labs(y="Global Amyloid Accumulation (Centiloids)", x= "Age")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=14), axis.title=element_text(size=16))+
  annotation_custom(grob)
dev.off()

# ChBF vol Age
#Dnseg_left
model <- lm (Dnseg_left ~ dems_age +samseg_sbtiv, data = full_df)
summary(model)

grob <-grobTree(textGrob(expression(atop(beta == 0.159, "p = 0.949")), x=0.85, y=0.1, hjust = 0,
                         gp=gpar(col="red", fontsize =14)))

png(file="Age vs L nbM trc-ds.png", units="in", width=7, height=5, res=300)
ggplot(full_df,aes(x=dems_age, y=Dnseg_left)) +
  geom_point(size=2, color="blue") + 
  geom_smooth(method='lm', color="black")+
  labs(y=expression("Left nbM volume, mm"^3), x= "Age")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=14), axis.title=element_text(size=16))+
  annotation_custom(grob)
dev.off()

#DnSeg Right
model <- lm (Dnseg_right ~ dems_age +samseg_sbtiv, data = full_df)
summary(model)

grob <-grobTree(textGrob(expression(atop(beta == -0.416, "p = 0.871")), x=0.85, y=0.1, hjust = 0,
                         gp=gpar(col="red", fontsize =14)))

png(file="Age vs R nbM trc-ds.png", units="in", width=7, height=5, res=300)
ggplot(full_df,aes(x=Centiloid, y=Dnseg_right)) +
  geom_point(size=2, color="blue") + 
  geom_smooth(method='lm', color="black")+
  labs(y=expression("Right nbM volume, mm"^3), x= "Age")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=14), axis.title=element_text(size=16))+
  annotation_custom(grob)
dev.off()

#ChBF left
model <- lm (Left_Basal_Forebrain ~ dems_age +samseg_sbtiv, data = full_df)
summary(model)

grob <-grobTree(textGrob(expression(atop(beta == 0.604, "p = 0.788")), x=0.85, y=0.1, hjust = 0,
                         gp=gpar(col="red", fontsize =14)))

png(file="Age vs L BF trc-ds.png", units="in", width=7, height=5, res=300)
ggplot(full_df,aes(x=dems_age, y=Left_Basal_Forebrain)) +
  geom_point(size=2, color="blue") + 
  geom_smooth(method='lm', color="black")+
  labs(y=expression("Left BF volume, mm"^3), x= "Age")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=14), axis.title=element_text(size=16))+
  annotation_custom(grob)
dev.off()

#ChBF Right
model <- lm (Right_Basal_Forebrain ~ dems_age +samseg_sbtiv, data = full_df)
summary(model)

grob <-grobTree(textGrob(expression(atop(beta == 0.320, "p = 0.853")), x=0.85, y=0.1, hjust = 0,
                         gp=gpar(col="red", fontsize =14)))
png(file="Age vs R BF trc-ds.png", units="in", width=7, height=5, res=300)
ggplot(full_df,aes(x=dems_age, y=Right_Basal_Forebrain)) +
  geom_point(size=2, color="blue") + 
  geom_smooth(method='lm', color="black")+
  labs(y=expression("Right BF volume, mm"^3), x= "Age")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=14), axis.title=element_text(size=16))+
  annotation_custom(grob)
dev.off()

names(cog_df)[names(cog_df) == "record_id"] <- 
  "SUBJECT"


new_cog_df<- aggregate(.~SUBJECT, data = cog_df, FUN = max, na.rm=TRUE, na.action=NULL)

new_cog_df <- new_cog_df %>%
  select(-contains(c("SUBJECT")))

names(new_cog_df)[names(new_cog_df) == "id"] <- 
  "SUBJECT"

full_df <- merge(full_df, new_cog_df, by=c("SUBJECT"))

full_df<- transform(full_df, iq_composite = as.numeric(iq_composite))

model <- lm (Dnseg_left ~ Centiloid +iq_composite+samseg_sbtiv, data = full_df)
summary(model)

model <- lm (Dnseg_right ~ Centiloid +iq_composite+samseg_sbtiv, data = full_df)
summary(model)

model <- lm (Right_Basal_Forebrain ~ Centiloid +iq_composite+samseg_sbtiv, data = full_df)
summary(model)

model <- lm (Left_Basal_Forebrain ~ Centiloid +iq_composite+samseg_sbtiv, data = full_df)
summary(model)

model <- lm (Dnseg_left ~ dems_age +iq_composite+samseg_sbtiv, data = full_df)
summary(model)

model <- lm (Dnseg_right ~ dems_age +iq_composite+samseg_sbtiv+samseg_sbtiv, data = full_df)
summary(model)

model <- lm (Right_Basal_Forebrain ~ dems_age +iq_composite+samseg_sbtiv, data = full_df)
summary(model)

model <- lm (Left_Basal_Forebrain ~ dems_age +iq_composite+samseg_sbtiv, data = full_df)
summary(model)

