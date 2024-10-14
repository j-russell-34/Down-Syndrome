library(tidyverse)
library(dplyr)
library(purrr)
library(glue)
library(readxl)
library(ggpp)
library(grid)

in_dir <- '/Users/jasonrussell/Documents/Study_data/Down Syndrome/ABCDS'

dnsegdata_df <- read.csv(glue("{in_dir}/ABCDS_DnSeg_volumes.csv"))
sclimbicdata_df <- read.csv(glue("{in_dir}/ABCDS_sclimbic_volumes.csv"))
samseg_df <- read.csv(glue("{in_dir}/ABCDS_SAMSEG.csv"))
dnseg_abcds_df <- read.csv(glue("{in_dir}/ABC-DS DnSeg DF for R 20231221.csv"))
sclimbic_abcds_df <- read.csv(glue("{in_dir}/ABC-DS ScLimbic DF for R 20231221.csv"))
ad_diagnosis_df <-read.csv(glue("{in_dir}/Diagnostic_Consensus_Classification.csv"))

#clean ad_diagnosis_df to just include baseline diagnosis
ad_diagnosis_df <- ad_diagnosis_df[ad_diagnosis_df$event_code =="bl",]
colnames(ad_diagnosis_df)[which(names(ad_diagnosis_df) == "subject_label")] <- 
  "SUBJECT"


#function to calculate p-value of linear model
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

#join dataframes by subject ID

sclimbic_abcds_df <- merge(sclimbic_abcds_df, sclimbicdata_df, by="SUBJECT")
sclimbic_abcds_df <- merge(sclimbic_abcds_df, samseg_df, by="SUBJECT")
sclimbic_abcds_df <- merge(sclimbic_abcds_df, ad_diagnosis_df, by="SUBJECT")
sclimbic_abcds_df <- sclimbic_abcds_df[!duplicated(sclimbic_abcds_df$SUBJECT),]

dnseg_abcds_df <- merge(dnseg_abcds_df, dnsegdata_df, by="SUBJECT")
dnseg_abcds_df <- merge(dnseg_abcds_df, samseg_df, by="SUBJECT")
dnseg_abcds_df <- merge(dnseg_abcds_df, ad_diagnosis_df, by="SUBJECT")
dnseg_abcds_df <- dnseg_abcds_df[!duplicated(dnseg_abcds_df$SUBJECT),]

colnames(sclimbic_abcds_df)[which(names(sclimbic_abcds_df) == "Gender")] <- "Sex"
colnames(dnseg_abcds_df)[which(names(dnseg_abcds_df) == "Gender")] <- "Sex"

#calculate R2 and p then graph for each volume vs centiloid. Then split into 
#male/female

#Dnseg Left
model <- lm (Left_DnSeg ~ Amyloid..centiloids.+ samseg_sbtiv,
             data = dnseg_abcds_df)
summary(model)



#generate custom annotation to be placed at absolute 0.8,0.92
grob <-grobTree(textGrob(expression(atop(beta == -0.141, "p = 0.0415")), x=0.8, y=0.92, hjust = 0,
                         gp=gpar(col="red", fontsize =18, fontface="bold")))

png(file="Left nbM vs Amyloid.png", units="in", width=7, height=5, res=300)
ggplot(dnseg_abcds_df,aes(Amyloid..centiloids., Left_DnSeg)) +
  geom_point(size=2, color="blue") + 
  geom_smooth(method='lm', color="black")+
  labs(y=expression("Left nbM Volume, mm"^3), x= "Global Amyloid Accumulation (Centiloids)")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=18), axis.title=element_text(size=20))+
  annotation_custom(grob)
dev.off()

model <- lm (Left_DnSeg ~ Amyloid..centiloids. * Sex + samseg_sbtiv ,
             data = dnseg_abcds_df)
summary(model)


#Dnseg Right
model <- lm (Right_DnSeg ~ Amyloid..centiloids. + samseg_sbtiv ,
             data = dnseg_abcds_df)
summary(model)


#generate custom annotation to be placed at absolute 0.7,0.9
grob <-grobTree(textGrob(expression(atop(beta == -0.203, "p = 0.0115")), x=0.8, y=0.88, hjust = 0,
                         gp=gpar(col="red", fontsize =18, fontface="bold")))

png(file="Right nbM vs Amyloid.png", units="in", width=7, height=5, res=300)
ggplot(dnseg_abcds_df,aes(Amyloid..centiloids., Right_DnSeg)) +
  geom_point(size=2, color="blue") + 
  geom_smooth(method='lm', color="black")+
  labs(y=expression("Right nbM Volume, mm"^3), x= "Global Amyloid Accumulation (Centiloids)")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=18), axis.title=element_text(size=20))+
  annotation_custom(grob)
dev.off()

model <- lm (Right_DnSeg ~ Amyloid..centiloids. * Sex + samseg_sbtiv ,
             data = dnseg_abcds_df)
summary(model)

#ScLimbic Left
model <- lm (Left.Basal.Forebrain ~ Amyloid..centiloids. + samseg_sbtiv, 
             data = sclimbic_abcds_df)
summary(model)


#generate custom annotation to be placed at absolute 0.7,0.1
grob <-grobTree(textGrob(expression(atop(beta == -0.374, "p < 0.0001")), x=0.7, y=0.1, hjust = 0,
                         gp=gpar(col="red", fontsize =18, fontface="bold")))

png(file="Left ChBF vs Amyloid.png", units="in", width=7, height=5, res=300)
ggplot(sclimbic_abcds_df,aes(Amyloid..centiloids., Left.Basal.Forebrain)) +
  geom_point(size=2, color="blue") + 
  geom_smooth(method='lm', color="black")+
  labs(y=expression("Left ChBF Volume, mm"^3), x= "Global Amyloid Accumulation (Centiloids)")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=18), axis.title=element_text(size=20))+
  annotation_custom(grob)
dev.off()

model <- lm (Left.Basal.Forebrain ~ Amyloid..centiloids. * Sex + samseg_sbtiv, 
             data = sclimbic_abcds_df)
summary(model)

#ScLimbic Right
model <- lm (Right.Basal.Forebrain ~ Amyloid..centiloids. + samseg_sbtiv, 
             data = sclimbic_abcds_df)
summary(model)

#generate custom annotation to be placed at absolute 0.7,0.1
grob <-grobTree(textGrob(expression(atop(beta == -0.238, "p = 0.0011")), x=0.7, y=0.1, hjust = 0,
                         gp=gpar(col="red", fontsize =18, fontface="bold")))

png(file="Right ChBF vs Amyloid.png", units="in", width=7, height=5, res=300)
ggplot(sclimbic_abcds_df,aes(Amyloid..centiloids., Right.Basal.Forebrain)) +
  geom_point(size=2, color="blue") + 
  geom_smooth(method='lm', color="black")+
  labs(y=expression("Right ChBF Volume, mm"^3), x= "Global Amyloid Accumulation (Centiloids)")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=18), axis.title=element_text(size=20))+
  annotation_custom(grob)
dev.off()

model <- lm (Right.Basal.Forebrain ~ Amyloid..centiloids. * Sex + samseg_sbtiv, 
             data = sclimbic_abcds_df)
summary(model)


#generate custom annotation to be placed at absolute 0.7,0.1
grob <-grobTree(textGrob('sex x centiloid  interaction = 0.0126', x=0.01, y=0.03, hjust = 0,
                         gp=gpar(col="red", fontsize =18, fontface="bold")))

png(file="Right ChBF vs Age.png SEX.png", units="in", width=7, height=5, res=300)
ggplot(sclimbic_abcds_df,aes(Amyloid..centiloids., Right.Basal.Forebrain, col=Sex)) +
  geom_point(size=2) + 
  geom_smooth(method='lm')+
  labs(y=expression("Right ChBF Volume, mm"^3), x= "Global Amyloid Accumulation (Centiloids)")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=18), axis.title=element_text(size=20))+
  annotation_custom(grob)
dev.off()

# plot age as predictor for centiloid

model <- lm (Amyloid..centiloids. ~ Age, 
             data = dnseg_abcds_df)

R2 <- summary(model)$r.squared
R2 <- round(R2,digits= 4)
pval <- lmp(model)
pval <- round(pval,digits = 4)
R2 <- sprintf("%.4f", R2)
if(pval == 0){
  pval <- "<0.0001"
}

#generate custom annotation to be placed at absolute 0.1,0.9
grob <-grobTree(textGrob(glue('R² = {R2}\np = {pval}'), x=0.1, y=0.9, hjust = 0,
                         gp=gpar(col="red", fontsize =18)))

png(file="Amyloid vs Age.png", units="in", width=7, height=5, res=300)
ggplot(dnseg_abcds_df,aes(Age, Amyloid..centiloids.)) +
  geom_point(size=2, color="blue") + 
  geom_smooth(method='lm', color="black")+
  labs(y="Global Amyloid Accumulation (Centiloids)", x= "Age")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=18), axis.title=element_text(size=20))+
  annotation_custom(grob)
dev.off()

model <- lm (Amyloid..centiloids. ~ Age * Sex, 
             data = dnseg_abcds_df)
summary(model)

#plot age vs BF vol

#DnSeg Left

model <- lm (Left_DnSeg ~ Age + samseg_sbtiv, 
             data = dnseg_abcds_df)
summary(model)

#generate custom annotation to be placed at absolute 0.7,0.9
grob <-grobTree(textGrob(expression(atop(beta == -0.737, "p = 0.0068")), x=0.75, y=0.9, hjust = 0,
                         gp=gpar(col="red", fontsize =18, fontface="bold")))

png(file="Left nbM vs Age.png", units="in", width=7, height=5, res=300)
ggplot(dnseg_abcds_df,aes(Age, Left_DnSeg)) +
  geom_point(size=2, color="blue") + 
  geom_smooth(method='lm', color="black")+
  labs(y= expression("Left nbM Volume, mm"^3), x= "Age")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=18), axis.title=element_text(size=20))+
  annotation_custom(grob)
dev.off()

model <- lm (Left_DnSeg ~ Age * Sex + samseg_sbtiv, 
             data = dnseg_abcds_df)
summary(model)


#Dnseg Right
model <- lm (Right_DnSeg ~ Age + samseg_sbtiv, 
             data = dnseg_abcds_df)
summary(model)

#generate custom annotation to be placed at absolute 0.7,0.9
grob <-grobTree(textGrob(expression(atop(beta == -1.069, "p = 0.0007")), x=0.75, y=0.9, hjust = 0,
                         gp=gpar(col="red", fontsize =18, fontface="bold")))

png(file="Right nbM vs Age.png", units="in", width=7, height=5, res=300)
ggplot(dnseg_abcds_df,aes(Age, Right_DnSeg)) +
  geom_point(size=2, color="blue") + 
  geom_smooth(method='lm', color="black")+
  labs(y=expression("Right nbM Volume, mm"^3), x= "Age")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=18), axis.title=element_text(size=20))+
  annotation_custom(grob)
dev.off()

model <- lm (Right_DnSeg ~ Age * Sex +samseg_sbtiv, 
             data = dnseg_abcds_df)
summary(model)

#ScLimbic Left
model <- lm (Left.Basal.Forebrain ~ Age + samseg_sbtiv, 
             data = sclimbic_abcds_df)
summary(model)

#generate custom annotation to be placed at absolute 0.7,0.1
grob <-grobTree(textGrob(expression(atop(beta == -1.662, "p < 0.0001")), x=0.1, y=0.1, hjust = 0,
                         gp=gpar(col="red", fontsize =18, fontface="bold")))


png(file="Left ChBF vs Age.png", units="in", width=7, height=5, res=300)
ggplot(sclimbic_abcds_df,aes(Age, Left.Basal.Forebrain)) +
  geom_point(size=2, color="blue") + 
  geom_smooth(method='lm', color="black")+
  labs(y=expression("Left ChBF Volume, mm"^3), x= "Age")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=18), axis.title=element_text(size=20))+
  annotation_custom(grob)
dev.off()

model <- lm (Left.Basal.Forebrain ~ Age * Sex, 
             data = sclimbic_abcds_df)
summary(model)

#ScLimbic Right
model <- lm (Right.Basal.Forebrain ~ Age + samseg_sbtiv, 
             data = sclimbic_abcds_df)
summary(model)

#generate custom annotation to be placed at absolute 0.7,0.1
grob <-grobTree(textGrob(expression(atop(beta == -1.005, "p = 0.0005")), x=0.1, y=0.1, hjust = 0,
                         gp=gpar(col="red", fontsize =18, fontface="bold")))


png(file="Right ChBF vs Age.png", units="in", width=7, height=5, res=300)
ggplot(sclimbic_abcds_df,aes(Age, Right.Basal.Forebrain)) +
  geom_point(size=2, color="blue") + 
  geom_smooth(method='lm', color="black")+
  labs(y=expression("Right ChBF Volume, mm"^3), x= "Age")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=18), axis.title=element_text(size=20))+
  annotation_custom(grob)
dev.off()

model <- lm (Right.Basal.Forebrain ~ Age * Sex, 
             data = sclimbic_abcds_df)
summary(model)


age_avg <-mean(sclimbic_abcds_df$Age, na.rm = TRUE)
print(age_avg)

min_age <-min(sclimbic_abcds_df$Age, na.rm =TRUE)
print(min_age)

max_age <- max(sclimbic_abcds_df$Age, na.rm = TRUE)
print(max_age)

