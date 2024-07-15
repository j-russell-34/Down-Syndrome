library(tidyverse)
library(dplyr)
library(purrr)
library(glue)
library(readxl)
library(ggpp)
library(grid)
library(tibble)

#import table/excels to df
feobv_suvr_df <- read.table("DSCHOL_FEOBV_SUVR.txt",
  sep="\t", header=FALSE)
pib_suvr_df <- read.table("DSCHOL_PiB_SUVR.txt",
                            sep="\t", header=FALSE)
centiloid_df <- read_excel("TRC-DS centiloids 20231215.xlsx")

#apply headers to columns
colnames(feobv_suvr_df) <- c("SUBJECT", "REF_REGION", "PVC", "ROI", "SUVR")
colnames(pib_suvr_df) <- c("SUBJECT", "REF_REGION", "PVC", "ROI", "SUVR")

#change from long to wide
feobv_suvr_df_wide <- reshape(data = feobv_suvr_df,
                              idvar = "SUBJECT",
                              v.names = c("SUVR"),
                              timevar = "ROI",
                              direction = "wide")

pib_suvr_df_wide <- reshape(data = pib_suvr_df,
                              idvar = "SUBJECT",
                              v.names = c("SUVR"),
                              timevar = "ROI",
                              direction = "wide")

#add column to feobv_df for timepoint of scan relative to trc_ds BL
feobv_suvr_df_wide["VISIT"] <- c("Baseline", "Baseline", "M16","Baseline",
                                 "Baseline", "Baseline", "Baseline", "Baseline",
                                 "Baseline", "Baseline")

#merge dataframes based on visit and subject name
#when participants without a m16 amyloid scn will have to edit
feobv_suvr_df_wide <- merge(feobv_suvr_df_wide, centiloid_df, by=c("SUBJECT",
                                                                   "VISIT"))

#remove ROIs that are not of interest
feobv_df <- feobv_suvr_df_wide %>%
  select(-contains(c("White-Matter", "Cerebellum", "Brain-Stem", "CSF", 
        "choroid", "Air", "Skull", "Vermis", "Pons", "Head")))

#generate blank lists to add correlation analysis too
correlation_coef <- list()
correlation_pvalues <- list()
ROI <- list()

#function to calculate p-value of linear model
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

#convert - to . in column headers
names(feobv_df) <- make.names(names(feobv_df))

#iterate through list of ROIs, lm with centiloid value
for(i in 5:88){
  model <- lm (feobv_df[,c(i)] ~ Centiloid, data = feobv_df)
  correlation_coef <- append(correlation_coef, list(summary(model)$r.squared))
  correlation_pvalues <- append(correlation_pvalues, list(lmp(model)))
  ROI <- append(ROI, list(colnames(feobv_df)[i]))
}

#combine 3 lists to individual list and convert to single dataframe for correlations
df_list <- list(ROI, correlation_coef, correlation_pvalues)

correlation_df <- as.data.frame(do.call(cbind, df_list))

#apply headers to df
colnames(correlation_df)[1] <- "ROI"
colnames(correlation_df)[2] <- "R2"
colnames(correlation_df)[3] <- "pvalue"



#plot any regions with p value less than 0.1


#pull R2 and p from df and round
R2 <- correlation_df[which(correlation_df$ROI == "SUVR.Left.Thalamus"), "R2" ]
R2 <- round(R2[[1]],digits= 4)
pval <- correlation_df[which(correlation_df$ROI == "SUVR.Left.Thalamus"), "pvalue"]
pval <- round(pval[[1]],digits = 4)
R2 <- sprintf("%.4f", R2)
pval <- sprintf("%.4f", pval)
#generate custom annotation to be placed at absolute 0.1,0.1
grob <-grobTree(textGrob(glue('R² = {R2}\np = {pval}'), x=0.1, y=0.1, hjust = 0,
                         gp=gpar(col="red", fontsize =18, fontface="bold")))

png(file="Amyloid vs L Thalamus FEOBV.png", units="in", width=7, height=5, res=300)
ggplot(feobv_df,aes(Centiloid, SUVR.Left.Thalamus)) +
  geom_point(size=6, color="blue") + 
  geom_smooth(method='lm', color="black")+
  labs(y="Left Thalamus FEOBV SUVR", x= "Global Amyloid Accumulation (Centiloids)")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
  axis.text=element_text(size=16), axis.title=element_text(size=18, face="bold"))+
  annotation_custom(grob)
dev.off()


R2 <- correlation_df[which(correlation_df$ROI == "SUVR.Left.Amygdala"), "R2" ]
R2 <- round(R2[[1]],digits= 4)
pval <- correlation_df[which(correlation_df$ROI == "SUVR.Left.Amygdala"), "pvalue"]
pval <- round(pval[[1]],digits = 4)  
R2 <- sprintf("%.4f", R2)
pval <- sprintf("%.4f", pval)
grob <-grobTree(textGrob(glue('R² = {R2}\np = {pval}'), x=0.1, y=0.1, hjust = 0,
                         gp=gpar(col="red", fontsize =18, fontface="bold")))

png(file="Amyloid vs L Amygdala FEOBV.png", units="in", width=7, height=5, res=300)
ggplot(feobv_df,aes(Centiloid, SUVR.Left.Amygdala)) +
  geom_point(size=6, color="blue") + 
  geom_smooth(method='lm', color="black")+
  labs(y="Left Amygdala FEOBV SUVR", x= "Global Amyloid Accumulation (Centiloids)")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
  axis.text=element_text(size=16), axis.title=element_text(size=18, face="bold"))+
  annotation_custom(grob)
dev.off()
  
  
R2 <- correlation_df[which(correlation_df$ROI == "SUVR.Right.Thalamus"), "R2" ]
R2 <- round(R2[[1]],digits= 4)
pval <- correlation_df[which(correlation_df$ROI == "SUVR.Right.Thalamus"), "pvalue"]
pval <- round(pval[[1]],digits = 4) 
R2 <- sprintf("%.4f", R2)
pval <- sprintf("%.4f", pval)
grob <-grobTree(textGrob(glue('R² = {R2}\np = {pval}'), x=0.1, y=0.1, hjust = 0,
                         gp=gpar(col="red", fontsize =18, fontface="bold")))

png(file="Amyloid vs R Thal FEOBV.png", units="in", width=7, height=5, res=300)
ggplot(feobv_df,aes(Centiloid, SUVR.Right.Thalamus)) +
  geom_point(size=6, color="blue") + 
  geom_smooth(method='lm', color="black")+
  labs(y="Right Thalamus FEOBV SUVR", x= "Global Amyloid Accumulation (Centiloids)")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
  axis.text=element_text(size=16), axis.title=element_text(size=18, face="bold"))+
  annotation_custom(grob)
dev.off()
  
R2 <- correlation_df[which(correlation_df$ROI == "SUVR.Right.Amygdala"), "R2" ]
R2 <- round(R2[[1]],digits= 4)
pval <- correlation_df[which(correlation_df$ROI == "SUVR.Right.Amygdala"), "pvalue"]
pval <- round(pval[[1]],digits = 4)  
R2 <- sprintf("%.4f", R2)
pval <- sprintf("%.4f", pval)
grob <-grobTree(textGrob(glue('R² = {R2}\np = {pval}'), x=0.1, y=0.1, hjust = 0,
                         gp=gpar(col="red", fontsize =18, fontface="bold")))

png(file="Amyloid vs R amygdala.png", units="in", width=7, height=5, res=300)
ggplot(feobv_df,aes(Centiloid, SUVR.Right.Amygdala)) +
  geom_point(size=6, color="blue") + 
  geom_smooth(method='lm', color="black")+
  labs(y="Right Amygdala FEOBV SUVR", x= "Global Amyloid Accumulation (Centiloids)")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
  axis.text=element_text(size=16), axis.title=element_text(size=18, face="bold"))+
  annotation_custom(grob)
dev.off()

R2 <- correlation_df[which(correlation_df$ROI == "SUVR.ctx.lh.entorhinal"), "R2" ]
R2 <- round(R2[[1]],digits= 4)
pval <- correlation_df[which(correlation_df$ROI == "SUVR.ctx.lh.entorhinal"), "pvalue"]
pval <- round(pval[[1]],digits = 4)
R2 <- sprintf("%.4f", R2)
pval <- sprintf("%.4f", pval)
grob <-grobTree(textGrob(glue('R² = {R2}\np = {pval}'), x=0.1, y=0.1, hjust = 0,
                         gp=gpar(col="red", fontsize =18, fontface="bold")))

png(file="Amyloid vs l entorhinal FEOBV.png", units="in", width=7, height=5, res=300)
ggplot(feobv_df,aes(Centiloid, SUVR.ctx.lh.entorhinal)) +
  geom_point(size=6, color="blue") + 
  geom_smooth(method='lm', color="black")+
  labs(y="Left Entorhinal FEOBV SUVR", x= "Global Amyloid Accumulation (Centiloids)")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=16), axis.title=element_text(size=18, face="bold"))+
  annotation_custom(grob)
dev.off()



#associations between amyloid SUVR and FEOBV SUVR
#convert - to . in column headers
names(pib_suvr_df_wide) <- make.names(names(pib_suvr_df_wide))

#remove ROIs that are not of interest
pib_suvr_df_wide <- pib_suvr_df_wide %>%
  select(-contains(c("White.Matter", "Cerebellum", "Brain.Stem", "CSF", 
                     "choroid", "Air", "Skull", "Vermis", "Pons", "Head")))

visit_list = list("Baseline","M16", "Baseline", "Baseline","Baseline",
                  "Baseline", "Baseline", "Baseline", "Baseline")

pib_suvr_df <- add_column(pib_suvr_df_wide, Visit = visit_list, .before = "REF_REGION")

#remove participants that do not have an FEOBV scan
pib_suvr_df <- pib_suvr_df[ pib_suvr_df$SUBJECT %in% feobv_df$SUBJECT, ]

#generate blank lists to add correlation analysis too
correlation_coef <- list()
correlation_pvalues <- list()
ROI <- list()

#iterate through list of ROIs, lm with centiloid value
for(i in 5:88){
  model <- lm (feobv_df[,c(i)] ~ pib_suvr_df[,c(i)])
  correlation_coef <- append(correlation_coef, list(summary(model)$r.squared))
  correlation_pvalues <- append(correlation_pvalues, list(lmp(model)))
  ROI <- append(ROI, list(colnames(feobv_df)[i]))
}

#combine 3 lists to individual list and convert to single dataframe for correlations
df_list <- list(ROI, correlation_coef, correlation_pvalues)

correlation_df_PiB <- as.data.frame(do.call(cbind, df_list))

#apply headers to df
colnames(correlation_df_PiB)[1] <- "ROI"
colnames(correlation_df_PiB)[2] <- "R2"
colnames(correlation_df_PiB)[3] <- "pvalue"

#plot any regions with p value less than 0.1
#Right Amygdala, SUVR.ctx.lh.fusiform, L inferior temporal, L preicalcarine, R preicalcarine

#pull R2 and p from df and round
R2 <- correlation_df_PiB[which(correlation_df_PiB$ROI == "SUVR.Right.Amygdala"), "R2" ]
R2 <- round(R2[[1]],digits= 4)
pval <- correlation_df_PiB[which(correlation_df_PiB$ROI == "SUVR.Right.Amygdala"), "pvalue"]
pval <- round(pval[[1]],digits = 4)
R2 <- sprintf("%.4f", R2)
pval <- sprintf("%.4f", pval)
#generate custom annotation to be placed at absolute 0.1,0.1
grob <-grobTree(textGrob(glue('R² = {R2}\np = {pval}'), x=0.1, y=0.1, hjust = 0,
                         gp=gpar(col="red", fontsize =18, fontface="bold")))

png(file="R Amygdala PiB vs FEOBV.png", units="in", width=7, height=5, res=300)
ggplot(data=data.frame(x=pib_suvr_df$"SUVR.Right.Amygdala",
                       y=feobv_df$"SUVR.Right.Amygdala"),
       aes(x=x, y=y)) +
  geom_point(size=6, color="blue") + 
  geom_smooth(method='lm', color="black")+
  labs(y="Right Amygdala FEOBV SUVR", x= "Right Amygdala PiB SUVR")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=16), axis.title=element_text(size=18, face="bold"))+
  annotation_custom(grob)
dev.off()

R2 <- correlation_df_PiB[which(correlation_df_PiB$ROI == "SUVR.ctx.lh.fusiform"), "R2" ]
R2 <- round(R2[[1]],digits= 4)
pval <- correlation_df_PiB[which(correlation_df_PiB$ROI == "SUVR.ctx.lh.fusiform"), "pvalue"]
pval <- round(pval[[1]],digits = 4)
R2 <- sprintf("%.4f", R2)
pval <- sprintf("%.4f", pval)
#generate custom annotation to be placed at absolute 0.1,0.1
grob <-grobTree(textGrob(glue('R² = {R2}\np = {pval}'), x=0.1, y=0.1, hjust = 0,
                         gp=gpar(col="red", fontsize =18, fontface="bold")))

png(file="L Fusiform PiB vs FEOBV.png", units="in", width=7, height=5, res=300)
ggplot(data=data.frame(x=pib_suvr_df$"SUVR.ctx.lh.fusiform",
                       y=feobv_df$"SUVR.ctx.lh.fusiform"),
       aes(x=x, y=y)) +
  geom_point(size=6, color="blue") + 
  geom_smooth(method='lm', color="black")+
  labs(y="Left Fusiform FEOBV SUVR", x= "Left Fusiform PiB SUVR")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=16), axis.title=element_text(size=18, face="bold"))+
  annotation_custom(grob)
dev.off()

R2 <- correlation_df_PiB[which(correlation_df_PiB$ROI == "SUVR.ctx.lh.inferiortemporal"), "R2" ]
R2 <- round(R2[[1]],digits= 4)
pval <- correlation_df_PiB[which(correlation_df_PiB$ROI == "SUVR.ctx.lh.inferiortemporal"), "pvalue"]
pval <- round(pval[[1]],digits = 4)
R2 <- sprintf("%.4f", R2)
pval <- sprintf("%.4f", pval)
#generate custom annotation to be placed at absolute 0.1,0.1
grob <-grobTree(textGrob(glue('R² = {R2}\np = {pval}'), x=0.1, y=0.1, hjust = 0,
                         gp=gpar(col="red", fontsize =18, fontface="bold")))

png(file="L Inferiortemporal PiB vs FEOBV.png", units="in", width=7, height=5, res=300)
ggplot(data=data.frame(x=pib_suvr_df$"SUVR.ctx.lh.inferiortemporal",
                       y=feobv_df$"SUVR.ctx.lh.inferiortemporal"),
       aes(x=x, y=y)) +
  geom_point(size=6, color="blue") + 
  geom_smooth(method='lm', color="black")+
  labs(y="Left Inferior Temporal FEOBV SUVR", x= "Left Inferior Temporal PiB SUVR")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=16), axis.title=element_text(size=18, face="bold"))+
  annotation_custom(grob)
dev.off()

R2 <- correlation_df_PiB[which(correlation_df_PiB$ROI == "SUVR.ctx.lh.pericalcarine"), "R2" ]
R2 <- round(R2[[1]],digits= 4)
pval <- correlation_df_PiB[which(correlation_df_PiB$ROI == "SUVR.ctx.lh.pericalcarine"), "pvalue"]
pval <- round(pval[[1]],digits = 4)
R2 <- sprintf("%.4f", R2)
pval <- sprintf("%.4f", pval)
#generate custom annotation to be placed at absolute 0.1,0.1
grob <-grobTree(textGrob(glue('R² = {R2}\np = {pval}'), x=0.1, y=0.1, hjust = 0,
                         gp=gpar(col="red", fontsize =18, fontface="bold")))

png(file="L Pericalcarine PiB vs FEOBV.png", units="in", width=7, height=5, res=300)
ggplot(data=data.frame(x=pib_suvr_df$"SUVR.ctx.lh.pericalcarine",
                       y=feobv_df$"SUVR.ctx.lh.pericalcarine"),
       aes(x=x, y=y)) +
  geom_point(size=6, color="blue") + 
  geom_smooth(method='lm', color="black")+
  labs(y="Left Pericalcarine FEOBV SUVR", x= "Left Pericalcarine PiB SUVR")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=16), axis.title=element_text(size=18, face="bold"))+
  annotation_custom(grob)
dev.off()

R2 <- correlation_df_PiB[which(correlation_df_PiB$ROI == "SUVR.ctx.rh.pericalcarine"), "R2" ]
R2 <- round(R2[[1]],digits= 4)
pval <- correlation_df_PiB[which(correlation_df_PiB$ROI == "SUVR.ctx.rh.pericalcarine"), "pvalue"]
pval <- round(pval[[1]],digits = 4)
R2 <- sprintf("%.4f", R2)
pval <- sprintf("%.4f", pval)
#generate custom annotation to be placed at absolute 0.1,0.1
grob <-grobTree(textGrob(glue('R² = {R2}\np = {pval}'), x=0.1, y=0.1, hjust = 0,
                         gp=gpar(col="red", fontsize =18, fontface="bold")))

png(file="R Pericalcarine PiB vs FEOBV.png", units="in", width=7, height=5, res=300)
ggplot(data=data.frame(x=pib_suvr_df$"SUVR.ctx.rh.pericalcarine",
                       y=feobv_df$"SUVR.ctx.rh.pericalcarine"),
       aes(x=x, y=y)) +
  geom_point(size=6, color="blue") + 
  geom_smooth(method='lm', color="black")+
  labs(y="Right Pericalcarine FEOBV SUVR", x= "Right Pericalcarine PiB SUVR")+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=16), axis.title=element_text(size=18, face="bold"))+
  annotation_custom(grob)
dev.off()
