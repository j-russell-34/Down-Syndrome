library(tidyverse)
library(dplyr)
library(purrr)
library(glue)
library(readxl)
library(ggpp)
library(grid)
library(glmnet)

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

dnseg_abcds_df <- merge(dnseg_abcds_df, dnsegdata_df, by="SUBJECT")
dnseg_abcds_df <- merge(dnseg_abcds_df, samseg_df, by="SUBJECT")
dnseg_abcds_df <- dnseg_abcds_df[!duplicated(dnseg_abcds_df$SUBJECT),]

dnseg_abcds_df$totalch4avg <- (dnseg_abcds_df$Left_DnSeg + 
                                 dnseg_abcds_df$Left_DnSeg)/2

dnseg_abcds_df$Sex <- factor(dnseg_abcds_df$Gender)
dnseg_abcds_df$avgch4divtiv <- dnseg_abcds_df$totalch4avg/dnseg_abcds_df$samseg_sbtiv

genderencoded <- model.matrix(~ Gender - 1, data = sclimbic_abcds_df)
genderencoded_df <- as.data.frame(genderencoded)
genderencoded_df$SUBJECT <- sclimbic_abcds_df$SUBJECT

sclimbic_abcds_df <- merge(sclimbic_abcds_df, genderencoded_df, by="SUBJECT")

#LASSO for sclimbic

y <- sclimbic_abcds_df$totalBFavg
x <- data.matrix(sclimbic_abcds_df[,c('Age', 'Amyloid..centiloids.', 'Sex',
                                  'samseg_sbtiv')])


#perform k-fold cross-validation to find optimal lambda value
cv_model <- cv.glmnet(x, y, alpha = 1)
best_lambda <- cv_model$lambda.min
best_lambda

#produce plot of test MSE by lambda value
#plot(cv_model) 

best_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)
coef(best_model)

#use fitted best model to make predictions
y_predicted <- predict(best_model, s = best_lambda, newx = x)


#find SST and SSE
sst <- sum((y - mean(y))^2)
sse <- sum((y_predicted - y)^2)

#find R-Squared
rsq <- 1 - sse/sst
rsq




#LASSO for dnseg

y1 <- dnseg_abcds_df$totalch4avg
x1 <- data.matrix(dnseg_abcds_df[,c('Age', 'Amyloid..centiloids.', 'Sex',
                                      'samseg_sbtiv')])


#perform k-fold cross-validation to find optimal lambda value
cv_model <- cv.glmnet(x1, y1, alpha = 1)
best_lambda1 <- cv_model$lambda.min
best_lambda1

#produce plot of test MSE by lambda value
#plot(cv_model) 

best_model1 <- glmnet(x1, y1, alpha = 1, lambda = best_lambda1)
coef(best_model1)

#use fitted best model to make predictions
y_predicted1 <- predict(best_model1, s = best_lambda1, newx = x1)


#find SST and SSE
sst1 <- sum((y1 - mean(y1))^2)
sse1 <- sum((y_predicted1 - y1)^2)

#find R-Squared
rsq <- 1 - sse1/sst1
rsq
