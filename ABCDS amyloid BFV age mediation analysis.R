
dnseg_abcds_df <- read.csv("ABC-DS DnSeg DF for R 20231221.csv")
sclimbic_abcds_df <- read.csv("ABC-DS ScLimbic DF for R 20231221.csv")

library(mediation)

#mediation for left nBM

model.age.lnbm <- lm(Left.nbM.Volume.normalized.to.ICV ~ Age, data = dnseg_abcds_df)
summary(model.age.lnbm)

model.age.amyloid <- lm(Amyloid..centiloids. ~ Age, data = dnseg_abcds_df)
summary(model.age.amyloid)

model.ageamyloid.lnbm <- lm(Left.nbM.Volume.normalized.to.ICV ~ Age + 
                             Amyloid..centiloids., data = dnseg_abcds_df)
summary(model.ageamyloid.lnbm)

lnbm.mediation.results <- mediate(model.age.amyloid, model.ageamyloid.lnbm, treat="Age",  
                             mediator="Amyloid..centiloids.", boot=TRUE, sims=500)

summary(lnbm.mediation.results)


#repeat for right nBM

model.age.rnbm <- lm(Right.nbM.Volume.normalized.to.ICV ~ Age, data = dnseg_abcds_df)
summary(model.age.rnbm)

model.age.amyloid <- lm(Amyloid..centiloids. ~ Age, data = dnseg_abcds_df)
summary(model.age.amyloid)

model.ageamyloid.rnbm <- lm(Right.nbM.Volume.normalized.to.ICV ~ Age + 
                              Amyloid..centiloids., data = dnseg_abcds_df)
summary(model.ageamyloid.rnbm)

rnbm.mediation.results <- mediate(model.age.amyloid, model.ageamyloid.rnbm, treat="Age",  
                                  mediator="Amyloid..centiloids.", boot=TRUE, sims=500)

summary(rnbm.mediation.results)

#repeat for left chbf

model.age.lchbf <- lm(Left.ChBF.Volume.normalized.to.ICV ~ Age, data = sclimbic_abcds_df)
summary(model.age.lchbf)

model.age.amyloid <- lm(Amyloid..centiloids. ~ Age, data = sclimbic_abcds_df)
summary(model.age.amyloid)

model.ageamyloid.lchbf <- lm(Left.ChBF.Volume.normalized.to.ICV ~ Age + 
                               Amyloid..centiloids., data = sclimbic_abcds_df)
summary(model.ageamyloid.lchbf)

lchbf.mediation.results <- mediate(model.age.amyloid, model.ageamyloid.lchbf, treat="Age",  
                                   mediator="Amyloid..centiloids.", boot=TRUE, sims=500)

summary(lchbf.mediation.results)

#repeat for right chbf

model.age.rchbf <- lm(Right.ChBF.Volume.normalized.to.ICV ~ Age, data = sclimbic_abcds_df)
summary(model.age.rchbf)

model.age.amyloid <- lm(Amyloid..centiloids. ~ Age, data = sclimbic_abcds_df)
summary(model.age.amyloid)

model.ageamyloid.rchbf <- lm(Right.ChBF.Volume.normalized.to.ICV ~ Age + 
                              Amyloid..centiloids., data = sclimbic_abcds_df)
summary(model.ageamyloid.rchbf)

rchbf.mediation.results <- mediate(model.age.amyloid, model.ageamyloid.rchbf, treat="Age",  
                                  mediator="Amyloid..centiloids.", boot=TRUE, sims=500)

summary(rchbf.mediation.results)