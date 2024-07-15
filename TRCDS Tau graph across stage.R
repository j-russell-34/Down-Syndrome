library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(latex2exp)
library(ragg)
library(ggtext)
library(ggrepel)


braak_df <- read.csv("DSCHOL_MK6240_Braak6.csv")

braak_df$Centiloid <- c(-4.17, 8.71, 46.78,-0.92, 1.72,
                        38.46, -3.68, 2.28)

braak_long_df <- pivot_longer(
  data=braak_df,
  cols=c("stage1", "stage2", "stage3", "stage4", "stage5", "stage6"),
  names_to = "Braak_Stage",
  values_to = "SUVR"
)

pal <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00")


ggplot(braak_long_df, aes(x=Braak_Stage, y=SUVR, fill = Braak_Stage))+ 
  scale_fill_brewer(palette = "Greys") +
  stat_summary(geom="col", fun = "mean", width=0.7, color = "black", show.legend = FALSE) +
  geom_point(position = position_jitter(width = 0.17), aes(size = Centiloid, 
                                                          color = SUBJECT))+
  scale_size_continuous(range = c(2,6))+
  scale_color_brewer(palette = "Paired") +
  theme_hc() +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=12,face="bold")) +
  theme(axis.title.x = element_text(size=14,  margin=margin(t=20, r=0, b=0, l=0)),
        axis.title.y = element_text(size=14, margin=margin(t=0, r=20, b=0, l=0)))+
  stat_summary(fun.data=mean_sdl, 
               fun.args = list(mult=1), 
               geom="errorbar", width=.2)+
  scale_x_discrete(labels=c("stage1" = "Stage 1", "stage2" = "Stage 2",
                            "stage3" = "Stage 3", "stage4" = "Stage 4",
                            "stage5" = "Stage 5", "stage6" = "Stage 6")) +
  labs(x="Braak Stage", color = "Subject Amyloid\n (Centiloids)\n")+
  ylab(TeX(r"($\lbrack ^{18}F\rbrack $MK-6240 SUVR)", bold = TRUE)) +
  theme(legend.position = "right") + guides(fill=FALSE) + guides(size=FALSE) +
  scale_color_manual(labels = c("-4.17", "8.71", "46.78","-0.92", "1.72",
                                "38.46", "-3.68", "2.28"), values = pal)
  
 


ggsave("Tau_in_TRCDS.png", device=agg_png, width = 18, height = 15, units="cm", res=600)





  