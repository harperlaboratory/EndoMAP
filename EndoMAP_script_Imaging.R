## Title: EndoMAP
## Author: Miguel A. Gonzalez-Lozano
## Date: 08/14/2024
  
library(viridis)
library(gridExtra)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(lme4)

###
### Co-localization and statistical analysis
###

file_Expression = ".../IF_sum159.csv"
t = read.csv(file_Expression, header=T, stringsAsFactors=F, skipNul = TRUE) 
t$Prots = paste0(t$Prot1, "-", t$Prot2)
t$Reps = as.character(t$Rep)

#Filter data based on area after threshold
t = t[t$Thresholded.M2 != 0,]
t = t[t$Area.A >5 & t$Area.B>5,]

#Stats
t.summary.f <- t %>%
  group_by(Prots) %>%
  summarise(
    sd.M2 = sd(Thresholded.M2, na.rm = TRUE),
    Num = length(na.omit(Thresholded.M2)),
    sem.M2 = sd.M2 / Num,
    Thresholded.M2 = mean(Thresholded.M2, na.rm = TRUE)
  )

P = combn(unique(t$Prots),2,simplify = T) #all pair comparisons
t$Well = paste0(t$Prots,t$Sample,t$Rep) # To separate wells as independent

Pval.n.w = data.frame()
for (i in 1:ncol(P)) { # Stats for all pairwise comparisons
  
  full.lmer <- lmer(Thresholded.M2 ~ Prots + (1|Rep/Well), 
                    data = t[t$Prots == P[1,i] | t$Prots == P[2,i],], REML = FALSE) #Nested wells
  reduced.lmer <- lmer(Thresholded.M2 ~ 1 + (1|Rep/Well), 
                       data = t[t$Prots == P[1,i] | t$Prots == P[2,i],], REML = FALSE)
  sig.M2 = anova(reduced.lmer, full.lmer)
  Pval.n.w = rbind(Pval.n.w, c(P[1,i],P[2,i],sig.M2$`Pr(>Chisq)`[2]))
}
colnames(Pval.n.w) = c("Prot1", "Prot2", "pval")
Pval.n.w$pval = as.numeric(Pval.n.w$pval)

plot = ggplot(aes(x = factor(Prots, level = c("TMEM9-CLCN3", "CLCN3-TMEM9", "TMEM9-EEA1", "CLCN3-EEA1", "TMEM9-LAMP1", "CLCN3-LAMP1")), y = Thresholded.M2), data = t) + 
  stat_summary(
    aes(color = Reps, width = 0.4), fun.data="mean_se", fun.args = list(mult=1), 
    linewidth = 0.4, position = position_dodge(0.6), geom = "errorbar")+
  stat_summary(
    aes(color = Reps, width = 0.6), fun.y="mean", fun.args = list(mult=1),
    linewidth = 0.5, position = position_dodge(0.6), geom = "crossbar")+
  scale_color_manual(values=viridis(4))+
  ylab("M2")+ 
  xlab("")+ 
  theme_bw() + 
  theme(axis.text.x = element_text(color="black",angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(color="black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
  geom_text(aes(label= paste0("n=",Num)), data = t.summary.f, angle = 0, vjust = 0, position = position_fill(vjust = 0)) +
  geom_signif(comparisons = list(c(P[1,1], P[2,1]),c(P[1,2], P[2,2]),c(P[1,6], P[2,6])),
              annotations= c(round(Pval.n.w[Pval.n.w$Prot1 == P[1,1] & Pval.n.w$Prot2 == P[2,1],"pval"],3),
                             round(Pval.n.w[Pval.n.w$Prot1 == P[1,2] & Pval.n.w$Prot2 == P[2,2],"pval"],3),
                             round(Pval.n.w[Pval.n.w$Prot1 == P[1,6] & Pval.n.w$Prot2 == P[2,6],"pval"],3)),
              tip_length = 0, vjust=0, y_position = c(1, 1.1, 0.9))+
  ggtitle("MixedModel nested wells indep.")

