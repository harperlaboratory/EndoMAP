## Title: EndoMAP
## Author: Miguel A. Gonzalez-Lozano
## Date: 08/14/2024
  

library (RColorBrewer)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(FactoMineR)
library(factoextra)
library(pROC)
library(reshape2)

###
### Endosomal Meta-analysis
###


### Multiple Correspondence Analysis (MCA) - Extended Data Figure 1b
file_Expression = ".../Endosome_MetaAnalysis.csv"
t = read.csv(file_Expression, header=T, stringsAsFactors=F, skipNul = TRUE)
t = t[,2:ncol(t)] #Remove gene names

file_Expression = ".../Endosome_MetaAnalysis_Metadata.csv"
t.meta = read.csv(file_Expression, header=T, stringsAsFactors=F, skipNul = TRUE)

res.mca <- MCA(t, graph = FALSE)
name <- list(name = paste0(colnames(t),"_TRUE"))

b = fviz_mca_var(res.mca, col.var = "contrib", select.var = name,
             repel = TRUE,
             ggtheme = theme_minimal()
)
b = b$data
b$name = gsub('_TRUE','',b$name) #Remove "_TRUE" from dataset names

Nprots = as.data.frame(colSums(t,na.rm = TRUE)) #Count and add total number of proteins in each dataset
Nprots = cbind(Nprots, rownames(Nprots))
colnames(Nprots) = c("Nproteins","name")
b = merge(b, Nprots, by = 'name', all.x = T)

colnames(t.meta)[1] = "name"
b = merge(b, t.meta, by = 'name', all.x = T) #Add metadata

ggplot(b, aes(x=x, y=y)) + 
  xlab("Dim1.(11.9%)") + ylab("Dim2.(11.8%)")+
  scale_size(range = c(10,40))+
  geom_point(aes(size=Nproteins, fill=Isolation.method), shape=21)+ 
  scale_fill_brewer(palette="Set3")+
  geom_text_repel(label=b$name) +
  geom_vline(xintercept = 0, linetype="dashed", color = "black", size=0.5)+
  geom_hline(yintercept = 0, linetype="dashed", color = "black", size=0.5)+
  theme_minimal()


#### Receiver operating characteristic (ROC) curves - Figure 1b
file_Expression = "C:/Users/migue/Desktop/USA/projects/Endosome data analysis/Database/Figures/Endosome_Database_Merged_ROC_PPIs.csv"
t = read.csv(file_Expression, header=T, stringsAsFactors=F, skipNul = TRUE)

# Binomial logistic regression
glm.fit = glm(as.numeric(t$Endo) ~ as.numeric(t$Count), family = binomial) # Dataset count

t.crop = t
t.crop[is.na(t.crop$Abundance.EndoIP),"Abundance.EndoIP"] = 0 # Convert NAs to ceros
glm.fit2 = glm(as.numeric(t.crop$Endo) ~ as.numeric(t.crop$Abundance.EndoIP), family = binomial) #Abundance in Endo-IP

t.comb = t.crop
glm.fit.S = glm(as.numeric(t.comb$Endo) ~ as.numeric(t.comb$Count) + as.numeric(t.comb$Abundance.EndoIP) + as.numeric(t.comb$PPIs), family = binomial) # Combined score
glm.fit.xl = glm(as.numeric(t$Endo) ~ as.numeric(t$PPIs), family = binomial) # Endosomal protein interactions

# Plot
par(pty = "s")
roc(t$Endo, glm.fit$fitted.values, plot = TRUE, legacy.axes = T, percent = T,
    xlab = "False Positive percentage", ylab="True Positive percentage",
    col=brewer.pal(4,"Dark2"), lwd=4, print.auc = T, partial.auc = c(100,90))
roc(t.crop$Endo, glm.fit2$fitted.values, plot = TRUE, legacy.axes = T, percent = T,
    xlab = "False Positive percentage", ylab="True Positive percentage",
    col=brewer.pal(4,"Dark2")[2], lwd=4, print.auc = T, print.auc.y = 45, partial.auc = c(100,90), add=T)
roc(t.comb$Endo, glm.fit.S$fitted.values, plot = TRUE, legacy.axes = T, percent = T,
    xlab = "False Positive percentage", ylab="True Positive percentage",
    col=brewer.pal(4,"Dark2")[3], lwd=4, print.auc = T, print.auc.y = 35, partial.auc = c(100,90),add=T)
roc(t$Endo, glm.fit.xl$fitted.values, plot = TRUE, legacy.axes = T, percent = T,
    xlab = "False Positive percentage", ylab="True Positive percentage",
    col=brewer.pal(4,"Dark2")[4], lwd=4, print.auc = T, print.auc.y = 40, partial.auc = c(100,90),add=T)

legend("bottomright", legend=c("Paper count", "Abundance",  "PPIs", "Combined Score"),
       col=brewer.pal(4,"Dark2")[c(1,2,4,3)], lwd=4)

par(pty = "m")

score = data.frame(t.comb, score = glm.fit.S$fitted.values) # Make table with results
score = score[order(score$score, decreasing = F),]
score$rank = 1:nrow(score)

ggplot(data=score, aes(x = rank, y = score))+
  geom_point(aes(color=Endo))+
  # geom_text_repel(data = subset(score, Group > 0 & score > 0.0734),
  #                 aes(label = Gene), size = 4,
  #                 max.overlaps = 20, min.segment.length = 0) +
  geom_hline(yintercept = 0.0734, linetype="dashed", color = "black", size=0.5)+
  ylab("Combined score")+ xlab("Rank")+
  theme_minimal()
