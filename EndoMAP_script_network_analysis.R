## Title: EndoMAP
## Author: Miguel A. Gonzalez-Lozano
## Date: 08/14/2024
  

library(igraph)
library(reshape2)
library(dplyr)
library (RColorBrewer)

###
### EndoMAP network analysis
###


# Load data
file_Expression = ".../Network.csv"
t = read.csv(file_Expression, header=T, stringsAsFactors=F, skipNul = TRUE)

# Filtering
t = t[t$Final_annot != "Nucleus" & t$Final_annot.1 != "Nucleus",]# Remove dubious interactions

t = t[!(t$Gene1.order %in% c("EEA1", "UBC", "FTH1", "FTL", "LDHA", "LDHB", "ALB")) & # Exclude contaminants
        !(t$Gene2.order %in% c("EEA1", "UBC", "FTH1", "FTL", "LDHA", "LDHB", "ALB")),]
t = t[!grepl("KRT",t$Gene1.order) & !grepl("KRT",t$Gene1.order),]

# Inclusion list
l.incl = t[t$Final_annot %in% c("Endosome") | t$Final_annot.1 %in% c("Endosome"), c("Gene1.order","Gene2.order")]
l.incl = unique(c(l.incl$Gene1.order, l.incl$Gene2.order)) #Include all interactions between proteins
l.incl1 = l.incl
# Second interactor
l.incl = t[t$Gene1.order %in% l.incl | t$Gene2.order %in% l.incl, c("Gene1.order","Gene2.order")]
l.incl = unique(c(l.incl$Gene1.order, l.incl$Gene2.order))
l.incl2 = l.incl
# Prune second interactors
l.excl = setdiff(l.incl2, l.incl1) # Second layer of proteins
t.prune = t[t$Gene1.order %in% l.incl & t$Gene2.order %in% l.incl,]
t.prune$dint = t.prune$Gene1.order %in% l.incl1 | t.prune$Gene2.order %in% l.incl1 # Direct endosomal interactors
t.prune = t.prune[t.prune$Gene1.order %in% l.excl | t.prune$Gene2.order %in% l.excl,]
t.prune = t.prune[t.prune$dint, ] # Second layer of proteins connected to direct endo interactors

l.incl3 = data.frame()
for (i in l.excl) {
  tmp = t.prune[t.prune$Gene1.order == i | t.prune$Gene2.order == i,]
  n.inter = sum(tmp$Crosslink.Type == "Inter", na.rm = TRUE) > 0
  n.bn = sum(tmp$Crosslink.Type == "BN", na.rm = TRUE) > 1
  l.incl3 = rbind(l.incl3, cbind(i, n.inter | n.bn))
}
l.incl3 = l.incl3[l.incl3$V2 == TRUE, 1]
l.incl = c(l.incl1, l.incl3)

t = t[t$Gene1.order %in% l.incl & t$Gene2.order %in% l.incl,]

# Generate network
t.filter = t[, c("Gene1.order","Gene2.order")]
tg= graph_from_data_frame(t.filter, directed = F)
tg = simplify(tg, remove.multiple = T, remove.loops = T)
tg = delete.vertices(tg, which(degree(tg)<1)) # remove nodes that have no edges


### Degree distribution main component - Extended Figure 3a
vert_ids = V(tg)[components(tg)$membership == 1] # Select only main component of the network
tg.main = induced_subgraph(tg, vert_ids)

deg = igraph::degree(tg.main, mode="all")
deg.hist = deg
deg.hist[deg.hist > 10] = 10
p1 = hist(deg.hist, breaks=c(0:10), main="Histogram of node degree",
          xlim = c(0,10), xaxt = "n",
          xlab = "Node Degree") # Histogram until degree 10
axis(side = 1, at = 1:10 - 0.5, labels=c(1:9,">10"))

### Extended Figure 3b
deg.dist = degree_distribution(tg.main, cumulative=F, mode="all")
plot(x = 0:(length(deg.dist)-1), # degree_distribution starts at degree 0 (no edges)
     y = deg.dist,
     log = "xy",
     xlab = "Node Degree",
     ylab = "Probability",
     main = "Main component of the network")
pf = power.law.fit(deg.dist, impelementation = "plfit") # Power law fit
lines(0:(length(deg.dist)-1), 
      (0:(length(deg.dist)-1))^-pf$alpha, col="#b00606")


### Mean distance main component - Extended Figure 3c
dist = data.frame(dist = distance_table(tg.main, directed = F)$res,
                  Short = 1:length(distance_table(tg.main, directed = F)$res))
barplot(dist$dist, names.arg = dist$Short,
        xlab = "Shortest path distance", ylab = "Frequency",
        main = "Main component of the network")
abline(v = mean_distance(tg.main,directed = FALSE), 
       col="black", lwd=2, lty=2)
legend("topright", inset = 0.02,
       legend=c(paste0("Average shortest path distance: ", round(mean_distance(tg.main,directed = FALSE),1))),
       cex=0.8, box.lty=0)


### Distances between annotations - Figure 2d
dist = distances(tg.main, weights = NA , mode = "all")
dist[upper.tri(dist, diag = FALSE)] = 0 # Get only half of the corr matrix
dist = melt(dist)
dist = dist[dist$value != 0,]

Annot1 = subset(t, select = c("Gene1.order", "Complex1")) #Get annotations
Annot2 = subset(t, select = c("Gene2.order", "Complex1.1"))
colnames(Annot2) = colnames(Annot1)
Annot = rbind(Annot1, Annot2)
Annot = Annot[!duplicated(Annot),]
rm(Annot1, Annot2)
dist = left_join(dist, Annot, by = c("Var1" = "Gene1.order"))
dist = left_join(dist, Annot, by = c("Var2" = "Gene1.order"))

h.dist.s = dist[dist$Complex1.x == dist$Complex1.y & dist$Complex1.x != "",] # Distance between annotated
h.dist.d = dist[dist$Complex1.x != dist$Complex1.y & dist$Complex1.x != "",] # Distance between different annotations
h.dist.b = dist[dist$Complex1.x == dist$Complex1.y & dist$Complex1.x == "",] # Distance between not annotated

h.dist.s[h.dist.s$value > 10, "value"] = 10 # Cap max distance
h.dist.d[h.dist.d$value > 10, "value"] = 10
h.dist.b[h.dist.b$value > 10, "value"] = 10

h1 = hist(h.dist.b$value, breaks = seq(0, 10, 1)) # Change histogram to density
h1$density = h1$counts/sum(h1$counts)*100
h2 = hist(h.dist.d$value, breaks = seq(0, 10, 1)) # Change histogram to density
h2$density = h2$counts/sum(h2$counts)*100
h3 = hist(h.dist.s$value, breaks = seq(0, 10, 1)) # Change histogram to density
h3$density = h3$counts/sum(h3$counts)*100

plot(h1,freq = FALSE,
     col = adjustcolor(brewer.pal(n = 3, name = "Blues")[1], 1),
     xlab = "Path distance", ylab = "Proportion", xaxt = "n",
     ylim = c(0, max(c(h1$density, h2$density, h3$density))),
     main = "Parth distance with annotations (main component)")
plot(h2, freq = FALSE,
     col = adjustcolor(brewer.pal(n = 3, name = "Blues")[2], 0.5),lwd = 3, add = T)
plot(h3, freq = FALSE,
     col = adjustcolor(brewer.pal(n = 3, name = "Blues")[3], 0.5),lwd = 3, add = T)
axis(side = 1, at = 1:max(h1$breaks) - 0.5, labels=c(1:9,">10"))
legend("topright", legend=c(paste0("Not annotated (", length(h.dist.b$value), ")"), 
                            paste0("Between annotations (", length(h.dist.d$value), ")"),
                            paste0("Within annotation (", length(h.dist.s$value), ")")),
       fill=c(brewer.pal(n = 3, name = "Blues")), cex=0.8, box.lty=0)
rm(h1, h2, h3, h.dist.b, h.dist.d, h.dist.s)


### Direct neighbors - Figure 2e
dn.p = data.frame() # Fraction of (annotated) direct neighbors that are in the same complex
for (i in unique(Annot[Annot$Complex1 != "", "Gene1.order"]) ) { 
  tmp = dist[dist$Var1 == i | dist$Var2  == i,] # Each protein
  tmp = tmp[tmp$value == 1, ] # Direct neighbors
  tmp = tmp[tmp$Complex1.x != "" & tmp$Complex1.y != "", ] # Only annotated neighbors
  a = nrow(tmp[tmp$Complex1.x == tmp$Complex1.y,  ]) / nrow(tmp)
  dn.p = rbind(dn.p, c(i, a))
}
colnames(dn.p) = c("Complex", "Fraction.Dnei")
dn.p$Fraction.Dnei = as.numeric(dn.p$Fraction.Dnei)

# Randomized network control
tg.main.rewired <- rewire(tg.main, with = keeping_degseq(niter = vcount(tg.main) * 100, loops = FALSE))
tg.main.rewired = simplify(tg.main.rewired, remove.multiple = T, remove.loops = T)
tg.main.rewired = delete.vertices(tg.main.rewired, which(degree(tg.main.rewired)<1)) # remove nodes that have no edges

vert_ids = V(tg.main.rewired)[components(tg.main.rewired)$membership == 1] # Select only main component of the network
tg.main.rewired = induced_subgraph(tg.main.rewired, vert_ids)

# Distances between annotations control
dist.r = distances(tg.main.rewired, weights = NA , mode = "all")
dist.r[upper.tri(dist.r, diag = FALSE)] = 0 # Get only half of the corr matrix
dist.r = melt(dist.r)
dist.r = dist.r[dist.r$value != 0,]
dist.r = left_join(dist.r, Annot, by = c("Var1" = "Gene1.order"))
dist.r = left_join(dist.r, Annot, by = c("Var2" = "Gene1.order"))

# Direct neighbors
dn.p.r = data.frame() # Fraction of (annotated) direct neighbors that are in the same complex
for (i in unique(Annot[Annot$Complex1 != "", "Gene1.order"]) ) { 
  tmp = dist.r[dist.r$Var1 == i | dist.r$Var2  == i,] # Each protein
  tmp = tmp[tmp$value == 1, ] # Direct neighbors
  tmp = tmp[tmp$Complex1.x != "" & tmp$Complex1.y != "", ] # Only annotated neighbors
  a = nrow(tmp[tmp$Complex1.x == tmp$Complex1.y,  ]) / nrow(tmp)
  dn.p.r = rbind(dn.p.r, c(i, a))
}
colnames(dn.p.r) = c("Complex", "Fraction.Dnei")
dn.p.r$Fraction.Dnei = as.numeric(dn.p.r$Fraction.Dnei)

#Plot
h3 = hist(dn.p$Fraction.Dnei, breaks = seq(0, 1, 0.1)) # Change histogram to density
h3$density = h3$counts/sum(h3$counts)*100
h3.r = hist(dn.p.r$Fraction.Dnei, breaks = seq(0, 1, 0.1)) # Change histogram to density
h3.r$density = h3.r$counts/sum(h3.r$counts)*100

plot(h3.r,freq = FALSE,
     col = adjustcolor(brewer.pal(n = 3, name = "Blues")[1], 1),
     xlab = "Fraction direct neighbors", ylab = "Density", xaxt = "n",
     ylim = c(0, max(c(h3$density, h3.r$density))),
     main = "Fraction of (annotated) direct neighbors that are in the same complex (main component)")
plot(h3, freq = FALSE,
     col = adjustcolor(brewer.pal(n = 3, name = "Blues")[3], 0.5),lwd = 3, add = T)
axis(side = 1, at = seq(0.1, 1, 0.1) - 0.05 , labels=seq(0.1, 1, 0.1))
legend("topright",legend=c(paste0("Randomized network (", sum(!is.na(dn.p.r$Fraction.Dnei)), ")"), 
                           paste0("Endosomal network (", sum(!is.na(dn.p$Fraction.Dnei)), ")")),
       fill=c(brewer.pal(n = 3, name = "Blues"))[c(1,3)], cex=0.8, box.lty=0)
rm(h3,h3.r)

### Network and communities
# Cluster and color
C = cluster_edge_betweenness(tg.main, weights = NULL)# Choose clustering
mod = round(modularity(C),2)
modularity(tg.main, C$membership)
col = rainbow(max(C$membership))
V(tg.main)$color = col[C$membership]

# Add weight based on clusters
weight.community=function(row, membership, weigth.within, weight.between){
  if(as.numeric(membership[which(names(membership)==row[1])])==as.numeric(membership[which(names(membership)==row[2])])){
    weight=weigth.within
  }else{
    weight=weight.between
  }
  return(weight)
}

E(tg.main)$weight=apply(get.edgelist(tg.main), 1, weight.community, membership(C), 100,1) #Adjust distance between and within clusters

set.seed(10)
l = layout_with_fr(tg.main, weights=E(tg.main)$weight)

par(mar=c(0,0,0,0))
plot(tg.main, vertex.frame.color = NA, vertex.label.cex=0.1, 
     vertex.label.color='black', vertex.size = 1,
     vertex.label.dist = 0, edge.color = adjustcolor("grey40", 0.5), edge.width = 0.3, layout = l, rescale=T)
title(paste0("Layout cluster weighted fr modularity:", mod ,
             " Nodes = ",vcount(tg.main)," Edges = ", ecount(tg.main)), cex.main=1, line=-1, xpd=T)
dev.off()
