## Title: EndoMAP
## Author: Miguel A. Gonzalez-Lozano
## Date: 08/14/2024
  
library(dplyr)
library(tidyr)


###
### PCprophet results - Protein complex to protein interactions
###

# Import PCprophet output files
file_Expression = ".../rf.txt"
t = read.delim(file_Expression, header=T, stringsAsFactors=F, skipNul = TRUE)

file_Expression = ".../cmplx_combined.txt"
d = read.delim(file_Expression, header=T, stringsAsFactors=F, skipNul = TRUE) 

file_Expression = ".../peak_list.txt"
e = read.delim(file_Expression, header=T, stringsAsFactors=F, skipNul = TRUE) 

ppi = data.frame()
for (i in unique(t$ID)) {
  tryCatch({ #Skip and print errors
    
    x = t[t$ID ==i, "POS"] # Get score
    y = strsplit(d[d$ID ==i, "MB"], split = "#")[[1]] # Get members and make pairs
    y = sort(y) # Sort alphabetical (relevant later for selecting ppis)
    y.comb = t(combn(y,2)) # ppis combinations
    
    z.tmp = e[e$ID == i,] # Get peaks for complex
    z = y.comb
    for(id in 1:nrow(z.tmp)){ # Get peaks for proteins in ppi
      z[z %in% z.tmp$MB[id]] = z.tmp$SEL[id] #Change proteins for peak
    }
    
    z = paste(z[,1],z[,2], sep = ";") 
    
    x.y.z = cbind(y.comb,x,length(y),z)
    
    ppi = rbind(ppi,x.y.z) # Paste in df
    
  }, error=function(e){cat("ERROR :", i, conditionMessage(e), "\n")})
}
colnames(ppi) = c("Gene1","Gene2","POS","Complex_Size", "Peak")
ppi$POS = as.numeric(ppi$POS)
ppi$Complex_Size = as.numeric(ppi$Complex_Size)


ppi$PPI = paste(ppi$Gene1,ppi$Gene2,sep="-")

# Filter by fraction
Max.fr = 63 # BSA peak
ppi = separate(data = ppi, col = Peak, into = c("Peak1", "Peak2"), sep = ";")
ppi = filter(ppi, Peak1<Max.fr, Peak2<Max.fr)

# Filter by complex members
Max.size = 25
ppi = filter(ppi, Complex_Size<Max.size)

# Filter repeated ppis: get best score
ppi_best = data.frame()
for (j in unique(ppi$PPI)){
  
  tmp = ppi[ppi$PPI ==j,]  
  
  ppi_best = rbind(ppi_best, cbind(j, max(tmp$POS)))
}
colnames(ppi_best) = c("PPI","Best_Score")

write.table(ppi_best,
            sprintf(paste(file_Expression2,"_ppis_best_filter.tsv",sep = "")), sep="\t", col.names=T, row.names=F)
