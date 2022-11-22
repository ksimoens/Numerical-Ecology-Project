library(tidyverse)

# load the necessary functions
source('functions.R')

# plot the results of the PCoA
# input = plot data.frame (see functions.R)
#         eigenvalues Evalues
#         fileName to be given to output png
#         option: =3 -> coloured by country
#                 =4 -> coloured by region
plotPCoA <- function(df.plot, Evalues, fileName, option){
  # proportion of variation explained (for axis labels)
  lambda <- round(Evalues[1:2] / sum(Evalues) *100, 2)  

  p <- df.plot %>%
  # plot site scores
  ggplot() + geom_point(aes(x=PCo1,y=PCo2,col=df.plot[,option],shape=df.plot[,option]),size=1.5) +
             theme_bw() +
             # colour scheme (see functions.R) 
             scale_colour_manual(
             	name = names(df.plot)[option],
             	labels = sort(unique(df.plot[,option])),
             	values=getColours(df.plot[,option])$colour) +
             scale_shape_manual(
             	name = names(df.plot)[option],
             	labels =sort(unique(df.plot[,option])),
             	values=getColours(df.plot[,option])$shape) +
             xlab(paste0("PCo1 (",lambda[1]," %)")) + ylab(paste0("PCo2 (",lambda[2]," %)")) +
             theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
             geom_hline(yintercept=0,linetype="dashed") + geom_vline(xintercept=0,linetype="dashed")

  p %>% ggsave(fileName,.,device='png',width=15,height=10,units='cm')
}

# plot PCoA for all sites and loci
Fst.mat <- read.csv("Output/distancesFst_tot.csv",header=T,row.names=1) %>%
           as.matrix()

result.tot <- PCoA(Fst.mat)
Evalues.tot <- result.tot$values
Evectors.tot <- result.tot$vectors

df.plot.tot <- makePlotDF(Evectors.tot,rownames(Fst.mat))

plotPCoA(df.plot.tot,Evalues.tot,"PCoA.tot.region.png",4)
plotPCoA(df.plot.tot,Evalues.tot,"PCoA.tot.country.png",3)
# write PCo axes to csv file for use in RDA
write.csv(Evectors.tot,"Output/PCo/PCo.tot.csv")

# plot PCoA for Atlantic sites and loci
Fst.mat.atl <- Fst.mat[-c(32:38),-c(32:38)]

result.atl <- PCoA(Fst.mat.atl)
Evalues.atl <- result.atl$values
Evectors.atl <- result.atl$vectors
# manually reverse PC1 (for consistency)
Evectors.atl[,1] <- -Evectors.atl[,1]

df.plot.atl <- makePlotDF(Evectors.atl,rownames(Fst.mat.atl))

plotPCoA(df.plot.atl,Evalues.atl,"PCoA.atl.region.png",4)
plotPCoA(df.plot.atl,Evalues.atl,"PCoA.atl.country.png",3)
write.csv(Evectors.atl,"Output/PCo/PCo.atl.csv")

# plot PCoA for all sites and neutral loci
Fst.neut <- read.csv("Output/distancesFst_neut.csv",header=T,row.names=1) %>%
           as.matrix()

result.neut.tot <- PCoA(Fst.neut)
Evalues.neut.tot <- result.neut.tot$values
Evectors.neut.tot <- result.neut.tot$vectors
# manually reverse PC2
Evectors.neut.tot[,2] <- -Evectors.neut.tot[,2]

df.plot.neut.tot <- makePlotDF(Evectors.neut.tot,rownames(Fst.neut))

plotPCoA(df.plot.neut.tot,Evalues.neut.tot,"PCoA.neut.tot.region.png",4)
plotPCoA(df.plot.neut.tot,Evalues.neut.tot,"PCoA.neut.tot.country.png",3)
write.csv(Evectors.neut.tot,"Output/PCo/PCo.neut.tot.csv")

# plot PCoA for Atlantic sites and neutral loci
Fst.neut <- as.data.frame(Fst.neut)
Fst.neut.atl <- Fst.neut[!(row.names(Fst.neut) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]
Fst.neut.atl <- as.matrix(Fst.neut.atl %>% dplyr::select(-Laz,-Tar,-Sar,-Ale,-The,-Tor,-Sky))

result.neut.atl <- PCoA(Fst.neut.atl)
Evalues.neut.atl <- result.neut.atl$values
Evectors.neut.atl <- result.neut.atl$vectors
Evectors.neut.atl[,2] <- -Evectors.neut.atl[,2]

df.plot.neut.atl <- makePlotDF(Evectors.neut.atl,rownames(Fst.neut.atl))

plotPCoA(df.plot.neut.atl,Evalues.neut.atl,"PCoA.neut.atl.region.png",4)
plotPCoA(df.plot.neut.atl,Evalues.neut.atl,"PCoA.neut.atl.country.png",3)
write.csv(Evectors.neut.atl,"Output/PCo/PCo.neut.atl.csv")

# plot PCoA for all sites and outlier loci
Fst.sel <- read.csv("Output/distancesFst_sel.csv",header=T,row.names=1) %>%
           as.matrix()

result.sel.tot <- PCoA(Fst.sel)
Evalues.sel.tot <- result.sel.tot$values
Evectors.sel.tot <- result.sel.tot$vectors
# manually reverse PC2
Evectors.sel.tot[,2] <- -Evectors.sel.tot[,2]

df.plot.sel.tot <- makePlotDF(Evectors.sel.tot,rownames(Fst.sel))

plotPCoA(df.plot.sel.tot,Evalues.sel.tot,"PCoA.sel.tot.region.png",4)
plotPCoA(df.plot.sel.tot,Evalues.sel.tot,"PCoA.sel.tot.country.png",3)
write.csv(Evectors.sel.tot,"Output/PCo/PCo.sel.tot.csv")

# plot PCoA for Atlantic sites and outlier loci
Fst.sel <- as.data.frame(Fst.sel)
Fst.sel.atl <- Fst.sel[!(row.names(Fst.sel) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]
Fst.sel.atl <- as.matrix(Fst.sel.atl %>% dplyr::select(-Laz,-Tar,-Sar,-Ale,-The,-Tor,-Sky))

result.sel.atl <- PCoA(Fst.sel.atl)
Evalues.sel.atl <- result.sel.atl$values
Evectors.sel.atl <- result.sel.atl$vectors
# manually reverse PC1 and PC2
Evectors.sel.atl[,1:2] <- -Evectors.sel.atl[,1:2]

df.plot.sel.atl <- makePlotDF(Evectors.sel.atl,rownames(Fst.sel.atl))

plotPCoA(df.plot.sel.atl,Evalues.sel.atl,"PCoA.sel.atl.region.png",4)
plotPCoA(df.plot.sel.atl,Evalues.sel.atl,"PCoA.sel.atl.country.png",3)
write.csv(Evectors.sel.atl,"Output/PCo/PCo.sel.atl.csv")