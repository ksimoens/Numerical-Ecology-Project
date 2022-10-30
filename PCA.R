library(vegan)
library(tidyverse)
library(ggforce)

countryList <- function(populations){
  reference <- data.frame(
                  pop = c("Brd","Cro","Eye","Heb","Iom","Ios","Loo","Lyn","Ork","Pad","Pem","She","Sbs","Sul",
                          "Jer",
                          "Idr",
                          "Hel",
                          "Ale","Sky","The","Tor",
                          "Cor","Hoo","Kil","Mul","Ven",
                          "Laz","Tar","Sar",
                          "Oos",
                          "Ber","Flo","Sin","Tro",
                          "Vig",
                          "Gul","Kav","Lys"),
                  country = c(rep("GRB",14),"CHA","FRA","DEU",rep("HEL",4),rep("IRL",5),rep("ITA",3),"NDL",rep("NOR",4),"ESP",rep("SVE",3))
                )
  for(i in 1:nrow(reference)){
    populations[populations == reference$pop[i]] <- reference$country[i]
  }

  return(populations)
}

regionList <- function(populations){
  reference <- data.frame(
                  pop = c("Brd","Cro","Eye","Heb","Iom","Ios","Loo","Lyn","Ork","Pad","Pem","She","Sbs","Sul",
                          "Jer",
                          "Idr",
                          "Hel",
                          "Ale","Sky","The","Tor",
                          "Cor","Hoo","Kil","Mul","Ven",
                          "Laz","Tar","Sar",
                          "Oos",
                          "Ber","Tro",
                          "Flo","Sin",
                          "Vig",
                          "Gul","Kav","Lys"),
                  region = c(rep("ATL",17),rep("MED",4),rep("ATL",5),rep("MED",3),rep("ATL",3),rep("SKA",2),"ATL",rep("SKA",3))
                )
  for(i in 1:nrow(reference)){
    populations[populations == reference$pop[i]] <- reference$region[i]
  }

  return(populations)
}

getColours <- function(pop){
	col_key <- data.frame(
			population = c("CHA","DEU","ESP","FRA","GRB","HEL","IRL","ITA","NDL","NOR","SVE","ATL","MED","SKA"),
			colour = c(rep(c("#C7E020FF","#24868EFF","#440154FF"),3),c("#C7E020FF","#440154FF"),c("#C7E020FF","#24868EFF","#440154FF")),
			shape = c(rep(15,3),rep(17,3),rep(19,3),18,18,15,17,19))

	col_res <- col_key[col_key$population %in% unique(pop),]$colour
	sha_res <- col_key[col_key$population %in% unique(pop),]$shape
	return(list(colour = col_res,shape = sha_res))
}

makePlotDF <- function(Evectors,names){
  df.plot <- Evectors[,1:2] %>% as.data.frame()
  rownames(df.plot) <- names
  df.plot$country <- rownames(df.plot)
  df.plot$country <- countryList(df.plot$country)
  df.plot$region <- rownames(df.plot)
  df.plot$region <- regionList(df.plot$region)

  print(df.plot)
  return(df.plot)
}

plotPCA <- function(df.plot, Evalues, U, scale, fileName, option){
  lambda <- sprintf("%.2f",round(Evalues[1:2] / sum(Evalues) *100, 2))
  rad <- scale*sqrt(2/length(Evalues))
  p <- ggplot() + geom_segment(data=U,x=rep(0,nrow(U)),y=rep(0,nrow(U)),xend=U[,1],yend=U[,2], 
             	lineend='round', arrow = arrow(length = unit(0.1, "inches"))) + 
  			 geom_circle(aes(x0=0,y0=0,r=rad)) +
  			 geom_text(data=U,aes(x=1.12*PC1,y=1.07*PC2,label=rownames(U)),size=3) +
  			 geom_point(data=df.plot,aes(x=PC1,y=PC2,col=df.plot[,option],shape=df.plot[,option]),size=1.5) +
             theme_bw() + 
			 scale_colour_manual(
             	name = names(df.plot)[option],
             	labels = sort(unique(df.plot[,option])),
             	values=getColours(df.plot[,option])$colour) +
             scale_shape_manual(
             	name = names(df.plot)[option],
             	labels =sort(unique(df.plot[,option])),
             	values=getColours(df.plot[,option])$shape) +
             xlab(paste0("PC1 (",lambda[1]," %)")) + ylab(paste0("PC2 (",lambda[2]," %)")) +
             theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
             geom_hline(yintercept=0,linetype="dashed") + geom_vline(xintercept=0,linetype="dashed") +
             guides(col=guide_legend(title=colnames(df.plot)[option])) +
             xlim(1.15*min(df.plot$PC1,U$PC1,-rad),1.15*max(df.plot$PC1,U$PC1,rad)) +
             ylim(1.15*min(df.plot$PC2,U$PC2,-rad),1.15*max(df.plot$PC2,U$PC2,rad)) 

  p %>% ggsave(fileName,.,device='png',width=15,height=10,units='cm')
}

PCAOutput <- function(allele.freq,thresh){
	PCA <- rda(allele.freq)
	PCA_sum <- summary(PCA)
	PCA.good <- as.data.frame(goodness(PCA, model = "CA", choices = 1:3))
	print(PCA.good[order(-PCA.good$PC2),])

	SNP_list <- rownames(PCA.good[PCA.good$PC2 > thresh,])
	site.scores <- scores(PCA, scaling=1, display="sites")
	site.scores[,1] <- -site.scores[,1]
	species.scores <- 0.5*scores(PCA, scaling=1, display="species")
	scale <- 0.5*attributes(species.scores)$const
	species.scores <- as.data.frame(species.scores[rownames(species.scores) %in% SNP_list,])
	rownames(species.scores) <- substr(rownames(species.scores),start=2,stop=nchar(rownames(species.scores))-2)
	species.scores$PC1 <- -species.scores$PC1

	Evalues <- as.vector(PCA_sum$cont$importance[2,])

	PCAOut <- list(Evectors=site.scores, Umatrix=species.scores, Evalues=Evalues, Scale=scale)

}

al.freq.tot <- read.csv("Output/allele_freq_total.csv", header=T,row.names=1)
PCA.tot <- PCAOutput(al.freq.tot,0.90)
df.plot.tot <- makePlotDF(PCA.tot$Evectors, rownames(PCA.tot$Evectors))
plotPCA(df.plot.tot,PCA.tot$Evalues,PCA.tot$Umatrix,PCA.tot$Scale,"PCA_tot_country.png",3)
plotPCA(df.plot.tot,PCA.tot$Evalues,PCA.tot$Umatrix,PCA.tot$Scale,"PCA_tot_region.png",4)

al.freq.atl <- al.freq.tot[!(row.names(al.freq.tot) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]
PCA.atl <- PCAOutput(al.freq.atl,0.7)
PCA.atl$Evectors[,2] <- -PCA.atl$Evectors[,2]
PCA.atl$Umatrix[,2] <- -PCA.atl$Umatrix[,2]
df.plot.atl <- makePlotDF(PCA.atl$Evectors, rownames(PCA.atl$Evectors))
plotPCA(df.plot.atl,PCA.atl$Evalues,PCA.atl$Umatrix,PCA.atl$Scale,"PCA_atl_country.png",3)
plotPCA(df.plot.atl,PCA.atl$Evalues,PCA.atl$Umatrix,PCA.atl$Scale,"PCA_atl_region.png",4)

al.freq.sel.tot <- read.csv("Output/allele_freq_sel.csv", header=T, row.names=1)
PCA.sel.tot <- PCAOutput(al.freq.sel.tot,0.01)
df.plot.sel.tot <- makePlotDF(PCA.sel.tot$Evectors, rownames(PCA.sel.tot$Evectors))
plotPCA(df.plot.sel.tot,PCA.sel.tot$Evalues,PCA.sel.tot$Umatrix,PCA.sel.tot$Scale,"PCA_sel_tot_country.png",3)
plotPCA(df.plot.sel.tot,PCA.sel.tot$Evalues,PCA.sel.tot$Umatrix,PCA.sel.tot$Scale,"PCA_sel_tot_region.png",4)

al.freq.sel.atl <- al.freq.sel.tot[!(row.names(al.freq.sel.tot) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]
PCA.sel.atl <- PCAOutput(al.freq.sel.atl,0.01)
df.plot.sel.atl <- makePlotDF(PCA.sel.atl$Evectors, rownames(PCA.sel.atl$Evectors))
plotPCA(df.plot.sel.atl,PCA.sel.atl$Evalues,PCA.sel.atl$Umatrix,PCA.sel.atl$Scale,"PCA_sel_atl_country.png",3)
plotPCA(df.plot.sel.atl,PCA.sel.atl$Evalues,PCA.sel.atl$Umatrix,PCA.sel.atl$Scale,"PCA_sel_atl_region.png",4)

al.freq.neut.tot <- read.csv("Output/allele_freq_neut.csv", header=T, row.names=1)
PCA.neut.tot <- PCAOutput(al.freq.neut.tot,0.80)
df.plot.neut.tot <- makePlotDF(PCA.neut.tot$Evectors, rownames(PCA.neut.tot$Evectors))
plotPCA(df.plot.neut.tot,PCA.neut.tot$Evalues,PCA.neut.tot$Umatrix,PCA.neut.tot$Scale,"PCA_neut_tot_country.png",3)
plotPCA(df.plot.neut.tot,PCA.neut.tot$Evalues,PCA.neut.tot$Umatrix,PCA.neut.tot$Scale,"PCA_neut_tot_region.png",4)

al.freq.neut.atl <- al.freq.neut.tot[!(row.names(al.freq.neut.tot) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]
PCA.neut.atl <- PCAOutput(al.freq.neut.atl,0.70)
df.plot.neut.atl <- makePlotDF(PCA.neut.atl$Evectors, rownames(PCA.neut.atl$Evectors))
plotPCA(df.plot.neut.atl,PCA.neut.atl$Evalues,PCA.neut.atl$Umatrix,PCA.neut.atl$Scale,"PCA_neut_atl_country.png",3)
plotPCA(df.plot.neut.atl,PCA.neut.atl$Evalues,PCA.neut.atl$Umatrix,PCA.neut.atl$Scale,"PCA_neut_atl_region.png",4)
