library(tidyverse)

PCoA <- function(Fst){
	A <- -0.5 * Fst * Fst

	identity <- diag(ncol(A))
	ones <- rep(1,ncol(A))
	side <- identity - ones %*% t(ones) / ncol(A)

	Delta <- side %*% A %*% side

	Evalues <- eigen(Delta)$values
	print(Evalues)

	c <- abs(Evalues[length(Evalues)])
	print(c)
	c.mat <- matrix(rep(c,ncol(A)*ncol(A)),ncol=ncol(A))
	c.mat <- c.mat - c*diag(ncol(A))

	Fst.corr <- sqrt(Fst * Fst + 2*c.mat)
	A <- -0.5 * Fst.corr * Fst.corr
	Delta <- side %*% A %*% side

	Evalues <- eigen(Delta)$values
	Evectors <- eigen(Delta)$vectors
	names <- paste0("PCo",1:ncol(Evectors))
	colnames(Evectors) <- names
	rownames(Evectors) <- rownames(Fst)

  Evalues <- Evalues[which(Evalues > 1e-10)]
  Evectors <- Evectors[,which(Evalues > 1e-10)]

  for(i in 1:ncol(Evectors)){
    Evectors[,i] <- Evectors[,i]*Evalues[i]
  }

	return(list(values=Evalues,vectors=Evectors))
}

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

makePlotDF <- function(Evalues,Evectors,names){
  df.plot <- Evectors[,1:2] %>% as.data.frame()
  rownames(df.plot) <- names
  df.plot$country <- rownames(df.plot)
  df.plot$country <- countryList(df.plot$country)
  df.plot$region <- rownames(df.plot)
  df.plot$region <- regionList(df.plot$region)

  print(df.plot)
  return(df.plot)
  
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

plotPCoA <- function(df.plot, Evalues, fileName, option){
  lambda <- round(Evalues[1:2] / sum(Evalues) *100, 2)  

  p <- df.plot %>%
  ggplot() + geom_point(aes(x=PCo1,y=PCo2,col=df.plot[,option],shape=df.plot[,option]),size=1.5) +
             theme_bw() + 
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

Fst.mat <- read.csv("Output/distancesFst_tot.csv",header=T,row.names=1) %>%
           as.matrix()

result.tot <- PCoA(Fst.mat)
Evalues.tot <- result.tot$values
Evectors.tot <- result.tot$vectors

df.plot.tot <- makePlotDF(Evalues.tot,Evectors.tot,rownames(Fst.mat))

plotPCoA(df.plot.tot,Evalues.tot,"PCoA.tot.region.png",4)
plotPCoA(df.plot.tot,Evalues.tot,"PCoA.tot.country.png",3)
write.csv(Evectors.tot,"Output/PCo/PCo.tot.csv")

Fst.mat.atl <- Fst.mat[-c(32:38),-c(32:38)]

result.atl <- PCoA(Fst.mat.atl)
Evalues.atl <- result.atl$values
Evectors.atl <- result.atl$vectors
Evectors.atl[,1] <- -Evectors.atl[,1]

df.plot.atl <- makePlotDF(Evalues.atl,Evectors.atl,rownames(Fst.mat.atl))

plotPCoA(df.plot.atl,Evalues.atl,"PCoA.atl.region.png",4)
plotPCoA(df.plot.atl,Evalues.atl,"PCoA.atl.country.png",3)
write.csv(Evectors.atl,"Output/PCo/PCo.atl.csv")

Fst.neut <- read.csv("Output/distancesFst_neut.csv",header=T,row.names=1) %>%
           as.matrix()

result.neut.tot <- PCoA(Fst.neut)
Evalues.neut.tot <- result.neut.tot$values
Evectors.neut.tot <- result.neut.tot$vectors
Evectors.neut.tot[,2] <- -Evectors.neut.tot[,2]

df.plot.neut.tot <- makePlotDF(Evalues.neut.tot,Evectors.neut.tot,rownames(Fst.neut))

plotPCoA(df.plot.neut.tot,Evalues.neut.tot,"PCoA.neut.tot.region.png",4)
plotPCoA(df.plot.neut.tot,Evalues.neut.tot,"PCoA.neut.tot.country.png",3)
write.csv(Evectors.neut.tot,"Output/PCo/PCo.neut.tot.csv")

Fst.neut <- as.data.frame(Fst.neut)
Fst.neut.atl <- Fst.neut[!(row.names(Fst.neut) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]
Fst.neut.atl <- as.matrix(Fst.neut.atl %>% dplyr::select(-Laz,-Tar,-Sar,-Ale,-The,-Tor,-Sky))

result.neut.atl <- PCoA(Fst.neut.atl)
Evalues.neut.atl <- result.neut.atl$values
Evectors.neut.atl <- result.neut.atl$vectors
Evectors.neut.atl[,2] <- -Evectors.neut.atl[,2]

df.plot.neut.atl <- makePlotDF(Evalues.neut.atl,Evectors.neut.atl,rownames(Fst.neut.atl))

plotPCoA(df.plot.neut.atl,Evalues.neut.atl,"PCoA.neut.atl.region.png",4)
plotPCoA(df.plot.neut.atl,Evalues.neut.atl,"PCoA.neut.atl.country.png",3)
write.csv(Evectors.neut.atl,"Output/PCo/PCo.neut.atl.csv")

Fst.sel <- read.csv("Output/distancesFst_sel.csv",header=T,row.names=1) %>%
           as.matrix()

result.sel.tot <- PCoA(Fst.sel)
Evalues.sel.tot <- result.sel.tot$values
Evectors.sel.tot <- result.sel.tot$vectors

df.plot.sel.tot <- makePlotDF(Evalues.sel.tot,Evectors.sel.tot,rownames(Fst.sel))

plotPCoA(df.plot.sel.tot,Evalues.sel.tot,"PCoA.sel.tot.region.png",4)
plotPCoA(df.plot.sel.tot,Evalues.sel.tot,"PCoA.sel.tot.country.png",3)
write.csv(Evectors.sel.tot,"Output/PCo/PCo.sel.tot.csv")

Fst.sel <- as.data.frame(Fst.sel)
Fst.sel.atl <- Fst.sel[!(row.names(Fst.sel) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]
Fst.sel.atl <- as.matrix(Fst.sel.atl %>% dplyr::select(-Laz,-Tar,-Sar,-Ale,-The,-Tor,-Sky))

result.sel.atl <- PCoA(Fst.sel.atl)
Evalues.sel.atl <- result.sel.atl$values
Evectors.sel.atl <- result.sel.atl$vectors
Evectors.sel.atl[,1:2] <- -Evectors.sel.atl[,1:2]

df.plot.sel.atl <- makePlotDF(Evalues.sel.atl,Evectors.sel.atl,rownames(Fst.sel.atl))

plotPCoA(df.plot.sel.atl,Evalues.sel.atl,"PCoA.sel.atl.region.png",4)
plotPCoA(df.plot.sel.atl,Evalues.sel.atl,"PCoA.sel.atl.country.png",3)
write.csv(Evectors.sel.atl,"Output/PCo/PCo.sel.atl.csv")
