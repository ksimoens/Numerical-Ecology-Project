library(vegan)
library(adespatial)
library(tidyverse)
library(ggforce)
library(ggnewscale)
library(rdacca.hp)

Fst.atl <- read.csv("Output/PCo/PCo.neut.atl.csv",header=T,row.names=1)

linear <- read.csv("Output/PCoSpatial.csv",header=T,row.names=1)[,1:2]

dbMEM <- read.csv("Output/dbMEM.csv",header=T,row.names=1)

AEM <- read.csv("Output/AEM.csv",header=T,row.names=1)

Fst.atl <- Fst.atl[match(rownames(AEM),rownames(Fst.atl)),]
linear <- linear[match(rownames(AEM),rownames(linear)),]
dbMEM <- dbMEM[match(rownames(AEM),rownames(dbMEM)),]

rda.linear <- rda(Fst.atl,linear)
R2adj.linear <- RsquareAdj(rda.linear)
print(anova(rda.linear,permutations = how(nperm=9999)))

rda.dbMEM <- rda(Fst.atl, dbMEM)
R2adj.dbMEM <- RsquareAdj(rda.dbMEM)
print(anova(rda.dbMEM,permutations = how(nperm=9999)))
sel.dbMEM <- forward.sel(Fst.atl, dbMEM, adjR2thresh=R2adj.dbMEM, alpha=0.05, nperm=9999)
print(sel.dbMEM)

dbMEM.sel <- dbMEM[,sel.dbMEM$order]

rda.AEM <- rda(Fst.atl, AEM)
R2adj.AEM <- RsquareAdj(rda.AEM)
print(anova(rda.AEM,permutations = how(nperm=9999)))
sel.AEM <- forward.sel(Fst.atl, AEM, adjR2thresh=R2adj.AEM, alpha=0.05, nperm=9999)
print(sel.AEM)

AEM.sel <- AEM[,sel.AEM$order]

fractions <- vector()
R2adj <- vector()
degrees <- vector()
Fvalues <- vector()
pvalues <- vector()

rda.part <- varpart(Fst.atl, linear, dbMEM.sel, AEM.sel)
rda.part.df <- as.data.frame(rbind(rda.part$part$fract,rda.part$part$indfract,rda.part$part$contr1))
rda.part.df <- rda.part.df[rda.part.df$Testable==T,]
rownames(rda.part.df) <- c("[lin]","[MEM]","[AEM]","[lin] + [MEM]","[lin] + [AEM]","[MEM] + [AEM]",
	"[lin] + [MEM] + [AEM]","[lin] | [MEM] + [AEM]","[MEM] | [lin] + [AEM]","[AEM] | [lin] + [MEM]",
	"[lin] | [AEM]","[lin] | [MEM]","[MEM] | [AEM]","[MEM] | [lin]","[AEM] | [lin]","[AEM] | [MEM]")

key <- list("[lin]"=linear, "[MEM]"=dbMEM.sel, "[AEM]"=AEM.sel)

makeRDA <- function(rowname){
	str_res <- strsplit(rowname, " +")[[1]]
	print(str_res)
	str_res <- str_res[str_res!="+"]
	n <- length(str_res)
	index <- which(str_res=="|")
	if(length(index) == 0){
		index <- n+1
	} else {
		index <- index[1]
	}
	variable <- as.data.frame(key[names(key)==str_res[1]])
	if(index > 2){
		for(i in 2:(index-1)){
			var <- str_res[i]
			variable <- cbind(variable,key[names(key)==var])
		}
	}
	print(names(variable))
	if(index != n+1){
		conditional <- as.data.frame(key[names(key)==str_res[index+1]])
		if(index + 1 < n){
			for(i in (index+2):(n)){
				var <- str_res[i]
				conditional <- cbind(conditional,key[names(key)==var])
			}
		}
		rda.out <- rda(Fst.atl, variable, conditional)
		print(names(conditional))
	} else {
		rda.out <- rda(Fst.atl, variable)
	}
	return(rda.out)
}

for(i in 1:nrow(rda.part.df)){
	fractions <- c(fractions, rownames(rda.part.df)[i])
	R2adj <- c(R2adj, rda.part.df$Adj.R.square[i])

	rda <- makeRDA(rownames(rda.part.df)[i])
	print(summary(rda))
	an <- anova(rda, permutations = how(nperm=9999))
	deg <- paste0(an$Df[1],"(",an$Df[1]+an$Df[2],")")
	degrees <- c(degrees, deg)
	Fvalues <- c(Fvalues, an$F[1])
	pvalues <- c(pvalues, an$"Pr(>F)"[1])
}

df.parts <- data.frame(fraction=fractions, RsquareAdj=R2adj, dof=degrees, F=Fvalues, "p-value"=pvalues)
print(df.parts)

plotRDA <- function(df.plot, Evalues, VAR, scale, fileName, option, scaling){
  lambda <- sprintf("%.2f",round(Evalues[1:2] / sum(Evalues) *100, 2))
  rad <- scale*sqrt(2/length(Evalues))
  p <- ggplot() + geom_segment(data=VAR,aes(x=rep(0,nrow(VAR)),y=rep(0,nrow(VAR)),xend=VAR[,1],yend=VAR[,2],col=variable), 
             	lineend='round', arrow = arrow(length = unit(0.1, "inches"))) + 
             scale_color_viridis_d(option='rocket') +
             new_scale_colour() +
             geom_circle(aes(x0=0,y0=0,r=rad),linetype=scaling) +  			 
  			 geom_text(data=VAR,aes(x=1.12*RDA1,y=1.07*RDA2,label=rownames(VAR)),size=3) +
  			 geom_point(data=df.plot,aes(x=RDA1,y=RDA2,col=df.plot[,option]),size=1.5) +
             theme_bw() + scale_color_viridis_d(option='magma') +
             xlab(paste0("RDA1 (",lambda[1]," %)")) + ylab(paste0("RDA2 (",lambda[2]," %)")) +
             theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
             geom_hline(yintercept=0,linetype="dashed") + geom_vline(xintercept=0,linetype="dashed") +
             guides(col=guide_legend(title=colnames(df.plot)[option])) +
             xlim(1.15*min(df.plot$RDA1,-rad,VAR$RDA1),1.15*max(df.plot$RDA1,rad,VAR$RDA1)) +
             ylim(1.15*min(df.plot$RDA2,-rad,VAR$RDA2),1.15*max(df.plot$RDA2,rad,VAR$RDA2)) 

  p %>% ggsave(fileName,.,device='png',width=15,height=10,units='cm')
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

rda.final <- rda(Fst.atl, cbind(linear,dbMEM.sel,AEM.sel))
rda.sum <- summary(rda.final)
Evalues <- as.vector(rda.sum$cont$importance[2,])
site.constraints <- scores(rda.final, scaling=1, display=c("lc"))
site.constraints[,2] <- -site.constraints[,2]
variables <- 0.05*scores(rda.final, scaling=1, display=c("bp"))
variables[,2] <- -variables[,2]
scale.fac <- 0.05*attributes(variables)$const
variables <- as.data.frame(variables)
rownames(variables)[3] <- "MEM4"
variables$variable <- c("linear","linear","MEM",rep("AEM",9))

df.plot.s1 <- makePlotDF(site.constraints, rownames(site.constraints))
plotRDA(df.plot.s1,Evalues,variables,scale.fac,"RDA_Fst_neut_s1_country.png",3,1)

site.constraints.s2 <- scores(rda.final, scaling=2, display=c("lc"))
site.constraints.s2[,2] <- -site.constraints.s2[,2]
variables.s2 <- 0.05*scores(rda.final, scaling=2, display=c("bp"))
variables.s2[,2] <- -variables.s2[,2]
scale.fac.s2 <- 0.05*attributes(variables.s2)$const
variables.s2 <- as.data.frame(variables.s2)
rownames(variables.s2)[3] <- "MEM4"
variables.s2$variable <- c("linear","linear","MEM",rep("AEM",9))

df.plot.s2 <- makePlotDF(site.constraints.s2, rownames(site.constraints.s2))
plotRDA(df.plot.s2,Evalues,variables.s2,scale.fac.s2,"RDA_Fst_neut_s2_country.png",3,0)

site.unconstrained <- as.data.frame(scores(rda.final, scaling=1, choices=c(10,11), display=c("sites")))
site.unconstrained$country <- countryList(rownames(site.unconstrained))
Evalues <- as.vector(rda.sum$cont$importance[2,10:11])
lambda <- sprintf("%.2f",round(Evalues[1:2]*100, 2))

p <- ggplot() + geom_point(data=site.unconstrained,aes(x=PC1,y=PC2,col=country),size=1.5) +
             theme_bw() + scale_color_viridis_d(option='magma') +
             xlab(paste0("PC1 (",lambda[1]," %)")) + ylab(paste0("PC2 (",lambda[2]," %)")) +
             theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
             geom_hline(yintercept=0,linetype="dashed") + geom_vline(xintercept=0,linetype="dashed")              

p %>% ggsave("resid_PCA_Fst_neut.png",.,device='png',width=15,height=10,units='cm')

df.dbMEM <- as.data.frame(dbMEM.sel)
names(df.dbMEM) <- "MEM4"
rdacca <- rdacca.hp(Fst.atl,list(as.data.frame(linear), df.dbMEM, as.data.frame(AEM.sel)), method='RDA')
rdacca.sep <- rdacca.hp(Fst.atl,cbind(linear,df.dbMEM,AEM.sel), method='RDA')

permutest <- permu.hp(Fst.atl,list(as.data.frame(env.atl.norm.sel), as.data.frame(linear), as.data.frame(dbMEM.sel)), method='RDA', permutations=9999)
permutest.sep <- permu.hp(Fst.atl,cbind(env.atl.norm.sel,linear,dbMEM.sel), method='RDA', permutations=9999)

p <- plot.rdaccahp(rdacca.sep) + theme_bw()
p %>% ggsave("rdacca_sep_Fst_neut.png",.,device='png',width=15,height=7.5,units='cm')