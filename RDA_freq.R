library(vegan)
library(adespatial)
library(tidyverse)
library(ggforce)

freq <- read.csv("Output/allele_freq_sel.csv",header=T,row.names=1)
freq.atl <- freq[!(row.names(freq) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]

env <- read.table("Output/EnvMatrix.txt",header=T,row.names=1)
env.atl <- env[!(row.names(env) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]

linear <- read.csv("Output/PCoSpatial.csv",header=T,row.names=1)[,1:2]

dbMEM <- read.csv("Output/dbMEM.csv",header=T,row.names=1)

AEM <- read.csv("Output/AEM.csv",header=T,row.names=1)

freq.atl <- freq.atl[match(rownames(AEM),rownames(freq.atl)),]
env.atl <- env.atl[match(rownames(AEM),rownames(env.atl)),]
linear <- linear[match(rownames(AEM),rownames(linear)),]
dbMEM <- dbMEM[match(rownames(AEM),rownames(dbMEM)),]

env.atl.norm <- scale(env.atl, center=T, scale=T)

rda.env <- rda(freq.atl, env.atl.norm)
R2adj.env <- RsquareAdj(rda.env)
print(anova(rda.env,permutations = how(nperm=9999)))
sel.env <- forward.sel(freq.atl, env.atl.norm, adjR2thresh=R2adj.env, alpha=0.05, nperm=9999)
print(sel.env)

env.atl.norm.sel <- env.atl.norm[,sel.env$order]

rda.linear <- rda(freq.atl,linear)
R2adj.linear <- RsquareAdj(rda.linear)
print(anova(rda.linear,permutations = how(nperm=9999)))

rda.dbMEM <- rda(freq.atl, dbMEM)
R2adj.dbMEM <- RsquareAdj(rda.dbMEM)
print(anova(rda.dbMEM,permutations = how(nperm=9999)))
sel.dbMEM <- forward.sel(freq.atl, dbMEM, adjR2thresh=R2adj.dbMEM, alpha=0.05, nperm=9999)
print(sel.dbMEM)

dbMEM.sel <- dbMEM[,sel.dbMEM$order]

rda.AEM <- rda(freq.atl, AEM)
R2adj.AEM <- RsquareAdj(rda.AEM)
print(anova(rda.AEM,permutations = how(nperm=9999)))
sel.AEM <- forward.sel(freq.atl, AEM, adjR2thresh=R2adj.AEM, alpha=0.05, nperm=9999)
print(sel.AEM)

fractions <- vector()
R2adj <- vector()
degrees <- vector()
Fvalues <- vector()
pvalues <- vector()

rda.part <- varpart(freq.atl, env.atl.norm.sel, linear, dbMEM.sel)
rda.part.df <- as.data.frame(rbind(rda.part$part$fract,rda.part$part$indfract,rda.part$part$contr1))
rda.part.df <- rda.part.df[rda.part.df$Testable==T,]
rownames(rda.part.df) <- c("[env]","[lin]","[MEM]","[env] + [lin]","[env] + [MEM]","[lin] + [MEM]",
	"[env] + [lin] + [MEM]","[env] | [lin] + [MEM]","[lin] | [env] + [MEM]","[MEM] | [env] + [lin]",
	"[env] | [MEM]","[env] | [lin]","[lin] | [MEM]","[lin] | [env]","[MEM] | [env]","[MEM] | [lin]")

key <- list("[env]"=env.atl.norm.sel, "[lin]"=linear, "[MEM]"=dbMEM.sel)

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
		rda.out <- rda(freq.atl, variable, conditional)
		print(names(conditional))
	} else {
		rda.out <- rda(freq.atl, variable)
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

plotRDA <- function(df.plot, Evalues, U, VAR, scale, fileName, option, scaling){
  lambda <- sprintf("%.2f",round(Evalues[1:2] / sum(Evalues) *100, 2))
  rad <- scale*sqrt(2/length(Evalues))
  p <- ggplot() + geom_segment(data=U,x=rep(0,nrow(U)),y=rep(0,nrow(U)),xend=U[,1],yend=U[,2], 
             	lineend='round', arrow = arrow(length = unit(0.1, "inches"))) +
             geom_segment(data=VAR,aes(x=rep(0,nrow(VAR)),y=rep(0,nrow(VAR)),xend=VAR[,1],yend=VAR[,2],col=variable), 
             	linetype='dashed',lineend='round', arrow = arrow(length = unit(0.1, "inches"))) + 
             scale_color_viridis_d(option='rocket') +
             new_scale_colour() +
             geom_circle(aes(x0=0,y0=0,r=rad),linetype=scaling) +  			 
  			 geom_text(data=U,aes(x=1.12*RDA1,y=1.07*RDA2,label=rownames(U)),size=3) +
  			 geom_text(data=VAR,aes(x=1.12*RDA1,y=1.07*RDA2,label=rownames(VAR)),size=3) +
  			 geom_point(data=df.plot,aes(x=RDA1,y=RDA2,col=df.plot[,option]),size=1.5) +
             theme_bw() + scale_color_viridis_d(option='magma') +
             xlab(paste0("RDA1 (",lambda[1]," %)")) + ylab(paste0("RDA2 (",lambda[2]," %)")) +
             theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
             geom_hline(yintercept=0,linetype="dashed") + geom_vline(xintercept=0,linetype="dashed") +
             guides(col=guide_legend(title=colnames(df.plot)[option])) +
             xlim(1.15*min(df.plot$RDA1,U$RDA1,-rad,VAR$RDA1),1.15*max(df.plot$RDA1,U$RDA1,rad,VAR$RDA1)) +
             ylim(1.15*min(df.plot$RDA2,U$RDA2,-rad,VAR$RDA2),1.15*max(df.plot$RDA2,U$RDA2,rad,VAR$RDA2)) 

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

rda.final <- rda(freq.atl, cbind(env.atl.norm.sel,linear,dbMEM.sel))
rda.sum <- summary(rda.final)
Evalues <- as.vector(rda.sum$cont$importance[2,])
site.constraints <- scores(rda.final, scaling=1, display=c("lc"))
species.scores <- 0.4*scores(rda.final, scaling=1, display=c("species"))
rownames(species.scores) <- substr(rownames(species.scores),start=2,stop=nchar(rownames(species.scores))-2)
variables <- as.data.frame(0.4*scores(rda.final, scaling=1, display=c("bp")))
variables$variable <- c("environment","environment","linear","linear","MEM","MEM","MEM")
scale.fac <- 0.4*attributes(species.scores)$const
species.scores <- as.data.frame(species.scores)

df.plot.s1 <- makePlotDF(site.constraints, rownames(site.constraints))
plotRDA(df.plot.s1,Evalues,species.scores,variables,scale.fac,"RDA_s1_country.png",3,1)

site.constraints.s2 <- scores(rda.final, scaling=2, display=c("lc"))
species.scores.s2 <- scores(rda.final, scaling=2, display=c("species"))
rownames(species.scores.s2) <- substr(rownames(species.scores.s2),start=2,stop=nchar(rownames(species.scores.s2))-2)
variables.s2 <- as.data.frame(scores(rda.final, scaling=2, display=c("bp")))
variables.s2$variable <- c("environment","environment","linear","linear","MEM","MEM","MEM")
scale.fac.s2 <- attributes(species.scores.s2)$const
species.scores.s2 <- as.data.frame(species.scores.s2)

df.plot.s2 <- makePlotDF(site.constraints.s2, rownames(site.constraints.s2))
plotRDA(df.plot.s2,Evalues,species.scores.s2,variables.s2,scale.fac.s2,"RDA_s2_country.png",3,0)