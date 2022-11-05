library(vegan)
library(adespatial)
library(tidyverse)
library(ggforce)
library(ggnewscale)
library(rdacca.hp)

source('functions.R')

freq <- read.csv("Output/allele_freq_sel.csv",header=T,row.names=1)
freq.atl <- freq[!(row.names(freq) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]

SSspec <- vector()
SStot <- vector()
for(i in 1:length(freq.atl[1,])){
  aver <- mean(freq.atl[,i])
  SSi <- 0
  for(j in 1:length(freq.atl[,1])){
    SSi <- SSi + (freq.atl[j,i]-aver)^2
  }
  SSspec <- append(SSspec,SSi/length(freq.atl[,1]))
  SStot <- append(SStot,SSi)
}

beta <- sum(SSspec)/(length(SSspec)-1)
SStotFin <- sum(SStot) 

env <- read.csv("Output/EnvMatrix.csv",header=T,row.names=1)
env.atl <- env[!(row.names(env) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]
env.atl <- env.atl[,-c(6,7)]

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
             geom_text(data=VAR,aes(x=1.12*RDA1,y=1.07*RDA2,label=rownames(VAR)),size=3) +
             scale_color_viridis_d(option='rocket') +
             new_scale_colour() +
             geom_circle(aes(x0=0,y0=0,r=rad),linetype=scaling) +  			 
  			 geom_text(data=U,aes(x=1.12*RDA1,y=1.07*RDA2,label=rownames(U)),size=3) +
  			 geom_point(data=df.plot,aes(x=RDA1,y=RDA2,col=df.plot[,option],shape=df.plot[,option]),size=1.5) +
             theme_bw() + 
			 scale_colour_manual(
             	name = names(df.plot)[option],
             	labels = sort(unique(df.plot[,option])),
             	values=getColours(df.plot[,option])$colour) +
             scale_shape_manual(
             	name = names(df.plot)[option],
             	labels =sort(unique(df.plot[,option])),
             	values=getColours(df.plot[,option])$shape) +
             xlab(paste0("RDA1 (",lambda[1]," %)")) + ylab(paste0("RDA2 (",lambda[2]," %)")) +
             theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
             geom_hline(yintercept=0,linetype="dashed") + geom_vline(xintercept=0,linetype="dashed") +
             guides(col=guide_legend(title=colnames(df.plot)[option])) +
             xlim(1.15*min(df.plot$RDA1,U$RDA1,-rad,VAR$RDA1),1.3*max(df.plot$RDA1,U$RDA1,rad,VAR$RDA1)) +
             ylim(1.15*min(df.plot$RDA2,U$RDA2,-rad,VAR$RDA2),1.15*max(df.plot$RDA2,U$RDA2,rad,VAR$RDA2)) 

  p %>% ggsave(fileName,.,device='png',width=15,height=10,units='cm')
}

rda.final <- rda(freq.atl, cbind(env.atl.norm.sel,linear,dbMEM.sel))
rda.sum <- summary(rda.final)
Evalues <- as.vector(rda.sum$cont$importance[2,])
site.constraints <- scores(rda.final, scaling=1, display=c("lc"))
site.constraints[,1] <- -site.constraints[,1]
species.scores <- 0.4*scores(rda.final, scaling=1, display=c("species"))
rownames(species.scores) <- substr(rownames(species.scores),start=2,stop=nchar(rownames(species.scores))-2)
species.scores[,1] <- -species.scores[,1]
variables <- as.data.frame(0.4*scores(rda.final, scaling=1, display=c("bp")))
variables$variable <- c("environment","environment","linear","linear","MEM","MEM","MEM")
variables[,1] <- -variables[,1]
scale.fac <- 0.4*attributes(species.scores)$const
species.scores <- as.data.frame(species.scores)

df.plot.s1 <- makePlotDF(site.constraints, rownames(site.constraints))
plotRDA(df.plot.s1,Evalues,species.scores,variables,scale.fac,"RDA_freq_sel_s1_country.png",3,1)

site.constraints.s2 <- scores(rda.final, scaling=2, display=c("lc"))
site.constraints.s2[,1] <- -site.constraints.s2[,1]
species.scores.s2 <- scores(rda.final, scaling=2, display=c("species"))
rownames(species.scores.s2) <- substr(rownames(species.scores.s2),start=2,stop=nchar(rownames(species.scores.s2))-2)
species.scores.s2[,1] <- -species.scores.s2[,1]
variables.s2 <- as.data.frame(scores(rda.final, scaling=2, display=c("bp")))
variables.s2$variable <- c("environment","environment","linear","linear","MEM","MEM","MEM")
variables.s2[,1] <- -variables.s2[,1]
scale.fac.s2 <- attributes(species.scores.s2)$const
species.scores.s2 <- as.data.frame(species.scores.s2)

df.plot.s2 <- makePlotDF(site.constraints.s2, rownames(site.constraints.s2))
plotRDA(df.plot.s2,Evalues,species.scores.s2,variables.s2,scale.fac.s2,"RDA_freq_sel_s2_country.png",3,0)

site.unconstrained <- as.data.frame(scores(rda.final, scaling=1, choices=c(8,9), display=c("sites")))
site.unconstrained$country <- countryList(rownames(site.unconstrained))
Evalues <- as.vector(rda.sum$cont$importance[2,8:9])
lambda <- sprintf("%.2f",round(Evalues[1:2]*100, 2))

p <- ggplot() + geom_point(data=site.unconstrained,aes(x=PC1,y=PC2,col=country),size=1.5) +
             theme_bw() + scale_color_viridis_d(option='magma') +
             xlab(paste0("PC1 (",lambda[1]," %)")) + ylab(paste0("PC2 (",lambda[2]," %)")) +
             theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
             geom_hline(yintercept=0,linetype="dashed") + geom_vline(xintercept=0,linetype="dashed")              

p %>% ggsave("resid_PCA_freq_sel.png",.,device='png',width=15,height=10,units='cm')

rdacca <- rdacca.hp(freq.atl,list(as.data.frame(env.atl.norm.sel), as.data.frame(linear), as.data.frame(dbMEM.sel)), method='RDA')
rdacca.sep <- rdacca.hp(freq.atl,cbind(env.atl.norm.sel,linear,dbMEM.sel), method='RDA')

permutest <- permu.hp(freq.atl,list(as.data.frame(env.atl.norm.sel), as.data.frame(linear), as.data.frame(dbMEM.sel)), method='RDA', permutations=9999)
permutest.sep <- permu.hp(freq.atl,cbind(env.atl.norm.sel,linear,dbMEM.sel), method='RDA', permutations=9999)

p <- plot.rdaccahp(rdacca.sep) + theme_bw()
p %>% ggsave("rdacca_sep_freq_sel.png",.,device='png',width=15,height=7.5,units='cm')