library(vegan)
library(adegenet)
library(hierfstat)
library(tidyverse)
library(adespatial)

source('functions.R')

rdaSelectionFst <- function(gen.df,name){
	dat.new <- df2genind(gen.df[,-1],ncode=2,pop=gen.df[,1])
	print(dat.new)
	dat.hierf <- genind2hierfstat(dat.new)
	dist <- calcFst(dat.hierf)

	res <- PCoA(dist)
	Fst <- res$vectors
	Fst <- Fst[match(rownames(AEM),rownames(Fst)),]

	rda.env <- rda(Fst, env.atl.norm)
	R2adj.env <- RsquareAdj(rda.env)
	sel.env <- forward.sel(Fst, env.atl.norm, adjR2thresh=R2adj.env, alpha=0.05, nperm=9999)
	env.sel <- env.atl.norm[,sel.env$order]

	rda.dbMEM <- rda(Fst, dbMEM)
	R2adj.dbMEM <- RsquareAdj(rda.dbMEM)
	sel.dbMEM <- forward.sel(Fst, dbMEM, adjR2thresh=R2adj.dbMEM, alpha=0.05, nperm=9999)
	dbMEM.sel <- dbMEM[,sel.dbMEM$order]

	rda.AEM <- rda(Fst, AEM)
	R2adj.AEM <- RsquareAdj(rda.AEM)
	sel.AEM <- forward.sel(Fst, AEM, adjR2thresh=R2adj.AEM, alpha=0.05, nperm=9999)
	AEM.sel <- AEM[,sel.AEM$order]

	rda.env.lin.MEM.AEM <- rda(Fst,cbind(env.sel,linear,dbMEM.sel,AEM.sel))
	R2adj.env.lin.MEM.AEM <- RsquareAdj(rda.env.lin.MEM.AEM)
	An.env.lin.MEM.AEM <- anova(rda.env.lin.MEM.AEM,permutations = how(nperm=9999))

	rda.env_lin.MEM.AEM <- rda(Fst,env.sel,cbind(linear,dbMEM.sel,AEM.sel))
	R2adj.env_lin.MEM.AEM <- RsquareAdj(rda.env_lin.MEM.AEM)
	An.env_lin.MEM.AEM <- anova(rda.env_lin.MEM.AEM,permutations = how(nperm=9999))

	out <- list(SNP=name, R2_tot=R2adj.env.lin.MEM.AEM$r.squared, R2adj_tot=R2adj.env.lin.MEM.AEM$adj.r.squared,
				p_tot=An.env.lin.MEM.AEM[[4]][1], R2_mar=R2adj.env_lin.MEM.AEM$r.squared,
				R2adj_mar=R2adj.env_lin.MEM.AEM$adj.r.squared,p_mar=An.env_lin.MEM.AEM[[4]][1])
	print(as.data.frame(out))

	return(out)
}

dat.neut <- read.genepop("Input/lobster_1278ind_71snps_40pop.gen")
dat.df.neut <- genind2df(dat.neut)

rownames(dat.df.neut)[substr(rownames(dat.df.neut),1,2)=="SK"] <- paste0("Sky",substr(rownames(dat.df.neut[substr(rownames(dat.df.neut),1,2)=="SK",]),3,4))
rownames(dat.df.neut)[substr(rownames(dat.df.neut),1,2)=="TH"] <- paste0("The",substr(rownames(dat.df.neut[substr(rownames(dat.df.neut),1,2)=="TH",]),3,4))
rownames(dat.df.neut)[substr(rownames(dat.df.neut),1,2)=="AL"] <- paste0("Ale",substr(rownames(dat.df.neut[substr(rownames(dat.df.neut),1,2)=="AL",]),3,4))
rownames(dat.df.neut)[substr(rownames(dat.df.neut),1,2)=="CH"] <- paste0("Tor",substr(rownames(dat.df.neut[substr(rownames(dat.df.neut),1,2)=="CH",]),3,4))
rownames(dat.df.neut)[substr(rownames(dat.df.neut),1,2)=="Sr"] <- paste0("Sar17_",16:22)
rownames(dat.df.neut)[substr(rownames(dat.df.neut),1,3)=="Lro"] <- paste0("Idr17_",30:61)

levels(dat.df.neut$pop) <- c(levels(dat.df.neut$pop),"Sky49","The44","Ale89","Tor74")
dat.df.neut$pop[dat.df.neut$pop=="CH74"] <- "Tor74"
dat.df.neut$pop[dat.df.neut$pop=="SK49"] <- "Sky49"
dat.df.neut$pop[dat.df.neut$pop=="TH44"] <- "The44"
dat.df.neut$pop[dat.df.neut$pop=="AL89"] <- "Ale89"
dat.df.neut$pop[dat.df.neut$pop=="Sr7"] <- "Sar17_15"
dat.df.neut$pop[dat.df.neut$pop=="Lro33"] <- "Idr17_29"
dat.df.neut$pop <- droplevels(dat.df.neut$pop)

dat.df.neut <- dat.df.neut[!(dat.df.neut$pop %in% c("Laz7","Tar7","Sar17_15","Ale89","The44","Tor74","Sky49")), ]
dat.df.neut$pop <- droplevels(dat.df.neut$pop)

env <- read.csv("Output/EnvMatrix.csv",header=T,row.names=1)
env.atl <- env[!(row.names(env) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]
env.atl <- env.atl[,-c(6,7)]

linear <- read.csv("Output/PCoSpatial.csv",header=T,row.names=1)[,1:2]

dbMEM <- read.csv("Output/dbMEM.csv",header=T,row.names=1)

AEM <- read.csv("Output/AEM.csv",header=T,row.names=1)

env.atl <- env.atl[match(rownames(AEM),rownames(env.atl)),]
linear <- linear[match(rownames(AEM),rownames(linear)),]
dbMEM <- dbMEM[match(rownames(AEM),rownames(dbMEM)),]

env.atl.norm <- scale(env.atl, center=T, scale=T)

output <- rdaSelectionFst(dat.df.neut,"total")

for(i in 15:ncol(dat.df.neut)){
	dat.sel <- dat.df.neut[,-i]
	output <- rbind(output, rdaSelectionFst(dat.sel,names(dat.df.neut)[i]))
}

dat.sel <- dat.df.neut[,!(names(dat.df.neut) %in% c("39107"))]
dat.new <- df2genind(dat.sel[,-1],ncode=2,pop=gen.df[,1])
print(dat.new)
dat.hierf <- genind2hierfstat(dat.new)

dist <- calcFst(dat.hierf)

Fst <- PCoA(dist)$vectors
Fst <- Fst[match(rownames(AEM),rownames(Fst)),]

rda.env <- rda(Fst, env.atl.norm)
R2adj.env <- RsquareAdj(rda.env)
print(anova(rda.env,permutations = how(nperm=9999)))
sel.env <- forward.sel(Fst, env.atl.norm, adjR2thresh=R2adj.env, alpha=0.05, nperm=9999)
env.sel <- env.atl.norm[,sel.env$order]

rda.linear <- rda(Fst,linear)
R2adj.linear <- RsquareAdj(rda.linear)
print(anova(rda.linear,permutations = how(nperm=9999)))

rda.dbMEM <- rda(Fst, dbMEM)
R2adj.dbMEM <- RsquareAdj(rda.dbMEM)
print(anova(rda.dbMEM,permutations = how(nperm=9999)))
sel.dbMEM <- forward.sel(Fst, dbMEM, adjR2thresh=R2adj.env, alpha=0.05, nperm=9999)
dbMEM.sel <- dbMEM[,sel.dbMEM$order]

rda.AEM <- rda(Fst, AEM)
R2adj.AEM <- RsquareAdj(rda.AEM)
print(anova(rda.AEM,permutations = how(nperm=9999)))
sel.AEM <- forward.sel(Fst, AEM, adjR2thresh=R2adj.AEM, alpha=0.05, nperm=9999)
AEM.sel <- AEM[,sel.AEM$order]

rda.env_all <- rda(Fst, env.sel, cbind(linear,dbMEM.sel,AEM.sel))
print(anova(rda.env_all,permutations = how(nperm=9999)))

#####################################################################################################
freq.in <- read.csv("Output/allele_freq_neut.csv",header=T,row.names=1)
freq.atl <- freq[!(row.names(freq) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]
freq.atl <- freq.atl[match(rownames(AEM),rownames(freq.atl)),]



rdaSelectionFreq <- function(gen.df,name){
	freq <- gen.df

	rda.env <- rda(freq, env.atl.norm)
	R2adj.env <- RsquareAdj(rda.env)
	sel.env <- forward.sel(freq, env.atl.norm, adjR2thresh=R2adj.env, alpha=0.05, nperm=9999)
	env.sel <- env.atl.norm[,sel.env$order]

	rda.dbMEM <- rda(freq, dbMEM)
	R2adj.dbMEM <- RsquareAdj(rda.dbMEM)
	sel.dbMEM <- forward.sel(freq, dbMEM, adjR2thresh=R2adj.dbMEM, alpha=0.05, nperm=9999)
	dbMEM.sel <- dbMEM[,sel.dbMEM$order]

	rda.AEM <- rda(freq, AEM)
	R2adj.AEM <- RsquareAdj(rda.AEM)
	sel.AEM <- forward.sel(freq, AEM, adjR2thresh=R2adj.AEM, alpha=0.05, nperm=9999)
	AEM.sel <- AEM[,sel.AEM$order]

	rda.env.lin.MEM.AEM <- rda(freq,cbind(env.sel,linear,dbMEM.sel,AEM.sel))
	R2adj.env.lin.MEM.AEM <- RsquareAdj(rda.env.lin.MEM.AEM)
	An.env.lin.MEM.AEM <- anova(rda.env.lin.MEM.AEM,permutations = how(nperm=9999))

	rda.env_lin.MEM.AEM <- rda(freq,env.sel,cbind(linear,dbMEM.sel,AEM.sel))
	R2adj.env_lin.MEM.AEM <- RsquareAdj(rda.env_lin.MEM.AEM)
	An.env_lin.MEM.AEM <- anova(rda.env_lin.MEM.AEM,permutations = how(nperm=9999))

	out <- list(SNP=name, R2_tot=R2adj.env.lin.MEM.AEM$r.squared, R2adj_tot=R2adj.env.lin.MEM.AEM$adj.r.squared,
				p_tot=An.env.lin.MEM.AEM[[4]][1], R2_mar=R2adj.env_lin.MEM.AEM$r.squared,
				R2adj_mar=R2adj.env_lin.MEM.AEM$adj.r.squared,p_mar=An.env_lin.MEM.AEM[[4]][1])
	print(as.data.frame(out))

	return(out)
}

env.atl <- env.atl[,-c(6,7)]
env.atl.norm <- scale(env.atl, center=T, scale=T)

output <- rdaSelectionFreq(freq.atl,"total")

for(i in 1:ncol(freq.atl)){
	dat.sel <- freq.atl[,-i]
	output <- rbind(output, rdaSelectionFreq(dat.sel,names(freq.atl)[i]))
}

neutral_freq <- read.csv("Output/neutral_freq.csv", row.names=1, header=T)
name <- substr(rownames(neutral_freq)[-1],start=2,stop=nchar(rownames(neutral_freq)[-1])-2)
df.neut.freq <- data.frame(SNP=factor(name,levels=name),
							p=neutral_freq$p_mar[-1])

p_freq <- ggplot(data=df.neut.freq) + 
	geom_point(aes(x=SNP,y=p)) + geom_hline(yintercept=neutral_freq[1,]$p_mar,linetype="dashed") + theme_bw() +
	theme(axis.text.x = element_text(angle = 90, size = 7, hjust = 0, vjust=0.5), panel.grid.major.y = element_blank()) + 
	ylab("conditional p-value") + xlab("removed SNP") + ylim(0,1) + theme(plot.margin = margin(0.7,0.1,0.1,0.1, "cm"))

neutral_Fst <- read.csv("Output/neutral_Fst.csv", row.names=1, header=T)

df.neut.Fst <- data.frame(SNP=factor(rownames(neutral_Fst)[-1],levels=rownames(neutral_Fst)[-1]),
							p=neutral_Fst$p_mar[-1])

p_Fst <- ggplot(data=df.neut.Fst) + 
	geom_point(aes(x=SNP,y=p)) + geom_hline(yintercept=neutral_Fst[1,]$p_mar,linetype="dashed") + theme_bw() +
	theme(axis.text.x = element_text(angle = 90, size = 7, hjust = 0, vjust=0.5), panel.grid.major.y = element_blank()) + 
	ylab("conditional p-value") + xlab("removed SNP") + ylim(0,1)+ theme(plot.margin = margin(0.7,0.1,0.1,0.1, "cm"))

q <- ggarrange(p_Fst, p_freq,
				labels=c("A","B"),
				nrow=2, hjust=-2, vjust=1.5)

ggsave("neutral.png", q, device='png', width=18, height=10, units='cm')