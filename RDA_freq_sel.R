library(vegan)
library(adespatial)
library(tidyverse)
library(ggforce)
library(ggnewscale)
library(rdacca.hp)

# load the necessary functions
source('functions.R')

# load the allele frequencies with outlier loci for Atlantic sites
freq <- read.csv("Output/allele_freq_sel.csv",header=T,row.names=1)
freq.atl <- freq[!(row.names(freq) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]

# calculate the total variance of the allele frequency matrix
SSspec <- vector()
SStot <- vector()
for(i in 1:length(freq.atl[1,])){
  # mean of column i
  aver <- mean(freq.atl[,i])
  SSi <- 0
  # Sum of Squares of column i
  for(j in 1:length(freq.atl[,1])){
    SSi <- SSi + (freq.atl[j,i]-aver)^2
  }
  # list of variances of columns
  SSspec <- append(SSspec,SSi/length(freq.atl[,1]))
  # list of Sum of Squares of columns
  SStot <- append(SStot,SSi)
}

# beta diversity as total variance of matrix
# = total Sum of Squares divided by (number of columns - 1)
beta <- sum(SSspec)/(length(SSspec)-1)
SStotFin <- sum(SStot) 

# load the environmental data
env <- read.csv("Output/EnvMatrix.csv",header=T,row.names=1)
# only for Atlantic sites
env.atl <- env[!(row.names(env) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]
# remove Productivity variables (see report)
env.atl <- env.atl[,-c(6,7)]

# load the linear coordinates (PCo)
linear <- read.csv("Output/PCoSpatial.csv",header=T,row.names=1)[,1:2]

# load the dbMEMs
dbMEM <- read.csv("Output/dbMEM.csv",header=T,row.names=1)

# load the AEMs
AEM <- read.csv("Output/AEM.csv",header=T,row.names=1)

# ensure they all have the same order of rows
freq.atl <- freq.atl[match(rownames(AEM),rownames(freq.atl)),]
env.atl <- env.atl[match(rownames(AEM),rownames(env.atl)),]
linear <- linear[match(rownames(AEM),rownames(linear)),]
dbMEM <- dbMEM[match(rownames(AEM),rownames(dbMEM)),]

# scale the environmental variables (standardise)
env.atl.norm <- scale(env.atl, center=T, scale=T)

# RDA with only environment
rda.env <- rda(freq.atl, env.atl.norm)
R2adj.env <- RsquareAdj(rda.env)
print(anova(rda.env,permutations = how(nperm=9999)))
# forward selection of the environmental variables
# R2thresh: R2 < R2 of the full model
# alpha < 0.05 for each individual variable
sel.env <- forward.sel(freq.atl, env.atl.norm, adjR2thresh=R2adj.env, alpha=0.05, nperm=9999)
print(sel.env)
# matrix of selected environmental variables
env.atl.norm.sel <- env.atl.norm[,sel.env$order]

# RDA with only linear coordinates (both are kept)
rda.linear <- rda(freq.atl,linear)
R2adj.linear <- RsquareAdj(rda.linear)
print(anova(rda.linear,permutations = how(nperm=9999)))

# RDA with only dbMEMs
rda.dbMEM <- rda(freq.atl, dbMEM)
R2adj.dbMEM <- RsquareAdj(rda.dbMEM)
print(anova(rda.dbMEM,permutations = how(nperm=9999)))
# forward selection of the dbMEMs (full model is not considered)
sel.dbMEM <- forward.sel(freq.atl, dbMEM, adjR2thresh=R2adj.dbMEM, alpha=0.05, nperm=9999)
print(sel.dbMEM)
# matrix of selected dbMEMs
dbMEM.sel <- dbMEM[,sel.dbMEM$order]

# RDA with only AEMs
rda.AEM <- rda(freq.atl, AEM)
R2adj.AEM <- RsquareAdj(rda.AEM)
print(anova(rda.AEM,permutations = how(nperm=9999)))
# forward selection of the AEMs (full model is not considered)
sel.AEM <- forward.sel(freq.atl, AEM, adjR2thresh=R2adj.AEM, alpha=0.05, nperm=9999)
print(sel.AEM)
# no AEMs are selected (full model significance > 0.05)

# initialise containers
fractions <- vector()
R2adj <- vector()
degrees <- vector()
Fvalues <- vector()
pvalues <- vector()

# do variation partitioning with the selected functions
# variation partitioning with varpart
rda.part <- varpart(freq.atl, env.atl.norm.sel, linear, dbMEM.sel)
# extract the parts with R2 and testable
rda.part.df <- as.data.frame(rbind(rda.part$part$fract,rda.part$part$indfract,rda.part$part$contr1))
# only keep the testable parts
rda.part.df <- rda.part.df[rda.part.df$Testable==T,]
# manually add row names
rownames(rda.part.df) <- c("[env]","[lin]","[MEM]","[env] + [lin]","[env] + [MEM]","[lin] + [MEM]",
	"[env] + [lin] + [MEM]","[env] | [lin] + [MEM]","[lin] | [env] + [MEM]","[MEM] | [env] + [lin]",
	"[env] | [MEM]","[env] | [lin]","[lin] | [MEM]","[lin] | [env]","[MEM] | [env]","[MEM] | [lin]")

# use key list
key <- list("[env]"=env.atl.norm.sel, "[lin]"=linear, "[MEM]"=dbMEM.sel)

# test the testable fractions
for(i in 1:nrow(rda.part.df)){
  # name of the fraction
	fractions <- c(fractions, rownames(rda.part.df)[i])
  # R2 value
	R2adj <- c(R2adj, rda.part.df$Adj.R.square[i])

  # construct the RDA corresponding to the fraction
	rda <- makeRDA(rownames(rda.part.df)[i])
	print(summary(rda))
  # test the fraction (permutation)
	an <- anova(rda, permutations = how(nperm=9999))
  # number of degrees of freedom
	deg <- paste0(an$Df[1],"(",an$Df[1]+an$Df[2],")")
	degrees <- c(degrees, deg)
  # F-value of the test
	Fvalues <- c(Fvalues, an$F[1])
  # p-value of the test
	pvalues <- c(pvalues, an$"Pr(>F)"[1])
}

# put everything together in output data.frame
df.parts <- data.frame(fraction=fractions, RsquareAdj=R2adj, dof=degrees, F=Fvalues, "p-value"=pvalues)
print(df.parts)

# plot the output of the RDA
# input = plot data.frame (see functions.R)
#         eigenvalues Evalues
#         species scores U
#         triplot scores (variables) VAR
#         scaling factor
#         fileName to give to the output png
#         option: =3 colour according to country
#                 =4 colour according to region
#         scaling: =0 scaling 2 in vegan
#                  =1 scaling 1 in vegan
plotRDA <- function(df.plot, Evalues, U, VAR, scale, fileName, option, scaling){
  # proportion of the variance explained by the axis
  lambda <- sprintf("%.2f",round(Evalues[1:2] / sum(Evalues) *100, 2))
  # radius of the circle of equal variance
  rad <- scale*sqrt(2/length(Evalues))
  p <- ggplot() + 
              # plot the species scores
              geom_segment(data=U,x=rep(0,nrow(U)),y=rep(0,nrow(U)),xend=U[,1],yend=U[,2], 
             	lineend='round', arrow = arrow(length = unit(0.1, "inches"))) +
              # plot the variable scores
             geom_segment(data=VAR,aes(x=rep(0,nrow(VAR)),y=rep(0,nrow(VAR)),xend=VAR[,1],yend=VAR[,2],col=variable), 
             	linetype='dashed',lineend='round', arrow = arrow(length = unit(0.1, "inches"))) + 
             # plot the variable labels
             geom_text(data=VAR,aes(x=1.12*RDA1,y=1.07*RDA2,label=rownames(VAR)),size=3) +
             # viridis scheme for the variables
             scale_color_viridis_d(option='rocket') +
             # reset the colour scale (ggnewscale)
             new_scale_colour() +
             # circle of equal variance (if scaling = 1)
             geom_circle(aes(x0=0,y0=0,r=rad),linetype=scaling) +  			 
             # labels for the species scores
  			 geom_text(data=U,aes(x=1.12*RDA1,y=1.07*RDA2,label=rownames(U)),size=3) +
         # plot the site scores
  			 geom_point(data=df.plot,aes(x=RDA1,y=RDA2,col=df.plot[,option],shape=df.plot[,option]),size=1.5) +
             theme_bw() + 
             # colour scheme for the sites (see functions.R)
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
             # legend title
             guides(col=guide_legend(title=colnames(df.plot)[option])) +
             # manually adjust limits
             xlim(1.15*min(df.plot$RDA1,U$RDA1,-rad,VAR$RDA1),1.3*max(df.plot$RDA1,U$RDA1,rad,VAR$RDA1)) +
             ylim(1.15*min(df.plot$RDA2,U$RDA2,-rad,VAR$RDA2),1.15*max(df.plot$RDA2,U$RDA2,rad,VAR$RDA2)) 

  p %>% ggsave(fileName,.,device='png',width=15,height=10,units='cm')
}

# do the total final RDA with the selected variables
rda.final <- rda(freq.atl, cbind(env.atl.norm.sel,linear,dbMEM.sel))
rda.sum <- summary(rda.final)
Evalues <- as.vector(rda.sum$cont$importance[2,])
# site scores = linear combinations of explanatory variables
site.constraints <- scores(rda.final, scaling=1, display=c("lc"))
# manually reverse PC1 (for consistency)
site.constraints[,1] <- -site.constraints[,1]
# manually scale species scores
species.scores <- 0.4*scores(rda.final, scaling=1, display=c("species"))
# adjust locus names
rownames(species.scores) <- substr(rownames(species.scores),start=2,stop=nchar(rownames(species.scores))-2)
# manually reverse PC1
species.scores[,1] <- -species.scores[,1]
# manually scale variable scores
variables <- as.data.frame(0.4*scores(rda.final, scaling=1, display=c("bp")))
# name variables
variables$variable <- c("environment","environment","linear","linear","MEM","MEM","MEM")
# manually reverse PC1
variables[,1] <- -variables[,1]
# extract scale factor ! before making data.frame of species.scores
scale.fac <- 0.4*attributes(species.scores)$const
species.scores <- as.data.frame(species.scores)

df.plot.s1 <- makePlotDF(site.constraints, rownames(site.constraints))
plotRDA(df.plot.s1,Evalues,species.scores,variables,scale.fac,"RDA_freq_sel_s1_country.png",3,1)

# also for scaling 2
site.constraints.s2 <- scores(rda.final, scaling=2, display=c("lc"))
# manually reverse PC1
site.constraints.s2[,1] <- -site.constraints.s2[,1]
species.scores.s2 <- scores(rda.final, scaling=2, display=c("species"))
rownames(species.scores.s2) <- substr(rownames(species.scores.s2),start=2,stop=nchar(rownames(species.scores.s2))-2)
# manually reverse PC1
species.scores.s2[,1] <- -species.scores.s2[,1]
variables.s2 <- as.data.frame(scores(rda.final, scaling=2, display=c("bp")))
variables.s2$variable <- c("environment","environment","linear","linear","MEM","MEM","MEM")
# manually reverse PC1
variables.s2[,1] <- -variables.s2[,1]
scale.fac.s2 <- attributes(species.scores.s2)$const
species.scores.s2 <- as.data.frame(species.scores.s2)

df.plot.s2 <- makePlotDF(site.constraints.s2, rownames(site.constraints.s2))
plotRDA(df.plot.s2,Evalues,species.scores.s2,variables.s2,scale.fac.s2,"RDA_freq_sel_s2_country.png",3,0)

# PCA of the unconstrained variation
# PCA on residuals
# axes 8 and 9 are the first unconstrained axes
site.unconstrained <- as.data.frame(scores(rda.final, scaling=1, choices=c(8,9), display=c("sites")))
site.unconstrained$country <- countryList(rownames(site.unconstrained))
Evalues <- as.vector(rda.sum$cont$importance[2,8:9])
lambda <- sprintf("%.2f",round(Evalues[1:2]*100, 2))

# plot this unconstrained PCA
p <- ggplot() + geom_point(data=site.unconstrained,aes(x=PC1,y=PC2,col=country),size=1.5) +
             theme_bw() + scale_color_viridis_d(option='magma') +
             xlab(paste0("PC1 (",lambda[1]," %)")) + ylab(paste0("PC2 (",lambda[2]," %)")) +
             theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
             geom_hline(yintercept=0,linetype="dashed") + geom_vline(xintercept=0,linetype="dashed")              

p %>% ggsave("resid_PCA_freq_sel.png",.,device='png',width=15,height=10,units='cm')

# hierarchical variation partitioning
rdacca <- rdacca.hp(freq.atl,list(as.data.frame(env.atl.norm.sel), as.data.frame(linear), as.data.frame(dbMEM.sel)), method='RDA')
rdacca.sep <- rdacca.hp(freq.atl,cbind(env.atl.norm.sel,linear,dbMEM.sel), method='RDA')

permutest <- permu.hp(freq.atl,list(as.data.frame(env.atl.norm.sel), as.data.frame(linear), as.data.frame(dbMEM.sel)), method='RDA', permutations=9999)
# !!!!!!!!!!!!!!!! permu.hp(individual variables) takes a long time to run
permutest.sep <- permu.hp(freq.atl,cbind(env.atl.norm.sel,linear,dbMEM.sel), method='RDA', permutations=9999)

p <- plot.rdaccahp(rdacca.sep) + theme_bw()
p %>% ggsave("rdacca_sep_freq_sel.png",.,device='png',width=15,height=7.5,units='cm')