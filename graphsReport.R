library(vegan)
library(tidyverse)
library(ggforce)
library(ggpubr)
library(ggnewscale)

source('functions.R')

plotPCoA <- function(df.plot, Evalues, option){
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

  return(p)
}

Fst.mat <- read.csv("Output/distancesFst_tot.csv",header=T,row.names=1) %>%
           as.matrix()

result.tot <- PCoA(Fst.mat)
Evalues.tot <- result.tot$values
Evectors.tot <- result.tot$vectors

df.plot.tot <- makePlotDF(Evectors.tot,rownames(Fst.mat))

p_Fst <- plotPCoA(df.plot.tot,Evalues.tot,3)

plotPCA <- function(df.plot, Evalues, U, text, scale, option){
  lambda <- sprintf("%.2f",round(Evalues[1:2] / sum(Evalues) *100, 2))
  rad <- scale*sqrt(2/length(Evalues))
  p <- ggplot() + geom_segment(data=U,x=rep(0,nrow(U)),y=rep(0,nrow(U)),xend=U[,1],yend=U[,2], 
              lineend='round', arrow = arrow(length = unit(0.1, "inches"))) + 
         geom_circle(aes(x0=0,y0=0,r=rad)) +
         geom_text(data=text,aes(x=x,y=y,label=names),size=3) +
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

  return(p)
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
df.text <- data.frame(x= c(0.2, 0.44, 0.51, 0.51, 0.45, -0.16, -0.21, -0.3, -0.28, -0.6) ,
                      y= c(0.35, 0.07, 0.03, -0.02, -0.09, -0.47, -0.37, -0.43, -0.23, -0.08),
                      names = c(15128, 29889, 81462, 15581, 11291, 65576, 6157, 65064, 53314, 58053))
p_freq <- plotPCA(df.plot.tot,PCA.tot$Evalues,PCA.tot$Umatrix,df.text,PCA.tot$Scale,3)

q <- ggarrange(p_Fst, p_freq, 
          labels = c("A", "B"),
          nrow = 2)
ggsave("PCoAvsPCA_tot.png",q,device="png",width=15,height=20,units='cm')

Fst.mat.atl <- Fst.mat[-c(32:38),-c(32:38)]

result.atl <- PCoA(Fst.mat.atl)
Evalues.atl <- result.atl$values
Evectors.atl <- result.atl$vectors
Evectors.atl[,1] <- -Evectors.atl[,1]

df.plot.atl <- makePlotDF(Evectors.atl,rownames(Fst.mat.atl))

p_Fst <- plotPCoA(df.plot.atl,Evalues.atl,3)

al.freq.atl <- al.freq.tot[!(row.names(al.freq.tot) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]
PCA.atl <- PCAOutput(al.freq.atl,0.7)
PCA.atl$Evectors[,2] <- -PCA.atl$Evectors[,2]
PCA.atl$Umatrix[,2] <- -PCA.atl$Umatrix[,2]
df.plot.atl <- makePlotDF(PCA.atl$Evectors, rownames(PCA.atl$Evectors))
df.text <- data.frame(x= c(-0.01, 0.22, 0.31, 0.25, 0.24, 0.74, 0.43) ,
                      y= c(0.28, 0.12, 0.1, 0.0505, 0.01, 0.09, -0.16),
                      names = c(39107, 6157, 65064, 39876, 65576, 53935, 42395))
p_freq <- plotPCA(df.plot.atl,PCA.atl$Evalues,PCA.atl$Umatrix,df.text,PCA.atl$Scale,3)

q <- ggarrange(p_Fst, p_freq, 
          labels = c("A", "B"),
          nrow = 2)
ggsave("PCoAvsPCA_atl.png",q,device="png",width=15,height=20,units='cm')

Fst.neut <- read.csv("Output/distancesFst_neut.csv",header=T,row.names=1) %>%
           as.matrix()

result.neut.tot <- PCoA(Fst.neut)
Evalues.neut.tot <- result.neut.tot$values
Evectors.neut.tot <- result.neut.tot$vectors
Evectors.neut.tot[,2] <- -Evectors.neut.tot[,2]

df.plot.neut.tot <- makePlotDF(Evectors.neut.tot,rownames(Fst.neut))

p_Fst_neut <- plotPCoA(df.plot.neut.tot,Evalues.neut.tot,3)

Fst.neut <- as.data.frame(Fst.neut)
Fst.neut.atl <- Fst.neut[!(row.names(Fst.neut) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]
Fst.neut.atl <- as.matrix(Fst.neut.atl %>% dplyr::select(-Laz,-Tar,-Sar,-Ale,-The,-Tor,-Sky))

result.neut.atl <- PCoA(Fst.neut.atl)
Evalues.neut.atl <- result.neut.atl$values
Evectors.neut.atl <- result.neut.atl$vectors
Evectors.neut.atl[,2] <- -Evectors.neut.atl[,2]

df.plot.neut.atl <- makePlotDF(Evectors.neut.atl,rownames(Fst.neut.atl))

p_Fst_neut_atl <- plotPCoA(df.plot.neut.atl,Evalues.neut.atl,3)

Fst.sel <- read.csv("Output/distancesFst_sel.csv",header=T,row.names=1) %>%
           as.matrix()

result.sel.tot <- PCoA(Fst.sel)
Evalues.sel.tot <- result.sel.tot$values
Evectors.sel.tot <- result.sel.tot$vectors
Evectors.sel.tot[,2] <- -Evectors.sel.tot[,2]

df.plot.sel.tot <- makePlotDF(Evectors.sel.tot,rownames(Fst.sel))

p_Fst_sel <- plotPCoA(df.plot.sel.tot,Evalues.sel.tot,3)

Fst.sel <- as.data.frame(Fst.sel)
Fst.sel.atl <- Fst.sel[!(row.names(Fst.sel) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]
Fst.sel.atl <- as.matrix(Fst.sel.atl %>% dplyr::select(-Laz,-Tar,-Sar,-Ale,-The,-Tor,-Sky))

result.sel.atl <- PCoA(Fst.sel.atl)
Evalues.sel.atl <- result.sel.atl$values
Evectors.sel.atl <- result.sel.atl$vectors
Evectors.sel.atl[,1:2] <- -Evectors.sel.atl[,1:2]

df.plot.sel.atl <- makePlotDF(Evectors.sel.atl,rownames(Fst.sel.atl))

p_Fst_sel_atl <- plotPCoA(df.plot.sel.atl,Evalues.sel.atl,3)

plotPCA <- function(df.plot, Evalues, U, scale, option){
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

  return(p )
}

al.freq.sel.tot <- read.csv("Output/allele_freq_sel.csv", header=T, row.names=1)
PCA.sel.tot <- PCAOutput(al.freq.sel.tot,0.01)
df.plot.sel.tot <- makePlotDF(PCA.sel.tot$Evectors, rownames(PCA.sel.tot$Evectors))
p_freq_sel <- plotPCA(df.plot.sel.tot,PCA.sel.tot$Evalues,PCA.sel.tot$Umatrix,PCA.sel.tot$Scale,3)

al.freq.sel.atl <- al.freq.sel.tot[!(row.names(al.freq.sel.tot) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]
PCA.sel.atl <- PCAOutput(al.freq.sel.atl,0.01)
df.plot.sel.atl <- makePlotDF(PCA.sel.atl$Evectors, rownames(PCA.sel.atl$Evectors))
p_freq_sel_atl <- plotPCA(df.plot.sel.atl,PCA.sel.atl$Evalues,PCA.sel.atl$Umatrix,PCA.sel.atl$Scale,3)

al.freq.neut.tot <- read.csv("Output/allele_freq_neut.csv", header=T, row.names=1)
PCA.neut.tot <- PCAOutput(al.freq.neut.tot,0.80)
df.plot.neut.tot <- makePlotDF(PCA.neut.tot$Evectors, rownames(PCA.neut.tot$Evectors))
p_freq_neut <- plotPCA(df.plot.neut.tot,PCA.neut.tot$Evalues,PCA.neut.tot$Umatrix,PCA.neut.tot$Scale,3)

al.freq.neut.atl <- al.freq.neut.tot[!(row.names(al.freq.neut.tot) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]
PCA.neut.atl <- PCAOutput(al.freq.neut.atl,0.70)
df.plot.neut.atl <- makePlotDF(PCA.neut.atl$Evectors, rownames(PCA.neut.atl$Evectors))
p_freq_neut_atl <- plotPCA(df.plot.neut.atl,PCA.neut.atl$Evalues,PCA.neut.atl$Umatrix,PCA.neut.atl$Scale,3)

q <- ggarrange(p_Fst_sel, p_freq_sel, p_Fst_sel_atl, p_freq_sel_atl, p_Fst_neut, p_freq_neut, p_Fst_neut_atl, p_freq_neut_atl,
        labels=c("A","B","C","D","E","F","G","H"),
        ncol=2,nrow=4)

ggsave("append_PCoAvsPCA.png",q,device='png',width=30,height=40,units='cm')

plotRDA <- function(df.plot, Evalues, VAR, scale, text, option, scaling){
  lambda <- sprintf("%.2f",round(Evalues[1:2] / sum(Evalues) *100, 2))
  rad <- scale*sqrt(2/length(Evalues))
  p <- ggplot() + geom_segment(data=VAR,aes(x=rep(0,nrow(VAR)),y=rep(0,nrow(VAR)),xend=VAR[,1],yend=VAR[,2],col=variable), 
              lineend='round', arrow = arrow(length = unit(0.1, "inches"))) + 
             scale_color_viridis_d(option='rocket') +
             new_scale_colour() +
             geom_circle(aes(x0=0,y0=0,r=rad),linetype=scaling) +        
         geom_text(data=text,aes(x=x,y=y,label=names),size=3) +
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
             xlim(1.3*min(df.plot$RDA1,-rad,VAR$RDA1),1.05*max(df.plot$RDA1,rad,VAR$RDA1)) +
             ylim(1.15*min(df.plot$RDA2,-rad,VAR$RDA2),1.15*max(df.plot$RDA2,rad,VAR$RDA2)) 

  return(p)
}

Fst.atl <- read.csv("Output/PCo/PCo.sel.atl.csv",header=T,row.names=1)

env <- read.csv("Output/EnvMatrix.csv",header=T,row.names=1)
env.atl <- env[!(row.names(env) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]

linear <- read.csv("Output/PCoSpatial.csv",header=T,row.names=1)[,1:2]

dbMEM <- read.csv("Output/dbMEM.csv",header=T,row.names=1)

AEM <- read.csv("Output/AEM.csv",header=T,row.names=1)

Fst.atl <- Fst.atl[match(rownames(AEM),rownames(Fst.atl)),]
env.atl <- env.atl[match(rownames(AEM),rownames(env.atl)),]
linear <- linear[match(rownames(AEM),rownames(linear)),]
dbMEM <- dbMEM[match(rownames(AEM),rownames(dbMEM)),]

env.atl.norm <- scale(env.atl, center=T, scale=T)

env.atl.norm.sel <- env.atl.norm[,1:3]

dbMEM.sel <- dbMEM[,c(1,3,5)]

rda.final <- rda(Fst.atl, cbind(env.atl.norm.sel,linear,dbMEM.sel))
rda.sum <- summary(rda.final)
Evalues <- as.vector(rda.sum$cont$importance[2,])
site.constraints <- scores(rda.final, scaling=1, display=c("lc"))
site.constraints[,1] <- -site.constraints[,1]
site.constraints[,2] <- -site.constraints[,2]
variables <- 0.4*scores(rda.final, scaling=1, display=c("bp"))
variables[,1] <- -variables[,1]
variables[,2] <- -variables[,2]
scale.fac <- 0.4*attributes(variables)$const
variables <- as.data.frame(variables)
variables$variable <- c("environment","environment","environment","linear","linear","MEM","MEM","MEM")

df.text <- data.frame(x=c(0.15, 0.29, 0.16, -0.35, -0.22, -0.32, -0.30, -0.15),
                      y=c(0.015, 0.005, -0.022, -0.01, -0.012, -0.002, 0.01, 0.022),
                      names=c("MEM5", "SST range", "PCo2", "PCo1", "SST mean", "MEM1", "SAL mean", "MEM3"))

#adjust xlim in plotRDA before plotting

df.plot.s1 <- makePlotDF(site.constraints, rownames(site.constraints))
p_s1 <- plotRDA(df.plot.s1,Evalues,variables,scale.fac,df.text,3,1)

site.constraints.s2 <- scores(rda.final, scaling=2, display=c("lc"))
site.constraints.s2[,1] <- -site.constraints.s2[,1]
site.constraints.s2[,2] <- -site.constraints.s2[,2]
variables.s2 <- scores(rda.final, scaling=2, display=c("bp"))
variables.s2[,1] <- -variables.s2[,1]
variables.s2[,2] <- -variables.s2[,2]
scale.fac.s2 <- attributes(variables.s2)$const
variables.s2 <- as.data.frame(variables.s2)
variables.s2$variable <- c("environment","environment","environment","linear","linear","MEM","MEM","MEM")

df.text <- data.frame(x=c(0.3, 0.82, 0.4, -0.97, -0.8, -0.92, -0.9, -0.35),
                      y=c(0.415, 0.125, -0.48, -0.24, -0.14, -0.11, 0.07, 0.45),
                      names=c("MEM5", "SST range", "PCo2", "PCo1", "SST mean", "MEM1", "SAL mean", "MEM3"))

df.plot.s2 <- makePlotDF(site.constraints.s2, rownames(site.constraints.s2))
p_s2 <- plotRDA(df.plot.s2,Evalues,variables.s2,scale.fac.s2,df.text,3,0)

q <- ggarrange(p_s1, p_s2, 
          labels = c("A", "B"),
          nrow = 2)
ggsave("RDA_Fst_sel.png",q,device="png",width=15,height=20,units='cm')

Fst.atl <- read.csv("Output/PCo/PCo.neut.atl.csv",header=T,row.names=1)
Fst.atl <- Fst.atl[match(rownames(AEM),rownames(Fst.atl)),]

dbMEM.sel <- dbMEM[,c(4)]
AEM.sel <- AEM[,c(1,10,12,15,17,18,19,22,27)]

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

df.text <- data.frame(x=c(0.005, 0.008, 0.005, 0.016, 0.0125, 0.0125, 0.0045, 0.011, 0.010, -0.01, -0.019, -0.013),
                      y=c(0.0075, 0.0055, 0.0014, 0.0016, -0.0006, -0.0015, -0.0017, -0.0023, -0.004, -0.012, -0.0008, 0.0006),
                      names=c("AEM1", "AEM27", "PCo2", "MEM4", "AEM10", "AEM17", "AEM19", "AEM22", "AEM18", "PCo1", "AEM15", "AEM12"))

df.plot.s1 <- makePlotDF(site.constraints, rownames(site.constraints))
p_s1 <- plotRDA(df.plot.s1,Evalues,variables,scale.fac,df.text,3,1)
ggsave("test.png", p_s1, device='png', height=10, width=15, units='cm')

site.constraints.s2 <- scores(rda.final, scaling=2, display=c("lc"))
site.constraints.s2[,2] <- -site.constraints.s2[,2]
variables.s2 <- 0.05*scores(rda.final, scaling=2, display=c("bp"))
variables.s2[,2] <- -variables.s2[,2]
scale.fac.s2 <- 0.05*attributes(variables.s2)$const
variables.s2 <- as.data.frame(variables.s2)
variables.s2$variable <- c("linear","linear","MEM",rep("AEM",9))

df.text <- data.frame(x=c(0.009, 0.013, 0.008, 0.027, 0.023, 0.023, 0.008, 0.019, 0.018, -0.0185, -0.035, -0.023),
                      y=c(0.028, 0.02, 0.005, 0.007, -0.002, -0.0055, -0.006, -0.0085, -0.0145, -0.045, -0.003, 0.002),
                      names=c("AEM1", "AEM27", "PCo2", "MEM4", "AEM10", "AEM17", "AEM19", "AEM22", "AEM18", "PCo1", "AEM15", "AEM12"))

df.plot.s2 <- makePlotDF(site.constraints.s2, rownames(site.constraints.s2))
p_s2 <- plotRDA(df.plot.s2,Evalues,variables.s2,scale.fac.s2,df.text,3,0)

q <- ggarrange(p_s1, p_s2, 
          labels = c("A", "B"),
          nrow = 2)
ggsave("RDA_Fst_neut.png",q,device="png",width=15,height=20,units='cm')


plotRDA <- function(df.plot, Evalues, U, VAR, scale, option, scaling){
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
             xlim(1.15*min(df.plot$RDA1,U$RDA1,-rad,VAR$RDA1),1.15*max(df.plot$RDA1,U$RDA1,rad,VAR$RDA1)) +
             ylim(1.15*min(df.plot$RDA2,U$RDA2,-rad,VAR$RDA2),1.15*max(df.plot$RDA2,U$RDA2,rad,VAR$RDA2)) 

  return(p)
}

freq <- read.csv("Output/allele_freq_sel.csv",header=T,row.names=1)
freq.atl <- freq[!(row.names(freq) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]
freq.atl <- freq.atl[match(rownames(AEM),rownames(freq.atl)),]

env.atl.norm.sel <- env.atl.norm[,c(1,2)]
dbMEM.sel <- dbMEM[,c(1,3,5)]

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
p_sel_s1 <- plotRDA(df.plot.s1,Evalues,species.scores,variables,scale.fac,3,1)

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
p_sel_s2 <- plotRDA(df.plot.s2,Evalues,species.scores.s2,variables.s2,scale.fac.s2,3,0)

freq <- read.csv("Output/allele_freq_neut.csv",header=T,row.names=1)
freq.atl <- freq[!(row.names(freq) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]
freq.atl <- freq.atl[match(rownames(AEM),rownames(freq.atl)),]

dbMEM.sel <- dbMEM[,c(1,3,4,5,7)]
AEM.sel <- AEM[,c(1,15,27)]

rda.final <- rda(freq.atl, cbind(linear,dbMEM.sel,AEM.sel))
rda.good <- as.data.frame(goodness(rda.final, model = "CCA", choices = 1:3))
SNP_list <- rownames(rda.good[rda.good$RDA2 > 0.5,])
rda.sum <- summary(rda.final)
Evalues <- as.vector(rda.sum$cont$importance[2,])
site.constraints <- scores(rda.final, scaling=1, display=c("lc"))
site.constraints[,1] <- -site.constraints[,1]
species.scores <- 0.5*scores(rda.final, scaling=1, display=c("species"))
scale.fac <- 0.5*attributes(species.scores)$const
species.scores <- as.data.frame(species.scores[rownames(species.scores) %in% SNP_list,])
rownames(species.scores) <- substr(rownames(species.scores),start=2,stop=nchar(rownames(species.scores))-2)
species.scores[,1] <- -species.scores[,1]
variables <- as.data.frame(0.4*scores(rda.final, scaling=1, display=c("bp")))
variables$variable <- c("linear","linear",rep("MEM",5),rep("AEM",3))
variables[,1] <- -variables[,1]

df.plot.s1 <- makePlotDF(site.constraints, rownames(site.constraints))
p_neut_s1 <- plotRDA(df.plot.s1,Evalues,species.scores,variables,scale.fac,3,1)

site.constraints.s2 <- scores(rda.final, scaling=2, display=c("lc"))
site.constraints.s2[,1] <- -site.constraints.s2[,1]
species.scores.s2 <- scores(rda.final, scaling=2, display=c("species"))
scale.fac.s2 <- attributes(species.scores.s2)$const
species.scores.s2 <- as.data.frame(species.scores.s2[rownames(species.scores.s2) %in% SNP_list,])
rownames(species.scores.s2) <- substr(rownames(species.scores.s2),start=2,stop=nchar(rownames(species.scores.s2))-2)
species.scores.s2[,1] <- -species.scores.s2[,1]
variables.s2 <- as.data.frame(scores(rda.final, scaling=2, display=c("bp")))
variables.s2$variable <- c("linear","linear",rep("MEM",5),rep("AEM",3))
variables.s2[,1] <- -variables.s2[,1]

df.plot.s2 <- makePlotDF(site.constraints.s2, rownames(site.constraints.s2))
p_neut_s2 <- plotRDA(df.plot.s2,Evalues,species.scores.s2,variables.s2,scale.fac.s2,3,0)

q <- ggarrange(p_sel_s1, p_sel_s2, p_neut_s1, p_neut_s2, 
          labels = c("A", "B", "C", "D"),
          nrow = 2, ncol=2)
ggsave("RDA_freq.png",q,device="png",width=30,height=20,units='cm')