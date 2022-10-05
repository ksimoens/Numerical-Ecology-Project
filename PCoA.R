library(tidyverse)

Fst.mat <- read.csv("Output/distancesFst_tot.csv",header=T,row.names=1) %>%
           as.matrix()

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
	return(list(values=Evalues,vectors=Evectors))
}

result.tot <- PCoA(Fst.mat)
Evalues.tot <- result.tot$values
Evectors.tot <- result.tot$vectors


Evectors.tot <- Evectors.tot[,-c(37:38)]
Evalues.tot <- Evalues.tot[-c(37:38)]
for(i in 1:ncol(Evectors.tot)){
  Evectors.tot[,i] <- Evectors.tot[,i]*Evalues.tot[i]
}
print(Evectors.tot)

df.plot.tot <- Evectors.tot[,1:2] %>%
               as.data.frame()
rownames(df.plot.tot) <- rownames(Fst.mat)
countryList <- c("NOR","NOR","NOR","SVE","SVE","SVE","NOR","DEU","NDL","GRB","GRB","GRB","GRB","GRB","GRB","GRB","GRB","CHA","GRB","GRB","GRB","GRB","GRB","GRB","IRL","IRL","IRL","IRL","IRL","FRA","ESP","ITA","ITA","ITA","HEL","HEL","HEL","HEL")
df.plot.tot$country <- countryList
regionList <- c("ATL","ATL",rep("SKA",5),rep("ATL",24),rep("MED",7))
df.plot.tot$region <- regionList

print(df.plot.tot)

lambda.tot <- Evalues.tot[1:2] / sum(Evalues.tot)

print(lambda.tot)

p.tot.region <- df.plot.tot %>%
  ggplot() + geom_point(aes(x=PCo1,y=PCo2,col=region),size=1.5) +
             theme_bw() + scale_color_viridis_d() +
             xlab("PCo1 (27.74 %)") + ylab("PCo2 (8.17 %)") +
             theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
             geom_hline(yintercept=0,linetype="dashed") + geom_vline(xintercept=0,linetype="dashed") 

p.tot.region %>% ggsave("PCoA.tot.region.png",.,device='png',width=15,height=10,units='cm')

p.tot.country <- df.plot.tot %>%
  ggplot() + geom_point(aes(x=PCo1,y=PCo2,col=country),size=1.5) +
             theme_bw() + scale_color_viridis_d(option='magma') +
             xlab("PCo1 (27.74 %)") + ylab("PCo2 (8.17 %)") +
             theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
             geom_hline(yintercept=0,linetype="dashed") + geom_vline(xintercept=0,linetype="dashed") 

p.tot.country %>% ggsave("PCoA.tot.country.png",.,device='png',width=15,height=10,units='cm')

Fst.mat.atl <- Fst.mat[-c(32:38),-c(32:38)]

result.atl <- PCoA(Fst.mat.atl)
Evalues.atl <- result.atl$values
Evectors.atl <- result.atl$vectors

Evectors.atl <- Evectors.atl[,-c(30:31)]
Evalues.atl <- Evalues.atl[-c(30:31)]
for(i in 1:ncol(Evectors.atl)){
  Evectors.atl[,i] <- Evectors.atl[,i]*Evalues.atl[i]
}

df.plot.atl <- Evectors.atl[,1:2] %>%
               as.data.frame()
rownames(df.plot.atl) <- rownames(Fst.mat.atl)
countryList <- c("NOR","NOR","NOR","SVE","SVE","SVE","NOR","DEU","NDL","GRB","GRB","GRB","GRB","GRB","GRB","GRB","GRB","CHA","GRB","GRB","GRB","GRB","GRB","GRB","IRL","IRL","IRL","IRL","IRL","FRA","ESP")
df.plot.atl$country <- countryList
regionList <- c("ATL","ATL",rep("SKA",5),rep("ATL",24))
df.plot.atl$region <- regionList

print(df.plot.atl)

lambda.atl <- Evalues.atl[1:2] / sum(Evalues.atl)

print(lambda.atl)

p.atl.region <- df.plot.atl %>%
  ggplot() + geom_point(aes(x=PCo1,y=PCo2,col=region),size=1.5) +
             theme_bw() + scale_color_viridis_d() +
             xlab("PCo1 (17.51 %)") + ylab("PCo2 (5.22 %)") +
             theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
             geom_hline(yintercept=0,linetype="dashed") + geom_vline(xintercept=0,linetype="dashed") 

p.atl.region %>% ggsave("PCoA.atl.region.png",.,device='png',width=15,height=10,units='cm')

p.atl.country <- df.plot.atl %>%
  ggplot() + geom_point(aes(x=PCo1,y=PCo2,col=country),size=1.5) +
             theme_bw() + scale_color_viridis_d(option='magma') +
             xlab("PCo1 (17.51 %)") + ylab("PCo2 (5.22 %)") +
             theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
             geom_hline(yintercept=0,linetype="dashed") + geom_vline(xintercept=0,linetype="dashed") 

p.atl.country %>% ggsave("PCoA.atl.country.png",.,device='png',width=15,height=10,units='cm')

Fst.neut <- read.csv("Output/distancesFst_neut.csv",header=T,row.names=1) %>%
           as.matrix()

result.neut.tot <- PCoA(Fst.neut)
Evalues.neut.tot <- result.neut.tot$values
Evectors.neut.tot <- result.neut.tot$vectors

Evectors.neut.tot <- Evectors.neut.tot[,-c(37:38)]
Evalues.neut.tot <- Evalues.neut.tot[-c(37:38)]
for(i in 1:ncol(Evectors.neut.tot)){
  Evectors.neut.tot[,i] <- Evectors.neut.tot[,i]*Evalues.neut.tot[i]
}
print(Evectors.neut.tot)

df.plot.neut.tot <- Evectors.neut.tot[,1:2] %>%
               as.data.frame()
rownames(df.plot.neut.tot) <- rownames(Fst.neut)
print(df.plot.neut.tot)

countryList <- c("HEL","NOR","GRB","IRL","GRB","GRB","NOR","SVE","GRB","DEU","IRL","FRA","GRB","GRB","CHA","SVE","IRL","ITA","GRB","GRB","SVE","IRL","NDL","GRB","GRB","GRB","ITA","GRB","GRB","NOR","HEL","GRB","ITA","HEL","HEL","NOR","IRL","ESP")
df.plot.neut.tot$country <- countryList
regionList <- c("MED","ATL",rep("ATL",4),"SKA","SKA",rep("ATL",7),"SKA","ATL","MED","ATL","ATL","SKA",rep("ATL",5),"MED","ATL","ATL","SKA","MED","ATL",rep("MED",3),"ATL","ATL","ATL")
df.plot.neut.tot$region <- regionList

print(df.plot.neut.tot)

lambda.neut.tot <- Evalues.neut.tot[1:2] / sum(Evalues.neut.tot)

print(lambda.neut.tot)

p.neut.tot.region <- df.plot.neut.tot %>%
  ggplot() + geom_point(aes(x=PCo1,y=PCo2,col=region),size=1.5) +
             theme_bw() + scale_color_viridis_d() +
             xlab("PCo1 (12.76 %)") + ylab("PCo2 (6.59 %)") +
             theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
             geom_hline(yintercept=0,linetype="dashed") + geom_vline(xintercept=0,linetype="dashed") 

p.neut.tot.region %>% ggsave("PCoA.neut.tot.region.png",.,device='png',width=15,height=10,units='cm')

p.neut.tot.country <- df.plot.neut.tot %>%
  ggplot() + geom_point(aes(x=PCo1,y=PCo2,col=country),size=1.5) +
             theme_bw() + scale_color_viridis_d(option='magma') +
             xlab("PCo1 (12.76 %)") + ylab("PCo2 (6.59 %)") +
             theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
             geom_hline(yintercept=0,linetype="dashed") + geom_vline(xintercept=0,linetype="dashed") 

p.neut.tot.country %>% ggsave("PCoA.neut.tot.country.png",.,device='png',width=15,height=10,units='cm')

Fst.neut <- as.data.frame(Fst.neut)
Fst.neut.atl <- Fst.neut[!(row.names(Fst.neut) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]
Fst.neut.atl <- as.matrix(Fst.neut.atl %>% select(-Laz,-Tar,-Sar,-Ale,-The,-Tor,-Sky))

result.neut.atl <- PCoA(Fst.neut.atl)
Evalues.neut.atl <- result.neut.atl$values
Evectors.neut.atl <- result.neut.atl$vectors

print(Evalues.neut.atl)

Evectors.neut.atl <- Evectors.neut.atl[,-c(30:31)]
Evalues.neut.atl <- Evalues.neut.atl[-c(30:31)]
for(i in 1:ncol(Evectors.neut.atl)){
  Evectors.neut.atl[,i] <- Evectors.neut.atl[,i]*Evalues.neut.atl[i]
}
print(Evectors.neut.atl)

df.plot.neut.atl <- Evectors.neut.atl[,1:2] %>%
               as.data.frame()
rownames(df.plot.neut.atl) <- rownames(Fst.neut.atl)
print(df.plot.neut.atl)

countryList <- c("NOR","GRB","IRL","GRB","GRB","NOR","SVE","GRB","DEU","IRL","FRA","GRB","GRB","CHA","SVE","IRL","GRB","GRB","SVE","IRL","NDL","GRB","GRB","GRB","GRB","GRB","NOR","GRB","NOR","IRL","ESP")
df.plot.neut.atl$country <- countryList
regionList <- c("ATL",rep("ATL",4),"SKA","SKA",rep("ATL",7),"SKA","ATL","ATL","ATL","SKA",rep("ATL",5),"ATL","ATL","SKA","ATL","ATL","ATL","ATL")
df.plot.neut.atl$region <- regionList

print(df.plot.neut.atl)

lambda.neut.atl <- Evalues.neut.atl[1:2] / sum(Evalues.neut.atl)

print(lambda.neut.atl)

p.neut.atl.region <- df.plot.neut.atl %>%
  ggplot() + geom_point(aes(x=PCo1,y=PCo2,col=region),size=1.5) +
             theme_bw() + scale_color_viridis_d() +
             xlab("PCo1 (13.04 %)") + ylab("PCo2 (6.41 %)") +
             theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
             geom_hline(yintercept=0,linetype="dashed") + geom_vline(xintercept=0,linetype="dashed") 

p.neut.atl.region %>% ggsave("PCoA.neut.atl.region.png",.,device='png',width=15,height=10,units='cm')

p.neut.atl.country <- df.plot.neut.atl %>%
  ggplot() + geom_point(aes(x=PCo1,y=PCo2,col=country),size=1.5) +
             theme_bw() + scale_color_viridis_d(option='magma') +
             xlab("PCo1 (13.04 %)") + ylab("PCo2 (6.41 %)") +
             theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
             geom_hline(yintercept=0,linetype="dashed") + geom_vline(xintercept=0,linetype="dashed") 

p.neut.atl.country %>% ggsave("PCoA.neut.atl.country.png",.,device='png',width=15,height=10,units='cm')

Fst.sel <- read.csv("Output/distancesFst_sel.csv",header=T,row.names=1) %>%
           as.matrix()

result.sel.tot <- PCoA(Fst.sel)
Evalues.sel.tot <- result.sel.tot$values
Evectors.sel.tot <- result.sel.tot$vectors

Evectors.sel.tot <- Evectors.sel.tot[,-c(37:38)]
Evalues.sel.tot <- Evalues.sel.tot[-c(37:38)]
for(i in 1:ncol(Evectors.sel.tot)){
  Evectors.sel.tot[,i] <- Evectors.sel.tot[,i]*Evalues.sel.tot[i]
}
print(Evectors.sel.tot)

df.plot.sel.tot <- Evectors.sel.tot[,1:2] %>%
               as.data.frame()
rownames(df.plot.sel.tot) <- rownames(Fst.sel)
print(df.plot.sel.tot)

countryList <- c("HEL","NOR","GRB","IRL","GRB","GRB","NOR","SVE","GRB","DEU","IRL","FRA","GRB","GRB","CHA","SVE","IRL","ITA","GRB","GRB","SVE","IRL","NDL","GRB","GRB","GRB","ITA","GRB","GRB","NOR","HEL","GRB","ITA","HEL","HEL","NOR","IRL","ESP")
df.plot.sel.tot$country <- countryList
regionList <- c("MED","ATL",rep("ATL",4),"SKA","SKA",rep("ATL",7),"SKA","ATL","MED","ATL","ATL","SKA",rep("ATL",5),"MED","ATL","ATL","SKA","MED","ATL",rep("MED",3),"ATL","ATL","ATL")
df.plot.sel.tot$region <- regionList

print(df.plot.sel.tot)

lambda.sel.tot <- Evalues.sel.tot[1:2] / sum(Evalues.sel.tot)

print(lambda.sel.tot)

p.sel.tot.region <- df.plot.sel.tot %>%
  ggplot() + geom_point(aes(x=PCo1,y=PCo2,col=region),size=1.5) +
             theme_bw() + scale_color_viridis_d() +
             xlab("PCo1 (41.76 %)") + ylab("PCo2 (14.23 %)") +
             theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
             geom_hline(yintercept=0,linetype="dashed") + geom_vline(xintercept=0,linetype="dashed") 

p.sel.tot.region %>% ggsave("PCoA.sel.tot.region.png",.,device='png',width=15,height=10,units='cm')

p.sel.tot.country <- df.plot.sel.tot %>%
  ggplot() + geom_point(aes(x=PCo1,y=PCo2,col=country),size=1.5) +
             theme_bw() + scale_color_viridis_d(option='magma') +
             xlab("PCo1 (41.76 %)") + ylab("PCo2 (14.23 %)") +
             theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
             geom_hline(yintercept=0,linetype="dashed") + geom_vline(xintercept=0,linetype="dashed") 

p.sel.tot.country %>% ggsave("PCoA.sel.tot.country.png",.,device='png',width=15,height=10,units='cm')

Fst.sel <- as.data.frame(Fst.sel)
Fst.sel.atl <- Fst.sel[!(row.names(Fst.sel) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]
Fst.sel.atl <- as.matrix(Fst.sel.atl %>% select(-Laz,-Tar,-Sar,-Ale,-The,-Tor,-Sky))

result.sel.atl <- PCoA(Fst.sel.atl)
Evalues.sel.atl <- result.sel.atl$values
Evectors.sel.atl <- result.sel.atl$vectors

print(Evalues.sel.atl)

Evectors.sel.atl <- Evectors.sel.atl[,-c(30:31)]
Evalues.sel.atl <- Evalues.sel.atl[-c(30:31)]
for(i in 1:ncol(Evectors.sel.atl)){
  Evectors.sel.atl[,i] <- Evectors.sel.atl[,i]*Evalues.sel.atl[i]
}
print(Evectors.sel.atl)

df.plot.sel.atl <- Evectors.sel.atl[,1:2] %>%
               as.data.frame()
rownames(df.plot.sel.atl) <- rownames(Fst.sel.atl)
print(df.plot.sel.atl)

countryList <- c("NOR","GRB","IRL","GRB","GRB","NOR","SVE","GRB","DEU","IRL","FRA","GRB","GRB","CHA","SVE","IRL","GRB","GRB","SVE","IRL","NDL","GRB","GRB","GRB","GRB","GRB","NOR","GRB","NOR","IRL","ESP")
df.plot.sel.atl$country <- countryList
regionList <- c("ATL",rep("ATL",4),"SKA","SKA",rep("ATL",7),"SKA","ATL","ATL","ATL","SKA",rep("ATL",5),"ATL","ATL","SKA","ATL","ATL","ATL","ATL")
df.plot.sel.atl$region <- regionList

print(df.plot.sel.atl)

lambda.sel.atl <- Evalues.sel.atl[1:2] / sum(Evalues.sel.atl)

print(lambda.sel.atl)

p.sel.atl.region <- df.plot.sel.atl %>%
  ggplot() + geom_point(aes(x=PCo1,y=PCo2,col=region),size=1.5) +
             theme_bw() + scale_color_viridis_d() +
             xlab("PCo1 (32.74 %)") + ylab("PCo2 (3.85 %)") +
             theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
             geom_hline(yintercept=0,linetype="dashed") + geom_vline(xintercept=0,linetype="dashed") 

p.sel.atl.region %>% ggsave("PCoA.sel.atl.region.png",.,device='png',width=15,height=10,units='cm')

p.sel.atl.country <- df.plot.sel.atl %>%
  ggplot() + geom_point(aes(x=PCo1,y=PCo2,col=country),size=1.5) +
             theme_bw() + scale_color_viridis_d(option='magma') +
             xlab("PCo1 (32.74 %)") + ylab("PCo2 (3.85 %)") +
             theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
             geom_hline(yintercept=0,linetype="dashed") + geom_vline(xintercept=0,linetype="dashed") 

p.sel.atl.country %>% ggsave("PCoA.sel.atl.country.png",.,device='png',width=15,height=10,units='cm')