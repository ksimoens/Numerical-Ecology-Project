library(ade4)
library(tidyverse)

Fst.dist <- read.csv("Output/distancesFst_neut.csv",header=T, row.names=1)

inwater.dist <- read.csv("Output/distances.csv",header=T, row.names=1)
Fst.dist <- Fst.dist[match(rownames(inwater.dist),rownames(Fst.dist)),]
Fst.dist <- Fst.dist[,match(names(inwater.dist),names(Fst.dist))]
print(Fst.dist)

df.plot <- expand.grid(as.matrix(Fst.dist))

inwater.dist <- as.matrix(as.dist(t(inwater.dist),upper=T, diag=T))
inwater.vec <- expand.grid(inwater.dist)$Var1
df.plot$inwater <- inwater.vec


df.plot$SITE1 <- rep(colnames(Fst.dist), each = ncol(Fst.dist))
df.plot$SITE2 <- rep(colnames(Fst.dist), ncol(Fst.dist))

index <- vector()
for(i in 0:37){
	for(j in 1:(i+1)){
		index <- append(index, i*38+j)
	}
}

df.plot <- df.plot[-index,]
rownames(df.plot) <- 1:nrow(df.plot)

df.plot$REGION1 <- rep("",nrow(df.plot))
df.plot[df.plot$SITE1 %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky"),]$REGION1 <- "MED"
df.plot[!(df.plot$SITE1 %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]$REGION1 <- "ATL"
df.plot[df.plot$SITE1 %in% c("Oos"),]$REGION1 <- "OOS"
df.plot$REGION2 <- rep("",nrow(df.plot))
df.plot[df.plot$SITE2 %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky"),]$REGION2 <- "MED"
df.plot[!(df.plot$SITE2 %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]$REGION2 <- "ATL"
df.plot[df.plot$SITE2 %in% c("Oos"),]$REGION2 <- "OOS"

df.plot$PAIR <- paste(df.plot$REGION1, df.plot$REGION2, sep = "-")

df.plot[df.plot$PAIR == "MED-ATL",]$PAIR <- "ATL-MED"
df.plot[df.plot$PAIR == "OOS-ATL",]$PAIR <- "ATL-OOS"

names(df.plot)[1] <- "Fst"

print(summary(lm(data=df.plot, Fst~inwater )))

Fst.mantel <- as.dist(Fst.dist)
inwater.mantel <- as.dist(inwater.dist)
print(ade4::mantel.rtest(Fst.mantel, inwater.mantel, nrepet=9999))

p <- ggplot(data=df.plot,aes(x=inwater,y=Fst,col=PAIR)) + geom_point(size=1.5) +
	theme_bw() + scale_color_viridis_d(option='magma') +
	xlab("in-water distance (km)") + ylab(expression(F[st])) + 
	geom_abline(intercept = -5.470e-3, slope = 1.578e-5, color="red", linetype="dashed", size=1) + 
	geom_label(x=5400,y=0.0125,label=expression(paste(R^2 == 0.8781, " (p < 0.0001)")),show.legend=F,size=4)

p %>% ggsave("Mantel_total.png",.,device='png',width=15,height=10,units='cm')

df.plot.atl <- df.plot[df.plot$REGION1!="MED" ,]
df.plot.atl <- df.plot.atl[df.plot.atl$REGION2!="MED",]
rownames(df.plot.atl) <- 1:nrow(df.plot.atl)

Fst.dist.atl <- Fst.dist[!(rownames(Fst.dist) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),
	!(colnames(Fst.dist) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky"))]

inwater.dist.atl <- inwater.dist[!(rownames(inwater.dist) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),
	!(colnames(inwater.dist) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky"))]

Fst.mantel.atl <- as.dist(Fst.dist.atl)
inwater.mantel.atl <- as.dist(inwater.dist.atl)
print(ade4::mantel.rtest(Fst.mantel.atl,inwater.mantel.atl,nrepet=9999))

print(summary(lm(data=df.plot.atl, Fst~inwater)))

p <- ggplot(data=df.plot.atl,aes(x=inwater,y=Fst,col=PAIR)) + geom_point(size=1.5) +
	theme_bw() + scale_color_viridis_d(option='magma') + xlim(0,3000) +
	xlab("in-water distance (km)") + ylab(expression(F[st])) + 
	geom_abline(intercept = 3.386e-3, slope = 3.576e-6, color="red", linetype="dashed", size=1) +
	geom_label(x=2300,y=0.045,label=expression(paste(R^2 == 0.1795, " (p = 0.06)")),show.legend=F,size=4)

p %>% ggsave("Mantel_atl.png",.,device='png',width=15,height=10,units='cm')

df.plot.oos <- df.plot.atl[df.plot.atl$SITE1 != "Oos",]
df.plot.oos <- df.plot.oos[df.plot.oos$SITE2 != "Oos",]
rownames(df.plot.oos) <- 1:nrow(df.plot.oos)

df.plot.oos$COUNTRY1 <- ""
df.plot.oos$COUNTRY2 <- ""

df.plot.oos[df.plot.oos$SITE1 %in% c("Brd","Cro","Eye","Heb","Iom","Ios","Loo","Lyn","Ork","Pad","Pem","She","Sbs","Sul","Cor","Hoo","Kil","Mul","Ven","Jer"),]$COUNTRY1 <- "BRT"
df.plot.oos[df.plot.oos$SITE2 %in% c("Brd","Cro","Eye","Heb","Iom","Ios","Loo","Lyn","Ork","Pad","Pem","She","Sbs","Sul","Cor","Hoo","Kil","Mul","Ven","Jer"),]$COUNTRY2 <- "BRT"
df.plot.oos[df.plot.oos$SITE1 %in% c("Idr","Vig"),]$COUNTRY1 <- "BIS"
df.plot.oos[df.plot.oos$SITE2 %in% c("Idr","Vig"),]$COUNTRY2 <- "BIS"
df.plot.oos[df.plot.oos$SITE1 %in% c("Hel","Ber","Flo","Sin","Tro","Gul","Kav","Lys"),]$COUNTRY1 <- "SKA"
df.plot.oos[df.plot.oos$SITE2 %in% c("Hel","Ber","Flo","Sin","Tro","Gul","Kav","Lys"),]$COUNTRY2 <- "SKA"

df.plot.oos$PAIR <- paste(df.plot.oos$COUNTRY1, df.plot.oos$COUNTRY2, sep = "-")
df.plot.oos[df.plot.oos$PAIR == "BIS-BRT",]$PAIR <- "BRT-BIS"
df.plot.oos[df.plot.oos$PAIR == "BIS-SKA",]$PAIR <- "SKA-BIS"
df.plot.oos[df.plot.oos$PAIR == "SKA-BRT",]$PAIR <- "BRT-SKA"

Fst.dist.oos <- Fst.dist.atl[!(rownames(Fst.dist.atl) %in% c("Oos")),
	!(colnames(Fst.dist.atl) %in% c("Oos"))]

inwater.dist.oos <- inwater.dist.atl[!(rownames(inwater.dist.atl) %in% c("Oos")),
	!(colnames(inwater.dist.atl) %in% c("Oos"))]

Fst.mantel.oos <- as.dist(Fst.dist.oos)
inwater.mantel.oos <- as.dist(inwater.dist.oos)
print(ade4::mantel.rtest(Fst.mantel.oos,inwater.mantel.oos,nrepet=9999))

print(summary(lm(data=df.plot.oos, Fst~inwater)))

p <- ggplot(data=df.plot.oos,aes(x=inwater,y=Fst)) + geom_point(size=1.5) +
	theme_bw() + 
	xlab("in-water distance (km)") + ylab(expression(F[st])) +
	geom_abline(intercept = -1.530e-3, slope = 4.619e-6, color="red", linetype="dashed", size=1) +
	geom_label(x=500,y=0.03,label=expression(paste(R^2 == 0.4497, " (p < 0.0001)")),show.legend=F,size=4)

p %>% ggsave("Mantel_oos.png",.,device='png',width=15,height=10,units='cm')
