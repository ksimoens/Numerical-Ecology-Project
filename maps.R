library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(sf)
library(ggspatial)
library(cowplot)
library(gdistance)
library(geosphere)
library(maptools)

world <- ne_countries(scale = "large", returnclass = "sf")

sf_use_s2(FALSE)

map <- ggplot(data = world) +
    geom_sf() +
    coord_sf(xlim = c(-15,15), ylim = c(40,65), expand = FALSE) + 
    theme_bw() +
    theme(panel.background = element_rect(fill = "#d0dfff"),
    	  panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "grey")) +
    xlab("longitude") + ylab("latitude") +
    annotation_scale(location = "tl", width_hint = 0.5) + 
    annotation_north_arrow(location="tl", which_north = "true", 
        pad_x = unit(0.3, "cm"), pad_y = unit(1, "cm"),
        style = north_arrow_fancy_orienteering)

inset <- ggplot() +
	geom_sf(data=world,fill="grey") +
    coord_sf(xlim = c(-17,42), ylim = c(33,72), expand = FALSE) + 
    geom_rect(aes(xmin=-15,xmax=15,ymin=40,ymax=65),col='red',fill=rgb(1,0,0,0.25),size=2) +
    theme_void() + 
    theme(panel.background = element_rect(fill = "white") )

map_inset <- ggdraw() +
  draw_plot(map) +
  draw_plot(inset, x = 0.577, y = -0.273,width=0.3)

coord <- read.csv("Input/coordinates.csv",header=T,row.names=1)
coord <- coord[!(row.names(coord) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]

dbMEM <- read.csv("Output/dbMEM.csv",header=T,row.names=1)
dbMEM <- dbMEM[match(rownames(coord),rownames(dbMEM)),]
dbMEM.sel <- dbMEM[,c(1,3,4,5,7)]

dbMEM.coord <- cbind(coord,dbMEM.sel)

mapMEM1 <- map + geom_point(data=dbMEM.coord,aes(x=lon,y=lat,col=MEM1),size=3) + scale_color_viridis_c(option='magma')
mapMEM1 <- ggdraw() +
  draw_plot(mapMEM1) +
  draw_plot(inset, x = 0.504, y = -0.261,width=0.3)

ggsave("mapMEM1.png",mapMEM1,device="png",width=15,height=15,units='cm')

mapMEM3 <- map + geom_point(data=dbMEM.coord,aes(x=lon,y=lat,col=MEM3),size=3) + scale_color_viridis_c(option='magma')
mapMEM3 <- ggdraw() +
  draw_plot(mapMEM3) +
  draw_plot(inset, x = 0.513, y = -0.261,width=0.3)

ggsave("mapMEM3.png",mapMEM3,device="png",width=15,height=15,units='cm')

mapMEM3 <- map + geom_point(data=dbMEM.coord,aes(x=lon,y=lat,col=MEM3),size=3) + scale_color_viridis_c(option='magma')
mapMEM3 <- ggdraw() +
  draw_plot(mapMEM3) +
  draw_plot(inset, x = 0.513, y = -0.261,width=0.3)

ggsave("mapMEM3.png",mapMEM3,device="png",width=15,height=15,units='cm')

mapMEM4 <- map + geom_point(data=dbMEM.coord,aes(x=lon,y=lat,col=MEM4),size=3) + scale_color_viridis_c(option='magma')
mapMEM4 <- ggdraw() +
  draw_plot(mapMEM4) +
  draw_plot(inset, x = 0.513, y = -0.261,width=0.3)

ggsave("mapMEM4.png",mapMEM4,device="png",width=15,height=15,units='cm')

mapMEM5 <- map + geom_point(data=dbMEM.coord,aes(x=lon,y=lat,col=MEM5),size=3) + scale_color_viridis_c(option='magma')
mapMEM5 <- ggdraw() +
  draw_plot(mapMEM5) +
  draw_plot(inset, x = 0.513, y = -0.261,width=0.3)

ggsave("mapMEM5.png",mapMEM5,device="png",width=15,height=15,units='cm')

mapMEM7 <- map + geom_point(data=dbMEM.coord,aes(x=lon,y=lat,col=MEM7),size=3) + scale_color_viridis_c(option='magma')
mapMEM7 <- ggdraw() +
  draw_plot(mapMEM7) +
  draw_plot(inset, x = 0.513, y = -0.261,width=0.3)

ggsave("mapMEM7.png",mapMEM7,device="png",width=15,height=15,units='cm')

env <- read.csv("Output/EnvMatrix.csv",header=T,row.names=1)
env <- env[!(row.names(coord) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]
env <- env[match(rownames(coord),rownames(env)),]

env.sel <- env[,c(1,2,3)]

env.coord <- cbind(coord,env.sel)

mapSSTMean <- map + geom_point(data=env.coord,aes(x=lon,y=lat,col=SST_mean),size=3) + scale_color_viridis_c(option='magma') +
              guides(col=guide_colourbar(title="SST mean (°C)"))
mapSSTMean <- ggdraw() +
  draw_plot(mapSSTMean) +
  draw_plot(inset, x = 0.460, y = -0.255,width=0.3)

ggsave("mapSSTMean.png",mapSSTMean,device="png",width=15,height=15,units='cm')

mapSSTRange <- map + geom_point(data=env.coord,aes(x=lon,y=lat,col=SST_range),size=3) + scale_color_viridis_c(option='magma') +
              guides(col=guide_colourbar(title="SST range (°C)"))
mapSSTRange <- ggdraw() +
  draw_plot(mapSSTRange) +
  draw_plot(inset, x = 0.459, y = -0.255,width=0.3)

ggsave("mapSSTRange.png",mapSSTRange,device="png",width=15,height=15,units='cm')

mapSALMean <- map + geom_point(data=env.coord,aes(x=lon,y=lat,col=SAL_mean),size=3) + scale_color_viridis_c(option='magma') +
              guides(col=guide_colourbar(title="mean salinity"))
mapSALMean <- ggdraw() +
  draw_plot(mapSALMean) +
  draw_plot(inset, x = 0.4765, y = -0.261,width=0.3)

#ggsave("mapSALMean.png",mapSALMean,device="png",width=15,height=15,units='cm')

lin <- read.csv("Output/PCoSpatial.csv",header=T,row.names=1)[,1:2]
lin <- lin[match(rownames(coord),rownames(env)),]

lin.coord <- cbind(coord,lin)

mapLin1 <- map + geom_point(data=lin.coord,aes(x=lon,y=lat,col=PCo1),size=3) + scale_color_viridis_c(option='magma') +
              guides(col=guide_colourbar(title="spatial PCo1"))
mapLin1 <- ggdraw() +
  draw_plot(mapLin1) +
  draw_plot(inset, x = 0.477, y = -0.261,width=0.3)

#ggsave("mapLin1.png",mapLin1,device="png",width=15,height=15,units='cm')

mapLin2 <- map + geom_point(data=lin.coord,aes(x=lon,y=lat,col=PCo2),size=3) + scale_color_viridis_c(option='magma') +
              guides(col=guide_colourbar(title="spatial PCo2"))
mapLin2 <- ggdraw() +
  draw_plot(mapLin2) +
  draw_plot(inset, x = 0.478, y = -0.261,width=0.3)

#ggsave("mapLin2.png",mapLin2,device="png",width=15,height=15,units='cm')

AEM <- read.csv("Output/AEM.csv",header=T,row.names=1)
AEM <- AEM[match(rownames(coord),rownames(AEM)),]
AEM.sel <- AEM[,c(1,10,12,15,17,18,19,22,27)]

AEM.coord <- cbind(coord,AEM)

mapAEM1 <- map + geom_point(data=AEM.coord,aes(x=lon,y=lat,col=AEM1),size=3) + scale_color_viridis_c(option='magma')
mapAEM1 <- ggdraw() +
  draw_plot(mapAEM1) +
  draw_plot(inset, x = 0.5107, y = -0.261,width=0.3)

#ggsave("mapAEM1.png",mapAEM1,device="png",width=15,height=15,units='cm')

mapAEM10 <- map + geom_point(data=AEM.coord,aes(x=lon,y=lat,col=AEM10),size=3) + scale_color_viridis_c(option='magma')
mapAEM10 <- ggdraw() +
  draw_plot(mapAEM10) +
  draw_plot(inset, x = 0.505, y = -0.261,width=0.3)

#ggsave("mapAEM10.png",mapAEM10,device="png",width=15,height=15,units='cm')

mapAEM12 <- map + geom_point(data=AEM.coord,aes(x=lon,y=lat,col=AEM12),size=3) + scale_color_viridis_c(option='magma')
mapAEM12 <- ggdraw() +
  draw_plot(mapAEM12) +
  draw_plot(inset, x = 0.505, y = -0.261,width=0.3)

#ggsave("mapAEM12.png",mapAEM12,device="png",width=15,height=15,units='cm')

mapAEM15 <- map + geom_point(data=AEM.coord,aes(x=lon,y=lat,col=AEM15),size=3) + scale_color_viridis_c(option='magma')
mapAEM15 <- ggdraw() +
  draw_plot(mapAEM15) +
  draw_plot(inset, x = 0.505, y = -0.261,width=0.3)

#ggsave("mapAEM15.png",mapAEM15,device="png",width=15,height=15,units='cm')

mapAEM17 <- map + geom_point(data=AEM.coord,aes(x=lon,y=lat,col=AEM17),size=3) + scale_color_viridis_c(option='magma')
mapAEM17 <- ggdraw() +
  draw_plot(mapAEM17) +
  draw_plot(inset, x = 0.505, y = -0.261,width=0.3)

#ggsave("mapAEM17.png",mapAEM17,device="png",width=15,height=15,units='cm')

mapAEM18 <- map + geom_point(data=AEM.coord,aes(x=lon,y=lat,col=AEM18),size=3) + scale_color_viridis_c(option='magma')
mapAEM18 <- ggdraw() +
  draw_plot(mapAEM18) +
  draw_plot(inset, x = 0.498, y = -0.261,width=0.3)

#ggsave("mapAEM18.png",mapAEM18,device="png",width=15,height=15,units='cm')

mapAEM19 <- map + geom_point(data=AEM.coord,aes(x=lon,y=lat,col=AEM19),size=3) + scale_color_viridis_c(option='magma')
mapAEM19 <- ggdraw() +
  draw_plot(mapAEM19) +
  draw_plot(inset, x = 0.505, y = -0.261,width=0.3)

#ggsave("mapAEM19.png",mapAEM19,device="png",width=15,height=15,units='cm')

mapAEM22 <- map + geom_point(data=AEM.coord,aes(x=lon,y=lat,col=AEM22),size=3) + scale_color_viridis_c(option='magma')
mapAEM22 <- ggdraw() +
  draw_plot(mapAEM22) +
  draw_plot(inset, x = 0.505, y = -0.261,width=0.3)

#ggsave("mapAEM22.png",mapAEM22,device="png",width=15,height=15,units='cm')

mapAEM27 <- map + geom_point(data=AEM.coord,aes(x=lon,y=lat,col=AEM27),size=3) + scale_color_viridis_c(option='magma')
mapAEM27 <- ggdraw() +
  draw_plot(mapAEM27) +
  draw_plot(inset, x = 0.505, y = -0.261,width=0.3)

#ggsave("mapAEM27.png",mapAEM27,device="png",width=15,height=15,units='cm')

coord <- read.csv("Input/coordinates.csv",header=T,row.names=1)

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

coord$country <- countryList(rownames(coord))

data(wrld_simpl)
shp <- wrld_simpl
shp <- crop(shp, extent(-15.788,31.0281,33.227,71.0976))
r <- raster(nrow=1000, ncol=1000)
extent(r) <- extent(-15.788,31.0281,33.227,71.0976)
r <- rasterize(shp,r,output='text')
r[is.na(r)] <- -999
r[r>-999] <- NA
r[r==-999] <- 1
coord[1,1] <- 0
coord[7,2] <- 50.34
coord[10,1] <- -5.03
coord[10,2] <- 50.58
coord[18,2] <- 40.80
coord[19,1] <- 24.52
coord[21,1] <- 23.53
coord[21,2] <- 40.15
coord[22,1] <- -8.25
coord[28,1] <- 11.59
coord[28,2] <- 42.16
trans <- transition(r, mean, directions=16)
trans <- geoCorrection(trans,'c')
Tro <- as.numeric(coord[rownames(coord)=="Tro",1:2])
Sky <- as.numeric(coord[rownames(coord)=="Sky",1:2])
TrotoSky <- shortestPath(trans, Tro, Sky, output="SpatialLines")
TrotoSky <- SpatialLinesDataFrame(TrotoSky, data = data.frame(ID = 1))

getColours <- function(pop){
  col_key <- data.frame(
      population = c("CHA","DEU","ESP","FRA","GRB","HEL","IRL","ITA","NDL","NOR","SVE","ATL","MED","SKA"),
      colour = c(rep(c("#C7E020FF","#24868EFF","#440154FF"),3),c("#C7E020FF","#440154FF"),c("#C7E020FF","#24868EFF","#440154FF")),
      shape = c(rep(15,3),rep(17,3),rep(19,3),18,18,15,17,19))

  col_res <- col_key[col_key$population %in% unique(pop),]$colour
  sha_res <- col_key[col_key$population %in% unique(pop),]$shape
  return(list(colour = col_res,shape = sha_res))
}

lines <- ggplot() +
	geom_sf(data=world) +
    coord_sf(xlim = c(-17,42), ylim = c(33,72), expand = FALSE) + 
    theme_bw() +
    theme(panel.background = element_rect(fill = "#d0dfff"),
    	  panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "grey")) +
    xlab("longitude") + ylab("latitude") +
    annotation_scale(location = "tl", width_hint = 0.5) + 
    annotation_north_arrow(location="tl", which_north = "true", 
        pad_x = unit(0.5, "cm"), pad_y = unit(1, "cm"),
        style = north_arrow_fancy_orienteering) +
    geom_point(data=coord,aes(x=lon,y=lat,col=country,shape=country),size=3) + 
    scale_colour_manual(
              name = "country",
              labels = sort(unique(coord$country)),
              values=getColours(coord$country)$colour) +
             scale_shape_manual(
              name = "country",
              labels =sort(unique(coord$country)),
              values=getColours(coord$country)$shape) + 
    geom_path(data=TrotoSky,aes(x=long,y=lat),size=1.5)

ggsave("mapLine.png",lines,device='png',width=15,height=15,units='cm')