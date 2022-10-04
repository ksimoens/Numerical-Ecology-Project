library(maptools)
library(gdistance)
library(geosphere)

data(wrld_simpl)

coord <- read.csv("Input/coordinates.csv",header=T, row.names=1)

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

dists <- matrix(0,nrow=nrow(coord),ncol=nrow(coord))
rownames(dists) <- rownames(coord)
colnames(dists) <- rownames(coord)

trans <- transition(r, mean, directions=16)
trans <- geoCorrection(trans,'c')

for(i in 1:nrow(dists)){
	A <- as.numeric(coord[i,])
	for(j in (i+1):nrow(dists) ){
		B <- as.numeric(coord[j,])
		AtoB_line <- shortestPath(trans, A, B, output="SpatialLines")
		AtoB_dist <- lengthLine(AtoB_line)/1000
		out <- paste(rownames(coord[i,]), rownames(coord[j,]), AtoB_dist, sep="  ")
		print(out)
		dists[i,j] <- AtoB_dist
	}
}

print(dists)
