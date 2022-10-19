library(maptools)
library(gdistance)
library(geosphere)

# load the shapefile data from maptools
data(wrld_simpl)

# load the coordinates of the sampling sites: 2 columns = lon and lat
coord <- read.csv("Input/coordinates.csv",header=T, row.names=1)

# crop the world to my problem: Europe
shp <- wrld_simpl
shp <- crop(shp, extent(-15.788,31.0281,33.227,71.0976))

# make a dummy raster with a number of pixels:
# larger number of pixels -> better resolution BUT takes a lot more from the computer
# what number? -> trial and error (start with this one)
r <- raster(nrow=1000, ncol=1000)
# set geographic extent of dummy raster
extent(r) <- extent(-15.788,31.0281,33.227,71.0976)
# rasterise the map shapefile using the dummy raster
r <- rasterize(shp,r,output='text')

# penalise going over land
# first set pixel values in water to -999
r[is.na(r)] <- -999
# secondly set pixel values on land to NA
r[r>-999] <- NA
# thirdly reset pixel values in water to 1
r[r==-999] <- 1

# cheat: move sites that fall on land pixels to water pixels
# you will have to check this for your own sites. This will happen due to resolution issues.
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

#initialise distance matrix (upper, no diagonal)
dists <- matrix(0,nrow=nrow(coord),ncol=nrow(coord))
rownames(dists) <- rownames(coord)
colnames(dists) <- rownames(coord)

# make transition matrix (from gdistance)
# depending on the number of pixels in the raster, this can take a while and will take memory
# for each pixel the connection to all the other pixels is calculated
# directions parameter fixes the possible connections you can have to neighbouring pixels
trans <- transition(r, mean, directions=16)
# make some correction
trans <- geoCorrection(trans,'c')

# for each row i
for(i in 1:nrow(dists)){
	# get the coordinate of the site at row i
	A <- as.numeric(coord[i,])
	# for each column j starting from i+1 (only upper and no diagonal)
	for(j in (i+1):nrow(dists) ){
		# get the coordinate of the site at column j
		B <- as.numeric(coord[j,])
		# calculate the shortest path between coord i and coord j: output = spatialLine
		AtoB_line <- shortestPath(trans, A, B, output="SpatialLines")
		# calculate the length of the line (in kilometers) -> will be the geographic distance (geosphere)
		AtoB_dist <- lengthLine(AtoB_line)/1000
		# make some output to follow the progress (can be slow)
		out <- paste(rownames(coord[i,]), rownames(coord[j,]), AtoB_dist, sep="  ")
		print(out)
		# save the output to the distance matrix
		dists[i,j] <- AtoB_dist
	}
}

# print the results
print(dists)
# save the results
write.csv(dists,"distances.csv", quote=F, row.names=T, col.names=T)
