library(raster) # raster library: necessity for GIS in R
				# it can do almost anything

# read the coordinates from a csv-file
# use read.table if you have a txt-file
# format: site		lon		lat 
#		  Brd      -0.17   54.07			
#           
# longitude and latitude in WGS 84
# row.names = 1 will take the site column as row names
# header = T will take the first row as column names
coord <- read.csv('Input/coordinates.csv',row.names=1,header=T)

# I have .tif files for each environmental variable.
# e.g. SST_mean.tif is a raster file which contains the mean SST for each point on the planet.
# use raster("filename") to read raster files (.tif e.g.)
# I place all the rasters in a list to extract everything in one go. You can do one by one. 
rasterList <- list( SST_mean = raster('Input/EnvData/SST_mean.tif'),
					SST_range = raster('Input/EnvData/SST_range.tif'),
					SAL_mean = raster('Input/EnvData/SAL_mean.tif'),
					SAL_range = raster('Input/EnvData/SAL_range.tif'),
					VEL_mean = raster('Input/EnvData/VEL_mean.tif'),
					PP_mean = raster('Input/EnvData/PP_mean.tif'),
					PP_range = raster('Input/EnvData/PP_range.tif'),
					bathy = raster('Input/EnvData/bathymetry.tif')
					)

# res <- extract(ras,coord) is the function you are after:
# it takes the value from the raster (ras) - e.g. mean SST - at each coordinate in the dataframe (coord)
# the result (res) is a vector of the same length as nrow in coord with the environmental variable at each coordinate
# the lapply function does this for each raster in the list and adds it as columns to a final matrix.
env <- lapply(rasterList, function(ras){res <- extract(ras,coord); return(res)})
# transform to dataframe
env <- as.data.frame(env)
# rownames = site names (from the coordinate file)
rownames(env) <- rownames(coord)
# cheating a bit:
# the bathymetry at site 10 is positive (+1). This can happen due to resolution issues.
# I just change it to -1.
env$bathy[10] <- -1*env$bathy[10]
print(env)
# write the environmental table to a file.
# write.table for .txt; write.csv for .csv...
# col,names = row.names = T: to get the names in the file
# this will be matrix X in the RDA.
write.table(env, file='Output/EnvMatrix.txt', quote=F, col.names = T, row.names = T)
