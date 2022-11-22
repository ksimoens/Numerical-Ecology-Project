library(adespatial)

# Read in-water distances
distances <- read.csv("Output/distances.csv",header=T,row.names=1)

# Select Atlantic sites
distances <- distances[!(row.names(distances) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]
distances <- distances[,!(names(distances) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky"))]

# Transform to distance object
distances_dist <- as.dist(t(distances))

# Build dbMEMs, only positive eigenvalues
dbMEM <- adespatial::dbmem(distances_dist, MEM.autocor="positive", store.listw=T, silent=F)
write.csv(dbMEM,"Output/dbMEM.csv")