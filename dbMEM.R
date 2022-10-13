library(adespatial)

distances <- read.csv("Output/distances.csv",header=T,row.names=1)

distances <- distances[!(row.names(distances) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]
distances <- distances[,!(names(distances) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky"))]

distances_dist <- as.dist(t(distances))

dbMEM <- adespatial::dbmem(distances_dist, MEM.autocor="positive", store.listw=T, silent=F)
write.csv(dbMEM,"Output/dbMEM.csv")