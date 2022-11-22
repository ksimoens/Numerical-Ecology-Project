# PCoA with the in-water distances (see PCoA.R and functions.R for details)
# load in-water distances
distances <- read.csv("Output/distances.csv", header=T, row.names=1)

# only for Atlantic sites
dist_mat <- as.matrix(as.dist(t(distances), upper=T, diag=T))
dist_mat <- dist_mat[!(row.names(dist_mat) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]
dist_mat <- dist_mat[,!(colnames(dist_mat) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky"))]

A <- -0.5 * dist_mat * dist_mat

ones <- rep(1,ncol(dist_mat))
identity <- diag(ncol(dist_mat))
side <- identity - ones %*% t(ones) / ncol(dist_mat)

Delta <- side %*% A %*% side

Evalues <- eigen(Delta)$values

c <- abs(Evalues[length(Evalues)])

c.mat <- matrix(rep(c,ncol(A)*ncol(A)),ncol=ncol(A))
c.mat <- c.mat - c*diag(ncol(A))

dist.corr <- sqrt(dist_mat * dist_mat + 2*c.mat)
A <- -0.5 * dist.corr * dist.corr
Delta <- side %*% A %*% side

Evalues <- eigen(Delta)$values
Evectors <- eigen(Delta)$vectors

print(Evalues)

Evectors <- Evectors[,-c(30,31)]
Evalues <- Evalues[-c(30,31)]
for(i in 1:ncol(Evectors)){
	Evectors[,i] <- Evectors[,i] * Evalues[i]
}
colnames(Evectors) <- paste0("PCo",1:ncol(Evectors))
rownames(Evectors) <- rownames(dist_mat)

# output the linear PCo for use in RDA
write.csv(Evectors,"Output/PCoSpatial.csv")