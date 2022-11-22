library(tidyverse)
library(adespatial)

# Based on https://rpubs.com/lbenestan/667437

# Read connectivity matrix (based on ocean currents)
con_mat <- read.csv("Input/AssymDist.csv", header=T, row.names=1)

con_mat <- as.matrix(con_mat)

# Transform to long format
melted_con <- melt(t(con_mat))
colnames(melted_con) <- c("SITE1","SITE2","Dispersal")

# Plot the connectivity matrix
p <- ggplot(data = melted_con, aes(SITE1, SITE2, fill = factor(Dispersal)))+
  geom_tile(color = "black")+
  scale_fill_manual(values=c("white","#f8dcdc","#f89696","red"), 
                       name="Dispersal probabilities") +
  ylab("Source location")+
  xlab("Destination location")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 8, hjust = 1))+ 
  theme(legend.position="top",legend.text = element_text(color = "black", size = 6)) +
  coord_fixed()

upper <- con_mat[upper.tri(con_mat, diag = TRUE)]
lower <- con_mat[lower.tri(con_mat, diag = TRUE)]

# Find highest probability of two between two sites
highest <- vector()
for (i in 1 : length(upper)) {
  if (upper[i] < lower[i]) {
    highest[i] <- lower[i]
  } else {
    highest[i] <- upper[i]
  }
}

# Replace the upper by the highest value
con_mat2 <- con_mat
con_mat2[upper.tri(con_mat2, diag = TRUE)] <- highest

# Set the lower to zero
con_mat2[lower.tri(con_mat2, diag = T)] <- 0
# Diagonal values remain the same
diag(con_mat2) <- diag(con_mat)

# Make dataframe with all pairs of sites with highest values
con_mat_ex <- expand.grid(con_mat2)

# Make unique identifiers
con_mat_ex$SITE1 <- rep(colnames(con_mat2), each = 31)
con_mat_ex$SITE2 <- rep(colnames(con_mat2), 31)
con_mat_ex$POP1id <- rep(seq(1,31,1), each = 31)
con_mat_ex$POP2id <- rep(seq(1,31,1), 31)

# Make pair of sites
con_mat_ex$PAIR <- paste(con_mat_ex$SITE1, con_mat_ex$SITE2, sep = "-")
# Remove the single sites
con_mat_ex <- con_mat_ex[,-c(2,3)]

colnames(con_mat_ex) <- c("PROB", "SITE1", "SITE2", "PAIR")

# Remove the pairs with no connection
con.mat.highest.prob.no0 <- dplyr::filter(con_mat_ex, PROB > 0)
print(nrow(con.mat.highest.prob.no0))

# Create matrix of identifiers for edges between pairs of sites
edges.highest.prob <- con.mat.highest.prob.no0[,c(2,3)]
edges <- edges.highest.prob

# No self-connections
edges <- dplyr::filter(edges, SITE1 != SITE2)

# Coordinate information 
xy <- read.csv("Input/coordinates.csv",header=T, row.names=1)
xy <- xy[!(row.names(xy) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]
xy$Nb <- 1:31
xy <- xy[,c(3,1,2)]

# Build binary site-by-edge matrix and plot
bin.mat <- aem.build.binary(coords=xy,link=edges, plot.connexions = TRUE)

# Weigh the edges by the connection probability
weights <- dplyr::filter(con.mat.highest.prob.no0, SITE1 != SITE2)
weight.vec <- as.vector(weights[,1])

# Create the AEM vectors
# rm.link0: remove the artificial upstream site
cuke.aem.wt <- adespatial::aem(aem.build.binary = bin.mat, weight = weight.vec, rm.link0 = TRUE)

# Clean up and export
AEM.vectors.wt <- as.data.frame(cuke.aem.wt$vectors)
rownames(AEM.vectors.wt) <- rownames(xy)
names(AEM.vectors.wt) <- paste0("AEM",1:ncol(AEM.vectors.wt))
write.csv(AEM.vectors.wt,"Output/AEM.csv")
