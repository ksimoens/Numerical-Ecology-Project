library(tidyverse)
library(adespatial)

# https://rpubs.com/lbenestan/667437

con_mat <- read.csv("Input/AssymDist.csv", header=T, row.names=1)

con_mat <- as.matrix(con_mat)

melted_con <- melt(t(con_mat))
colnames(melted_con) <- c("SITE1","SITE2","Dispersal")

p <- ggplot(data = melted_con, aes(SITE1, SITE2, fill = factor(Dispersal)))+
  geom_tile(color = "black")+
  scale_fill_manual(values=c("white","#f8dcdc","#f89696","red"), 
                       name="Dispersal probabilities") +
  ylab("Sampling location A")+
  xlab("Sampling location B")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 8, hjust = 1))+ 
  theme(legend.position="top",legend.text = element_text(color = "black", size = 6)) +
  coord_fixed()

upper <- con_mat[upper.tri(con_mat, diag = TRUE)]
lower <- con_mat[lower.tri(con_mat, diag = TRUE)]

highest <- vector()
for (i in 1 : length(upper)) {
  if (upper[i] < lower[i]) {
    highest[i] <- lower[i]
  } else {
    highest[i] <- upper[i]
  }
}

con_mat2 <- con_mat
con_mat2[upper.tri(con_mat2, diag = TRUE)] <- highest

con_mat2[lower.tri(con_mat2, diag = T)] <- 0
diag(con_mat2) <- diag(con_mat)

con_mat_ex <- expand.grid(con_mat2)

con_mat_ex$SITE1 <- rep(colnames(con_mat2), each = 31)
con_mat_ex$SITE2 <- rep(colnames(con_mat2), 31)
con_mat_ex$POP1id <- rep(seq(1,31,1), each = 31)
con_mat_ex$POP2id <- rep(seq(1,31,1), 31)

con_mat_ex$PAIR <- paste(con_mat_ex$SITE1, con_mat_ex$SITE2, sep = "-")
con_mat_ex <- con_mat_ex[,-c(2,3)]

colnames(con_mat_ex) <- c("PROB", "SITE1", "SITE2", "PAIR")

con.mat.highest.prob.no0 <- dplyr::filter(con_mat_ex, PROB > 0)
print(nrow(con.mat.highest.prob.no0))

edges.highest.prob <- con.mat.highest.prob.no0[,c(2,3)]
edges <- edges.highest.prob

edges <- dplyr::filter(edges, SITE1 != SITE2)

xy <- read.csv("Input/coordinates.csv",header=T, row.names=1)
xy <- xy[!(row.names(xy) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]
xy$Nb <- 1:31
xy <- xy[,c(3,1,2)]

bin.mat <- aem.build.binary(coords=xy,link=edges, plot.connexions = TRUE)

weights <- dplyr::filter(con.mat.highest.prob.no0, SITE1 != SITE2)
weight.vec <- as.vector(weights[,1])

cuke.aem.wt <- adespatial::aem(aem.build.binary = bin.mat, weight = weight.vec, rm.link0 = TRUE)

AEM.vectors.wt <- as.data.frame(cuke.aem.wt$vectors)
rownames(AEM.vectors.wt) <- rownames(xy)
names(AEM.vectors.wt) <- paste0("AEM",1:ncol(AEM.vectors.wt))
write.csv(AEM.vectors.wt,"Output/AEM.csv")
