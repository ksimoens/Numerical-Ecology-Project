library(adegenet)
library(diveRsity)

data <- read.genepop("Input/lobster_1278ind_79snps_40pop.gen")

diffstats <- diffCalc("Input/lobster_1278ind_79snps_40pop.gen", fst=T, pairwise=T, bs_locus=F, bs_pairwise=F)

Fst <- diffstats$pairwise$Fst

Fst.dist <- as.dist(Fst, diag=T, upper=T)

Fst.mat <- as.matrix(Fst.dist)
rownames(Fst.mat) <- substr(rownames(Fst.mat),1,3)
colnames(Fst.mat) <- substr(colnames(Fst.mat),1,3)
Fst.mat[Fst.mat < 0] <- 0

write.csv(Fst.mat,"Output/distancesFst.csv")
