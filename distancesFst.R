library(adegenet)
library(diveRsity)

calcFst <- function(gendata){

	diffstats <- diffCalc(gendata, fst=T, pairwise=T, bs_locus=F, bs_pairwise=F)

	Fst <- diffstats$pairwise$Fst

	Fst.dist <- as.dist(Fst, diag=T, upper=T)

	Fst.mat <- as.matrix(Fst.dist)
	rownames(Fst.mat) <- substr(rownames(Fst.mat),1,3)
	colnames(Fst.mat) <- substr(colnames(Fst.mat),1,3)
	Fst.mat[Fst.mat < 0] <- 0

	return(Fst.mat)

}

Fst.tot <- calcFst("Input/lobster_1278ind_79snps_40pop.gen")
write.csv(Fst.tot,"Output/distancesFst_tot.csv")

Fst.sel <- calcFst("Input/lobster_1278ind_8snps_40pop.gen")
write.csv(Fst.sel,"Output/distancesFst_sel.csv")

Fst.neut <- calcFst("Input/lobster_1278ind_71snps_40pop.gen")
write.csv(Fst.neut,"Output/distancesFst_neut.csv")
