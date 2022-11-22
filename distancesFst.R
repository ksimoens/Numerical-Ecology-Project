library(adegenet)
library(hierfstat)

# load necessary functions
source('functions.R')

# read gendata with all SNPs
dat.tot <- read.genepop("Input/lobster_1278ind_79snps_40pop.gen")
# transform to df to manipulate names
dat.df.tot <- genind2df(dat.tot)

# consistent naming + combine temporal samples
# Sr = Sar13 -> Sar17; Lro = Idr16 -> Idr17
# change individual names
rownames(dat.df.tot)[substr(rownames(dat.df.tot),1,2)=="Sr"] <- paste0("Sar17_",16:22)
rownames(dat.df.tot)[substr(rownames(dat.df.tot),1,3)=="Lro"] <- paste0("Idr17_",30:61)

# change population names
dat.df.tot$pop[dat.df.tot$pop=="Sr7"] <- "Sar17_15"
dat.df.tot$pop[dat.df.tot$pop=="Lro33"] <- "Idr17_29"
# remove unused population levels
dat.df.tot$pop <- droplevels(dat.df.tot$pop)

# transform back to genind
# column 1 = population
# ncode = 2 -> two characters for one genotype
dat.new.tot <- df2genind(dat.df.tot[,-1],ncode=2,pop=dat.df.tot[,1])
# transform to hierfstat
dat.hierf.tot <- genind2hierfstat(dat.new.tot)

# calculate Fst; in functions.R
# can take a minute or two
dist.tot <- calcFst(dat.hierf.tot)
write.csv(dist.tot,"Output/distancesFst_tot.csv")

# outlier data
dat.sel <- read.genepop("Input/lobster_1278ind_8snps_40pop.gen")
dat.df.sel <- genind2df(dat.sel)

# consistent naming
rownames(dat.df.sel)
rownames(dat.df.sel)[substr(rownames(dat.df.sel),1,2)=="SK"] <- paste0("Sky",substr(rownames(dat.df.sel[substr(rownames(dat.df.sel),1,2)=="SK",]),3,4))
rownames(dat.df.sel)[substr(rownames(dat.df.sel),1,2)=="TH"] <- paste0("The",substr(rownames(dat.df.sel[substr(rownames(dat.df.sel),1,2)=="TH",]),3,4))
rownames(dat.df.sel)[substr(rownames(dat.df.sel),1,2)=="AL"] <- paste0("Ale",substr(rownames(dat.df.sel[substr(rownames(dat.df.sel),1,2)=="AL",]),3,4))
rownames(dat.df.sel)[substr(rownames(dat.df.sel),1,2)=="CH"] <- paste0("Tor",substr(rownames(dat.df.sel[substr(rownames(dat.df.sel),1,2)=="CH",]),3,4))
rownames(dat.df.sel)[substr(rownames(dat.df.sel),1,2)=="Sr"] <- paste0("Sar17_",16:22)
rownames(dat.df.sel)[substr(rownames(dat.df.sel),1,3)=="Lro"] <- paste0("Idr17_",30:61)

# add population levels
levels(dat.df.sel$pop) <- c(levels(dat.df.sel$pop),"Sky49","The44","Ale89","Tor74")
dat.df.sel$pop[dat.df.sel$pop=="CH74"] <- "Tor74"
dat.df.sel$pop[dat.df.sel$pop=="SK49"] <- "Sky49"
dat.df.sel$pop[dat.df.sel$pop=="TH44"] <- "The44"
dat.df.sel$pop[dat.df.sel$pop=="AL89"] <- "Ale89"
dat.df.sel$pop[dat.df.sel$pop=="Sr7"] <- "Sar17_15"
dat.df.sel$pop[dat.df.sel$pop=="Lro33"] <- "Idr17_29"
dat.df.sel$pop <- droplevels(dat.df.sel$pop)

dat.new.sel <- df2genind(dat.df.sel[,-1],ncode=2,pop=dat.df.sel[,1])
dat.hierf.sel <- genind2hierfstat(dat.new.sel)

# can take a minute or two
dist.sel <- calcFst(dat.hierf.sel)
write.csv(dist.sel,"Output/distancesFst_sel.csv")

# neutral data
dat.neut <- read.genepop("Input/lobster_1278ind_71snps_40pop.gen")
dat.df.neut <- genind2df(dat.neut)

# consistent naming
rownames(dat.df.neut)[substr(rownames(dat.df.neut),1,2)=="SK"] <- paste0("Sky",substr(rownames(dat.df.neut[substr(rownames(dat.df.neut),1,2)=="SK",]),3,4))
rownames(dat.df.neut)[substr(rownames(dat.df.neut),1,2)=="TH"] <- paste0("The",substr(rownames(dat.df.neut[substr(rownames(dat.df.neut),1,2)=="TH",]),3,4))
rownames(dat.df.neut)[substr(rownames(dat.df.neut),1,2)=="AL"] <- paste0("Ale",substr(rownames(dat.df.neut[substr(rownames(dat.df.neut),1,2)=="AL",]),3,4))
rownames(dat.df.neut)[substr(rownames(dat.df.neut),1,2)=="CH"] <- paste0("Tor",substr(rownames(dat.df.neut[substr(rownames(dat.df.neut),1,2)=="CH",]),3,4))
rownames(dat.df.neut)[substr(rownames(dat.df.neut),1,2)=="Sr"] <- paste0("Sar17_",16:22)
rownames(dat.df.neut)[substr(rownames(dat.df.neut),1,3)=="Lro"] <- paste0("Idr17_",30:61)

levels(dat.df.neut$pop) <- c(levels(dat.df.neut$pop),"Sky49","The44","Ale89","Tor74")
dat.df.neut$pop[dat.df.neut$pop=="CH74"] <- "Tor74"
dat.df.neut$pop[dat.df.neut$pop=="SK49"] <- "Sky49"
dat.df.neut$pop[dat.df.neut$pop=="TH44"] <- "The44"
dat.df.neut$pop[dat.df.neut$pop=="AL89"] <- "Ale89"
dat.df.neut$pop[dat.df.neut$pop=="Sr7"] <- "Sar17_15"
dat.df.neut$pop[dat.df.neut$pop=="Lro33"] <- "Idr17_29"
dat.df.neut$pop <- droplevels(dat.df.neut$pop)

dat.new.neut <- df2genind(dat.df.neut[,-1],ncode=2,pop=dat.df.neut[,1])
dat.hierf.neut <- genind2hierfstat(dat.new.neut)

# can take a minute or two
dist.neut <- calcFst(dat.hierf.neut)
write.csv(dist.sel,"Output/distancesFst_neut.csv")