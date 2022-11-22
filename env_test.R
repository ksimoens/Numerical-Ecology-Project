library(vegan)
library(gridExtra)
library(tidyverse)

# read environmental data
env <- read.csv("Output/EnvMatrix.csv", row.names=1, header=T)
# only Atlantic sites
env.atl <- env[!(row.names(env) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]

# standardise variables
env.atl.norm <- as.data.frame(scale(env.atl, scale=T, center=T))

# function to calculate VIF from variable matrix
VIF_analysis <- function(x){
	# lm() requires df
	x <- as.data.frame(x)

	varname <- vector()
	Rsquared <- vector()
	VIF <- vector()

	for(i in 1:ncol(x)){
		# get variable name
		varname <- c(varname, colnames(x)[i])
		# make model that regresses variable i to all other variables
		mod <- lm(data=x[,-i], x[,i]~.)

		# calculate R2 of the model
		R2 <- summary(mod)$r.squared
		Rsquared <- c(Rsquared,R2)

		# calculate VIF
		VIF <- c(VIF,1/(1-R2))
	}
	output <- data.frame(variable=varname, Rsquared=Rsquared, VIF=VIF)
}

VIF <- VIF_analysis(env.atl.norm)
print(VIF)

# plot variable distributions
p_SST_m <- ggplot(env.atl.norm,aes(x=SST_mean)) + geom_histogram(col='black',fill='grey') + theme_bw() +
	theme(panel.grid.minor = element_blank()) + xlab("mean SST") 
p_SST_r <- ggplot(env.atl.norm,aes(x=SST_range)) + geom_histogram(col='black',fill='grey') + theme_bw() +
	theme(panel.grid.minor = element_blank()) + xlab("range SST") 

p_SAL_m <- ggplot(env.atl.norm,aes(x=SAL_mean)) + geom_histogram(col='black',fill='grey') + theme_bw() +
	theme(panel.grid.minor = element_blank()) + xlab("mean salinity") 
p_SAL_r <- ggplot(env.atl.norm,aes(x=SAL_range)) + geom_histogram(col='black',fill='grey') + theme_bw() +
	theme(panel.grid.minor = element_blank()) + xlab("range salinity") 

p_VEL_m <- ggplot(env.atl.norm,aes(x=VEL_mean)) + geom_histogram(col='black',fill='grey') + theme_bw() +
	theme(panel.grid.minor = element_blank()) + xlab("mean current velocity") 

p_PP_m <- ggplot(env.atl.norm,aes(x=PP_mean)) + geom_histogram(col='black',fill='grey') + theme_bw() +
	theme(panel.grid.minor = element_blank()) + xlab("mean productivity") 
p_PP_r <- ggplot(env.atl.norm,aes(x=PP_range)) + geom_histogram(col='black',fill='grey') + theme_bw() +
	theme(panel.grid.minor = element_blank()) + xlab("range productivity") 

p_BAT <- ggplot(env.atl.norm,aes(x=bathy)) + geom_histogram(col='black',fill='grey') + theme_bw() +
	theme(panel.grid.minor = element_blank()) + xlab("bathymetry") 

g <- arrangeGrob(p_SST_m, p_SST_r, p_SAL_m, p_SAL_r, p_VEL_m, p_PP_m, p_PP_r, p_BAT,nrow=4,ncol=2)

g %>% ggsave("normality.png",.,device='png',width=20,height=40,units='cm')



