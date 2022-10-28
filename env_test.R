library(vegan)
library(gridExtra)
library(tidyverse)

env <- read.csv("Output/EnvMatrix.csv", row.names=1, header=T)
env.atl <- env[!(row.names(env) %in% c("Laz","Tar","Sar","Ale","The","Tor","Sky")),]

env.atl.norm <- as.data.frame(scale(env.atl, scale=T, center=T))

VIF_analysis <- function(x){
	x <- as.data.frame(x)

	varname <- vector()
	Rsquared <- vector()
	VIF <- vector()

	for(i in 1:ncol(x)){
		varname <- c(varname, colnames(x)[i])
		mod <- lm(data=x[,-i], x[,i]~.)

		R2 <- summary(mod)$r.squared
		Rsquared <- c(Rsquared,R2)

		VIF <- c(VIF,1/(1-R2))
	}
	output <- data.frame(variable=varname, Rsquared=Rsquared, VIF=VIF)
}

VIF <- VIF_analysis(env.atl.norm)
print(VIF)


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

g <- arrangeGrob(p_SST_m, p_SST_r, p_SAL_m, p_SAL_r, p_VEL_m, p_PP_m, p_PP_r, p_BAT,ncol=4)

g %>% ggsave("normality.png",.,device='png',width=40,height=20,units='cm')



